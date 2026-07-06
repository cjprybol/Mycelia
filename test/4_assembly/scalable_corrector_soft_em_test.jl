# Stage 2 (td-e70t): soft-EM v1 edge-memory wiring. After each successful Viterbi
# decode the ML path's edges accumulate responsibility (1.0 for the single argmax
# path) into a SoftEdgeWeightAccumulator; the NEXT iteration's graph is re-weighted
# from that accumulator so the M-step consumes probability-weighted evidence
# instead of raw coverage counts. Rare (error) edges — traversed by few paths —
# accrue little soft weight and are non-increasing across iterations, which is the
# "probability memory" that makes emergent cleaning possible.
#
# CAVEAT (v1): with a single argmax path per read, responsibility is always 1.0
# (soft weight == hard count). Full emergent coalescence to ~1 contig needs the v2
# competing-paths E-step (several candidate paths per read); this test asserts the
# v1 mechanism (accumulation + consumption + non-increasing rare edges), not full
# coalescence.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/scalable_corrector_soft_em_test.jl")'

import Test
import Mycelia
import FASTX
import Random

const _BASES = ['A', 'C', 'G', 'T']

function _toy_fastq_records(rng; reflen = 1000, n_reads = 120, readlen = 80, err = 0.0)
    ref = join(rand(rng, _BASES, reflen))
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        for j in 1:readlen
            err > 0 && rand(rng) < err && (seq[j] = rand(rng, filter(!=(seq[j]), _BASES)))
        end
        push!(records, FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
    end
    return records
end

Test.@testset "scalable corrector soft-EM v1 (td-e70t)" begin
    R = Mycelia.Rhizomorph

    Test.@testset "registry consumption + byte-identical revert" begin
        # A soft-registered edge returns its soft weight from compute_edge_weight;
        # after clear!, the raw coverage count is restored EXACTLY (the invariant
        # that keeps every non-soft-EM path unchanged).
        reads = _toy_fastq_records(Random.MersenneTwister(3); n_reads = 60, err = 0.0)
        graph = R.build_qualmer_graph(reads, 13; mode = :canonical)
        edge_labels = collect(Mycelia.Rhizomorph.MetaGraphsNext.edge_labels(graph))
        Test.@test !isempty(edge_labels)
        (src, dst) = first(edge_labels)
        edge_data = graph[src, dst]
        raw = R.compute_edge_weight(edge_data)
        Test.@test raw >= 1.0

        acc = R.SoftEdgeWeightAccumulator()
        R.accumulate_path_probability!(acc, [(src, dst)], 0.37)
        R.register_soft_edge_weights!(graph, acc)
        Test.@test R.compute_edge_weight(edge_data) == 0.37   # soft weight consumed
        # An edge NOT in the accumulator keeps its raw count.
        if length(edge_labels) > 1
            (s2, d2) = edge_labels[2]
            Test.@test R.compute_edge_weight(graph[s2, d2]) ==
                       R.count_total_observations(graph[s2, d2])
        end
        R.clear_soft_edge_weights!()
        Test.@test R.compute_edge_weight(edge_data) == raw    # byte-identical revert
    end

    Test.@testset "decode pass populates the accumulator (responsibility 1.0)" begin
        reads = _toy_fastq_records(Random.MersenneTwister(5); n_reads = 80, err = 0.0)
        graph = R.build_qualmer_graph(reads, 13; mode = :canonical)
        acc = R.SoftEdgeWeightAccumulator()
        _updated, _n,
        _skip = Mycelia.improve_read_set_likelihood(
            reads, graph, 13; graph_mode = :canonical, soft_weights = acc)
        # The E-step accumulated ML-path edges: non-empty, and every soft weight is
        # a sum of unit (1.0) responsibilities ⇒ integer-valued and >= 1.
        Test.@test !isempty(acc.weights)
        Test.@test all(w -> w >= 1.0 && isapprox(w, round(w)), values(acc.weights))
    end

    Test.@testset "rare (error) edge soft weight is non-increasing across iterations" begin
        # Two manual EM iterations. Iteration 2's graph is re-weighted from
        # iteration 1's accumulator (the M-step consumption), then decoded on the
        # CORRECTED reads. The probability-memory / decay property is AGGREGATE in
        # v1 (responsibility is always 1.0, so a single edge's soft weight is just
        # its ML-path coverage and can rise if a read is corrected onto it): the
        # rare-edge (error) POPULATION shrinks. We assert both the count of rare
        # edges and their total tail mass are non-increasing.
        #
        # v1 CAVEAT: full coalescence to ~1 contig (and strict per-edge decay)
        # needs the v2 competing-paths E-step; this asserts the v1 aggregate trend.
        reads = _toy_fastq_records(Random.MersenneTwister(9); n_reads = 150, err = 0.008)
        k = 13
        rare_threshold = 2.0   # <=2 supporting ML paths ⇒ error-like edge

        graph1 = R.build_qualmer_graph(reads, k; mode = :canonical)
        acc1 = R.SoftEdgeWeightAccumulator()
        corrected1, _n1,
        _s1 = Mycelia.improve_read_set_likelihood(
            reads, graph1, k; graph_mode = :canonical, soft_weights = acc1)

        graph2 = R.build_qualmer_graph(corrected1, k; mode = :canonical)
        R.register_soft_edge_weights!(graph2, acc1)   # consume iter-1 soft memory
        acc2 = R.SoftEdgeWeightAccumulator()
        _corrected2, _n2,
        _s2 = Mycelia.improve_read_set_likelihood(
            corrected1, graph2, k; graph_mode = :canonical, soft_weights = acc2)
        R.clear_soft_edge_weights!()

        Test.@test !isempty(acc1.weights)
        # The distribution is right-skewed: consensus edges accrue high weight,
        # error edges little (the probability-memory signal).
        Test.@test maximum(values(acc1.weights)) > rare_threshold

        rare_count1 = count(w -> w <= rare_threshold, values(acc1.weights))
        rare_count2 = count(w -> w <= rare_threshold, values(acc2.weights))
        tail_mass1 = sum((w for w in values(acc1.weights) if w <= rare_threshold); init = 0.0)
        tail_mass2 = sum((w for w in values(acc2.weights) if w <= rare_threshold); init = 0.0)

        Test.@test rare_count1 > 0
        # Error-edge population is non-increasing (aggregate decay).
        Test.@test rare_count2 <= rare_count1
        Test.@test tail_mass2 <= tail_mass1 + 1e-9
    end
end
