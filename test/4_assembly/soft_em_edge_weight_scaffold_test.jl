# Soft-EM edge weighting tests (td-e70t) — v2 competing-paths ACTIVE.
#
# The graph-as-HMM correction redesign replaces the hard-count M-step
# (`compute_edge_weight` == raw coverage) with probability-weighted evidence:
# each competing candidate path's RESPONSIBILITY (posterior over the read's
# candidate set) is accumulated onto its edges. Error edges, traversed only by
# low-probability paths, accrue little soft weight and decay below the
# emergent-cleaning gate WITHOUT tip-clipping.
#
# This file exercises the accumulation PRIMITIVES (softmax responsibility,
# accumulator, single-path hook) at the unit level AND the v2 activation: the
# final testset drives the competing-paths E-step + the M-step consumption
# feedback loop end-to-end (register soft weights -> compute_edge_weight decays ->
# next iteration's responsibility sharpens) and asserts a tracked error edge's
# soft weight strictly decreases across EM iterations below the 0.01 gate.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/soft_em_edge_weight_scaffold_test.jl")'

import BioSequences
import FASTX
import Kmers
import Mycelia
import Random
import Test

Test.@testset "Soft-EM edge weighting scaffold (td-e70t)" begin
    Test.@testset "path_responsibility softmax normalization" begin
        # A single candidate path owns all the responsibility.
        Test.@test Mycelia.Rhizomorph.path_responsibility(-3.2, Float64[-3.2]) ≈ 1.0
        # Competing paths split responsibility by relative likelihood; the sum is 1.
        logps = Float64[-1.0, -2.0, -5.0]
        rs = [Mycelia.Rhizomorph.path_responsibility(lp, logps) for lp in logps]
        Test.@test sum(rs) ≈ 1.0
        Test.@test rs[1] > rs[2] > rs[3]          # higher logp => more responsibility
        # A far-lower-likelihood (error) path gets vanishing responsibility.
        Test.@test Mycelia.Rhizomorph.path_responsibility(-20.0, Float64[-1.0, -20.0]) < 1e-6
    end

    Test.@testset "accumulator sums path responsibilities onto edges" begin
        acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        true_edge = ("ATG", "TGC")
        error_edge = ("ATG", "TGA")
        Mycelia.Rhizomorph.accumulate_path_probability!(acc, [true_edge], 0.9)
        Mycelia.Rhizomorph.accumulate_path_probability!(acc, [error_edge], 0.1)
        Mycelia.Rhizomorph.accumulate_path_probability!(acc, [true_edge], 0.95)

        Test.@test Mycelia.Rhizomorph.soft_edge_weight(acc, true_edge) ≈ 1.85
        Test.@test Mycelia.Rhizomorph.soft_edge_weight(acc, error_edge) ≈ 0.1
        # Unseen edge falls back to the prior.
        Test.@test Mycelia.Rhizomorph.soft_edge_weight(acc, ("X", "Y"); prior = 0.0) == 0.0
        # The error edge is already far below the emergent-cleaning gate (0.01 per
        # read of responsibility); the true edge dominates.
        Test.@test Mycelia.Rhizomorph.soft_edge_weight(acc, error_edge) <
                   Mycelia.Rhizomorph.soft_edge_weight(acc, true_edge)
    end

    Test.@testset "error-edge soft weight decreases across (simulated) EM iterations" begin
        # Simulate the E->M loop at the primitive level: each iteration a read
        # offers a TRUE path and an ERROR path. As correction sharpens the graph
        # across iterations, the error path's relative likelihood falls (the gap
        # widens), so its per-iteration responsibility -- the soft weight it
        # deposits -- monotonically decreases. This is the decay mechanism that,
        # once the M-step consumes it, drives emergent cleaning.
        true_edge = ("ATG", "TGC")
        error_edge = ("ATG", "TGA")
        true_logp = -1.0
        error_contributions = Float64[]
        for iteration in 1:5
            error_logp = true_logp - Float64(iteration)   # gap widens each iteration
            candidate_logps = Float64[true_logp, error_logp]
            acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
            Mycelia.Rhizomorph.accumulate_path_probability!(
                acc, [true_edge],
                Mycelia.Rhizomorph.path_responsibility(true_logp, candidate_logps))
            Mycelia.Rhizomorph.accumulate_path_probability!(
                acc, [error_edge],
                Mycelia.Rhizomorph.path_responsibility(error_logp, candidate_logps))
            push!(error_contributions, Mycelia.Rhizomorph.soft_edge_weight(acc, error_edge))
        end
        # Strictly monotone decrease of the error edge's accumulated soft weight.
        Test.@test all(diff(error_contributions) .< 0)
        Test.@test error_contributions[end] < error_contributions[1]
    end

    Test.@testset "hook accumulates real decoded paths (responsibility 1 for single path)" begin
        records = FASTX.FASTQ.Record[]
        for index in 1:20
            push!(records, FASTX.FASTQ.Record(
                "truth_$index", "ATGCA", String(fill('I', 5))))
        end
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; dataset_id = "soft_em_scaffold", mode = :singlestrand)

        observation = [
            Mycelia.QualityObservation(Kmers.DNAKmer{3}("ATG"), UInt8[40, 40, 40]),
            Mycelia.QualityObservation(Kmers.DNAKmer{3}("TGC"), UInt8[40, 40, 40]),
            Mycelia.QualityObservation(Kmers.DNAKmer{3}("GCA"), UInt8[40, 40, 40])
        ]
        result = Mycelia.correct_observations(graph, [observation])

        edges = Mycelia._decoded_path_edges(only(result.paths))
        Test.@test !isempty(edges)

        acc = Mycelia.Rhizomorph.SoftEdgeWeightAccumulator()
        Mycelia.accumulate_soft_em_edge_weights!(acc, result.paths)
        # A single argmax path per read => responsibility 1.0 on every traversed
        # edge (soft weight == hard count in the degenerate single-candidate case).
        for edge in edges
            Test.@test Mycelia.Rhizomorph.soft_edge_weight(acc, edge) ≈ 1.0
        end
    end

    Test.@testset "ERROR-EDGE DECAY across EM iterations (v2 competing paths, td-e70t)" begin
        # td-e70t acceptance, now ACTIVE (v2). Controlled fixture: a random backbone
        # at high coverage plus ONE coverage-1 error read. The single substitution
        # creates coverage-1 error k-mers/edges while every consensus edge is
        # high-coverage. The v2 competing-paths E-step
        # (`accumulate_competing_paths!`) gives the error branch a small
        # responsibility (a genuine split — v1 was always 1.0), and REGISTERING that
        # soft weight (the M-step) makes `compute_edge_weight` return the decayed
        # value, so the next iteration's transition scoring sharpens the split
        # further — a positive feedback loop. We isolate the loop by re-accumulating
        # on the SAME reads across iterations (stable edge identity) and assert the
        # tracked error edge's soft weight decreases strictly and falls below the
        # 0.01 emergent-cleaning gate, WITHOUT any explicit tip/bubble removal.
        R = Mycelia.Rhizomorph
        bases = ['A', 'C', 'G', 'T']
        rng = Random.MersenneTwister(7)
        L, cov, readlen, errpos, kk = 300, 25, 100, 150, 13
        backbone = collect(join(rand(rng, bases, L)))
        reads = FASTX.FASTQ.Record[]
        qual = String(fill('I', readlen))
        for i in 1:cov
            s = rand(rng, 1:(L - readlen + 1))
            push!(reads, FASTX.FASTQ.Record("c$i", String(backbone[s:(s + readlen - 1)]), qual))
        end
        errseq = copy(backbone)
        errseq[errpos] = rand(rng, filter(!=(errseq[errpos]), bases))
        s0 = errpos - readlen ÷ 2
        push!(reads, FASTX.FASTQ.Record("err", String(errseq[s0:(s0 + readlen - 1)]), qual))

        graph = R.build_qualmer_graph(reads, kk; mode = :canonical)
        err_edges = [(a, b) for (a, b) in R.MetaGraphsNext.edge_labels(graph)
                     if R.count_total_observations(graph[a, b]) == 1]
        Test.@test !isempty(err_edges)

        # Iteration 1 E-step (raw weights).
        acc1 = R.SoftEdgeWeightAccumulator()
        for rd in reads
            Mycelia.accumulate_competing_paths!(acc1, rd, graph, kk; graph_mode = :canonical)
        end
        # v2 genuinely SPLITS: at least one error edge received a fractional (<1)
        # responsibility because a consensus alternative outscored the observed
        # error branch. (v1's single-path decode made every responsibility 1.0.)
        competed = [(e, acc1.weights[e]) for e in err_edges
                    if haskey(acc1.weights, e) && acc1.weights[e] < 1.0]
        Test.@test !isempty(competed)
        tracked = last(sort(competed; by = x -> x[2]))[1]   # largest still-competed edge

        # Feedback iterations: register the previous accumulator (M-step), then
        # re-accumulate; `compute_edge_weight` now reads the decayed soft weight.
        trace = Float64[acc1.weights[tracked]]
        prev = acc1
        for _ in 1:4
            acc = R.SoftEdgeWeightAccumulator()
            try
                R.register_soft_edge_weights!(graph, prev)
                for rd in reads
                    Mycelia.accumulate_competing_paths!(acc, rd, graph, kk; graph_mode = :canonical)
                end
            finally
                R.clear_soft_edge_weights!()   # never leak the process-global registry
            end
            push!(trace, get(acc.weights, tracked, 0.0))
            prev = acc
        end
        # Strictly decreasing across EM iterations, and below the 0.01 gate.
        Test.@test all(diff(trace) .< 0)
        Test.@test trace[end] < 0.01
        # Registry cleared ⇒ the raw coverage count is restored byte-for-byte.
        Test.@test R.compute_edge_weight(graph[tracked...]) ==
                   R.count_total_observations(graph[tracked...])
    end
end
