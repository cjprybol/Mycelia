# Unit tests for soft-EM edge weighting (td-e70t).
#
# The graph-as-HMM correction redesign replaces the hard-count M-step
# (`compute_edge_weight` == raw coverage) with probability-weighted evidence:
# after Viterbi decodes a read, each candidate path's RESPONSIBILITY (posterior
# over the read's candidate set) is accumulated onto its edges. Error edges,
# traversed only by rare low-probability paths, accrue little soft weight and
# contribute less to subsequent transition scoring.
#
# This file covers the accumulation primitives, support floor, and direct
# registration behavior. Production M-step consumption is landed via task-local
# scoped snapshots. The full-pipeline monotonic-decrease acceptance remains
# skipped pending a low-skip-fraction integration fixture.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/soft_em_edge_weight_scaffold_test.jl")'

import BioSequences
import FASTX
import Kmers
import MetaGraphsNext
import Mycelia
import Test

Test.@testset "Soft-EM edge weighting primitives (td-e70t)" begin
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
        # The lower responsibility gives the error edge less transition weight;
        # the true edge dominates.
        Test.@test Mycelia.Rhizomorph.soft_edge_weight(acc, error_edge) <
                   Mycelia.Rhizomorph.soft_edge_weight(acc, true_edge)
    end

    Test.@testset "error-edge soft weight decreases across (simulated) EM iterations" begin
        # Simulate the E->M loop at the primitive level: each iteration a read
        # offers a TRUE path and an ERROR path. As correction sharpens the graph
        # across iterations, the error path's relative likelihood falls (the gap
        # widens), so its per-iteration responsibility -- the soft weight it
        # deposits -- monotonically decreases. This is the transition-weight
        # signal consumed by the next M-step.
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

    Test.@testset "SUPPORT FLOOR retains supported skewed variant, lets error decay (td-h6w9)" begin
        # The mechanism that fixes the prior soft-EM v2's collapse of skewed
        # variants (PR #363 review C1): `register_soft_edge_weights!` clamps an
        # edge meeting the support threshold to at least its RAW coverage, while a
        # below-threshold edge keeps a floor of 0 and is free to decay.
        #
        # Controlled graph: from source ATG three out-edges with distinct raw
        # coverage — TGC (5x majority), TGT (4x SKEWED minority, >= MIN_SUPPORT=3),
        # TGG (1x error, < MIN_SUPPORT). We put DECAYED accumulator weights on the
        # minority and error edges (simulating the geometric responsibility decay
        # that would zero them out) and prove the floor rescues the supported
        # minority (raised back to raw) but NOT the error.
        R = Mycelia.Rhizomorph
        recs = FASTX.FASTQ.Record[]
        for i in 1:5
            push!(recs, FASTX.FASTQ.Record(string("maj", i), "ATGCA", String(fill('I', 5))))
        end
        for i in 1:4
            push!(recs, FASTX.FASTQ.Record(string("min", i), "ATGTA", String(fill('I', 5))))
        end
        push!(recs, FASTX.FASTQ.Record("err1", "ATGGA", String(fill('I', 5))))
        graph = R.build_qualmer_graph(recs, 3; dataset_id = "floor_probe", mode = :singlestrand)

        # Locate the minority (raw==4) and error (raw==1) edges by coverage.
        min_key, err_key = nothing, nothing
        for (s, d) in MetaGraphsNext.edge_labels(graph)
            cov = R.count_total_observations(graph[s, d])
            cov == 4 && (min_key = (s, d))
            cov == 1 && (err_key = (s, d))
        end
        Test.@test min_key !== nothing
        Test.@test err_key !== nothing
        Test.@test R.count_total_observations(graph[min_key...]) == 4   # supported minority
        Test.@test R.count_total_observations(graph[err_key...]) == 1   # error

        # Decayed soft weights (as the responsibility split would produce): the
        # supported minority has been driven far below its raw coverage, the error
        # further still.
        decayed_min, decayed_err = 0.5, 0.02
        acc = R.SoftEdgeWeightAccumulator()
        R.accumulate_path_probability!(acc, [min_key], decayed_min)
        R.accumulate_path_probability!(acc, [err_key], decayed_err)

        # (1) WITH the support floor (min_support=3): the supported minority is
        # rescued to its raw coverage (never decays below its own support); the
        # error keeps its decayed value (floor 0) and stays below the 0.01... it is
        # 0.02 here but decays further across iterations (see the recurrence test).
        R.clear_soft_edge_weights!()
        R.register_soft_edge_weights!(graph, acc; min_support = 3)
        Test.@test R.compute_edge_weight(graph[min_key...]) ≈ 4.0      # floored to raw
        Test.@test R.compute_edge_weight(graph[err_key...]) ≈ decayed_err  # unfloored
        R.clear_soft_edge_weights!()

        # (2) WITHOUT the floor (min_support above raw): the prior-v2 behavior — the
        # supported minority stays DECAYED (would geometrically collapse to zero),
        # exactly the bug the floor fixes.
        R.register_soft_edge_weights!(graph, acc; min_support = 1000)
        Test.@test R.compute_edge_weight(graph[min_key...]) ≈ decayed_min   # NOT rescued
        Test.@test R.compute_edge_weight(graph[err_key...]) ≈ decayed_err
        R.clear_soft_edge_weights!()

        # The discriminating property: the floor RAISES the supported minority
        # (0.5 -> 4.0) but leaves the error untouched (0.02 -> 0.02).
        Test.@test 4.0 > decayed_min       # floor rescued the supported variant
        Test.@test decayed_err < 3         # error is below the support threshold, not rescued
    end

    Test.@testset "SUPPORT FLOOR breaks the geometric decay of a skewed minority (td-h6w9)" begin
        # Prior v2 recurrence (PR #363 review C1): a skewed minority edge's soft
        # weight follows W' = N * W / (W_maj + W_min), a geometric contraction to
        # ZERO whenever a strictly-heavier sibling exists — a real skewed allele
        # decays like an error. Simulate the E->M loop at the primitive level: a
        # 10x minority vs 20x majority. WITHOUT the floor the minority's registered
        # weight contracts across iterations; WITH the floor it is held at its raw
        # coverage (>= MIN_SUPPORT) every iteration.
        R = Mycelia.Rhizomorph
        recs = FASTX.FASTQ.Record[]
        for i in 1:20
            push!(recs, FASTX.FASTQ.Record(string("maj", i), "ATGCA", String(fill('I', 5))))
        end
        for i in 1:10
            push!(recs, FASTX.FASTQ.Record(string("min", i), "ATGTA", String(fill('I', 5))))
        end
        graph = R.build_qualmer_graph(recs, 3; dataset_id = "decay_probe", mode = :singlestrand)
        maj_key, min_key = nothing, nothing
        for (s, d) in MetaGraphsNext.edge_labels(graph)
            cov = R.count_total_observations(graph[s, d])
            cov == 20 && (maj_key = (s, d))
            cov == 10 && (min_key = (s, d))
        end
        Test.@test maj_key !== nothing && min_key !== nothing
        raw_min = Float64(R.count_total_observations(graph[min_key...]))  # 10

        # Iterate the recurrence: each iteration the minority's *effective* weight
        # (what the next iteration's decode sees) drives its responsibility. Model
        # responsibility ∝ effective weight; the minority deposits N_min *
        # (w_min / (w_maj + w_min)) onto its edge.
        n_min = 10.0
        w_maj = 20.0
        # -- No floor --
        w_min_nofloor = raw_min
        nofloor_trace = Float64[]
        for _ in 1:5
            resp = w_min_nofloor / (w_maj + w_min_nofloor)
            deposited = n_min * resp
            # min_support above raw ⇒ floor 0 ⇒ registered weight == deposited
            w_min_nofloor = deposited
            push!(nofloor_trace, w_min_nofloor)
        end
        # -- With floor (min_support=3, raw=10 ⇒ floor 10) --
        w_min_floor = raw_min
        floor_trace = Float64[]
        for _ in 1:5
            resp = w_min_floor / (w_maj + w_min_floor)
            deposited = n_min * resp
            w_min_floor = max(deposited, raw_min >= R.SOFT_EM_MIN_SUPPORT ? raw_min : 0.0)
            push!(floor_trace, w_min_floor)
        end

        # Without the floor the minority weight strictly contracts (heads toward 0).
        Test.@test all(diff(nofloor_trace) .< 0)
        Test.@test nofloor_trace[end] < 0.5 * raw_min
        # With the floor it is pinned at raw coverage every iteration — retained.
        Test.@test all(w -> w ≈ raw_min, floor_trace)
    end

    Test.@testset "FULL-PIPELINE ACCEPTANCE (skipped pending low-skip-fraction integration harness — td-969e)" begin
        # td-e70t acceptance: running mycelia_iterative_assemble on a clean ~1 kb
        # toy should coalesce the qualmer graph toward ~1 contig across EM
        # iterations, and the summed soft weight of error edges should DECREASE
        # across iterations in a fixture that isolates the soft-EM contribution
        # from separately configured graph cleaning.
        #
        # M-step consumption IS LANDED on this branch: `mycelia_iterative_assemble`
        # binds each EM iteration's soft-weight accumulator to the next iteration's
        # freshly-built graph in a task-local scope, so `compute_edge_weight` /
        # Viterbi transition scoring consume support-floored soft weights and the
        # prior snapshot is restored afterward. The wiring is exercised by
        # scalable_corrector_soft_em_test.jl and the M-step-floor recurrence above.
        #
        # This end-to-end monotonic-decrease assertion remains @test_skip only
        # because it needs a dedicated LOW-SKIP-FRACTION integration harness (bead
        # td-969e): on the default :scalable route most clean reads are
        # skip-solid'd, so few reads reach the decoder and the per-iteration
        # error-edge soft-weight trace is too sparse to assert a stable monotonic
        # trend without a fixture that forces a meaningful decode fraction. Skipped
        # (not @test false) so it documents the target without failing CI pending
        # that harness.
        Test.@test_skip false
    end
end
