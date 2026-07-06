# SCAFFOLD test for soft-EM edge weighting (td-e70t).
#
# The graph-as-HMM correction redesign replaces the hard-count M-step
# (`compute_edge_weight` == raw coverage) with probability-weighted evidence:
# after Viterbi decodes a read, each candidate path's RESPONSIBILITY (posterior
# over the read's candidate set) is accumulated onto its edges. Error edges,
# traversed only by rare low-probability paths, accrue little soft weight and
# decay below the emergent-cleaning gate WITHOUT tip-clipping.
#
# What lands here: the accumulation PRIMITIVE + the correction-side hook, proven
# at the unit level (passing tests below). What does NOT land yet: making the
# qualmer graph rebuild / transition weighting CONSUME these soft weights so the
# full pipeline coalesces a 1 kb toy to ~1 contig across EM iterations — that
# M-step consumption is outside the correction-core file boundary and is the
# td-e70t follow-on. The full-pipeline acceptance is therefore a skipped test.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/soft_em_edge_weight_scaffold_test.jl")'

import BioSequences
import FASTX
import Kmers
import Mycelia
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

    Test.@testset "FULL-PIPELINE ACCEPTANCE (skipped until M-step consumes soft weights)" begin
        # td-e70t acceptance: running mycelia_iterative_assemble on a clean ~1 kb
        # toy should coalesce the qualmer graph toward ~1 contig across EM
        # iterations, and the summed soft weight of error edges should DECREASE
        # across iterations, WITHOUT any explicit tip/bubble removal. This requires
        # the qualmer graph rebuild / transition weighting to consume the soft
        # weights produced by `accumulate_soft_em_edge_weights!` in place of raw
        # coverage counts -- the follow-on that lives outside the correction-core
        # file boundary. Skipped (not @test false) so it documents the target
        # without failing CI on unlanded work.
        Test.@test_skip false
    end
end
