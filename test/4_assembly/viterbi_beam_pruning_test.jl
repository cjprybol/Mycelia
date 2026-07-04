# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_beam_pruning_test.jl")'
# ```
#
# Regression for td-63qy: the exact Viterbi corrector (`_viterbi_correct_observation`)
# never beam-capped its (vertex, strand) frontier, so the reachable-state set grew
# ~unboundedly with read-length depth on real branchy graphs — 21 billion
# allocations and a hard crash on a 48 kb phage. `ViterbiCorrectionConfig` now
# carries `beam_width` (default `typemax(Int)` = exact/unbounded); a finite value
# caps the frontier to the top-K states by score each depth.

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

function beam_decoded_label_strings(
    result::Mycelia.ViterbiCorrectionResult
)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

Test.@testset "Viterbi corrector beam pruning (td-63qy)" begin
    # A simple linear graph: the frontier is small, so the exact and bounded
    # decoders must agree, and pruning must never fire.
    linear_records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
    linear_graph = Mycelia.Rhizomorph.build_kmer_graph(
        linear_records, 3; dataset_id="beam_linear", mode=:singlestrand
    )
    linear_observed = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGA"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
    ]

    Test.@testset "default beam_width is unbounded (exact, no regression)" begin
        result = Mycelia.correct_observations(linear_graph, [linear_observed])
        # The exact answer from the existing B8 fixture — must be byte-identical.
        Test.@test beam_decoded_label_strings(result) == ["ATG", "TGC", "GCG", "CGT"]
        # No pruning on the default (unbounded) path.
        Test.@test get(only(result.paths).diagnostics, :beam_pruned, 0) == 0
    end

    Test.@testset "a finite beam still recovers the linear correction" begin
        config = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, beam_width=1)
        result = Mycelia.correct_observations(linear_graph, [linear_observed]; config=config)
        # Anti-crash guarantee: a bounded corrector completes and returns a
        # valid, complete correction on a linear graph.
        Test.@test beam_decoded_label_strings(result) == ["ATG", "TGC", "GCG", "CGT"]
    end

    # A bubble graph (two paths sharing a start/end k-mer) forces the frontier
    # above one state at the branch, so a beam_width=1 corrector MUST prune —
    # and must still complete rather than explode.
    bubble_records = [
        FASTX.FASTA.Record("path_a", BioSequences.dna"ATGCGTA"),
        FASTX.FASTA.Record("path_b", BioSequences.dna"ATGAGTA"),
    ]
    bubble_graph = Mycelia.Rhizomorph.build_kmer_graph(
        bubble_records, 3; dataset_id="beam_bubble", mode=:singlestrand
    )
    bubble_observed = [
        Kmers.DNAKmer{3}("ATG"),
        Kmers.DNAKmer{3}("TGC"),
        Kmers.DNAKmer{3}("GCG"),
        Kmers.DNAKmer{3}("CGT"),
        Kmers.DNAKmer{3}("GTA"),
    ]

    Test.@testset "beam pruning activates on a branch and still completes" begin
        config = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, beam_width=1)
        result = Mycelia.correct_observations(bubble_graph, [bubble_observed]; config=config)
        decoded = only(result.corrected_observations)
        # Completes with a non-empty correction (did not explode / crash).
        Test.@test !isempty(decoded)
        # Pruning fired at least once (the branch pushed the frontier past 1).
        Test.@test get(only(result.paths).diagnostics, :beam_pruned, 0) >= 1
    end

    Test.@testset "unbounded decoder on the bubble is at least as good as beamed" begin
        exact = Mycelia.correct_observations(bubble_graph, [bubble_observed])
        beamed_config = Mycelia.ViterbiCorrectionConfig(alphabet=:DNA, beam_width=1)
        beamed = Mycelia.correct_observations(bubble_graph, [bubble_observed]; config=beamed_config)
        # Beam search is an approximation: the exact score is an upper bound.
        Test.@test only(exact.paths).score >= only(beamed.paths).score - 1e-9
    end

    Test.@testset "beam_width must be positive" begin
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_width=0)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_width=-4)
    end
end
