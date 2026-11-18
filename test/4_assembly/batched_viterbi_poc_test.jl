# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/batched_viterbi_poc_test.jl")'
# ```
#
# Equivalence-oracle acceptance gate for the batched array-frontier corrector PoC
# (bead td-qoo3; ADR docs/design/2026-07-06-gpu-simd-corrector-acceleration.md).
#
# The batched decoder reformulates the scalar per-read decoder
# (`_viterbi_correct_observation`) into the layout a SIMD/GPU kernel can consume:
# a dense array frontier, a depth-outer / read-inner loop, length-binned batches,
# and a graph structure-of-arrays built once per batch. This test is the ADR's
# acceptance criterion (§Prototype scope item 3): on toy graphs + reads (including
# reads carrying an injected error — the "hard window" the residual decode
# targets), the batched decoder MUST return byte-identical corrected paths and
# log-probability scores to the scalar `correct_observations`, per read. It is
# also the regression guard for the future Phase A (SIMD) and Phase B (GPU) kernels.

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

if !isdefined(Main, :test_throws_message)
    include(joinpath(dirname(@__DIR__), "test_helpers.jl"))
end

# Corrected label strings for read `index` of a ViterbiCorrectionResult.
function poc_decoded_label_strings(
    result::Mycelia.ViterbiCorrectionResult, index::Integer
)::Vector{String}
    return [string(label) for label in result.corrected_observations[index]]
end

Test.@testset "Batched array-frontier Viterbi corrector PoC (td-qoo3)" begin
    Test.@testset "n-gram (TEXT, singlestrand) batch — equivalence oracle" begin
        # Reference "ABCDEFGH" -> 3-grams ABC,BCD,CDE,DEF,EFG,FGH (a linear chain).
        graph = Mycelia.Rhizomorph.build_ngram_graph(
            ["ABCDEFGH"], 3; dataset_id="batched_poc_ngram"
        )

        # A batch mixing lengths (exercises length-binning) and clean vs. error
        # reads (the error reads are the residual "hard windows" the decode fixes):
        #   o1: len-4, substitution error "BCX" at position 2 -> corrects to BCD
        #   o2: len-3, clean
        #   o3: len-4, substitution error "CDX" at position 3 -> corrects to CDE
        #   o4: len-3, error "EFZ" -> corrects to EFG
        #   o5: len-4, clean tail
        observations = [
            ["ABC", "BCX", "CDE", "DEF"],
            ["DEF", "EFG", "FGH"],
            ["ABC", "BCD", "CDX", "DEF"],
            ["DEF", "EFZ", "FGH"],
            ["CDE", "DEF", "EFG", "FGH"],
        ]

        oracle = Mycelia.batched_equivalence_oracle(graph, observations)

        # THE acceptance gate: batched == scalar, byte-identical, every read.
        Test.@test oracle.passes
        Test.@test all(oracle.per_read)
        Test.@test length(oracle.per_read) == length(observations)

        # Sanity: the error reads were actually corrected onto the reference chain.
        Test.@test poc_decoded_label_strings(oracle.batched, 1) ==
            ["ABC", "BCD", "CDE", "DEF"]
        Test.@test poc_decoded_label_strings(oracle.batched, 3) ==
            ["ABC", "BCD", "CDE", "DEF"]
        Test.@test poc_decoded_label_strings(oracle.batched, 4) ==
            ["DEF", "EFG", "FGH"]

        # Scores are compared with `===` inside the oracle (bit-identical Float64);
        # assert the corrected labels match the scalar decoder read-for-read here too.
        for index in 1:length(observations)
            Test.@test poc_decoded_label_strings(oracle.batched, index) ==
                poc_decoded_label_strings(oracle.scalar, index)
            Test.@test oracle.batched.paths[index].score ===
                oracle.scalar.paths[index].score
        end

        # The batch really did length-bin (len-3 and len-4 present) and ran the
        # batched interface.
        Test.@test oracle.batched.diagnostics[:interface] == :metagraphs_next_batched
        Test.@test oracle.batched.diagnostics[:length_bins] == [3, 4]
        Test.@test oracle.batched.diagnostics[:strand_mode] == :singlestrand
    end

    Test.@testset "DNA k-mer (singlestrand) batch — equivalence oracle" begin
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGTAC")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id="batched_poc_kmer", mode=:singlestrand
        )
        # "TGA" is an injected substitution for the true "TGC"; the ML path
        # corrects it back onto the reference walk.
        observations = [
            [
                Kmers.DNAKmer{3}("ATG"),
                Kmers.DNAKmer{3}("TGA"),
                Kmers.DNAKmer{3}("GCG"),
                Kmers.DNAKmer{3}("CGT"),
            ],
            [
                Kmers.DNAKmer{3}("GCG"),
                Kmers.DNAKmer{3}("CGT"),
                Kmers.DNAKmer{3}("GTA"),
            ],
        ]

        oracle = Mycelia.batched_equivalence_oracle(graph, observations)

        Test.@test oracle.passes
        Test.@test all(oracle.per_read)
        Test.@test poc_decoded_label_strings(oracle.batched, 1) ==
            ["ATG", "TGC", "GCG", "CGT"]
        for index in 1:length(observations)
            Test.@test poc_decoded_label_strings(oracle.batched, index) ==
                poc_decoded_label_strings(oracle.scalar, index)
            Test.@test oracle.batched.paths[index].score ===
                oracle.scalar.paths[index].score
        end
    end

    Test.@testset "single-read batch reproduces the scalar decoder" begin
        graph = Mycelia.Rhizomorph.build_ngram_graph(
            ["ABCDEF"], 3; dataset_id="batched_poc_single"
        )
        observations = [["ABC", "BCX", "CDE", "DEF"]]
        oracle = Mycelia.batched_equivalence_oracle(graph, observations)
        Test.@test oracle.passes
        Test.@test poc_decoded_label_strings(oracle.batched, 1) ==
            ["ABC", "BCD", "CDE", "DEF"]
    end

    Test.@testset "PoC guards reject out-of-scope configs" begin
        graph = Mycelia.Rhizomorph.build_ngram_graph(
            ["ABCDEF"], 3; dataset_id="batched_poc_guard"
        )
        observations = [["ABC", "BCD", "CDE", "DEF"]]

        # Beam-pruned (inexact) decode is out of PoC scope.
        beam_config = Mycelia.ViterbiCorrectionConfig(beam_width=2)
        Test.@test_throws ArgumentError Mycelia.batched_correct_observations(
            graph, observations; config=beam_config
        )
    end

    # Default `strand_mode=:auto` on a directed DNA graph with reverse-strand
    # evidence resolves to `:doublestrand`. That used to throw even though the
    # caller never asked for multi-strand decoding. Explicit `:doublestrand`
    # stays out of PoC scope.
    Test.@testset ":auto on a directed doublestrand graph does not throw" begin
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGTAC")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id="batched_poc_auto_ds", mode=:doublestrand
        )
        Test.@test Mycelia._resolve_viterbi_strand_mode(graph, :DNA, :auto) ==
                   :doublestrand

        observations = [
            [
                Kmers.DNAKmer{3}("ATG"),
                Kmers.DNAKmer{3}("TGC"),
                Kmers.DNAKmer{3}("GCG"),
                Kmers.DNAKmer{3}("CGT"),
            ],
        ]

        result = Mycelia.batched_correct_observations(graph, observations)
        Test.@test result.diagnostics[:strand_mode] == :doublestrand
        scalar = Mycelia.correct_observations(graph, observations)
        Test.@test result.diagnostics[:strand_mode] == scalar.diagnostics[:strand_mode]
        Test.@test poc_decoded_label_strings(result, 1) ==
                   poc_decoded_label_strings(scalar, 1)
        Test.@test result.paths[1].score === scalar.paths[1].score

        singlestrand_config = Mycelia.ViterbiCorrectionConfig(strand_mode=:singlestrand)
        forced = Mycelia.batched_correct_observations(
            graph, observations; config=singlestrand_config)
        Test.@test forced.diagnostics[:strand_mode] == :singlestrand

        explicit = Mycelia.ViterbiCorrectionConfig(strand_mode=:doublestrand)
        test_throws_message(
            ArgumentError,
            "explicit multi-strand"
        ) do
            Mycelia.batched_correct_observations(
                graph, observations; config=explicit)
        end
    end
end
