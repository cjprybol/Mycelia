# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_fixed_length_graph_correction_test.jl")'
# ```

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

struct FixedLengthQualityKmerObservation{K}
    Kmer::K
    quality_scores::Vector{UInt8}
end

function fixed_length_quality_record(
    identifier::AbstractString, sequence::AbstractString, phred::AbstractVector{<:Integer}
)::FASTX.FASTQ.Record
    quality = String([Char(score + 33) for score in phred])
    return FASTX.FASTQ.Record(identifier, sequence, quality)
end

function fixed_length_decoded_label_strings(
    result::Mycelia.ViterbiCorrectionResult
)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

Test.@testset "Fixed-length graph Viterbi correction" begin
    Test.@testset "k-mer graph correction recovers injected substitution" begin
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id="viterbi_b5_kmer", mode=:singlestrand
        )
        observed = [
            Kmers.DNAKmer{3}("ATG"),
            Kmers.DNAKmer{3}("TGA"),
            Kmers.DNAKmer{3}("GCG"),
            Kmers.DNAKmer{3}("CGT"),
        ]

        result = Mycelia.correct_observations(graph, [observed])
        beam_result = Mycelia.correct_observations(
            graph,
            [observed];
            config = Mycelia.ViterbiCorrectionConfig(beam_width = 1)
        )

        Test.@test result.diagnostics[:alphabet] == :DNA
        Test.@test result.diagnostics[:emission_model] == :alphabet_parameterized
        Test.@test fixed_length_decoded_label_strings(result) ==
            ["ATG", "TGC", "GCG", "CGT"]
        Test.@test beam_result.diagnostics[:beam_width] == 1
        Test.@test only(beam_result.paths).diagnostics[:max_retained_states] <= 1
        Test.@test fixed_length_decoded_label_strings(beam_result) ==
            fixed_length_decoded_label_strings(result)
    end

    Test.@testset "beam width validation rejects non-positive beams" begin
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_width = 0)
        Test.@test_throws ArgumentError Mycelia.ViterbiCorrectionConfig(beam_width = -1)
    end

    Test.@testset "qualmer graph correction reuses quality-aware emissions" begin
        records = FASTX.FASTQ.Record[]
        for index in 1:12
            push!(
                records,
                fixed_length_quality_record(
                    "truth_$index", "ATGCGT", [40, 40, 40, 40, 40, 40]
                ),
            )
        end
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; dataset_id="viterbi_b5_qualmer", mode=:singlestrand
        )
        observed = [
            FixedLengthQualityKmerObservation(Kmers.DNAKmer{3}("ATG"), UInt8[40, 40, 40]),
            FixedLengthQualityKmerObservation(Kmers.DNAKmer{3}("TGA"), UInt8[2, 2, 2]),
            FixedLengthQualityKmerObservation(Kmers.DNAKmer{3}("GCG"), UInt8[40, 40, 40]),
            FixedLengthQualityKmerObservation(Kmers.DNAKmer{3}("CGT"), UInt8[40, 40, 40]),
        ]

        result = Mycelia.correct_observations(graph, [observed])
        path = only(result.paths)

        Test.@test result.diagnostics[:alphabet] == :DNA
        Test.@test result.diagnostics[:emission_model] == :quality_aware
        Test.@test path.diagnostics[:emission_scoring] == :quality_aware
        Test.@test fixed_length_decoded_label_strings(result) ==
            ["ATG", "TGC", "GCG", "CGT"]
    end

    Test.@testset "n-gram graph correction treats text labels as fixed-length text" begin
        graph = Mycelia.Rhizomorph.build_ngram_graph(
            ["ABCDEF"], 3; dataset_id="viterbi_b5_ngram"
        )
        observed = ["ABC", "BCX", "CDE", "DEF"]

        result = Mycelia.correct_observations(graph, [observed])
        path = only(result.paths)

        Test.@test result.diagnostics[:alphabet] == :TEXT
        Test.@test result.diagnostics[:emission_model] == :alphabet_parameterized
        Test.@test path.diagnostics[:emission_scoring] == :alphabet_parameterized
        Test.@test fixed_length_decoded_label_strings(result) ==
            ["ABC", "BCD", "CDE", "DEF"]
    end
end
