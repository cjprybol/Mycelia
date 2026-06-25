# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_quality_emission_test.jl")'
# ```

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

struct QualityKmerObservation{K}
    Kmer::K
    quality_scores::Vector{UInt8}
end

function flat_viterbi_emission_logp(
        observed_unit::Any,
        node::Any,
        alphabet::Symbol;
        quality = nothing,
        error_rate::Float64 = 0.01
)::Float64
    return Mycelia.default_viterbi_emission_logp(
        observed_unit,
        node,
        alphabet;
        error_rate = error_rate
    )
end

function quality_record(
        identifier::AbstractString,
        sequence::AbstractString,
        phred::AbstractVector{<:Integer}
)::FASTX.FASTQ.Record
    quality = String([Char(score + 33) for score in phred])
    return FASTX.FASTQ.Record(identifier, sequence, quality)
end

function decoded_label_strings(result::Mycelia.ViterbiCorrectionResult)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

Test.@testset "Quality-aware Viterbi emissions" begin
    Test.@testset "per-base Phred changes mismatch likelihood relative to flat emissions" begin
        observed = Kmers.DNAKmer{3}("TGA")
        expected = Kmers.DNAKmer{3}("TGC")
        flat = Mycelia.default_viterbi_emission_logp(
            observed,
            expected,
            :DNA;
            error_rate = 0.01
        )
        low_quality = Mycelia.default_viterbi_emission_logp(
            observed,
            expected,
            :DNA;
            quality = UInt8[2, 2, 2],
            error_rate = 0.01
        )
        high_quality = Mycelia.default_viterbi_emission_logp(
            observed,
            expected,
            :DNA;
            quality = UInt8[40, 40, 40],
            error_rate = 0.01
        )

        Test.@test low_quality > flat
        Test.@test high_quality < flat
    end

    Test.@testset "qualmer evidence combines Phred across coverage" begin
        records = [
            quality_record("read_1", "ATGC", [20, 20, 20, 20]),
            quality_record("read_2", "ATGC", [30, 30, 30, 30])
        ]
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records,
            3;
            dataset_id = "viterbi_quality_combined",
            mode = :singlestrand
        )
        observation = [Kmers.DNAKmer{3}("ATG")]

        result = Mycelia.correct_observations(graph, [observation])
        path = only(result.paths)
        expected_score = Mycelia.default_viterbi_emission_logp(
            Kmers.DNAKmer{3}("ATG"),
            Kmers.DNAKmer{3}("ATG"),
            :DNA;
            quality = Float64[50.0, 50.0, 50.0]
        )

        Test.@test result.diagnostics[:emission_model] == :quality_aware
        Test.@test path.diagnostics[:emission_scoring] == :quality_aware
        Test.@test path.score ≈ expected_score
    end

    Test.@testset "quality-aware correction can differ from flat-rate correction" begin
        records = FASTX.FASTQ.Record[]
        for index in 1:20
            push!(records, quality_record("truth_$index", "ATGC", [40, 40, 40, 40]))
        end
        push!(records, quality_record("error", "ATGA", [40, 2, 2, 2]))
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records,
            3;
            dataset_id = "viterbi_quality_choice",
            mode = :singlestrand
        )
        observed = [
            QualityKmerObservation(Kmers.DNAKmer{3}("ATG"), UInt8[40, 40, 40]),
            QualityKmerObservation(Kmers.DNAKmer{3}("TGA"), UInt8[2, 2, 2])
        ]

        flat_config = Mycelia.ViterbiCorrectionConfig(
            emission_logp = flat_viterbi_emission_logp
        )
        flat_result = Mycelia.correct_observations(graph, [observed]; config = flat_config)
        quality_result = Mycelia.correct_observations(graph, [observed])

        Test.@test decoded_label_strings(flat_result) == ["ATG", "TGA"]
        Test.@test decoded_label_strings(quality_result) == ["ATG", "TGC"]
        Test.@test only(quality_result.paths).diagnostics[:emission_scoring] == :quality_aware
    end
end
