# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_variable_length_graph_correction_test.jl")'
# ```

import BioSequences
import FASTX
import MetaGraphsNext
import Mycelia
import Test

struct VariableLengthQualityObservation{S}
    sequence::S
    quality_scores::Vector{UInt8}
end

function variable_length_quality_record(
        identifier::AbstractString,
        sequence::AbstractString,
        phred::AbstractVector{<:Integer}
)::FASTX.FASTQ.Record
    quality = String([Char(score + 33) for score in phred])
    return FASTX.FASTQ.Record(identifier, sequence, quality)
end

function variable_length_decoded_label_strings(
        result::Mycelia.ViterbiCorrectionResult
)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

function variable_length_edge_overlaps(graph::MetaGraphsNext.MetaGraph)::Vector{Int}
    return sort([
        getproperty(graph[edge_label...], :overlap_length)
        for edge_label in MetaGraphsNext.edge_labels(graph)
    ])
end

Test.@testset "Variable-length graph Viterbi correction" begin
    Test.@testset "FASTA OLC graph correction recovers injected substitution" begin
        records = [
            FASTX.FASTA.Record("read_1", BioSequences.dna"ATGCG"),
            FASTX.FASTA.Record("read_2", BioSequences.dna"GCGTA"),
            FASTX.FASTA.Record("read_3", BioSequences.dna"GTACC")
        ]
        graph = Mycelia.Rhizomorph.build_fasta_graph_olc(records; min_overlap = 3)
        observed = [
            BioSequences.dna"ATGCG",
            BioSequences.dna"GCGTT",
            BioSequences.dna"GTACC"
        ]

        result = Mycelia.correct_observations(graph, [observed])
        path = only(result.paths)

        Test.@test variable_length_edge_overlaps(graph) == [3, 3]
        Test.@test result.diagnostics[:alphabet] == :DNA
        Test.@test result.diagnostics[:emission_model] == :alphabet_parameterized
        Test.@test result.diagnostics[:transition_model] == :normalized_overlap_length
        Test.@test path.diagnostics[:transition_scoring] == :normalized_overlap_length
        Test.@test variable_length_decoded_label_strings(result) == [
            "ATGCG", "GCGTA", "GTACC"
        ]
    end

    Test.@testset "FASTQ OLC graph correction reuses quality-aware emissions" begin
        records = [
            variable_length_quality_record("read_1", "ATGCG", [40, 40, 40, 40, 40]),
            variable_length_quality_record("read_2", "GCGTA", [40, 40, 40, 40, 40]),
            variable_length_quality_record("read_3", "GTACC", [40, 40, 40, 40, 40])
        ]
        graph = Mycelia.Rhizomorph.build_fastq_graph_olc(records; min_overlap = 3)
        observed = [
            BioSequences.dna"ATGCG",
            VariableLengthQualityObservation(
                BioSequences.dna"GCGTT", UInt8[2, 2, 2, 2, 2]
            ),
            BioSequences.dna"GTACC"
        ]

        result = Mycelia.correct_observations(graph, [observed])
        path = only(result.paths)

        Test.@test variable_length_edge_overlaps(graph) == [3, 3]
        Test.@test result.diagnostics[:alphabet] == :DNA
        Test.@test result.diagnostics[:emission_model] == :quality_aware
        Test.@test result.diagnostics[:transition_model] == :normalized_overlap_length
        Test.@test path.diagnostics[:emission_scoring] == :quality_aware
        Test.@test path.diagnostics[:transition_scoring] == :normalized_overlap_length
        Test.@test variable_length_decoded_label_strings(result) == [
            "ATGCG", "GCGTA", "GTACC"
        ]
    end

    Test.@testset "string OLC graph correction recovers injected text error" begin
        graph = Mycelia.Rhizomorph.build_string_graph(
            ["ABCDE", "CDEFG", "EFGHI"];
            dataset_id = "viterbi_b6_string",
            min_overlap = 3
        )
        observed = ["ABCDE", "CDXFG", "EFGHI"]

        result = Mycelia.correct_observations(graph, [observed])
        path = only(result.paths)

        Test.@test variable_length_edge_overlaps(graph) == [3, 3]
        Test.@test result.diagnostics[:alphabet] == :TEXT
        Test.@test result.diagnostics[:emission_model] == :alphabet_parameterized
        Test.@test result.diagnostics[:transition_model] == :normalized_overlap_length
        Test.@test path.diagnostics[:transition_scoring] == :normalized_overlap_length
        Test.@test variable_length_decoded_label_strings(result) == [
            "ABCDE", "CDEFG", "EFGHI"
        ]
    end
end
