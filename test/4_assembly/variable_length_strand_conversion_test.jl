# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/variable_length_strand_conversion_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/variable_length_strand_conversion_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Variable-length FASTA/FASTQ strand conversion tests

import Test
import Mycelia
import FASTX
import BioSequences
import MetaGraphsNext
import Graphs

function make_dna_variable_graph(vertex_data_type, edge_data_type)
    return MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = BioSequences.LongDNA{4},
        vertex_data_type = vertex_data_type,
        edge_data_type = edge_data_type,
        weight_function = Mycelia.Rhizomorph.compute_edge_weight
    )
end

function get_single_evidence_entry(data, dataset_id, observation_id)
    evidence = Mycelia.Rhizomorph.get_observation_evidence(data, dataset_id, observation_id)
    Test.@test !isnothing(evidence)
    Test.@test length(evidence) == 1
    return only(evidence)
end

function test_evidence_contains_strand(data, dataset_id, observation_id, strand)
    evidence = Mycelia.Rhizomorph.get_observation_evidence(data, dataset_id, observation_id)
    Test.@test !isnothing(evidence)
    Test.@test any(entry -> entry.strand == strand, evidence)
    return evidence
end

Test.@testset "Variable-length strand conversions" begin
    Test.@testset "FASTA and FASTQ builders convert basic singlestrand graphs" begin
        fasta_records = [FASTX.FASTA.Record("contig1", "ATGCAG")]
        ss_fasta = Mycelia.Rhizomorph.build_fasta_graph(
            fasta_records; dataset_id = "fasta", min_overlap = 3)
        ds_fasta = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(ss_fasta)
        Test.@test Mycelia.Rhizomorph.vertex_count(ds_fasta) >= 2
        canon_fasta = Mycelia.Rhizomorph.convert_variable_length_to_canonical(ss_fasta)
        Test.@test Mycelia.Rhizomorph.vertex_count(canon_fasta) == 1

        fastq_records = [FASTX.FASTQ.Record("read1", "ATGCAG", "IIIIII")]
        ss_fastq = Mycelia.Rhizomorph.build_fastq_graph(
            fastq_records; dataset_id = "fastq", min_overlap = 3)
        ds_fastq = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(ss_fastq)
        Test.@test Mycelia.Rhizomorph.vertex_count(ds_fastq) >= 2
        canon_fastq = Mycelia.Rhizomorph.convert_variable_length_to_canonical(ss_fastq)
        Test.@test Mycelia.Rhizomorph.vertex_count(canon_fastq) == 1
    end

    Test.@testset "Empty graphs return unchanged" begin
        empty_vertex = Mycelia.Rhizomorph.BioSequenceVertexData(BioSequences.LongDNA{4}("AT"))
        empty_graph = make_dna_variable_graph(
            typeof(empty_vertex), Mycelia.Rhizomorph.BioSequenceEdgeData)

        Test.@test Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(empty_graph) ===
                   empty_graph
        Test.@test Mycelia.Rhizomorph.convert_variable_length_to_canonical(empty_graph) ===
                   empty_graph
    end

    Test.@testset "Doublestrand conversion adds reverse-complement vertices and edges" begin
        seq_a = BioSequences.LongDNA{4}("ATGCC")
        seq_b = BioSequences.LongDNA{4}("CCAAT")
        rc_a = BioSequences.reverse_complement(seq_a)
        rc_b = BioSequences.reverse_complement(seq_b)

        graph = make_dna_variable_graph(
            typeof(Mycelia.Rhizomorph.BioSequenceVertexData(seq_a)),
            Mycelia.Rhizomorph.BioSequenceEdgeData
        )

        vdata_a = Mycelia.Rhizomorph.BioSequenceVertexData(seq_a)
        Mycelia.Rhizomorph.add_evidence!(
            vdata_a, "dataset", "seq_a", Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
        graph[seq_a] = vdata_a

        vdata_b = Mycelia.Rhizomorph.BioSequenceVertexData(seq_b)
        Mycelia.Rhizomorph.add_evidence!(
            vdata_b, "dataset", "seq_b", Mycelia.Rhizomorph.EvidenceEntry(2, Mycelia.Rhizomorph.Forward))
        graph[seq_b] = vdata_b

        edge_data = Mycelia.Rhizomorph.BioSequenceEdgeData(3)
        Mycelia.Rhizomorph.add_evidence!(
            edge_data,
            "dataset",
            "seq_a_to_seq_b",
            Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Forward)
        )
        graph[seq_a, seq_b] = edge_data

        ds_graph = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(graph)

        Test.@test MetaGraphsNext.nv(ds_graph) == 4
        Test.@test MetaGraphsNext.ne(ds_graph) == 2
        Test.@test haskey(ds_graph, rc_a)
        Test.@test haskey(ds_graph, rc_b)
        Test.@test haskey(ds_graph, rc_b, rc_a)
        Test.@test ds_graph[rc_b, rc_a].overlap_length == 3

        rc_edge_entry = get_single_evidence_entry(ds_graph[rc_b, rc_a], "dataset", "seq_a_to_seq_b")
        Test.@test rc_edge_entry.strand == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "Palindromic vertices do not duplicate during strand conversion" begin
        palindromic = BioSequences.LongDNA{4}("ATAT")
        graph = make_dna_variable_graph(
            typeof(Mycelia.Rhizomorph.BioSequenceVertexData(palindromic)),
            Mycelia.Rhizomorph.BioSequenceEdgeData
        )

        vertex_data = Mycelia.Rhizomorph.BioSequenceVertexData(palindromic)
        Mycelia.Rhizomorph.add_evidence!(
            vertex_data,
            "dataset",
            "palindrome",
            Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward)
        )
        graph[palindromic] = vertex_data

        ds_graph = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(graph)
        canon_graph = Mycelia.Rhizomorph.convert_variable_length_to_canonical(graph)

        Test.@test MetaGraphsNext.nv(ds_graph) == 1
        Test.@test MetaGraphsNext.ne(ds_graph) == 0
        Test.@test haskey(ds_graph, palindromic)
        Test.@test MetaGraphsNext.nv(canon_graph) == 1
        Test.@test haskey(canon_graph, palindromic)
    end

    Test.@testset "Canonical conversion merges reverse-complement evidence and edges" begin
        seq_a = BioSequences.LongDNA{4}("ATGCC")
        seq_b = BioSequences.LongDNA{4}("CCAAT")
        rc_a = BioSequences.reverse_complement(seq_a)
        rc_b = BioSequences.reverse_complement(seq_b)
        canon_a = BioSequences.canonical(seq_a)
        canon_b = BioSequences.canonical(seq_b)

        graph = make_dna_variable_graph(
            typeof(Mycelia.Rhizomorph.BioSequenceVertexData(seq_a)),
            Mycelia.Rhizomorph.BioSequenceEdgeData
        )

        vertex_a = Mycelia.Rhizomorph.BioSequenceVertexData(seq_a)
        Mycelia.Rhizomorph.add_evidence!(
            vertex_a,
            "dataset",
            "forward_vertex",
            Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward)
        )
        graph[seq_a] = vertex_a

        vertex_rc_a = Mycelia.Rhizomorph.BioSequenceVertexData(rc_a)
        Mycelia.Rhizomorph.add_evidence!(
            vertex_rc_a,
            "dataset",
            "reverse_vertex",
            Mycelia.Rhizomorph.EvidenceEntry(11, Mycelia.Rhizomorph.Forward)
        )
        graph[rc_a] = vertex_rc_a

        vertex_b = Mycelia.Rhizomorph.BioSequenceVertexData(seq_b)
        Mycelia.Rhizomorph.add_evidence!(
            vertex_b,
            "dataset",
            "forward_target",
            Mycelia.Rhizomorph.EvidenceEntry(20, Mycelia.Rhizomorph.Forward)
        )
        graph[seq_b] = vertex_b

        vertex_rc_b = Mycelia.Rhizomorph.BioSequenceVertexData(rc_b)
        Mycelia.Rhizomorph.add_evidence!(
            vertex_rc_b,
            "dataset",
            "reverse_target",
            Mycelia.Rhizomorph.EvidenceEntry(21, Mycelia.Rhizomorph.Forward)
        )
        graph[rc_b] = vertex_rc_b

        forward_edge = Mycelia.Rhizomorph.BioSequenceEdgeData(3)
        Mycelia.Rhizomorph.add_evidence!(
            forward_edge,
            "dataset",
            "forward_edge",
            Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Forward)
        )
        graph[seq_a, seq_b] = forward_edge

        mixed_edge = Mycelia.Rhizomorph.BioSequenceEdgeData(3)
        Mycelia.Rhizomorph.add_evidence!(
            mixed_edge,
            "dataset",
            "mixed_edge",
            Mycelia.Rhizomorph.EdgeEvidenceEntry(3, 4, Mycelia.Rhizomorph.Forward)
        )
        graph[rc_a, seq_b] = mixed_edge

        canon_graph = Mycelia.Rhizomorph.convert_variable_length_to_canonical(graph)

        Test.@test MetaGraphsNext.nv(canon_graph) == 2
        Test.@test MetaGraphsNext.ne(canon_graph) == 1
        Test.@test haskey(canon_graph, canon_a)
        Test.@test haskey(canon_graph, canon_b)
        Test.@test haskey(canon_graph, canon_a, canon_b)

        merged_vertex = canon_graph[canon_a]
        Test.@test Mycelia.Rhizomorph.count_total_observations(merged_vertex) == 2
        Test.@test get_single_evidence_entry(
            merged_vertex, "dataset", "forward_vertex").strand == Mycelia.Rhizomorph.Forward
        Test.@test get_single_evidence_entry(
            merged_vertex, "dataset", "reverse_vertex").strand == Mycelia.Rhizomorph.Reverse

        merged_edge = canon_graph[canon_a, canon_b]
        Test.@test Mycelia.Rhizomorph.count_total_observations(merged_edge) == 2
        test_evidence_contains_strand(
            merged_edge, "dataset", "forward_edge", Mycelia.Rhizomorph.Forward)
        test_evidence_contains_strand(
            merged_edge, "dataset", "mixed_edge", Mycelia.Rhizomorph.Reverse)
    end

    Test.@testset "Quality-aware doublestrand conversion infers quality edge data" begin
        seq = BioSequences.LongDNA{4}("ATGC")
        rc_seq = BioSequences.reverse_complement(seq)
        quality_scores = UInt8[31, 32, 33, 34]
        graph = make_dna_variable_graph(
            typeof(Mycelia.Rhizomorph.QualityBioSequenceVertexData(seq, quality_scores)),
            Mycelia.Rhizomorph.QualityBioSequenceEdgeData
        )

        vertex_data = Mycelia.Rhizomorph.QualityBioSequenceVertexData(seq, quality_scores)
        Mycelia.Rhizomorph.add_evidence!(
            vertex_data,
            "dataset",
            "quality_read",
            Mycelia.Rhizomorph.QualityEvidenceEntry(
                1, Mycelia.Rhizomorph.Forward, quality_scores)
        )
        graph[seq] = vertex_data

        ds_graph = Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(graph)

        Test.@test ds_graph[rc_seq] isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
        Test.@test ds_graph[rc_seq].quality_scores == quality_scores

        quality_edge = Mycelia.Rhizomorph.QualityBioSequenceEdgeData(2)
        Mycelia.Rhizomorph.add_evidence!(
            quality_edge,
            "dataset",
            "quality_edge",
            Mycelia.Rhizomorph.EdgeQualityEvidenceEntry(
                1,
                2,
                Mycelia.Rhizomorph.Forward,
                UInt8[31, 32],
                UInt8[33, 34]
            )
        )
        ds_graph[seq, rc_seq] = quality_edge

        stored_edge = ds_graph[seq, rc_seq]
        Test.@test stored_edge isa Mycelia.Rhizomorph.QualityBioSequenceEdgeData
        Test.@test stored_edge.overlap_length == 2
        Test.@test get_single_evidence_entry(
            stored_edge, "dataset", "quality_edge").strand == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "Quality-aware palindromic vertices preserve quality data" begin
        palindromic = BioSequences.LongDNA{4}("ATAT")
        quality_scores = UInt8[35, 36, 37, 38]
        graph = make_dna_variable_graph(
            typeof(Mycelia.Rhizomorph.QualityBioSequenceVertexData(palindromic, quality_scores)),
            Mycelia.Rhizomorph.QualityBioSequenceEdgeData
        )

        vertex_data = Mycelia.Rhizomorph.QualityBioSequenceVertexData(palindromic, quality_scores)
        Mycelia.Rhizomorph.add_evidence!(
            vertex_data,
            "dataset",
            "palindrome_fastq",
            Mycelia.Rhizomorph.QualityEvidenceEntry(
                1, Mycelia.Rhizomorph.Forward, quality_scores)
        )
        graph[palindromic] = vertex_data

        canon_graph = Mycelia.Rhizomorph.convert_variable_length_to_canonical(graph)
        converted_vertex = canon_graph[palindromic]

        Test.@test converted_vertex isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
        Test.@test converted_vertex.quality_scores == quality_scores
        Test.@test get_single_evidence_entry(
            converted_vertex, "dataset", "palindrome_fastq").quality_scores == quality_scores
    end

    Test.@testset "Non-nucleotide graphs reject strand conversions" begin
        string_records = ["hello world"]
        string_graph = Mycelia.Rhizomorph.build_string_graph(
            string_records; dataset_id = "str", min_overlap = 3)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(
            string_graph)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.convert_variable_length_to_canonical(
            string_graph)
    end
end

println("✓ Variable-length strand conversion tests completed")
