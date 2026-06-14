import Test
import BioSequences
import FASTX
import Graphs
import Kmers
import MetaGraphsNext
import Mycelia

function assembly_core_test_error_contains(f, expected::AbstractString)
    did_throw = false
    try
        f()
    catch err
        did_throw = true
        Test.@test occursin(expected, sprint(showerror, err))
    end
    Test.@test did_throw
end

function assembly_core_test_write_fasta(path::AbstractString, records)
    open(path, "w") do io
        for (identifier, sequence) in records
            println(io, ">$identifier")
            println(io, sequence)
        end
    end
    return path
end

function assembly_core_test_write_fastq(path::AbstractString, records)
    open(path, "w") do io
        for (identifier, sequence, quality) in records
            println(io, "@$identifier")
            println(io, sequence)
            println(io, "+")
            println(io, quality)
        end
    end
    return path
end

function assembly_core_test_empty_kmer_graph()
    return MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = Kmers.DNAKmer{3},
        vertex_data_type = Mycelia.Rhizomorph.KmerVertexData{Kmers.DNAKmer{3}},
        edge_data_type = Mycelia.Rhizomorph.KmerEdgeData,
        weight_function = Mycelia.Rhizomorph.compute_edge_weight
    )
end

function assembly_core_test_empty_qualmer_graph()
    return MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = Kmers.DNAKmer{3},
        vertex_data_type = Mycelia.Rhizomorph.QualmerVertexData{Kmers.DNAKmer{3}},
        edge_data_type = Mycelia.Rhizomorph.QualmerEdgeData,
        weight_function = Mycelia.Rhizomorph.compute_edge_weight
    )
end

function assembly_core_test_empty_fastq_graph()
    return MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = BioSequences.LongDNA{4},
        vertex_data_type = Mycelia.Rhizomorph.QualityBioSequenceVertexData{
            BioSequences.LongDNA{4}
        },
        edge_data_type = Mycelia.Rhizomorph.QualityBioSequenceEdgeData,
        weight_function = Mycelia.Rhizomorph.compute_edge_weight
    )
end

Test.@testset "Assembly Core Helper Branch Coverage" begin
    Test.@testset "K-mer graph statistics and filters" begin
        records = [
            FASTX.FASTA.Record("read1", "ATGCAT"),
            FASTX.FASTA.Record("read2", "ATGCAT")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; dataset_id = "reads")
        stats = Mycelia.Rhizomorph.get_kmer_statistics(graph)

        Test.@test stats[:k] == 3
        Test.@test stats[:sequence_type] == "DNA"
        Test.@test stats[:num_datasets] == 1
        Test.@test stats[:num_vertices] == MetaGraphsNext.nv(graph)

        high_global = Mycelia.Rhizomorph.find_high_coverage_kmers(graph, 2)
        high_dataset = Mycelia.Rhizomorph.find_high_coverage_kmers(
            graph, 2; dataset_id = "reads")
        low_dataset = Mycelia.Rhizomorph.find_low_coverage_kmers(
            graph, 1; dataset_id = "reads")

        Test.@test Set(high_global) == Set(MetaGraphsNext.labels(graph))
        Test.@test Set(high_dataset) == Set(MetaGraphsNext.labels(graph))
        Test.@test isempty(low_dataset)

        empty_stats = Mycelia.Rhizomorph.get_kmer_statistics(
            assembly_core_test_empty_kmer_graph())
        Test.@test empty_stats[:num_vertices] == 0
        Test.@test empty_stats[:k] === nothing

        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph.build_kmer_graph(
                FASTX.FASTQ.Record[], 3; memory_profile = :ultralight_quality),
            "Cannot build graph from empty record set"
        )
        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph.build_kmer_graph(
                records, 3; memory_profile = :compact),
            "Invalid memory_profile"
        )
        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph.build_kmer_graph(records, 3; mode = :triplestrand),
            "Invalid mode"
        )
    end

    Test.@testset "K-mer graph file helper validation" begin
        mktempdir() do dir
            fasta_a = assembly_core_test_write_fasta(
                joinpath(dir, "sample_a.fna"), [("seq1", "ATGCAT")])
            fasta_b = assembly_core_test_write_fasta(
                joinpath(dir, "sample_b.fna"), [("seq2", "ATGCAT")])

            graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [fasta_a, fasta_b],
                3;
                memory_profile = :lightweight,
                mode = :singlestrand
            )
            dataset_ids = Set{String}()
            for label in MetaGraphsNext.labels(graph)
                union!(dataset_ids, keys(graph[label].dataset_counts))
            end
            Test.@test dataset_ids == Set(["sample_a", "sample_b"])

            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_kmer_graph_from_file(
                    joinpath(dir, "missing.fasta"), 3),
                "File not found"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_kmer_graph_from_files(String[], 3),
                "No files provided"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_kmer_graph_from_files(
                    [fasta_a], 3; mode = :triplestrand),
                "Invalid mode"
            )
        end
    end

    Test.@testset "Qualmer graph statistics and quality filters" begin
        records = [
            FASTX.FASTQ.Record("read1", "ATGCAT", "IIIIII"),
            FASTX.FASTQ.Record("read2", "ATGCAT", "IIIIII")
        ]
        graph = Mycelia.Rhizomorph.build_qualmer_graph(records, 3; dataset_id = "reads")
        stats = Mycelia.Rhizomorph.get_qualmer_statistics(graph)

        Test.@test stats[:k] == 3
        Test.@test stats[:sequence_type] == "DNA"
        Test.@test stats[:num_datasets] == 1
        Test.@test stats[:mean_joint_quality] !== nothing
        Test.@test stats[:min_joint_quality] <= stats[:max_joint_quality]

        high_quality = Mycelia.Rhizomorph.find_high_quality_kmers(
            graph, 80; dataset_id = "reads")
        missing_dataset = Mycelia.Rhizomorph.find_high_quality_kmers(
            graph, 1; dataset_id = "missing")
        too_high = Mycelia.Rhizomorph.find_high_quality_kmers(graph, 81)

        Test.@test Set(high_quality) == Set(MetaGraphsNext.labels(graph))
        Test.@test isempty(missing_dataset)
        Test.@test isempty(too_high)

        empty_stats = Mycelia.Rhizomorph.get_qualmer_statistics(
            assembly_core_test_empty_qualmer_graph())
        Test.@test empty_stats[:num_vertices] == 0
        Test.@test empty_stats[:mean_joint_quality] === nothing

        ultralight = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; memory_profile = :ultralight_quality)
        lightweight = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; memory_profile = :lightweight_quality)
        Test.@test MetaGraphsNext.nv(ultralight) > 0
        Test.@test MetaGraphsNext.nv(lightweight) > 0

        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph.build_qualmer_graph(
                records, 3; memory_profile = :compact),
            "Invalid memory_profile"
        )
        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph.build_qualmer_graph(records, 3; mode = :triplestrand),
            "Invalid mode"
        )
    end

    Test.@testset "Qualmer graph file helper validation" begin
        mktempdir() do dir
            fastq_a = assembly_core_test_write_fastq(
                joinpath(dir, "sample_a.fastq"), [("read1", "ATGCAT", "IIIIII")])
            fastq_b = assembly_core_test_write_fastq(
                joinpath(dir, "sample_b.fastq"), [("read2", "ATGCAT", "IIIIII")])
            fasta = assembly_core_test_write_fasta(
                joinpath(dir, "not_fastq.fasta"), [("seq1", "ATGCAT")])
            empty = joinpath(dir, "empty.fastq")
            touch(empty)

            graph = Mycelia.Rhizomorph.build_qualmer_graph_from_files(
                [fastq_a, fastq_b], 3)
            Test.@test MetaGraphsNext.nv(graph) > 0
            Test.@test Set(keys(graph[first(MetaGraphsNext.labels(graph))].evidence)) ==
                       Set(["sample_a", "sample_b"])

            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_qualmer_graph_from_file(
                    joinpath(dir, "missing.fastq"), 3),
                "File not found"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_qualmer_graph_from_file(empty, 3),
                "No records found"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_qualmer_graph_from_file(fasta, 3),
                "Qualmer graphs require FASTQ input"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_qualmer_graph_from_files(
                    [fastq_a, fasta], 3),
                "Qualmer graphs require FASTQ input"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_qualmer_graph_from_files(String[], 3),
                "No files provided"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_qualmer_graph_from_files(
                    [fastq_a], 3; mode = :triplestrand),
                "Invalid mode"
            )
        end
    end

    Test.@testset "FASTQ graph statistics and quality filters" begin
        records = [
            FASTX.FASTQ.Record("high", "ATGCA", "IIIII"),
            FASTX.FASTQ.Record("low", "GCATG", "!!!!!")
        ]
        graph = Mycelia.Rhizomorph.build_fastq_graph(
            records; dataset_id = "reads", min_overlap = 3)
        stats = Mycelia.Rhizomorph.get_fastq_graph_statistics(
            graph; dataset_id = "reads")

        Test.@test stats[:sequence_type] == "DNA"
        Test.@test stats[:num_vertices] == 2
        Test.@test stats[:num_edges] >= 1
        Test.@test stats[:mean_overlap_length] !== nothing
        Test.@test stats[:mean_sequence_length] == 5.0
        Test.@test stats[:mean_quality] !== nothing

        high_quality = Mycelia.Rhizomorph.find_high_quality_sequences(
            graph, 40; dataset_id = "reads")
        missing_dataset = Mycelia.Rhizomorph.find_high_quality_sequences(
            graph, 1; dataset_id = "missing")
        too_high = Mycelia.Rhizomorph.find_high_quality_sequences(graph, 41)

        Test.@test length(high_quality) == 1
        Test.@test string(only(high_quality)) == "ATGCA"
        Test.@test isempty(missing_dataset)
        Test.@test isempty(too_high)

        empty_stats = Mycelia.Rhizomorph.get_fastq_graph_statistics(
            assembly_core_test_empty_fastq_graph())
        Test.@test empty_stats[:num_vertices] == 0
        Test.@test empty_stats[:mean_quality] === nothing

        qualityless_graph = assembly_core_test_empty_fastq_graph()
        qualityless_sequence = BioSequences.LongDNA{4}("ATGC")
        qualityless_graph[qualityless_sequence] =
            Mycelia.Rhizomorph.QualityBioSequenceVertexData(qualityless_sequence)
        qualityless_records = Mycelia.Rhizomorph.fastq_graph_to_records(
            qualityless_graph, "contig")
        Test.@test length(qualityless_records) == 1
        Test.@test FASTX.FASTQ.sequence(qualityless_records[1]) == "ATGC"
        Test.@test FASTX.FASTQ.quality(qualityless_records[1]) == "IIII"
    end

    Test.@testset "FASTQ graph file helper validation" begin
        mktempdir() do dir
            fastq_a = assembly_core_test_write_fastq(
                joinpath(dir, "sample_a.fastq"), [("read1", "ATGCA", "IIIII")])
            fastq_b = assembly_core_test_write_fastq(
                joinpath(dir, "sample_b.fastq"), [("read2", "GCATG", "IIIII")])
            fasta = assembly_core_test_write_fasta(
                joinpath(dir, "not_fastq.fasta"), [("seq1", "ATGCA")])
            empty = joinpath(dir, "empty.fastq")
            touch(empty)

            graph = Mycelia.Rhizomorph.build_fastq_graph_from_files(
                [fastq_a, fastq_b]; min_overlap = 3)
            Test.@test MetaGraphsNext.nv(graph) > 0

            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_fastq_graph_from_file(
                    joinpath(dir, "missing.fastq")),
                "File not found"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_fastq_graph_from_file(empty),
                "No records found"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_fastq_graph_from_file(fasta),
                "FASTQ graphs require FASTQ input"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_fastq_graph_from_files([fastq_a, fasta]),
                "FASTQ graphs require FASTQ input"
            )
            assembly_core_test_error_contains(
                () -> Mycelia.Rhizomorph.build_fastq_graph_from_files(String[]),
                "No files provided"
            )
        end
    end

    Test.@testset "Variable-length strand conversion validation" begin
        string_graph = Mycelia.Rhizomorph.build_string_graph(["alpha", "lphab"]; min_overlap = 3)

        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph.convert_variable_length_to_doublestrand(string_graph),
            "Doublestrand conversion only supported for DNA/RNA"
        )
        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph.convert_variable_length_to_canonical(string_graph),
            "Canonical conversion only supported for DNA/RNA"
        )
    end

    Test.@testset "Rhizomorph assembly result helpers" begin
        Test.@test Mycelia.Rhizomorph._graph_mode_symbol(
            Mycelia.Rhizomorph.SingleStrand) == :singlestrand
        Test.@test Mycelia.Rhizomorph._graph_mode_symbol(
            Mycelia.Rhizomorph.DoubleStrand) == :doublestrand

        qualmer_graph = Mycelia.Rhizomorph.build_qualmer_graph(
            [FASTX.FASTQ.Record("read1", "ATGCAT", "IIIIII")],
            3;
            dataset_id = "reads"
        )
        label = first(MetaGraphsNext.labels(qualmer_graph))
        Test.@test Mycelia.Rhizomorph._qualmer_vertex_score(qualmer_graph[label]) > 0.0
        Test.@test length(Mycelia.Rhizomorph._qualmer_vertex_quality_scores(
            qualmer_graph[label])) == 3

        edge_label = first(MetaGraphsNext.edge_labels(qualmer_graph))
        Test.@test Mycelia.Rhizomorph._qualmer_edge_score(
            qualmer_graph[edge_label...]) >= 1.0
        Test.@test Mycelia.Rhizomorph._qualmer_edge_score(
            Mycelia.Rhizomorph.QualmerEdgeData()) == 1.0

        inconsistent = Mycelia.Rhizomorph.AssemblyResult(
            ["ATGC"], ["contig_1", "contig_2"])
        report = Mycelia.Rhizomorph.validate_assembly_structure(inconsistent)
        Test.@test !report["valid"]
        Test.@test "Contigs and contig_names have different lengths" in report["issues"]
        Test.@test "Marked as GFA compatible but no graph structure" in report["issues"]
        Test.@test !Mycelia.Rhizomorph.has_graph_structure(inconsistent)
        Test.@test !Mycelia.Rhizomorph.has_simplified_graph(inconsistent)
        Test.@test !Mycelia.Rhizomorph.has_paths(inconsistent)

        fastq_result = Mycelia.Rhizomorph.AssemblyResult(
            ["ATGC"],
            ["contig_1"];
            assembly_stats = Dict("quality_preserved" => true),
            fastq_contigs = [FASTX.FASTQ.Record("contig_1", "ATGC", "IIII")]
        )
        Test.@test Mycelia.Rhizomorph.has_quality_information(fastq_result)
        Test.@test length(Mycelia.Rhizomorph.get_fastq_contigs(fastq_result)) == 1

        mktempdir() do dir
            fastq_path = joinpath(dir, "contigs.fastq")
            saved_path = Mycelia.Rhizomorph.write_fastq_contigs(fastq_result, fastq_path)
            Test.@test saved_path == fastq_path
            Test.@test occursin("@contig_1", read(fastq_path, String))

            jld2_path = joinpath(dir, "assembly.jld2")
            Mycelia.Rhizomorph.save_assembly(fastq_result, jld2_path)
            loaded = Mycelia.Rhizomorph.load_assembly(jld2_path)
            Test.@test loaded.contigs == fastq_result.contigs
            Test.@test loaded.contig_names == fastq_result.contig_names
        end

        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph.write_gfa(inconsistent, "unused.gfa"),
            "AssemblyResult contains no graph"
        )
        assembly_core_test_error_contains(
            () -> Mycelia.Rhizomorph._prepare_observations([1, 2, 3]),
            "Unsupported reads format"
        )
    end

    Test.@testset "Iterative assembly decision helpers" begin
        Test.@test Mycelia._vertex_joint_probability((joint_probability = 0.25,)) == 0.25
        Test.@test Mycelia._vertex_mean_quality((mean_quality = 37.5,)) == 37.5
        Test.@test Mycelia._vertex_coverage((coverage = 4,)) == 4
        Test.@test Mycelia._vertex_coverage((coverage = [1, 2, 3],)) == 3

        empty_vertex = Mycelia.Rhizomorph.QualmerVertexData(Kmers.DNAKmer{3}("ATG"))
        Test.@test Mycelia._rhizomorph_first_dataset_id(empty_vertex) === nothing
        Test.@test Mycelia._rhizomorph_joint_probability(empty_vertex) == 0.0
        Test.@test Mycelia._vertex_mean_quality(empty_vertex) == 0.0
        Test.@test Mycelia._vertex_coverage(empty_vertex) == 0

        Test.@test Mycelia.calculate_unseen_kmer_penalty(Int8[40, 40, 40]) >
                   Mycelia.calculate_unseen_kmer_penalty(Int8[0, 0, 0])

        decreasing_history = [
            Dict{Symbol, Any}(:improvement_rate => 0.09),
            Dict{Symbol, Any}(:improvement_rate => 0.08),
            Dict{Symbol, Any}(:improvement_rate => 0.07)
        ]
        plateau_history = [
            Dict{Symbol, Any}(:improvement_rate => 0.01),
            Dict{Symbol, Any}(:improvement_rate => 0.02),
            Dict{Symbol, Any}(:improvement_rate => 0.01)
        ]
        Test.@test !Mycelia.sufficient_improvements(
            10, 100, 0.05; iteration_history = decreasing_history)
        Test.@test !Mycelia.sufficient_improvements(
            10, 100, 0.05; iteration_history = plateau_history)

        short_history = [Dict{Symbol, Any}(
            :total_improvements => 1,
            :total_reads => 10,
            :runtime_seconds => 0.1
        )]
        slow_converged_history = [
            Dict{Symbol, Any}(
                :total_improvements => 0,
                :total_reads => 10,
                :runtime_seconds => 2.0
            )
            for _ in 1:3
        ]
        productive_history = [
            Dict{Symbol, Any}(
                :total_improvements => 2,
                :total_reads => 10,
                :runtime_seconds => 0.1
            )
            for _ in 1:3
        ]

        Test.@test !Mycelia.should_continue_k_progression(productive_history, 11, 11)
        Test.@test Mycelia.should_continue_k_progression(short_history, 5, 11)
        Test.@test !Mycelia.should_continue_k_progression(slow_converged_history, 5, 11)
        Test.@test Mycelia.should_continue_k_progression(productive_history, 5, 11)
    end
end
