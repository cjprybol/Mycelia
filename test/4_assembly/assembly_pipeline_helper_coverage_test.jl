import Test
import Mycelia
import FASTX
import BioSequences
import MetaGraphsNext
import Graphs

struct MockAssemblyQualityVertexData
    evidence::Int
    dataset_ids::Vector{String}
    mean_quality::Union{Nothing, Vector{Float64}}
end

struct MockAssemblyEdgeData
    evidence::Int
end

Mycelia.Rhizomorph.count_evidence(data::MockAssemblyQualityVertexData) = data.evidence
Mycelia.Rhizomorph.count_evidence(data::MockAssemblyEdgeData) = data.evidence
Mycelia.Rhizomorph.get_all_dataset_ids(data::MockAssemblyQualityVertexData) = data.dataset_ids
Mycelia.Rhizomorph.get_vertex_mean_quality(
    data::MockAssemblyQualityVertexData,
    ::String
) = data.mean_quality

function build_mock_quality_graph()
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph(),
        label_type = String,
        vertex_data_type = MockAssemblyQualityVertexData,
        edge_data_type = MockAssemblyEdgeData,
        default_weight = 1.0
    )

    graph["AAA"] = MockAssemblyQualityVertexData(4, ["mock"], [35.2, 36.8, 42.1])
    graph["AAT"] = MockAssemblyQualityVertexData(3, ["mock"], [33.0, 38.0, 41.0])
    graph["ATC"] = MockAssemblyQualityVertexData(2, ["mock"], [31.0, 37.0, 39.0])
    graph["AAA", "AAT"] = MockAssemblyEdgeData(5)
    graph["AAT", "ATC"] = MockAssemblyEdgeData(3)

    return graph
end

function build_zero_quality_graph()
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph(),
        label_type = String,
        vertex_data_type = MockAssemblyQualityVertexData,
        edge_data_type = MockAssemblyEdgeData,
        default_weight = 1.0
    )

    graph["AAA"] = MockAssemblyQualityVertexData(0, ["mock"], nothing)
    graph["AAT"] = MockAssemblyQualityVertexData(0, ["mock"], nothing)
    graph["AAA", "AAT"] = MockAssemblyEdgeData(1)

    return graph
end

Test.@testset "Assembly Pipeline Helper Coverage" begin
    Test.@testset "Qualmer Helper Utilities" begin
        Test.@test Mycelia.Rhizomorph._graph_mode_symbol(Mycelia.Rhizomorph.SingleStrand) ==
                   :singlestrand
        Test.@test Mycelia.Rhizomorph._graph_mode_symbol(Mycelia.Rhizomorph.DoubleStrand) ==
                   :doublestrand

        mock_vertex = MockAssemblyQualityVertexData(7, ["mock"], [12.4, 61.1, -2.0])
        empty_dataset_vertex = MockAssemblyQualityVertexData(2, String[], [20.0, 21.0, 22.0])
        missing_quality_vertex = MockAssemblyQualityVertexData(2, ["mock"], nothing)
        mock_edge = MockAssemblyEdgeData(0)

        Test.@test Mycelia.Rhizomorph._qualmer_vertex_score(mock_vertex) == 7.0
        Test.@test Mycelia.Rhizomorph._qualmer_edge_score(mock_edge) == 1.0
        Test.@test Mycelia.Rhizomorph._qualmer_vertex_quality_scores(mock_vertex) ==
                   UInt8[12, 60, 0]
        Test.@test isempty(
            Mycelia.Rhizomorph._qualmer_vertex_quality_scores(empty_dataset_vertex)
        )
        Test.@test isempty(
            Mycelia.Rhizomorph._qualmer_vertex_quality_scores(missing_quality_vertex)
        )

        graph = build_mock_quality_graph()
        path = ["AAA", "AAT", "ATC"]

        Test.@test Mycelia.Rhizomorph._qualmer_path_to_sequence(path, graph) == "AAATC"
        Test.@test Mycelia.Rhizomorph._qualmer_path_to_sequence(String[], graph) == ""

        consensus = Mycelia.Rhizomorph._qualmer_path_to_consensus_fastq(
            path,
            graph,
            "mock_contig"
        )
        Test.@test String(FASTX.identifier(consensus)) == "mock_contig"
        Test.@test String(FASTX.sequence(consensus)) == "AAATC"
        Test.@test collect(UInt8, FASTX.quality_scores(consensus)) == UInt8[35, 37, 42, 41, 39]

        empty_consensus = Mycelia.Rhizomorph._qualmer_path_to_consensus_fastq(
            String[],
            graph,
            "empty"
        )
        Test.@test isempty(FASTX.sequence(empty_consensus))
        Test.@test isempty(collect(UInt8, FASTX.quality_scores(empty_consensus)))

        heaviest_path = Mycelia.Rhizomorph._find_heaviest_eulerian_path(graph, "AAA", Set())
        Test.@test heaviest_path == path

        viterbi_path = Mycelia.Rhizomorph._viterbi_optimal_path(graph, "AAA", Set(), 10)
        Test.@test viterbi_path == path

        weighted_walk = Mycelia.Rhizomorph._quality_weighted_walk(graph, "AAA", 10)
        Test.@test weighted_walk == path

        config = Mycelia.Rhizomorph.AssemblyConfig(
            k = 2,
            sequence_type = BioSequences.LongDNA{4},
            graph_mode = Mycelia.Rhizomorph.SingleStrand,
            use_quality_scores = true
        )
        discovered_paths = Mycelia.Rhizomorph._find_qualmer_paths(graph, config)
        Test.@test path in discovered_paths

        empty_graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph(),
            label_type = String,
            vertex_data_type = MockAssemblyQualityVertexData,
            edge_data_type = MockAssemblyEdgeData,
            default_weight = 1.0
        )
        Test.@test isempty(Mycelia.Rhizomorph._find_qualmer_paths(empty_graph, config))
        Test.@test isempty(
            Mycelia.Rhizomorph._generate_fastq_contigs_from_qualmer_graph(empty_graph, config)
        )

        fastq_contigs = Mycelia.Rhizomorph._generate_fastq_contigs_from_qualmer_graph(
            graph,
            config
        )
        Test.@test length(fastq_contigs) == 1
        Test.@test String(FASTX.sequence(only(fastq_contigs))) == "AAATC"

        zero_quality_contigs = Mycelia.Rhizomorph._generate_fastq_contigs_from_qualmer_graph(
            build_zero_quality_graph(),
            config
        )
        Test.@test length(zero_quality_contigs) == 1
        Test.@test collect(UInt8, FASTX.quality_scores(only(zero_quality_contigs))) ==
                   fill(UInt8(2), 4)

        contigs = Mycelia.Rhizomorph._generate_contigs_from_qualmer_graph(graph, config)
        Test.@test "AAATC" in contigs
    end

    Test.@testset "Observation and FASTQ Preparation Helpers" begin
        fastq_record = FASTX.FASTQ.Record("read_fastq", "ATCG", "IIII")
        fasta_record = FASTX.FASTA.Record("read_fasta", "TCGA")

        prepared_fastq = Test.@test_logs (
            :warn,
            r"Unsupported record type in tuple"
        ) (
            :warn,
            r"Unsupported observation type"
        ) Mycelia.Rhizomorph._prepare_fastq_observations([
            fastq_record,
            fasta_record,
            (fasta_record, 1),
            (fastq_record, 2),
            ("bad", 3),
            7
        ])

        Test.@test length(prepared_fastq) == 4
        Test.@test all(record -> record isa FASTX.FASTQ.Record, prepared_fastq)
        Test.@test String(FASTX.sequence(prepared_fastq[2])) == "TCGA"
        Test.@test collect(UInt8, FASTX.quality_scores(prepared_fastq[2])) ==
                   fill(UInt8(40), 4)

        Test.@test Mycelia.Rhizomorph._determine_kmer_type([fastq_record], 3) ===
                   Mycelia.Kmers.DNAKmer{3}
        Test.@test Mycelia.Rhizomorph._determine_kmer_type(FASTX.FASTQ.Record[], 4) ===
                   Mycelia.Kmers.DNAKmer{4}

        fasta_reads = [FASTX.FASTA.Record("dna1", "ATCG")]
        Test.@test Mycelia.Rhizomorph._prepare_observations(fasta_reads) === fasta_reads
        auto_config = Mycelia.Rhizomorph._auto_configure_assembly(
            fasta_reads;
            min_overlap = 2
        )
        Test.@test auto_config.k === nothing
        Test.@test auto_config.graph_mode == Mycelia.Rhizomorph.DoubleStrand
        Test.@test !auto_config.use_quality_scores
    end

    Test.@testset "Assembly Strategy Helpers" begin
        observations = [
            FASTX.FASTA.Record("read1", "ATCG"),
            FASTX.FASTA.Record("read2", "TCGA")
        ]

        ngram_config = Mycelia.Rhizomorph.AssemblyConfig(
            k = 2,
            sequence_type = String,
            graph_mode = Mycelia.Rhizomorph.SingleStrand
        )
        ngram_result = Mycelia.Rhizomorph._assemble_ngram_graph(observations, ngram_config)
        Test.@test ngram_result.assembly_stats["method"] == "NgramGraph"
        Test.@test Mycelia.Rhizomorph.has_graph_structure(ngram_result)

        string_config = Mycelia.Rhizomorph.AssemblyConfig(
            min_overlap = 2,
            sequence_type = String,
            graph_mode = Mycelia.Rhizomorph.SingleStrand
        )
        string_result = Mycelia.Rhizomorph._assemble_string_graph(observations, string_config)
        Test.@test string_result.assembly_stats["method"] == "StringGraph"
        Test.@test !isempty(string_result.contigs)

        mktempdir() do dir
            gfa_path = joinpath(dir, "string_result.gfa")
            Test.@test_logs (
                :warn,
                "AssemblyResult is marked as not GFA compatible - output may be invalid"
            ) (
                :info,
                r"Assembly written to GFA format:"
            ) Mycelia.Rhizomorph.write_gfa(
                Mycelia.Rhizomorph.AssemblyResult(
                    string_result.contigs,
                    string_result.contig_names;
                    graph = string_result.graph,
                    simplified_graph = string_result.simplified_graph,
                    gfa_compatible = false
                ),
                gfa_path
            ) == gfa_path
        end

        kmer_config = Mycelia.Rhizomorph.AssemblyConfig(
            k = 3,
            sequence_type = BioSequences.LongDNA{4},
            graph_mode = Mycelia.Rhizomorph.SingleStrand,
            use_quality_scores = false
        )
        kmer_result = Mycelia.Rhizomorph._assemble_kmer_graph(observations, kmer_config)
        Test.@test kmer_result.assembly_stats["method"] == "KmerGraph"

        probabilistic_contigs = Mycelia.Rhizomorph._generate_contigs_probabilistic(
            kmer_result.graph,
            kmer_config
        )
        Test.@test probabilistic_contigs isa Vector{String}

        biosequence_config = Mycelia.Rhizomorph.AssemblyConfig(
            min_overlap = 2,
            sequence_type = BioSequences.LongDNA{4},
            graph_mode = Mycelia.Rhizomorph.DoubleStrand,
            use_quality_scores = false
        )
        biosequence_result = Mycelia.Rhizomorph._assemble_biosequence_graph(
            observations,
            biosequence_config
        )
        Test.@test biosequence_result.assembly_stats["method"] == "BioSequenceGraph"
        Test.@test !isempty(biosequence_result.contigs)

        quality_biosequence_config = Mycelia.Rhizomorph.AssemblyConfig(
            min_overlap = 2,
            sequence_type = BioSequences.LongDNA{4},
            graph_mode = Mycelia.Rhizomorph.SingleStrand,
            use_quality_scores = true
        )
        quality_biosequence_result = Mycelia.Rhizomorph._assemble_quality_biosequence_graph(
            observations,
            quality_biosequence_config
        )
        Test.@test quality_biosequence_result.assembly_stats["method"] ==
                   "QualityBioSequenceGraph"

        qualmer_config = Mycelia.Rhizomorph.AssemblyConfig(
            k = 2,
            sequence_type = BioSequences.LongDNA{4},
            graph_mode = Mycelia.Rhizomorph.SingleStrand,
            use_quality_scores = true
        )
        qualmer_result = Mycelia.Rhizomorph._assemble_qualmer_graph(observations, qualmer_config)
        Test.@test qualmer_result.assembly_stats["method"] == "QualmerGraph"
        Test.@test Mycelia.Rhizomorph.has_quality_information(qualmer_result)

        hybrid_result = Test.@test_logs (:warn, r"Hybrid OLC not fully implemented") Mycelia.Rhizomorph._assemble_hybrid_olc(
            observations,
            kmer_config
        )
        Test.@test hybrid_result.assembly_stats["method"] == "KmerGraph"

        multi_k_result = Test.@test_logs (:warn, r"Multi-k assembly not fully implemented") Mycelia.Rhizomorph._assemble_multi_k(
            observations,
            kmer_config
        )
        Test.@test multi_k_result.assembly_stats["method"] == "KmerGraph"
    end

    Test.@testset "Miscellaneous Assembly Helpers" begin
        verbose_config = Mycelia.Rhizomorph.AssemblyConfig(
            k = 3,
            sequence_type = BioSequences.LongDNA{4},
            graph_mode = Mycelia.Rhizomorph.SingleStrand,
            verbose = true
        )
        Test.@test_logs (:info, r"helper log") Mycelia.Rhizomorph._log_info(
            verbose_config,
            "helper log"
        )
        Test.@test isnothing(
            Mycelia.Rhizomorph._log_info(
                Mycelia.Rhizomorph.AssemblyConfig(
                    k = 3,
                    sequence_type = BioSequences.LongDNA{4},
                    graph_mode = Mycelia.Rhizomorph.SingleStrand,
                    verbose = false
                ),
                "silent"
            )
        )

        records = Mycelia.Rhizomorph._contigs_to_records(["ATCG", "GG"])
        Test.@test length(records) == 2
        Test.@test String(FASTX.FASTA.identifier(records[1])) == "contig_1"

        Test.@test Mycelia.Rhizomorph._calculate_n_statistic(Int[], 0.5) == 0
        Test.@test Mycelia.Rhizomorph._calculate_l_statistic(Int[], 0.5) == 0

        reference_metrics = Mycelia.Rhizomorph._validate_against_reference(["ATCG"], "ATCGA")
        Test.@test reference_metrics["reference_provided"]
        Test.@test reference_metrics["reference_length"] == 5

        passthrough_graph = build_mock_quality_graph()
        Test.@test Mycelia.Rhizomorph._simplify_string_graph(passthrough_graph) ===
                   passthrough_graph
        Test.@test Mycelia.Rhizomorph._simplify_ngram_graph(passthrough_graph) ===
                   passthrough_graph
        Test.@test Mycelia.Rhizomorph._simplify_ngram_to_string_graph(passthrough_graph) ===
                   passthrough_graph

        long_contig = repeat("A", 31)
        Test.@test_logs (:warn, r"Viterbi polishing failed") Mycelia.Rhizomorph._polish_contig_viterbi(
            long_contig,
            passthrough_graph,
            FASTX.FASTA.Record[]
        ) == long_contig
    end
end
