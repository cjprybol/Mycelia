import Test
import Mycelia
import MetaGraphsNext
import Graphs

struct MockLabelGraph{K, V}
    data::Dict{K, V}
end

Base.getindex(graph::MockLabelGraph, key) = graph.data[key]
Base.haskey(graph::MockLabelGraph, key) = haskey(graph.data, key)
MetaGraphsNext.labels(graph::MockLabelGraph) = collect(keys(graph.data))

struct NoDatasetVertexData end
struct MissingQualityVertexData end
struct EvidenceOnlyVertexData end
struct MockImprovedGraph end
struct MockViterbiPreferredGraph end
struct MockStatPreferredGraph end
struct MockStatDecodeGraph end

Mycelia.Rhizomorph.get_all_dataset_ids(::NoDatasetVertexData) = String[]
Mycelia.Rhizomorph.get_all_dataset_ids(::MissingQualityVertexData) = ["mock"]
Mycelia.Rhizomorph.get_all_dataset_ids(::EvidenceOnlyVertexData) = ["mock"]

Mycelia.Rhizomorph.get_vertex_joint_quality(::NoDatasetVertexData, ::String) = nothing
Mycelia.Rhizomorph.get_vertex_joint_quality(::MissingQualityVertexData, ::String) = nothing
Mycelia.Rhizomorph.get_vertex_joint_quality(::EvidenceOnlyVertexData, ::String) = UInt8[
    35, 37, 39]

Mycelia.Rhizomorph.get_vertex_mean_quality(::NoDatasetVertexData, ::String) = nothing
Mycelia.Rhizomorph.get_vertex_mean_quality(::MissingQualityVertexData, ::String) = Float64[]
Mycelia.Rhizomorph.get_vertex_mean_quality(::EvidenceOnlyVertexData, ::String) = Float64[
    21.0, 24.0, 27.0]

Mycelia.Rhizomorph.count_evidence(::NoDatasetVertexData) = 0
Mycelia.Rhizomorph.count_evidence(::MissingQualityVertexData) = 2
Mycelia.Rhizomorph.count_evidence(::EvidenceOnlyVertexData) = 4

function Mycelia.find_optimal_sequence_path(
        read::Mycelia.FASTX.FASTQ.Record, graph::MockImprovedGraph, k::Int;
        graph_mode::Symbol = :canonical)
    improved = Mycelia.FASTX.FASTQ.Record(
        Mycelia.FASTX.identifier(read),
        replace(Mycelia.FASTX.sequence(String, read), "A" => "T", count = 1),
        Mycelia.FASTX.quality(read)
    )
    return improved, 0.5
end

Mycelia.calculate_read_likelihood(read::Mycelia.FASTX.FASTQ.Record, ::MockViterbiPreferredGraph, k::Int;
    graph_mode::Symbol = :canonical) = 0.0
Mycelia.try_viterbi_path_improvement(read::Mycelia.FASTX.FASTQ.Record, ::MockViterbiPreferredGraph, k::Int;
    graph_mode::Symbol = :canonical) = (read, 0.5)
Mycelia.try_statistical_path_resampling(read::Mycelia.FASTX.FASTQ.Record, ::MockViterbiPreferredGraph, k::Int;
    graph_mode::Symbol = :canonical) = nothing
Mycelia.try_local_path_improvements(read::Mycelia.FASTX.FASTQ.Record, ::MockViterbiPreferredGraph, k::Int, original_likelihood::Float64;
    graph_mode::Symbol = :canonical) = (read, 0.0)

Mycelia.calculate_read_likelihood(read::Mycelia.FASTX.FASTQ.Record, ::MockStatPreferredGraph, k::Int;
    graph_mode::Symbol = :canonical) = 0.0
Mycelia.try_viterbi_path_improvement(read::Mycelia.FASTX.FASTQ.Record, ::MockStatPreferredGraph, k::Int;
    graph_mode::Symbol = :canonical) = nothing
Mycelia.try_statistical_path_resampling(read::Mycelia.FASTX.FASTQ.Record, ::MockStatPreferredGraph, k::Int;
    graph_mode::Symbol = :canonical) = (read, 0.5)
Mycelia.try_local_path_improvements(read::Mycelia.FASTX.FASTQ.Record, ::MockStatPreferredGraph, k::Int, original_likelihood::Float64;
    graph_mode::Symbol = :canonical) = (read, 0.0)

function Mycelia.generate_alternative_qualmer_paths(
        read::Mycelia.FASTX.FASTQ.Record, ::MockStatDecodeGraph, k::Int;
        num_samples::Int = 5, graph_mode::Symbol = :canonical)
    return [[Mycelia.Qualmer(Mycelia.Kmers.DNAKmer{3}("TTT"), UInt8[30, 30, 30])]]
end

function Mycelia.calculate_sequence_likelihood(
        sequence::Mycelia.BioSequences.BioSequence, quality::Vector{Int8},
        ::MockStatDecodeGraph, k::Int; graph_mode::Symbol = :canonical)
    return string(sequence) == "TTT" ? 1.0 : 0.1
end

Test.@testset "Iterative Assembly Helpers" begin
    Test.@testset "Sparsity and singleton detection" begin
        reads = [
            Mycelia.fastq_record(identifier = "r1", sequence = "AAAAAA", quality_scores = fill(30, 6)),
            Mycelia.fastq_record(identifier = "r2", sequence = "AT", quality_scores = fill(30, 2))
        ]
        sparsity = Mycelia.calculate_sparsity(reads, 2)
        Test.@test isapprox(sparsity, 0.875; atol = 1e-6)
        Test.@test Mycelia.errors_are_singletons(reads, 2)
        Test.@test !Mycelia.errors_are_singletons(Mycelia.FASTX.FASTQ.Record[], 3)
    end

    Test.@testset "K selection and memory estimates" begin
        reads = [Mycelia.fastq_record(identifier = "r1", sequence = "AAAAAA", quality_scores = fill(30, 6))]
        Test.@test Mycelia.find_initial_k(reads; k_range = [2, 3], sparsity_threshold = 0.0) ==
                   2
        Test.@test Mycelia.find_initial_k(Mycelia.FASTX.FASTQ.Record[]; k_range = Int[]) == 3
        Test.@test Mycelia.next_prime_k(5; max_k = 10) == 7
        Test.@test Mycelia.next_prime_k(11; max_k = 12) == 11

        estimate = Mycelia.estimate_memory_usage(10, 3)
        Test.@test estimate == 2436
    end

    Test.@testset "Label helpers" begin
        dna_qmer = Mycelia.Qualmer(Mycelia.Kmers.DNAKmer{3}("ATG"), UInt8[30, 31, 32])
        Test.@test Mycelia._label_k_length(dna_qmer) == 3
        Test.@test Mycelia._label_k_length("ATG") == 3
        Test.@test Mycelia._label_k_length(Mycelia.Kmers.DNAKmer{3}("ATG")) == 3
        Test.@test Mycelia._label_k_length(Mycelia.BioSequences.LongDNA{4}("ATG")) == 3
        Test.@test Mycelia._label_k_length(1234) == 4
    end

    Test.@testset "Memory limit checks" begin
        empty_graph = Graphs.SimpleGraph(0)
        Test.@test Mycelia.check_memory_limits(empty_graph, 1)

        records = [Mycelia.FASTX.FASTA.Record("seq1", "ATGC")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; mode = :singlestrand)
        labels = collect(MetaGraphsNext.labels(graph))
        expected = Mycelia.estimate_memory_usage(length(labels), 3)
        Test.@test Mycelia.check_memory_limits(graph, expected + 1)
        Test.@test !Mycelia.check_memory_limits(graph, expected - 1)
    end

    Test.@testset "Graph resolution and likelihood helpers" begin
        canonical_kmer = Mycelia.Kmers.DNAKmer{3}("ATG")
        reverse_kmer = Mycelia.Kmers.DNAKmer{3}("CAT")
        higher_prob_kmer = Mycelia.Kmers.DNAKmer{3}("ATA")
        rna_kmer = Mycelia.Kmers.RNAKmer{3}("AUG")
        aa_kmer = Mycelia.Kmers.AAKmer{2}("AR")

        dna_graph = MockLabelGraph(Dict(
            canonical_kmer => (joint_probability = 0.2, mean_quality = 20.0, coverage = 4),
            higher_prob_kmer => (joint_probability = 0.95, mean_quality = 35.0, coverage = 8)
        ))
        rna_graph = MockLabelGraph(Dict(
            rna_kmer => (joint_probability = 0.8, mean_quality = 31.0, coverage = 5)
        ))
        aa_graph = MockLabelGraph(Dict(
            aa_kmer => (joint_probability = 0.75, mean_quality = 28.0, coverage = 3)
        ))

        Test.@test Mycelia._resolve_kmer_label(
            dna_graph, reverse_kmer; graph_mode = :canonical) == canonical_kmer
        Test.@test isnothing(Mycelia._resolve_kmer_label(
            dna_graph, Mycelia.Kmers.DNAKmer{3}("CCC"); graph_mode = :singlestrand))

        reverse_qmer = Mycelia.Qualmer(reverse_kmer, UInt8[30, 31, 32])
        resolved_qmer = Mycelia._resolve_qualmer_for_graph(
            dna_graph, reverse_qmer; graph_mode = :canonical)
        Test.@test !isnothing(resolved_qmer)
        Test.@test resolved_qmer.kmer == canonical_kmer
        Test.@test isnothing(Mycelia._resolve_qualmer_for_graph(
            dna_graph,
            Mycelia.Qualmer(Mycelia.Kmers.DNAKmer{3}("CCC"), UInt8[20, 20, 20]);
            graph_mode = :canonical))

        short_seq = Mycelia.BioSequences.LongDNA{4}("AT")
        Test.@test Mycelia.calculate_sequence_likelihood(
            short_seq, Int8[30, 30], dna_graph, 3; graph_mode = :canonical) == 0.0

        rna_likelihood = Mycelia.calculate_sequence_likelihood(
            "AUGA", Int8[40, 40, 40, 40], rna_graph, 3; graph_mode = :singlestrand)
        aa_likelihood = Mycelia.calculate_sequence_likelihood(
            "ARNA", Int8[35, 35, 35, 35], aa_graph, 2; graph_mode = :singlestrand)
        Test.@test isa(rna_likelihood, Float64)
        Test.@test isa(aa_likelihood, Float64)

        penalty_low = Mycelia.calculate_unseen_kmer_penalty(Int8[-5, -2, 0])
        penalty_high = Mycelia.calculate_unseen_kmer_penalty(Int8[40, 40, 40])
        Test.@test penalty_high > penalty_low > 0.0

        qualmer_likelihood = Mycelia.calculate_qualmer_likelihood(
            canonical_kmer, Int8[-1, 30, 35], dna_graph[canonical_kmer])
        Test.@test qualmer_likelihood > 0.0
    end

    Test.@testset "Vertex helper fallbacks" begin
        Test.@test Mycelia._rhizomorph_joint_probability(NoDatasetVertexData()) == 0.0
        Test.@test Mycelia._rhizomorph_joint_probability(MissingQualityVertexData()) == 0.0
        Test.@test Mycelia._vertex_mean_quality(NoDatasetVertexData()) == 0.0
        Test.@test Mycelia._vertex_mean_quality(MissingQualityVertexData()) == 0.0
        Test.@test Mycelia._vertex_coverage(EvidenceOnlyVertexData()) == 4

        fallback_prob = Mycelia._rhizomorph_joint_probability(EvidenceOnlyVertexData())
        Test.@test 0.0 < fallback_prob <= 1.0
        Test.@test isapprox(
            Mycelia._vertex_mean_quality(EvidenceOnlyVertexData()), 24.0; atol = 1e-6)

        property_data = (joint_probability = 0.7, mean_quality = 12.5, coverage = [1, 2, 3])
        Test.@test Mycelia._vertex_joint_probability(property_data) == 0.7
        Test.@test Mycelia._vertex_mean_quality(property_data) == 12.5
        Test.@test Mycelia._vertex_coverage(property_data) == 3

        property_int_coverage = (joint_probability = 0.6, mean_quality = 10.0, coverage = 9)
        Test.@test Mycelia._vertex_coverage(property_int_coverage) == 9
    end

    Test.@testset "Alternative generation and path conversions" begin
        dna_graph = MockLabelGraph(Dict(
            Mycelia.Kmers.DNAKmer{3}("ATG") => (joint_probability = 0.1, mean_quality = 20.0, coverage = 3),
            Mycelia.Kmers.DNAKmer{3}("ATA") => (joint_probability = 0.95, mean_quality = 35.0, coverage = 5)
        ))
        rna_graph = MockLabelGraph(Dict(
            Mycelia.Kmers.RNAKmer{3}("AUG") => (joint_probability = 0.2, mean_quality = 20.0, coverage = 3),
            Mycelia.Kmers.RNAKmer{3}("UGA") => (joint_probability = 0.8, mean_quality = 33.0, coverage = 4),
            Mycelia.Kmers.RNAKmer{3}("AUA") => (joint_probability = 0.9, mean_quality = 32.0, coverage = 4)
        ))
        aa_graph = MockLabelGraph(Dict(
            Mycelia.Kmers.AAKmer{2}("AR") => (joint_probability = 0.1, mean_quality = 19.0, coverage = 2),
            Mycelia.Kmers.AAKmer{2}("AA") => (joint_probability = 0.92, mean_quality = 40.0, coverage = 6)
        ))

        dna_alternatives = Mycelia.generate_kmer_alternatives(
            Mycelia.Kmers.DNAKmer{3}("ATG"), dna_graph; graph_mode = :singlestrand)
        rna_object_alternatives = Mycelia.generate_kmer_alternatives(
            Mycelia.Kmers.RNAKmer{3}("AUG"), rna_graph; graph_mode = :singlestrand)
        aa_object_alternatives = Mycelia.generate_kmer_alternatives(
            Mycelia.Kmers.AAKmer{2}("AR"), aa_graph; graph_mode = :singlestrand)
        rna_alternatives = Mycelia.generate_kmer_alternatives(
            "AUG", rna_graph; graph_mode = :singlestrand)
        aa_alternatives = Mycelia.generate_kmer_alternatives(
            "AR", aa_graph; graph_mode = :singlestrand)

        Test.@test Mycelia.Kmers.DNAKmer{3}("ATA") in dna_alternatives
        Test.@test Mycelia.Kmers.RNAKmer{3}("AUA") in rna_object_alternatives
        Test.@test Mycelia.Kmers.AAKmer{2}("AA") in aa_object_alternatives
        Test.@test Mycelia.Kmers.RNAKmer{3}("AUA") in rna_alternatives
        Test.@test Mycelia.Kmers.AAKmer{2}("AA") in aa_alternatives

        Test.@test isempty(Mycelia.generate_alternative_paths("ATG", dna_graph, 3; num_samples = 0))
        Test.@test Mycelia.kmer_path_to_sequence(Any["", "A"], 3) == ""

        sequence_path = Mycelia.sequence_to_qualmer_path("AUGA", "IIII", rna_graph, 3;
            graph_mode = :singlestrand)
        Test.@test length(sequence_path) == 2
        Test.@test isempty(Mycelia.sequence_to_qualmer_path(
            "AT", "II", MockLabelGraph(Dict{Mycelia.Kmers.DNAKmer{3}, Any}()), 3))

        empty_qualmer_result = Mycelia.qualmer_path_to_biosequence(Mycelia.Qualmer[])
        Test.@test isempty(empty_qualmer_result.sequence)
        Test.@test isempty(empty_qualmer_result.quality_scores)

        rna_qualmer_result = Mycelia.qualmer_path_to_biosequence([
            Mycelia.Qualmer(Mycelia.Kmers.RNAKmer{3}("AUG"), UInt8[30, 31, 32]),
            Mycelia.Qualmer(Mycelia.Kmers.RNAKmer{3}("UGA"), UInt8[31, 32, 33])
        ])
        Test.@test rna_qualmer_result.sequence isa Mycelia.BioSequences.LongRNA{4}
        Test.@test length(rna_qualmer_result.quality_scores) == 4

        aa_qualmer_result = Mycelia.qualmer_path_to_biosequence([
            Mycelia.Qualmer(Mycelia.Kmers.AAKmer{2}("AR"), UInt8[20, 21]),
            Mycelia.Qualmer(Mycelia.Kmers.AAKmer{2}("RN"), UInt8[21, 22])
        ])
        Test.@test aa_qualmer_result.sequence isa Mycelia.BioSequences.LongAA

        rna_kmer_result = Mycelia.kmer_path_to_biosequence([
            Mycelia.Kmers.RNAKmer{3}("AUG"),
            Mycelia.Kmers.RNAKmer{3}("UGA")
        ])
        aa_kmer_result = Mycelia.kmer_path_to_biosequence([
            Mycelia.Kmers.AAKmer{2}("AR"),
            Mycelia.Kmers.AAKmer{2}("RN")
        ])
        empty_dna_kmer_result = Mycelia.kmer_path_to_biosequence(Mycelia.Kmers.DNAKmer{3}[])
        dna_kmer_result = Mycelia.kmer_path_to_biosequence([
            Mycelia.Kmers.DNAKmer{3}("ATG"),
            Mycelia.Kmers.DNAKmer{3}("TGC")
        ])
        Test.@test rna_kmer_result isa Mycelia.BioSequences.LongRNA{4}
        Test.@test aa_kmer_result isa Mycelia.BioSequences.LongAA
        Test.@test empty_dna_kmer_result == Mycelia.BioSequences.LongDNA{4}("")
        Test.@test dna_kmer_result isa Mycelia.BioSequences.LongDNA{4}
    end

    Test.@testset "Decision and finalization helpers" begin
        plateau_history = [
            Dict{Symbol, Any}(:improvement_rate => 0.02),
            Dict{Symbol, Any}(:improvement_rate => 0.01),
            Dict{Symbol, Any}(:improvement_rate => 0.02)
        ]
        decreasing_history = [
            Dict{Symbol, Any}(:improvement_rate => 0.08),
            Dict{Symbol, Any}(:improvement_rate => 0.07),
            Dict{Symbol, Any}(:improvement_rate => 0.06)
        ]
        Test.@test !Mycelia.sufficient_improvements(
            8, 100, 0.05; iteration_history = decreasing_history)
        Test.@test !Mycelia.sufficient_improvements(
            8, 100, 0.05; iteration_history = plateau_history)
        Test.@test Mycelia.sufficient_improvements(
            12, 100, 0.05;
            iteration_history = [Dict{Symbol, Any}(:improvement_rate => 0.12)])

        k_history = [
            Dict{Symbol, Any}(:total_improvements => 0, :total_reads => 10, :runtime_seconds => 2.0),
            Dict{Symbol, Any}(:total_improvements => 0, :total_reads => 10, :runtime_seconds => 2.0),
            Dict{Symbol, Any}(:total_improvements => 0, :total_reads => 10, :runtime_seconds => 2.0)
        ]
        progressing_history = [
            Dict{Symbol, Any}(:total_improvements => 5, :total_reads => 10, :runtime_seconds => 0.2),
            Dict{Symbol, Any}(:total_improvements => 6, :total_reads => 10, :runtime_seconds => 0.2),
            Dict{Symbol, Any}(:total_improvements => 7, :total_reads => 10, :runtime_seconds => 0.2)
        ]

        Test.@test !Mycelia.should_continue_k_progression(progressing_history, 11, 11)
        Test.@test Mycelia.should_continue_k_progression([progressing_history[1]], 5, 11)
        Test.@test !Mycelia.should_continue_k_progression(k_history, 5, 11)
        Test.@test Mycelia.should_continue_k_progression(progressing_history, 5, 11)

        temp_dir = mktempdir()
        iteration_history = Dict(
            7 => [Dict{Symbol, Any}(:timestamp => "missing", :improvements_made => 0)]
        )
        summary_result = redirect_stdout(devnull) do
            Mycelia.finalize_iterative_assembly(temp_dir, [7], iteration_history, 1.5; verbose = true)
        end
        Test.@test isempty(summary_result[:final_assembly])
        Test.@test summary_result[:metadata][:final_fastq_file] ==
                   joinpath(temp_dir, "reads_k7_iter1_missing.fastq")
        rm(temp_dir; recursive = true, force = true)
    end

    Test.@testset "Quality adapters and batch processing" begin
        records = [
            Mycelia.fastq_record(identifier = "r1", sequence = "ATCGATCG", quality_scores = fill(40, 8)),
            Mycelia.fastq_record(identifier = "r2", sequence = "ATCGATCA", quality_scores = fill(40, 8))
        ]
        graph = Mycelia.Rhizomorph.build_qualmer_graph(records, 3; mode = :canonical)
        quality_vector = fill(Int8(40), 8)

        string_result = Mycelia.find_optimal_sequence_path(
            "ATCGATCG", quality_vector, graph, 3; graph_mode = :canonical)
        viterbi_result = Mycelia.try_viterbi_path_improvement(
            "ATCGATCG", quality_vector, graph, 3; graph_mode = :canonical)
        statistical_result = Mycelia.try_statistical_path_resampling(
            "ATCGATCG", quality_vector, graph, 3; graph_mode = :canonical)
        local_result = Mycelia.try_local_path_improvements(
            "ATCGATCG", quality_vector, graph, 3, 0.0; graph_mode = :canonical)

        Test.@test isa(string_result, Tuple{String, Float64})
        Test.@test isa(viterbi_result, Union{Tuple{String, Float64}, Nothing})
        Test.@test isa(statistical_result, Union{Tuple{String, Float64}, Nothing})
        Test.@test isa(local_result, Tuple{String, Float64})

        short_read = Mycelia.fastq_record(identifier = "short", sequence = "AT", quality_scores = fill(30, 2))
        returned_read, improved = Mycelia.improve_read_likelihood(short_read, graph, 3)
        Test.@test Mycelia.FASTX.identifier(returned_read) == "short"
        Test.@test !improved

        batch_reads, batch_improvements = redirect_stdout(devnull) do
            Mycelia.improve_read_set_likelihood(
                records, graph, 3;
                verbose = true,
                batch_size = 1,
                enable_parallel = Threads.nthreads() > 1,
                graph_mode = :canonical
            )
        end
        Test.@test length(batch_reads) == length(records)
        Test.@test batch_improvements >= 0

        empty_extended = Mycelia.adjust_quality_scores("", "ATCG", 0.5)
        Test.@test empty_extended == "IIII"
    end

    Test.@testset "Mock-driven branch coverage" begin
        read = Mycelia.fastq_record(
            identifier = "mock",
            sequence = "AAAT",
            quality_scores = fill(35, 4)
        )

        improved_read, was_improved = Mycelia.improve_read_likelihood(
            read, MockImprovedGraph(), 3)
        Test.@test was_improved
        Test.@test Mycelia.FASTX.sequence(String, improved_read) !=
                   Mycelia.FASTX.sequence(String, read)

        viterbi_choice = Mycelia.find_optimal_sequence_path(
            read, MockViterbiPreferredGraph(), 3)
        statistical_choice = Mycelia.find_optimal_sequence_path(
            read, MockStatPreferredGraph(), 3)
        Test.@test isa(viterbi_choice, Tuple{Mycelia.FASTX.FASTQ.Record, Float64})
        Test.@test isa(statistical_choice, Tuple{Mycelia.FASTX.FASTQ.Record, Float64})

        wrapped_viterbi = Mycelia.try_viterbi_path_improvement(
            "AAAT", fill(Int8(35), 4), MockViterbiPreferredGraph(), 3)
        Test.@test isa(wrapped_viterbi, Tuple{String, Float64})

        direct_stat = Mycelia.try_statistical_path_resampling(
            Mycelia.fastq_record(identifier = "stat", sequence = "AAA", quality_scores = fill(30, 3)),
            MockStatDecodeGraph(),
            3
        )
        wrapped_stat = Mycelia.try_statistical_path_resampling(
            "AAA", fill(Int8(30), 3), MockStatPreferredGraph(), 3)
        Test.@test isa(direct_stat, Tuple{Mycelia.FASTX.FASTQ.Record, Float64})
        Test.@test isa(wrapped_stat, Tuple{String, Float64})
    end
end
