import Test
import Mycelia
import FASTX
import BioSequences
import MetaGraphsNext
import StableRNGs

const BUILD_READ_LENGTH = 256
const BUILD_K = 15
const BASES = UInt8[0x41, 0x43, 0x47, 0x54]
const BENCHMARK_SINK = Base.RefValue{Any}(nothing)

function synthetic_dna_sequence(sequence_length::Int, seed::UInt32)
    data = Vector{UInt8}(undef, sequence_length)
    rng = StableRNGs.StableRNG(seed)

    for i in 1:sequence_length
        data[i] = rand(rng, BASES)
    end

    return BioSequences.LongDNA{4}(String(data))
end

function make_fasta_records(record_count::Int, record_length::Int)
    return [FASTX.FASTA.Record(
                "record_$(index)",
                synthetic_dna_sequence(record_length, UInt32(index))
            )
            for index in 1:record_count]
end

function build_unique_path_graph(sequence_length::Int; k::Int = BUILD_K, max_attempts::Int = 256)
    expected_vertex_count = sequence_length - k + 1
    base_seed = UInt32(sequence_length + 42)

    for attempt in 0:(max_attempts - 1)
        record = FASTX.FASTA.Record(
            "path_sequence_$(sequence_length)_$(attempt)",
            synthetic_dna_sequence(sequence_length, base_seed + UInt32(attempt))
        )

        graph = Mycelia.Rhizomorph.build_kmer_graph(
            [record],
            k;
            dataset_id = "performance_path_graph_$(sequence_length)_$(attempt)",
            mode = :singlestrand
        )

        vertex_count = length(MetaGraphsNext.labels(graph))
        if vertex_count == expected_vertex_count
            return graph, vertex_count
        end
    end

    error("Unable to synthesize a collision-free path graph for sequence length $(sequence_length)")
end

function best_elapsed_seconds(operation::Function; samples::Int = 3, repetitions::Int = 1)
    best = Inf

    for _ in 1:samples
        GC.gc()
        start_ns = Base.time_ns()
        for _ in 1:repetitions
            BENCHMARK_SINK[] = operation()
        end
        best = min(best, (Base.time_ns() - start_ns) / repetitions / 1e9)
    end

    return best
end

function measure_kmer_graph_build(record_count::Int; record_length::Int = BUILD_READ_LENGTH, k::Int = BUILD_K)
    records = make_fasta_records(record_count, record_length)

    Mycelia.Rhizomorph.build_kmer_graph(
        records,
        k;
        dataset_id = "performance_build_warmup",
        mode = :singlestrand
    )

    elapsed_seconds = best_elapsed_seconds(
        () -> Mycelia.Rhizomorph.build_kmer_graph(
            records,
            k;
            dataset_id = "performance_build_measure",
            mode = :singlestrand
        );
        samples = 3,
        repetitions = 3
    )

    graph = nothing
    GC.gc()
    allocated_bytes = @allocated begin
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records,
            k;
            dataset_id = "performance_build_allocations",
            mode = :singlestrand
        )
    end

    return (
        record_count = record_count,
        total_bases = record_count * record_length,
        elapsed_seconds = elapsed_seconds,
        allocated_bytes = allocated_bytes,
        vertex_count = length(MetaGraphsNext.labels(graph)),
        edge_count = length(MetaGraphsNext.edge_labels(graph))
    )
end

function measure_eulerian_path_finding(sequence_length::Int; k::Int = BUILD_K)
    graph, vertex_count = build_unique_path_graph(sequence_length; k = k)

    Mycelia.Rhizomorph.find_eulerian_paths_next(graph)

    elapsed_seconds = best_elapsed_seconds(
        () -> Mycelia.Rhizomorph.find_eulerian_paths_next(graph);
        samples = 3,
        repetitions = 10
    )

    paths = nothing
    GC.gc()
    allocated_bytes = @allocated begin
        paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
    end

    return (
        sequence_length = sequence_length,
        elapsed_seconds = elapsed_seconds,
        allocated_bytes = allocated_bytes,
        vertex_count = vertex_count,
        edge_count = length(MetaGraphsNext.edge_labels(graph)),
        path_count = length(paths),
        first_path_length = isempty(paths) ? 0 : length(first(paths))
    )
end

function normalized_metric_range(metrics, numerator::Symbol, denominator::Symbol)
    normalized = [getfield(metric, numerator) / getfield(metric, denominator)
                  for metric in metrics
                  if getfield(metric, numerator) > 0 &&
                     getfield(metric, denominator) > 0]

    if isempty(normalized)
        return 1.0
    end

    return maximum(normalized) / minimum(normalized)
end

Test.@testset "Assembly Performance Bounds" begin
    Test.@testset "K-mer graph construction scales linearly with input bases" begin
        build_metrics = [measure_kmer_graph_build(record_count)
                         for record_count in (8, 16, 32)]

        for metric in build_metrics
            Test.@test metric.vertex_count > 0
            Test.@test metric.edge_count > 0
            Test.@test metric.allocated_bytes / metric.total_bases < 6_000
        end

        Test.@test issorted([metric.vertex_count for metric in build_metrics])
        Test.@test issorted([metric.edge_count for metric in build_metrics])
        Test.@test normalized_metric_range(build_metrics, :elapsed_seconds, :total_bases) <
                   4.0
        Test.@test normalized_metric_range(build_metrics, :allocated_bytes, :total_bases) <
                   2.0
    end

    Test.@testset "Eulerian path finding stays within linear graph-size bounds" begin
        path_metrics = [measure_eulerian_path_finding(sequence_length)
                        for sequence_length in (1024, 2048, 4096)]

        for metric in path_metrics
            Test.@test metric.vertex_count > 0
            Test.@test metric.edge_count > 0
            Test.@test metric.path_count == 1
            Test.@test metric.first_path_length == metric.vertex_count
            Test.@test metric.allocated_bytes / metric.vertex_count < 1_500
        end

        Test.@test issorted([metric.vertex_count for metric in path_metrics])
        Test.@test normalized_metric_range(path_metrics, :elapsed_seconds, :vertex_count) <
                   3.0
        Test.@test normalized_metric_range(path_metrics, :allocated_bytes, :vertex_count) <
                   3.0
    end
end
