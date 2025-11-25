# # Data Processing Benchmark
#
# This benchmark evaluates the performance of data processing operations
# using realistic datasets and high-throughput scenarios.
#
# ## Benchmark Overview
#
# This benchmark tests:
# - Large-scale FASTQ processing
# - File format conversion performance
# - Memory usage with big datasets
# - I/O performance optimization
# - Parallel processing scalability

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Test
import Mycelia
import FASTX
import BenchmarkTools
import Random
import Statistics
import Dates
import BioSequences
import CSV
import DataFrames

# Include benchmark utilities
include("benchmark_utils.jl")

println("=== Data Processing Benchmark ===")
println("Start time: $(now())")

# ## Benchmark Configuration
#
# Configure benchmark parameters for different scales

# Small scale (development testing)
small_config = Dict(
    "n_files" => 5,
    "reads_per_file" => 10000,
    "read_length" => 150,
    "description" => "Small scale - quick validation"
)

# Medium scale (realistic testing)
medium_config = Dict(
    "n_files" => 50,
    "reads_per_file" => 100000,
    "read_length" => 15000,
    "description" => "Medium scale - realistic datasets"
)

# Large scale (scalability testing)
large_config = Dict(
    "n_files" => 100,
    "reads_per_file" => 1000000,
    "read_length" => 20000,
    "description" => "Large scale - scalability testing"
)

# Select configuration based on environment
config = if haskey(ENV, "BENCHMARK_SCALE")
    if ENV["BENCHMARK_SCALE"] == "small"
        small_config
    elseif ENV["BENCHMARK_SCALE"] == "medium"
        medium_config
    else
        large_config
    end
else
    small_config  # Default to small for safety
end

println("Benchmark Configuration: $(config["description"])")
println("Files: $(config["n_files"]), Reads per file: $(config["reads_per_file"])")

# ## Data Generation
#
# Generate realistic test datasets using NCBI reference genomes and read simulation

println("\n--- Generating Test Data ---")

# Generate test datasets using proper read simulation from reference genomes
test_data_dir = "test_data"
mkpath(test_data_dir)

# Create a temporary reference genome for read simulation
# Using random_fasta_record for reproducible test data generation
ref_fasta_path = joinpath(test_data_dir, "reference_genome.fasta")
ref_record = Mycelia.random_fasta_record(moltype=:DNA, seed=42, L=10000)
Mycelia.write_fasta(outfile=ref_fasta_path, records=[ref_record])

test_files = []
for i in 1:min(5, config["n_files"])  # Limit to 5 files for testing
    # Use simulate_illumina_reads for realistic FASTQ generation
    result = Mycelia.simulate_illumina_reads(
        fasta=ref_fasta_path, 
        read_count=config["reads_per_file"],
        len=config["read_length"],
        quiet=true
    )
    # Rename to our expected pattern
    filename = joinpath(test_data_dir, "test_reads_$(i).fastq")
    if isfile(result.fastq1)
        cp(result.fastq1, filename, force=true)
        push!(test_files, filename)
    end
end

println("Generated $(length(test_files)) test FASTQ files")

# ## FASTQ Processing Benchmarks
#
# Test FASTQ reading and processing performance

println("\n--- FASTQ Processing Benchmarks ---")

# Initialize benchmark suite
benchmark_suite = BenchmarkSuite("Data Processing Performance")

# Benchmark FASTQ record counting
for (i, test_file) in enumerate(test_files[1:3])  # Test first 3 files
    println("\nBenchmarking FASTQ record counting: $(basename(test_file))")
    
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.count_records, test_file;
        samples=3, seconds=10
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "count_records_file$(i)", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# Benchmark FASTQ to normalized table conversion
if length(test_files) > 0
    println("\nBenchmarking FASTQ to normalized table conversion")
    
    test_file = test_files[1]
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.fastx2normalized_table, test_file;
        samples=2, seconds=15
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "fastq_to_table_conversion", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# ## Format Conversion Benchmarks
#
# Test file format conversion performance

println("\n--- Format Conversion Benchmarks ---")

# Benchmark read length determination
if length(test_files) > 0
    println("\nBenchmarking read length determination")
    
    test_file = test_files[1]
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.determine_read_lengths, test_file;
        samples=3, seconds=10
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "determine_read_lengths", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# Benchmark duplication rate assessment
if length(test_files) > 0
    println("\nBenchmarking duplication rate assessment")
    
    test_file = test_files[1]
    
    # This function processes the file twice, so it's more intensive
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.assess_duplication_rates, test_file;
        samples=2, seconds=20
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "assess_duplication_rates", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# ## Memory Usage Benchmarks
#
# Test memory efficiency with large datasets

println("\n--- Memory Usage Benchmarks ---")

# Memory scaling analysis
println("\n--- Memory Scaling Analysis ---")

# Test how memory usage scales with file size
file_sizes = []
memory_usage = []

for test_file in test_files[1:3]
    file_size = filesize(test_file)
    push!(file_sizes, file_size)
    
    # Estimate memory requirement for processing
    memory_estimate = estimate_memory_requirement(Mycelia.fastx2normalized_table, test_file)
    push!(memory_usage, memory_estimate)
    
    println("File $(basename(test_file)):")
    println("  Size: $(round(file_size / 1e6, digits=2)) MB")
    println("  Estimated memory: $(round(memory_estimate / 1e6, digits=2)) MB")
    println("  Memory ratio: $(round(memory_estimate / file_size, digits=2))")
end

# ## I/O Performance Benchmarks
#
# Test input/output performance

println("\n--- I/O Performance Benchmarks ---")

# I/O Performance Analysis
println("\n--- I/O Performance Analysis ---")

# Test sequential file reading performance
if length(test_files) > 0
    test_file = test_files[1]
    
    println("Testing sequential FASTQ reading performance")
    
    # Benchmark raw file reading
    function read_fastq_sequential(filename)
        record_count = 0
        open(FASTX.FASTQ.Reader, filename) do reader
            for record in reader
                record_count += 1
            end
        end
        return record_count
    end
    
    benchmark_result, memory_stats = run_benchmark_with_memory(
        read_fastq_sequential, test_file;
        samples=3, seconds=10
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "sequential_fastq_reading", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# Save benchmark results
results_dir = "results"
mkpath(results_dir)
results_file = joinpath(results_dir, "data_processing_benchmark_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")
save_benchmark_results(benchmark_suite, results_file)

# Display summary
format_benchmark_summary(benchmark_suite)

# Cleanup test data
for test_file in test_files
    rm(test_file, force=true)
end
rm(test_data_dir, force=true, recursive=true)

println("\n=== Data Processing Benchmark Complete ===")
println("Results saved to: $results_file")
println("End time: $(Dates.now())")
# - Network vs local storage
# - Concurrent I/O operations

# ## Parallel Processing Benchmarks
#
# Test scalability with multiple threads/processes

println("\n--- Parallel Processing Benchmarks ---")

# Test thread scaling for parallel FASTQ processing
if length(test_files) >= 3
    println("\nBenchmarking parallel FASTQ processing")
    
    # Function to process multiple files in parallel
    function process_files_parallel(files)
        results = []
        Threads.@threads for file in files
            result = Mycelia.count_records(file)
            push!(results, result)
        end
        return results
    end
    
    # Test with different thread counts
    test_files_subset = test_files[1:min(3, length(test_files))]
    
    for thread_count in [1, min(2, Sys.CPU_THREADS), min(4, Sys.CPU_THREADS)]
        if thread_count <= Sys.CPU_THREADS
            println("Testing with $thread_count threads")
            
            # Set thread limit (approximation - Julia doesn't allow runtime thread changes)
            ENV["JULIA_NUM_THREADS"] = string(thread_count)
            
            benchmark_result, memory_stats = run_benchmark_with_memory(
                process_files_parallel, test_files_subset;
                samples=2, seconds=15
            )
            
            add_benchmark_result!(
                benchmark_suite, 
                "parallel_processing_$(thread_count)_threads", 
                benchmark_result, 
                memory_stats
            )
            
            println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
            println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
        end
    end
end

# Test load balancing with uneven file sizes
if length(test_files) >= 2
    println("\nBenchmarking load balancing with uneven workloads")
    
    # Create files of different sizes for load balancing test
    small_file = test_files[1]
    large_file_path = joinpath(test_data_dir, "large_test.fastq")
    
    # Generate a larger test file using simulate_illumina_reads
    large_result = Mycelia.simulate_illumina_reads(
        fasta=ref_fasta_path,
        read_count=config["reads_per_file"] * 3,
        len=config["read_length"],
        quiet=true
    )
    if isfile(large_result.fastq1)
        cp(large_result.fastq1, large_file_path, force=true)
    end
    
    mixed_files = [small_file, large_file_path]
    
    benchmark_result, memory_stats = run_benchmark_with_memory(
        process_files_parallel, mixed_files;
        samples=2, seconds=20
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "load_balancing_mixed_sizes", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Load balancing - Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    
    # Cleanup
    rm(large_file_path, force=true)
end

# ## Results Collection
#
# Collect and save benchmark results

println("\n--- Collecting Results ---")

# Calculate aggregate performance metrics
total_files_processed = length(test_files)
total_data_size = sum(filesize(f) for f in test_files if isfile(f))

# Calculate throughput metrics
if haskey(benchmark_suite.results, "count_records_file1")
    sample_time = benchmark_suite.results["count_records_file1"]["median_time"] / 1e9  # Convert to seconds
    sample_file_size = filesize(test_files[1])
    throughput_mbps = (sample_file_size / 1e6) / sample_time  # MB/s
else
    throughput_mbps = 0.0
end

# Enhanced results structure with comprehensive metrics
enhanced_results = Dict(
    "benchmark_name" => "data_processing",
    "timestamp" => Dates.now(),
    "configuration" => config,
    "system_info" => Dict(
        "julia_version" => VERSION,
        "cpu_threads" => Sys.CPU_THREADS,
        "total_memory" => Sys.total_memory(),
        "hostname" => gethostname(),
        "platform" => Sys.MACHINE
    ),
    "performance_summary" => Dict(
        "total_files_processed" => total_files_processed,
        "total_data_size_mb" => round(total_data_size / 1e6, digits=2),
        "throughput_mbps" => round(throughput_mbps, digits=2),
        "benchmark_duration_minutes" => round((Dates.now() - benchmark_suite.metadata["timestamp"]).value / 60000, digits=2)
    ),
    "detailed_results" => benchmark_suite.results,
    "resource_utilization" => Dict(
        "peak_memory_estimate_mb" => round(maximum([get(result, "memory", 0) for result in values(benchmark_suite.results)]) / 1e6, digits=2),
        "total_allocations" => sum([get(result, "allocations", 0) for result in values(benchmark_suite.results)])
    )
)

# Save enhanced results (this replaces the duplicate save at the end)
enhanced_results_file = joinpath(results_dir, "data_processing_enhanced_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")

open(enhanced_results_file, "w") do f
    JSON.print(f, enhanced_results, 2)
end

println("\n--- Performance Summary ---")
println("Files processed: $total_files_processed")
println("Total data size: $(round(total_data_size / 1e6, digits=2)) MB")
println("Estimated throughput: $(round(throughput_mbps, digits=2)) MB/s")
println("Peak memory usage: $(round(maximum([get(result, "memory", 0) for result in values(benchmark_suite.results)]) / 1e6, digits=2)) MB")

println("\nDetailed results saved to: $enhanced_results_file")
println("End time: $(Dates.now())")

nothing