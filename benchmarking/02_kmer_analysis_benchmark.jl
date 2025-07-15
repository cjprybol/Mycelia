# # K-mer Analysis Benchmark
#
# This benchmark evaluates k-mer analysis performance across different
# scales and parameters, testing both accuracy and computational efficiency.
#
# ## Benchmark Overview
#
# This benchmark tests:
# - K-mer counting performance at different k values
# - Dense vs sparse counting efficiency
# - Memory usage scaling with k-mer space
# - Parallel k-mer processing
# - K-mer spectrum analysis performance

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

# Include benchmark utilities
include("benchmark_utils.jl")

println("=== K-mer Analysis Benchmark ===")
println("Start time: $(now())")

# ## Benchmark Configuration
#
# Configure k-mer analysis parameters for different scales

# Small scale configuration
small_config = Dict(
    "genome_sizes" => [10000, 50000, 100000],
    "k_values" => [11, 15, 19, 21],
    "n_genomes" => 5,
    "coverage" => 10,
    "description" => "Small scale - quick validation"
)

# Medium scale configuration
medium_config = Dict(
    "genome_sizes" => [100000, 500000, 1000000],
    "k_values" => [15, 19, 21, 25, 31],
    "n_genomes" => 20,
    "coverage" => 30,
    "description" => "Medium scale - realistic datasets"
)

# Large scale configuration
large_config = Dict(
    "genome_sizes" => [1000000, 5000000, 10000000],
    "k_values" => [19, 21, 25, 31, 35],
    "n_genomes" => 50,
    "coverage" => 50,
    "description" => "Large scale - scalability testing"
)

# Select configuration
config = if haskey(ENV, "BENCHMARK_SCALE")
    if ENV["BENCHMARK_SCALE"] == "small"
        small_config
    elseif ENV["BENCHMARK_SCALE"] == "medium"
        medium_config
    else
        large_config
    end
else
    small_config
end

println("Benchmark Configuration: $(config["description"])")
println("Genome sizes: $(config["genome_sizes"])")
println("K values: $(config["k_values"])")
println("Number of genomes: $(config["n_genomes"])")

# ## Test Data Generation
#
# Generate test datasets for k-mer analysis

println("\n--- Generating Test Data ---")

# Generate test data using existing Mycelia simulation functions
test_sequences = Dict()
for genome_size in config["genome_sizes"]
    test_sequences[genome_size] = Mycelia.generate_test_sequences(genome_size, config["n_genomes"])
end

println("Generated test sequences for genome sizes: $(config["genome_sizes"])")

# ## Dense K-mer Counting Benchmarks
#
# Test dense k-mer counting performance

println("\n--- Dense K-mer Counting Benchmarks ---")

# Initialize benchmark suite
benchmark_suite = BenchmarkSuite("K-mer Analysis Performance")

# Test k-mer counting performance across different k values
for k in config["k_values"]
    println("\nBenchmarking k-mer counting with k=$k")
    
    # Test with medium-sized genome
    medium_genome_size = config["genome_sizes"][2]
    test_seq = test_sequences[medium_genome_size][1]
    
    # Benchmark k-mer counting
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.count_kmers, test_seq, k;
        samples=3, seconds=15
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "kmer_counting_k$(k)_size$(medium_genome_size)", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# ## Sparse K-mer Counting Benchmarks
#
# Test sparse k-mer counting performance

println("\n--- Sparse K-mer Counting Benchmarks ---")

# Benchmark canonical k-mer counting
for k in config["k_values"]
    println("\nBenchmarking canonical k-mer counting with k=$k")
    
    medium_genome_size = config["genome_sizes"][2]
    test_seq = test_sequences[medium_genome_size][1]
    
    # Benchmark canonical k-mer counting
    benchmark_result, memory_stats = run_benchmark_with_memory(
        Mycelia.count_canonical_kmers, test_seq, k;
        samples=3, seconds=15
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "canonical_kmer_counting_k$(k)_size$(medium_genome_size)", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# ## K-mer Spectrum Analysis Benchmarks
#
# Test k-mer spectrum analysis performance

println("\n--- K-mer Spectrum Analysis Benchmarks ---")

# Test k-mer frequency histogram generation
for k in [15, 21]  # Test with smaller k values for spectrum analysis
    println("\nBenchmarking k-mer frequency histogram with k=$k")
    
    medium_genome_size = config["genome_sizes"][2]
    test_seq = test_sequences[medium_genome_size][1]
    
    # Function to generate k-mer frequency histogram
    function generate_kmer_histogram(sequence, k)
        kmer_counts = Mycelia.count_kmers(sequence, k)
        frequencies = collect(values(kmer_counts))
        return Dict(freq => count(==(freq), frequencies) for freq in unique(frequencies))
    end
    
    benchmark_result, memory_stats = run_benchmark_with_memory(
        generate_kmer_histogram, test_seq, k;
        samples=3, seconds=15
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "kmer_histogram_k$(k)_size$(medium_genome_size)", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# Test genome size estimation accuracy
println("\nBenchmarking genome size estimation from k-mer spectrum")

k_estimate = 21
medium_genome_size = config["genome_sizes"][2]
test_seq = test_sequences[medium_genome_size][1]

# Function to estimate genome size from k-mer spectrum
function estimate_genome_size_from_kmers(sequence, k)
    kmer_counts = Mycelia.count_kmers(sequence, k)
    unique_kmers = length(kmer_counts)
    total_kmers = sum(values(kmer_counts))
    
    # Simple genome size estimation: total_kmers - k + 1 ≈ genome_size
    # This is a basic estimation; more sophisticated methods exist
    estimated_size = total_kmers - k + 1
    
    return Dict(
        "unique_kmers" => unique_kmers,
        "total_kmers" => total_kmers,
        "estimated_genome_size" => estimated_size,
        "actual_genome_size" => length(sequence)
    )
end

benchmark_result, memory_stats = run_benchmark_with_memory(
    estimate_genome_size_from_kmers, test_seq, k_estimate;
    samples=3, seconds=15
)

add_benchmark_result!(
    benchmark_suite, 
    "genome_size_estimation_k$(k_estimate)", 
    benchmark_result, 
    memory_stats
)

println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")

# Test the actual estimation for accuracy
estimation_result = estimate_genome_size_from_kmers(test_seq, k_estimate)
accuracy = 1.0 - abs(estimation_result["estimated_genome_size"] - estimation_result["actual_genome_size"]) / estimation_result["actual_genome_size"]
println("  Genome size estimation accuracy: $(round(accuracy * 100, digits=1))%")

# ## Scalability Benchmarks
#
# Test scalability with increasing data sizes

println("\n--- Scalability Benchmarks ---")

# Scaling benchmarks - test performance vs genome size
k_test = 21  # Use k=21 for scaling tests
println("\n--- Scaling Benchmarks (k=$k_test) ---")

# Function to generate test arguments for scaling
args_generator = function(genome_size)
    return (test_sequences[genome_size][1], k_test)
end

# Run scaling benchmark
scaling_results = run_scaling_benchmark(
    Mycelia.count_kmers,
    config["genome_sizes"],
    args_generator;
    samples=2,
    seconds=10
)

# Add scaling results to suite
for (size, result) in scaling_results
    add_benchmark_result!(
        benchmark_suite, 
        "scaling_kmer_size$(size)", 
        Dict("median_time" => result["median_time"], "memory" => result["memory"], "allocations" => result["allocations"]),
        result["memory_stats"]
    )
end

println("\nScaling Results:")
for (size, result) in scaling_results
    time_ms = result["median_time"] / 1e6
    memory_mb = result["memory"] / 1e6
    println("  Size $(size): $(round(time_ms, digits=2)) ms, $(round(memory_mb, digits=2)) MB")
end

# ## Accuracy Benchmarks
#
# Test accuracy of k-mer analysis algorithms

println("\n--- Accuracy Benchmarks ---")

# Memory efficiency analysis
println("\n--- Memory Efficiency Analysis ---")

for k in [15, 21, 31]  # Test key k values
    medium_genome_size = config["genome_sizes"][2]
    
    # Estimate theoretical memory requirement
    theoretical_kmers = 4^k  # Maximum possible k-mers
    actual_genome_size = medium_genome_size
    
    # Check if matrix would fit in memory
    memory_check = Mycelia.check_matrix_fits_in_memory(theoretical_kmers, 1, Float64)
    
    println("k=$k:")
    println("  Theoretical k-mer space: $(theoretical_kmers)")
    println("  Genome size: $(actual_genome_size)")
    println("  Dense matrix fits in memory: $memory_check")
    
    # Estimate sparse vs dense memory usage
    sparse_estimate = Mycelia.estimate_sparse_matrix_memory(actual_genome_size, 1, Float64)
    dense_estimate = Mycelia.estimate_dense_matrix_memory(theoretical_kmers, 1, Float64)
    
    println("  Sparse estimate: $(round(sparse_estimate / 1e6, digits=2)) MB")
    println("  Dense estimate: $(round(dense_estimate / 1e6, digits=2)) MB")
    println("  Memory ratio (sparse/dense): $(round(sparse_estimate / dense_estimate, digits=4))")
end

# Save benchmark results
results_dir = "results"
mkpath(results_dir)
results_file = joinpath(results_dir, "kmer_analysis_benchmark_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")
save_benchmark_results(benchmark_suite, results_file)

# Display summary
format_benchmark_summary(benchmark_suite)

println("\n=== K-mer Analysis Benchmark Complete ===")
println("Results saved to: $results_file")
println("End time: $(Dates.now())")

# ## Memory Efficiency Benchmarks
#
# Test memory usage optimization

println("\n--- Memory Efficiency Benchmarks ---")

# Test peak memory usage during k-mer counting
println("\nBenchmarking peak memory usage during k-mer operations")

k_memory_test = 21
large_genome_size = config["genome_sizes"][end]  # Use largest genome size

# Test memory-intensive operations
if haskey(test_sequences, large_genome_size)
    test_seq = test_sequences[large_genome_size][1]
    
    # Monitor memory during large k-mer counting
    function memory_intensive_kmer_count(sequence, k)
        GC.gc()  # Start with clean memory
        initial_memory = Base.gc_num()
        
        result = Mycelia.count_kmers(sequence, k)
        
        GC.gc()
        final_memory = Base.gc_num()
        
        return result, (final_memory.allocd - initial_memory.allocd)
    end
    
    benchmark_result, memory_stats = run_benchmark_with_memory(
        memory_intensive_kmer_count, test_seq, k_memory_test;
        samples=2, seconds=30
    )
    
    add_benchmark_result!(
        benchmark_suite, 
        "memory_intensive_kmer_k$(k_memory_test)_size$(large_genome_size)", 
        benchmark_result, 
        memory_stats
    )
    
    println("  Large genome k-mer counting:")
    println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
    println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")
end

# Test cache efficiency with repeated k-mer operations
println("\nBenchmarking cache efficiency with repeated operations")

medium_genome_size = config["genome_sizes"][2]
test_seq = test_sequences[medium_genome_size][1]

# Function to test cache performance
function repeated_kmer_operations(sequence, k, iterations=3)
    results = []
    for i in 1:iterations
        result = Mycelia.count_canonical_kmers(sequence, k)
        push!(results, length(result))
    end
    return results
end

benchmark_result, memory_stats = run_benchmark_with_memory(
    repeated_kmer_operations, test_seq, 15, 3;
    samples=3, seconds=20
)

add_benchmark_result!(
    benchmark_suite, 
    "cache_efficiency_repeated_ops", 
    benchmark_result, 
    memory_stats
)

println("  Cache efficiency test:")
println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")

# Test memory fragmentation patterns
println("\nAnalyzing memory fragmentation patterns")

# Function to create and destroy k-mer dictionaries rapidly
function memory_fragmentation_test(sequence, k_values)
    results = Dict()
    for k in k_values
        GC.gc()  # Force garbage collection
        result = Mycelia.count_kmers(sequence, k)
        results[k] = length(result)
        # Allow dictionary to be garbage collected
        result = nothing
    end
    GC.gc()
    return results
end

k_frag_values = [11, 15, 19]
benchmark_result, memory_stats = run_benchmark_with_memory(
    memory_fragmentation_test, test_seq, k_frag_values;
    samples=2, seconds=25
)

add_benchmark_result!(
    benchmark_suite, 
    "memory_fragmentation_test", 
    benchmark_result, 
    memory_stats
)

println("  Memory fragmentation test:")
println("  Median time: $(round(BenchmarkTools.median(benchmark_result).time / 1e6, digits=2)) ms")
println("  Memory: $(round(BenchmarkTools.median(benchmark_result).memory / 1e6, digits=2)) MB")

# ## Results Collection
#
# Collect comprehensive benchmark results

println("\n--- Collecting Results ---")

# Calculate k-mer analysis performance metrics
total_sequences_tested = sum(length(seqs) for seqs in values(test_sequences))
total_sequence_length = sum(sum(length(seq) for seq in seqs) for seqs in values(test_sequences))

# Calculate k-mer counting throughput
k_throughput_estimates = Dict()
for k in config["k_values"][1:3]  # Test first 3 k values
    test_key = "kmer_counting_k$(k)_size$(config["genome_sizes"][2])"
    if haskey(benchmark_suite.results, test_key)
        time_seconds = benchmark_suite.results[test_key]["median_time"] / 1e9
        genome_size = config["genome_sizes"][2]
        kmers_per_second = (genome_size - k + 1) / time_seconds
        k_throughput_estimates[k] = round(kmers_per_second, digits=0)
    end
end

# Analyze scaling performance
scaling_efficiency = Dict()
if length(config["genome_sizes"]) >= 2
    for i in 2:length(config["genome_sizes"])
        size_ratio = config["genome_sizes"][i] / config["genome_sizes"][i-1]
        
        # Find corresponding benchmark results
        prev_key = "scaling_kmer_size$(config["genome_sizes"][i-1])"
        curr_key = "scaling_kmer_size$(config["genome_sizes"][i])"
        
        if haskey(benchmark_suite.results, prev_key) && haskey(benchmark_suite.results, curr_key)
            prev_time = get(benchmark_suite.results[prev_key], "median_time", 1.0)
            curr_time = get(benchmark_suite.results[curr_key], "median_time", 1.0)
            time_ratio = curr_time / prev_time
            efficiency = size_ratio / time_ratio  # Ideal is 1.0 (linear scaling)
            scaling_efficiency[config["genome_sizes"][i]] = round(efficiency, digits=3)
        end
    end
end

# Comprehensive results structure
comprehensive_results = Dict(
    "benchmark_name" => "kmer_analysis",
    "timestamp" => Dates.now(),
    "configuration" => config,
    "system_info" => Dict(
        "julia_version" => VERSION,
        "cpu_threads" => Sys.CPU_THREADS,
        "total_memory" => Sys.total_memory(),
        "hostname" => gethostname(),
        "platform" => Sys.MACHINE
    ),
    "dataset_summary" => Dict(
        "total_sequences_tested" => total_sequences_tested,
        "total_sequence_length" => total_sequence_length,
        "genome_sizes_tested" => config["genome_sizes"],
        "k_values_tested" => config["k_values"]
    ),
    "performance_metrics" => Dict(
        "k_mer_throughput_per_second" => k_throughput_estimates,
        "scaling_efficiency" => scaling_efficiency,
        "benchmark_duration_minutes" => round((Dates.now() - benchmark_suite.metadata["timestamp"]).value / 60000, digits=2)
    ),
    "resource_utilization" => Dict(
        "peak_memory_mb" => round(maximum([get(result, "memory", 0) for result in values(benchmark_suite.results)]) / 1e6, digits=2),
        "total_allocations" => sum([get(result, "allocations", 0) for result in values(benchmark_suite.results)]),
        "memory_efficiency_score" => round(total_sequence_length / maximum([get(result, "memory", 1) for result in values(benchmark_suite.results)]), digits=2)
    ),
    "detailed_results" => benchmark_suite.results
)

# Save comprehensive results 
comprehensive_results_file = joinpath(results_dir, "kmer_analysis_comprehensive_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")

open(comprehensive_results_file, "w") do f
    JSON.print(f, comprehensive_results, 2)
end

println("\n--- K-mer Analysis Performance Summary ---")
println("Sequences tested: $total_sequences_tested")
println("Total sequence length: $(round(total_sequence_length / 1e6, digits=2)) Mbp")
println("K-mer throughput estimates:")
for (k, throughput) in k_throughput_estimates
    println("  k=$k: $(Int(throughput)) k-mers/second")
end
println("Scaling efficiency:")
for (size, efficiency) in scaling_efficiency
    println("  $(size)bp: $(efficiency) (1.0 = perfect linear scaling)")
end
println("Peak memory usage: $(round(maximum([get(result, "memory", 0) for result in values(benchmark_suite.results)]) / 1e6, digits=2)) MB")

println("\nComprehensive results saved to: $comprehensive_results_file")
println("End time: $(Dates.now())")

nothing