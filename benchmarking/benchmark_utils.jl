# Benchmarking Utilities for Mycelia
#
# This module provides common utilities for performance benchmarking
# including memory profiling, result serialization, and performance tracking.

import BenchmarkTools
import Statistics
import Dates
import JSON
import JLD2

"""
    BenchmarkSuite

A comprehensive benchmarking suite that tracks performance metrics,
memory usage, and provides regression analysis capabilities.
"""
struct BenchmarkSuite
    name::String
    results::Dict{String, Any}
    metadata::Dict{String, Any}
    
    function BenchmarkSuite(name::String)
        new(name, Dict{String, Any}(), Dict{String, Any}(
            "timestamp" => Dates.now(),
            "julia_version" => VERSION,
            "cpu_info" => Sys.cpu_info()[1],
            "memory_total" => Sys.total_memory()
        ))
    end
end

"""
    run_benchmark_with_memory(func, args...; kwargs...)

Run a function with BenchmarkTools while tracking memory allocation.
Returns both benchmark results and memory statistics.
"""
function run_benchmark_with_memory(func, args...; samples=5, seconds=30, kwargs...)
    # Get initial memory state
    GC.gc()
    initial_memory = Base.gc_num()
    
    # Run benchmark
    benchmark_result = BenchmarkTools.@benchmark func(args...; kwargs...) samples=samples seconds=seconds
    
    # Get final memory state
    GC.gc()
    final_memory = Base.gc_num()
    
    # Calculate memory statistics
    memory_stats = Dict(
        "allocations" => final_memory.allocd - initial_memory.allocd,
        "gc_time" => (final_memory.total_time - initial_memory.total_time) / 1e9,
        "gc_count" => final_memory.pause - initial_memory.pause
    )
    
    return benchmark_result, memory_stats
end

"""
    estimate_memory_requirement(func, args...; scale_factor=1.0)

Estimate memory requirements for a function call based on argument sizes
and apply an optional scale factor for safety margins.
"""
function estimate_memory_requirement(func, args...; scale_factor=1.0)
    total_size = 0
    
    for arg in args
        if isa(arg, AbstractArray)
            total_size += sizeof(arg)
        elseif isa(arg, AbstractString)
            total_size += sizeof(arg)
        elseif isa(arg, Dict)
            total_size += Base.summarysize(arg)
        else
            total_size += Base.summarysize(arg)
        end
    end
    
    return total_size * scale_factor
end

"""
    profile_memory_usage(func, args...; kwargs...)

Profile memory usage patterns during function execution.
Returns detailed memory allocation information.
"""
function profile_memory_usage(func, args...; kwargs...)
    # Clear previous allocations
    GC.gc()
    
    # Enable allocation tracking
    Profile.clear_malloc_data()
    
    # Run function with allocation tracking
    result = @time func(args...; kwargs...)
    
    # Get allocation profile
    allocation_data = Profile.fetch()
    
    return result, allocation_data
end

"""
    save_benchmark_results(suite::BenchmarkSuite, filepath::String)

Save benchmark results to a JSON file with metadata.
"""
function save_benchmark_results(suite::BenchmarkSuite, filepath::String)
    output = Dict(
        "suite_name" => suite.name,
        "metadata" => suite.metadata,
        "results" => suite.results
    )
    
    open(filepath, "w") do f
        JSON.print(f, output, 2)
    end
    
    println("Benchmark results saved to: $filepath")
end

"""
    load_benchmark_results(filepath::String)

Load benchmark results from a JSON file.
"""
function load_benchmark_results(filepath::String)
    return JSON.parsefile(filepath)
end

"""
    compare_benchmark_results(current_results, reference_results; threshold=0.05)

Compare current benchmark results against reference results.
Returns performance regression warnings if performance degrades beyond threshold.
"""
function compare_benchmark_results(current_results, reference_results; threshold=0.05)
    regressions = []
    improvements = []
    
    for (test_name, current_result) in current_results["results"]
        if haskey(reference_results["results"], test_name)
            ref_result = reference_results["results"][test_name]
            
            if haskey(current_result, "median_time") && haskey(ref_result, "median_time")
                current_time = current_result["median_time"]
                ref_time = ref_result["median_time"]
                
                change_ratio = (current_time - ref_time) / ref_time
                
                if change_ratio > threshold
                    push!(regressions, (test_name, change_ratio))
                elseif change_ratio < -threshold
                    push!(improvements, (test_name, abs(change_ratio)))
                end
            end
        end
    end
    
    return regressions, improvements
end

"""
    add_benchmark_result!(suite::BenchmarkSuite, name::String, benchmark_result, memory_stats=nothing)

Add a benchmark result to the suite with optional memory statistics.
"""
function add_benchmark_result!(suite::BenchmarkSuite, name::String, benchmark_result, memory_stats=nothing)
    result_data = Dict(
        "median_time" => BenchmarkTools.median(benchmark_result).time,
        "mean_time" => BenchmarkTools.mean(benchmark_result).time,
        "min_time" => BenchmarkTools.minimum(benchmark_result).time,
        "max_time" => BenchmarkTools.maximum(benchmark_result).time,
        "std_time" => BenchmarkTools.std(benchmark_result).time,
        "allocations" => BenchmarkTools.median(benchmark_result).allocs,
        "memory" => BenchmarkTools.median(benchmark_result).memory
    )
    
    if memory_stats !== nothing
        result_data["memory_stats"] = memory_stats
    end
    
    suite.results[name] = result_data
end

"""
    format_benchmark_summary(suite::BenchmarkSuite)

Format benchmark results into a human-readable summary.
"""
function format_benchmark_summary(suite::BenchmarkSuite)
    println("\n" * "="^60)
    println("BENCHMARK SUMMARY: $(suite.name)")
    println("="^60)
    println("Timestamp: $(suite.metadata["timestamp"])")
    println("Julia Version: $(suite.metadata["julia_version"])")
    println("Total Memory: $(round(suite.metadata["memory_total"] / 1e9, digits=2)) GB")
    println("-"^60)
    
    for (test_name, result) in suite.results
        median_time_ms = result["median_time"] / 1e6  # Convert to milliseconds
        memory_mb = result["memory"] / 1e6  # Convert to megabytes
        
        println("$test_name:")
        println("  Time (median): $(round(median_time_ms, digits=3)) ms")
        println("  Memory: $(round(memory_mb, digits=2)) MB")
        println("  Allocations: $(result["allocations"])")
        println()
    end
    
    println("="^60)
end

"""
    run_scaling_benchmark(func, input_sizes, args_generator; kwargs...)

Run benchmarks across different input sizes to analyze scaling behavior.
"""
function run_scaling_benchmark(func, input_sizes, args_generator; samples=3, seconds=10, kwargs...)
    results = Dict()
    
    for size in input_sizes
        println("Benchmarking with input size: $size")
        
        # Generate arguments for this size
        args = args_generator(size)
        
        # Run benchmark
        benchmark_result, memory_stats = run_benchmark_with_memory(
            func, args...; 
            samples=samples, 
            seconds=seconds, 
            kwargs...
        )
        
        results[size] = Dict(
            "median_time" => BenchmarkTools.median(benchmark_result).time,
            "memory" => BenchmarkTools.median(benchmark_result).memory,
            "allocations" => BenchmarkTools.median(benchmark_result).allocs,
            "memory_stats" => memory_stats
        )
        
        # Force garbage collection between runs
        GC.gc()
    end
    
    return results
end

"""
    check_performance_regression(current_file::String, reference_file::String; threshold=0.10)

Check for performance regressions by comparing current results with reference.
Prints warnings for any performance degradations beyond the threshold.
"""
function check_performance_regression(current_file::String, reference_file::String; threshold=0.10)
    if !isfile(reference_file)
        println("⚠️  No reference benchmark file found at: $reference_file")
        println("   Current results will serve as the new baseline.")
        return
    end
    
    current_results = load_benchmark_results(current_file)
    reference_results = load_benchmark_results(reference_file)
    
    regressions, improvements = compare_benchmark_results(
        current_results, reference_results, threshold=threshold
    )
    
    if !isempty(regressions)
        println("\n🚨 PERFORMANCE REGRESSIONS DETECTED:")
        for (test_name, change_ratio) in regressions
            println("  $test_name: $(round(change_ratio * 100, digits=1))% slower")
        end
    end
    
    if !isempty(improvements)
        println("\n✅ PERFORMANCE IMPROVEMENTS:")
        for (test_name, improvement_ratio) in improvements
            println("  $test_name: $(round(improvement_ratio * 100, digits=1))% faster")
        end
    end
    
    if isempty(regressions) && isempty(improvements)
        println("\n✅ No significant performance changes detected.")
    end
end