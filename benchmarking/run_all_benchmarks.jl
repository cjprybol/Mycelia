#!/usr/bin/env julia
# Benchmarking Suite Runner for Mycelia
# This runs performance benchmarks - suitable for HPC environments

println("""
================================================================================
Mycelia Benchmarking Suite
================================================================================
This will run performance benchmarks which may take significant time and resources.
Ensure you have:
- Sufficient memory (>32GB recommended)
- Multiple CPU cores for parallel tests
- Long timeout settings
- Time (this may take hours)
================================================================================
""")

# Activate the project
import Pkg
Pkg.activate(dirname(@__DIR__))

# Import required packages
import Mycelia
import BenchmarkTools
import Dates

println("\nStarting benchmarks at $(Dates.now())")
start_time = time()

try
    # Run benchmarks in order
    println("\n=== Data Processing Benchmark ===")
    include("01_data_processing_benchmark.jl")
    
    println("\n=== K-mer Analysis Benchmark ===")
    include("02_kmer_analysis_benchmark.jl")
    
    println("\n=== Assembly Benchmark ===")
    include("03_assembly_benchmark.jl")
    
    println("\n=== Annotation Benchmark ===")
    include("04_annotation_benchmark.jl")
    
    println("\n=== Comparative Analysis Benchmark ===")
    include("05_comparative_benchmark.jl")
    
    println("\n=== CoverM Coverage Benchmark ===")
    include("06_coverm_benchmark.jl")
    
    elapsed = round(time() - start_time, digits=1)
    println("\n✅ All benchmarks completed in $elapsed seconds!")
    
catch e
    elapsed = round(time() - start_time, digits=1)
    println("\n❌ Benchmarks failed after $elapsed seconds")
    rethrow(e)
end
