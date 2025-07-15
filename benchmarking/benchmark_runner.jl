#!/usr/bin/env julia

# Benchmark Runner with Performance Regression Tracking
#
# This script runs the complete Mycelia benchmark suite and checks for
# performance regressions against baseline results.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Dates
import JSON

include("benchmark_utils.jl")

"""
    run_benchmark_suite(scale="small"; check_regression=true, baseline_dir="baselines")

Run the complete Mycelia benchmark suite with performance regression checking.

# Arguments
- `scale`: Benchmark scale ("small", "medium", "large")
- `check_regression`: Whether to check for performance regressions
- `baseline_dir`: Directory containing baseline benchmark results
"""
function run_benchmark_suite(scale="small"; check_regression=true, baseline_dir="baselines")
    println("="^60)
    println("MYCELIA PERFORMANCE BENCHMARK SUITE")
    println("="^60)
    println("Scale: $scale")
    println("Start time: $(Dates.now())")
    println("Regression checking: $check_regression")
    
    # Set environment variable for benchmark scale
    ENV["BENCHMARK_SCALE"] = scale
    
    # Create results directory
    results_dir = "results"
    mkpath(results_dir)
    mkpath(baseline_dir)
    
    # Benchmark scripts to run
    benchmark_scripts = [
        ("01_data_processing_benchmark.jl", "Data Processing"),
        ("02_kmer_analysis_benchmark.jl", "K-mer Analysis"),
        ("03_assembly_benchmark.jl", "Assembly"),
        ("04_annotation_benchmark.jl", "Annotation"),
        # Add other benchmarks as they are implemented
        # ("05_comparative_benchmark.jl", "Comparative Genomics")
    ]
    
    # Run each benchmark
    benchmark_results = Dict{String, Any}()\n    
    for (script, description) in benchmark_scripts
        println("\\n" * "-"^40)
        println("Running: $description")
        println("-"^40)
        
        try
            # Run benchmark script
            script_path = joinpath(@__DIR__, script)
            if isfile(script_path)
                include(script_path)
                println("✅ $description benchmark completed successfully")
            else
                println("⚠️  Benchmark script not found: $script")
                continue
            end
            
            # Find the most recent results file for this benchmark
            pattern = replace(script, ".jl" => "_")
            result_files = filter(f -> startswith(f, pattern), readdir(results_dir))
            
            if !isempty(result_files)
                latest_result = joinpath(results_dir, sort(result_files)[end])
                benchmark_results[description] = latest_result
            end
            
        catch e
            println("❌ Error running $description benchmark:")
            println("   $(typeof(e)): $e")
            continue
        end
    end
    
    # Performance regression checking
    if check_regression
        println("\\n" * "="^40)
        println("PERFORMANCE REGRESSION ANALYSIS")
        println("="^40)
        
        for (description, results_file) in benchmark_results
            baseline_file = joinpath(baseline_dir, "$(replace(description, " " => "_"))_baseline.json")
            
            println("\\nChecking regressions for: $description")
            check_performance_regression(results_file, baseline_file, threshold=0.10)
            
            # Update baseline if this is the first run or if explicitly requested
            if !isfile(baseline_file) || get(ENV, "UPDATE_BASELINES", "false") == "true"
                cp(results_file, baseline_file, force=true)
                println("📝 Updated baseline: $baseline_file")
            end
        end
    end
    
    # Generate summary report
    generate_benchmark_report(benchmark_results, scale)
    
    println("\\n" * "="^60)
    println("BENCHMARK SUITE COMPLETE")
    println("="^60)
    println("End time: $(Dates.now())")
    println("Results directory: $results_dir")
    println("Baselines directory: $baseline_dir")
end

"""
    generate_benchmark_report(benchmark_results, scale)

Generate a comprehensive HTML report of benchmark results.
"""
function generate_benchmark_report(benchmark_results, scale)
    println("\\n--- Generating Benchmark Report ---")
    
    # Load all benchmark results
    all_results = Dict{String, Any}()
    for (description, results_file) in benchmark_results
        try
            all_results[description] = load_benchmark_results(results_file)
        catch e
            println("⚠️  Could not load results for $description: $e")
        end
    end
    
    # Generate HTML report
    html_content = generate_html_report(all_results, scale)
    
    # Save report
    report_file = "results/benchmark_report_$(scale)_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).html"
    open(report_file, "w") do f
        write(f, html_content)
    end
    
    println("📊 Benchmark report generated: $report_file")
end

"""
    generate_html_report(results, scale)

Generate HTML content for the benchmark report.
"""
function generate_html_report(results, scale)
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Mycelia Performance Benchmark Report - $scale Scale</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
            .benchmark-section { margin: 20px 0; border: 1px solid #ddd; border-radius: 5px; }
            .benchmark-header { background-color: #e8f4f8; padding: 15px; font-weight: bold; }
            .benchmark-content { padding: 15px; }
            table { width: 100%; border-collapse: collapse; margin: 10px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            .metric { margin: 10px 0; }
            .time { color: #2e8b57; }
            .memory { color: #4169e1; }
            .regression { color: #dc143c; font-weight: bold; }
            .improvement { color: #228b22; font-weight: bold; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Mycelia Performance Benchmark Report</h1>
            <p><strong>Scale:</strong> $scale</p>
            <p><strong>Generated:</strong> $(Dates.now())</p>
            <p><strong>Julia Version:</strong> $(VERSION)</p>
        </div>
    """
    
    for (benchmark_name, benchmark_data) in results
        html *= """
        <div class="benchmark-section">
            <div class="benchmark-header">$benchmark_name Benchmark</div>
            <div class="benchmark-content">
                <p><strong>Timestamp:</strong> $(get(benchmark_data, "timestamp", "N/A"))</p>
                <table>
                    <tr>
                        <th>Test Name</th>
                        <th>Median Time (ms)</th>
                        <th>Memory (MB)</th>
                        <th>Allocations</th>
                    </tr>
        """
        
        if haskey(benchmark_data, "results")
            for (test_name, test_result) in benchmark_data["results"]
                median_time_ms = get(test_result, "median_time", 0) / 1e6
                memory_mb = get(test_result, "memory", 0) / 1e6
                allocations = get(test_result, "allocations", 0)
                
                html *= """
                    <tr>
                        <td>$test_name</td>
                        <td class="time">$(round(median_time_ms, digits=2))</td>
                        <td class="memory">$(round(memory_mb, digits=2))</td>
                        <td>$allocations</td>
                    </tr>
                """
            end
        end
        
        html *= """
                </table>
            </div>
        </div>
        """
    end
    
    html *= """
    </body>
    </html>
    """
    
    return html
end

"""
    main()

Main entry point for the benchmark runner.
"""
function main()
    # Parse command line arguments
    scale = get(ARGS, 1, "small")
    check_regression = get(ENV, "CHECK_REGRESSION", "true") == "true"
    
    # Validate scale argument
    if !(scale in ["small", "medium", "large"])
        println("❌ Invalid scale: $scale")
        println("   Valid options: small, medium, large")
        exit(1)
    end
    
    # Run benchmark suite
    try
        run_benchmark_suite(scale, check_regression=check_regression)
        println("\\n✅ Benchmark suite completed successfully!")
    catch e
        println("\\n❌ Benchmark suite failed:")
        println("   $(typeof(e)): $e")
        rethrow(e)
    end
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end