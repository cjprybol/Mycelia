# # Mycelia Extended Tests and Tutorials Runner
#
# This script provides a single command interface for running extended tests,
# tutorials, and benchmarks. It supports both local execution and HPC submission.
#
# ## Usage
#
# ```bash
# # Run all tutorials locally
# julia --project=. run_extended_tests.jl tutorials
#
# # Run all benchmarks locally (use with caution)
# julia --project=. run_extended_tests.jl benchmarks
#
# # Submit benchmarks to SLURM HPC
# julia --project=. run_extended_tests.jl benchmarks --hpc
#
# # Convert tutorials to notebooks
# julia --project=. run_extended_tests.jl notebooks
#
# # Run specific tutorial
# julia --project=. run_extended_tests.jl tutorial 01_data_acquisition
#
# # Show help
# julia --project=. run_extended_tests.jl help
# ```

import Pkg
Pkg.activate(".")

using Dates
import Literate

# Configuration
const TUTORIALS_DIR = "tutorials"
const BENCHMARKS_DIR = "benchmarking"
const NOTEBOOKS_DIR = "tutorials/notebooks"
const RESULTS_DIR = "results"

"""
    show_help()

Display usage information and available commands.
"""
function show_help()
    println("""
    Mycelia Extended Tests and Tutorials Runner
    
    Usage: julia --project=. run_extended_tests.jl <command> [options]
    
    Commands:
        tutorials               Run all tutorials locally
        benchmarks              Run all benchmarks locally (resource intensive!)
        benchmarks --hpc        Submit benchmarks to SLURM HPC
        notebooks               Convert all tutorials to Jupyter notebooks
        reports                 Generate HTML reports from executed tutorials/benchmarks
        tutorial <name>         Run specific tutorial (e.g., 01_data_acquisition)
        help                    Show this help message
    
    Examples:
        julia --project=. run_extended_tests.jl tutorials
        julia --project=. run_extended_tests.jl benchmarks --hpc
        julia --project=. run_extended_tests.jl tutorial 01_data_acquisition
        julia --project=. run_extended_tests.jl notebooks
        julia --project=. run_extended_tests.jl reports
    
    Notes:
        - Tutorials use small datasets and run quickly
        - Benchmarks use large datasets and require significant resources
        - Use --hpc flag for benchmark submission to SLURM clusters
        - Reports generate HTML with executed code, outputs, and visualizations
        - Extended tests are NOT run in CI/CD - they're for manual execution
    """)
end

"""
    run_tutorials()

Execute all tutorials in sequence.
"""
function run_tutorials()
    println("=== Running All Mycelia Tutorials ===")
    println("Start time: $(now())")
    println()
    
    if !isdir(TUTORIALS_DIR)
        println("Error: Tutorials directory not found: $TUTORIALS_DIR")
        return false
    end
    
    # Find all tutorial files
    tutorial_files = filter(f -> endswith(f, ".jl") && f != "run_all_tutorials.jl", 
                           readdir(TUTORIALS_DIR))
    sort!(tutorial_files)
    
    if isempty(tutorial_files)
        println("No tutorial files found in $TUTORIALS_DIR")
        return false
    end
    
    println("Found $(length(tutorial_files)) tutorials:")
    for file in tutorial_files
        println("  - $file")
    end
    println()
    
    # Run tutorials
    results = []
    total_start = time()
    
    for tutorial_file in tutorial_files
        tutorial_path = joinpath(TUTORIALS_DIR, tutorial_file)
        println("Running $tutorial_file...")
        
        start_time = time()
        try
            include(tutorial_path)
            elapsed = time() - start_time
            push!(results, (tutorial_file, :success, elapsed))
            println("✓ $tutorial_file completed in $(round(elapsed, digits=2))s")
        catch e
            elapsed = time() - start_time
            push!(results, (tutorial_file, :error, elapsed, e))
            println("✗ $tutorial_file failed after $(round(elapsed, digits=2))s")
            println("  Error: $e")
        end
        println()
    end
    
    total_elapsed = time() - total_start
    
    # Summary
    println("=== Tutorial Summary ===")
    println("Total time: $(round(total_elapsed, digits=2))s")
    
    successful = count(r -> r[2] == :success, results)
    println("Successful: $successful/$(length(results))")
    println("End time: $(now())")
    
    return successful == length(results)
end

"""
    run_benchmarks(use_hpc=false)

Execute benchmarks either locally or submit to HPC.
"""
function run_benchmarks(use_hpc=false)
    if use_hpc
        println("=== Submitting Benchmarks to SLURM HPC ===")
        return submit_benchmarks_to_hpc()
    else
        println("=== Running Benchmarks Locally ===")
        println("⚠ WARNING: Benchmarks are resource-intensive!")
        println("   - May consume significant CPU, memory, and disk space")
        println("   - May take hours or days to complete")
        println("   - Consider using --hpc flag for SLURM submission")
        println()
        
        print("Continue with local execution? (y/N): ")
        response = readline()
        if lowercase(strip(response)) != "y"
            println("Benchmark execution cancelled")
            return false
        end
        
        return run_benchmarks_locally()
    end
end

"""
    submit_benchmarks_to_hpc()

Submit benchmark job to SLURM HPC system.
"""
function submit_benchmarks_to_hpc()
    benchmark_script = joinpath(BENCHMARKS_DIR, "run_all_benchmarks.sh")
    
    if !isfile(benchmark_script)
        println("Error: Benchmark script not found: $benchmark_script")
        return false
    end
    
    try
        # Submit job to SLURM
        cmd = `sbatch $benchmark_script`
        output = read(cmd, String)
        
        println("✓ Benchmark job submitted successfully")
        println(output)
        
        # Extract job ID if possible
        if contains(output, "Submitted batch job")
            job_id = match(r"Submitted batch job (\d+)", output)
            if job_id !== nothing
                println("Job ID: $(job_id.captures[1])")
                println("Monitor with: squeue -j $(job_id.captures[1])")
                println("View output: tail -f benchmarks_$(job_id.captures[1]).out")
            end
        end
        
        return true
    catch e
        println("Error submitting benchmark job: $e")
        return false
    end
end

"""
    run_benchmarks_locally()

Run benchmarks on local machine.
"""
function run_benchmarks_locally()
    println("Starting local benchmark execution...")
    
    if !isdir(BENCHMARKS_DIR)
        println("Error: Benchmarks directory not found: $BENCHMARKS_DIR")
        return false
    end
    
    # Create results directory
    if !isdir(RESULTS_DIR)
        mkdir(RESULTS_DIR)
    end
    
    # Find benchmark files
    benchmark_files = filter(f -> endswith(f, ".jl"), readdir(BENCHMARKS_DIR))
    sort!(benchmark_files)
    
    if isempty(benchmark_files)
        println("No benchmark files found in $BENCHMARKS_DIR")
        return false
    end
    
    println("Found $(length(benchmark_files)) benchmarks:")
    for file in benchmark_files
        println("  - $file")
    end
    println()
    
    # Run benchmarks
    results = []
    total_start = time()
    
    for benchmark_file in benchmark_files
        benchmark_path = joinpath(BENCHMARKS_DIR, benchmark_file)
        println("Running $benchmark_file...")
        
        start_time = time()
        try
            include(benchmark_path)
            elapsed = time() - start_time
            push!(results, (benchmark_file, :success, elapsed))
            println("✓ $benchmark_file completed in $(round(elapsed, digits=2))s")
        catch e
            elapsed = time() - start_time
            push!(results, (benchmark_file, :error, elapsed, e))
            println("✗ $benchmark_file failed after $(round(elapsed, digits=2))s")
            println("  Error: $e")
        end
        println()
    end
    
    total_elapsed = time() - total_start
    
    # Summary
    println("=== Benchmark Summary ===")
    println("Total time: $(round(total_elapsed, digits=2))s")
    
    successful = count(r -> r[2] == :success, results)
    println("Successful: $successful/$(length(results))")
    println("Results saved to: $RESULTS_DIR")
    println("End time: $(now())")
    
    return successful == length(results)
end

"""
    convert_to_notebooks()

Convert all tutorials to Jupyter notebooks.
"""
function convert_to_notebooks()
    println("=== Converting Tutorials to Notebooks ===")
    
    if !isdir(TUTORIALS_DIR)
        println("Error: Tutorials directory not found: $TUTORIALS_DIR")
        return false
    end
    
    # Create notebooks directory
    if !isdir(NOTEBOOKS_DIR)
        mkpath(NOTEBOOKS_DIR)
    end
    
    # Find tutorial files
    tutorial_files = filter(f -> endswith(f, ".jl") && f != "run_all_tutorials.jl", 
                           readdir(TUTORIALS_DIR))
    sort!(tutorial_files)
    
    if isempty(tutorial_files)
        println("No tutorial files found in $TUTORIALS_DIR")
        return false
    end
    
    println("Converting $(length(tutorial_files)) tutorials:")
    
    successful = 0
    for tutorial_file in tutorial_files
        tutorial_path = joinpath(TUTORIALS_DIR, tutorial_file)
        
        try
            println("  Converting $tutorial_file...")
            Literate.notebook(tutorial_path, NOTEBOOKS_DIR, execute=false)
            println("  ✓ $tutorial_file converted")
            successful += 1
        catch e
            println("  ✗ Failed to convert $tutorial_file: $e")
        end
    end
    
    println()
    println("Conversion complete: $successful/$(length(tutorial_files)) successful")
    println("Notebooks saved to: $NOTEBOOKS_DIR")
    
    return successful == length(tutorial_files)
end

"""
    generate_html_reports()

Generate HTML reports from tutorials and benchmarks by executing them as notebooks
and converting to HTML using nbconvert.
"""
function generate_html_reports()
    println("=== Generating HTML Reports ===")
    
    # Create output directories
    reports_dir = "reports"
    notebooks_dir = joinpath(reports_dir, "notebooks")
    html_dir = joinpath(reports_dir, "html")
    
    for dir in [reports_dir, notebooks_dir, html_dir]
        if !isdir(dir)
            mkpath(dir)
        end
    end
    
    # Check if jupyter/nbconvert is available
    if !success(`which jupyter`)
        println("Warning: jupyter not found. Install with: pip install jupyter nbconvert")
        println("Continuing with notebook generation only...")
        return convert_to_notebooks()
    end
    
    # Find all tutorial and benchmark files
    all_files = String[]
    
    # Add tutorials
    if isdir(TUTORIALS_DIR)
        for file in readdir(TUTORIALS_DIR)
            if endswith(file, ".jl") && !startswith(file, "run_") && !startswith(file, "generate_")
                push!(all_files, joinpath(TUTORIALS_DIR, file))
            end
        end
    end
    
    # Add benchmarks
    if isdir(BENCHMARKS_DIR)
        for file in readdir(BENCHMARKS_DIR)
            if endswith(file, ".jl") && !startswith(file, "run_") && !startswith(file, "generate_")
                push!(all_files, joinpath(BENCHMARKS_DIR, file))
            end
        end
    end
    
    if isempty(all_files)
        println("No tutorial or benchmark files found")
        return false
    end
    
    println("Found $(length(all_files)) files to process")
    
    successful = 0
    total_start = time()
    
    for file_path in all_files
        file_name = basename(file_path)
        println("Processing $file_name...")
        
        try
            # Step 1: Convert to executed notebook
            println("  Converting to executed notebook...")
            start_time = time()
            
            Literate.notebook(file_path, notebooks_dir, execute=true)
            
            notebook_file = joinpath(notebooks_dir, replace(file_name, ".jl" => ".ipynb"))
            
            if !isfile(notebook_file)
                println("  ✗ Failed to create notebook: $notebook_file")
                continue
            end
            
            elapsed = time() - start_time
            println("  ✓ Notebook created in $(round(elapsed, digits=1))s")
            
            # Step 2: Convert notebook to HTML
            println("  Converting notebook to HTML...")
            start_time = time()
            
            html_file = joinpath(html_dir, replace(file_name, ".jl" => ".html"))
            
            # Use nbconvert to create HTML
            cmd = `jupyter nbconvert --to html $notebook_file --output-dir $html_dir`
            
            if success(cmd)
                elapsed = time() - start_time
                println("  ✓ HTML report created in $(round(elapsed, digits=1))s")
                successful += 1
            else
                println("  ✗ Failed to convert notebook to HTML")
            end
            
        catch e
            println("  ✗ Error processing $file_name: $e")
        end
        
        println()
    end
    
    total_elapsed = time() - total_start
    
    # Generate index page
    println("Creating index page...")
    create_reports_index(html_dir, all_files)
    
    # Summary
    println("=== HTML Report Generation Summary ===")
    println("Total time: $(round(total_elapsed, digits=1))s")
    println("Successful: $successful/$(length(all_files))")
    println("Reports saved to: $html_dir")
    
    if successful > 0
        println("\nTo view reports:")
        println("  Open: $html_dir/index.html")
        println("  Or serve with: python -m http.server 8000 --directory $html_dir")
    end
    
    return successful == length(all_files)
end

"""
    create_reports_index(html_dir, source_files)

Create an index.html page linking to all generated reports.
"""
function create_reports_index(html_dir, source_files)
    # Separate tutorials and benchmarks
    tutorials = filter(f -> contains(f, "tutorials"), source_files)
    benchmarks = filter(f -> contains(f, "benchmarking"), source_files)
    
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Mycelia Documentation and Benchmark Reports</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
            .header { background-color: #f0f0f0; padding: 20px; border-radius: 10px; margin-bottom: 20px; }
            .section { margin: 20px 0; }
            .file-list { list-style-type: none; padding: 0; }
            .file-list li { margin: 10px 0; padding: 10px; background-color: #f9f9f9; border-radius: 5px; }
            .file-list a { text-decoration: none; color: #007acc; font-weight: bold; }
            .file-list a:hover { text-decoration: underline; }
            .description { color: #666; margin-top: 5px; }
            .usage { background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin: 10px 0; }
            .usage code { background-color: #e0e0e0; padding: 2px 4px; border-radius: 3px; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Mycelia Documentation and Benchmark Reports</h1>
            <p>Generated: $(now())</p>
            <p>Interactive reports with executed code, outputs, and visualizations</p>
        </div>
        
        <div class="section">
            <h2>📚 Tutorials</h2>
            <p>Step-by-step tutorials for learning Mycelia with executed examples:</p>
            <ul class="file-list">
    """
    
    for tutorial in tutorials
        file_name = basename(tutorial)
        html_name = replace(file_name, ".jl" => ".html")
        title = replace(replace(file_name, ".jl" => ""), "_" => " ")
        title = titlecase(title)
        
        description = get_file_description(tutorial)
        
        html_content *= """
                <li>
                    <a href="$html_name">$title</a>
                    <div class="description">$description</div>
                </li>
        """
    end
    
    html_content *= """
            </ul>
        </div>
        
        <div class="section">
            <h2>⚡ Benchmarks</h2>
            <p>Performance benchmarks and analysis with timing data and metrics:</p>
            <ul class="file-list">
    """
    
    for benchmark in benchmarks
        file_name = basename(benchmark)
        html_name = replace(file_name, ".jl" => ".html")
        title = replace(replace(file_name, ".jl" => ""), "_" => " ")
        title = titlecase(title)
        
        description = get_file_description(benchmark)
        
        html_content *= """
                <li>
                    <a href="$html_name">$title</a>
                    <div class="description">$description</div>
                </li>
        """
    end
    
    html_content *= """
            </ul>
        </div>
        
        <div class="section">
            <h2>🚀 Usage</h2>
            
            <div class="usage">
                <h3>Running Tutorials</h3>
                <code>julia --project=. tutorials/01_data_acquisition.jl</code><br>
                <code>julia --project=. run_extended_tests.jl tutorials</code>
            </div>
            
            <div class="usage">
                <h3>Running Benchmarks</h3>
                <code>julia --project=. run_extended_tests.jl benchmarks</code><br>
                <code>julia --project=. run_extended_tests.jl benchmarks --hpc</code>
            </div>
            
            <div class="usage">
                <h3>Generating Reports</h3>
                <code>julia --project=. run_extended_tests.jl reports</code><br>
                <code>python -m http.server 8000 --directory reports/html</code>
            </div>
        </div>
        
        <div class="section">
            <h2>📋 About</h2>
            <p>These reports are generated automatically from executable Julia scripts using 
            <a href="https://github.com/fredrikekre/Literate.jl">Literate.jl</a> and 
            <a href="https://nbconvert.readthedocs.io/">nbconvert</a>. Each report contains:</p>
            <ul>
                <li>✅ Executed code with real outputs</li>
                <li>📊 Generated plots and visualizations</li>
                <li>📈 Performance metrics and timing data</li>
                <li>📝 Comprehensive documentation and explanations</li>
            </ul>
            <p>For more information, see the <a href="https://cjprybol.github.io/Mycelia/dev/">Mycelia documentation</a>.</p>
        </div>
    </body>
    </html>
    """
    
    index_file = joinpath(html_dir, "index.html")
    open(index_file, "w") do io
        write(io, html_content)
    end
    
    println("✓ Index page created: $index_file")
end

"""
    get_file_description(file_path)

Extract description from file's docstring or comments.
"""
function get_file_description(file_path)
    descriptions = Dict(
        "01_data_acquisition.jl" => "Download genomic data from NCBI and simulate synthetic datasets",
        "02_quality_control.jl" => "Assess and improve sequencing data quality",
        "03_kmer_analysis.jl" => "K-mer counting, analysis, and genome size estimation",
        "04_genome_assembly.jl" => "HiFi genome assembly with hifiasm and quality assessment",
        "05_assembly_validation.jl" => "Comprehensive assembly validation and quality metrics",
        "06_gene_annotation.jl" => "Gene prediction and functional annotation",
        "07_comparative_genomics.jl" => "Pangenome analysis and phylogenetic reconstruction",
        "08_tool_integration.jl" => "External tool integration and workflow management",
        "01_data_processing_benchmark.jl" => "Performance benchmarking of data processing operations",
        "02_kmer_analysis_benchmark.jl" => "K-mer analysis performance across different scales",
        "03_assembly_benchmark.jl" => "Genome assembly performance and quality benchmarking",
        "04_annotation_benchmark.jl" => "Gene annotation performance and accuracy benchmarking",
        "05_comparative_benchmark.jl" => "Comparative genomics and pangenome analysis benchmarking"
    )
    
    file_name = basename(file_path)
    return get(descriptions, file_name, "Mycelia bioinformatics analysis")
end

"""
    run_specific_tutorial(tutorial_name)

Run a specific tutorial by name.
"""
function run_specific_tutorial(tutorial_name)
    # Add .jl extension if missing
    if !endswith(tutorial_name, ".jl")
        tutorial_name *= ".jl"
    end
    
    # Add number prefix if missing (e.g., "data_acquisition" -> "01_data_acquisition.jl")
    if !occursin(r"^\d+_", tutorial_name)
        # Try to find matching tutorial
        tutorial_files = filter(f -> endswith(f, ".jl"), readdir(TUTORIALS_DIR))
        matching = filter(f -> contains(f, tutorial_name), tutorial_files)
        
        if length(matching) == 1
            tutorial_name = matching[1]
        elseif length(matching) > 1
            println("Multiple tutorials match '$tutorial_name':")
            for file in matching
                println("  - $file")
            end
            return false
        else
            println("No tutorial found matching '$tutorial_name'")
            return false
        end
    end
    
    tutorial_path = joinpath(TUTORIALS_DIR, tutorial_name)
    
    if !isfile(tutorial_path)
        println("Tutorial not found: $tutorial_path")
        return false
    end
    
    println("Running tutorial: $tutorial_name")
    start_time = time()
    
    try
        include(tutorial_path)
        elapsed = time() - start_time
        println("✓ Tutorial completed in $(round(elapsed, digits=2))s")
        return true
    catch e
        elapsed = time() - start_time
        println("✗ Tutorial failed after $(round(elapsed, digits=2))s")
        println("Error: $e")
        return false
    end
end

"""
    main()

Main entry point for command-line interface.
"""
function main()
    if length(ARGS) == 0
        show_help()
        return
    end
    
    command = ARGS[1]
    
    if command == "help" || command == "--help" || command == "-h"
        show_help()
        
    elseif command == "tutorials"
        success = run_tutorials()
        exit(success ? 0 : 1)
        
    elseif command == "benchmarks"
        use_hpc = length(ARGS) > 1 && ARGS[2] == "--hpc"
        success = run_benchmarks(use_hpc)
        exit(success ? 0 : 1)
        
    elseif command == "notebooks"
        success = convert_to_notebooks()
        exit(success ? 0 : 1)
        
    elseif command == "reports"
        success = generate_html_reports()
        exit(success ? 0 : 1)
        
    elseif command == "tutorial"
        if length(ARGS) < 2
            println("Error: Please specify tutorial name")
            println("Usage: julia run_extended_tests.jl tutorial <name>")
            exit(1)
        end
        
        tutorial_name = ARGS[2]
        success = run_specific_tutorial(tutorial_name)
        exit(success ? 0 : 1)
        
    else
        println("Error: Unknown command '$command'")
        println("Run 'julia run_extended_tests.jl help' for usage information")
        exit(1)
    end
end

# Execute main function if script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end