# # Run All Mycelia Tutorials
#
# This script executes all Mycelia tutorials in sequence.
# Use with caution - may take significant time and resources.
#
# ## Usage
#
# ```bash
# # Run all tutorials
# julia --project=. tutorials/run_all_tutorials.jl
#
# # Run specific tutorial
# julia --project=. tutorials/01_data_acquisition.jl
#
# # Convert to notebooks
# julia --project=. -e 'include("tutorials/run_all_tutorials.jl"); convert_all_to_notebooks()'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Dates
import Literate

# List of all tutorials in execution order
TUTORIALS = [
    "01_data_acquisition.jl",
    "02_quality_control.jl",
    "03_kmer_analysis.jl",
    "04_genome_assembly.jl",
    "05_assembly_validation.jl",
    "06_gene_annotation.jl",
    "07_comparative_genomics.jl",
    "08_tool_integration.jl",
]

"""
    run_all_tutorials()

Execute all tutorials in sequence, collecting timing and status information.
"""
function run_all_tutorials()
    println("=== Running All Mycelia Tutorials ===")
    println("Start time: $(Dates.now())")
    println()
    
    results = []
    total_start = time()
    
    for tutorial in TUTORIALS
        if isfile(tutorial)
            println("Running $tutorial...")
            start_time = time()
            
            try
                ## Execute tutorial
                include(tutorial)
                
                elapsed = time() - start_time
                push!(results, (tutorial, :success, elapsed))
                println("✓ $tutorial completed in $(round(elapsed, digits=2))s")
                
            catch e
                elapsed = time() - start_time
                push!(results, (tutorial, :error, elapsed, e))
                println("✗ $tutorial failed after $(round(elapsed, digits=2))s")
                println("  Error: $e")
            end
            
            println()
        else
            println("⚠ $tutorial not found, skipping...")
            push!(results, (tutorial, :not_found, 0.0))
        end
    end
    
    total_elapsed = time() - total_start
    
    ## Print summary
    println("=== Tutorial Execution Summary ===")
    println("Total time: $(round(total_elapsed, digits=2))s")
    println()
    
    successful = 0
    for (tutorial, status, elapsed, error...) in results
        if status == :success
            println("✓ $tutorial ($(round(elapsed, digits=2))s)")
            successful += 1
        elseif status == :error
            println("✗ $tutorial ($(round(elapsed, digits=2))s) - $(error[1])")
        else
            println("⚠ $tutorial (not found)")
        end
    end
    
    println()
    println("Successful: $successful/$(length(TUTORIALS))")
    println("End time: $(Dates.now())")
    
    return results
end

"""
    convert_all_to_notebooks()

Convert all tutorials to Jupyter notebooks using Literate.jl.
"""
function convert_all_to_notebooks()
    println("=== Converting Tutorials to Notebooks ===")
    
    ## Create notebooks directory if it doesn't exist
    notebooks_dir = "notebooks"
    if !isdir(notebooks_dir)
        mkdir(notebooks_dir)
    end
    
    for tutorial in TUTORIALS
        if isfile(tutorial)
            println("Converting $tutorial to notebook...")
            
            try
                ## Convert to notebook
                Literate.notebook(tutorial, notebooks_dir, execute=false)
                println("✓ Converted $tutorial")
                
            catch e
                println("✗ Failed to convert $tutorial: $e")
            end
        else
            println("⚠ $tutorial not found, skipping...")
        end
    end
    
    println()
    println("Notebooks saved to: $notebooks_dir")
end

"""
    run_tutorial(tutorial_name)

Run a specific tutorial by name.
"""
function run_tutorial(tutorial_name)
    if !endswith(tutorial_name, ".jl")
        tutorial_name *= ".jl"
    end
    
    if isfile(tutorial_name)
        println("Running $tutorial_name...")
        start_time = time()
        
        try
            include(tutorial_name)
            elapsed = time() - start_time
            println("✓ $tutorial_name completed in $(round(elapsed, digits=2))s")
            return true
        catch e
            elapsed = time() - start_time
            println("✗ $tutorial_name failed after $(round(elapsed, digits=2))s")
            println("  Error: $e")
            return false
        end
    else
        println("✗ $tutorial_name not found")
        return false
    end
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    ## Only run if this script is executed directly
    if length(ARGS) == 0
        ## Run all tutorials
        results = run_all_tutorials()
    elseif ARGS[1] == "convert"
        ## Convert to notebooks
        convert_all_to_notebooks()
    else
        ## Run specific tutorial
        tutorial_name = ARGS[1]
        run_tutorial(tutorial_name)
    end
end