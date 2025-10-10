#!/usr/bin/env julia
# Documentation build verification script
# Usage: julia docs/verify_build.jl

println("=== Mycelia Documentation Build Verification ===\n")

# Change to docs directory
cd(dirname(@__FILE__))
project_root = abspath(joinpath(@__DIR__, ".."))

println("1. Checking environment...")
println("   Project root: $project_root")
println("   Docs directory: $(pwd())")

# Check for required files
println("\n2. Checking required files...")
required_files = [
    "make.jl",
    "src/index.md",
    "src/installation.md",
    "src/getting-started.md",
    "src/api-reference.md"
]

all_exist = true
for file in required_files
    if isfile(file)
        println("   ✓ $file")
    else
        println("   ✗ MISSING: $file")
        all_exist = false
    end
end

if !all_exist
    println("\n❌ Some required files are missing!")
    exit(1)
end

# Check tutorials
println("\n3. Checking tutorials...")
tutorials_dir = joinpath(project_root, "tutorials")
if isdir(tutorials_dir)
    tutorial_files = filter(f -> endswith(f, ".jl"), readdir(tutorials_dir))
    println("   Found $(length(tutorial_files)) tutorial files")
    for (i, f) in enumerate(tutorial_files[1:min(5, length(tutorial_files))])
        println("     - $f")
    end
    if length(tutorial_files) > 5
        println("     ... and $(length(tutorial_files) - 5) more")
    end
else
    println("   ⚠ Tutorials directory not found")
end

# Try to build
println("\n4. Building documentation...")
println("   This may take 2-3 minutes...\n")

build_start = time()
try
    include("make.jl")
    build_time = round(time() - build_start, digits=1)
    
    println("\n✅ Documentation build SUCCESSFUL!")
    println("   Build time: $(build_time)s")
    
    # Check output
    if isfile("build/index.html")
        file_size = filesize("build/index.html")
        println("   Output: build/index.html ($(round(file_size/1024, digits=1)) KB)")
    end
    
    if isdir("build")
        html_files = length(filter(f -> endswith(f, ".html"), 
                                   [joinpath(root, file) 
                                    for (root, dirs, files) in walkdir("build") 
                                    for file in files]))
        println("   Total HTML pages: $html_files")
    end
    
    println("\n📖 View documentation:")
    println("   file://$(pwd())/build/index.html")
    
catch e
    build_time = round(time() - build_start, digits=1)
    println("\n❌ Documentation build FAILED after $(build_time)s")
    println("Error: $e")
    exit(1)
end

println("\n=== Verification Complete ===")
