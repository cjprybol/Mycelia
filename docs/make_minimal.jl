# Minimal documentation build that doesn't load Mycelia module
# This allows us to build docs despite network restrictions

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

# Create a mock Mycelia module for Documenter
module Mycelia
    # Minimal module interface for documentation
end

# Install minimal deps in a temporary environment
import Pkg
Pkg.activate(mktempdir())
Pkg.add("Documenter")
Pkg.add("Literate")

using Documenter
import Literate

const PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
const TUTORIALS_DIR = joinpath(PROJECT_ROOT, "tutorials")
const GENERATED_DOCS_DIR = joinpath(@__DIR__, "src", "generated")

# Process tutorials
mkpath(GENERATED_DOCS_DIR)
tutorials_output_dir = joinpath(GENERATED_DOCS_DIR, "tutorials")
mkpath(tutorials_output_dir)

tutorial_files = []

println("Processing tutorials from: $TUTORIALS_DIR")
if isdir(TUTORIALS_DIR)
    for file in sort(readdir(TUTORIALS_DIR))
        if endswith(file, ".jl")
            file_path = joinpath(TUTORIALS_DIR, file)
            if isfile(file_path)
                println("  Processing: $file")
                try
                    Literate.markdown(
                        file_path, 
                        tutorials_output_dir;
                        credit = false,
                        flavor = Literate.DocumenterFlavor(),
                        execute = false
                    )
                    
                    md_filename = replace(file, ".jl" => ".md")
                    md_path = joinpath(tutorials_output_dir, md_filename)
                    if isfile(md_path)
                        content = read(md_path, String)
                        content = replace(content, r"```@example [^\n]*\n" => "```julia\n")
                        write(md_path, content)
                    end
                    
                    push!(tutorial_files, joinpath("generated", "tutorials", md_filename))
                catch e
                    println("    ERROR: $e")
                end
            end
        end
    end
end

println("\n✓ Processed $(length(tutorial_files)) tutorial files")

# Build minimal docs
println("\nBuilding documentation...")
makedocs(
    sitename = "Mycelia",
    modules = [Mycelia],
    authors = "Cameron Prybol <cameron.prybol@gmail.com> and contributors",
    repo = "https://github.com/cjprybol/Mycelia/blob/{commit}{path}#L{line}",
    format = Documenter.HTMLWriter.HTML(size_threshold = 1_000_000),
    build = joinpath(@__DIR__, "build"),
    source = joinpath(@__DIR__, "src"),
    doctest = false,
    checkdocs = :none,
    warnonly = true,  # Convert all errors to warnings
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Getting Started" => "getting-started.md",
    ]
)

println("\n✓ Documentation build complete!")
println("Output directory: $(joinpath(@__DIR__, "build"))")
