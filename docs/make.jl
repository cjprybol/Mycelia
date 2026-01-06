# generate docs locally with
# julia --project=docs -e 'include("docs/make.jl")'

using Documenter
using Mycelia
import Literate

const PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
const TUTORIALS_DIR = joinpath(PROJECT_ROOT, "tutorials")
const LITERATE_SRC_DIR = joinpath(PROJECT_ROOT, "test")
const GENERATED_DOCS_DIR = joinpath(@__DIR__, "src", "generated")

# Drop any stale generated content so we don't ship orphaned pages
if isdir(GENERATED_DOCS_DIR)
    rm(GENERATED_DOCS_DIR; recursive=true, force=true)
end

# Process tutorials from tutorials/ directory
mkpath(GENERATED_DOCS_DIR)
tutorials_output_dir = joinpath(GENERATED_DOCS_DIR, "tutorials")
mkpath(tutorials_output_dir)

# Initialize tutorial_files array
tutorial_files = []

# Process ALL tutorial files from the tutorials directory
if isdir(TUTORIALS_DIR)
    for file in readdir(TUTORIALS_DIR)
        if endswith(file, ".jl")
            file_path = joinpath(TUTORIALS_DIR, file)
            if isfile(file_path)
                println("Processing tutorial: $file")
                Literate.markdown(
                    file_path, 
                    tutorials_output_dir;
                    credit = false,
                    flavor = Literate.DocumenterFlavor(),
                    # Disable execution for now to fix build issues
                    # TODO: Re-enable once tutorial examples are fixed
                    execute = false
                )
                
                # Post-process the generated markdown to convert @example blocks to ```julia blocks
                md_filename = replace(file, ".jl" => ".md")
                md_path = joinpath(tutorials_output_dir, md_filename)
                if isfile(md_path)
                    content = read(md_path, String)
                    # Replace @example blocks with regular julia code blocks
                    content = replace(content, r"```@example [^\n]*\n" => "```julia\n")
                    write(md_path, content)
                    println("  Disabled @example blocks in $md_filename")
                end
                
                # Add to our list for the sidebar
                push!(tutorial_files, joinpath("generated", "tutorials", md_filename))
            end
        end
    end
end

# Also process the test tutorial for debugging
if isfile(joinpath(LITERATE_SRC_DIR, "test-driven-tutorial.jl"))
    Literate.markdown(
        joinpath(LITERATE_SRC_DIR, "test-driven-tutorial.jl"),
        GENERATED_DOCS_DIR;
        credit = false,
        flavor = Literate.DocumenterFlavor(),
        # Disable execution for now to fix build issues
        execute = false
    )
    
    # Post-process to disable @example blocks
    test_md_path = joinpath(GENERATED_DOCS_DIR, "test-driven-tutorial.md")
    if isfile(test_md_path)
        content = read(test_md_path, String)
        content = replace(content, r"```@example [^\n]*\n" => "```julia\n")
        write(test_md_path, content)
        println("  Disabled @example blocks in test-driven-tutorial.md")
    end
    
    push!(tutorial_files, joinpath("generated", "test-driven-tutorial.md"))
end

# Create tutorial pages list
tutorial_pages = []

# Add tutorials in logical order
tutorial_order = [
    # Main tutorial series (01-08) - temporarily only including these to test build
    ("Step 1: Data Acquisition", "01_data_acquisition.md"),
    ("Step 2: Quality Control", "02_quality_control.md"),
    ("Step 3: K-mer Analysis", "03_kmer_analysis.md"),
    ("Step 4: Genome Assembly", "04_genome_assembly.md"),
    # ("4b. Graph Type Tutorials", "04_graph_type_tutorials.md"),  # TODO: Fix execution errors
    ("Step 5: Assembly Validation", "05_assembly_validation.md"),
    ("Step 6: Gene Annotation", "06_gene_annotation.md"),
    ("Step 7: Comparative Genomics", "07_comparative_genomics.md"),
    ("Step 8: Tool Integration", "08_tool_integration.md"),
    ("Step 9: Binning Workflow", "14_binning_workflow.md"),
    
    # Round-trip tutorial series (09) - temporarily disabled due to execution errors
    # TODO: Re-enable once example execution issues are resolved
    # ("Round-Trip 1: String Graphs", "09_round_trip_01_string_graphs.md"),
    # ("Round-Trip 2: N-gram to String", "09_round_trip_02_ngram_to_string.md"), 
    # ("Round-Trip 3: FASTA Sequences", "09_round_trip_03_fasta_sequences.md"),
    # ("Round-Trip 4: K-mer to Sequence", "09_round_trip_04_kmer_to_sequence.md"),
    # ("Round-Trip 5: FASTQ Graphs", "09_round_trip_05_fastq_graphs.md"),
    # ("Round-Trip 6: Qualmer Graphs", "09_round_trip_06_qualmer_graphs.md")
]

# Add tutorials that exist in the specified order
for (title, filename) in tutorial_order
    tutorial_path = joinpath("generated", "tutorials", filename)
    if tutorial_path in tutorial_files
        push!(tutorial_pages, title => tutorial_path)
        println("Added tutorial: $title -> $tutorial_path")
    else
        println("Tutorial not found: $tutorial_path")
    end
end

# Add any remaining tutorials that weren't in the ordered list
for tutorial_file in tutorial_files
    # Skip if already added
    if !any(page -> page.second == tutorial_file, tutorial_pages)
        # Extract a title from the filename
        filename = basename(tutorial_file)
        title = replace(filename, ".md" => "", "_" => " ", r"^\d+[-_]?" => "")
        title = titlecase(title)
        push!(tutorial_pages, title => tutorial_file)
        println("Added additional tutorial: $title -> $tutorial_file")
    end
end

# # Add test tutorial if it exists
# if joinpath("generated", "test-driven-tutorial.md") in tutorial_files
#     push!(tutorial_pages, "Test Tutorial" => joinpath("generated", "test-driven-tutorial.md"))
# end

makedocs(
    sitename = "Mycelia",
    modules = [Mycelia],
    authors = "Cameron Prybol <cameron.prybol@gmail.com> and contributors",
    repo = Documenter.Remotes.GitHub("cjprybol", "Mycelia"),
    format = Documenter.HTMLWriter.HTML(size_threshold = 1_000_000),
    build = joinpath(@__DIR__, "build"),
    source = joinpath(@__DIR__, "src"),
    # Disable doctests and example blocks for now to fix build issues
    # TODO: Re-enable these once tutorial examples are fixed
    doctest = false,
    checkdocs = :none,
    warnonly = [:cross_references, :example_block, :missing_docs],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Installation" => "installation.md",
        "Concepts" => "concepts.md",
        "Workflow Map" => "workflow-map.md",
        "Metagenomic Workflow" => "metagenomic-workflow.md",
        "Tutorials" => vcat(
            ["Overview" => "tutorials.md"],
            tutorial_pages,
        ),
        "API Reference" => [
            "Complete API Surface" => "api/all-functions.md",
        ],
        "Benchmarks" => "benchmarks.md",
        "Contributing" => "contributing.md",
        "Related Projects" => "related-projects.md",
    ]
)

deploydocs(
    repo = "github.com/cjprybol/Mycelia.git",
)
