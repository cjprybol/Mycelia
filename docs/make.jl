using Mycelia
using Documenter
import Weave
import DataStructures

# chapters = String[]

PKG_BASE = dirname(dirname(pathof(Mycelia)))

# copy readme from repo to being main index file in documentation
cp("$(PKG_BASE)/README.md", "$(PKG_BASE)/docs/src/README.md", force=true)
# for svg in filter(x -> occursin(r"\.svg$", x), readdir(PKG_BASE, join=true))
#     cp(svg, "$(PKG_BASE)/docs/src/$(basename(svg))", force=true)
# end

pages = DataStructures.SortedDict()

for documentation_group in filter(x -> isdir(x) && !occursin(r"^\.", basename(x)), readdir("$(PKG_BASE)/docs/src", join=true))
    
    @show documentation_group
    # cleanup all current .ipynb exports
    for converted_notebook_file in filter(x -> occursin(r"\.ipynb\.md", x), readdir(documentation_group, join=true))
        @show converted_notebook_file
        rm(converted_notebook_file)
    end
    
    # convert all current .ipynb notebooks to markdown files
    for notebook in filter(x -> occursin(r".ipynb$", x), readdir(documentation_group, join=true))
        @show notebook
        out_path = notebook * ".md" 
        Weave.weave(notebook, out_path=out_path, doctype="github")
#         Weave.weave(notebook, doctype="github")
    end
    
    md_pages = filter(x -> occursin(r"\.md", x), readdir(documentation_group))
    md_pages = map(page -> "$(basename(documentation_group))/$(page)", md_pages)
    
    @show md_pages
    
#     pages[basename(documentation_group)] = 
end

# makedocs(;
#     modules=[Mycelia],
#     authors="Cameron Prybol <cameron.prybol@gmail.com> and contributors",
#     repo="https://github.com/cjprybol/Mycelia.jl/blob/{commit}{path}#L{line}",
#     sitename="Mycelia.jl",
#     format=Documenter.HTML(;
#         prettyurls=get(ENV, "CI", "false") == "true",
#         canonical="https://cjprybol.github.io/Mycelia.jl",
#         assets=String[],
#         collapselevel=1,
#     ),
#     pages=[
#         "Home" => "README.md",
#         collect(pages)...
#     ],
# )

# makedocs(;
#     modules=[Mycelia],
#     authors="Cameron Prybol <cameron.prybol@gmail.com> and contributors",
#     repo="https://github.com/cjprybol/Mycelia.jl/blob/{commit}{path}#L{line}",
#     sitename="Mycelia.jl",
#     format=Documenter.HTML(;
#         prettyurls=get(ENV, "CI", "false") == "true",
#         canonical="https://cjprybol.github.io/Mycelia.jl",
#         assets=String[],
#         collapselevel=1,
#     ),
#     pages=[
#         "Home" => "README.md",
#         "Tutorial" => [
#             "Constructing a graph from reference sequences",
#             "Constructing a graph from observations",
#             "Merging graphs",
#             "Subsetting graphs",
#             "Adding annotations to graphs",
#             "Adding taxonomy to graphs"
#         ],
#         "Applications" => [
#         ],
#         "Discussion" => chapters,
#         "Reference" => [
#             "docstrings.md",
#             "neo4j-notes.md"
#         ],
#         "Nesting" => [
#             "nested" => []
#         ]
#     ],
# )

# deploydocs(
#     repo = "github.com/cjprybol/Mycelia.git",
#     push_preview = true,
#     deps = nothing,
#     make = nothing,
#     devurl = "docs"
# )
