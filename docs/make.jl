using Mycelia
using Documenter
import Weave

chapters = String[]

PKG_BASE = dirname(dirname(pathof(Mycelia)))

# cleanup all current .ipynb exports
for converted_notebook_file in filter(x -> occursin(r"\.ipynb\.md", x), readdir("$(PKG_BASE)/docs/src", join=true))
    rm(converted_notebook_file)
end

# weave all current chapter notebooks into markdown files
for notebook in filter(x -> occursin(r"\.ipynb$", x), readdir("$(PKG_BASE)/docs/chapters/", join=true))
    out = basename(notebook) * ".md"
    Weave.weave(notebook, out_path="$(PKG_BASE)/docs/src/$(out)", doctype="github")
    push!(chapters, out)
end

# copy readme from repo to being main index file in documentation
cp("$(PKG_BASE)/README.md", "$(PKG_BASE)/docs/src/README.md", force=true)
for svg in filter(x -> occursin(r"\.svg$", x), readdir(PKG_BASE, join=true))
    cp(svg, "$(PKG_BASE)/docs/src/$(basename(svg))", force=true)
end

makedocs(;
    modules=[Mycelia],
    authors="Cameron Prybol <cameron.prybol@gmail.com> and contributors",
    repo="https://github.com/cjprybol/Mycelia.jl/blob/{commit}{path}#L{line}",
    sitename="Mycelia.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cjprybol.github.io/Mycelia.jl",
        assets=String[],
        collapselevel=1,
    ),
    pages=[
        "Home" => "README.md",
        "Tutorial" => [
        ],
        "Chapters" => chapters,
        "Applications" => [
        ],
        "Reference" => [
            "docstrings.md",
            "neo4j-notes.md"
        ],
        "Nesting" => [
            "nested" => []
        ]
    ],
)

deploydocs(
    repo = "github.com/cjprybol/Mycelia.git",
    push_preview = true,
    deps = nothing,
    make = nothing,
    devurl = "docs"
)
