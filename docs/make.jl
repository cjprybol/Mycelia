using Mycelia
using Documenter
import Weave

# run literate or pandoc .ipynb -> .md here

chapters = String[]

PKG_BASE = dirname(dirname(pathof(Mycelia)))

for notebook in filter(x -> occursin(r"\.ipynb$", x), readdir("$(PKG_BASE)/docs/chapters/", join=true))
    out = replace(basename(notebook), ".ipynb" => ".md")
#     out = replace(basename(notebook), ".ipynb" => ".html")
#     Weave.convert_doc(notebook, out, outformat="markdown")
    Weave.weave(notebook, out_path="$(PKG_BASE)/docs/src/$(out)", doctype="github")
    push!(chapters, out)
end

cp("$(PKG_BASE)/README.md", "$(PKG_BASE)/docs/src/index.md", force=true)
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
    ),
    pages=[
        "Home" => "index.md",
        "Documentation" => [
            "docstrings.md",
        ],
        "Chapters" => chapters
    ],
)

deploydocs(
    repo = "github.com/cjprybol/Mycelia.git",
    push_preview = true,
    deps = nothing,
    make = nothing,
    devurl = "docs"
)
