using Mycelia
using Documenter
import Weave

# run literate or pandoc .ipynb -> .md here

chapters = String[]

println(pwd())
for dir in readdir(pwd())
    @show dir
end

for notebook in filter(x -> occursin(r"\.ipynb$", x), readdir("chapters/", join=true))
    out = replace(basename(notebook), ".ipynb" => ".md")
#     out = replace(basename(notebook), ".ipynb" => ".html")
#     out = "src/$(html_out)"
#     Weave.weave(notebook, out_path=out)
#     Weave.convert_doc(notebook, out, outformat="markdown")
    Weave.weave(notebook, out_path="src/$(out)", doctype="github")
    push!(chapters, out)
end

cp("../README.md", "src/index.md", force=true)
for svg in filter(x -> occursin(r"\.svg$", x), readdir("..", join=true))
    cp(svg, "src/$(basename(svg))", force=true)
end

@show "here"
# Literate.markdown(inputfile, outputdir=pwd(); config::Dict=Dict(), kwargs...)

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

# deploydocs(
#     repo = "github.com/cjprybol/Mycelia.git",
#     push_preview = true,
#     deps = nothing,
#     make = nothing,
#     devurl = "docs"
# )
