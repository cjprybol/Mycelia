using Mycelia
using Documenter

# run literate or pandoc .ipynb -> .md here

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
        "Docuentation" => [
            "page1.md",
            "page2.md"
        ],
        "Chapters" => [
            "chapter1.md",
            "chapter2.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/cjprybol/Mycelia.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)
