using Eisenia
using Documenter

makedocs(;
    modules=[Eisenia],
    authors="Cameron Prybol <cameron.prybol@gmail.com> and contributors",
    repo="https://github.com/cjprybol/Eisenia.jl/blob/{commit}{path}#L{line}",
    sitename="Eisenia.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cjprybol.github.io/Eisenia.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
