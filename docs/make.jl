using Eisenia
using Documenter

makedocs(;
    modules=[Eisenia],
    authors="Cameron Prybol <cameron.prybol@gmail.com> and contributors",
    repo="https://github.com/Cameron Prybol/Eisenia.jl/blob/{commit}{path}#L{line}",
    sitename="Eisenia.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Cameron Prybol.gitlab.io/Eisenia.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
