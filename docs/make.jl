using Documenter
using Mycelia

# format=Documenter.HTML(;
# prettyurls=get(ENV, "CI", "false") == "true",
# canonical="https://cjprybol.github.io/Mycelia",
# assets=String[],
# collapselevel=1,
# ),
makedocs(
    sitename = "Mycelia",
    modules = [Mycelia],
    authors="Cameron Prybol <cameron.prybol@gmail.com> and contributors",
    repo="https://github.com/cjprybol/Mycelia/blob/{commit}{path}#L{line}",
    format = Documenter.HTMLWriter.HTML(size_threshold = 1_000_000),  # Increase the size threshold
    pages = [
        "Home" => "index.md",
        # Add other pages here
    ]
)

deploydocs(
    repo = "github.com/cjprybol/Mycelia.git",
    # push_preview = true,
    # deps = nothing,
    # make = nothing,
    # devurl = "docs"
)