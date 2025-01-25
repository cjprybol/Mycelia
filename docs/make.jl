using Documenter
using Mycelia

# julia --project=. docs/make.jl

makedocs(
    sitename = "Mycelia Documentation",
    modules = [Mycelia],
    format = Documenter.HTMLWriter.HTML(size_threshold = 300_000),  # Increase the size threshold
    pages = [
        "Home" => "index.md",
        # Add other pages here
    ]
)