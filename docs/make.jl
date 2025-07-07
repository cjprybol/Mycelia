# generate docs locally with
# julia --project=docs -e 'include("docs/make.jl")'

using Documenter
using Mycelia
import Literate

const PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
const LITERATE_SRC_DIR = joinpath(PROJECT_ROOT, "test")
const GENERATED_DOCS_DIR = joinpath(@__DIR__, "src", "generated")

# --- TEMPORARY: Only process the example tutorial for debugging Literate integration ---
mkpath(GENERATED_DOCS_DIR)
Literate.markdown(
    joinpath(LITERATE_SRC_DIR, "test-driven-tutorial.jl"),
    GENERATED_DOCS_DIR;
    credit = false,
    flavor = Literate.DocumenterFlavor()
)

# # --- ORIGINAL: Process all .jl files in test/ (reactivate later) ---
# for (root, _, files) in walkdir(LITERATE_SRC_DIR)
#     relative_path = relpath(root, LITERATE_SRC_DIR)
#     for file in filter(f -> endswith(f, ".jl") && f != "runtests.jl", files)
#         input_file = joinpath(root, file)
#         output_dir = joinpath(GENERATED_DOCS_DIR, relative_path)
#         mkpath(output_dir)
#         Literate.markdown(input_file, output_dir, credit = false, flavor = Literate.DocumenterFlavor())
#     end
# end

# Only include the generated tutorial for now
tutorial_pages = [
    joinpath("generated", "test-driven-tutorial.md")
]

makedocs(
    sitename = "Mycelia",
    modules = [Mycelia],
    authors = "Cameron Prybol <cameron.prybol@gmail.com> and contributors",
    repo = "https://github.com/cjprybol/Mycelia/blob/{commit}{path}#L{line}",
    format = Documenter.HTMLWriter.HTML(size_threshold = 1_000_000),
    pages = [
        "Home" => "index.md",
        "Tutorials" => tutorial_pages,
    ]
)

deploydocs(
    repo = "github.com/cjprybol/Mycelia.git",
)