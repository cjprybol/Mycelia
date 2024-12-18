using Mycelia
using Documenter
import Weave
import DataStructures

PKG_BASE = dirname(dirname(pathof(Mycelia)))

mkpath("$(PKG_BASE)/docs/src")
# copy readme from repo to being main index file in documentation
cp("$(PKG_BASE)/README.md", "$(PKG_BASE)/docs/src/index.md", force=true)

for documentation_group in filter(x -> isdir(x) && !occursin(r"^\.", basename(x)), readdir("$(PKG_BASE)/docs/_src", join=true))
    # convert all current .ipynb notebooks to markdown files
#     @show documentation_group
#     @show basename(documentation_group)
    outdir = joinpath("$(PKG_BASE)/docs/src", basename(documentation_group))
    mkpath(outdir)
    if basename(documentation_group) == "5.Development"
        for notebook in filter(x -> occursin(r".ipynb$", x), readdir(documentation_group, join=true))
            out_path = joinpath("$(PKG_BASE)/docs/src", basename(documentation_group), basename(notebook) * ".md")
            Weave.convert_doc(notebook, out_path)
        end
    #     copy non notebook files over
        for non_notebook in filter(x -> !occursin(r".ipynb$", x), readdir(documentation_group, join=true))
            out_path = joinpath("$(PKG_BASE)/docs/src", basename(documentation_group), basename(non_notebook))
            cp(non_notebook, out_path, force=true)
        end
    else
        for notebook in filter(x -> occursin(r".ipynb$", x), readdir(documentation_group, join=true))
            out_path = joinpath("$(PKG_BASE)/docs/src", basename(documentation_group), basename(notebook) * ".md")
            Weave.weave(notebook, out_path=out_path, doctype="github")
        end
    #     copy non notebook files over
        for non_notebook in filter(x -> !occursin(r".ipynb$", x), readdir(documentation_group, join=true))
            out_path = joinpath("$(PKG_BASE)/docs/src", basename(documentation_group), basename(non_notebook))
            cp(non_notebook, out_path, force=true)
        end
    end
end

pages = DataStructures.SortedDict()
for documentation_group in filter(x -> isdir(x) && !occursin(r"^\.", basename(x)), readdir("$(PKG_BASE)/docs/src", join=true))
    md_pages = filter(x -> occursin(r"\.md", x), readdir(documentation_group))
    md_pages = map(page -> "$(basename(documentation_group))/$(page)", md_pages)    
    pages[basename(documentation_group)] = md_pages
end

pages=[
    "Home" => "index.md",
    collect(pages)...
]

makedocs(;
    modules=[Mycelia],
    authors="Cameron Prybol <cameron.prybol@gmail.com> and contributors",
    repo="https://github.com/cjprybol/Mycelia/blob/{commit}{path}#L{line}",
    sitename="Mycelia",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cjprybol.github.io/Mycelia",
        assets=String[],
        collapselevel=1,
    ),
    pages=pages
)

deploydocs(
    repo = "github.com/cjprybol/Mycelia.git",
    push_preview = true,
    deps = nothing,
    make = nothing,
    devurl = "docs"
)

# rm("$(PKG_BASE)/docs/src", recursive=true)
