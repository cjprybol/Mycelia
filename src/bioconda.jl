"""
$(DocStringExtensions.TYPEDSIGNATURES)

Check whether a named Bioconda environment already exists.

# Arguments
- `pkg::String`: Name of the environment.

# Returns
`Bool` indicating if the environment is present.
"""
function check_bioconda_env_is_installed(pkg)
    # ensure conda environment is available
    if !isfile(CONDA_RUNNER)
        if (basename(CONDA_RUNNER) == "mamba")
            Conda.add("mamba")
        elseif (basename(CONDA_RUNNER) == "conda")
            Conda.update()
        end
    end
    # try
    current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
    if pkg in current_environments
        return true
    else
        return false
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a new Conda environment with a specified Bioconda package.

# Arguments
- `pkg::String`: Package name to install. Can include channel specification using 
the format "channel::package"

# Keywords
- `force::Bool=false`: If true, recreates the environment even if it already exists

# Details
The function creates a new Conda environment named after the package and installs
the package into it. It uses channel priority: conda-forge > bioconda > defaults.
If CONDA_RUNNER is set to 'mamba', it will ensure mamba is installed first.

# Examples
```julia
# Install basic package
add_bioconda_env("blast")

# Install from specific channel
add_bioconda_env("bioconda::blast")

# Force reinstallation
add_bioconda_env("blast", force=true)
```
# Notes
- Requires Conda.jl to be installed and configured
- Uses CONDA_RUNNER global variable to determine whether to use conda or mamba
- Cleans conda cache after installation
"""
function add_bioconda_env(pkg; force=false)
channel = nothing
if occursin("::", pkg)
    println("splitting $(pkg)")
    channel, pkg = split(pkg, "::")
    println("into channel:$(channel) pkg:$(pkg)")
end
already_installed = check_bioconda_env_is_installed(pkg)
if !already_installed || force
    @info "installing conda environment $(pkg)"
    if isnothing(channel)
        run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
    else
        run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(channel)::$(pkg) -y`)
    end
    run(`$(CONDA_RUNNER) clean --all -y`)
end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Update a package and its dependencies in its dedicated Conda environment.

# Arguments
- `pkg::String`: Name of the package/environment to update
"""
function update_bioconda_env(pkg)
run(`$(CONDA_RUNNER) update -n $(pkg) $(pkg) -y`)
# conda update --all -n <env_name>
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)
# """
# function add_bioconda_envs(;all=false, force=false)
#     if !isfile(CONDA_RUNNER) && (basename(CONDA_RUNNER) == "mamba")
#         Conda.add("mamba")
#     end
#     if !isfile(joinpath(Conda.BINDIR, "pigz"))
#         run(`$(CONDA_RUNNER) install pigz -y`)
#     end
#     current_environments = Set(first.(filter(x -> length(x) == 2, split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
#     # https://github.com/JuliaPy/Conda.jl/issues/185#issuecomment-1145149905
#     if all
#         for pkg in [
#             "art",
#             # "bioconvert",
#             "badread",
#             "bcftools",
#             "bedtools",
#             "blast",
#             "clair3-illumina",
#             "clair3",    
#             # "bwa",
#             # "bwa-mem2",
#             # "deepvariant",
#             "emboss",
#             "filtlong",
#             # "freebayes",
#             "flye",
#             "gatk4",
#             # "gffread",
#             "htslib",
#             "megahit",
#             "medaka",
#             "minimap2",
#             "mmseqs2",
#             "nanocaller",
#             "nanovar",
#             # "nanoq",
#             # "nanosim",
#             # "nanosim-h",
#             "ncbi-datasets-cli",
#             "pggb",
#             "picard",
#             # "polypolish",
#             "prodigal",
#             "raven-assembler",
#             "rtg-tools",
#             "samtools",
#             "sniffles",
#             "sourmash",
#             "spades",
#             "tabix",
#             "transtermhp",
#             "trim-galore",
#             "vcftools",
#             "vg"
#             ]
#             if !(pkg in current_environments) || force
#                 @info "installing conda environment $(pkg)"
#                 add_bioconda_env(pkg)
#             else
#                 @info "conda environment $(pkg) already present; set force=true to update/re-install"
#             end
#         end
#     end
#     run(`$(CONDA_RUNNER) clean --all -y`)
# end