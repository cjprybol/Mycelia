function _ensure_conda_env_vars!()
    conda_exe = CONDA_RUNNER
    if isfile(conda_exe)
        # Prefer Conda.jl's conda to avoid stale module paths.
        if !isfile(get(ENV, "CONDA_EXE", "")) || get(ENV, "CONDA_EXE", "") != conda_exe
            ENV["CONDA_EXE"] = conda_exe
        end
        conda_python = joinpath(dirname(conda_exe), "python")
        if isfile(conda_python)
            if !isfile(get(ENV, "CONDA_PYTHON_EXE", "")) ||
               get(ENV, "CONDA_PYTHON_EXE", "") != conda_python
                ENV["CONDA_PYTHON_EXE"] = conda_python
            end
        end
    end
end

"""
    _conda_env_prefix(env_name::AbstractString) -> Union{Nothing, String}

Return the filesystem prefix for a named conda environment, or `nothing` when
it cannot be determined.
"""
function _conda_env_prefix(env_name::AbstractString)
    _ensure_conda_env_vars!()
    try
        for line in
            filter(x -> !occursin(r"^#", x), readlines(`$(Mycelia.CONDA_RUNNER) env list`))
            m = match(r"^(\S+)\s+(\*)?\s*(.+?)\s*$", line)
            if !isnothing(m) && m.captures[1] == env_name
                return m.captures[3]
            end
        end
    catch e
        @debug "conda env prefix lookup failed" exception = e
        return nothing
    end
    return nothing
end

"""
    _conda_env_variable(env_name::AbstractString, variable_name::AbstractString) -> Union{Nothing, String}

Read a single environment variable from a named conda environment, returning
`nothing` if the environment or variable is unavailable.
"""
function _conda_env_variable(env_name::AbstractString, variable_name::AbstractString)
    _ensure_conda_env_vars!()
    if !isfile(Mycelia.CONDA_RUNNER)
        return nothing
    end
    try
        cmd = Cmd([
            Mycelia.CONDA_RUNNER, "run", "-n", env_name, "python", "-c",
            "import os; print(os.environ.get($(repr(variable_name)), ''))"
        ])
        value = strip(read(cmd, String))
        return isempty(value) ? nothing : value
    catch e
        @debug "conda env variable lookup failed" exception = e
        return nothing
    end
end

"""
    _vibrant_data_path_candidates_from_prefix(env_prefix::AbstractString) -> Vector{String}

Return candidate VIBRANT data directories under a conda environment prefix.
"""
function _vibrant_data_path_candidates_from_prefix(env_prefix::AbstractString)
    share_dir = joinpath(env_prefix, "share")
    if !isdir(share_dir)
        return String[]
    end
    candidates = joinpath.(
        Ref(share_dir),
        filter(name -> startswith(name, "vibrant"), readdir(share_dir)),
        Ref("db")
    )
    sort!(candidates)
    return candidates
end

const VIBRANT_REQUIRED_DATA_FILES = (
    joinpath("databases", "KEGG_profiles_prokaryotes.HMM"),
    joinpath("databases", "KEGG_profiles_prokaryotes.HMM.h3f"),
    joinpath("databases", "KEGG_profiles_prokaryotes.HMM.h3i"),
    joinpath("databases", "KEGG_profiles_prokaryotes.HMM.h3m"),
    joinpath("databases", "KEGG_profiles_prokaryotes.HMM.h3p"),
    joinpath("databases", "Pfam-A_v32.HMM"),
    joinpath("databases", "Pfam-A_v32.HMM.h3f"),
    joinpath("databases", "Pfam-A_v32.HMM.h3i"),
    joinpath("databases", "Pfam-A_v32.HMM.h3m"),
    joinpath("databases", "Pfam-A_v32.HMM.h3p"),
    joinpath("databases", "VOGDB94_phage.HMM"),
    joinpath("databases", "VOGDB94_phage.HMM.h3f"),
    joinpath("databases", "VOGDB94_phage.HMM.h3i"),
    joinpath("databases", "VOGDB94_phage.HMM.h3m"),
    joinpath("databases", "VOGDB94_phage.HMM.h3p"),
    joinpath("files", "VIBRANT_machine_model.sav"),
    joinpath("files", "VIBRANT_names.tsv")
)

"""
    _vibrant_databases_exist(data_path::AbstractString) -> Bool

Return `true` when a VIBRANT data directory contains the downloaded and pressed
database artifacts required by the upstream setup script.
"""
function _vibrant_databases_exist(data_path::AbstractString)
    if !isdir(data_path)
        return false
    end
    return all(isfile(joinpath(data_path, relpath))
    for relpath in VIBRANT_REQUIRED_DATA_FILES)
end

"""
    _vibrant_data_path_candidates(env_name::AbstractString="vibrant") -> Vector{String}

Return possible VIBRANT data directories for a named conda environment.
"""
function _vibrant_data_path_candidates(env_name::AbstractString = "vibrant")
    candidates = String[]
    data_path = _conda_env_variable(env_name, "VIBRANT_DATA_PATH")
    if !isnothing(data_path)
        push!(candidates, data_path)
    end
    if haskey(ENV, "VIBRANT_DATA_PATH") && !isempty(strip(ENV["VIBRANT_DATA_PATH"]))
        push!(candidates, strip(ENV["VIBRANT_DATA_PATH"]))
    end
    env_prefix = _conda_env_prefix(env_name)
    if !isnothing(env_prefix)
        append!(candidates, _vibrant_data_path_candidates_from_prefix(env_prefix))
    end
    return unique(candidates)
end

"""
    _vibrant_database_path(env_name::AbstractString="vibrant") -> Union{Nothing, String}

Return the first detected VIBRANT data directory containing downloaded
databases, or `nothing` when none are found.
"""
function _vibrant_database_path(env_name::AbstractString = "vibrant")
    for candidate in _vibrant_data_path_candidates(env_name)
        if _vibrant_databases_exist(candidate)
            return candidate
        end
    end
    return nothing
end

"""
    _install_vibrant() -> Union{Nothing, String}

Install VIBRANT (Virus Identification By iteRative ANnoTation) from Bioconda and download its required databases.

# Details
This function performs two main tasks:
1. Creates a new conda environment named 'vibrant' and installs the VIBRANT package
2. Downloads the required HMM databases using VIBRANT's download-db.sh script
   when they are not already present

# Implementation Notes
- The database location is fixed to VIBRANT's default location
- Uses CONDA_RUNNER to execute commands within the vibrant environment
- Returns the detected VIBRANT database path when one is available

# Internal Use
This function is called internally by `run_vibrant` to ensure VIBRANT is properly installed
before attempting to run analyses.
"""
function _install_vibrant()
    # install VIBRANT from Bioconda
    # run(`conda install -y -c bioconda vibrant`)
    Mycelia.add_bioconda_env("vibrant")

    existing_data_path = _vibrant_database_path("vibrant")
    if !isnothing(existing_data_path)
        @info "VIBRANT databases already present at $(existing_data_path); skipping download."
        return existing_data_path
    end

    # download required HMM databases
    # run(`bash -lc "download-db.sh"`)
    # don't love the default location, but that's ok for now - no obvious way to set to a different location
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n vibrant download-db.sh`)
    downloaded_data_path = _vibrant_database_path("vibrant")
    if isnothing(downloaded_data_path)
        error("VIBRANT database download completed, but the expected data files were not detected automatically.")
    end
    return downloaded_data_path
    # Likely unused optional arguments

    # -d: specify the location of the databases/ directory if moved from its default location.
    # -m: specify the location of the files/ directory if moved from its default location.
    # For -d and -m please specify the full path to the new location of the files. For example, -d new_location/databases/.
end

"""
    _setup_phageboost_environment(force_reinstall::Bool=false)

Set up the PhageBoost conda environment if it doesn't exist or if force_reinstall is true.
"""
function _setup_phageboost_environment(force_reinstall::Bool = false)
    env_name = "phageboost_env"

    if force_reinstall && _check_conda_env_exists(env_name)
        println("Removing existing PhageBoost environment...")
        try
            Base.run(`$(Mycelia.CONDA_RUNNER) env remove -n $env_name -y`)
        catch e
            @warn "Failed to remove existing environment: $e"
        end
    end

    if !_check_conda_env_exists(env_name) || force_reinstall
        println("Creating PhageBoost conda environment...")
        try
            Base.run(`$(Mycelia.CONDA_RUNNER) create -y -n $env_name python=3.7`)
            println("Installing PhageBoost...")
            # Base.run(`bash -lc "$(Mycelia.CONDA_RUNNER) activate $env_name && pip install PhageBoost"`)
            Base.run(`$(Mycelia.CONDA_RUNNER)  run --live-stream -n $env_name pip install PhageBoost`)
            println("PhageBoost environment setup completed")
        catch e
            throw(ErrorException("Failed to setup PhageBoost environment: $e"))
        end
    else
        println("PhageBoost environment already exists")
    end
end

# """
#     _validate_phageboost_installation()

# Validate that PhageBoost is properly installed and accessible.
# """
# function _validate_phageboost_installation()
#     println("Validating PhageBoost installation...")
#     try
#         # Test if PhageBoost command is available and shows help
#         result = Base.read(`bash -lc "conda activate phageboost_env && PhageBoost --help"`, String)
#         if Base.occursin("PhageBoost", result)
#             println("PhageBoost installation validated successfully")
#         else
#             throw(ErrorException("PhageBoost command available but output unexpected"))
#         end
#     catch e
#         throw(ErrorException("PhageBoost validation failed. Please check installation: $e"))
#     end
# end

# """
#     cleanup_phageboost_environment()

# Remove the PhageBoost conda environment.
# """
# function cleanup_phageboost_environment()
#     env_name = "phageboost_env"
#     if _check_conda_env_exists(env_name)
#         println("Removing PhageBoost environment...")
#         try
#             Base.run(`conda env remove -n $env_name -y`)
#             println("PhageBoost environment removed successfully")
#         catch e
#             @warn "Failed to remove PhageBoost environment: $e"
#         end
#     else
#         println("PhageBoost environment does not exist")
#     end
# end

"""
    _check_conda_env_exists(env_name::AbstractString) -> Bool

Check if a conda environment exists.
"""
function _check_conda_env_exists(env_name::AbstractString)
    _ensure_conda_env_vars!()
    try
        result = Base.read(`$(Mycelia.CONDA_RUNNER) env list`, String)
        return Base.occursin(env_name, result)
    catch e
        @debug "conda env existence check failed" exception = e
        return false
    end
end

"""
    conda_tool_version(env_name::AbstractString, cmd_parts::Vector{String})

Return the first line of a tool version string from a conda environment, or
missing when the tool or environment is unavailable.
"""
function conda_tool_version(env_name::AbstractString, cmd_parts::Vector{String})
    _ensure_conda_env_vars!()
    if !isfile(Mycelia.CONDA_RUNNER)
        return missing
    end
    try
        cmd = Cmd(vcat([Mycelia.CONDA_RUNNER, "run", "-n", env_name], cmd_parts))
        output = read(cmd, String)
        line = strip(first(split(output, '\n')))
        return isempty(line) ? missing : line
    catch e
        @debug "conda tool version check failed" exception = e
        return missing
    end
end

"""
    TOOL_VERSION_REGISTRY

Registry mapping tool environment names to their version-check commands and
expected major versions. Used by [`verify_tool_version`](@ref) for pre-flight
version validation.

Each entry maps `env_name => (cmd_parts, expected_major)`.
"""
const TOOL_VERSION_REGISTRY = Dict{String, Tuple{Vector{String}, Int}}(
    "metaphlan" => (["metaphlan", "--version"], 4),
    "kraken2" => (["kraken2", "--version"], 2),
    "megahit" => (["megahit", "--version"], 1),
    "spades" => (["spades.py", "--version"], 4),
    "busco" => (["busco", "--version"], 5),
    "checkm2" => (["checkm2", "--version"], 1),
    "diamond" => (["diamond", "version"], 2),
    "blast" => (["blastn", "-version"], 2),
    "minimap2" => (["minimap2", "--version"], 2),
    "samtools" => (["samtools", "--version"], 1),
    "fastp" => (["fastp", "--version"], 0),
    "prodigal" => (["prodigal", "-v"], 2),
    "bakta" => (["bakta", "--version"], 1)
)

"""
    verify_tool_version(env_name::AbstractString; expected_major::Union{Nothing, Integer} = nothing, throw_on_mismatch::Bool = true)

Check the installed version of a conda tool against an expected major version.
Returns a named tuple `(; version_string, major, minor, patch, ok)`.

When `expected_major` is `nothing`, looks up the tool in [`TOOL_VERSION_REGISTRY`](@ref).
If the tool is not in the registry and no `expected_major` is given, returns the
parsed version without comparison (`ok = true`).

Throws `ErrorException` on major-version mismatch when `throw_on_mismatch` is true;
otherwise logs a warning and returns `ok = false`.

# Examples
```julia
# Using registry (metaphlan expects major=4)
result = verify_tool_version("metaphlan")

# Explicit expected version
result = verify_tool_version("metaphlan"; expected_major = 4)

# Just parse, don't enforce
result = verify_tool_version("some_tool"; throw_on_mismatch = false)
```
"""
function verify_tool_version(env_name::AbstractString;
        expected_major::Union{Nothing, Integer} = nothing,
        throw_on_mismatch::Bool = true)
    # Resolve version check command and expected major from registry if needed
    if haskey(TOOL_VERSION_REGISTRY, env_name)
        registry_cmd, registry_major = TOOL_VERSION_REGISTRY[env_name]
        cmd_parts = registry_cmd
        if expected_major === nothing
            expected_major = registry_major
        end
    else
        # Unknown tool — try generic --version
        cmd_parts = [env_name, "--version"]
    end

    version_string = conda_tool_version(env_name, cmd_parts)
    if version_string === missing
        msg = "Could not determine version for '$env_name' (environment may not exist)"
        if throw_on_mismatch
            error(msg)
        else
            @warn msg
            return (; version_string = missing, major = nothing, minor = nothing,
                patch = nothing, ok = false)
        end
    end

    # Parse first version-like pattern (handles diverse formats like
    # "MetaPhlAn version 4.1.0", "blastn: 2.14.0+", "4.0.0rc1")
    m = match(r"(\d+)\.(\d+)(?:\.(\d+))?", version_string)
    if m === nothing
        msg = "Could not parse version from '$env_name' output: $version_string"
        if throw_on_mismatch
            error(msg)
        else
            @warn msg
            return (; version_string, major = nothing, minor = nothing,
                patch = nothing, ok = false)
        end
    end

    major = parse(Int, m.captures[1])
    minor = parse(Int, m.captures[2])
    patch = m.captures[3] !== nothing ? parse(Int, m.captures[3]) : nothing

    ok = true
    if expected_major !== nothing && major != expected_major
        ok = false
        msg = "Version mismatch for '$env_name': installed $version_string " *
              "(major=$major), expected major=$expected_major"
        if throw_on_mismatch
            error(msg)
        else
            @warn msg
        end
    else
        @info "Verified '$env_name' version: $version_string (major=$major)"
    end

    return (; version_string, major, minor, patch, ok)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Check whether a named Bioconda environment already exists.

# Arguments
- `pkg::String`: Name of the environment.

# Returns
`Bool` indicating if the environment is present.
"""
function check_bioconda_env_is_installed(pkg)
    _ensure_conda_env_vars!()
    # ensure conda environment is available
    if !isfile(CONDA_RUNNER)
        if (basename(CONDA_RUNNER) == "mamba")
            Conda.add("mamba")
        elseif (basename(CONDA_RUNNER) == "conda")
            Conda.update()
        end
    end
    # try
    current_environments = Set(first.(filter(x -> length(x) == 2,
        split.(filter(x -> !occursin(r"^#", x), readlines(`$(CONDA_RUNNER) env list`))))))
    if pkg in current_environments
        return true
    else
        return false
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a Conda environment from a YAML file.

# Arguments
- `yaml_file::AbstractString`: Path to the environment.yml file
- `env_name::AbstractString`: Name of the environment to create

# Keywords
- `force::Bool=false`: If true, remove existing environment before creation
"""
function create_conda_env_from_yaml(yaml_file::AbstractString, env_name::AbstractString; force::Bool = false)
    _ensure_conda_env_vars!()
    if !isfile(yaml_file)
        throw(ArgumentError("YAML file does not exist: $(yaml_file)"))
    end

    if check_bioconda_env_is_installed(env_name)
        if force
            @info "Removing existing environment '$env_name'..."
            run(`$(CONDA_RUNNER) env remove -n $(env_name) -y`)
        else
            @info "Environment '$env_name' already exists. Use force=true to recreate."
            return env_name
        end
    end

    @info "Creating environment '$env_name' from $(yaml_file)..."
    run(`$(CONDA_RUNNER) env create -f $(yaml_file) -n $(env_name)`)
    run(`$(CONDA_RUNNER) clean --all -y`)

    return env_name
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a new Conda environment with a specified Bioconda package.

# Arguments
- `pkg::String`: Package name to install. Can include channel specification using 
the format "channel::package"

# Keywords
- `force::Bool=false`: If true, recreates the environment even if it already exists
- `quiet::Bool=false`: If true, suppress conda output using --quiet flag

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
function add_bioconda_env(pkg; force = false, quiet = false)
    _ensure_conda_env_vars!()
    channel = nothing
    if occursin("::", pkg)
        if !quiet
            println("splitting $(pkg)")
        end
        channel, pkg = split(pkg, "::")
        if !quiet
            println("into channel:$(channel) pkg:$(pkg)")
        end
    end
    already_installed = check_bioconda_env_is_installed(pkg)
    if !already_installed || force
        if !quiet
            @info "installing conda environment $(pkg)"
        end
        if isnothing(channel)
            if quiet
                run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y --quiet`)
            else
                run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(pkg) -y`)
            end
        else
            if quiet
                run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(channel)::$(pkg) -y --quiet`)
            else
                run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda -c defaults --strict-channel-priority -n $(pkg) $(channel)::$(pkg) -y`)
            end
        end
        if quiet
            run(`$(CONDA_RUNNER) clean --all -y --quiet`)
        else
            run(`$(CONDA_RUNNER) clean --all -y`)
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Update a package and its dependencies in its dedicated Conda environment.

# Arguments
- `pkg::String`: Name of the package/environment to update
"""
function update_bioconda_env(pkg)
    _ensure_conda_env_vars!()
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
#             "augustus",
#             "megahit",
#             "medaka",
#             "metaeuk",
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
#             "prodigal-gv",
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
