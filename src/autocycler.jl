"""
Wrapper for Autocycler and the automated Autocycler bash pipeline.
Reference: https://github.com/rrwick/Autocycler
"""

const AUTOCYCLER_ENV_NAME = "autocycler"
const AUTOCYCLER_SCRIPT_URL = "https://raw.githubusercontent.com/rrwick/Autocycler/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh"
const AUTOCYCLER_ENV_URL = "https://raw.githubusercontent.com/rrwick/Autocycler/main/pipelines/Conda_environment_file_by_Ryan_Wick/environment.yml"

function _autocycler_paths()
    install_dir = joinpath(dirname(dirname(pathof(Mycelia))), "deps", "autocycler")
    script_path = joinpath(install_dir, "autocycler_full.sh")
    env_file_path = joinpath(install_dir, "environment.yml")
    return install_dir, script_path, env_file_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install Autocycler by creating the conda environment and downloading the automation script.

# Keywords
- `force::Bool=false`: Force recreation of the environment and re-download of the script.
"""
function install_autocycler(; force::Bool=false)
    install_dir, script_path, env_file_path = _autocycler_paths()
    mkpath(install_dir)

    if !isfile(env_file_path) || force
        @info "Downloading Autocycler environment.yml..."
        Downloads.download(AUTOCYCLER_ENV_URL, env_file_path)
    end

    create_conda_env_from_yaml(env_file_path, AUTOCYCLER_ENV_NAME; force=force)

    if !isfile(script_path) || force
        @info "Downloading autocycler_full.sh script..."
        Downloads.download(AUTOCYCLER_SCRIPT_URL, script_path)
        Base.chmod(script_path, 0o755)
    end

    @info "Autocycler installed successfully."
    return script_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the automated Autocycler pipeline.

# Keywords
- `long_reads::AbstractString`: Path to long reads (FASTQ)
- `out_dir::AbstractString`: Output directory
- `short_reads_1::Union{AbstractString, Nothing}=nothing`: Optional short reads R1 (FASTQ)
- `short_reads_2::Union{AbstractString, Nothing}=nothing`: Optional short reads R2 (FASTQ)
- `extra_args::Vector{String}=String[]`: Extra arguments passed to the script
"""
function run_autocycler(;
    long_reads::AbstractString,
    out_dir::AbstractString,
    short_reads_1::Union{AbstractString, Nothing}=nothing,
    short_reads_2::Union{AbstractString, Nothing}=nothing,
    extra_args::Vector{String}=String[]
)
    if !isfile(long_reads)
        throw(ArgumentError("Long-read FASTQ not found: $(long_reads)"))
    end

    if (isnothing(short_reads_1) && !isnothing(short_reads_2)) || (!isnothing(short_reads_1) && isnothing(short_reads_2))
        throw(ArgumentError("If providing short reads, both R1 and R2 must be provided."))
    end

    if !isnothing(short_reads_1) && !isfile(short_reads_1)
        throw(ArgumentError("Short-read R1 FASTQ not found: $(short_reads_1)"))
    end

    if !isnothing(short_reads_2) && !isfile(short_reads_2)
        throw(ArgumentError("Short-read R2 FASTQ not found: $(short_reads_2)"))
    end

    _, script_path, _ = _autocycler_paths()

    if !isfile(script_path) || !check_bioconda_env_is_installed(AUTOCYCLER_ENV_NAME)
        @warn "Autocycler not fully installed. Installing now..."
        install_autocycler()
    end

    cmd_args = String[script_path, long_reads, out_dir]

    if !isnothing(short_reads_1)
        push!(cmd_args, short_reads_1)
        push!(cmd_args, short_reads_2)
    end

    if !isempty(extra_args)
        append!(cmd_args, extra_args)
    end

    run_cmd = `$(CONDA_RUNNER) run --live-stream -n $(AUTOCYCLER_ENV_NAME) bash $(cmd_args)`

    @info "Starting Autocycler pipeline..."
    run(run_cmd)
    @info "Autocycler pipeline complete. Results in: $(out_dir)"

    return out_dir
end
