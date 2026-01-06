"""
Wrapper for bcalm (De Bruijn Graph construction) and GFA conversion.
References:
- https://github.com/GATB/bcalm
- https://github.com/GATB/bcalm/blob/master/scripts/convertToGFA.py
"""

const BCALM_ENV_NAME = "bcalm"
const BCALM_GFA_SCRIPT_URL = "https://raw.githubusercontent.com/GATB/bcalm/master/scripts/convertToGFA.py"

function _bcalm_paths()
    install_dir = joinpath(dirname(dirname(pathof(Mycelia))), "deps", "bcalm")
    script_path = joinpath(install_dir, "convertToGFA.py")
    return install_dir, script_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install bcalm via Bioconda and download the GFA conversion script.

# Keywords
- `force::Bool=false`: Force reinstallation and re-download of dependencies.
"""
function install_bcalm(; force::Bool=false)
    Mycelia.add_bioconda_env(BCALM_ENV_NAME; force=force)
    try
        run(`$(Mycelia.CONDA_RUNNER) run -n $(BCALM_ENV_NAME) python -c "import sys"`)
    catch
        @info "Installing python for bcalm conversion script..."
        run(`$(Mycelia.CONDA_RUNNER) install -y -n $(BCALM_ENV_NAME) python`)
    end

    install_dir, script_path = _bcalm_paths()
    mkpath(install_dir)

    if !isfile(script_path) || force
        @info "Downloading convertToGFA.py..."
        Downloads.download(BCALM_GFA_SCRIPT_URL, script_path)
    end

    return nothing
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run bcalm to construct a compacted De Bruijn Graph and convert it to GFA.

# Arguments
- `input_files`: Single FASTQ/FASTA path or vector of paths.
- `out_dir`: Directory for outputs.

# Keywords
- `kmer_size::Int=21`: Size of k-mers.
- `abundance_min::Int=2`: Minimum k-mer abundance.
- `threads::Int=Sys.CPU_THREADS`: Number of cores to use.
- `tmp_dir::Union{Nothing, AbstractString}=nothing`: Optional directory for manifest files.
- `force::Bool=false`: Force reinstallation and re-download of dependencies.

# Returns
NamedTuple with keys `:unitigs` (path to .unitigs.fa) and `:gfa` (path to .gfa).
"""
function run_bcalm(
    input_files::Union{AbstractString, AbstractVector{<:AbstractString}},
    out_dir::AbstractString;
    kmer_size::Int=21,
    abundance_min::Int=2,
    threads::Int=Sys.CPU_THREADS,
    tmp_dir::Union{Nothing, AbstractString}=nothing,
    force::Bool=false,
)
    install_bcalm(force=force)
    mkpath(out_dir)

    input_arg = ""
    if input_files isa AbstractVector{<:AbstractString}
        if isempty(input_files)
            throw(ArgumentError("input_files cannot be empty"))
        end
        for file in input_files
            if !isfile(file)
                throw(ArgumentError("Input file not found: $(file)"))
            end
        end

        manifest_dir = isnothing(tmp_dir) ? out_dir : tmp_dir
        mkpath(manifest_dir)
        manifest_path = joinpath(manifest_dir, "bcalm_input_manifest.txt")
        open(manifest_path, "w") do io
            for file in input_files
                println(io, file)
            end
        end
        input_arg = manifest_path
    else
        if !isfile(input_files)
            throw(ArgumentError("Input file not found: $(input_files)"))
        end
        input_arg = input_files
    end

    prefix = joinpath(out_dir, "bcalm_assembly")

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n $(BCALM_ENV_NAME) bcalm -in $input_arg -kmer-size $kmer_size -abundance-min $abundance_min -out $prefix -nb-cores $threads`
    @info "Running bcalm..."
    run(cmd)

    unitigs_fa = "$(prefix).unitigs.fa"
    if !isfile(unitigs_fa)
        error("bcalm failed to produce output at $(unitigs_fa)")
    end

    _, script_path = _bcalm_paths()
    if !isfile(script_path)
        error("Missing convertToGFA.py at $(script_path). Run install_bcalm() to fetch it.")
    end

    out_gfa = "$(prefix).gfa"
    convert_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n $(BCALM_ENV_NAME) python $script_path $unitigs_fa $out_gfa $kmer_size`
    @info "Converting to GFA..."
    run(convert_cmd)

    if !isfile(out_gfa)
        error("convertToGFA.py failed to produce output at $(out_gfa)")
    end

    return (unitigs=unitigs_fa, gfa=out_gfa)
end
