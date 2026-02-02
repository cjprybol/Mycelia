"""
Wrapper for ggcat (Graph Construction and Querying).
"""

const GGCAT_ENV_NAME = "ggcat"

function _ggcat_input_args(input_files::Union{
        AbstractString, AbstractVector{<:AbstractString}})
    if input_files isa AbstractVector{<:AbstractString}
        if isempty(input_files)
            throw(ArgumentError("input_files cannot be empty"))
        end
        for file in input_files
            if !isfile(file)
                throw(ArgumentError("Input file not found: $(file)"))
            end
        end
        return String.(input_files)
    end

    if !isfile(input_files)
        throw(ArgumentError("Input file not found: $(input_files)"))
    end
    return [String(input_files)]
end

function _ggcat_resolve_output_path(output_file::AbstractString)
    if isfile(output_file)
        return String(output_file)
    end

    alternates = String[]
    if !endswith(output_file, ".lz4")
        push!(alternates, "$(output_file).lz4")
    end
    if !endswith(output_file, ".gz")
        push!(alternates, "$(output_file).gz")
    end

    for candidate in alternates
        if isfile(candidate)
            return candidate
        end
    end

    return String(output_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Ensure the GGCAT environment is installed via Bioconda.
"""
function install_ggcat(; force::Bool = false)
    Mycelia.add_bioconda_env(GGCAT_ENV_NAME; force = force)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run `ggcat build` to construct a compacted or colored de Bruijn graph.

# Arguments
- `input_files::Union{String, Vector{String}}`: Path(s) to input FASTA/FASTQ files.
- `output_file::String`: Path for the output graph.
- `kmer_length::Int`: Length of k-mers (required).

# Keywords
- `threads::Int=1`: Number of threads to use.
- `min_multiplicity::Int=2`: Minimum k-mer abundance required to be included in the graph.
- `colors::Bool=false`: Whether to build a colored graph.
- `generate_links::Bool=false`: Generate maximal unitig connections (equivalent to `-e`).
- `matchtigs::Bool=false`: Compute greedy matchtigs (`-g`).
- `eulertigs::Bool=false`: Compute eulertigs (`--eulertigs`).
- `pathtigs::Bool=false`: Compute pathtigs (`--pathtigs`).
- `prefer_memory::Bool=false`: Use all given memory before writing to disk (`--prefer-memory`).
- `color_mapping::Union{Nothing,AbstractString}=nothing`: Path to a colored input list (`--colored-input-lists`).

# Returns
The path to the output graph file (resolved if ggcat appends `.lz4` or `.gz`).
"""
function ggcat_build(
        input_files::Union{AbstractString, AbstractVector{<:AbstractString}},
        output_file::AbstractString,
        kmer_length::Int;
        threads::Int = 1,
        min_multiplicity::Int = 2,
        colors::Bool = false,
        generate_links::Bool = false,
        matchtigs::Bool = false,
        eulertigs::Bool = false,
        pathtigs::Bool = false,
        prefer_memory::Bool = false,
        color_mapping::Union{Nothing, AbstractString} = nothing
)
    if kmer_length < 1
        throw(ArgumentError("kmer_length must be positive, got $(kmer_length)"))
    end
    if threads < 1
        throw(ArgumentError("threads must be positive, got $(threads)"))
    end
    if min_multiplicity < 1
        throw(ArgumentError("min_multiplicity must be positive, got $(min_multiplicity)"))
    end

    install_ggcat()
    mkpath(dirname(String(output_file)))

    input_args = _ggcat_input_args(input_files)
    cmd_flags = String[
    "build",
    "--kmer-length", string(kmer_length),
    "--threads-count", string(threads),
    "--min-multiplicity", string(min_multiplicity),
    "--output-file", String(output_file)
]

    if colors
        push!(cmd_flags, "--colors")
    end

    if !isnothing(color_mapping)
        if !isfile(color_mapping)
            throw(ArgumentError("color_mapping file not found: $(color_mapping)"))
        end
        push!(cmd_flags, "--colored-input-lists", String(color_mapping))
    end

    if generate_links
        push!(cmd_flags, "-e")
    end

    if matchtigs
        push!(cmd_flags, "-g")
    end

    if eulertigs
        push!(cmd_flags, "--eulertigs")
    end

    if pathtigs
        push!(cmd_flags, "--pathtigs")
    end

    if prefer_memory
        push!(cmd_flags, "--prefer-memory")
    end

    append!(cmd_flags, input_args)

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n $(GGCAT_ENV_NAME) ggcat $(cmd_flags)`
    run(cmd)

    resolved_output = _ggcat_resolve_output_path(output_file)
    if !isfile(resolved_output)
        error("ggcat build failed to produce output at $(output_file)")
    end

    return resolved_output
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run `ggcat query` to query k-mers in an existing graph.

# Arguments
- `graph_file::String`: Path to the existing graph file (output of `ggcat build`).
- `query_file::String`: Path to the FASTA/FASTQ file containing query sequences.
- `output_file::String`: Path where the query results will be saved.
- `kmer_length::Int`: k-mer length used during graph construction.

# Keywords
- `threads::Int=1`: Number of threads.
- `colors::Bool=false`: Set to true if querying a colored graph.
"""
function ggcat_query(
        graph_file::AbstractString,
        query_file::AbstractString,
        output_file::AbstractString,
        kmer_length::Int;
        threads::Int = 1,
        colors::Bool = false
)
    if kmer_length < 1
        throw(ArgumentError("kmer_length must be positive, got $(kmer_length)"))
    end
    if threads < 1
        throw(ArgumentError("threads must be positive, got $(threads)"))
    end
    if !isfile(graph_file)
        throw(ArgumentError("graph_file not found: $(graph_file)"))
    end
    if !isfile(query_file)
        throw(ArgumentError("query_file not found: $(query_file)"))
    end

    install_ggcat()
    mkpath(dirname(String(output_file)))

    cmd_flags = String[
    "query",
    "-k", string(kmer_length),
    "-j", string(threads)
]

    if colors
        push!(cmd_flags, "--colors")
    end

    append!(cmd_flags, [String(graph_file), String(query_file)])

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n $(GGCAT_ENV_NAME) ggcat $(cmd_flags)`
    run(pipeline(cmd, stdout = String(output_file)))

    if !isfile(output_file)
        error("ggcat query failed to produce output at $(output_file)")
    end

    return String(output_file)
end
