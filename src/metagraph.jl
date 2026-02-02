"""
Wrapper for MetaGraph (ratschlab/metagraph).
Reference: https://github.com/ratschlab/metagraph
"""

const METAGRAPH_CONDA_ENV = "metagraph"

function _metagraph_require_path(path::AbstractString; label::AbstractString = "path")
    path_str = String(path)
    if !ispath(path_str)
        throw(ArgumentError("$(label) not found: $(path_str)"))
    end
    return path_str
end

function _metagraph_input_args(inputs::AbstractVector{<:AbstractString}; label::AbstractString = "input")
    isempty(inputs) && throw(ArgumentError("$(label) list cannot be empty"))
    return [_metagraph_require_path(input_path; label = label) for input_path in inputs]
end

function _metagraph_output_base(outfile_base::AbstractString, out_dir::Union{
        Nothing, AbstractString})
    base = String(outfile_base)
    if isnothing(out_dir)
        output_dir = dirname(base)
        if output_dir != "." && !isempty(output_dir)
            mkpath(output_dir)
        end
        return base
    end
    out_dir_str = String(out_dir)
    mkpath(out_dir_str)
    return joinpath(out_dir_str, base)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get the MetaGraph command, optionally specifying the executable name (e.g., "metagraph", "metagraph_Protein", "metagraph_DNA5").
"""
function metagraph_cmd(args::Vector{String};
        env::AbstractString = METAGRAPH_CONDA_ENV,
        executable::AbstractString = "metagraph",
        live_stream::Bool = true)
    if live_stream
        return `$(CONDA_RUNNER) run --live-stream -n $(env) $(executable) $(args)`
    end
    return `$(CONDA_RUNNER) run -n $(env) $(executable) $(args)`
end

metagraph_cmd(args...; kwargs...) = metagraph_cmd(string.(collect(args)); kwargs...)

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install MetaGraph via Bioconda.
"""
function install_metagraph(; force::Bool = false, quiet::Bool = false)
    add_bioconda_env(METAGRAPH_CONDA_ENV; force = force, quiet = quiet)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run a generic MetaGraph command.

# Keywords
- `env`: Conda environment name (default: "metagraph")
- `executable`: Executable name (default: "metagraph", use "metagraph_Protein" for amino acids)
- `force_env`: Recreate the conda environment if true
- `quiet_env`: Suppress conda output during environment creation
- `mmap`: Enable memory-mapped loading for large graphs/annotations
- `workdir`: Working directory for execution
"""
function run_metagraph(args::Vector{String};
        env::AbstractString = METAGRAPH_CONDA_ENV,
        executable::AbstractString = "metagraph",
        live_stream::Bool = true,
        force_env::Bool = false,
        quiet_env::Bool = false,
        mmap::Bool = false,
        workdir::Union{Nothing, AbstractString} = nothing)
    add_bioconda_env(METAGRAPH_CONDA_ENV; force = force_env, quiet = quiet_env)
    args_final = args
    if mmap && !("--mmap" in args)
        args_final = copy(args)
        if isempty(args_final)
            push!(args_final, "--mmap")
        else
            insert!(args_final, 2, "--mmap")
        end
    end
    cmd = metagraph_cmd(args_final; env = env, executable = executable, live_stream = live_stream)
    if isnothing(workdir)
        run(cmd)
    else
        run(Cmd(cmd; dir = workdir))
    end
    return cmd
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a de Bruijn graph from sequences.

# Arguments
- `inputs`: Vector of input file paths (FASTA/FASTQ/KMC)

# Keywords
- `k`: k-mer size (default: 31)
- `outfile_base`: Output filename base (default: "graph")
- `out_dir`: Output directory
- `threads`: Number of parallel threads
- `mem_cap_gb`: Memory cap in GB
- `disk_swap`: Path for disk swap directory
- `mode`: Graph mode ("canonical", "primary", etc.)
- `count_kmers`: Build a counting de Bruijn graph
- `verbose`: Enable verbose output
- `executable`: Binary to use ("metagraph", "metagraph_Protein", etc.)
"""
function metagraph_build(inputs::AbstractVector{<:AbstractString};
        k::Integer = 31,
        outfile_base::String = "graph",
        out_dir::Union{Nothing, AbstractString} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        mem_cap_gb::Union{Nothing, Integer} = nothing,
        disk_swap::Union{Nothing, AbstractString} = nothing,
        mode::Union{Nothing, AbstractString} = nothing,
        count_kmers::Bool = false,
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    k > 0 || throw(ArgumentError("k must be positive, got $(k)"))
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end
    if !isnothing(mem_cap_gb)
        mem_cap_gb > 0 ||
            throw(ArgumentError("mem_cap_gb must be positive, got $(mem_cap_gb)"))
    end
    input_args = _metagraph_input_args(inputs; label = "input")

    args = String["build"]
    if verbose
        push!(args, "-v")
    end
    push!(args, "-k", string(k))
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end
    if !isnothing(mem_cap_gb)
        push!(args, "--mem-cap-gb", string(mem_cap_gb))
    end
    if !isnothing(disk_swap)
        push!(args, "--disk-swap", String(disk_swap))
    end
    if !isnothing(mode)
        push!(args, "--mode", String(mode))
    end
    if count_kmers
        push!(args, "--count-kmers")
    end

    final_output = _metagraph_output_base(outfile_base, out_dir)
    push!(args, "-o", final_output)

    append!(args, additional_args)
    append!(args, input_args)

    run_metagraph(args; executable = executable, kwargs...)
    return final_output
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Annotate a graph with labels.

# Arguments
- `graph`: Path to the input graph (.dbg)
- `annotations`: Vector of annotation file paths (FASTA etc.)

# Keywords
- `out_dir`: Output directory
- `outfile_base`: Output base name
- `anno_type`: Annotation type (e.g., "row", "column", "brwt")
- `count_kmers`: Count k-mers during annotation
- `count_width`: Bit width for count storage (e.g., 8, 16)
- `coordinates`: Store k-mer coordinates
- `anno_header`: Use sequence headers as labels
- `anno_filename`: Use input filenames as labels
- `anno_label`: Assign a custom label to the input file
- `separately`: Annotate each input file independently
- `threads_each`: Threads per annotation job when `separately=true`
- `fasta_anno`: Input is FASTA annotation
- `mem_cap_gb`: Memory cap in GB
- `disk_swap`: Path for disk swap directory
"""
function metagraph_annotate(graph::AbstractString,
        annotations::AbstractVector{<:AbstractString};
        out_dir::Union{Nothing, AbstractString} = nothing,
        outfile_base::String = "annotation",
        anno_type::Union{Nothing, AbstractString} = nothing,
        count_kmers::Bool = false,
        count_width::Union{Nothing, Integer} = nothing,
        coordinates::Bool = false,
        anno_header::Bool = false,
        anno_filename::Bool = false,
        anno_label::Union{Nothing, AbstractString} = nothing,
        separately::Bool = false,
        threads_each::Union{Nothing, Integer} = nothing,
        fasta_anno::Bool = false,
        threads::Union{Nothing, Integer} = get_default_threads(),
        mem_cap_gb::Union{Nothing, Integer} = nothing,
        disk_swap::Union{Nothing, AbstractString} = nothing,
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    graph_path = _metagraph_require_path(graph; label = "graph")
    anno_args = _metagraph_input_args(annotations; label = "annotation")
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end
    if !isnothing(threads_each)
        threads_each > 0 ||
            throw(ArgumentError("threads_each must be positive, got $(threads_each)"))
    end
    if !isnothing(count_width)
        count_width > 0 ||
            throw(ArgumentError("count_width must be positive, got $(count_width)"))
    end
    if !isnothing(mem_cap_gb)
        mem_cap_gb > 0 ||
            throw(ArgumentError("mem_cap_gb must be positive, got $(mem_cap_gb)"))
    end

    args = String["annotate"]
    if verbose
        push!(args, "-v")
    end
    push!(args, "-i", graph_path)
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end
    if !isnothing(mem_cap_gb)
        push!(args, "--mem-cap-gb", string(mem_cap_gb))
    end
    if !isnothing(disk_swap)
        push!(args, "--disk-swap", String(disk_swap))
    end
    if !isnothing(anno_type)
        push!(args, "--anno-type", String(anno_type))
    end
    if count_kmers
        push!(args, "--count-kmers")
    end
    if !isnothing(count_width)
        push!(args, "--count-width", string(count_width))
    end
    if coordinates
        push!(args, "--coordinates")
    end
    if anno_header
        push!(args, "--anno-header")
    end
    if anno_filename
        push!(args, "--anno-filename")
    end
    if !isnothing(anno_label)
        push!(args, "--anno-label", String(anno_label))
    end
    if separately
        push!(args, "--separately")
    end
    if !isnothing(threads_each)
        push!(args, "--threads-each", string(threads_each))
    end
    if fasta_anno
        push!(args, "--fasta-anno")
    end

    final_output = _metagraph_output_base(outfile_base, out_dir)
    push!(args, "-o", final_output)

    append!(args, additional_args)
    append!(args, anno_args)

    run_metagraph(args; executable = executable, kwargs...)
    return final_output
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Query the graph for sequences.

# Arguments
- `graph`: Input graph (.dbg)
- `queries`: Query sequences file (FASTA/FASTQ)

# Keywords
- `annotation`: Annotation file (.annodbg)
- `min_kmers_fraction`: Minimum fraction of k-mers required (0.0-1.0)
- `labels_delimiter`: Delimiter for labels in output
- `query_mode`: Query mode (e.g., "counts" or "coords")
"""
function metagraph_query(graph::AbstractString,
        queries::AbstractString;
        annotation::Union{Nothing, AbstractString} = nothing,
        min_kmers_fraction::Union{Nothing, Float64} = nothing,
        labels_delimiter::Union{Nothing, AbstractString} = nothing,
        query_mode::Union{Nothing, AbstractString} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    graph_path = _metagraph_require_path(graph; label = "graph")
    query_path = _metagraph_require_path(queries; label = "query")
    if !isnothing(min_kmers_fraction)
        if min_kmers_fraction < 0.0 || min_kmers_fraction > 1.0
            throw(ArgumentError("min_kmers_fraction must be between 0.0 and 1.0, got $(min_kmers_fraction)"))
        end
    end
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end

    args = String["query"]
    if verbose
        push!(args, "-v")
    end
    push!(args, "-i", graph_path)
    if !isnothing(annotation)
        push!(args, "-a", _metagraph_require_path(annotation; label = "annotation"))
    end
    if !isnothing(min_kmers_fraction)
        push!(args, "--min-kmers-fraction-label", string(min_kmers_fraction))
    end
    if !isnothing(labels_delimiter)
        push!(args, "--labels-delimiter", String(labels_delimiter))
    end
    if !isnothing(query_mode)
        push!(args, "--query-mode", String(query_mode))
    end
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end

    append!(args, additional_args)
    push!(args, query_path)

    run_metagraph(args; executable = executable, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Align sequences to the graph.

# Arguments
- `graph`: Input graph (.dbg)
- `queries`: Query sequences file

# Keywords
- `seeder`: Seeding method (e.g., "sketch")
- `batch_size`: Batch size for alignment
- `map`: Use pseudo-alignment with exact k-mer matches
- `query_presence`: Output presence/absence per sequence
- `count_kmers`: Output k-mer match counts per sequence
- `align_min_seed_length`: Minimum seed length for alignment
- `annotation`: Annotation file (.annodbg) to enable label-aware alignment
- `min_kmers_fraction`: Minimum fraction of k-mers required (0.0-1.0)
"""
function metagraph_align(graph::AbstractString,
        queries::AbstractString;
        seeder::Union{Nothing, AbstractString} = nothing,
        batch_size::Union{Nothing, Integer} = nothing,
        map::Bool = false,
        query_presence::Bool = false,
        count_kmers::Bool = false,
        align_min_seed_length::Union{Nothing, Integer} = nothing,
        annotation::Union{Nothing, AbstractString} = nothing,
        min_kmers_fraction::Union{Nothing, Float64} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    graph_path = _metagraph_require_path(graph; label = "graph")
    query_path = _metagraph_require_path(queries; label = "query")
    if !isnothing(batch_size)
        batch_size > 0 ||
            throw(ArgumentError("batch_size must be positive, got $(batch_size)"))
    end
    if !isnothing(align_min_seed_length)
        align_min_seed_length > 0 ||
            throw(ArgumentError("align_min_seed_length must be positive, got $(align_min_seed_length)"))
    end
    if !isnothing(min_kmers_fraction)
        if min_kmers_fraction < 0.0 || min_kmers_fraction > 1.0
            throw(ArgumentError("min_kmers_fraction must be between 0.0 and 1.0, got $(min_kmers_fraction)"))
        end
    end
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end

    args = String["align"]
    if verbose
        push!(args, "-v")
    end
    push!(args, "-i", graph_path)
    if !isnothing(annotation)
        push!(args, "-a", _metagraph_require_path(annotation; label = "annotation"))
    end
    if !isnothing(seeder)
        push!(args, "--seeder", String(seeder))
    end
    if !isnothing(batch_size)
        push!(args, "--batch-size", string(batch_size))
    end
    if map
        push!(args, "--map")
    end
    if query_presence
        push!(args, "--query-presence")
    end
    if count_kmers
        push!(args, "--count-kmers")
    end
    if !isnothing(align_min_seed_length)
        push!(args, "--align-min-seed-length", string(align_min_seed_length))
    end
    if !isnothing(min_kmers_fraction)
        push!(args, "--min-kmers-fraction-label", string(min_kmers_fraction))
    end
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end

    append!(args, additional_args)
    push!(args, query_path)

    run_metagraph(args; executable = executable, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Transform a graph representation.

# Keywords
- `out_file`: Output file path
- `to_fasta`: Emit contigs as FASTA
- `primary_kmers`: Extract primary contigs from canonical graphs
- `unitigs`: Emit unitigs instead of contigs
- `state`: Transform to specific state (e.g., "small" for compressed)
- `threads`: Number of parallel threads
"""
function metagraph_transform(input_graph::AbstractString;
        out_file::AbstractString,
        to_fasta::Bool = false,
        primary_kmers::Bool = false,
        unitigs::Bool = false,
        state::Union{Nothing, AbstractString} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    graph_path = _metagraph_require_path(input_graph; label = "graph")
    output_path = _metagraph_output_base(out_file, nothing)
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end

    args = String["transform"]
    if verbose
        push!(args, "-v")
    end
    if to_fasta
        push!(args, "--to-fasta")
    end
    if primary_kmers
        push!(args, "--primary-kmers")
    end
    if unitigs
        push!(args, "--unitigs")
    end
    if !isnothing(state)
        push!(args, "--state", String(state))
    end
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end
    push!(args, "-o", output_path)

    append!(args, additional_args)
    push!(args, graph_path)

    run_metagraph(args; executable = executable, kwargs...)
    return output_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Transform graph annotations (e.g., convert to Multi-BRWT).

# Arguments
- `input_annos`: Input annotation file or list of annotation files

# Keywords
- `out_file`: Output filename
- `anno_type`: Target annotation type
- `linkage`: Compute linkage
- `greedy`: Use greedy algorithm
- `subsample`: Subsample size
- `linkage_file`: Path to linkage file
- `row_diff_stage`: Stage for RowDiff transformation (0-2)
- `graph`: Graph path required for diff transforms
- `coordinates`: Transform coordinate-aware annotations
- `count_kmers`: Transform count-aware annotations
- `aggregate_columns`: Aggregate columns into a new annotation
- `min_count`: Minimum number of columns matched
- `max_count`: Maximum number of columns matched
- `min_fraction`: Minimum fraction of columns matched (0.0-1.0)
- `max_fraction`: Maximum fraction of columns matched (0.0-1.0)
- `min_value`: Minimum count value threshold
- `max_value`: Maximum count value threshold
- `rename_cols`: Rename columns using a rules file
- `parallel_nodes`: Concurrent nodes merged when building BRWT
- `mem_cap_gb`: Memory cap in GB
- `disk_swap`: Path for disk swap directory
"""
function metagraph_transform_anno(
        input_annos::Union{AbstractString, AbstractVector{<:AbstractString}};
        out_file::AbstractString,
        anno_type::Union{Nothing, AbstractString} = nothing,
        linkage::Bool = false,
        greedy::Bool = false,
        subsample::Union{Nothing, Integer} = nothing,
        linkage_file::Union{Nothing, AbstractString} = nothing,
        row_diff_stage::Union{Nothing, Integer} = nothing,
        graph::Union{Nothing, AbstractString} = nothing,
        coordinates::Bool = false,
        count_kmers::Bool = false,
        aggregate_columns::Bool = false,
        min_count::Union{Nothing, Integer} = nothing,
        max_count::Union{Nothing, Integer} = nothing,
        min_fraction::Union{Nothing, Float64} = nothing,
        max_fraction::Union{Nothing, Float64} = nothing,
        min_value::Union{Nothing, Integer} = nothing,
        max_value::Union{Nothing, Integer} = nothing,
        rename_cols::Union{Nothing, AbstractString} = nothing,
        parallel_nodes::Union{Nothing, Integer} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        mem_cap_gb::Union{Nothing, Integer} = nothing,
        disk_swap::Union{Nothing, AbstractString} = nothing,
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    input_args = if input_annos isa AbstractVector{<:AbstractString}
        _metagraph_input_args(input_annos; label = "annotation")
    else
        [_metagraph_require_path(input_annos; label = "annotation")]
    end
    output_path = _metagraph_output_base(out_file, nothing)
    if !isnothing(subsample)
        subsample > 0 ||
            throw(ArgumentError("subsample must be positive, got $(subsample)"))
    end
    if !isnothing(row_diff_stage)
        row_diff_stage >= 0 ||
            throw(ArgumentError("row_diff_stage must be non-negative, got $(row_diff_stage)"))
    end
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end
    if !isnothing(parallel_nodes)
        parallel_nodes > 0 ||
            throw(ArgumentError("parallel_nodes must be positive, got $(parallel_nodes)"))
    end
    if !isnothing(mem_cap_gb)
        mem_cap_gb > 0 ||
            throw(ArgumentError("mem_cap_gb must be positive, got $(mem_cap_gb)"))
    end
    if !isnothing(min_fraction)
        if min_fraction < 0.0 || min_fraction > 1.0
            throw(ArgumentError("min_fraction must be between 0.0 and 1.0, got $(min_fraction)"))
        end
    end
    if !isnothing(max_fraction)
        if max_fraction < 0.0 || max_fraction > 1.0
            throw(ArgumentError("max_fraction must be between 0.0 and 1.0, got $(max_fraction)"))
        end
    end
    if !isnothing(min_count)
        min_count >= 0 ||
            throw(ArgumentError("min_count must be non-negative, got $(min_count)"))
    end
    if !isnothing(max_count)
        max_count >= 0 ||
            throw(ArgumentError("max_count must be non-negative, got $(max_count)"))
    end
    if !isnothing(min_value)
        min_value >= 0 ||
            throw(ArgumentError("min_value must be non-negative, got $(min_value)"))
    end
    if !isnothing(max_value)
        max_value >= 0 ||
            throw(ArgumentError("max_value must be non-negative, got $(max_value)"))
    end

    args = String["transform_anno"]
    if verbose
        push!(args, "-v")
    end
    if aggregate_columns
        push!(args, "--aggregate-columns")
    end
    push!(args, "-o", output_path)
    if !isnothing(anno_type)
        push!(args, "--anno-type", String(anno_type))
    end
    if linkage
        push!(args, "--linkage")
    end
    if greedy
        push!(args, "--greedy")
    end
    if !isnothing(subsample)
        push!(args, "--subsample", string(subsample))
    end
    if !isnothing(linkage_file)
        push!(args, "--linkage-file", _metagraph_require_path(linkage_file; label = "linkage file"))
    end
    if !isnothing(row_diff_stage)
        push!(args, "--row-diff-stage", string(row_diff_stage))
    end
    if !isnothing(graph)
        push!(args, "-i", _metagraph_require_path(graph; label = "graph"))
    end
    if coordinates
        push!(args, "--coordinates")
    end
    if count_kmers
        push!(args, "--count-kmers")
    end
    if !isnothing(rename_cols)
        push!(args, "--rename-cols", _metagraph_require_path(rename_cols; label = "rename rules"))
    end
    if !isnothing(min_count)
        push!(args, "--min-count", string(min_count))
    end
    if !isnothing(max_count)
        push!(args, "--max-count", string(max_count))
    end
    if !isnothing(min_fraction)
        push!(args, "--min-fraction", string(min_fraction))
    end
    if !isnothing(max_fraction)
        push!(args, "--max-fraction", string(max_fraction))
    end
    if !isnothing(min_value)
        push!(args, "--min-value", string(min_value))
    end
    if !isnothing(max_value)
        push!(args, "--max-value", string(max_value))
    end
    if !isnothing(parallel_nodes)
        push!(args, "--parallel-nodes", string(parallel_nodes))
    end
    if !isnothing(mem_cap_gb)
        push!(args, "--mem-cap-gb", string(mem_cap_gb))
    end
    if !isnothing(disk_swap)
        push!(args, "--disk-swap", String(disk_swap))
    end
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end

    append!(args, additional_args)
    append!(args, input_args)

    run_metagraph(args; executable = executable, kwargs...)
    return output_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Relax a Multi-BRWT annotation to improve compression.
"""
function metagraph_relax_brwt(input_anno::AbstractString;
        out_file::AbstractString,
        relax_arity::Union{Nothing, Integer} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    input_path = _metagraph_require_path(input_anno; label = "annotation")
    output_path = _metagraph_output_base(out_file, nothing)
    if !isnothing(relax_arity)
        relax_arity > 0 ||
            throw(ArgumentError("relax_arity must be positive, got $(relax_arity)"))
    end
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end

    args = String["relax_brwt"]
    if verbose
        push!(args, "-v")
    end
    push!(args, "-o", output_path)
    if !isnothing(relax_arity)
        push!(args, "--relax-arity", string(relax_arity))
    end
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end

    append!(args, additional_args)
    push!(args, input_path)

    run_metagraph(args; executable = executable, kwargs...)
    return output_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run MetaGraph in server mode for querying via the Python API or HTTP.
"""
function metagraph_server_query(graph::AbstractString;
        annotation::Union{Nothing, AbstractString} = nothing,
        port::Union{Nothing, Integer} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    graph_path = _metagraph_require_path(graph; label = "graph")
    if !isnothing(port)
        port > 0 || throw(ArgumentError("port must be positive, got $(port)"))
    end
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end

    args = String["server_query"]
    if verbose
        push!(args, "-v")
    end
    push!(args, "-i", graph_path)
    if !isnothing(annotation)
        push!(args, "-a", _metagraph_require_path(annotation; label = "annotation"))
    end
    if !isnothing(port)
        push!(args, "--port", string(port))
    end
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end

    append!(args, additional_args)

    run_metagraph(args; executable = executable, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assemble sequences from the graph.

# Keywords
- `unitigs`: Assemble unitigs
- `diff_assembly_rules`: Path to JSON rules file
"""
function metagraph_assemble(graph::AbstractString;
        out_file::AbstractString,
        unitigs::Bool = false,
        annotation::Union{Nothing, AbstractString} = nothing,
        diff_assembly_rules::Union{Nothing, AbstractString} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    graph_path = _metagraph_require_path(graph; label = "graph")
    output_path = _metagraph_output_base(out_file, nothing)
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end

    args = String["assemble"]
    if verbose
        push!(args, "-v")
    end
    push!(args, "-o", output_path)
    if unitigs
        push!(args, "--unitigs")
    end
    if !isnothing(annotation)
        push!(args, "-a", _metagraph_require_path(annotation; label = "annotation"))
    end
    if !isnothing(diff_assembly_rules)
        push!(args, "--diff-assembly-rules",
            _metagraph_require_path(diff_assembly_rules; label = "diff assembly rules"))
    end
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end

    append!(args, additional_args)
    push!(args, graph_path)

    run_metagraph(args; executable = executable, kwargs...)
    return output_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get statistics for a graph or annotation.
"""
function metagraph_stats(input_file::AbstractString;
        annotation::Union{Nothing, AbstractString} = nothing,
        print_col_names::Bool = false,
        print_counts_hist::Bool = false,
        count_quantiles::Union{Nothing, String} = nothing,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    input_path = _metagraph_require_path(input_file; label = "input")

    args = String["stats"]
    if !isnothing(annotation)
        push!(args, "-a", _metagraph_require_path(annotation; label = "annotation"))
    end
    if print_col_names
        push!(args, "--print-col-names")
    end
    if print_counts_hist
        push!(args, "--print-counts-hist")
    end
    if !isnothing(count_quantiles)
        push!(args, "--count-quantiles", String(count_quantiles))
    end

    append!(args, additional_args)
    push!(args, input_path)

    run_metagraph(args; executable = executable, kwargs...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Clean the graph (remove sequencing errors).

# Keywords
- `to_fasta`: Emit cleaned contigs as FASTA
- `prune_tips`: Prune tips shorter than this threshold
- `prune_unitigs`: Prune unitigs below this count threshold
- `fallback`: Fallback abundance threshold when auto-estimation fails
"""
function metagraph_clean(graph::AbstractString;
        out_file::AbstractString,
        to_fasta::Bool = false,
        prune_tips::Union{Nothing, Integer} = nothing,
        prune_unitigs::Union{Nothing, Integer} = nothing,
        fallback::Union{Nothing, Integer} = nothing,
        threads::Union{Nothing, Integer} = get_default_threads(),
        verbose::Bool = true,
        executable::AbstractString = "metagraph",
        additional_args::Vector{String} = String[],
        kwargs...)
    graph_path = _metagraph_require_path(graph; label = "graph")
    output_path = _metagraph_output_base(out_file, nothing)
    if !isnothing(prune_tips)
        prune_tips >= 0 ||
            throw(ArgumentError("prune_tips must be non-negative, got $(prune_tips)"))
    end
    if !isnothing(prune_unitigs)
        prune_unitigs >= 0 ||
            throw(ArgumentError("prune_unitigs must be non-negative, got $(prune_unitigs)"))
    end
    if !isnothing(fallback)
        fallback >= 0 ||
            throw(ArgumentError("fallback must be non-negative, got $(fallback)"))
    end
    if !isnothing(threads)
        threads > 0 || throw(ArgumentError("threads must be positive, got $(threads)"))
    end

    args = String["clean"]
    if verbose
        push!(args, "-v")
    end
    if to_fasta
        push!(args, "--to-fasta")
    end
    if !isnothing(prune_tips)
        push!(args, "--prune-tips", string(prune_tips))
    end
    if !isnothing(prune_unitigs)
        push!(args, "--prune-unitigs", string(prune_unitigs))
    end
    if !isnothing(fallback)
        push!(args, "--fallback", string(fallback))
    end
    push!(args, "-o", output_path)
    if !isnothing(threads)
        push!(args, "--parallel", string(threads))
    end

    append!(args, additional_args)
    push!(args, graph_path)

    run_metagraph(args; executable = executable, kwargs...)
    return output_path
end
