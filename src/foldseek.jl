"""
$(DocStringExtensions.TYPEDSIGNATURES)

Ensure the Foldseek environment is installed via Bioconda.
"""
function install_foldseek(; force=false)
    Mycelia.add_bioconda_env("foldseek"; force=force)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run `foldseek easy-search` to compare protein structures.

# Arguments
- `query_file::String`: Path to PDB/mmCIF file or directory of structures.
- `target_database::String`: Path to target database or PDB/mmCIF directory.
- `output_file::String`: Path to the output file (usually .m8 or .tsv).

# Keywords
- `tmp_dir::String=mktempdir()`: Temporary directory for Foldseek intermediates.
- `sensitivity::Float64=9.5`: Adjust sensitivity-speed trade-off. Lower is faster.
- `alignment_type::Int=2`:
  - 0: 3Di Gotoh-Smith-Waterman (local)
  - 1: TMalign (global, slow)
  - 2: 3Di+AA Gotoh-Smith-Waterman (local, default)
- `format_output::String`: Columns for output format.
- `format_mode::Int=0`: 0 for BLAST-tab, 3 for HTML, 4 for BLAST-tab custom.
- `threads::Int=1`: Number of CPU threads.
- `exhaustive_search::Bool=false`: Skip prefilter for all-vs-all alignment (slower).
- `gpu::Bool=false`: Enable GPU acceleration (requires CUDA-enabled system).
"""
function foldseek_easy_search(
    query_file::String,
    target_database::String,
    output_file::String;
    tmp_dir::String=mktempdir(),
    sensitivity::Float64=9.5,
    alignment_type::Int=2,
    format_output::String="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
    format_mode::Int=0,
    threads::Int=1,
    exhaustive_search::Bool=false,
    gpu::Bool=false,
)
    install_foldseek()

    cmd_flags = [
        "-s", string(sensitivity),
        "--alignment-type", string(alignment_type),
        "--format-output", format_output,
        "--format-mode", string(format_mode),
        "--threads", string(threads),
    ]

    if exhaustive_search
        push!(cmd_flags, "--exhaustive-search", "1")
    end

    if gpu
        push!(cmd_flags, "--gpu", "1")
    end

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n foldseek foldseek easy-search $query_file $target_database $output_file $tmp_dir $cmd_flags`
    run(cmd)

    if startswith(tmp_dir, tempdir()) && isdir(tmp_dir)
        rm(tmp_dir; recursive=true, force=true)
    end

    return output_file
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run `foldseek easy-cluster` to cluster protein structures.

# Arguments
- `input_file::String`: Path to input PDB/mmCIF files or formatted database.
- `output_prefix::String`: Prefix for output files (e.g., resulting in `prefix_cluster.tsv`).

# Keywords
- `tmp_dir::String=mktempdir()`: Temporary directory for Foldseek intermediates.
- `min_seq_id::Float64=0.0`: Minimum sequence identity for clustering.
- `coverage::Float64=0.8`: Fraction of aligned residues required.
- `cov_mode::Int=0`:
  - 0: Coverage of query and target
  - 1: Coverage of target
  - 2: Coverage of query
- `alignment_type::Int=2`: 0 (3Di), 1 (TMalign), 2 (3Di+AA).
- `sensitivity::Float64=9.5`: Adjust sensitivity-speed trade-off. Lower is faster.
- `threads::Int=1`: Number of CPU threads.
"""
function foldseek_easy_cluster(
    input_file::String,
    output_prefix::String;
    tmp_dir::String=mktempdir(),
    min_seq_id::Float64=0.0,
    coverage::Float64=0.8,
    cov_mode::Int=0,
    alignment_type::Int=2,
    sensitivity::Float64=9.5,
    threads::Int=1,
)
    install_foldseek()

    cmd_flags = [
        "--min-seq-id", string(min_seq_id),
        "-c", string(coverage),
        "--cov-mode", string(cov_mode),
        "--alignment-type", string(alignment_type),
        "-s", string(sensitivity),
        "--threads", string(threads),
    ]

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n foldseek foldseek easy-cluster $input_file $output_prefix $tmp_dir $cmd_flags`
    run(cmd)

    if startswith(tmp_dir, tempdir()) && isdir(tmp_dir)
        rm(tmp_dir; recursive=true, force=true)
    end

    return "$(output_prefix)_cluster.tsv"
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert PDB/mmCIF files into a Foldseek database.
"""
function foldseek_createdb(
    input_files::String,
    output_db_path::String;
    threads::Int=1,
)
    install_foldseek()

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n foldseek foldseek createdb $input_files $output_db_path --threads $threads`
    run(cmd)

    return output_db_path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Download pre-assembled Foldseek databases (e.g., "PDB", "Alphafold-SwissProt").

# Common Databases
- "PDB"
- "Alphafold-SwissProt"
- "Alphafold-Proteome"
- "GMGCL"
"""
function foldseek_databases(
    database_name::String,
    destination_dir::String;
    tmp_dir::String=mktempdir(),
    threads::Int=1,
)
    install_foldseek()

    if !isdir(destination_dir)
        mkpath(destination_dir)
    end

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n foldseek foldseek databases $database_name $destination_dir $tmp_dir --threads $threads`
    run(cmd)

    if startswith(tmp_dir, tempdir()) && isdir(tmp_dir)
        rm(tmp_dir; recursive=true, force=true)
    end

    return destination_dir
end
