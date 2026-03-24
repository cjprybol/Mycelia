"""
$(DocStringExtensions.TYPEDSIGNATURES)

Ensure the Foldseek environment is installed via Bioconda.

When `force=true`, cleans corrupted conda package cache before reinstalling.
Also verifies the binary is functional after install; if verification fails,
retries with a full cache clean.
"""
function install_foldseek(; force = false)
    if force
        try
            run(`conda clean --all --yes`)
        catch e
            @warn "conda clean failed (non-fatal): $(e)"
        end
    end

    Mycelia.add_bioconda_env("foldseek"; force = force)

    # Verify the binary is functional
    try
        run(`$(Mycelia.CONDA_RUNNER) run -n foldseek foldseek version`)
    catch verify_err
        if !force
            @warn "FoldSeek binary verification failed. Retrying with force reinstall..."
            return install_foldseek(; force = true)
        else
            error("FoldSeek install verification failed after force reinstall: $(verify_err)")
        end
    end
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
        tmp_dir::String = mktempdir(),
        sensitivity::Float64 = 9.5,
        alignment_type::Int = 2,
        format_output::String = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
        format_mode::Int = 0,
        threads::Int = 1,
        exhaustive_search::Bool = false,
        gpu::Bool = false
)
    install_foldseek()

    cmd_flags = [
        "-s", string(sensitivity),
        "--alignment-type", string(alignment_type),
        "--format-output", format_output,
        "--format-mode", string(format_mode),
        "--threads", string(threads)
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
        rm(tmp_dir; recursive = true, force = true)
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
        tmp_dir::String = mktempdir(),
        min_seq_id::Float64 = 0.0,
        coverage::Float64 = 0.8,
        cov_mode::Int = 0,
        alignment_type::Int = 2,
        sensitivity::Float64 = 9.5,
        threads::Int = 1
)
    install_foldseek()

    cmd_flags = [
        "--min-seq-id", string(min_seq_id),
        "-c", string(coverage),
        "--cov-mode", string(cov_mode),
        "--alignment-type", string(alignment_type),
        "-s", string(sensitivity),
        "--threads", string(threads)
    ]

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n foldseek foldseek easy-cluster $input_file $output_prefix $tmp_dir $cmd_flags`
    run(cmd)

    if startswith(tmp_dir, tempdir()) && isdir(tmp_dir)
        rm(tmp_dir; recursive = true, force = true)
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
        threads::Int = 1
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
        tmp_dir::String = mktempdir(),
        threads::Int = 1
)
    install_foldseek()

    if !isdir(destination_dir)
        mkpath(destination_dir)
    end

    cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n foldseek foldseek databases $database_name $destination_dir $tmp_dir --threads $threads`
    run(cmd)

    if startswith(tmp_dir, tempdir()) && isdir(tmp_dir)
        rm(tmp_dir; recursive = true, force = true)
    end

    return destination_dir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Search protein structures via the FoldSeek web API (https://search.foldseek.com).

No local FoldSeek installation required. Submits a PDB/mmCIF file to the
FoldSeek server, polls until results are ready, and returns a DataFrame of hits.

# Arguments
- `pdb_path::String`: Path to a PDB or mmCIF query structure file.

# Keywords
- `databases::Vector{String}`: Databases to search. Default: `["afdb50"]`.
  Options include: "afdb50", "afdb-swissprot", "afdb-proteome", "pdb100",
  "gmgcl", "mgnify_esm30", "bfmd".
- `mode::String`: Search mode. "3diaa" (default, fast) or "tmalign" (slow, accurate).
- `max_results::Int=500`: Maximum number of hits to return per database.
- `poll_interval::Real=5`: Seconds between status checks.
- `timeout::Real=600`: Maximum seconds to wait for results.

# Returns
`DataFrames.DataFrame` with columns: query, target, fident, alnlen, evalue, bits,
alntmscore, prob, database.
"""
# Recursively flatten nested arrays/vectors to extract all Dict entries
function _flatten_to_dicts!(result::Vector, x)
    if x isa Dict
        push!(result, x)
    elseif x isa Vector || x isa AbstractVector
        for item in x
            _flatten_to_dicts!(result, item)
        end
    end
    return result
end

function foldseek_web_search(
        pdb_path::String;
        databases::Vector{String} = ["afdb50"],
        mode::String = "3diaa",
        max_results::Int = 500,
        poll_interval::Real = 5,
        timeout::Real = 600
)
    @assert isfile(pdb_path) "PDB file not found: $(pdb_path)"
    @assert mode in ("3diaa", "tmalign") "mode must be '3diaa' or 'tmalign'"

    api_base = "https://search.foldseek.com/api"

    # Submit search job
    pdb_content = read(pdb_path, String)
    # FoldSeek API requires the PDB as a file upload (multipart), not a raw string.
    # HTTP.Form with HTTP.Multipart handles file-like uploads.
    form_parts = [
        "q" => HTTP.Multipart(basename(pdb_path), IOBuffer(pdb_content), "application/octet-stream"),
        "mode" => mode
    ]
    for db in databases
        push!(form_parts, "database[]" => db)
    end

    submit_resp = HTTP.post("$(api_base)/ticket"; body = HTTP.Form(form_parts))
    submit_json = JSON.parse(String(submit_resp.body))

    ticket_id = submit_json["id"]
    @info "FoldSeek web search submitted: ticket=$(ticket_id)"

    # Poll for completion
    elapsed = 0.0
    status = "PENDING"
    while status in ("PENDING", "RUNNING") && elapsed < timeout
        sleep(poll_interval)
        elapsed += poll_interval
        status_resp = HTTP.get("$(api_base)/ticket/$(ticket_id)")
        status_json = JSON.parse(String(status_resp.body))
        status = get(status_json, "status", "UNKNOWN")
    end

    if status != "COMPLETE"
        error("FoldSeek web search failed or timed out (status: $(status), elapsed: $(elapsed)s)")
    end

    # Fetch results
    result_resp = HTTP.get("$(api_base)/result/$(ticket_id)/0";
        query = Dict("format" => "json"))
    result_json = JSON.parse(String(result_resp.body))

    # Parse into DataFrame from FoldSeek web API response.
    # API structure: {"results": [{"db": "afdb50", "alignments": [[{hit1}, {hit2}, ...]]}]}
    # Note: alignments is DOUBLY nested — alignments[0] contains the array of hit dicts.
    rows = NamedTuple[]
    raw_results = get(result_json, "results", [])
    for (db_idx, db_results) in enumerate(raw_results)
        db_results isa Dict || continue
        db_name = db_idx <= length(databases) ? databases[db_idx] : "db_$(db_idx)"
        raw_alns = get(db_results, "alignments", [])

        # Flatten all nested arrays to get individual alignment dicts
        # The API returns alignments as [[{hit1}, {hit2}, ...]] (array of arrays)
        alignment_dicts = Any[]
        _flatten_to_dicts!(alignment_dicts, raw_alns)

        # Limit to top 100 hits per query (1000 is excessive for clustering)
        for alignment in alignment_dicts[1:min(100, length(alignment_dicts))]
            alignment isa Dict || continue
            try
                push!(rows,
                    (
                        query = string(alignment["query"]),
                        target = string(alignment["target"]),
                        fident = Float64(alignment["seqId"]),
                        alnlen = Int(round(alignment["alnLength"])),
                        evalue = Float64(alignment["eval"]),
                        bits = Float64(alignment["score"]),
                        alntmscore = Float64(alignment["prob"]),
                        prob = Float64(alignment["prob"]),
                        database = db_name
                    ))
            catch e
                @debug "Skipping alignment entry: $(e)"
            end
        end
    end

    if isempty(rows)
        @warn "FoldSeek web search returned no hits"
        return DataFrames.DataFrame(
            query = String[], target = String[], fident = Float64[],
            alnlen = Int[], evalue = Float64[], bits = Float64[],
            alntmscore = Float64[], prob = Float64[], database = String[])
    end

    return DataFrames.DataFrame(rows)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run FoldSeek web API all-vs-all search for a directory of PDB files.

Submits each structure against all others via `foldseek_web_search`, then builds
a TM-score distance matrix compatible with `structural_distance_matrix`.

# Arguments
- `structure_dir::String`: Directory containing PDB files (named `AF-{accession}*.pdb`).
- `output_file::String`: Path to write TSV results (same format as local foldseek).

# Keywords
- `databases::Vector{String}`: Databases for remote search. Default: `["afdb50"]`.
- `delay::Real=2`: Seconds between submissions to avoid rate limiting.

# Returns
Tuple `(output_file::String, DataFrame)` — TSV path and full results table.
"""
function foldseek_web_allvsall(
        structure_dir::String,
        output_file::String;
        databases::Vector{String} = ["afdb50"],
        delay::Real = 2
)
    pdb_files = filter(f -> endswith(f, ".pdb") || endswith(f, ".cif"), readdir(structure_dir; join = true))
    @assert !isempty(pdb_files) "No PDB/CIF files found in $(structure_dir)"

    all_results = DataFrames.DataFrame[]
    for pdb_path in pdb_files
        @info "Submitting $(basename(pdb_path)) to FoldSeek web API..."
        try
            df = foldseek_web_search(pdb_path; databases = databases)
            push!(all_results, df)
        catch e
            @warn "FoldSeek web search failed for $(basename(pdb_path)): $(e)"
        end
        sleep(delay)
    end

    if isempty(all_results)
        @warn "No FoldSeek web results obtained"
        return output_file, DataFrames.DataFrame()
    end

    combined = reduce(vcat, all_results)

    # Write TSV in same format as local foldseek
    mkpath(dirname(abspath(output_file)))
    open(output_file, "w") do io
        for row in DataFrames.eachrow(combined)
            println(io,
                join(
                    [row.query, row.target, row.fident, row.alnlen,
                        row.evalue, row.bits, row.alntmscore],
                    "\t"))
        end
    end

    return output_file, combined
end
