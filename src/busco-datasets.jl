"""
    run_busco_download_datasets(datasets::Union{String,Vector{String}}="all"; outdir::String="busco_downloads", threads::Int=get_default_threads())

Pre-download BUSCO lineage datasets using the BUSCO CLI `--download` flag. This is useful to warm caches for CI or cluster environments.

# Arguments
- `datasets`: One of a dataset name, a vector of dataset names, or groups like "all", "prokaryota", "eukaryota", "virus".
- `outdir`: Destination directory for the downloads.
- `threads`: Number of CPU threads to use for downloads (passed to BUSCO `--cpu`).

# Notes
- BUSCO requires network access for downloads.
- Uses the same Bioconda environment as `run_busco`.
"""
function run_busco_download_datasets(datasets::Union{String, Vector{String}} = "all";
        outdir::String = "busco_downloads", threads::Int = get_default_threads())
    add_bioconda_env("busco")
    mkpath(outdir)

    # Normalize datasets to vector of strings
    ds_vec = datasets isa String ? [datasets] : datasets

    for ds in ds_vec
        cmd_parts = [
            "busco",
            "--download", ds,
            "--out_path", outdir,
            "--cpu", string(threads)
        ]
        try
            println("Downloading BUSCO dataset: $(ds)")
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n busco $cmd_parts`)
        catch e
            error("BUSCO dataset download failed for $(ds): $(e)")
        end
    end

    return outdir
end

"""
    list_busco_datasets(; parse_output::Bool=true)

List available BUSCO lineage datasets via the BUSCO CLI `--list-datasets` flag.

# Keywords
- `parse_output`: When true, return parsed dataset identifiers matching `*_odb*` in addition to the raw CLI output.

# Returns
- Named tuple with `datasets::Vector{String}` (may be empty if parsing is disabled or no matches are found) and `raw::String` containing the BUSCO stdout.

# Notes
- Shares the same Bioconda environment as other BUSCO wrappers.
"""
function list_busco_datasets(; parse_output::Bool = true)
    add_bioconda_env("busco")

    cmd_parts = ["busco", "--list-datasets"]

    output = try
        read(`$(Mycelia.CONDA_RUNNER) run -n busco $cmd_parts`, String)
    catch e
        error("BUSCO dataset listing failed: $(e)")
    end

    datasets = parse_output ? _parse_busco_dataset_list(output) : String[]

    return (datasets = datasets, raw = output)
end

function _parse_busco_dataset_list(output::AbstractString)
    datasets = String[]
    for line in split(output, '\n')
        for m in eachmatch(r"([A-Za-z0-9_]+_odb\d+)", line)
            push!(datasets, m.captures[1])
        end
    end
    return unique(datasets)
end
