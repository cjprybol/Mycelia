const NCBI_DATASETS_ENV = "ncbi-datasets-cli"

function _datasets_flag_parts(flags::Dict{String,Any})
    parts = String[]
    for (key, value) in flags
        if value === nothing || value === false
            continue
        end
        flag = startswith(key, "-") ? key : "--" * key
        if value === true
            push!(parts, flag)
        elseif value isa AbstractVector
            push!(parts, flag)
            push!(parts, join(string.(value), ","))
        else
            push!(parts, flag)
            push!(parts, string(value))
        end
    end
    return parts
end

function _datasets_merge_flags(flags::Dict{String,Any}; api_key::String="", debug::Bool=false, no_progressbar::Bool=true)
    merged = Dict{String,Any}()
    if !isempty(api_key)
        merged["api-key"] = api_key
    end
    debug && (merged["debug"] = true)
    no_progressbar && (merged["no-progressbar"] = true)
    for (key, value) in flags
        merged[key] = value
    end
    return merged
end

function _datasets_cmd(parts::Vector{String}; live_stream::Bool=true)
    conda_parts = ["run"]
    if live_stream
        push!(conda_parts, "--live-stream")
    end
    append!(conda_parts, ["-n", NCBI_DATASETS_ENV])
    full_parts = vcat([Mycelia.CONDA_RUNNER], conda_parts, parts)
    return Cmd(full_parts)
end

function _datasets_cmd_parts(cmd_type::String, subcmd::Union{String,Nothing}, args::Vector{String}, flags::Dict{String,Any})
    parts = ["datasets", cmd_type]
    if subcmd !== nothing && !isempty(subcmd)
        push!(parts, subcmd)
    end
    append!(parts, args)
    append!(parts, _datasets_flag_parts(flags))
    return parts
end

function _dataformat_cmd_parts(schema::String; format::String="tsv", fields::Vector{String}=String[], template::String="")
    parts = ["dataformat", format, schema]
    flags = Dict{String,Any}()
    if !isempty(template)
        flags["template"] = template
    elseif !isempty(fields)
        flags["fields"] = join(fields, ",")
    end
    append!(parts, _datasets_flag_parts(flags))
    return parts
end

function _datasets_validate_choice(value::String, valid::Vector{String}, label::String)
    if !(value in valid)
        error("Invalid $(label): $(value). Must be one of $(join(valid, ", "))")
    end
end

function _datasets_default_zip_name(input::AbstractString)
    safe = replace(input, r"\s+" => "_")
    return endswith(lowercase(safe), ".zip") ? safe : safe * ".zip"
end

"""
    run_datasets_cli(cmd_type::String, subcmd::Union{String,Nothing}, args::Vector{String};
        flags::Dict{String,Any}=Dict{String,Any}(),
        capture_output::Bool=false,
        api_key::String="",
        debug::Bool=false,
        no_progressbar::Bool=true,
        max_attempts::Int=1,
        initial_retry_delay::Float64=5.0
    )

Run a `datasets` CLI command inside the `ncbi-datasets-cli` Conda environment.
Returns the captured output when `capture_output=true`, otherwise returns `nothing`.
"""
function run_datasets_cli(cmd_type::String, subcmd::Union{String,Nothing}, args::Vector{String};
    flags::Dict{String,Any}=Dict{String,Any}(),
    capture_output::Bool=false,
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=true,
    max_attempts::Int=1,
    initial_retry_delay::Float64=5.0
)
    add_bioconda_env(NCBI_DATASETS_ENV)
    merged_flags = _datasets_merge_flags(flags; api_key=api_key, debug=debug, no_progressbar=no_progressbar)
    cmd_parts = _datasets_cmd_parts(cmd_type, subcmd, args, merged_flags)
    cmd = _datasets_cmd(cmd_parts; live_stream=!capture_output)
    if capture_output
        return with_retry(max_attempts=max_attempts, initial_delay=initial_retry_delay) do
            read(cmd, String)
        end
    end
    return with_retry(max_attempts=max_attempts, initial_delay=initial_retry_delay) do
        run(cmd)
    end
end

"""
    run_dataformat_cli(args::Vector{String};
        flags::Dict{String,Any}=Dict{String,Any}(),
        capture_output::Bool=false
    )

Run the `dataformat` CLI inside the `ncbi-datasets-cli` Conda environment.
"""
function run_dataformat_cli(args::Vector{String};
    flags::Dict{String,Any}=Dict{String,Any}(),
    capture_output::Bool=false
)
    add_bioconda_env(NCBI_DATASETS_ENV)
    cmd_parts = ["dataformat"]
    append!(cmd_parts, args)
    append!(cmd_parts, _datasets_flag_parts(flags))
    cmd = _datasets_cmd(cmd_parts; live_stream=!capture_output)
    return capture_output ? read(cmd, String) : run(cmd)
end

"""
    datasets_cli_version() -> String

Return the installed NCBI datasets CLI version string.
"""
function datasets_cli_version()
    add_bioconda_env(NCBI_DATASETS_ENV)
    cmd = _datasets_cmd(["datasets", "--version"]; live_stream=false)
    return strip(read(cmd, String))
end

"""
    parse_datasets_json(output::AbstractString)

Parse a JSON string produced by the datasets CLI into a dictionary or vector.
"""
function parse_datasets_json(output::AbstractString)
    stripped = strip(output)
    isempty(stripped) && error("datasets JSON output is empty")
    return JSON.parse(stripped)
end

"""
    parse_datasets_jsonl(output::AbstractString) -> Vector{Dict{String,Any}}

Parse JSON Lines output from the datasets CLI into a vector of dictionaries.
"""
function parse_datasets_jsonl(output::AbstractString)::Vector{Dict{String,Any}}
    results = Vector{Dict{String,Any}}()
    for line in split(output, '\n')
        stripped = strip(line)
        if !isempty(stripped)
            push!(results, JSON.parse(stripped))
        end
    end
    return results
end

"""
    datasets_dataformat(input_file::String;
        schema::String="genome",
        format::String="tsv",
        fields::Vector{String}=String[],
        template::String=""
    )

Run `dataformat` on a JSONL file and return a DataFrame for TSV output or a raw string otherwise.
"""
function datasets_dataformat(input_file::String;
    schema::String="genome",
    format::String="tsv",
    fields::Vector{String}=String[],
    template::String=""
)
    add_bioconda_env(NCBI_DATASETS_ENV)
    cmd_parts = _dataformat_cmd_parts(schema; format=format, fields=fields, template=template)
    push!(cmd_parts, "--inputfile")
    push!(cmd_parts, input_file)
    cmd = _datasets_cmd(cmd_parts; live_stream=false)
    io = open(cmd)
    try
        if lowercase(format) == "tsv"
            return CSV.read(io, DataFrames.DataFrame, delim='\t', header=1)
        end
        return read(io, String)
    finally
        close(io)
    end
end

function _datasets_summary_dataframe(summary_subcmd::String, args::Vector{String};
    schema::String,
    flags::Dict{String,Any}=Dict{String,Any}(),
    fields::Vector{String}=String[],
    template::String="",
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=false,
    max_attempts::Int=1,
    initial_retry_delay::Float64=5.0
)
    add_bioconda_env(NCBI_DATASETS_ENV)
    merged_flags = _datasets_merge_flags(flags; api_key=api_key, debug=debug, no_progressbar=no_progressbar)
    cmd_parts = _datasets_cmd_parts("summary", summary_subcmd, args, merged_flags)
    return with_retry(max_attempts=max_attempts, initial_delay=initial_retry_delay) do
        datasets_cmd = join(Base.shell_escape.(cmd_parts), " ")
        dataformat_parts = _dataformat_cmd_parts(schema; format="tsv", fields=fields, template=template)
        dataformat_cmd = join(Base.shell_escape.(dataformat_parts), " ")
        full_cmd = "$(datasets_cmd) | $(dataformat_cmd)"
        cmd = Cmd([Mycelia.CONDA_RUNNER, "run", "-n", NCBI_DATASETS_ENV, "bash", "-lc", full_cmd])
        io = open(cmd)
        try
            CSV.read(io, DataFrames.DataFrame, delim='\t', header=1)
        finally
            close(io)
        end
    end
end

function _datasets_summary(summary_subcmd::String, args::Vector{String};
    schema::String,
    as_dataframe::Bool=true,
    as_json_lines::Bool=true,
    flags::Dict{String,Any}=Dict{String,Any}(),
    fields::Vector{String}=String[],
    template::String="",
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=false,
    max_attempts::Int=1,
    initial_retry_delay::Float64=5.0
)
    summary_flags = copy(flags)
    as_json_lines && (summary_flags["as-json-lines"] = true)
    if as_dataframe
        summary_flags["as-json-lines"] = true
        return _datasets_summary_dataframe(summary_subcmd, args;
            schema=schema,
            flags=summary_flags,
            fields=fields,
            template=template,
            api_key=api_key,
            debug=debug,
            no_progressbar=no_progressbar,
            max_attempts=max_attempts,
            initial_retry_delay=initial_retry_delay
        )
    end
    output = run_datasets_cli("summary", summary_subcmd, args;
        flags=summary_flags,
        capture_output=true,
        api_key=api_key,
        debug=debug,
        no_progressbar=no_progressbar,
        max_attempts=max_attempts,
        initial_retry_delay=initial_retry_delay
    )
    return as_json_lines ? parse_datasets_jsonl(output) : parse_datasets_json(output)
end

"""
    datasets_genome_summary(;
        taxon::Union{String,Nothing}=nothing,
        accession::Union{String,Nothing}=nothing,
        assembly_source::String="all",
        limit::Union{String,Int}="all",
        as_json_lines::Bool=true,
        as_dataframe::Bool=true,
        fields::Vector{String}=["accession"],
        api_key::String="",
        debug::Bool=false,
        no_progressbar::Bool=false
    )

Fetch genome summary metadata via the datasets CLI.
Returns a DataFrame by default when `as_dataframe=true`.
"""
function datasets_genome_summary(;
    taxon::Union{String,Nothing}=nothing,
    accession::Union{String,Nothing}=nothing,
    assembly_source::String="all",
    limit::Union{String,Int}="all",
    as_json_lines::Bool=true,
    as_dataframe::Bool=true,
    fields::Vector{String}=["accession"],
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=false
)
    if (taxon === nothing && accession === nothing) || (taxon !== nothing && accession !== nothing)
        error("Provide exactly one of taxon or accession")
    end
    input_type = taxon !== nothing ? "taxon" : "accession"
    input_value = taxon !== nothing ? taxon : accession
    flags = Dict{String,Any}()
    if !isempty(assembly_source)
        flags["assembly-source"] = assembly_source
    end
    if !(limit === nothing)
        flags["limit"] = string(limit)
    end
    return _datasets_summary("genome", [input_type, input_value];
        schema="genome",
        as_dataframe=as_dataframe,
        as_json_lines=as_json_lines,
        flags=flags,
        fields=fields,
        api_key=api_key,
        debug=debug,
        no_progressbar=no_progressbar
    )
end

function datasets_genome_summary(taxon::Union{String,Nothing}; kwargs...)
    return datasets_genome_summary(; taxon=taxon, kwargs...)
end

"""
    datasets_download_genome(input::String;
        input_type::String="taxon",
        include::Vector{String}=["genome", "gff3"],
        dehydrated::Bool=false,
        outdir::String=pwd(),
        filename::String="",
        extract::Bool=false,
        api_key::String="",
        debug::Bool=false,
        no_progressbar::Bool=true
    )

Download a genome package using the datasets CLI. Returns a named tuple with zip_path and directory.
"""
function datasets_download_genome(input::String;
    input_type::String="taxon",
    include::Vector{String}=["genome", "gff3"],
    dehydrated::Bool=false,
    outdir::String=pwd(),
    filename::String="",
    extract::Bool=false,
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=true,
    max_attempts::Int=3,
    initial_retry_delay::Float64=10.0
)
    _datasets_validate_choice(input_type, ["taxon", "accession"], "input_type")
    zip_name = isempty(filename) ? _datasets_default_zip_name(input) : _datasets_default_zip_name(filename)
    outpath = joinpath(outdir, zip_name)
    mkpath(outdir)
    flags = Dict{String,Any}()
    if !isempty(include)
        flags["include"] = join(include, ",")
    end
    dehydrated && (flags["dehydrated"] = true)
    flags["filename"] = outpath
    if isfile(outpath)
        @info "$(outpath) already exists, skipping download..."
    else
        run_datasets_cli("download", "genome", [input_type, input];
            flags=flags,
            capture_output=false,
            api_key=api_key,
            debug=debug,
            no_progressbar=no_progressbar,
            max_attempts=max_attempts,
            initial_retry_delay=initial_retry_delay
        )
    end
    extract_dir = nothing
    if extract
        extract_dir = joinpath(outdir, replace(basename(outpath), r"\.zip$" => ""))
        if !isdir(extract_dir) || isempty(readdir(extract_dir))
            mkpath(extract_dir)
            run(`unzip -q -d $(extract_dir) $(outpath)`)
        end
    end
    return (zip_path = outpath, directory = extract_dir)
end

"""
    datasets_gene_summary(input::String;
        input_type::String="symbol",
        taxon::String="human",
        as_json_lines::Bool=true,
        as_dataframe::Bool=true,
        fields::Vector{String}=["gene-id"],
        api_key::String="",
        debug::Bool=false,
        no_progressbar::Bool=false
    )

Fetch gene summary metadata via the datasets CLI.
"""
function datasets_gene_summary(input::String;
    input_type::String="symbol",
    taxon::String="human",
    as_json_lines::Bool=true,
    as_dataframe::Bool=true,
    fields::Vector{String}=["gene-id"],
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=false
)
    _datasets_validate_choice(input_type, ["gene-id", "symbol", "accession"], "input_type")
    flags = Dict{String,Any}()
    if !isempty(taxon)
        flags["taxon"] = taxon
    end
    return _datasets_summary("gene", [input_type, input];
        schema="gene",
        as_dataframe=as_dataframe,
        as_json_lines=as_json_lines,
        flags=flags,
        fields=fields,
        api_key=api_key,
        debug=debug,
        no_progressbar=no_progressbar
    )
end

"""
    datasets_download_gene(input::String;
        input_type::String="symbol",
        taxon::String="human",
        include::Vector{String}=["gene", "protein", "cds"],
        dehydrated::Bool=false,
        outdir::String=pwd(),
        filename::String="",
        extract::Bool=false,
        api_key::String="",
        debug::Bool=false,
        no_progressbar::Bool=true
    )

Download gene data using the datasets CLI. Returns a named tuple with zip_path and directory.
"""
function datasets_download_gene(input::String;
    input_type::String="symbol",
    taxon::String="human",
    include::Vector{String}=["gene", "protein", "cds"],
    dehydrated::Bool=false,
    outdir::String=pwd(),
    filename::String="",
    extract::Bool=false,
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=true,
    max_attempts::Int=3,
    initial_retry_delay::Float64=10.0
)
    _datasets_validate_choice(input_type, ["gene-id", "symbol", "accession"], "input_type")
    zip_name = isempty(filename) ? _datasets_default_zip_name(input) : _datasets_default_zip_name(filename)
    outpath = joinpath(outdir, zip_name)
    mkpath(outdir)
    flags = Dict{String,Any}()
    if !isempty(include)
        flags["include"] = join(include, ",")
    end
    if !isempty(taxon)
        flags["taxon"] = taxon
    end
    dehydrated && (flags["dehydrated"] = true)
    flags["filename"] = outpath
    if isfile(outpath)
        @info "$(outpath) already exists, skipping download..."
    else
        run_datasets_cli("download", "gene", [input_type, input];
            flags=flags,
            capture_output=false,
            api_key=api_key,
            debug=debug,
            no_progressbar=no_progressbar,
            max_attempts=max_attempts,
            initial_retry_delay=initial_retry_delay
        )
    end
    extract_dir = nothing
    if extract
        extract_dir = joinpath(outdir, replace(basename(outpath), r"\.zip$" => ""))
        if !isdir(extract_dir) || isempty(readdir(extract_dir))
            mkpath(extract_dir)
            run(`unzip -q -d $(extract_dir) $(outpath)`)
        end
    end
    return (zip_path = outpath, directory = extract_dir)
end

"""
    datasets_download_orthologs(gene_id::Int;
        taxon_filter::String,
        include::Vector{String}=["gene"],
        outdir::String=pwd(),
        filename::String="",
        extract::Bool=false,
        api_key::String="",
        debug::Bool=false,
        no_progressbar::Bool=true
    )

Download orthologs for a gene ID using the datasets CLI.
"""
function datasets_download_orthologs(gene_id::Int;
    taxon_filter::String,
    include::Vector{String}=["gene"],
    outdir::String=pwd(),
    filename::String="",
    extract::Bool=false,
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=true,
    max_attempts::Int=3,
    initial_retry_delay::Float64=10.0
)
    isempty(taxon_filter) && error("taxon_filter must be provided")
    zip_name = isempty(filename) ? _datasets_default_zip_name("orthologs_$(gene_id)") : _datasets_default_zip_name(filename)
    outpath = joinpath(outdir, zip_name)
    mkpath(outdir)
    flags = Dict{String,Any}("taxon" => taxon_filter, "ortholog" => true, "filename" => outpath)
    if !isempty(include)
        flags["include"] = join(include, ",")
    end
    if isfile(outpath)
        @info "$(outpath) already exists, skipping download..."
    else
        run_datasets_cli("download", "gene", ["gene-id", string(gene_id)];
            flags=flags,
            capture_output=false,
            api_key=api_key,
            debug=debug,
            no_progressbar=no_progressbar,
            max_attempts=max_attempts,
            initial_retry_delay=initial_retry_delay
        )
    end
    extract_dir = nothing
    if extract
        extract_dir = joinpath(outdir, replace(basename(outpath), r"\.zip$" => ""))
        if !isdir(extract_dir) || isempty(readdir(extract_dir))
            mkpath(extract_dir)
            run(`unzip -q -d $(extract_dir) $(outpath)`)
        end
    end
    return (zip_path = outpath, directory = extract_dir)
end

"""
    datasets_virus_summary(; taxon::String="", lineage::String="", host::String="",
        as_json_lines::Bool=true, as_dataframe::Bool=true, fields::Vector{String}=["accession"],
        api_key::String="", debug::Bool=false, no_progressbar::Bool=false)

Fetch virus summary metadata via the datasets CLI.
"""
function datasets_virus_summary(;
    taxon::String="",
    lineage::String="",
    host::String="",
    as_json_lines::Bool=true,
    as_dataframe::Bool=true,
    fields::Vector{String}=["accession"],
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=false
)
    flags = Dict{String,Any}()
    !isempty(taxon) && (flags["taxon"] = taxon)
    !isempty(lineage) && (flags["lineage"] = lineage)
    !isempty(host) && (flags["host"] = host)
    isempty(flags) && error("Provide at least one of taxon, lineage, or host")
    return _datasets_summary("virus", String[];
        schema="virus-genome",
        as_dataframe=as_dataframe,
        as_json_lines=as_json_lines,
        flags=flags,
        fields=fields,
        api_key=api_key,
        debug=debug,
        no_progressbar=no_progressbar
    )
end

"""
    datasets_download_virus(; taxon::String="", lineage::String="", host::String="",
        include::Vector{String}=["genome", "cds", "protein"],
        dehydrated::Bool=false,
        outdir::String=pwd(),
        filename::String="",
        extract::Bool=false,
        api_key::String="",
        debug::Bool=false,
        no_progressbar::Bool=true
    )

Download virus datasets using the datasets CLI. Returns a named tuple with zip_path and directory.
"""
function datasets_download_virus(;
    taxon::String="",
    lineage::String="",
    host::String="",
    include::Vector{String}=["genome", "cds", "protein"],
    dehydrated::Bool=false,
    outdir::String=pwd(),
    filename::String="",
    extract::Bool=false,
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=true,
    max_attempts::Int=3,
    initial_retry_delay::Float64=10.0
)
    flags = Dict{String,Any}()
    !isempty(taxon) && (flags["taxon"] = taxon)
    !isempty(lineage) && (flags["lineage"] = lineage)
    !isempty(host) && (flags["host"] = host)
    isempty(flags) && error("Provide at least one of taxon, lineage, or host")
    if !isempty(include)
        flags["include"] = join(include, ",")
    end
    dehydrated && (flags["dehydrated"] = true)
    zip_name = isempty(filename) ? _datasets_default_zip_name("virus") : _datasets_default_zip_name(filename)
    outpath = joinpath(outdir, zip_name)
    mkpath(outdir)
    flags["filename"] = outpath
    if isfile(outpath)
        @info "$(outpath) already exists, skipping download..."
    else
        run_datasets_cli("download", "virus", String[];
            flags=flags,
            capture_output=false,
            api_key=api_key,
            debug=debug,
            no_progressbar=no_progressbar,
            max_attempts=max_attempts,
            initial_retry_delay=initial_retry_delay
        )
    end
    extract_dir = nothing
    if extract
        extract_dir = joinpath(outdir, replace(basename(outpath), r"\.zip$" => ""))
        if !isdir(extract_dir) || isempty(readdir(extract_dir))
            mkpath(extract_dir)
            run(`unzip -q -d $(extract_dir) $(outpath)`)
        end
    end
    return (zip_path = outpath, directory = extract_dir)
end

"""
    datasets_taxonomy_summary(taxon::String;
        as_json_lines::Bool=true,
        as_dataframe::Bool=true,
        fields::Vector{String}=String[],
        template::String="tax-summary",
        api_key::String="",
        debug::Bool=false,
        no_progressbar::Bool=false
    )

Fetch taxonomy summaries via the datasets CLI.
"""
function datasets_taxonomy_summary(taxon::String;
    as_json_lines::Bool=true,
    as_dataframe::Bool=true,
    fields::Vector{String}=String[],
    template::String="tax-summary",
    api_key::String="",
    debug::Bool=false,
    no_progressbar::Bool=false
)
    return _datasets_summary("taxonomy", ["taxon", taxon];
        schema="taxonomy",
        as_dataframe=as_dataframe,
        as_json_lines=as_json_lines,
        fields=fields,
        template=template,
        api_key=api_key,
        debug=debug,
        no_progressbar=no_progressbar
    )
end

"""
    datasets_taxonomy_tree(taxon::String; api_key::String="", debug::Bool=false, no_progressbar::Bool=false)

Return taxonomy tree information for a taxon using the datasets CLI.
Falls back to the summary output if tree support is unavailable.
"""
function datasets_taxonomy_tree(taxon::String; api_key::String="", debug::Bool=false, no_progressbar::Bool=false)
    flags = Dict{String,Any}("children" => true, "as-json-lines" => true)
    try
        output = run_datasets_cli("summary", "taxonomy", ["taxon", taxon];
            flags=flags,
            capture_output=true,
            api_key=api_key,
            debug=debug,
            no_progressbar=no_progressbar
        )
        return parse_datasets_jsonl(output)
    catch e
        @warn "datasets taxonomy tree failed; returning taxonomy summary output" exception=(e, catch_backtrace())
        return datasets_taxonomy_summary(taxon; as_dataframe=false, as_json_lines=true, api_key=api_key, debug=debug, no_progressbar=no_progressbar)
    end
end

"""
    datasets_rehydrate(directory::String; gzip::Bool=true, max_workers::Int=4)

Rehydrate a dehydrated datasets package in-place.
"""
function datasets_rehydrate(directory::String; gzip::Bool=true, max_workers::Int=4)
    flags = Dict{String,Any}("directory" => directory)
    gzip && (flags["gzip"] = true)
    max_workers > 0 && (flags["max-workers"] = max_workers)
    run_datasets_cli("rehydrate", nothing, String[];
        flags=flags,
        capture_output=false
    )
    return directory
end
