# Classification Tools for Metagenomic Analysis
#
# This file contains wrappers for metagenomic classification tools:
# - sourmash: MinHash-based sequence search and classification
# - mash: MinHash sketching and containment screening
# - metaphlan: Marker-based taxonomic profiling
# - metabuli: Fast metagenomic classification
#
# All tools use conda/mamba environments for dependency management.

# ============================================================================
# Sourmash - MinHash-based sequence analysis
# ============================================================================

"""
    run_sourmash_sketch(;input_files, outdir, k_sizes=[21, 31, 51], scaled=1000,
                        molecule="dna", singleton=false, name=nothing)

Create sourmash signatures (sketches) from input sequences.

# Arguments
- `input_files::Vector{String}`: Input FASTA/FASTQ files
- `outdir::String`: Output directory for signature files
- `k_sizes::Vector{Int}`: K-mer sizes to use (default: [21, 31, 51])
- `scaled::Int`: Scaled parameter for sketching (default: 1000)
- `molecule::String`: Molecule type - "dna", "protein", or "dayhoff" (default: "dna")
- `singleton::Bool`: Create individual signatures per sequence (default: false)
- `name::Union{Nothing, String}`: Optional name for the signature

# Returns
Named tuple with:
- `outdir`: Output directory path
- `signatures`: Vector of signature file paths

# Example
```julia
result = run_sourmash_sketch(
    input_files=["reads.fq.gz"],
    outdir="sourmash_output",
    k_sizes=[31],
    scaled=1000
)
```
"""
function run_sourmash_sketch(;
        input_files::Vector{String},
        outdir::String,
        k_sizes::Vector{Int} = [21, 31, 51],
        scaled::Int = 1000,
        molecule::String = "dna",
        singleton::Bool = false,
        name::Union{Nothing, String} = nothing)

    # Validate inputs
    for f in input_files
        isfile(f) || error("Input file not found: $(f)")
    end
    molecule in ["dna", "protein", "dayhoff"] || error("Invalid molecule type: $(molecule)")

    Mycelia.add_bioconda_env("sourmash")
    mkpath(outdir)

    signatures = String[]
    for input_file in input_files
        basename_clean = replace(basename(input_file), r"\.(fa|fasta|fq|fastq)(\.gz)?$"i => "")
        sig_file = joinpath(outdir, "$(basename_clean).sig")

        if !isfile(sig_file)
            k_param = join(["k=$(k)" for k in k_sizes], ",")

            cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n sourmash sourmash sketch $(molecule)
                -p $(k_param),scaled=$(scaled)
                -o $(sig_file)
                $(input_file)`

            if singleton
                cmd = `$(cmd) --singleton`
            end
            if !isnothing(name)
                cmd = `$(cmd) --name $(name)`
            end

            run(cmd)
        end
        push!(signatures, sig_file)
    end

    return (; outdir, signatures)
end

"""
    run_sourmash_search(;query_sig, database_sig, outdir, threshold=0.1,
                        k_size=31, best_only=false, num_results=nothing)

Search query signatures against a database of signatures.

# Arguments
- `query_sig::String`: Query signature file
- `database_sig::String`: Database signature file (or directory)
- `outdir::String`: Output directory
- `threshold::Float64`: Minimum containment threshold (default: 0.1)
- `k_size::Int`: K-mer size to use (default: 31)
- `best_only::Bool`: Return only the best match (default: false)
- `num_results::Union{Nothing, Int}`: Maximum number of results

# Returns
Named tuple with:
- `outdir`: Output directory path
- `results_csv`: Path to CSV results file
"""
function run_sourmash_search(;
        query_sig::String,
        database_sig::String,
        outdir::String,
        threshold::Float64 = 0.1,
        k_size::Int = 31,
        best_only::Bool = false,
        num_results::Union{Nothing, Int} = nothing)
    isfile(query_sig) || error("Query signature not found: $(query_sig)")
    (isfile(database_sig) || isdir(database_sig)) ||
        error("Database not found: $(database_sig)")

    Mycelia.add_bioconda_env("sourmash")
    mkpath(outdir)

    query_basename = replace(basename(query_sig), r"\.sig(\.gz)?$"i => "")
    results_csv = joinpath(outdir, "$(query_basename)_search_results.csv")

    if !isfile(results_csv)
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n sourmash sourmash search
            $(query_sig) $(database_sig)
            -k $(k_size)
            --threshold $(threshold)
            -o $(results_csv)`

        if best_only
            cmd = `$(cmd) --best-only`
        end
        if !isnothing(num_results)
            cmd = `$(cmd) -n $(num_results)`
        end

        run(cmd)
    end

    return (; outdir, results_csv)
end

"""
    run_sourmash_gather(;query_sig, database_sig, outdir,
                        k_size=31, threshold_bp=50000)

Find the best reference genomes matching a metagenome signature.

# Arguments
- `query_sig::String`: Query metagenome signature
- `database_sig::String`: Database of reference signatures
- `outdir::String`: Output directory
- `k_size::Int`: K-mer size (default: 31)
- `threshold_bp::Int`: Minimum base pairs for a match (default: 50000)

# Returns
Named tuple with:
- `outdir`: Output directory path
- `results_csv`: Path to CSV results file
- `results_matches`: Path to matching signatures output
"""
function run_sourmash_gather(;
        query_sig::String,
        database_sig::String,
        outdir::String,
        k_size::Int = 31,
        threshold_bp::Int = 50000)
    isfile(query_sig) || error("Query signature not found: $(query_sig)")
    (isfile(database_sig) || isdir(database_sig)) ||
        error("Database not found: $(database_sig)")

    Mycelia.add_bioconda_env("sourmash")
    mkpath(outdir)

    query_basename = replace(basename(query_sig), r"\.sig(\.gz)?$"i => "")
    results_csv = joinpath(outdir, "$(query_basename)_gather.csv")
    results_matches = joinpath(outdir, "$(query_basename)_gather_matches.sig")

    if !isfile(results_csv)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n sourmash sourmash gather
            $(query_sig) $(database_sig)
            -k $(k_size)
            --threshold-bp $(threshold_bp)
            -o $(results_csv)
            --save-matches $(results_matches)`)
    end

    return (; outdir, results_csv, results_matches)
end

"""
    parse_sourmash_gather_output(file)

Parse a Sourmash gather output CSV file into a DataFrame.

# Arguments
- `file::String`: Path to the sourmash gather CSV output file

# Returns
A DataFrame containing the gather results with columns including:
- `intersect_bp`: Base pairs in the intersection
- `f_orig_query`: Fraction of original query matched
- `f_match`: Fraction of the reference matched
- `name`: Reference genome name
- And other sourmash gather output columns
"""
function parse_sourmash_gather_output(file::String)
    isfile(file) || error("Sourmash gather output not found: $(file)")
    return CSV.read(file, DataFrames.DataFrame)
end

"""
    parse_sourmash_search_output(file)

Parse a Sourmash search output CSV file into a DataFrame.

# Arguments
- `file::String`: Path to the sourmash search CSV output file

# Returns
A DataFrame containing the search results with columns including:
- `similarity`: Jaccard similarity score
- `name`: Reference genome name
- `filename`: Path to signature file
- `md5`: MD5 hash of the signature
"""
function parse_sourmash_search_output(file::String)
    isfile(file) || error("Sourmash search output not found: $(file)")
    return CSV.read(file, DataFrames.DataFrame)
end

# ============================================================================
# Mash - MinHash sketching, screening, and distance
# ============================================================================

"""
    run_mash_sketch(; input_files, outdir, k=21, s=1000, r=false,
                    min_copies=nothing, threads=get_default_threads(),
                    output_prefix=nothing, additional_args=String[])

Create Mash sketches from input sequences or reads.

# Arguments
- `input_files::Vector{String}`: Input FASTA/FASTQ files
- `outdir::String`: Output directory for sketch files
- `k::Int`: K-mer size (default: 21)
- `s::Int`: Sketch size (default: 1000)
- `r::Bool`: Treat input as reads (`-r`, default: false)
- `min_copies::Union{Nothing,Int}`: Minimum copies for reads (`-m`, requires `r=true`)
- `threads::Int`: Thread count (`-p`)
- `output_prefix::Union{Nothing,String}`: Optional prefix for sketch outputs
- `additional_args::Vector{String}`: Extra CLI args appended to `mash sketch`

# Returns
Named tuple with:
- `outdir`: Output directory path
- `sketches`: Vector of sketch file paths
"""
function run_mash_sketch(;
        input_files::Vector{String},
        outdir::String,
        k::Int = 21,
        s::Int = 1000,
        r::Bool = false,
        min_copies::Union{Nothing, Int} = nothing,
        threads::Int = get_default_threads(),
        output_prefix::Union{Nothing, String} = nothing,
        additional_args::Vector{String} = String[])
    for f in input_files
        isfile(f) || error("Input file not found: $(f)")
    end
    k > 0 || error("k must be positive")
    s > 0 || error("s must be positive")
    threads > 0 || error("threads must be positive")
    if !isnothing(min_copies)
        min_copies > 0 || error("min_copies must be positive")
        r || error("min_copies requires r=true for read sketching")
    end

    Mycelia.add_bioconda_env("mash")
    mkpath(outdir)

    sketches = String[]
    for input_file in input_files
        basename_clean = replace(basename(input_file), r"\.(fa|fasta|fq|fastq)(\.gz)?$"i => "")
        prefix = if isnothing(output_prefix)
            joinpath(outdir, basename_clean)
        elseif length(input_files) == 1
            joinpath(outdir, output_prefix)
        else
            joinpath(outdir, "$(output_prefix)_$(basename_clean)")
        end
        sketch_file = prefix * ".msh"

        if !isfile(sketch_file)
            sketch_args = ["sketch", "-k", string(k), "-s", string(s),
                "-p", string(threads), "-o", prefix]
            if r
                push!(sketch_args, "-r")
            end
            if !isnothing(min_copies)
                append!(sketch_args, ["-m", string(min_copies)])
            end
            append!(sketch_args, additional_args)
            push!(sketch_args, input_file)

            sketch_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash $sketch_args`
            run(sketch_cmd)
        end

        push!(sketches, sketch_file)
    end

    return (; outdir, sketches)
end

"""
    run_mash_paste(; out_file, in_files)

Combine Mash sketch files into a single database sketch.

# Arguments
- `out_file::String`: Output sketch path (or prefix) for the combined database
- `in_files::Vector{String}`: Input `.msh` files to combine

# Returns
`String` path to the combined `.msh` file.
"""
function run_mash_paste(;
        out_file::String,
        in_files::Vector{String})
    for f in in_files
        isfile(f) || error("Input sketch not found: $(f)")
    end

    Mycelia.add_bioconda_env("mash")

    output_prefix = endswith(out_file, ".msh") ? out_file[1:(end - 4)] : out_file
    mkpath(dirname(output_prefix))
    output_sketch = output_prefix * ".msh"

    if !isfile(output_sketch)
        paste_args = ["paste", output_prefix]
        append!(paste_args, in_files)
        paste_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash $paste_args`
        run(paste_cmd)
    end

    return output_sketch
end

"""
    run_mash_dist(; reference, query, outdir, threads=get_default_threads(),
                  output_tsv=nothing, additional_args=String[])

Compute Mash distances between a reference and query.

# Arguments
- `reference::String`: Reference sketch or sequence
- `query::String`: Query sketch or sequence
- `outdir::String`: Output directory for results
- `threads::Int`: Thread count (`-p`)
- `output_tsv::Union{Nothing,String}`: Optional output file path
- `additional_args::Vector{String}`: Extra CLI args appended to `mash dist`

# Returns
Named tuple with:
- `outdir`: Output directory path
- `results_tsv`: Path to the distance output file
"""
function run_mash_dist(;
        reference::String,
        query::String,
        outdir::String,
        threads::Int = get_default_threads(),
        output_tsv::Union{Nothing, String} = nothing,
        additional_args::Vector{String} = String[])
    isfile(reference) || error("Reference not found: $(reference)")
    isfile(query) || error("Query not found: $(query)")
    threads > 0 || error("threads must be positive")

    Mycelia.add_bioconda_env("mash")
    mkpath(outdir)

    ref_base = replace(basename(reference), r"\.(msh|fa|fasta|fq|fastq)(\.gz)?$"i => "")
    query_base = replace(basename(query), r"\.(msh|fa|fasta|fq|fastq)(\.gz)?$"i => "")
    results_tsv = isnothing(output_tsv) ?
                  joinpath(outdir, "$(ref_base)_vs_$(query_base)_mash_dist.tsv") :
                  output_tsv

    if !isfile(results_tsv)
        dist_args = ["dist", "-p", string(threads)]
        append!(dist_args, additional_args)
        append!(dist_args, [reference, query])
        dist_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash $dist_args`
        open(results_tsv, "w") do io
            run(pipeline(dist_cmd, stdout = io))
        end
    end

    return (; outdir, results_tsv)
end

"""
    run_mash_screen(; reference, query, outdir, winner_takes_all=true,
                    threads=get_default_threads(), min_identity=nothing,
                    output_tsv=nothing, additional_args=String[])

Screen query sequences against a reference sketch database for containment.

# Arguments
- `reference::String`: Reference sketch file (`.msh`)
- `query::Union{String,Vector{String}}`: Query sequences (FASTA/FASTQ)
- `outdir::String`: Output directory for results
- `winner_takes_all::Bool`: Use Mash `-w` (default: true)
- `threads::Int`: Thread count (`-p`)
- `min_identity::Union{Nothing,Float64}`: Minimum identity threshold (`-i`)
- `output_tsv::Union{Nothing,String}`: Optional output file path
- `additional_args::Vector{String}`: Extra CLI args appended to `mash screen`

# Returns
Named tuple with:
- `outdir`: Output directory path
- `results_tsv`: Path to the screen output file
"""
function run_mash_screen(;
        reference::String,
        query::Union{String, Vector{String}},
        outdir::String,
        winner_takes_all::Bool = true,
        threads::Int = get_default_threads(),
        min_identity::Union{Nothing, Float64} = nothing,
        output_tsv::Union{Nothing, String} = nothing,
        additional_args::Vector{String} = String[])
    isfile(reference) || error("Reference not found: $(reference)")
    query_files = query isa String ? [query] : query
    for f in query_files
        isfile(f) || error("Query not found: $(f)")
    end
    threads > 0 || error("threads must be positive")
    if !isnothing(min_identity) && (min_identity < 0.0 || min_identity > 1.0)
        error("min_identity must be between 0.0 and 1.0")
    end

    Mycelia.add_bioconda_env("mash")
    mkpath(outdir)

    query_base = replace(basename(query_files[1]), r"\.(msh|fa|fasta|fq|fastq)(\.gz)?$"i => "")
    output_name = "$(basename(reference))_vs_$(query_base)_mash_screen.tsv"
    results_tsv = isnothing(output_tsv) ? joinpath(outdir, output_name) : output_tsv

    if !isfile(results_tsv)
        screen_args = ["screen", "-p", string(threads)]
        if winner_takes_all
            push!(screen_args, "-w")
        end
        if !isnothing(min_identity)
            append!(screen_args, ["-i", string(min_identity)])
        end
        append!(screen_args, additional_args)
        push!(screen_args, reference)
        append!(screen_args, query_files)

        screen_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash $screen_args`
        open(results_tsv, "w") do io
            run(pipeline(screen_cmd, stdout = io))
        end
    end

    return (; outdir, results_tsv)
end

"""
    parse_mash_dist_output(file)

Parse a Mash distance output file into a DataFrame.
"""
function parse_mash_dist_output(file::String)
    isfile(file) || error("Mash dist output not found: $(file)")
    return CSV.read(
        file,
        DataFrames.DataFrame;
        delim = '\t',
        header = ["reference", "query", "distance", "p_value", "matching_hashes"],
        comment = "#"
    )
end

"""
    parse_mash_screen_output(file)

Parse a Mash screen output file into a DataFrame.
"""
function parse_mash_screen_output(file::String)
    isfile(file) || error("Mash screen output not found: $(file)")
    return CSV.read(
        file,
        DataFrames.DataFrame;
        delim = '\t',
        header = ["identity", "shared_hashes", "median_multiplicity",
            "p_value", "query", "reference"],
        comment = "#"
    )
end

# ============================================================================
# MetaPhlAn - Marker-based taxonomic profiling
# ============================================================================

function resolve_metaphlan_db_index(; db_index::Union{Nothing, String} = nothing)
    if db_index !== nothing && !isempty(db_index)
        return db_index, true
    end
    value = get(ENV, "METAPHLAN_DB_INDEX", "")
    if !isempty(value)
        return value, true
    end
    return nothing, false
end

function _metaphlan_db_present(db_dir::AbstractString, db_index::Union{Nothing, String} = nothing)
    if !isdir(db_dir)
        return false
    end
    if db_index !== nothing && !isempty(db_index)
        return isfile(joinpath(db_dir, "$(db_index).pkl"))
    end
    return any(name -> endswith(name, ".pkl"), readdir(db_dir))
end

function download_metaphlan_db(;
        db_dir::AbstractString = Mycelia.DEFAULT_METAPHLAN_DB_PATH,
        db_index::Union{Nothing, String} = nothing,
        force::Bool = false)
    db_index_val, _ = resolve_metaphlan_db_index(db_index = db_index)
    if !force && _metaphlan_db_present(db_dir, db_index_val)
        return db_dir
    end
    mkpath(db_dir)
    Mycelia.add_bioconda_env("metaphlan")
    cmd_args = String["metaphlan", "--install", "--db_dir", db_dir]
    if db_index_val !== nothing && !isempty(db_index_val)
        append!(cmd_args, ["--index", db_index_val])
    end
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n metaphlan $(cmd_args)`)
    return db_dir
end

function get_metaphlan_db_path(;
        require::Bool = true,
        db_dir::Union{Nothing, String} = nothing,
        db_index::Union{Nothing, String} = nothing,
        download::Bool = true)
    resolved_dir = db_dir
    if resolved_dir === nothing || isempty(resolved_dir)
        value = get(ENV, "METAPHLAN_DB_DIR", "")
        if isempty(value)
            value = get(ENV, "METAPHLAN_DB_PATH", "")
        end
        if isempty(value)
            value = get(ENV, "METAPHLAN_BOWTIE2DB", "")
            if !isempty(value)
                @warn "METAPHLAN_BOWTIE2DB is deprecated; use METAPHLAN_DB_DIR instead."
            end
        end
        resolved_dir = isempty(value) ? Mycelia.DEFAULT_METAPHLAN_DB_PATH : value
    end

    db_index_val, _ = resolve_metaphlan_db_index(db_index = db_index)
    if _metaphlan_db_present(resolved_dir, db_index_val)
        return resolved_dir
    end

    if download
        return download_metaphlan_db(db_dir = resolved_dir, db_index = db_index_val)
    end

    if require
        error("MetaPhlAn database not found at $(resolved_dir). Set METAPHLAN_DB_DIR (preferred) or METAPHLAN_DB_PATH to override.")
    end
    return nothing
end

"""
    run_metaphlan(;input_file, outdir, input_type="fastq",
                  nprocs=get_default_threads(), db_dir=nothing,
                  db_index=nothing, unknown_estimation=true, stat_q=0.2,
                  long_reads=false)

Run MetaPhlAn for marker-based taxonomic profiling.

# Arguments
- `input_file::String`: Input file (FASTQ, FASTA, or BAM)
- `outdir::String`: Output directory
- `input_type::String`: Input type - "fastq", "fasta", "mapout", "sam" (default: "fastq")
- `nprocs::Int`: Number of threads (default: get_default_threads())
- `db_dir::Union{Nothing, String}`: Path to database directory (passed to MetaPhlAn `--db_dir`)
- `db_index::Union{Nothing, String}`: MetaPhlAn database index (optional; uses METAPHLAN_DB_INDEX when set; otherwise downloads MetaPhlAn's default database)
- `unknown_estimation::Bool`: Include unknown/unclassified fraction estimation (default: true)
- `stat_q::Float64`: Quantile for robust average (default: 0.2)
- `long_reads::Bool`: Use long-read mode (MetaPhlAn `--long_reads`, uses minimap2)

# Returns
Named tuple with:
- `outdir`: Output directory path
- `profile_txt`: Path to taxonomic profile output
- `mapout`: Path to mapping output (also returned as `bowtie2_out` for compatibility)

If `db_dir` is not provided, the MetaPhlAn database is resolved from
METAPHLAN_DB_DIR (preferred) or `$(Mycelia.DEFAULT_METAPHLAN_DB_PATH)` and
downloaded automatically when missing.

# Example
```julia
result = run_metaphlan(
    input_file="reads.fq.gz",
    outdir="metaphlan_output",
    nprocs=8
)
```
"""
function run_metaphlan(;
        input_file::String,
        outdir::String,
        input_type::String = "fastq",
        nprocs::Int = get_default_threads(),
        db_dir::Union{Nothing, String} = nothing,
        db_index::Union{Nothing, String} = nothing,
        unknown_estimation::Bool = true,
        stat_q::Float64 = 0.2,
        long_reads::Bool = false)
    isfile(input_file) || error("Input file not found: $(input_file)")
    input_type in ["fastq", "fasta", "mapout", "sam"] ||
        error("Invalid input_type: $(input_type)")

    resolved_index, index_explicit = resolve_metaphlan_db_index(db_index = db_index)
    if db_dir === nothing || isempty(db_dir)
        db_dir = get_metaphlan_db_path(db_index = resolved_index)
    end

    Mycelia.add_bioconda_env("metaphlan")
    mkpath(outdir)

    input_basename = replace(basename(input_file), r"\.(fa|fasta|fq|fastq|bam|sam)(\.gz)?$"i => "")
    profile_txt = joinpath(outdir, "$(input_basename)_profile.txt")
    mapout_file = joinpath(outdir, "$(input_basename)_mapout.bz2")

    if !isfile(profile_txt)
        cmd_args = String[
        "$(Mycelia.CONDA_RUNNER)", "run", "--live-stream", "-n", "metaphlan",
        "metaphlan",
        input_file,
        "--input_type", input_type,
        "--nproc", string(nprocs),
        "-o", profile_txt,
        "--mapout", mapout_file,
        "--stat_q", string(stat_q)
]

        if !isnothing(db_dir)
            push!(cmd_args, "--db_dir")
            push!(cmd_args, db_dir)
        end
        if index_explicit
            push!(cmd_args, "--index")
            push!(cmd_args, resolved_index)
        end

        if !unknown_estimation
            push!(cmd_args, "--skip_unclassified_estimation")
        end

        if long_reads
            push!(cmd_args, "--long_reads")
        end

        run(Cmd(cmd_args))
    end

    return (; outdir, profile_txt, mapout = mapout_file)
end

"""
    parse_metaphlan_profile(profile_file::String)::DataFrames.DataFrame

Parse a MetaPhlAn profile output into a DataFrame.

# Arguments
- `profile_file::String`: Path to MetaPhlAn profile.txt output

# Returns
DataFrame with columns:
- `clade_name`: Full taxonomic path
- `taxid`: NCBI taxonomy ID (if available)
- `relative_abundance`: Relative abundance percentage
- `additional_species`: Additional species info (if available)
"""
function parse_metaphlan_profile(profile_file::String)::DataFrames.DataFrame
    isfile(profile_file) || error("Profile file not found: $(profile_file)")

    # Read the file, skipping comment lines
    lines = filter(l -> !startswith(l, "#"), readlines(profile_file))

    if isempty(lines)
        return DataFrames.DataFrame(
            clade_name = String[],
            taxid = Union{Int, Missing}[],
            relative_abundance = Float64[],
            additional_species = Union{String, Missing}[]
        )
    end

    # Parse lines
    clade_names = String[]
    taxids = Union{Int, Missing}[]
    abundances = Float64[]
    additional_species = Union{String, Missing}[]

    for line in lines
        isempty(strip(line)) && continue
        parts = split(line, '\t')
        length(parts) >= 3 || continue

        rel_abundance_str = strip(parts[3])
        rel_abundance = tryparse(Float64, rel_abundance_str)
        rel_abundance === nothing && continue

        push!(clade_names, parts[1])
        push!(abundances, rel_abundance)

        if length(parts) >= 2 && !isempty(strip(parts[2]))
            taxid_val = tryparse(Int, strip(parts[2]))
            push!(taxids, isnothing(taxid_val) ? missing : taxid_val)
        else
            push!(taxids, missing)
        end

        if length(parts) >= 4 && !isempty(strip(parts[4]))
            push!(additional_species, strip(parts[4]))
        else
            push!(additional_species, missing)
        end
    end

    return DataFrames.DataFrame(
        clade_name = clade_names,
        taxid = taxids,
        relative_abundance = abundances,
        additional_species = additional_species
    )
end

# ============================================================================
# Metabuli - Fast metagenomic classification
# ============================================================================

function _metabuli_db_present(db_dir::AbstractString)
    if !isdir(db_dir)
        return false
    end
    for marker in ("db.info", "info")
        if isfile(joinpath(db_dir, marker))
            return true
        end
    end
    return false
end

function _resolve_metabuli_db_dir(db_root::AbstractString, db_name::AbstractString)
    target = joinpath(db_root, db_name)
    if _metabuli_db_present(target)
        return target
    end
    candidates = filter(path -> _metabuli_db_present(path), readdir(db_root; join = true))
    if !isempty(candidates)
        named = filter(path -> occursin(db_name, basename(path)), candidates)
        if length(named) == 1
            return first(named)
        elseif length(candidates) == 1
            return first(candidates)
        end
    end
    error("Metabuli database download completed but no valid db.info/info was found under $(db_root).")
end

function download_metabuli_db(;
        db_name::AbstractString,
        db_root::AbstractString = Mycelia.DEFAULT_METABULI_DB_PATH,
        force::Bool = false,
        keep_archive::Bool = false)
    url = get(Mycelia.METABULI_DB_URLS, db_name, nothing)
    if url === nothing
        valid = join(collect(keys(Mycelia.METABULI_DB_URLS)), ", ")
        error("Unknown Metabuli database '$(db_name)'. Available options: $(valid).")
    end
    mkpath(db_root)

    if !force && _metabuli_db_present(joinpath(db_root, db_name))
        return joinpath(db_root, db_name)
    end

    archive_name = basename(url)
    archive_path = joinpath(db_root, archive_name)
    if force || !isfile(archive_path) || filesize(archive_path) == 0
        Mycelia.with_retry() do
            Downloads.download(url, archive_path)
        end
    end
    Mycelia.tar_extract(tarchive = archive_path, directory = db_root)
    db_dir = _resolve_metabuli_db_dir(db_root, db_name)
    if !keep_archive && isfile(archive_path)
        rm(archive_path; force = true)
    end
    return db_dir
end

function get_metabuli_db_path(;
        require::Bool = true,
        db_name::Union{Nothing, String} = nothing,
        db_root::Union{Nothing, String} = nothing,
        download::Bool = true)
    for key in ("METABULI_DB", "METABULI_DB_PATH")
        value = get(ENV, key, "")
        if !isempty(value)
            if isdir(value)
                return value
            end
            if require
                error("Metabuli database directory not found at $(value). Set METABULI_DB/METABULI_DB_PATH to a valid path.")
            end
            return nothing
        end
    end

    root_value = get(ENV, "METABULI_DB_ROOT", "")
    resolved_root = db_root === nothing ?
                    (isempty(root_value) ? Mycelia.DEFAULT_METABULI_DB_PATH : root_value) :
                    db_root
    name_value = get(ENV, "METABULI_DB_NAME", "")
    resolved_name = db_name === nothing ?
                    (isempty(name_value) ? Mycelia.DEFAULT_METABULI_DB_NAME : name_value) :
                    db_name
    candidate = joinpath(resolved_root, resolved_name)
    if _metabuli_db_present(candidate)
        return candidate
    end

    if download
        return download_metabuli_db(db_name = resolved_name, db_root = resolved_root)
    end

    if require
        error("Metabuli database directory not found. Set METABULI_DB/METABULI_DB_PATH or METABULI_DB_NAME (default: $(Mycelia.DEFAULT_METABULI_DB_NAME)).")
    end
    return nothing
end

"""
    run_metabuli_classify(;input_files, outdir, database_path=nothing,
                          seq_mode="1", threads=get_default_threads(),
                          min_score=nothing, min_sp_score=nothing,
                          max_ram_gb=nothing)

Run Metabuli for fast metagenomic classification.

# Arguments
- `input_files::Vector{String}`: Input sequence files (FASTA/FASTQ)
- `outdir::String`: Output directory
- `database_path::Union{Nothing, String}`: Path to Metabuli database. If not provided,
  uses METABULI_DB/METABULI_DB_PATH or auto-downloads the named database under
  METABULI_DB_ROOT (default: `$(Mycelia.DEFAULT_METABULI_DB_PATH)`). Set
  METABULI_DB_NAME to choose a download (default: `$(Mycelia.DEFAULT_METABULI_DB_NAME)`).
- `seq_mode::String`: Sequence mode - "1" for single-end, "2" for paired-end (default: "1")
- `threads::Int`: Number of threads (default: get_default_threads())
- `min_score::Union{Nothing, Float64}`: Minimum score threshold
- `min_sp_score::Union{Nothing, Float64}`: Minimum species-level score
- `max_ram_gb::Union{Nothing, Int}`: Max RAM (GiB) for Metabuli. Defaults to an auto-detected
  value capped at 128 GiB to avoid overcommitting memory on smaller machines.

# Returns
Named tuple with:
- `outdir`: Output directory path
- `report_file`: Path to classification report
- `classifications_file`: Path to per-read classifications

# Example
```julia
result = run_metabuli_classify(
    input_files=["reads.fq.gz"],
    database_path="/path/to/metabuli_db",
    outdir="metabuli_output"
)
```
"""
function _default_metabuli_max_ram_gb()::Int
    total_memory = Int(Sys.total_memory())
    available_memory = try
        Int(Sys.free_memory())
    catch
        total_memory
    end

    slurm_threads = Mycelia._detect_slurm_cpu_allocation()
    slurm_gpus = Mycelia._detect_slurm_gpu_allocation()
    slurm_memory,
    _ = Mycelia._detect_slurm_memory_bytes(
        cpu_allocation = slurm_threads,
        gpu_allocation = slurm_gpus,
        total_memory = total_memory,
        available_memory = available_memory
    )

    usable_memory = slurm_memory === nothing ? min(available_memory, total_memory) :
                    min(slurm_memory, available_memory, total_memory)
    safe_bytes = Int(floor(usable_memory * 0.90))
    gib_bytes = 1024^3
    safe_gb = max(1, Int(floor(safe_bytes / gib_bytes)))
    return min(safe_gb, 128)
end

function run_metabuli_classify(;
        input_files::Vector{String},
        outdir::String,
        database_path::Union{Nothing, String} = nothing,
        seq_mode::String = "1",
        threads::Int = get_default_threads(),
        min_score::Union{Nothing, Float64} = nothing,
        min_sp_score::Union{Nothing, Float64} = nothing,
        max_ram_gb::Union{Nothing, Int} = nothing)

    # Validate inputs
    for f in input_files
        isfile(f) || error("Input file not found: $(f)")
    end

    # Determine database path
    if database_path === nothing || isempty(database_path)
        database_path = get_metabuli_db_path()
    end
    isdir(database_path) || error("Database directory not found: $(database_path)")
    seq_mode in ["1", "2", "3"] ||
        error("Invalid seq_mode: $(seq_mode). Use '1' for single-end, '2' for paired-end interleaved, '3' for paired-end separate")

    max_ram_value = max_ram_gb === nothing ? Mycelia._default_metabuli_max_ram_gb() :
                    max_ram_gb
    max_ram_value > 0 || error("max_ram_gb must be positive (got $(max_ram_value)).")

    Mycelia.add_bioconda_env("metabuli")
    mkpath(outdir)

    job_id = "metabuli_classification"
    report_file = joinpath(outdir, "$(job_id)_report.tsv")
    classifications_file = joinpath(outdir, "$(job_id)_classifications.tsv")

    if !isfile(report_file)
        # Build command
        cmd_args = String[
        "$(Mycelia.CONDA_RUNNER)", "run", "--live-stream", "-n", "metabuli",
        "metabuli", "classify",
        input_files...,
        database_path,
        outdir,
        job_id,
        "--seq-mode", seq_mode,
        "--threads", string(threads),
        "--max-ram", string(max_ram_value)
]

        if !isnothing(min_score)
            push!(cmd_args, "--min-score")
            push!(cmd_args, string(min_score))
        end

        if !isnothing(min_sp_score)
            push!(cmd_args, "--min-sp-score")
            push!(cmd_args, string(min_sp_score))
        end

        run(Cmd(cmd_args))
    end

    return (; outdir, report_file, classifications_file)
end

"""
    run_metabuli_build_db(;reference_fasta, taxonomy_dir, outdir,
                          threads=get_default_threads(), split_num=4096)

Build a Metabuli database from reference sequences.

# Arguments
- `reference_fasta::String`: Input reference sequences (FASTA)
- `taxonomy_dir::String`: Directory with NCBI taxonomy files (names.dmp, nodes.dmp)
- `outdir::String`: Output database directory
- `threads::Int`: Number of threads (default: get_default_threads())
- `split_num::Int`: Number of database splits (default: 4096)

# Returns
Named tuple with:
- `database_path`: Path to the built database
"""
function run_metabuli_build_db(;
        reference_fasta::String,
        taxonomy_dir::String,
        outdir::String,
        threads::Int = get_default_threads(),
        split_num::Int = 4096)
    isfile(reference_fasta) || error("Reference FASTA not found: $(reference_fasta)")
    isdir(taxonomy_dir) || error("Taxonomy directory not found: $(taxonomy_dir)")

    Mycelia.add_bioconda_env("metabuli")
    mkpath(outdir)

    # Check for required taxonomy files
    for required in ["names.dmp", "nodes.dmp"]
        f = joinpath(taxonomy_dir, required)
        isfile(f) || error("Required taxonomy file not found: $(f)")
    end

    if !_metabuli_db_present(outdir)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n metabuli metabuli build
            $(reference_fasta)
            $(taxonomy_dir)
            $(outdir)
            --threads $(threads)
            --split-num $(split_num)`)
    end

    if !_metabuli_db_present(outdir)
        error("Metabuli database build completed but no valid db.info/info was found under $(outdir).")
    end

    return (; database_path = outdir)
end

"""
    parse_metabuli_report(report_file::String)::DataFrames.DataFrame

Parse a Metabuli classification report into a DataFrame.

# Arguments
- `report_file::String`: Path to Metabuli report.tsv output

# Returns
DataFrame with columns: percentage, num_reads, num_direct_reads, rank, taxid, name

# Notes
Metabuli report files have no header row. Column order follows Kraken2 report format:
1. percentage - Percentage of reads classified to the clade rooted at this taxon
2. num_reads - Number of reads classified to the clade (clade_count)
3. num_direct_reads - Number of reads classified directly to this taxon (taxon_count)
4. rank - Taxonomic rank (no rank, superkingdom, phylum, class, order, family, genus, species, etc.)
5. taxid - NCBI taxonomy ID
6. name - Scientific name (may have leading spaces indicating hierarchy depth)
"""
function parse_metabuli_report(report_file::String)::DataFrames.DataFrame
    isfile(report_file) || error("Report file not found: $(report_file)")

    # Metabuli report is tab-separated with no header row
    # Column count may vary, but first 6 are standard Kraken2-style columns
    df = DataFrames.DataFrame(CSV.File(
        report_file;
        delim = '\t',
        header = false
    ))

    # Rename the standard 6 columns (if present)
    standard_names = [:percentage, :num_reads, :num_direct_reads, :rank, :taxid, :name]
    ncols = DataFrames.ncol(df)
    if ncols >= 6
        # Rename first 6 columns to standard names
        rename_pairs = [Symbol("Column", i) => standard_names[i] for i in 1:6]
        DataFrames.rename!(df, rename_pairs...)
    elseif ncols > 0
        # Partial rename for files with fewer columns
        rename_pairs = [Symbol("Column", i) => standard_names[i] for i in 1:ncols]
        DataFrames.rename!(df, rename_pairs...)
    end

    return df
end

"""
    parse_metabuli_classifications(classifications_file::String)::DataFrames.DataFrame

Parse Metabuli per-read classifications into a DataFrame.

# Arguments
- `classifications_file::String`: Path to Metabuli classifications.tsv output

# Returns
DataFrame with per-read classification assignments
"""
function parse_metabuli_classifications(classifications_file::String)::DataFrames.DataFrame
    isfile(classifications_file) ||
        error("Classifications file not found: $(classifications_file)")

    # Read the classifications file
    df = DataFrames.DataFrame(CSV.File(classifications_file; delim = '\t', header = false))

    # Rename columns based on typical metabuli output format
    if DataFrames.ncol(df) >= 3
        DataFrames.rename!(df, 1 => :classified, 2 => :read_id, 3 => :taxid)
    end

    return df
end

"""
    run_clamlst(genome_file::String; 
                db_dir::String=joinpath(homedir(), "workspace", "pymlst", "claMLSTDB"),
                species::String="Escherichia coli", 
                outdir::String=dirname(genome_file),
                threads::Int=get_default_threads(),
                force_db_update::Bool=false)

Run `claMLST` (from `pymlst`) to perform MLST typing on a genome.

# Arguments
- `genome_file`: Path to the input genome FASTA file (gzip supported).
- `db_path`: Path to store/look for the claMLST database. Defaults to `~/workspace/pymlst/claMLSTDB`.
- `species`: Species name to import/search if the database needs initialization (default: "Escherichia coli").
- `outdir`: Directory for output files.
- `threads`: Number of threads to use.
- `force_db_update`: If `true`, re-imports the database even if it exists.

# Details
This function automatically:
1. Installs the `pymlst` conda environment.
2. Checks if the `claMLST` database exists at `db_path`. If not, it imports it for the specified `species`.
3. Decompresses the input genome to a temporary file (required by `claMLST`).
4. Runs the search and returns the path to the output.

# Returns
- `String`: Path to the output directory.
"""
function run_clamlst(genome_file::String;
        db_path::String = joinpath(homedir(), "workspace", "pymlst", "claMLSTDB"),
        species::String = "Escherichia coli",
        outdir::String = dirname(genome_file),
        threads::Int = get_default_threads(),
        force_db_update::Bool = false)
    mkpath(outdir)
    mkpath(dirname(db_path))

    # Fix 1: naming output based on original file, not temp file
    base_name = replace(basename(genome_file), r"\.gz$" => "")
    base_name_no_ext = replace(base_name, r"\.(fasta|fa|fna)$" => "")
    output_tsv = joinpath(outdir, "$(base_name_no_ext)_clamlst.tsv")

    # Check for existence of final output
    if isfile(output_tsv)
        return output_tsv
    end

    # only install if we need to make the file - skip if it's already present
    Mycelia.add_bioconda_env("pymlst")

    # Fix 2: Simplified DB check. 
    # If the database exists, we assume the DB is already present
    db_exists = isfile(db_path) || isdir(db_path)
    if !db_exists || force_db_update
        @info "Initializing claMLST database for $species at $db_dir..."
        try
            # Prepare arguments
            cmd_args = ["run", "--live-stream", "-n", "pymlst",
                "claMLST", "import", db_dir, species]

            # If forcing an update on an existing DB, attempt to use --force if supported, 
            # or the user might need to manually clear the folder. 
            # The error message "use --force to override it" suggests the flag exists.
            if db_exists && force_db_update
                push!(cmd_args, "--force")
            end

            run(`$(Mycelia.CONDA_RUNNER) $cmd_args`)
        catch e
            @warn "Failed to import claMLST database: $e"
            # If strictly required, you might want to rethrow() here, 
            # but usually we proceed to check if search works.
        end
    end

    # Handle Decompression
    is_gzipped = endswith(genome_file, ".gz")
    input_path_to_use = genome_file
    temp_file = nothing

    try
        if is_gzipped
            # temp_file = joinpath(outdir, "tmp_" * base_name)
            temp_file = joinpath(outdir, base_name)
            open(temp_file, "w") do output_io
                run(pipeline(`gunzip -c $genome_file`, stdout = output_io))
            end
            input_path_to_use = temp_file
        end

        # Run Search
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n pymlst claMLST search $(db_dir) $(input_path_to_use)`

        open(output_tsv, "w") do io
            run(pipeline(cmd, stdout = io))
        end

        return output_tsv

    finally
        # Cleanup temp file
        if temp_file !== nothing && isfile(temp_file)
            rm(temp_file)
        end
    end
end

"""
    run_ectyper(genome_file::String; 
                outdir::String=joinpath(dirname(genome_file), "ectyper_output"),
                percent_identity::Int=90,
                percent_coverage::Int=50,
                verify_species::Bool=true,
                threads::Int=get_default_threads())

Run `ectyper` to predict *Escherichia coli* serotype and identify species.

# Arguments
- `genome_file`: Path to input FASTA file.
- `outdir`: Output directory.
- `percent_identity`: Minimum percent identity for antigen match (default: 90).
- `percent_coverage`: Minimum percent coverage for antigen match (default: 50).
- `verify_species`: Enable species verification module (default: true).

# Returns
- `String`: Path to the output directory containing `output.tsv`.
"""
function run_ectyper(genome_file::String;
        outdir::String = joinpath(dirname(genome_file), "ectyper_output"),
        percent_identity::Int = 90,
        percent_coverage::Int = 50,
        verify_species::Bool = true,
        threads::Int = get_default_threads())

    # ECTyper creates a file named 'output.tsv' inside the output directory
    expected_output = joinpath(outdir, "output.tsv")

    if isfile(expected_output)
        return expected_output
    end

    if !isfile(genome_file)
        error("Input genome file not found: $genome_file")
    end

    # Clean up directory if it exists but is incomplete (no output.tsv)
    if isdir(outdir)
        rm(outdir, recursive = true, force = true)
    end

    # Updated flags for new ectyper version:
    # Use specific flags for O-antigen (-opid/cov) and H-antigen (-hpid/cov)
    cmd_args = [
        "ectyper",
        "-i", genome_file,
        "-o", outdir,
        "-opid", string(percent_identity),  # O-antigen Identity
        "-opcov", string(percent_coverage), # O-antigen Coverage
        "-hpid", string(percent_identity),  # H-antigen Identity
        "-hpcov", string(percent_coverage), # H-antigen Coverage
        "-c", string(threads)
    ]

    if verify_species
        push!(cmd_args, "--verify")
    end

    Mycelia.add_bioconda_env("ectyper")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ectyper $cmd_args`)

    return outdir
end

"""
    run_ezclermont(genome_file::String; 
                   outdir::String=dirname(genome_file))

Run `ezclermont` to determine the Clermont phylogroup of an *E. coli* genome.

# Arguments
- `genome_file`: Path to input FASTA file.
- `outdir`: Output directory to save the results file.

# Details
EzClermont prints results to standard output. This function captures that output 
to a `.csv` file named after the input genome in the `outdir`.
It handles gzipped inputs by decompressing them temporarily.

# Returns
- `String`: Path to the generated CSV result file.
"""
function run_ezclermont(genome_file::String;
        outdir::String = dirname(genome_file))
    mkpath(outdir)

    base_name = replace(basename(genome_file), r"\.gz$" => "")
    base_name_no_ext = splitext(base_name)[1]
    result_file = joinpath(outdir, "$(base_name_no_ext)_ezclermont.csv")

    if isfile(result_file)
        return result_file
    end

    # Decompression handling
    is_gzipped = endswith(genome_file, ".gz")
    input_path_to_use = genome_file
    temp_file = nothing

    if is_gzipped
        temp_file = joinpath(outdir, "tmp_ez_" * base_name)
        open(temp_file, "w") do output_io
            run(pipeline(`gzip -dc $genome_file`, stdout = output_io))
        end
        input_path_to_use = temp_file
    end

    try
        Mycelia.add_bioconda_env("ezclermont")
        cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n ezclermont ezclermont $(input_path_to_use)`
        open(result_file, "w") do io
            run(pipeline(cmd, stdout = io))
        end
        return result_file
    finally
        if temp_file !== nothing && isfile(temp_file)
            rm(temp_file)
        end
    end
end

"""
    list_kleborate_modules()

Print the list of available Kleborate modules to standard output.
Useful for determining which modules to pass to `run_kleborate`.
"""
function list_kleborate_modules()
    Mycelia.add_bioconda_env("kleborate")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kleborate kleborate --list_modules`)
end

"""
    run_kleborate(assemblies::Vector{String}; 
                  outdir::String="kleborate_results",
                  preset::Union{String, Nothing}="kpsc",
                  modules::Union{Vector{String}, Nothing}=nothing,
                  all_amr::Bool=false,
                  trim_headers::Bool=false,
                  threads::Int=get_default_threads())

Run `kleborate` on a set of genome assemblies.

# Arguments
- `assemblies`: Vector of paths to assembly FASTA files.
- `outdir`: Output directory.
- `preset`: Analysis preset (default: "kpsc"). Options: "kpsc", "kosc", "escherichia". 
            Set to `nothing` if using `modules`.
- `modules`: Vector of specific module names to run (e.g., `["klebsiella_pneumo_complex__amr", "general__contig_stats"]`).
             If provided, this overrides `preset`.
- `trim_headers`: If true, adds `--trim_headers` to remove module prefixes from output columns.

# Returns
- `String`: Path to the output directory containing results.
"""
function run_kleborate(assemblies::Vector{String};
        outdir::String = "kleborate_results",
        preset::Union{String, Nothing} = "kpsc",
        modules::Union{Vector{String}, Nothing} = nothing,
        trim_headers::Bool = false,
        threads::Int = get_default_threads())
    mkpath(outdir)

    # # Check for existing output
    # expected_output = joinpath(outdir, "Kleborate_results.txt")
    # if isfile(expected_output)
    #     return expected_output
    # end
    # Heuristic check: if the directory contains any *_output.txt file, assume run is complete.
    # We can't predict the exact filename easily because it depends on the modules/species found.
    existing_outputs = filter(f -> endswith(f, "output.txt"), readdir(outdir))
    if !isempty(existing_outputs)
        # return (;outdir, existing_outputs)
        return outdir
    end

    if isempty(assemblies)
        error("No assembly files provided for Kleborate.")
    end

    cmd_args = ["kleborate", "--outdir", outdir]

    # Handle Modules vs Preset
    if modules !== nothing && !isempty(modules)
        # Use specific modules
        module_list = join(modules, ",")
        push!(cmd_args, "--modules", module_list)
    elseif preset !== nothing
        # Use preset
        push!(cmd_args, "--preset", preset)
    else
        error("You must provide either a `preset` or a list of `modules`.")
    end

    if trim_headers
        push!(cmd_args, "--trim_headers")
    end

    # Append assemblies
    append!(cmd_args, ["--assemblies"])
    append!(cmd_args, assemblies)

    Mycelia.add_bioconda_env("kleborate")
    @info "Running Kleborate on $(length(assemblies)) genomes..."
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n kleborate $cmd_args`)

    return outdir
end

"""
    run_kleborate(assembly_file::String; kwargs...)

Convenience wrapper for running Kleborate on a single file.
"""
function run_kleborate(assembly_file::String; kwargs...)
    return run_kleborate([assembly_file]; kwargs...)
end
