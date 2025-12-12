# Classification Tools for Metagenomic Analysis
#
# This file contains wrappers for metagenomic classification tools:
# - sourmash: MinHash-based sequence search and classification
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
        k_sizes::Vector{Int}=[21, 31, 51],
        scaled::Int=1000,
        molecule::String="dna",
        singleton::Bool=false,
        name::Union{Nothing, String}=nothing)

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
        threshold::Float64=0.1,
        k_size::Int=31,
        best_only::Bool=false,
        num_results::Union{Nothing, Int}=nothing)

    isfile(query_sig) || error("Query signature not found: $(query_sig)")
    (isfile(database_sig) || isdir(database_sig)) || error("Database not found: $(database_sig)")

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
        k_size::Int=31,
        threshold_bp::Int=50000)

    isfile(query_sig) || error("Query signature not found: $(query_sig)")
    (isfile(database_sig) || isdir(database_sig)) || error("Database not found: $(database_sig)")

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

# ============================================================================
# MetaPhlAn - Marker-based taxonomic profiling
# ============================================================================

"""
    run_metaphlan(;input_file, outdir, input_type="fastq",
                  nprocs=get_default_threads(), db_dir=nothing,
                  unknown_estimation=true, stat_q=0.2,
                  long_reads=false)

Run MetaPhlAn for marker-based taxonomic profiling.

# Arguments
- `input_file::String`: Input file (FASTQ, FASTA, or BAM)
- `outdir::String`: Output directory
- `input_type::String`: Input type - "fastq", "fasta", "mapout", "sam" (default: "fastq")
- `nprocs::Int`: Number of threads (default: get_default_threads())
- `db_dir::Union{Nothing, String}`: Path to database directory (passed to MetaPhlAn `--db_dir`)
- `unknown_estimation::Bool`: Include unknown/unclassified fraction estimation (default: true)
- `stat_q::Float64`: Quantile for robust average (default: 0.2)
- `long_reads::Bool`: Use long-read mode (MetaPhlAn `--long_reads`, uses minimap2)

# Returns
Named tuple with:
- `outdir`: Output directory path
- `profile_txt`: Path to taxonomic profile output
- `mapout`: Path to mapping output (also returned as `bowtie2_out` for compatibility)

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
        input_type::String="fastq",
        nprocs::Int=get_default_threads(),
        db_dir::Union{Nothing, String}=nothing,
        unknown_estimation::Bool=true,
        stat_q::Float64=0.2,
        long_reads::Bool=false)

    isfile(input_file) || error("Input file not found: $(input_file)")
    input_type in ["fastq", "fasta", "mapout", "sam"] || error("Invalid input_type: $(input_type)")

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

        if !unknown_estimation
            push!(cmd_args, "--skip_unclassified_estimation")
        end

        if long_reads
            push!(cmd_args, "--long_reads")
        end

        run(Cmd(cmd_args))
    end

    return (; outdir, profile_txt, mapout=mapout_file)
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
            clade_name=String[],
            taxid=Union{Int, Missing}[],
            relative_abundance=Float64[],
            additional_species=Union{String, Missing}[]
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
        clade_name=clade_names,
        taxid=taxids,
        relative_abundance=abundances,
        additional_species=additional_species
    )
end

# ============================================================================
# Metabuli - Fast metagenomic classification
# ============================================================================

"""
    run_metabuli_classify(;input_files, outdir, database_path=nothing,
                          seq_mode="1", threads=get_default_threads(),
                          min_score=nothing, min_sp_score=nothing)

Run Metabuli for fast metagenomic classification.

# Arguments
- `input_files::Vector{String}`: Input sequence files (FASTA/FASTQ)
- `outdir::String`: Output directory
- `database_path::Union{Nothing, String}`: Path to Metabuli database. If not provided,
  checks METABULI_DB environment variable, then falls back to `\$HOME/workspace/metabuli`
- `seq_mode::String`: Sequence mode - "1" for single-end, "2" for paired-end (default: "1")
- `threads::Int`: Number of threads (default: get_default_threads())
- `min_score::Union{Nothing, Float64}`: Minimum score threshold
- `min_sp_score::Union{Nothing, Float64}`: Minimum species-level score

# Returns
Named tuple with:
- `outdir`: Output directory path
- `report_file`: Path to classification report
- `classifications_file`: Path to per-read classifications

# Example
```julia
result = run_metabuli_classify(
    input_files=["reads.fq.gz"],
    outdir="metabuli_output"
)
```
"""
function run_metabuli_classify(;
        input_files::Vector{String},
        outdir::String,
        database_path::Union{Nothing, String}=nothing,
        seq_mode::String="1",
        threads::Int=get_default_threads(),
        min_score::Union{Nothing, Float64}=nothing,
        min_sp_score::Union{Nothing, Float64}=nothing)

    # Validate inputs
    for f in input_files
        isfile(f) || error("Input file not found: $(f)")
    end

    # Determine database path
    if database_path === nothing
        database_path = get(ENV, "METABULI_DB", nothing)
        if database_path === nothing
            home_db = joinpath(homedir(), "workspace", "metabuli")
            if isdir(home_db)
                database_path = home_db
            else
                error("No database_path provided and no database found at default location ($(home_db)). " *
                      "Set METABULI_DB environment variable or provide database_path parameter.")
            end
        end
    end
    isdir(database_path) || error("Database directory not found: $(database_path)")
    seq_mode in ["1", "2", "3"] || error("Invalid seq_mode: $(seq_mode). Use '1' for single-end, '2' for paired-end interleaved, '3' for paired-end separate")

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
            "--threads", string(threads)
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
        threads::Int=get_default_threads(),
        split_num::Int=4096)

    isfile(reference_fasta) || error("Reference FASTA not found: $(reference_fasta)")
    isdir(taxonomy_dir) || error("Taxonomy directory not found: $(taxonomy_dir)")

    Mycelia.add_bioconda_env("metabuli")
    mkpath(outdir)

    # Check for required taxonomy files
    for required in ["names.dmp", "nodes.dmp"]
        f = joinpath(taxonomy_dir, required)
        isfile(f) || error("Required taxonomy file not found: $(f)")
    end

    db_info_file = joinpath(outdir, "db.info")

    if !isfile(db_info_file)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n metabuli metabuli build
            $(reference_fasta)
            $(taxonomy_dir)
            $(outdir)
            --threads $(threads)
            --split-num $(split_num)`)
    end

    return (; database_path=outdir)
end

"""
    parse_metabuli_report(report_file::String)::DataFrames.DataFrame

Parse a Metabuli classification report into a DataFrame.

# Arguments
- `report_file::String`: Path to Metabuli report.tsv output

# Returns
DataFrame with classification statistics per taxon
"""
function parse_metabuli_report(report_file::String)::DataFrames.DataFrame
    isfile(report_file) || error("Report file not found: $(report_file)")

    # Metabuli report is tab-separated
    df = DataFrames.DataFrame(CSV.File(report_file; delim='\t', header=true))

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
    isfile(classifications_file) || error("Classifications file not found: $(classifications_file)")

    # Read the classifications file
    df = DataFrames.DataFrame(CSV.File(classifications_file; delim='\t', header=false))

    # Rename columns based on typical metabuli output format
    if DataFrames.ncol(df) >= 3
        DataFrames.rename!(df, 1 => :classified, 2 => :read_id, 3 => :taxid)
    end

    return df
end
