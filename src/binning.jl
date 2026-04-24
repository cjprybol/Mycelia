function _vamb_env_name()
    return "mycelia_vamb"
end

function _ensure_vamb_installed()
    env_name = _vamb_env_name()
    if !_check_conda_env_exists(env_name)
        run(`$(Mycelia.CONDA_RUNNER) create -y -n $(env_name) python=3.11 pip`)
    end

    vamb_ok = false
    try
        vamb_ok = success(`$(Mycelia.CONDA_RUNNER) run -n $(env_name) vamb --version`)
    catch
        vamb_ok = false
    end
    if !vamb_ok
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) python -m pip install --upgrade pip setuptools wheel`)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) python -m pip install vamb`)
    end

    vamb_ok = success(`$(Mycelia.CONDA_RUNNER) run -n $(env_name) vamb --version`)
    vamb_ok || error("VAMB installation failed in conda env: $(env_name)")
    return env_name
end

"""
    run_vamb(; contigs_fasta, depth_file, outdir, minfasta::Int=2000,
             threads::Int=get_default_threads(), extra_args::Vector{String}=String[])

Run VAMB to bin contigs using sequence composition and coverage.
Uses `vamb bin default`.

# Arguments
- `contigs_fasta::String`: FASTA file with assembled contigs
- `depth_file::String`: JGI depth table or VAMB abundance TSV
- `outdir::String`: Output directory for VAMB results
- `minfasta::Int`: Minimum contig length to include (passed to `-m`, default: 2000)
- `threads::Int`: Number of threads to use (default: all CPUs)
- `extra_args::Vector{String}`: Additional command-line arguments

# Returns
Named tuple with:
- `outdir`: Output directory
- `clusters_tsv`: Path to the VAMB clusters table (contig → bin), if produced
"""
function _vamb_abundance_tsv(depth_file::String; output_dir::Union{Nothing, String} = nothing)
    header = open(depth_file, "r") do io
        readline(io)
    end
    header_cols = split(chomp(header), '\t')
    header_lc = lowercase.(header_cols)
    is_vamb = !isempty(header_lc) && header_lc[1] == "contigname" &&
              !("contiglen" in header_lc) && !("totalavgdepth" in header_lc)
    if is_vamb
        return depth_file
    end

    target_dir = output_dir === nothing ? dirname(depth_file) : output_dir
    mkpath(target_dir)
    abundance_tsv = joinpath(target_dir, "vamb_abundance.tsv")
    if isfile(abundance_tsv)
        return abundance_tsv
    end

    open(depth_file, "r") do io
        header_cols = split(chomp(readline(io)), '\t')
        sample_cols = Int[]
        sample_names = String[]
        for (idx, name) in enumerate(header_cols)
            if idx <= 3
                continue
            end
            name_lc = lowercase(name)
            if occursin("var", name_lc)
                continue
            end
            push!(sample_cols, idx)
            push!(sample_names, name)
        end
        isempty(sample_cols) &&
            error("No sample depth columns found in JGI depth file: $(depth_file)")

        open(abundance_tsv, "w") do out
            write(out, "contigname")
            for name in sample_names
                write(out, '\t', name)
            end
            write(out, '\n')
            for line in eachline(io)
                isempty(line) && continue
                fields = split(line, '\t')
                write(out, fields[1])
                for idx in sample_cols
                    value = idx <= length(fields) ? fields[idx] : "0"
                    write(out, '\t', value)
                end
                write(out, '\n')
            end
        end
    end

    return abundance_tsv
end

function _find_first_matching_file(dir::String, patterns::Vector{Regex}; recursive::Bool = false)
    if recursive
        for (dirpath, _, filenames) in walkdir(dir)
            for name in sort(filenames)
                for pattern in patterns
                    if occursin(pattern, name)
                        return joinpath(dirpath, name)
                    end
                end
            end
        end
    else
        for name in sort(readdir(dir))
            path = joinpath(dir, name)
            isfile(path) || continue
            for pattern in patterns
                if occursin(pattern, name)
                    return path
                end
            end
        end
    end
    return nothing
end

function _find_first_matching_dir(dir::String, patterns::Vector{Regex}; recursive::Bool = false)
    if recursive
        for (dirpath, dirnames, _) in walkdir(dir)
            for name in sort(dirnames)
                for pattern in patterns
                    if occursin(pattern, name)
                        return joinpath(dirpath, name)
                    end
                end
            end
        end
    else
        for name in sort(readdir(dir))
            path = joinpath(dir, name)
            isdir(path) || continue
            for pattern in patterns
                if occursin(pattern, name)
                    return path
                end
            end
        end
    end
    return nothing
end

function run_vamb(; contigs_fasta::String, depth_file::String, outdir::String,
        minfasta::Int = 2000, threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[])
    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(depth_file) || error("Depth file not found: $(depth_file)")
    ispath(outdir) && error("VAMB output directory already exists: $(outdir)")
    mkpath(dirname(outdir))

    env_name = _ensure_vamb_installed()
    abundance_tsv = _vamb_abundance_tsv(depth_file)
    cmd_args = String[
    "vamb", "bin", "default",
    "--fasta", contigs_fasta,
    "--abundance_tsv", abundance_tsv,
    "--outdir", outdir,
    "-m", string(minfasta)
]
    if !any(arg -> arg == "-p" || startswith(arg, "--threads"), extra_args)
        append!(cmd_args, ["-p", string(threads)])
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) $(cmd_args)`)

    clusters_tsv = _find_first_matching_file(outdir, [
        r"vae_clusters.*\.tsv$", r"clusters.*\.tsv$"])
    return (; outdir, clusters_tsv)
end

"""
    run_metabat2(; contigs_fasta, depth_file, outdir,
                 min_contig::Int=1500, threads::Int=get_default_threads(),
                 seed::Int=42, extra_args::Vector{String}=String[])

Run MetaBAT2 for metagenomic binning.

# Arguments
- `contigs_fasta::String`: Assembled contigs FASTA
- `depth_file::String`: Depth/coverage table from jgi_summarize_bam_contig_depths
- `outdir::String`: Output directory (will contain bin FASTAs with prefix `bin`)
- `min_contig::Int`: Minimum contig length (default: 1500)
- `threads::Int`: Number of threads (default: all CPUs)
- `seed::Int`: Random seed (default: 42)
- `extra_args::Vector{String}`: Additional MetaBAT2 arguments

# Returns
Named tuple with:
- `outdir`: Output directory
- `bins_prefix`: Output prefix used for generated bins
"""
function run_metabat2(; contigs_fasta::String, depth_file::String, outdir::String,
        min_contig::Int = 1500, threads::Int = get_default_threads(), seed::Int = 42,
        extra_args::Vector{String} = String[])
    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(depth_file) || error("Depth file not found: $(depth_file)")
    mkpath(outdir)

    add_bioconda_env("metabat2")

    bins_prefix = joinpath(outdir, "bin")
    cmd_args = String[
    "metabat2",
    "-i", contigs_fasta,
    "-a", depth_file,
    "-o", bins_prefix,
    "-m", string(min_contig),
    "-t", string(threads),
    "-s", string(seed)
]
    append!(cmd_args, extra_args)

    run(`$(CONDA_RUNNER) run --live-stream -n metabat2 $(cmd_args)`)
    return (; outdir, bins_prefix)
end

"""
    run_metacoag(; contigs_fasta, assembly_graph, mapping_file, outdir,
                 assembler::String="custom", threads::Int=get_default_threads(),
                 extra_args::Vector{String}=String[])

Run MetaCoAG, which integrates assembly graph structure and coverage for binning.

# Arguments
- `contigs_fasta::String`: Assembled contigs FASTA
- `assembly_graph::String`: Assembly graph file (GFA/GFA1)
- `mapping_file::String`: Read mapping/coverage file
- `outdir::String`: Output directory
- `assembler::String`: Assembler name (default: "custom")
- `threads::Int`: Number of threads (default: all CPUs)
- `extra_args::Vector{String}`: Additional MetaCoAG arguments

# Returns
Named tuple with:
- `outdir`: Output directory
- `bins_dir`: Directory containing bin FASTAs (if produced)
- `bins_tsv`: Contig → bin assignments table (if produced)
"""
function run_metacoag(; contigs_fasta::String, assembly_graph::String,
        mapping_file::String, outdir::String, assembler::String = "custom",
        threads::Int = get_default_threads(), extra_args::Vector{String} = String[])
    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(assembly_graph) || error("Assembly graph not found: $(assembly_graph)")
    isfile(mapping_file) || error("Mapping/coverage file not found: $(mapping_file)")
    mkpath(outdir)

    add_bioconda_env("metacoag")

    bins_dir = joinpath(outdir, "bins")
    bins_tsv = _find_first_matching_file(outdir, [r"bin.*\.tsv$"])
    cmd_args = String[
    "metacoag",
    "--assembler", assembler,
    "--graph", assembly_graph,
    "--contigs", contigs_fasta,
    "--abundance", mapping_file,
    "--output", outdir,
    "--nthreads", string(threads)
]
    append!(cmd_args, extra_args)

    run(`$(CONDA_RUNNER) run --live-stream -n metacoag $(cmd_args)`)
    return (; outdir, bins_dir, bins_tsv)
end

"""
    run_comebin(; contigs_fasta, bam_path, outdir,
                 views::Int=6, threads::Int=get_default_threads(),
                 temperature::Union{Nothing,Float64}=nothing,
                 embedding_size::Int=2048, coverage_embedding_size::Int=2048,
                 batch_size::Int=1024, extra_args::Vector{String}=String[])

Run COMEBin (constraint-based binning).

# Arguments
- `contigs_fasta::String`: Contig FASTA
- `bam_path::String`: Path to BAM file or directory containing BAMs
- `outdir::String`: Output directory
- `views::Int`: Number of contrastive learning views
- `threads::Int`: Number of threads (default: all CPUs)
- `temperature::Union{Nothing,Float64}`: Optional loss temperature (defaults to COMEBin heuristic)
- `embedding_size::Int`: Embedding size for combine network
- `coverage_embedding_size::Int`: Embedding size for coverage network
- `batch_size::Int`: Training batch size
- `extra_args::Vector{String}`: Additional COMEBin arguments

# Returns
Named tuple with:
- `outdir`: Output directory
- `bins_dir`: Directory containing bin FASTAs (if produced)
- `bins_tsv`: Contig → bin assignments table (if produced)
"""
function run_comebin(; contigs_fasta::String, bam_path::String, outdir::String,
        views::Int = 6, threads::Int = get_default_threads(),
        temperature::Union{Nothing, Float64} = nothing,
        embedding_size::Int = 2048, coverage_embedding_size::Int = 2048,
        batch_size::Int = 1024, extra_args::Vector{String} = String[])
    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    if !(isfile(bam_path) || isdir(bam_path))
        error("BAM path not found: $(bam_path)")
    end
    if isdir(bam_path)
        bam_files = filter(name -> endswith(lowercase(name), ".bam"), readdir(bam_path))
        isempty(bam_files) && error("No BAM files found in directory: $(bam_path)")
    end
    mkpath(outdir)

    add_bioconda_env("comebin")

    cmd_args = String[
    "run_comebin.sh",
    "-a", contigs_fasta,
    "-o", outdir,
    "-p", bam_path,
    "-n", string(views),
    "-t", string(threads),
    "-e", string(embedding_size),
    "-c", string(coverage_embedding_size),
    "-b", string(batch_size)
]
    if temperature !== nothing
        push!(cmd_args, "-l")
        push!(cmd_args, string(temperature))
    end
    append!(cmd_args, extra_args)

    run(`$(CONDA_RUNNER) run --live-stream -n comebin $(cmd_args)`)
    bins_dir = _find_first_matching_dir(outdir, [r"bins$"]; recursive = true)
    bins_tsv = _find_first_matching_file(outdir, [r"bins.*\.tsv$"]; recursive = true)
    return (; outdir, bins_dir, bins_tsv)
end

"""
    run_drep_dereplicate(; genomes, outdir, completeness_threshold=nothing,
                         contamination_threshold=nothing, ani_threshold=0.99,
                         threads::Int=get_default_threads(), extra_args::Vector{String}=String[])

Run dRep for MAG dereplication.

# Arguments
- `genomes::Vector{String}`: Genome/MAG FASTA files to dereplicate
- `outdir::String`: Output directory
- `completeness_threshold::Union{Nothing, Float64}`: Minimum completeness to keep a genome
- `contamination_threshold::Union{Nothing, Float64}`: Maximum contamination allowed
- `ani_threshold::Float64`: ANI threshold for secondary clustering (default: 0.99)
- `threads::Int`: Number of threads (default: all CPUs)
- `extra_args::Vector{String}`: Additional dRep arguments

# Returns
Named tuple with:
- `outdir`: Output directory
- `winning_genomes`: Path to winning genomes CSV (if produced)
"""
function run_drep_dereplicate(; genomes::Vector{String}, outdir::String,
        completeness_threshold::Union{Nothing, Float64} = nothing,
        contamination_threshold::Union{Nothing, Float64} = nothing,
        ani_threshold::Float64 = 0.99,
        threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[])
    isempty(genomes) && error("No genomes provided to dRep")
    for genome in genomes
        isfile(genome) || error("Genome file not found: $(genome)")
    end
    mkpath(outdir)

    add_bioconda_env("drep")

    cmd_args = String[
    "dRep", "dereplicate", outdir,
    "-g"
]
    append!(cmd_args, genomes)
    append!(cmd_args, ["-p", string(threads), "-sa", string(ani_threshold)])

    if !isnothing(completeness_threshold)
        push!(cmd_args, "-comp")
        push!(cmd_args, string(completeness_threshold))
    end

    if !isnothing(contamination_threshold)
        push!(cmd_args, "-con")
        push!(cmd_args, string(contamination_threshold))
    end

    append!(cmd_args, extra_args)

    run(`$(CONDA_RUNNER) run --live-stream -n drep $(cmd_args)`)

    winning_genomes = _find_first_matching_file(
        outdir,
        [r"Widb\.csv$", r"dereplicated_genomes\.csv$"];
        recursive = true
    )
    return (; outdir, winning_genomes)
end

"""
    run_skder(; genomes, outdir, ani_threshold=99.5, af_threshold=50.0,
              mode=:dynamic, filter_mge=false, determine_clusters=true,
              symlink=true, threads::Int=get_default_threads(),
              extra_args::Vector{String}=String[])

Run skDER for ANI-based reference dereplication.

skDER is purpose-built for reducing redundancy in reference panels used for
metagenomic competitive mapping. `dynamic` mode approximates single-linkage
clustering and minimizes the representative count (best for mapping panels);
`greedy` mode strictly enforces the user ANI/AF cutoffs (best for comparative
genomics). See Salamzade & Kalan, *Microbial Genomics* 2025
(https://doi.org/10.1099/mgen.0.001438).

The wrapper stages input genome paths into a temporary directory of symlinks
before invoking skDER, which sidesteps the `ARG_MAX` command-line limit when
dereplicating thousands of genomes. It also sanitizes `PYTHONPATH` and sets
`PYTHONNOUSERSITE=1` in the subprocess env to prevent user-site shadowing of
the conda env's numpy (which breaks skDER's Python 3.13 install on hosts
where `~/.local/lib/python3.9/site-packages/numpy` exists).

# Arguments
- `genomes::Vector{String}`: Genome FASTA files to dereplicate. Must be
  non-empty; every path must exist.
- `outdir::String`: Output directory; created if missing. If it already
  exists skDER may fail (it refuses to overwrite); callers are responsible
  for removing stale directories on retry.
- `ani_threshold::Float64=99.5`: ANI percentage cutoff (skDER `-i`).
  Rodriguez-R et al. 2024 within-species knee for E. coli.
- `af_threshold::Float64=50.0`: Aligned-fraction percentage cutoff (skDER
  `-f`). skDER's default; use 95.0 for tighter strain-level dereplication
  tied to gene-content similarity.
- `mode::Symbol=:dynamic`: `:dynamic` (metagenomic mapping panel) or
  `:greedy` (comparative genomics).
- `filter_mge::Bool=false`: Pass `--filter-mge` (`-fm`) to mask mobile
  genetic elements (plasmids, prophages) before ANI/AF computation using
  PhiSpy. Useful for UPEC and other MGE-heavy taxa.
- `determine_clusters::Bool=true`: Pass `--determine-clusters` (`-n`) to
  produce the secondary cluster assignment TSV that `parse_skder_clusters`
  expects. Keep enabled unless you only need the representative set.
- `symlink::Bool=true`: Pass `--symlink` (`-l`) so the representative
  directory contains symlinks to the input files rather than copies.
- `threads::Int`: Thread count (skDER `-c`; default `get_default_threads()`).
- `extra_args::Vector{String}`: Additional CLI arguments forwarded to skDER.

# Returns
Named tuple with:
- `outdir`: The output directory.
- `representatives_dir`: Path to the directory containing representative
  FASTAs (or `nothing` if not produced).
- `representatives`: `Vector{String}` of representative FASTA paths
  (resolved to absolute paths when `symlink=true`).
- `cluster_info`: Path to the skDER cluster-info TSV mapping redundant
  genomes to their representative (`nothing` if not produced).
"""
function run_skder(; genomes::Vector{String}, outdir::String,
        ani_threshold::Float64 = 99.5,
        af_threshold::Float64 = 50.0,
        mode::Symbol = :dynamic,
        filter_mge::Bool = false,
        determine_clusters::Bool = true,
        symlink::Bool = true,
        threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[])
    isempty(genomes) && error("No genomes provided to skDER")
    seen_basenames = Set{String}()
    for genome in genomes
        isfile(genome) || error("Genome file not found: $(genome)")
        bn = basename(genome)
        bn in seen_basenames &&
            error("Duplicate input basename: $(bn); skDER staging requires unique basenames")
        push!(seen_basenames, bn)
    end
    mode in (:dynamic, :greedy) ||
        error("Invalid skDER mode: $(mode); expected :dynamic or :greedy")
    0.0 <= ani_threshold <= 100.0 ||
        error("Invalid skDER ani_threshold: $(ani_threshold); expected in [0, 100]")
    0.0 <= af_threshold <= 100.0 ||
        error("Invalid skDER af_threshold: $(af_threshold); expected in [0, 100]")
    mkpath(outdir)

    _ensure_skder_env()

    # Stage input genomes as symlinks in a single directory to avoid ARG_MAX
    # overflow when dereplicating thousands of genomes.
    staging_dir = mktempdir(prefix = "skder_inputs_")
    try
        for genome in genomes
            Base.Filesystem.symlink(abspath(genome),
                joinpath(staging_dir, basename(genome)))
        end

        cmd_args = String["skder",
        "-g", staging_dir,
        "-o", outdir,
        "-i", string(ani_threshold),
        "-f", string(af_threshold),
        "-d", string(mode),
        "-c", string(threads)]
        if filter_mge
            push!(cmd_args, "--filter-mge")
        end
        if determine_clusters
            push!(cmd_args, "--determine-clusters")
        end
        if symlink
            push!(cmd_args, "--symlink")
        end
        append!(cmd_args, extra_args)

        # Sanitize Python env: PYTHONNOUSERSITE=1 disables ~/.local/lib user
        # site-packages which shadow the conda env's numpy; PYTHONPATH= drops
        # any inherited PYTHONPATH that would poison `import numpy`.
        cmd = addenv(`$(CONDA_RUNNER) run --live-stream -n skder $(cmd_args)`,
            "PYTHONNOUSERSITE" => "1",
            "PYTHONPATH" => "")
        run(cmd)
    finally
        rm(staging_dir; recursive = true, force = true)
    end

    representatives_dir = _find_first_matching_dir(
        outdir,
        [r"Dereplicated_Representative_Genomes$", r"Representative_Genomes$",
            r"representatives?$"i];
        recursive = true
    )
    representatives = String[]
    if !isnothing(representatives_dir) && isdir(representatives_dir)
        for f in readdir(representatives_dir; join = true)
            if isfile(f) && (endswith(f, ".fa") || endswith(f, ".fna") ||
                endswith(f, ".fasta") || endswith(f, ".fa.gz") ||
                endswith(f, ".fna.gz") || endswith(f, ".fasta.gz"))
                push!(representatives, f)
            end
        end
        sort!(representatives)
    end
    cluster_info = _find_first_matching_file(
        outdir,
        [r"skDER_Clustering\.txt$", r"skDER_Cluster_Info\.tsv$",
            r"Cluster_Info\.tsv$"];
        recursive = true
    )

    return (; outdir, representatives_dir, representatives, cluster_info)
end

"""
    _ensure_skder_env()

Ensure the `skder` conda env exists and is pinned to `numpy<2`. skDER 1.3.4
crashes in its multiprocessing N50 step under numpy 2.x with
`only 0-dimensional arrays can be converted to Python scalars`. Upstream fix
is expected but unreleased as of April 2026; pin aggressively here and relax
once skDER tags a numpy-2-compatible release.
"""
function _ensure_skder_env()
    env_name = "skder"
    if check_bioconda_env_is_installed(env_name)
        numpy_ok = _conda_env_exec_ok(env_name,
            ["python",
                "-c",
                "import numpy as n; import sys; " *
                "sys.exit(0 if int(n.__version__.split('.')[0]) < 2 else 1)"])
        numpy_ok && return env_name
        @info "Existing skder env has incompatible numpy; rebuilding with numpy<2 pin."
        run(`$(CONDA_RUNNER) env remove -n $(env_name) -y`)
    end
    run(`$(CONDA_RUNNER) create -c conda-forge -c bioconda --strict-channel-priority -y -n $(env_name) skder "numpy<2"`)
    return env_name
end

"""
    parse_skder_clusters(file::String)::DataFrames.DataFrame

Parse a skDER cluster-info TSV into a DataFrame with columns
`genome`, `representative`, and `cluster_id`. Column names in the input file
may vary by skDER version; the parser accepts common aliases.
"""
function parse_skder_clusters(file::String)::DataFrames.DataFrame
    isfile(file) || error("skDER cluster-info file not found: $(file)")
    df = DataFrames.DataFrame(CSV.File(file; delim = '\t', ignorerepeated = true))
    names_map = Dict(string(col) => col for col in DataFrames.names(df))

    genome_key = findfirst(k -> haskey(names_map, k),
        ("genome", "Genome", "query", "Query"))
    rep_aliases = ("nearest_representative_genome",
        "representative", "Representative",
        "Representative_Genome", "representative_genome")
    rep_key = findfirst(k -> haskey(names_map, k), rep_aliases)
    cluster_key = findfirst(k -> haskey(names_map, k),
        ("cluster_id", "Cluster_ID", "cluster", "Cluster"))

    isnothing(genome_key) &&
        error("No genome column (genome/Genome/query) found in $(file)")
    isnothing(rep_key) &&
        error("No representative column found in $(file)")

    selection = [
        names_map[("genome", "Genome", "query", "Query")[genome_key]] => :genome,
        names_map[rep_aliases[rep_key]] => :representative
    ]
    if !isnothing(cluster_key)
        push!(selection,
            names_map[("cluster_id", "Cluster_ID", "cluster", "Cluster")[cluster_key]] => :cluster_id)
    end
    df_out = DataFrames.select(df, selection...)
    df_out.genome = string.(df_out.genome)
    df_out.representative = string.(df_out.representative)
    if "cluster_id" in names(df_out)
        df_out.cluster_id = string.(df_out.cluster_id)
    end
    return df_out
end

# -----------------------------------------------------------------------------
# Parsers
# -----------------------------------------------------------------------------

"""
    parse_bin_assignments(file::String; contig_col::AbstractString=\"contig\", bin_col::AbstractString=\"bin\")

Parse a generic contig → bin assignment table into a DataFrame.

The table is expected to be tab-delimited and contain at least the columns
named by `contig_col` and `bin_col`.
"""
function parse_bin_assignments(file::String; contig_col::AbstractString = "contig",
        bin_col::AbstractString = "bin")::DataFrames.DataFrame
    isfile(file) || error("Assignment file not found: $(file)")
    df = DataFrames.DataFrame(CSV.File(file; delim = '\t', ignorerepeated = true))
    names_map = Dict(string(col) => col for col in DataFrames.names(df))

    haskey(names_map, contig_col) || error("Column $(contig_col) not found in $(file)")
    haskey(names_map, bin_col) || error("Column $(bin_col) not found in $(file)")

    contig_sym = names_map[contig_col]
    bin_sym = names_map[bin_col]

    return DataFrames.select(df, contig_sym => :contig, bin_sym => :bin)
end

"""
    parse_drep_clusters(file::String)::DataFrames.DataFrame

Parse dRep dereplication cluster assignments (CSV) into a DataFrame
with columns `genome`, `secondary_cluster`, and `representative`.
"""
function parse_drep_clusters(file::String)::DataFrames.DataFrame
    isfile(file) || error("dRep cluster file not found: $(file)")
    df = DataFrames.DataFrame(CSV.File(file; delim = ',', ignorerepeated = true))
    names_map = Dict(string(col) => col for col in DataFrames.names(df))

    required = ("genome", "secondary_cluster", "representative")
    for col in required
        haskey(names_map, col) || error("Column $(col) not found in $(file)")
    end

    df_out = DataFrames.select(df,
        names_map["genome"] => :genome,
        names_map["secondary_cluster"] => :secondary_cluster,
        names_map["representative"] => :representative)

    # Normalize to String columns for consistent downstream handling
    df_out.genome = string.(df_out.genome)
    df_out.secondary_cluster = string.(df_out.secondary_cluster)
    df_out.representative = string.(df_out.representative)

    return df_out
end

# -----------------------------------------------------------------------------
# Additional binning tool wrappers (Taxometer, TaxVAMB, GenomeFace) and MAG merging (MAGmax)
# -----------------------------------------------------------------------------

"""
    run_taxometer(; contigs_fasta, depth_file, taxonomy_file, outdir,
                  threads=get_default_threads(), extra_args=String[])

Run Taxometer for contig binning using coverage and taxonomy priors.
Uses `vamb taxometer`.

# Arguments
- `contigs_fasta::String`: FASTA file with assembled contigs
- `depth_file::String`: JGI depth table or VAMB abundance TSV
- `taxonomy_file::String`: Taxonomy TSV
- `outdir::String`: Output directory
- `threads::Int`: Number of threads (default: all CPUs)
- `extra_args::Vector{String}`: Additional command-line arguments
"""
function run_taxometer(; contigs_fasta::String, depth_file::String, taxonomy_file::String,
        outdir::String, threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[])
    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(depth_file) || error("Depth file not found: $(depth_file)")
    isfile(taxonomy_file) || error("Taxonomy file not found: $(taxonomy_file)")
    ispath(outdir) && error("Taxometer output directory already exists: $(outdir)")
    mkpath(dirname(outdir))

    env_name = _ensure_vamb_installed()
    abundance_tsv = _vamb_abundance_tsv(depth_file)
    cmd_args = String[
    "vamb", "taxometer",
    "--fasta", contigs_fasta,
    "--abundance_tsv", abundance_tsv,
    "--taxonomy", taxonomy_file,
    "--outdir", outdir
]
    if !any(arg -> arg == "-p" || startswith(arg, "--threads"), extra_args)
        append!(cmd_args, ["-p", string(threads)])
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) $(cmd_args)`)
    return (; outdir)
end

"""
    run_taxvamb(; contigs_fasta, depth_file, taxonomy_file, outdir,
                threads=get_default_threads(), extra_args=String[])

Run TaxVAMB hybrid binning (VAMB with taxonomy priors).
Uses `vamb bin taxvamb`.

# Arguments
- `contigs_fasta::String`: FASTA file with assembled contigs
- `depth_file::String`: JGI depth table or VAMB abundance TSV
- `taxonomy_file::String`: Taxonomy TSV
- `outdir::String`: Output directory
- `threads::Int`: Number of threads (default: all CPUs)
- `extra_args::Vector{String}`: Additional command-line arguments
"""
function run_taxvamb(; contigs_fasta::String, depth_file::String, taxonomy_file::String,
        outdir::String, threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[])
    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(depth_file) || error("Depth file not found: $(depth_file)")
    isfile(taxonomy_file) || error("Taxonomy file not found: $(taxonomy_file)")
    ispath(outdir) && error("TaxVAMB output directory already exists: $(outdir)")
    mkpath(dirname(outdir))

    env_name = _ensure_vamb_installed()
    abundance_tsv = _vamb_abundance_tsv(depth_file)
    cmd_args = String[
    "vamb", "bin", "taxvamb",
    "--fasta", contigs_fasta,
    "--abundance_tsv", abundance_tsv,
    "--taxonomy", taxonomy_file,
    "--outdir", outdir
]
    if !any(arg -> arg == "-p" || startswith(arg, "--threads"), extra_args)
        append!(cmd_args, ["-p", string(threads)])
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) $(cmd_args)`)
    clusters_tsv = _find_first_matching_file(outdir, [
        r"vae_clusters.*\.tsv$", r"clusters.*\.tsv$"])
    return (; outdir, clusters_tsv)
end

"""
    run_genomeface(; contigs_fasta, coverage_table, outdir, threads=get_default_threads(), extra_args=String[])

Run GenomeFace for binning with marker- and coverage-aware clustering.

Note: This wrapper is disabled until the `genomeface` executable can be located again.
"""
function run_genomeface(; contigs_fasta::String, coverage_table::String, outdir::String,
        threads::Int = get_default_threads(), extra_args::Vector{String} = String[])
    error("GenomeFace wrapper disabled: genomeface executable unavailable. Re-enable when the binary is found.")
end

"""
    run_magmax_merge(; bins_dirs, outdir, threads=get_default_threads(), extra_args=String[])

Merge MAGs/bins from multiple binners using MAGmax.
"""
function run_magmax_merge(; bins_dirs::Vector{String}, outdir::String,
        threads::Int = get_default_threads(), extra_args::Vector{String} = String[])
    isempty(bins_dirs) && error("No bins directories provided to MAGmax")
    for dir in bins_dirs
        isdir(dir) || error("Bins directory not found: $(dir)")
    end
    mkpath(outdir)

    add_bioconda_env("magmax")

    bins_input_dir = joinpath(outdir, "magmax_bins_input")
    mkpath(bins_input_dir)
    for dir in bins_dirs
        for name in readdir(dir)
            src = joinpath(dir, name)
            isfile(src) || continue
            dest = joinpath(bins_input_dir, "$(basename(dir))__$(name)")
            cp(src, dest; force = true)
        end
    end

    cmd_args = String["magmax", "--bindir", bins_input_dir, "--threads", string(threads)]
    if !any(arg -> arg in ("-f", "--format"), extra_args)
        exts = String[]
        for name in readdir(bins_input_dir)
            ext = lowercase(replace(splitext(name)[2], "." => ""))
            isempty(ext) && continue
            push!(exts, ext)
        end
        unique_exts = sort(unique(exts))
        if !isempty(unique_exts)
            append!(cmd_args, ["--format", unique_exts[1]])
        end
    end
    append!(cmd_args, extra_args)

    cmd = Cmd(`$(CONDA_RUNNER) run --live-stream -n magmax $(cmd_args)`; dir = outdir)
    run(cmd)
    return (; outdir, bins_input_dir)
end
