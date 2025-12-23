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
- `depth_file::String`: JGI-style depth/coverage table
- `outdir::String`: Output directory for VAMB results
- `minfasta::Int`: Minimum contig length to include (default: 2000)
- `threads::Int`: Number of threads to use (default: all CPUs)
- `extra_args::Vector{String}`: Additional command-line arguments

# Returns
Named tuple with:
- `outdir`: Output directory
- `clusters_tsv`: Path to the VAMB clusters table (contig → bin)
"""
function run_vamb(; contigs_fasta::String, depth_file::String, outdir::String,
        minfasta::Int=2000, threads::Int=get_default_threads(),
        extra_args::Vector{String}=String[])

    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(depth_file) || error("Depth file not found: $(depth_file)")
    mkpath(outdir)

    env_name = _ensure_vamb_installed()
    cmd_args = String[
        "vamb", "bin", "default",
        "--fasta", contigs_fasta,
        "--jgi", depth_file,
        "--outdir", outdir,
        "--minfasta", string(minfasta),
    ]
    if !any(arg -> arg == "-p" || startswith(arg, "--threads"), extra_args)
        append!(cmd_args, ["-p", string(threads)])
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) $(cmd_args)`)

    clusters_tsv = joinpath(outdir, "clusters.tsv")
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
        min_contig::Int=1500, threads::Int=get_default_threads(), seed::Int=42,
        extra_args::Vector{String}=String[])

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
        "-s", string(seed),
    ]
    append!(cmd_args, extra_args)

    run(`$(CONDA_RUNNER) run --live-stream -n metabat2 $(cmd_args)`)
    return (; outdir, bins_prefix)
end

"""
    run_metacoag(; contigs_fasta, assembly_graph, mapping_file, outdir,
                 threads::Int=get_default_threads(), extra_args::Vector{String}=String[])

Run MetaCoAG, which integrates assembly graph structure and coverage for binning.

# Arguments
- `contigs_fasta::String`: Assembled contigs FASTA
- `assembly_graph::String`: Assembly graph file (GFA/GFA1)
- `mapping_file::String`: Read mapping/coverage file
- `outdir::String`: Output directory
- `threads::Int`: Number of threads (default: all CPUs)
- `extra_args::Vector{String}`: Additional MetaCoAG arguments

# Returns
Named tuple with:
- `outdir`: Output directory
- `bins_tsv`: Contig → bin assignments table
"""
function run_metacoag(; contigs_fasta::String, assembly_graph::String,
        mapping_file::String, outdir::String, threads::Int=get_default_threads(),
        extra_args::Vector{String}=String[])

    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(assembly_graph) || error("Assembly graph not found: $(assembly_graph)")
    isfile(mapping_file) || error("Mapping/coverage file not found: $(mapping_file)")
    mkpath(outdir)

    add_bioconda_env("metacoag")

    bins_tsv = joinpath(outdir, "bins.tsv")
    cmd_args = String[
        "metacoag",
        "-a", contigs_fasta,
        "-g", assembly_graph,
        "-m", mapping_file,
        "-o", outdir,
        "-t", string(threads),
    ]
    append!(cmd_args, extra_args)

    run(`$(CONDA_RUNNER) run --live-stream -n metacoag $(cmd_args)`)
    return (; outdir, bins_tsv)
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
- `bins_tsv`: Contig → bin assignments table
"""
function run_comebin(; contigs_fasta::String, bam_path::String, outdir::String,
        views::Int=6, threads::Int=get_default_threads(),
        temperature::Union{Nothing,Float64}=nothing,
        embedding_size::Int=2048, coverage_embedding_size::Int=2048,
        batch_size::Int=1024, extra_args::Vector{String}=String[])

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
        "-b", string(batch_size),
    ]
    if temperature !== nothing
        push!(cmd_args, "-l")
        push!(cmd_args, string(temperature))
    end
    append!(cmd_args, extra_args)

    run(`$(CONDA_RUNNER) run --live-stream -n comebin $(cmd_args)`)
    bins_tsv = joinpath(outdir, "bins.tsv")
    return (; outdir, bins_tsv)
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
        completeness_threshold::Union{Nothing, Float64}=nothing,
        contamination_threshold::Union{Nothing, Float64}=nothing,
        ani_threshold::Float64=0.99,
        threads::Int=get_default_threads(),
        extra_args::Vector{String}=String[])

    isempty(genomes) && error("No genomes provided to dRep")
    for genome in genomes
        isfile(genome) || error("Genome file not found: $(genome)")
    end
    mkpath(outdir)

    add_bioconda_env("drep")

    cmd_args = String[
        "dRep", "dereplicate", outdir,
        "-g",
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

    winning_genomes = joinpath(outdir, "dereplicated_genomes.csv")
    return (; outdir, winning_genomes)
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
function parse_bin_assignments(file::String; contig_col::AbstractString="contig",
        bin_col::AbstractString="bin")::DataFrames.DataFrame

    isfile(file) || error("Assignment file not found: $(file)")
    df = DataFrames.DataFrame(CSV.File(file; delim='\t', ignorerepeated=true))
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
    df = DataFrames.DataFrame(CSV.File(file; delim=',', ignorerepeated=true))
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
"""
function run_taxometer(; contigs_fasta::String, depth_file::String, taxonomy_file::String,
        outdir::String, threads::Int=get_default_threads(),
        extra_args::Vector{String}=String[])

    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(depth_file) || error("Depth file not found: $(depth_file)")
    isfile(taxonomy_file) || error("Taxonomy file not found: $(taxonomy_file)")
    mkpath(outdir)

    env_name = _ensure_vamb_installed()
    cmd_args = String[
        "vamb", "taxometer",
        "--fasta", contigs_fasta,
        "--jgi", depth_file,
        "--taxonomy", taxonomy_file,
        "--outdir", outdir,
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
"""
function run_taxvamb(; contigs_fasta::String, depth_file::String, taxonomy_file::String,
        outdir::String, threads::Int=get_default_threads(),
        extra_args::Vector{String}=String[])

    isfile(contigs_fasta) || error("Contigs FASTA not found: $(contigs_fasta)")
    isfile(depth_file) || error("Depth file not found: $(depth_file)")
    isfile(taxonomy_file) || error("Taxonomy file not found: $(taxonomy_file)")
    mkpath(outdir)

    env_name = _ensure_vamb_installed()
    cmd_args = String[
        "vamb", "bin", "taxvamb",
        "--fasta", contigs_fasta,
        "--jgi", depth_file,
        "--taxonomy", taxonomy_file,
        "--outdir", outdir,
    ]
    if !any(arg -> arg == "-p" || startswith(arg, "--threads"), extra_args)
        append!(cmd_args, ["-p", string(threads)])
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) $(cmd_args)`)
    clusters_tsv = joinpath(outdir, "clusters.tsv")
    return (; outdir, clusters_tsv)
end

"""
    run_genomeface(; contigs_fasta, coverage_table, outdir, threads=get_default_threads(), extra_args=String[])

Run GenomeFace for binning with marker- and coverage-aware clustering.

Note: This wrapper is disabled until the `genomeface` executable can be located again.
"""
function run_genomeface(; contigs_fasta::String, coverage_table::String, outdir::String,
        threads::Int=get_default_threads(), extra_args::Vector{String}=String[])
    error("GenomeFace wrapper disabled: genomeface executable unavailable. Re-enable when the binary is found.")
end

"""
    run_magmax_merge(; bins_dirs, outdir, threads=get_default_threads(), extra_args=String[])

Merge MAGs/bins from multiple binners using MAGmax.
"""
function run_magmax_merge(; bins_dirs::Vector{String}, outdir::String,
        threads::Int=get_default_threads(), extra_args::Vector{String}=String[])

    isempty(bins_dirs) && error("No bins directories provided to MAGmax")
    for dir in bins_dirs
        isdir(dir) || error("Bins directory not found: $(dir)")
    end
    mkpath(outdir)

    add_bioconda_env("magmax")

    cmd_args = String["magmax", "-o", outdir, "-t", string(threads)]
    for dir in bins_dirs
        push!(cmd_args, "-i")
        push!(cmd_args, dir)
    end
    append!(cmd_args, extra_args)

    run(`$(CONDA_RUNNER) run --live-stream -n magmax $(cmd_args)`)
    return (; outdir)
end
