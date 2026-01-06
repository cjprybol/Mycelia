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
function _vamb_abundance_tsv(depth_file::String; output_dir::Union{Nothing,String}=nothing)
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
        isempty(sample_cols) && error("No sample depth columns found in JGI depth file: $(depth_file)")

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

function _find_first_matching_file(dir::String, patterns::Vector{Regex}; recursive::Bool=false)
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

function _find_first_matching_dir(dir::String, patterns::Vector{Regex}; recursive::Bool=false)
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
        minfasta::Int=2000, threads::Int=get_default_threads(),
        extra_args::Vector{String}=String[])

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
        "-m", string(minfasta),
    ]
    if !any(arg -> arg == "-p" || startswith(arg, "--threads"), extra_args)
        append!(cmd_args, ["-p", string(threads)])
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) $(cmd_args)`)

    clusters_tsv = _find_first_matching_file(outdir, [r"vae_clusters.*\.tsv$", r"clusters.*\.tsv$"])
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
        mapping_file::String, outdir::String, assembler::String="custom",
        threads::Int=get_default_threads(), extra_args::Vector{String}=String[])

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
        "--nthreads", string(threads),
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
    bins_dir = _find_first_matching_dir(outdir, [r"bins$"]; recursive=true)
    bins_tsv = _find_first_matching_file(outdir, [r"bins.*\.tsv$"]; recursive=true)
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

    winning_genomes = _find_first_matching_file(
        outdir,
        [r"Widb\.csv$", r"dereplicated_genomes\.csv$"];
        recursive=true
    )
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

# Arguments
- `contigs_fasta::String`: FASTA file with assembled contigs
- `depth_file::String`: JGI depth table or VAMB abundance TSV
- `taxonomy_file::String`: Taxonomy TSV
- `outdir::String`: Output directory
- `threads::Int`: Number of threads (default: all CPUs)
- `extra_args::Vector{String}`: Additional command-line arguments
"""
function run_taxometer(; contigs_fasta::String, depth_file::String, taxonomy_file::String,
        outdir::String, threads::Int=get_default_threads(),
        extra_args::Vector{String}=String[])

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

# Arguments
- `contigs_fasta::String`: FASTA file with assembled contigs
- `depth_file::String`: JGI depth table or VAMB abundance TSV
- `taxonomy_file::String`: Taxonomy TSV
- `outdir::String`: Output directory
- `threads::Int`: Number of threads (default: all CPUs)
- `extra_args::Vector{String}`: Additional command-line arguments
"""
function run_taxvamb(; contigs_fasta::String, depth_file::String, taxonomy_file::String,
        outdir::String, threads::Int=get_default_threads(),
        extra_args::Vector{String}=String[])

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
        "--outdir", outdir,
    ]
    if !any(arg -> arg == "-p" || startswith(arg, "--threads"), extra_args)
        append!(cmd_args, ["-p", string(threads)])
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n $(env_name) $(cmd_args)`)
    clusters_tsv = _find_first_matching_file(outdir, [r"vae_clusters.*\.tsv$", r"clusters.*\.tsv$"])
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

    bins_input_dir = joinpath(outdir, "magmax_bins_input")
    mkpath(bins_input_dir)
    for dir in bins_dirs
        for name in readdir(dir)
            src = joinpath(dir, name)
            isfile(src) || continue
            dest = joinpath(bins_input_dir, "$(basename(dir))__$(name)")
            cp(src, dest; force=true)
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

    cmd = Cmd(`$(CONDA_RUNNER) run --live-stream -n magmax $(cmd_args)`; dir=outdir)
    run(cmd)
    return (; outdir, bins_input_dir)
end
