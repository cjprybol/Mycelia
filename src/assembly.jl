"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run MEGAHIT assembler for metagenomic short read assembly and downstream graph export.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::Union{Nothing,String}`: Output directory path (default: inferred from FASTQ prefix + "_megahit")
- `min_contig_len::Int`: Minimum contig length (default: 200)
- `k_list::String`: k-mer sizes to use (default: "21,29,39,59,79,99,119,141")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to final contigs file
- `fastg::String`: Path to MEGAHIT FASTG export
- `gfa::String`: Path to GFA converted with gfatools

# Details
- Automatically creates and uses a conda environment with megahit
- Exports both FASTG and GFA via gfatools
- Optimized for metagenomic assemblies with varying coverage
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_megahit(;fastq1, fastq2=nothing, outdir=nothing, min_contig_len=200, k_list="21,29,39,59,79,99,119,141", threads=get_default_threads())
    Mycelia.add_bioconda_env("megahit")
    # Default output directory derived from FASTQ prefix
    if isnothing(outdir)
        cleaned1 = replace(fastq1, Mycelia.FASTQ_REGEX => "")
        if isnothing(fastq2)
            prefix = cleaned1
        else
            cleaned2 = replace(fastq2, Mycelia.FASTQ_REGEX => "")
            prefix = Mycelia.find_matching_prefix(cleaned1, cleaned2)
            if isempty(prefix)
                prefix = cleaned1
            end
        end
        outdir = prefix * "_megahit"
    end
    # Infer the final k-mer size from contig identifiers produced by MEGAHIT
    infer_final_k(contigs_path) = begin
        ks = Int[]
        open(contigs_path) do io
            for record in FASTX.FASTA.Reader(io)
                id = FASTX.identifier(record)
                m = match(r"^k(\d+)", id)
                m === nothing && continue
                push!(ks, parse(Int, m.captures[1]))
            end
        end
        @assert !isempty(ks) "Could not infer k-mer size from contig identifiers in $(contigs_path)"
        unique_ks = unique(ks)
        @assert length(unique_ks) == 1 "Expected a single final k-mer size, found $(unique_ks)"
        return first(unique_ks)
    end
    
    # MEGAHIT requires the output directory to not exist, so check output file first
    contigs_path = joinpath(outdir, "final.contigs.fa")
    if !isfile(contigs_path)
        # Remove output directory if it exists (MEGAHIT will create it)
        if isdir(outdir)
            rm(outdir, recursive=true)
        end
        
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n megahit megahit -r $(fastq1) -o $(outdir) --min-contig-len $(min_contig_len) --k-list $(k_list) -t $(threads)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n megahit megahit -1 $(fastq1) -2 $(fastq2) -o $(outdir) --min-contig-len $(min_contig_len) --k-list $(k_list) -t $(threads)`)
        end
    else
        # If output already exists, ensure directory exists for return value
        mkpath(outdir)
    end

    # Derive FASTG and GFA outputs
    fastg_path = replace(contigs_path, ".fa" => ".fastg")
    if !isfile(fastg_path)
        final_k = infer_final_k(contigs_path)
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n megahit megahit_core contig2fastg $(final_k) $(contigs_path)`, fastg_path))
        @assert isfile(fastg_path)
        @assert filesize(fastg_path) > 0
    end

    gfa_path = fastg_path * ".gfa"
    if !isfile(gfa_path)
        # if isfile(fastg_path) && (filesize(fastg_path) > 0)
        Mycelia.add_bioconda_env("gfatools")
        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n gfatools gfatools view $(fastg_path)`, gfa_path))
        # else
        #     @warn "final.contigs.fastg not found, skipping gfatools step"
        #     gfa_path = missing
        # end
    end

    return (;outdir, contigs=contigs_path, fastg=fastg_path, gfa=gfa_path)
end

"""
    run_assembler(label, run_fn)

Run an assembler helper with timing and error reporting.
"""
function run_assembler(label, run_fn)
    result = nothing
    runtime = missing
    try
        runtime = @elapsed begin
            result = run_fn()
        end
        println("  $(label) completed in $(round(runtime, digits=2))s")
    catch e
        @warn "$(label) failed" exception=e
    end
    return result, runtime
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Build a bash command string to run MEGAHIT and downstream graph export.

This is intended for use with SLURM `sbatch --wrap` (e.g. via `scg_sbatch`).

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::Union{Nothing,String}`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::Union{Nothing,String}`: Output directory path (default: inferred from FASTQ prefix + "_megahit")
- `min_contig_len::Int`: Minimum contig length (default: 200)
- `k_list::String`: k-mer sizes to use (default: "21,29,39,59,79,99,119,141")
- `threads::Int`: Number of CPU threads

# Returns
Named tuple containing:
- `cmd::String`: Bash snippet that assembles with MEGAHIT and exports FASTG + GFA
- `outdir::String`: Path to output directory
- `contigs::String`: Expected path to final contigs file
- `fastg::String`: Expected path to MEGAHIT FASTG export
- `gfa::String`: Expected path to GFA converted with gfatools

# Notes
- The generated bash is idempotent:
  - Skips assembly if `final.contigs.fa` already exists.
  - Skips FASTG/GFA export if their files already exist.
- The MEGAHIT environment is ensured on the Julia side via `Mycelia.add_bioconda_env`.
"""
function megahit_cmd(;
    fastq1,
    fastq2::Union{Nothing,String}=nothing,
    outdir::Union{Nothing,String}=nothing,
    min_contig_len::Int=200,
    k_list::String="21,29,39,59,79,99,119,141",
    threads::Int=get_default_threads(),
)
    # Ensure environments exist (shared filesystem, so compute nodes will see them)
    Mycelia.add_bioconda_env("megahit")
    Mycelia.add_bioconda_env("gfatools")

    # Derive default outdir as in run_megahit
    if isnothing(outdir)
        cleaned1 = replace(fastq1, Mycelia.FASTQ_REGEX => "")
        if isnothing(fastq2)
            prefix = cleaned1
        else
            cleaned2 = replace(fastq2, Mycelia.FASTQ_REGEX => "")
            prefix = Mycelia.find_matching_prefix(cleaned1, cleaned2)
            if isempty(prefix)
                prefix = cleaned1
            end
        end
        outdir = prefix * "_megahit"
    end

    contigs_path = joinpath(outdir, "final.contigs.fa")
    fastg_path   = replace(contigs_path, ".fa" => ".fastg")
    gfa_path     = fastg_path * ".gfa"

    lines = String[]

    # Safe shell behaviour
    push!(lines, "set -euo pipefail")
    push!(lines, "")
    push!(lines, "# MEGAHIT assembly")
    push!(lines, "if [ ! -f \"$(contigs_path)\" ]; then")
    push!(lines, "  rm -rf \"$(outdir)\" || true")
    if isnothing(fastq2)
        push!(lines,
            "  $(Mycelia.CONDA_RUNNER) run --live-stream -n megahit " *
            "megahit -r \"$(fastq1)\" -o \"$(outdir)\" " *
            "--min-contig-len $(min_contig_len) --k-list $(k_list) -t $(threads)"
        )
    else
        push!(lines,
            "  $(Mycelia.CONDA_RUNNER) run --live-stream -n megahit " *
            "megahit -1 \"$(fastq1)\" -2 \"$(fastq2)\" -o \"$(outdir)\" " *
            "--min-contig-len $(min_contig_len) --k-list $(k_list) -t $(threads)"
        )
    end
    push!(lines, "fi")
    push!(lines, "")

    # Infer final k on the node (AWK version of infer_final_k)
    push!(lines, "# Convert contigs to FASTG")
    push!(lines, "if [ ! -f \"$(fastg_path)\" ]; then")
    push!(lines, "  final_k=\$(awk '")
    push!(lines, "    /^>k[0-9]+/ {")
    push!(lines, "      header = substr(\$0, 2)      # drop leading >")
    push!(lines, "      sub(/^k/, \"\", header)      # drop leading k")
    push!(lines, "      split(header, a, \"_\")      # e.g. 141_1 -> a[1]=141")
    push!(lines, "      ks[a[1]] = 1")
    push!(lines, "    }")
    push!(lines, "    END {")
    push!(lines, "      n = 0")
    push!(lines, "      for (k in ks) { last = k; n++ }")
    push!(lines,
        "      if (n == 0) { print \"ERROR: Could not infer k-mer size from contig identifiers in $(contigs_path)\" > \"/dev/stderr\"; exit 1 }"
    )
    push!(lines,
        "      if (n > 1)  { print \"ERROR: Expected a single final k-mer size, found multiple values\" > \"/dev/stderr\"; exit 1 }"
    )
    push!(lines, "      print last")
    push!(lines, "    }")
    push!(lines, "  ' \"$(contigs_path)\")")
    push!(lines,
        "  $(Mycelia.CONDA_RUNNER) run --live-stream -n megahit " *
        "megahit_core contig2fastg \"\$final_k\" \"$(contigs_path)\" > \"$(fastg_path)\""
    )
    push!(lines, "fi")
    push!(lines, "")

    push!(lines, "# Convert FASTG to GFA")
    push!(lines, "if [ ! -f \"$(gfa_path)\" ]; then")
    push!(lines,
        "  $(Mycelia.CONDA_RUNNER) run --live-stream -n gfatools " *
        "gfatools view \"$(fastg_path)\" > \"$(gfa_path)\""
    )
    push!(lines, "fi")

    cmd = join(lines, "\n")
    return (; cmd, outdir, contigs=contigs_path, fastg=fastg_path, gfa=gfa_path)
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaSPAdes assembler for metagenomic short read assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "metaspades_output")
- `k_list::String`: k-mer sizes to use (default: "21,33,55,77")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `scaffolds::String`: Path to scaffolds file
- `graph::String`: Path to assembly graph in GFA format

# Details
- Automatically creates and uses a conda environment with spades
- Designed for metagenomic datasets with uneven coverage
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_metaspades(;fastq1, fastq2=nothing, outdir="metaspades_output", k_list="21,33,55,77", threads = get_default_threads())
    Mycelia.add_bioconda_env("spades")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "contigs.fasta"))
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades metaspades.py -s $(fastq1) -o $(outdir) -k $(k_list) -t $(threads)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades metaspades.py -1 $(fastq1) -2 $(fastq2) -o $(outdir) -k $(k_list) -t $(threads)`)
        end
    end
    return (;outdir, contigs=joinpath(outdir, "contigs.fasta"), scaffolds=joinpath(outdir, "scaffolds.fasta"), graph=joinpath(outdir, "assembly_graph_with_scaffolds.gfa"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run SPAdes assembler for single genome isolate assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "spades_output")
- `k_list::String`: k-mer sizes to use (default: "21,33,55,77")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to assembled contigs file
- `scaffolds::String`: Path to scaffolds file
- `graph::String`: Path to assembly graph in GFA format

# Details
- Automatically creates and uses a conda environment with spades
- Designed for single bacterial/archaeal genome assembly
- Optimized for uniform coverage isolate data
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_spades(;fastq1, fastq2=nothing, outdir="spades_output", k_list="21,33,55,77", threads = get_default_threads())
    Mycelia.add_bioconda_env("spades")
    mkpath(outdir)
    
    if !isfile(joinpath(outdir, "contigs.fasta"))
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades spades.py -s $(fastq1) -o $(outdir) -k $(k_list) -t $(threads)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades spades.py -1 $(fastq1) -2 $(fastq2) -o $(outdir) -k $(k_list) -t $(threads)`)
        end
    end
    return (;outdir, contigs=joinpath(outdir, "contigs.fasta"), scaffolds=joinpath(outdir, "scaffolds.fasta"), graph=joinpath(outdir, "assembly_graph_with_scaffolds.gfa"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run SKESA assembler for high-accuracy bacterial genome assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "skesa_output")
- `min_contig_len::Int`: Minimum contig length (default: 200)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file

# Details
- Uses conservative assembly approach for high accuracy
- Optimized for bacterial genomes with uniform coverage
- Automatically creates and uses a conda environment with skesa
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_skesa(;fastq1, fastq2=nothing, outdir="skesa_output", min_contig_len=200, threads = get_default_threads())
    Mycelia.add_bioconda_env("skesa")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "contigs.fa")
    if !isfile(contigs_file)
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n skesa skesa --reads $(fastq1) --contigs_out $(contigs_file) --cores $(threads) --min_contig $(min_contig_len)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n skesa skesa --reads $(fastq1),$(fastq2) --contigs_out $(contigs_file) --cores $(threads) --min_contig $(min_contig_len)`)
        end
    end
    return (;outdir, contigs=contigs_file)
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Run IDBA-UD assembler for metagenomic assembly with uneven depth.

# # Arguments
# - `fastq1::String`: Path to first paired-end FASTQ file
# - `fastq2::String`: Path to second paired-end FASTQ file
# - `outdir::String`: Output directory path (default: "idba_ud_output")
# - `min_k::Int`: Minimum k-mer size (default: 20)
# - `max_k::Int`: Maximum k-mer size (default: 100)
# - `step::Int`: K-mer size increment step (default: 20)

# # Returns
# Named tuple containing:
# - `outdir::String`: Path to output directory
# - `contigs::String`: Path to contigs file
# - `scaffolds::String`: Path to scaffolds file

# # Details
# - Uses multi-k-mer iterative assembly approach
# - Specifically designed for metagenomic data with uneven coverage depths
# - Automatically creates and uses a conda environment with idba
# - Requires paired-end reads merged to fasta format
# - Skips assembly if output directory already exists
# - Utilizes all available CPU threads
# """
# function run_idba_ud(;fastq1, fastq2, outdir="idba_ud_output", min_k=20, max_k=100, step=20)
#     Mycelia.add_bioconda_env("idba")
#     mkpath(outdir)
    
#     # IDBA-UD requires paired reads in fasta format
#     merged_fasta = joinpath(outdir, "merged_reads.fa")
#     contig_file = joinpath(outdir, "contig.fa")
#     scaffold_file = joinpath(outdir, "scaffold.fa")
    
#     if !isfile(contig_file)
#         # Convert and merge FASTQ to FASTA format for IDBA-UD
#         if !isfile(merged_fasta)
#             run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n idba fq2fa --merge $(fastq1) $(fastq2) $(merged_fasta)`)
#         end
        
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n idba idba_ud -r $(merged_fasta) -o $(outdir) --mink $(min_k) --maxk $(max_k) --step $(step) --num_threads $(Sys.CPU_THREADS)`)
#     end
#     return (;outdir, contigs=contig_file, scaffolds=scaffold_file)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Flye assembler for long read assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "flye_output")
- `genome_size::String`: Estimated genome size (e.g., "5m", "1.2g") (default: nothing, auto-estimated)
- `read_type::String`: Type of reads ("pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq")
- `min_overlap::Int`: Minimum overlap between reads (default: nothing, auto-selected)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly FASTA file
- `graph::String`: Path to assembly graph in GFA format

# Details
- Automatically creates and uses a conda environment with flye
- Supports various long read technologies and quality levels
- Skips assembly if output directory already exists
- Thread count is determined by get_default_threads() and capped at 128 (Flye's maximum)
"""
function run_flye(;fastq, outdir="flye_output", genome_size=nothing, read_type="pacbio-hifi", min_overlap=nothing, threads = min(get_default_threads(), FLYE_MAX_THREADS))
    Mycelia.add_bioconda_env("flye")
    mkpath(outdir)

    if !isfile(joinpath(outdir, "assembly.fasta"))
        cmd_args = ["flye", "--$(read_type)", fastq, "--out-dir", outdir, "--threads", string(threads)]
        if !isnothing(genome_size)
            push!(cmd_args, "--genome-size", string(genome_size))
        end
        if !isnothing(min_overlap)
            push!(cmd_args, "--min-overlap", string(min_overlap))
        end
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n flye $(cmd_args)`)
    end
    return (;outdir, assembly=joinpath(outdir, "assembly.fasta"), graph=joinpath(outdir, "assembly_graph.gfa"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaFlye assembler for long-read metagenomic assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "metaflye_output")
- `genome_size::String`: Estimated genome size (e.g., "5m", "1.2g") (default: nothing, auto-estimated)
- `read_type::String`: Type of reads ("pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq")
- `meta::Bool`: Enable metagenome mode (default: true)
- `min_overlap::Int`: Minimum overlap between reads (default: auto-selected)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly FASTA file
- `graph::String`: Path to assembly graph in GFA format

# Details
- Uses metaFlye's repeat graph approach optimized for metagenomic data
- Implements solid k-mer selection combining global and local k-mer distributions
- Handles uneven coverage and strain variation in metagenomic samples
- Automatically creates and uses a conda environment with flye
- Skips assembly if output directory already exists
- Thread count is determined by get_default_threads() and capped at 128 (Flye's maximum)
"""
function run_metaflye(;fastq, outdir="metaflye_output", genome_size=nothing, read_type="pacbio-hifi", meta=true, min_overlap=nothing, threads = min(get_default_threads(), FLYE_MAX_THREADS))
    Mycelia.add_bioconda_env("flye")
    mkpath(outdir)

    if !isfile(joinpath(outdir, "assembly.fasta"))
        cmd_args = ["flye", "--$(read_type)", fastq, "--out-dir", outdir, "--threads", string(threads)]
        if !isnothing(genome_size)
            push!(cmd_args, "--genome-size", string(genome_size))
        end

        if meta
            push!(cmd_args, "--meta")
        end

        if !isnothing(min_overlap)
            push!(cmd_args, "--min-overlap", string(min_overlap))
        end

        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n flye $(cmd_args)`)
    end
    return (;outdir, assembly=joinpath(outdir, "assembly.fasta"), graph=joinpath(outdir, "assembly_graph.gfa"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Canu assembler for long read assembly.

# Arguments
- `fastq::String`: Path to input FASTQ file containing long reads
- `outdir::String`: Output directory path (default: "canu_output")
- `genome_size::String`: Estimated genome size (e.g., "5m", "1.2g") (default: nothing, auto-estimated)
- `read_type::String`: Type of reads ("pacbio", "nanopore")
- `threads::Integer`: Maximum threads to use for the run (default: get_default_threads())
- `cor_threads::Integer`: Threads to allocate to the correction stage (default: matches `threads`)
- `stopOnLowCoverage::Integer`: Minimum coverage required to continue assembly (default: 10)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly FASTA file
- `graph::String`: Path to assembly graph in GFA format

# Details
- Automatically creates and uses a conda environment with canu
- Includes error correction, trimming, and assembly stages
- Skips assembly if output directory already exists
- Thread count is determined by get_default_threads()
- Can reduce `stopOnLowCoverage` for CI environments or low-coverage datasets
- Sets `corThreads` to the requested thread count to avoid Canu's default (4) exceeding `maxThreads` on small runners
- Uses `saveReads=false` to skip saving intermediate corrected/trimmed reads (avoids gzip issues on some filesystems)
- Uses `useGrid=false` to disable automatic SLURM/grid job submission
"""
function run_canu(;fastq, outdir="canu_output", genome_size, read_type="pacbio", stopOnLowCoverage=10, threads = get_default_threads(), cor_threads = threads)
    Mycelia.add_bioconda_env("canu")
    mkpath(outdir)
    threads = max(threads, 1)
    cor_threads = clamp(cor_threads, 1, threads)

    prefix = splitext(basename(fastq))[1]
    if !isfile(joinpath(outdir, "$(prefix).contigs.fasta"))
        if read_type == "pacbio"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -pacbio $(fastq) maxThreads=$(threads) corThreads=$(cor_threads) stopOnLowCoverage=$(stopOnLowCoverage) saveReads=false useGrid=false`)
        elseif read_type == "nanopore"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -nanopore $(fastq) maxThreads=$(threads) corThreads=$(cor_threads) stopOnLowCoverage=$(stopOnLowCoverage) saveReads=false useGrid=false`)
        else
            error("Unsupported read type: $(read_type). Use 'pacbio' or 'nanopore'.")
        end
    end
    return (;outdir, assembly=joinpath(outdir, "$(prefix).contigs.fasta"), graph=joinpath(outdir, "$(prefix).contigs.gfa"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run the hifiasm genome assembler on PacBio HiFi reads.

# Arguments
- `fastq::String`: Path to input FASTQ file containing HiFi reads
- `outdir::String`: Output directory path (default: "\${basename(fastq)}_hifiasm")
- `bloom_filter::Int`: Bloom filter flag (default: -1 for automatic, 0 to disable 16GB filter)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `hifiasm_outprefix::String`: Prefix used for hifiasm output files

# Details
- Automatically creates and uses a conda environment with hifiasm
- Uses primary assembly mode (--primary) optimized for inbred samples
- Skips assembly if output files already exist at the specified prefix
- Utilizes all available CPU threads
- Bloom filter can be disabled (-f0) for small genomes to reduce memory usage from 16GB
"""
function run_hifiasm(;fastq, outdir=basename(fastq) * "_hifiasm", bloom_filter=-1, threads = get_default_threads())
    Mycelia.add_bioconda_env("hifiasm")
    hifiasm_outprefix = joinpath(outdir, basename(fastq) * ".hifiasm")
    # Check if output directory exists before trying to read it
    hifiasm_outputs = isdir(outdir) ? filter(x -> occursin(hifiasm_outprefix, x), readdir(outdir, join=true)) : String[]
    # https://hifiasm.readthedocs.io/en/latest/faq.html#are-inbred-homozygous-genomes-supported
    if isempty(hifiasm_outputs)
        mkpath(outdir)
        # Build command with optional bloom filter flag
        cmd_args = ["hifiasm", "--primary", "-l0", "-o", hifiasm_outprefix, "-t", string(threads)]
        if bloom_filter >= 0
            push!(cmd_args, "-f$(bloom_filter)")
        end
        push!(cmd_args, fastq)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hifiasm $(cmd_args)`)
    end
    return (;outdir, hifiasm_outprefix)
end

"""
    hifiasm_primary_contigs(hifiasm_result) -> Union{String,Nothing}

Return the primary contigs FASTA path produced by `run_hifiasm`, if present.
"""
function hifiasm_primary_contigs(hifiasm_result)
    if hifiasm_result === nothing
        return nothing
    end

    candidates = [
        hifiasm_result.hifiasm_outprefix * ".p_ctg.fa",
        hifiasm_result.hifiasm_outprefix * ".p_ctg.fasta"
    ]
    for candidate in candidates
        if isfile(candidate)
            return candidate
        end
    end
    return nothing
end

# disabled due to poor performance relative to other long read metagenomic assemblers in tests
# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Run hifiasm-meta assembler for metagenomic PacBio HiFi assembly.

# # Arguments
# - `fastq::String`: Path to input FASTQ file containing HiFi reads
# - `outdir::String`: Output directory path (default: "\${basename(fastq)}_hifiasm_meta")
# - `bloom_filter::Int`: Bloom filter flag (default: -1 for automatic, 0 to disable 16GB filter)
# - `read_selection::Bool`: Enable read selection for mock/small datasets (default: false)

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run hifiasm-meta for metagenomic long-read assembly.

# Arguments
- `fastq::String`: Input long-read FASTQ
- `outdir::String`: Output directory path (default: "`basename(fastq)_hifiasm_meta`")
- `bloom_filter::Int`: Bloom filter size; use -1 to leave default, 0 to disable for small genomes
- `read_selection::Bool`: Enable `-S` read selection for low-complexity/mock datasets

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `hifiasm_outprefix::String`: Prefix used for hifiasm-meta output files

# Details
- Uses hifiasm-meta's string graph approach for metagenomic strain resolution
- Automatically creates and uses a conda environment with `hifiasm_meta`
- Skips assembly if output files already exist at the specified prefix
- Utilizes all available CPU threads
- Bloom filter can be disabled (`-f0`) for small genomes to reduce memory usage from ~16GB
- Read selection (`-S`) can be enabled for mock/small datasets to handle low complexity data
"""
function run_hifiasm_meta(;fastq, outdir=basename(fastq) * "_hifiasm_meta", bloom_filter=-1, read_selection=false, threads = get_default_threads())
    Mycelia.add_bioconda_env("hifiasm_meta")
    hifiasm_outprefix = joinpath(outdir, basename(fastq) * ".hifiasm_meta")
    # Check if output directory exists before trying to read it
    hifiasm_outputs = isdir(outdir) ? filter(x -> occursin(hifiasm_outprefix, x), readdir(outdir, join=true)) : String[]

    if isempty(hifiasm_outputs)
        mkpath(outdir)
        # Build command with optional flags
        cmd_args = ["hifiasm_meta", "-t", string(threads), "-o", hifiasm_outprefix]
        if bloom_filter >= 0
            push!(cmd_args, "-f$(bloom_filter)")
        end
        if read_selection
            push!(cmd_args, "-S")
        end
        push!(cmd_args, fastq)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hifiasm_meta $(cmd_args)`)
    end
    return (;outdir, hifiasm_outprefix)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run hybrid assembly combining short and long reads using Unicycler.

# Arguments
- `short_1::String`: Path to first short read FASTQ file
- `short_2::String`: Path to second short read FASTQ file (optional)
- `long_reads::String`: Path to long read FASTQ file
- `outdir::String`: Output directory path (default: "unicycler_output")
- `threads::Int`: Number of threads to use (default: `get_default_threads()`)
- `spades_options::Union{Nothing,String}`: Extra SPAdes options (default: `nothing`)
- `kmers::Union{Nothing,String}`: Explicit SPAdes k-mers (e.g., "21,33,55")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `assembly::String`: Path to final assembly file
- `graph::String`: Path to assembly graph in GFA format

# Details
- Automatically creates and uses a conda environment with unicycler
- Combines short read accuracy with long read scaffolding
- Skips assembly if output directory already exists
- Utilizes all available CPU threads
"""
function run_unicycler(;short_1, short_2=nothing, long_reads, outdir="unicycler_output", threads = get_default_threads(), spades_options::Union{Nothing,String}=nothing, kmers::Union{Nothing,String}=nothing)
    Mycelia.add_bioconda_env("unicycler")
    
    # Unicycler requires the output directory to not exist, so check output file first
    if !isfile(joinpath(outdir, "assembly.fasta"))
        # Remove output directory if it exists (Unicycler will create it)
        if isdir(outdir)
            rm(outdir, recursive=true)
        end
        
        spades_args = isnothing(spades_options) ? String[] : ["--spades_options=$(spades_options)"]
        kmer_args = isnothing(kmers) ? String[] : ["--kmers=$(kmers)"]
        if isnothing(short_2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n unicycler unicycler -s $(short_1) -l $(long_reads) -o $(outdir) -t $(threads) $(kmer_args...) $(spades_args...)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n unicycler unicycler -1 $(short_1) -2 $(short_2) -l $(long_reads) -o $(outdir) -t $(threads) $(kmer_args...) $(spades_args...)`)
        end
    else
        # If output already exists, ensure directory exists for return value
        mkpath(outdir)
    end
    return (;outdir, assembly=joinpath(outdir, "assembly.fasta"), graph=joinpath(outdir, "assembly.gfa"))
end

# ============================================================================
# Hybrid/Plasmid Assembly and Contig Reorientation
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install the Hybracter databases in `db_dir` if missing.

# Arguments
- `db_dir::AbstractString`: Directory to store the databases.

# Keywords
- `force::Bool=false`: If true, re-download and overwrite the database directory.

# Returns
- `String`: Path to the database directory.
"""
function install_hybracter_db(db_dir::AbstractString; force::Bool=false)
    db_dir = abspath(String(db_dir))
    if isdir(db_dir) && !force && !isempty(readdir(db_dir))
        @info "Hybracter database already present; skipping download." db_dir
        return db_dir
    end
    if force && isdir(db_dir)
        rm(db_dir; recursive=true, force=true)
    end
    mkpath(db_dir)
    Mycelia.add_bioconda_env("hybracter")
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hybracter hybracter install --databases $(db_dir)`)
    return db_dir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Hybracter hybrid assembly for a single isolate with long and short reads.

# Arguments
- `long_reads::AbstractString`: Long read FASTQ file.
- `read1::AbstractString`: Short read R1 FASTQ file.
- `read2::AbstractString`: Short read R2 FASTQ file.
- `chrom_size::Int`: Estimated chromosome size (bp).

# Keywords
- `sample_name::AbstractString`: Sample identifier used in Hybracter outputs.
- `outdir::AbstractString`: Output directory for the run.
- `db_dir::Union{Nothing,AbstractString}=nothing`: Database directory (downloaded if missing).
- `threads::Int=get_default_threads()`: Thread count.
- `extra_args::Vector{String}=String[]`: Additional Hybracter CLI arguments.
- `force::Bool=false`: If true, re-run and overwrite existing outputs.

# Returns
Named tuple containing:
- `outdir::String`: Output directory.
- `sample_dir::String`: Sample subdirectory within the output directory.
- `summary_tsv::String`: Summary TSV path.
- `final_fasta::String`: Final assembly FASTA path.
"""
function run_hybracter_hybrid_single(long_reads::AbstractString,
                                     read1::AbstractString,
                                     read2::AbstractString,
                                     chrom_size::Int;
                                     sample_name::AbstractString,
                                     outdir::AbstractString,
                                     db_dir::Union{Nothing,AbstractString}=nothing,
                                     threads::Int=get_default_threads(),
                                     extra_args::Vector{String}=String[],
                                     force::Bool=false)
    isfile(long_reads) || error("Long read file not found: $(long_reads)")
    isfile(read1) || error("Read 1 file not found: $(read1)")
    isfile(read2) || error("Read 2 file not found: $(read2)")
    isempty(sample_name) && error("sample_name cannot be empty")

    outdir = abspath(String(outdir))
    sample_dir = joinpath(outdir, sample_name)
    complete_final = joinpath(outdir, "FINAL_OUTPUT", "complete", sample_name * "_final.fasta")
    incomplete_final = joinpath(outdir, "FINAL_OUTPUT", "incomplete", sample_name * "_final.fasta")
    legacy_final = joinpath(sample_dir, "final_assembly.fasta")
    complete_summary = joinpath(outdir, "FINAL_OUTPUT", "complete", sample_name * "_summary.tsv")
    incomplete_summary = joinpath(outdir, "FINAL_OUTPUT", "incomplete", sample_name * "_summary.tsv")
    legacy_summary = joinpath(sample_dir, "summary.tsv")

    if !force && (isfile(complete_final) || isfile(incomplete_final) || isfile(legacy_final))
        final_fasta = isfile(complete_final) ? complete_final : (isfile(incomplete_final) ? incomplete_final : legacy_final)
        summary_tsv = isfile(complete_summary) ? complete_summary : (isfile(incomplete_summary) ? incomplete_summary : legacy_summary)
        return (;outdir, sample_dir, summary_tsv, final_fasta)
    end

    if force && isdir(sample_dir)
        rm(sample_dir; recursive=true, force=true)
    end

    if db_dir === nothing || isempty(db_dir)
        db_dir = joinpath(homedir(), "workspace", ".hybracter")
    end
    if !isdir(db_dir) || isempty(readdir(db_dir))
        install_hybracter_db(db_dir; force=force)
    end

    Mycelia.add_bioconda_env("hybracter")

    cmd_args = String[
        "hybrid-single",
        "--longreads", long_reads,
        "--short_one", read1,
        "--short_two", read2,
        "--chromosome", string(chrom_size),
        "--sample", sample_name,
        "--output", outdir,
        "--databases", db_dir
    ]
    if !any(arg -> arg == "-t" || arg == "--threads" || startswith(arg, "--threads="), extra_args)
        push!(cmd_args, "--threads", string(max(threads, 1)))
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hybracter hybracter $(cmd_args)`)
    final_fasta = isfile(complete_final) ? complete_final : (isfile(incomplete_final) ? incomplete_final : legacy_final)
    summary_tsv = isfile(complete_summary) ? complete_summary : (isfile(incomplete_summary) ? incomplete_summary : legacy_summary)
    return (;outdir, sample_dir, summary_tsv, final_fasta)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Hybracter long-read-only assembly for a single isolate.

# Arguments
- `long_reads::AbstractString`: Long read FASTQ file.
- `chrom_size::Int`: Estimated chromosome size (bp).

# Keywords
- `sample_name::AbstractString`: Sample identifier used in Hybracter outputs.
- `outdir::AbstractString`: Output directory for the run.
- `db_dir::Union{Nothing,AbstractString}=nothing`: Database directory (downloaded if missing).
- `threads::Int=get_default_threads()`: Thread count.
- `extra_args::Vector{String}=String[]`: Additional Hybracter CLI arguments.
- `force::Bool=false`: If true, re-run and overwrite existing outputs.

# Returns
Named tuple containing:
- `outdir::String`: Output directory.
- `sample_dir::String`: Sample subdirectory within the output directory.
- `summary_tsv::String`: Summary TSV path.
- `final_fasta::String`: Final assembly FASTA path.
"""
function run_hybracter_long_single(long_reads::AbstractString,
                                   chrom_size::Int;
                                   sample_name::AbstractString,
                                   outdir::AbstractString,
                                   db_dir::Union{Nothing,AbstractString}=nothing,
                                   threads::Int=get_default_threads(),
                                   extra_args::Vector{String}=String[],
                                   force::Bool=false)
    isfile(long_reads) || error("Long read file not found: $(long_reads)")
    isempty(sample_name) && error("sample_name cannot be empty")

    outdir = abspath(String(outdir))
    sample_dir = joinpath(outdir, sample_name)
    complete_final = joinpath(outdir, "FINAL_OUTPUT", "complete", sample_name * "_final.fasta")
    incomplete_final = joinpath(outdir, "FINAL_OUTPUT", "incomplete", sample_name * "_final.fasta")
    legacy_final = joinpath(sample_dir, "final_assembly.fasta")
    complete_summary = joinpath(outdir, "FINAL_OUTPUT", "complete", sample_name * "_summary.tsv")
    incomplete_summary = joinpath(outdir, "FINAL_OUTPUT", "incomplete", sample_name * "_summary.tsv")
    legacy_summary = joinpath(sample_dir, "summary.tsv")

    if !force && (isfile(complete_final) || isfile(incomplete_final) || isfile(legacy_final))
        final_fasta = isfile(complete_final) ? complete_final : (isfile(incomplete_final) ? incomplete_final : legacy_final)
        summary_tsv = isfile(complete_summary) ? complete_summary : (isfile(incomplete_summary) ? incomplete_summary : legacy_summary)
        return (;outdir, sample_dir, summary_tsv, final_fasta)
    end

    if force && isdir(sample_dir)
        rm(sample_dir; recursive=true, force=true)
    end

    if db_dir === nothing || isempty(db_dir)
        db_dir = joinpath(homedir(), "workspace", ".hybracter")
    end
    if !isdir(db_dir) || isempty(readdir(db_dir))
        install_hybracter_db(db_dir; force=force)
    end

    Mycelia.add_bioconda_env("hybracter")

    cmd_args = String[
        "long-single",
        "--longreads", long_reads,
        "--chromosome", string(chrom_size),
        "--sample", sample_name,
        "--output", outdir,
        "--databases", db_dir
    ]
    if !any(arg -> arg == "-t" || arg == "--threads" || startswith(arg, "--threads="), extra_args)
        push!(cmd_args, "--threads", string(max(threads, 1)))
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hybracter hybracter $(cmd_args)`)
    final_fasta = isfile(complete_final) ? complete_final : (isfile(incomplete_final) ? incomplete_final : legacy_final)
    summary_tsv = isfile(complete_summary) ? complete_summary : (isfile(incomplete_summary) ? incomplete_summary : legacy_summary)
    return (;outdir, sample_dir, summary_tsv, final_fasta)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run DNAapler to reorient contigs in FASTA or GFA format.

# Arguments
- `input_file::AbstractString`: Contigs FASTA/GFA file.

# Keywords
- `outdir::AbstractString`: Output directory.
- `prefix::AbstractString="dnaapler"`: Output prefix.
- `threads::Int=get_default_threads()`: Thread count.
- `extra_args::Vector{String}=String[]`: Additional DNAapler CLI arguments.
- `force::Bool=false`: If true, re-run and overwrite existing outputs.

# Returns
Named tuple containing:
- `outdir::String`: Output directory.
- `reoriented::String`: Expected reoriented contig file path.
"""
function run_dnaapler_all(input_file::AbstractString;
                          outdir::AbstractString,
                          prefix::AbstractString="dnaapler",
                          threads::Int=get_default_threads(),
                          extra_args::Vector{String}=String[],
                          force::Bool=false)
    isfile(input_file) || error("Input file not found: $(input_file)")

    outdir = abspath(String(outdir))
    mkpath(outdir)

    lower_name = lowercase(input_file)
    ext = (endswith(lower_name, ".gfa") || endswith(lower_name, ".gfa.gz")) ? ".gfa" : ".fasta"
    reoriented = joinpath(outdir, "$(prefix)_reoriented$(ext)")

    if !force && isfile(reoriented)
        return (;outdir, reoriented)
    end

    Mycelia.add_bioconda_env("dnaapler")

    cmd_args = String[
        "all",
        "--input", input_file,
        "--output", outdir,
        "--prefix", prefix
    ]
    if !any(arg -> arg == "-t" || arg == "--threads" || startswith(arg, "--threads="), extra_args)
        push!(cmd_args, "--threads", string(max(threads, 1)))
    end
    if force
        push!(cmd_args, "--force")
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n dnaapler dnaapler $(cmd_args)`)
    return (;outdir, reoriented)
end

function _plassembler_db_ready(db_dir::AbstractString)
    return isdir(db_dir) && any(path -> endswith(path, ".msh"), readdir(db_dir; join=true))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Install the Plassembler database in `db_dir` if missing.

# Arguments
- `db_dir::AbstractString`: Directory to store the database.

# Keywords
- `force::Bool=false`: If true, re-download and overwrite the database directory.

# Returns
- `String`: Path to the database directory.
"""
function install_plassembler_db(db_dir::AbstractString; force::Bool=false)
    db_dir = abspath(String(db_dir))
    db_ready = _plassembler_db_ready(db_dir)
    if db_ready && !force
        @info "Plassembler database already present; skipping download." db_dir
        return db_dir
    end
    if isdir(db_dir) && (!db_ready || force)
        rm(db_dir; recursive=true, force=true)
    end
    mkpath(db_dir)
    Mycelia.add_bioconda_env("plassembler")
    try
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plassembler plassembler download -d $(db_dir) --force`)
    catch e
        tarballs = filter(path -> endswith(path, ".tar.gz"), readdir(db_dir; join=true))
        if isempty(tarballs)
            rethrow(e)
        end
        tarball = length(tarballs) == 1 ? tarballs[1] : sort(tarballs)[end]
        @warn "Plassembler download failed; attempting manual extraction." exception=(e, catch_backtrace())
        tmp_dir = mktempdir(dirname(db_dir))
        try
            Tar.extract(tarball, tmp_dir)
        catch extract_error
            rm(tmp_dir; recursive=true, force=true)
            rethrow(extract_error)
        end
        rm(db_dir; recursive=true, force=true)
        mv(tmp_dir, db_dir)
    end
    if !_plassembler_db_ready(db_dir)
        error("Plassembler database install failed: no .msh files found in $(db_dir).")
    end
    return db_dir
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Plassembler for hybrid plasmid assembly.

# Arguments
- `long_reads::AbstractString`: Long read FASTQ file.
- `read1::AbstractString`: Short read R1 FASTQ file.
- `read2::AbstractString`: Short read R2 FASTQ file.
- `chrom_size::Int`: Estimated chromosome size (bp).

# Keywords
- `outdir::AbstractString`: Output directory.
- `db_dir::Union{Nothing,AbstractString}=nothing`: Database directory (downloaded if missing).
- `threads::Int=get_default_threads()`: Thread count.
- `extra_args::Vector{String}=String[]`: Additional Plassembler CLI arguments.
- `force::Bool=false`: If true, re-run and overwrite existing outputs.

# Returns
Named tuple containing:
- `outdir::String`: Output directory.
- `plasmids_fasta::String`: Plasmid assemblies FASTA path.
- `chromosome_fasta::String`: Chromosome assembly FASTA path (if produced).
- `summary_tsv::String`: Summary TSV path.
"""
function run_plassembler(long_reads::AbstractString,
                         read1::AbstractString,
                         read2::AbstractString,
                         chrom_size::Int;
                         outdir::AbstractString,
                         db_dir::Union{Nothing,AbstractString}=nothing,
                         threads::Int=get_default_threads(),
                         extra_args::Vector{String}=String[],
                         force::Bool=false)
    isfile(long_reads) || error("Long read file not found: $(long_reads)")
    isfile(read1) || error("Read 1 file not found: $(read1)")
    isfile(read2) || error("Read 2 file not found: $(read2)")

    outdir = abspath(String(outdir))
    plasmids_fasta = joinpath(outdir, "plassembler_plasmids.fasta")
    chromosome_fasta = joinpath(outdir, "chromosome.fasta")
    summary_tsv = joinpath(outdir, "plassembler_summary.tsv")

    if !force && isfile(plasmids_fasta)
        return (;outdir, plasmids_fasta, chromosome_fasta, summary_tsv)
    end

    if force && isdir(outdir)
        rm(outdir; recursive=true, force=true)
    end
    mkpath(dirname(outdir))

    if db_dir === nothing || isempty(db_dir)
        db_dir = joinpath(homedir(), "workspace", ".plassembler")
    end
    if !_plassembler_db_ready(db_dir)
        install_plassembler_db(db_dir; force=force)
    end

    Mycelia.add_bioconda_env("plassembler")

    cmd_args = String[
        "run",
        "-l", long_reads,
        "-1", read1,
        "-2", read2,
        "-c", string(chrom_size),
        "-o", outdir,
        "-d", db_dir
    ]
    if !any(arg -> arg == "-t" || arg == "--threads" || startswith(arg, "--threads="), extra_args)
        push!(cmd_args, "-t", string(max(threads, 1)))
    end
    if force
        push!(cmd_args, "--force")
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plassembler plassembler $(cmd_args)`)
    return (;outdir, plasmids_fasta, chromosome_fasta, summary_tsv)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Plassembler for long-read-only plasmid assembly.

# Arguments
- `long_reads::AbstractString`: Long read FASTQ file.
- `chrom_size::Int`: Estimated chromosome size (bp).

# Keywords
- `outdir::AbstractString`: Output directory.
- `db_dir::Union{Nothing,AbstractString}=nothing`: Database directory (downloaded if missing).
- `threads::Int=get_default_threads()`: Thread count.
- `extra_args::Vector{String}=String[]`: Additional Plassembler CLI arguments.
- `force::Bool=false`: If true, re-run and overwrite existing outputs.

# Returns
Named tuple containing:
- `outdir::String`: Output directory.
- `plasmids_fasta::String`: Plasmid assemblies FASTA path.
- `chromosome_fasta::String`: Chromosome assembly FASTA path (if produced).
- `summary_tsv::String`: Summary TSV path.
"""
function run_plassembler_long(long_reads::AbstractString,
                              chrom_size::Int;
                              outdir::AbstractString,
                              db_dir::Union{Nothing,AbstractString}=nothing,
                              threads::Int=get_default_threads(),
                              extra_args::Vector{String}=String[],
                              force::Bool=false)
    isfile(long_reads) || error("Long read file not found: $(long_reads)")

    outdir = abspath(String(outdir))
    plasmids_fasta = joinpath(outdir, "plassembler_plasmids.fasta")
    chromosome_fasta = joinpath(outdir, "chromosome.fasta")
    summary_tsv = joinpath(outdir, "plassembler_summary.tsv")

    if !force && isfile(plasmids_fasta)
        return (;outdir, plasmids_fasta, chromosome_fasta, summary_tsv)
    end

    if force && isdir(outdir)
        rm(outdir; recursive=true, force=true)
    end
    mkpath(dirname(outdir))

    if db_dir === nothing || isempty(db_dir)
        db_dir = joinpath(homedir(), "workspace", ".plassembler")
    end
    if !_plassembler_db_ready(db_dir)
        install_plassembler_db(db_dir; force=force)
    end

    Mycelia.add_bioconda_env("plassembler")

    cmd_args = String[
        "long",
        "-l", long_reads,
        "-c", string(chrom_size),
        "-o", outdir,
        "-d", db_dir
    ]
    if !any(arg -> arg == "-t" || arg == "--threads" || startswith(arg, "--threads="), extra_args)
        push!(cmd_args, "-t", string(max(threads, 1)))
    end
    if force
        push!(cmd_args, "--force")
    end
    append!(cmd_args, extra_args)

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plassembler plassembler $(cmd_args)`)
    return (;outdir, plasmids_fasta, chromosome_fasta, summary_tsv)
end

# ============================================================================
# PLASS / PenguiN Assemblers (protein and nucleotide)
# ============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assemble proteins with PLASS.

# Arguments
- `reads1::String`: Path to first reads file (FASTQ/FASTA)
- `reads2::Union{String,Nothing}`: Optional paired-end mate
- `outdir::String`: Output directory (default: "plass_output")
- `min_seq_id::Union{Real,Nothing}`: Overlap identity threshold (e.g., 0.9)
- `min_length::Union{Int,Nothing}`: Minimum codon length for ORF prediction (default in tool: 40)
- `evalue::Union{Real,Nothing}`: E-value threshold for overlaps
- `num_iterations::Union{Int,Nothing}`: Number of assembly iterations
- `filter_proteins::Bool`: Whether to keep the neural network protein filter on (default true)
- `threads::Int`: Number of threads to use (default: `get_default_threads()`)

# Returns
Named tuple with:
- `outdir::String`: Output directory
- `assembly::String`: Protein assembly output path
- `tmpdir::String`: Temporary working directory used by PLASS
"""
function run_plass_assemble(;reads1::String, reads2::Union{String,Nothing}=nothing, outdir::String="plass_output",
    min_seq_id::Union{Real,Nothing}=nothing, min_length::Union{Int,Nothing}=nothing, evalue::Union{Real,Nothing}=nothing,
    num_iterations::Union{Int,Nothing}=nothing, filter_proteins::Bool=true, threads::Int=get_default_threads())

    Mycelia.add_bioconda_env("plass")

    if !isfile(reads1)
        error("reads1 not found: $reads1")
    end
    if !isnothing(reads2) && !isfile(reads2)
        error("reads2 not found: $reads2")
    end

    mkpath(outdir)
    tmpdir = joinpath(outdir, "tmp")
    mkpath(tmpdir)
    assembly_out = joinpath(outdir, "assembly.fas")

    cmd = ["plass", "assemble"]
    push!(cmd, reads1)
    if !isnothing(reads2)
        push!(cmd, reads2)
    end
    push!(cmd, assembly_out, tmpdir)

    if !isnothing(min_seq_id)
        push!(cmd, "--min-seq-id", string(min_seq_id))
    end
    if !isnothing(min_length)
        push!(cmd, "--min-length", string(min_length))
    end
    if !isnothing(evalue)
        push!(cmd, "-e", string(evalue))
    end
    if !isnothing(num_iterations)
        push!(cmd, "--num-iterations", string(num_iterations))
    end
    # Tool default is filter on; allow disabling
    if !filter_proteins
        push!(cmd, "--filter-proteins", "0")
    end
    push!(cmd, "--threads", string(threads))

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`)

    return (;outdir, assembly=assembly_out, tmpdir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assemble nucleotides with PenguiN using protein-guided mode.

# Arguments
- `reads1::String`: Path to first reads file (FASTQ/FASTA)
- `reads2::Union{String,Nothing}`: Optional paired-end mate
- `outdir::String`: Output directory (default: "penguin_guided_output")
- `threads::Int`: Number of threads to use (default: `get_default_threads()`)

# Returns
Named tuple with:
- `outdir::String`: Output directory
- `assembly::String`: Nucleotide assembly output path
- `tmpdir::String`: Temporary working directory
"""
function run_penguin_guided_nuclassemble(;reads1::String, reads2::Union{String,Nothing}=nothing, outdir::String="penguin_guided_output",
    threads::Int=get_default_threads())
    Mycelia.add_bioconda_env("plass")

    if !isfile(reads1)
        error("reads1 not found: $reads1")
    end
    if !isnothing(reads2) && !isfile(reads2)
        error("reads2 not found: $reads2")
    end

    mkpath(outdir)
    tmpdir = joinpath(outdir, "tmp")
    mkpath(tmpdir)
    assembly_out = joinpath(outdir, "assembly.fasta")

    cmd = ["penguin", "guided_nuclassemble"]
    push!(cmd, reads1)
    if !isnothing(reads2)
        push!(cmd, reads2)
    end
    push!(cmd, assembly_out, tmpdir)

    push!(cmd, "--threads", string(threads))

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`)

    return (;outdir, assembly=assembly_out, tmpdir)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Assemble nucleotides with PenguiN (nucleotide-only mode).

# Arguments
- `reads1::String`: Path to first reads file (FASTQ/FASTA)
- `reads2::Union{String,Nothing}`: Optional paired-end mate
- `outdir::String`: Output directory (default: "penguin_output")
- `threads::Int`: Number of threads to use (default: `get_default_threads()`)

# Returns
Named tuple with:
- `outdir::String`: Output directory
- `assembly::String`: Nucleotide assembly output path
- `tmpdir::String`: Temporary working directory
"""
function run_penguin_nuclassemble(;reads1::String, reads2::Union{String,Nothing}=nothing, outdir::String="penguin_output",
    threads::Int=get_default_threads())
    Mycelia.add_bioconda_env("plass")

    if !isfile(reads1)
        error("reads1 not found: $reads1")
    end
    if !isnothing(reads2) && !isfile(reads2)
        error("reads2 not found: $reads2")
    end

    mkpath(outdir)
    tmpdir = joinpath(outdir, "tmp")
    mkpath(tmpdir)
    assembly_out = joinpath(outdir, "assembly.fasta")

    cmd = ["penguin", "nuclassemble"]
    push!(cmd, reads1)
    if !isnothing(reads2)
        push!(cmd, reads2)
    end
    push!(cmd, assembly_out, tmpdir)

    push!(cmd, "--threads", string(threads))

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`)

    return (;outdir, assembly=assembly_out, tmpdir)
end

# ISN'T AVAILABLE IN BIOCONDA
# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Run QuickMerge for assembly merging and scaffolding.

# # Arguments
# - `self_assembly::String`: Path to self assembly FASTA file (primary assembly)
# - `ref_assembly::String`: Path to reference assembly FASTA file (hybrid assembly)
# - `outdir::String`: Output directory path (default: "quickmerge_output")
# - `hco::Float64`: HCO threshold for overlap detection (default: 5.0)
# - `c::Float64`: Coverage cutoff for contig filtering (default: 1.5)
# - `l::Int`: Minimum alignment length (default: 5000)

# # Returns
# Named tuple containing:
# - `outdir::String`: Path to output directory
# - `merged_assembly::String`: Path to merged assembly file

# # Details
# - Uses MUMmer-based approach for identifying high-confidence overlaps
# - Merges assemblies by splicing contigs at overlap boundaries
# - Optimized for merging complementary assemblies (e.g., short+long read)
# - Automatically creates and uses a conda environment with quickmerge
# - Skips analysis if output files already exist
# """
# function run_quickmerge(self_assembly::String, ref_assembly::String; outdir::String="quickmerge_output", hco::Float64=5.0, c::Float64=1.5, l::Int=5000)
#     Mycelia.add_bioconda_env("bioconda::quickmerge")
#     mkpath(outdir)
    
#     merged_assembly = joinpath(outdir, "merged.fasta")
    
#     if !isfile(merged_assembly)
#         # Run QuickMerge
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n quickmerge merge_wrapper.py $(self_assembly) $(ref_assembly) -pre $(outdir)/merged -hco $(hco) -c $(c) -l $(l)`)
#     end
    
#     return (;outdir, merged_assembly)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Apollo for HMM-based assembly polishing.

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `reads_file::String`: Path to reads FASTQ file (can be comma-separated for paired reads)
- `outdir::String`: Output directory path (default: "\${assembly_file}_apollo")

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `polished_assembly::String`: Path to polished assembly file

# Details
- Uses profile Hidden Markov Models for technology-independent polishing
- Constructs HMM graphs for consensus sequence correction
- Works with both short and long read technologies
- Automatically creates and uses a conda environment with apollo
- Falls back to Racon polishing when the Apollo CLI is unavailable
- Skips analysis if output files already exist
"""
function run_apollo(assembly_file::String, reads_file::String; outdir::String=assembly_file * "_apollo")
    Mycelia.add_bioconda_env("apollo")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    mkpath(outdir)
    
    basename_assembly = splitext(basename(assembly_file))[1]
    polished_assembly = joinpath(outdir, basename_assembly * "_polished.fasta")
    
    apollo_prefix = _conda_env_prefix("apollo")
    apollo_bin = apollo_prefix === nothing ? nothing : joinpath(apollo_prefix, "bin", "apollo")
    apollo_available = apollo_bin !== nothing && isfile(apollo_bin)

    if !isfile(polished_assembly)
        if apollo_available
            # Map reads to assembly first
            bam_file = joinpath(outdir, basename_assembly * ".bam")
            if !isfile(bam_file)
                threads = get_default_threads()
                minimap_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -ax map-pb $(assembly_file) $(reads_file)`
                samtools_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools sort -@ $(threads) -o $(bam_file)`
                run(pipeline(minimap_cmd, samtools_cmd))
                run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n samtools samtools index $(bam_file)`)
            end

            # Run Apollo polishing
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n apollo apollo polish -a $(assembly_file) -b $(bam_file) -o $(polished_assembly)`)
        else
            unicycler_prefix = _conda_env_prefix("unicycler")
            if unicycler_prefix !== nothing && isfile(joinpath(unicycler_prefix, "bin", "racon"))
                racon_env = "unicycler"
            else
                Mycelia.add_bioconda_env("racon")
                racon_env = "racon"
            end
            overlaps_file = joinpath(outdir, basename_assembly * ".paf")
            if !isfile(overlaps_file)
                minimap_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n minimap2 minimap2 -x map-pb $(assembly_file) $(reads_file)`
                run(pipeline(minimap_cmd, overlaps_file))
            end
            racon_cmd = `$(Mycelia.CONDA_RUNNER) run --live-stream -n $(racon_env) racon $(reads_file) $(overlaps_file) $(assembly_file)`
            open(polished_assembly, "w") do io
                run(pipeline(racon_cmd, stdout=io))
            end
        end
    end
    
    return (;outdir, polished_assembly)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Homopolish for reference-based homopolymer error correction.

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `local_db_path::String`: Path to local FASTA database of homologous genomes
- `outdir::String`: Output directory path (default: "\${assembly_file}_homopolish")
- `model_path::String`: Path to Homopolish model (default: auto-detected)
- `sketch_path::String`: Path to a Mash sketch (optional; generated from `local_db_path`)
- `threads::Int`: Number of threads to use for Homopolish (default: `get_default_threads()`)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `polished_assembly::String`: Path to polished assembly file

# Details
- Uses reference genomes and neural networks for homopolymer correction
- Specifically targets homopolymer run errors common in long-read sequencing
- Auto-selects the PacBio CLR model when available
- Builds a Mash sketch when `sketch_path` is not provided
- Automatically creates and uses a conda environment with homopolish
- Utilizes requested CPU threads for polishing
- Skips analysis if output files already exist
"""
function _conda_env_prefix(env_name::String)
    env_lines = try
        readlines(`$(Mycelia.CONDA_RUNNER) env list`)
    catch
        return nothing
    end

    for line in env_lines
        stripped = strip(line)
        if isempty(stripped) || startswith(stripped, "#")
            continue
        end
        parts = split(stripped)
        if !isempty(parts) && parts[1] == env_name && length(parts) >= 2
            return parts[end]
        end
    end

    return nothing
end

function _find_homopolish_model_path(model_name::String)
    env_prefix = _conda_env_prefix("homopolish")
    if env_prefix === nothing
        return nothing
    end

    for (root, _, files) in walkdir(env_prefix)
        if model_name in files
            return joinpath(root, model_name)
        end
    end

    return nothing
end

function run_homopolish(assembly_file::String, local_db_path::String; outdir::String=assembly_file * "_homopolish", model_path::String="", threads = get_default_threads(), sketch_path::String="")
    Mycelia.add_bioconda_env("homopolish")
    mkpath(outdir)
    
    basename_assembly = splitext(basename(assembly_file))[1]
    polished_assembly = joinpath(outdir, basename_assembly * "_homopolished.fasta")
    
    if isempty(model_path)
        model_path = _find_homopolish_model_path("pb.pkl")
        if model_path === nothing
            error("Homopolish model not found in homopolish environment; pass model_path explicitly.")
        end
    end

    if isempty(sketch_path)
        if isempty(local_db_path)
            error("Homopolish requires sketch_path or local_db_path pointing to a FASTA file.")
        end
        Mycelia.add_bioconda_env("mash")
        sketch_prefix = joinpath(outdir, basename_assembly * "_homopolish_db")
        sketch_path = sketch_prefix * ".msh"
        if !isfile(sketch_path)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash sketch -o $(sketch_prefix) $(local_db_path)`)
        end
    end

    if !isfile(polished_assembly)
        cmd_args = ["homopolish", "polish", "-a", assembly_file, "-s", sketch_path, "-m", model_path, "-o", outdir, "-t", string(threads)]
        
        # Run Homopolish
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n homopolish $(cmd_args)`)
        
        # Move output to expected location if needed
        homopolish_output = joinpath(outdir, "homopolished_" * basename_assembly)
        if isfile(homopolish_output) && !isfile(polished_assembly)
            mv(homopolish_output, polished_assembly)
        end
    end
    
    return (;outdir, polished_assembly)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run HyLight for hybrid strain-resolved metagenomic assembly.

# Arguments
- `short_reads_1::String`: Path to first short read FASTQ file
- `short_reads_2::String`: Path to second short read FASTQ file  
- `long_reads::String`: Path to long read FASTQ file
- `outdir::String`: Output directory path (default: "hylight_output")
- `threads::Int`: Number of threads to use (default: `get_default_threads()`)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `strain_assemblies::String`: Path to strain-resolved assemblies directory

# Details
- Uses hybrid overlap graphs combining short and long reads
- Implements "cross hybrid" mutual support strategy for strain resolution
- Separates closely related strains in metagenomic samples
- Automatically creates and uses a conda environment with hylight
- Utilizes requested CPU threads for HyLight run
- Skips assembly if output directory already exists
"""
function run_hylight(short_reads_1::String, short_reads_2::String, long_reads::String; outdir::String="hylight_output", threads = get_default_threads())
    isfile(short_reads_1) || error("Short read 1 file not found: $(short_reads_1)")
    isfile(short_reads_2) || error("Short read 2 file not found: $(short_reads_2)")
    isfile(long_reads) || error("Long read file not found: $(long_reads)")

    Mycelia.add_bioconda_env("hylight")
    mkpath(outdir)
    
    strain_assemblies = joinpath(outdir, "strain_assemblies")
    
    if !isdir(strain_assemblies)
        # Run HyLight assembly
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hylight hylight -1 $(short_reads_1) -2 $(short_reads_2) -l $(long_reads) -o $(outdir) -t $(threads)`)
    end
    
    return (;outdir, strain_assemblies)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run STRONG for strain resolution on assembly graphs.

# Arguments
- `assembly_graph::String`: Path to assembly graph file (GFA format)
- `reads_file::String`: Path to reads FASTQ file (can be comma-separated for paired reads)
- `outdir::String`: Output directory path (default: "strong_output")
- `nb_strains::Int`: Expected number of strains (default: auto-detect)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `strain_unitigs::String`: Path to strain-resolved unitigs

# Details
- Uses BayesPaths algorithm for haplotype resolution on assembly graphs
- Implements Bayesian approach to separate strain variants
- Supports multi-sample strain analysis
- Automatically creates and uses a conda environment with strong
- Skips analysis if output files already exist
"""
function run_strong(assembly_graph::String, reads_file::String; outdir::String="strong_output", nb_strains::Int=0)
    Mycelia.add_bioconda_env("strong")
    mkpath(outdir)
    
    strain_unitigs = joinpath(outdir, "strain_unitigs.fasta")
    
    if !isfile(strain_unitigs)
        cmd_args = ["STRONG", assembly_graph, reads_file, "-o", outdir]
        
        if nb_strains > 0
            push!(cmd_args, "-s", string(nb_strains))
        end
        
        # Run STRONG
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strong $(cmd_args)`)
    end
    
    return (;outdir, strain_unitigs)
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Run Ray assembler for de novo genome assembly.

# # Arguments
# - `reads_files::Vector{String}`: Vector of FASTQ file paths (supports paired-end and single-end)
# - `outdir::String`: Output directory path (default: "ray_output")
# - `k::Int`: K-mer size for assembly (default: 31)
# - `min_contig_length::Int`: Minimum contig length (default: 200)

# # Returns
# Named tuple containing:
# - `outdir::String`: Path to output directory
# - `contigs::String`: Path to contigs file
# - `scaffolds::String`: Path to scaffolds file

# # Details
# - Uses Ray's distributed de Bruijn graph approach
# - Supports both single-end and paired-end reads
# - Automatically creates and uses a conda environment with ray
# - Skips assembly if output directory already exists
# - Utilizes all available CPU threads
# """
# function run_ray(reads_files::Vector{String}; outdir::String="ray_output", k::Int=31, min_contig_length::Int=200)
#     Mycelia.add_bioconda_env("ray")
#     mkpath(outdir)
    
#     contigs_file = joinpath(outdir, "Contigs.fasta")
#     scaffolds_file = joinpath(outdir, "Scaffolds.fasta")
    
#     if !isfile(contigs_file)
#         # Build Ray command arguments
#         cmd_args = ["Ray", "-k", string(k), "-o", outdir, "-minimum-contig-length", string(min_contig_length)]
        
#         # Add reads files - Ray auto-detects paired vs single-end
#         if length(reads_files) == 2
#             # Paired-end reads
#             push!(cmd_args, "-p", reads_files[1], reads_files[2])
#         elseif length(reads_files) == 1
#             # Single-end reads  
#             push!(cmd_args, "-s", reads_files[1])
#         else
#             # Multiple libraries - treat as single-end
#             for reads_file in reads_files
#                 push!(cmd_args, "-s", reads_file)
#             end
#         end
        
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ray $(cmd_args)`)
#     end
    
#     return (;outdir, contigs=contigs_file, scaffolds=scaffolds_file)
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Run Ray Meta for metagenomic assembly.

# # Arguments
# - `reads_files::Vector{String}`: Vector of FASTQ file paths (supports paired-end and single-end)
# - `outdir::String`: Output directory path (default: "ray_meta_output")
# - `k::Int`: K-mer size for assembly (default: 31)
# - `min_contig_length::Int`: Minimum contig length (default: 200)
# - `enable_communities::Bool`: Enable community detection (default: true)

# # Returns
# Named tuple containing:
# - `outdir::String`: Path to output directory
# - `contigs::String`: Path to contigs file
# - `scaffolds::String`: Path to scaffolds file

# # Details
# - Uses Ray Meta's distributed approach for metagenomic data
# - Supports community detection for binning related sequences
# - Handles uneven coverage typical in metagenomic samples
# - Automatically creates and uses a conda environment with ray
# - Skips assembly if output directory already exists
# """
# function run_ray_meta(reads_files::Vector{String}; outdir::String="ray_meta_output", k::Int=31, min_contig_length::Int=200, enable_communities::Bool=true)
#     Mycelia.add_bioconda_env("ray")
#     mkpath(outdir)
    
#     contigs_file = joinpath(outdir, "Contigs.fasta")
#     scaffolds_file = joinpath(outdir, "Scaffolds.fasta")
    
#     if !isfile(contigs_file)
#         # Build Ray Meta command arguments
#         cmd_args = ["Ray", "-k", string(k), "-o", outdir, "-minimum-contig-length", string(min_contig_length)]
        
#         # Enable metagenomic mode
#         push!(cmd_args, "-enable-neighbourhoods")
        
#         if enable_communities
#             push!(cmd_args, "-enable-neighbourhoods")
#         end
        
#         # Add reads files
#         if length(reads_files) == 2
#             # Paired-end reads
#             push!(cmd_args, "-p", reads_files[1], reads_files[2])
#         elseif length(reads_files) == 1
#             # Single-end reads
#             push!(cmd_args, "-s", reads_files[1])
#         else
#             # Multiple libraries
#             for reads_file in reads_files
#                 push!(cmd_args, "-s", reads_file)
#             end
#         end
        
#         run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ray $(cmd_args)`)
#     end
    
#     return (;outdir, contigs=contigs_file, scaffolds=scaffolds_file)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Velvet assembler for short-read genome assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "velvet_output")
- `k::Int`: K-mer size for assembly (default: 31)
- `exp_cov::String`: Expected coverage (default: "auto")
- `min_contig_lgth::Int`: Minimum contig length (default: 200)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file

# Details
- Uses Velvet's de Bruijn graph approach for short reads
- Two-step process: velveth (indexing) + velvetg (assembly)
- Optimized for Illumina short reads
- Automatically creates and uses a conda environment with velvet
- Skips assembly if output directory already exists
"""
function run_velvet(fastq1::String; fastq2::Union{String,Nothing}=nothing, outdir::String="velvet_output", k::Int=31, exp_cov::String="auto", min_contig_lgth::Int=200)
    Mycelia.add_bioconda_env("velvet")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "contigs.fa")
    
    if !isfile(contigs_file)
        # Step 1: velveth (k-mer hashing and indexing)
        if isnothing(fastq2)
            # Single-end reads
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet velveth $(outdir) $(k) -short -fastq $(fastq1)`)
        else
            # Paired-end reads
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet velveth $(outdir) $(k) -shortPaired -fastq -separate $(fastq1) $(fastq2)`)
        end
        
        # Step 2: velvetg (graph construction and traversal)
        cmd_args = ["velvetg", outdir, "-min_contig_lgth", string(min_contig_lgth)]
        
        if exp_cov != "auto"
            push!(cmd_args, "-exp_cov", exp_cov)
        end
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet $(cmd_args)`)
    end
    
    return (;outdir, contigs=contigs_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaVelvet assembler for metagenomic short-read assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "metavelvet_output")
- `k::Int`: K-mer size for assembly (default: 31)
- `exp_cov::String`: Expected coverage (default: "auto")
- `min_contig_lgth::Int`: Minimum contig length (default: 200)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to assembled contigs file

# Details
- Uses metaVelvet's approach for metagenomic data with varying coverage
- Three-step process: velveth (indexing) + velvetg (read tracking) + meta-velvetg
- Designed for mixed community samples with uneven coverage
- Automatically creates and uses a conda environment with velvet
- Skips assembly if output directory already exists
"""
function run_metavelvet(fastq1::String; fastq2::Union{String,Nothing}=nothing, outdir::String="metavelvet_output", k::Int=31, exp_cov::String="auto", min_contig_lgth::Int=200)
    Mycelia.add_bioconda_env("velvet")
    Mycelia.add_bioconda_env("metavelvet")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "meta-velvetg.contigs.fa")
    
    if !isfile(contigs_file)
        # Step 1: velveth (k-mer hashing and indexing) - same as regular velvet
        if isnothing(fastq2)
            # Single-end reads
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet velveth $(outdir) $(k) -short -fastq $(fastq1)`)
        else
            # Paired-end reads
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet velveth $(outdir) $(k) -shortPaired -fastq -separate $(fastq1) $(fastq2)`)
        end
        
        # Step 2: velvetg with read tracking (required by meta-velvetg)
        velvetg_args = ["velvetg", outdir, "-read_trkg", "yes", "-min_contig_lgth", string(min_contig_lgth)]
        if exp_cov != "auto"
            push!(velvetg_args, "-exp_cov", exp_cov)
        end
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet $(velvetg_args)`)

        # Step 3: meta-velvetg (metagenomic graph construction and traversal)
        cmd_args = ["meta-velvetg", outdir, "-min_contig_lgth", string(min_contig_lgth)]
        
        if exp_cov != "auto"
            push!(cmd_args, "-exp_cov", exp_cov)
        end
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n metavelvet $(cmd_args)`)
    end
    
    return (;outdir, contigs=contigs_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run VelvetG with optimized parameters for genome assembly.

# Arguments
- `velvet_dir::String`: Path to velveth output directory
- `outdir::String`: Output directory path (default: velvet_dir)
- `exp_cov::Union{String,Float64}`: Expected coverage (default: "auto")
- `cov_cutoff::Union{String,Float64}`: Coverage cutoff (default: "auto")  
- `min_contig_lgth::Int`: Minimum contig length (default: 200)
- `scaffolding::Bool`: Enable scaffolding (default: true)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `stats::String`: Path to assembly statistics file

# Details
- Advanced velvetg run with parameter optimization
- Supports scaffolding and coverage optimization
- Can be run on existing velveth output
- Automatically creates and uses a conda environment with velvet
"""
function run_velvetg(velvet_dir::String; outdir::String=velvet_dir, exp_cov::Union{String,Float64}="auto", cov_cutoff::Union{String,Float64}="auto", min_contig_lgth::Int=200, scaffolding::Bool=true)
    Mycelia.add_bioconda_env("velvet")
    
    contigs_file = joinpath(outdir, "contigs.fa")
    stats_file = joinpath(outdir, "stats.txt")
    
    # Build velvetg command
    cmd_args = ["velvetg", velvet_dir, "-min_contig_lgth", string(min_contig_lgth)]
    
    if exp_cov != "auto"
        push!(cmd_args, "-exp_cov", string(exp_cov))
    end
    
    if cov_cutoff != "auto"
        push!(cmd_args, "-cov_cutoff", string(cov_cutoff))
    end
    
    if scaffolding
        push!(cmd_args, "-scaffolding", "yes")
    else
        push!(cmd_args, "-scaffolding", "no")
    end
    
    # Copy output to specified directory if different
    if outdir != velvet_dir
        mkpath(outdir)
        push!(cmd_args, "-read_trkg", "yes")  # Enable read tracking for copying
    end
    
    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n velvet $(cmd_args)`)
    
    # Copy results if needed
    if outdir != velvet_dir
        for file in ["contigs.fa", "stats.txt", "Graph", "LastGraph"]
            src_file = joinpath(velvet_dir, file)
            dst_file = joinpath(outdir, file)
            if isfile(src_file)
                cp(src_file, dst_file, force=true)
            end
        end
    end
    
    return (;outdir, contigs=contigs_file, stats=stats_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaMDBG assembler for metagenomic long-read assembly.

# Arguments
- `hifi_reads::Union{String,Vector{String},Nothing}`: Path(s) to HiFi (PacBio) read files (default: nothing)
- `ont_reads::Union{String,Vector{String},Nothing}`: Path(s) to ONT (Nanopore) read files (default: nothing)
- `outdir::String`: Output directory path (default: "metamdbg_output")
- `abundance_min::Int`: Minimum abundance threshold (default: 3)
- `threads::Int`: Number of threads to use (default: get_default_threads())
- `graph_k::Int`: K-mer resolution level for graph generation (default: 21)

Note: Must provide either `hifi_reads`, `ont_reads`, or both. Cannot be both nothing.
Graph generation: Automatically generates assembly graphs using the specified k-mer resolution.

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `contigs::String`: Path to contigs file
- `graph::String`: Path to assembly graph file

# Details
- Uses metaMDBG's minimizer-space de Bruijn graphs for metagenomic assembly
- Specifically designed for long-read metagenomic data
- Handles species abundance variation and strain diversity
- Produces both contigs and assembly graphs
- Automatically creates and uses a conda environment with metamdbg
- Skips assembly if output directory already exists
"""
function run_metamdbg(; hifi_reads::Union{String,Vector{String},Nothing}=nothing,
                      ont_reads::Union{String,Vector{String},Nothing}=nothing,
                      outdir::String="metamdbg_output",
                      abundance_min::Int=3,
                      threads::Int=get_default_threads(),
                      graph_k::Int=21)
    
    # Validate input - must have at least one read type
    if isnothing(hifi_reads) && isnothing(ont_reads)
        error("Must provide either hifi_reads, ont_reads, or both")
    end
    
    Mycelia.add_bioconda_env("metamdbg")
    mkpath(outdir)
    
    contigs_file = joinpath(outdir, "contigs.fasta")
    graph_file = joinpath(outdir, "graph.gfa")
    
    if !isfile(contigs_file)
        # Build command arguments
        cmd_args = ["metaMDBG", "asm", "--out-dir", outdir]
        
        # Add HiFi reads if provided
        if !isnothing(hifi_reads)
            push!(cmd_args, "--in-hifi")
            if isa(hifi_reads, String)
                push!(cmd_args, hifi_reads)
            else
                append!(cmd_args, hifi_reads)
            end
        end
        
        # Add ONT reads if provided
        if !isnothing(ont_reads)
            push!(cmd_args, "--in-ont")
            if isa(ont_reads, String)
                push!(cmd_args, ont_reads)
            else
                append!(cmd_args, ont_reads)
            end
        end
        
        # Add other parameters
        push!(cmd_args, "--min-abundance", string(abundance_min))
        push!(cmd_args, "--threads", string(threads))
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n metamdbg $(cmd_args)`)
        
        # Generate assembly graph using metaMDBG gfa command
        gfa_cmd_args = [
            "metaMDBG", "gfa",
            "--assembly-dir", outdir,
            "--k", string(graph_k),
            "--threads", string(threads)
        ]
        
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n metamdbg $(gfa_cmd_args)`)
        
        # Set expected output file names
        contigs_file = joinpath(outdir, "contigs.fasta.gz")  # metaMDBG compresses output
        
        # Find the generated graph file (format: assemblyGraph_k{k}_{length}bps.gfa)
        graph_files = filter(f -> startswith(basename(f), "assemblyGraph_k") && endswith(f, ".gfa"), 
                           readdir(outdir, join=true))
        
        if isempty(graph_files)
            error("metaMDBG failed to generate assembly graph file. Expected file pattern: assemblyGraph_k*.gfa in $outdir")
        end
        
        graph_file = first(graph_files)
    end
    
    return (;outdir, contigs=contigs_file, graph=graph_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run Strainy for strain phasing from long reads.

# Arguments
- `assembly_file::String`: Path to assembly FASTA file
- `long_reads::String`: Path to long read FASTQ file
- `outdir::String`: Output directory path (default: "strainy_output")
- `mode::String`: Analysis mode ("transform" or "phase", default: "phase")
- `threads::Int`: Number of threads to use for mapping and phasing (default: `get_default_threads()`)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `strain_assemblies::String`: Path to strain-phased assemblies

# Details
- Uses connection graph-based read clustering for strain separation
- Implements long-read phasing to resolve strain-level variants
- Generates strain unitigs and simplified assembly graphs
- Automatically creates and uses a conda environment with strainy
- Utilizes requested CPU threads for read mapping
- Skips analysis if output files already exist
"""
function run_strainy(assembly_file::String, long_reads::String; outdir::String="strainy_output", mode::String="phase", threads = get_default_threads())
    Mycelia.add_bioconda_env("strainy")
    mkpath(outdir)
    
    strain_assemblies = joinpath(outdir, "strain_assemblies.fasta")
    
    if !isfile(strain_assemblies)
        # First map reads to assembly
        bam_file = joinpath(outdir, "mapped_reads.bam")
        if !isfile(bam_file)
            run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy minimap2 -ax map-ont $(assembly_file) $(long_reads)`, `$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy samtools sort -@ $(threads) -o $(bam_file)`))
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy samtools index $(bam_file)`)
        end
        
        # Run Strainy
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy strainy --bam $(bam_file) --fasta $(assembly_file) --output $(outdir) --mode $(mode)`)
    end
    
    return (;outdir, strain_assemblies)
end
