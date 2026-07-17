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
function run_megahit(;
        fastq1,
        fastq2 = nothing,
        outdir = nothing,
        min_contig_len = 200,
        k_list = "21,29,39,59,79,99,119,141",
        threads::Int = get_default_threads(),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "megahit",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing,
        ensure_env::Bool = true)
    if ensure_env
        Mycelia.add_bioconda_env("megahit")
    end
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

    if executor !== nothing
        command_parts = megahit_cmd(
            fastq1 = fastq1,
            fastq2 = fastq2,
            outdir = outdir,
            min_contig_len = min_contig_len,
            k_list = k_list,
            threads = threads
        )
        job = Mycelia.build_execution_job(
            cmd = command_parts.cmd,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (;
            outdir = command_parts.outdir,
            contigs = command_parts.contigs,
            fastg = command_parts.fastg,
            gfa = command_parts.gfa
        )
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
            rm(outdir, recursive = true)
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
        run(pipeline(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n megahit megahit_core contig2fastg $(final_k) $(contigs_path)`,
            fastg_path))
        @assert isfile(fastg_path)
        @assert filesize(fastg_path) > 0
    end

    gfa_path = fastg_path * ".gfa"
    if !isfile(gfa_path)
        # if isfile(fastg_path) && (filesize(fastg_path) > 0)
        if ensure_env
            Mycelia.add_bioconda_env("gfatools")
        end
        run(pipeline(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n gfatools gfatools view $(fastg_path)`,
            gfa_path))
        # else
        #     @warn "final.contigs.fastg not found, skipping gfatools step"
        #     gfa_path = missing
        # end
    end

    return (; outdir, contigs = contigs_path, fastg = fastg_path, gfa = gfa_path)
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
        fastq2::Union{Nothing, String} = nothing,
        outdir::Union{Nothing, String} = nothing,
        min_contig_len::Int = 200,
        k_list::String = "21,29,39,59,79,99,119,141",
        threads::Int = get_default_threads()
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
    fastg_path = replace(contigs_path, ".fa" => ".fastg")
    gfa_path = fastg_path * ".gfa"

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
    return (; cmd, outdir, contigs = contigs_path, fastg = fastg_path, gfa = gfa_path)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaSPAdes assembler for metagenomic short read assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "metaspades_output")
- `k_list::String`: k-mer sizes to use (default: "21,33,55,77")
- `only_assembler::Bool`: Skip read correction (default: false)

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
function run_metaspades(;
        fastq1,
        fastq2 = nothing,
        outdir = "metaspades_output",
        k_list = "21,33,55,77",
        threads::Int = get_default_threads(),
        only_assembler::Bool = false,
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "metaspades",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("spades")
    mkpath(outdir)

    if executor !== nothing
        cmd_args = ["metaspades.py", "-o", outdir, "-k", k_list, "-t", string(threads)]
        only_assembler && push!(cmd_args, "--only-assembler")
        if isnothing(fastq2)
            push!(cmd_args, "-s", fastq1)
        else
            push!(cmd_args, "-1", fastq1, "-2", fastq2)
        end
        cmd = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n spades $(cmd_args)`
        )
        script = join([
            "set -euo pipefail",
            "mkdir -p \"$(outdir)\"",
            "if [ ! -f \"$(joinpath(outdir, "contigs.fasta"))\" ]; then",
            "  $(cmd)",
            "fi"
        ], "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, contigs = joinpath(outdir, "contigs.fasta"),
            scaffolds = joinpath(outdir, "scaffolds.fasta"),
            graph = joinpath(outdir, "assembly_graph_with_scaffolds.gfa"))
    end

    if !isfile(joinpath(outdir, "contigs.fasta"))
        cmd_args = ["metaspades.py", "-o", outdir, "-k", k_list, "-t", string(threads)]
        only_assembler && push!(cmd_args, "--only-assembler")
        if isnothing(fastq2)
            push!(cmd_args, "-s", fastq1)
        else
            push!(cmd_args, "-1", fastq1, "-2", fastq2)
        end
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades $(cmd_args)`)
    end
    return (; outdir, contigs = joinpath(outdir, "contigs.fasta"),
        scaffolds = joinpath(outdir, "scaffolds.fasta"),
        graph = joinpath(outdir, "assembly_graph_with_scaffolds.gfa"))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run SPAdes assembler for single genome isolate assembly.

# Arguments
- `fastq1::String`: Path to first paired-end FASTQ file
- `fastq2::String`: Path to second paired-end FASTQ file (optional for single-end)
- `outdir::String`: Output directory path (default: "spades_output")
- `k_list::String`: k-mer sizes to use (default: "21,33,55,77")
- `only_assembler::Bool`: Skip read correction (default: false)

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
function run_spades(;
        fastq1,
        fastq2 = nothing,
        outdir = "spades_output",
        k_list = "21,33,55,77",
        threads::Int = get_default_threads(),
        only_assembler::Bool = false,
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "spades",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("spades")
    mkpath(outdir)

    if executor !== nothing
        cmd_args = ["spades.py", "-o", outdir, "-k", k_list, "-t", string(threads)]
        only_assembler && push!(cmd_args, "--only-assembler")
        if isnothing(fastq2)
            push!(cmd_args, "-s", fastq1)
        else
            push!(cmd_args, "-1", fastq1, "-2", fastq2)
        end
        cmd = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n spades $(cmd_args)`
        )
        script = join([
            "set -euo pipefail",
            "mkdir -p \"$(outdir)\"",
            "if [ ! -f \"$(joinpath(outdir, "contigs.fasta"))\" ]; then",
            "  $(cmd)",
            "fi"
        ], "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, contigs = joinpath(outdir, "contigs.fasta"),
            scaffolds = joinpath(outdir, "scaffolds.fasta"),
            graph = joinpath(outdir, "assembly_graph_with_scaffolds.gfa"))
    end

    if !isfile(joinpath(outdir, "contigs.fasta"))
        cmd_args = ["spades.py", "-o", outdir, "-k", k_list, "-t", string(threads)]
        only_assembler && push!(cmd_args, "--only-assembler")
        if isnothing(fastq2)
            push!(cmd_args, "-s", fastq1)
        else
            push!(cmd_args, "-1", fastq1, "-2", fastq2)
        end
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n spades $(cmd_args)`)
    end
    return (; outdir, contigs = joinpath(outdir, "contigs.fasta"),
        scaffolds = joinpath(outdir, "scaffolds.fasta"),
        graph = joinpath(outdir, "assembly_graph_with_scaffolds.gfa"))
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
function run_skesa(;
        fastq1,
        fastq2 = nothing,
        outdir = "skesa_output",
        min_contig_len = 200,
        threads::Int = get_default_threads(),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "skesa",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("skesa")
    mkpath(outdir)

    contigs_file = joinpath(outdir, "contigs.fa")
    if executor !== nothing
        cmd = if isnothing(fastq2)
            Mycelia.command_string(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n skesa skesa --reads $(fastq1) --contigs_out $(contigs_file) --cores $(threads) --min_contig $(min_contig_len)`
            )
        else
            Mycelia.command_string(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n skesa skesa --reads $(fastq1),$(fastq2) --contigs_out $(contigs_file) --cores $(threads) --min_contig $(min_contig_len)`
            )
        end
        script = join([
            "set -euo pipefail",
            "mkdir -p \"$(outdir)\"",
            "if [ ! -f \"$(contigs_file)\" ]; then",
            "  $(cmd)",
            "fi"
        ], "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, contigs = contigs_file)
    end

    if !isfile(contigs_file)
        if isnothing(fastq2)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n skesa skesa --reads $(fastq1) --contigs_out $(contigs_file) --cores $(threads) --min_contig $(min_contig_len)`)
        else
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n skesa skesa --reads $(fastq1),$(fastq2) --contigs_out $(contigs_file) --cores $(threads) --min_contig $(min_contig_len)`)
        end
    end
    return (; outdir, contigs = contigs_file)
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
function run_flye(;
        fastq,
        outdir = "flye_output",
        genome_size = nothing,
        read_type = "pacbio-hifi",
        min_overlap = nothing,
        threads::Int = min(get_default_threads(), FLYE_MAX_THREADS),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "flye",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("flye")
    mkpath(outdir)

    if executor !== nothing
        cmd_args = ["flye", "--$(read_type)", fastq, "--out-dir",
            outdir, "--threads", string(threads)]
        if !isnothing(genome_size)
            push!(cmd_args, "--genome-size", string(genome_size))
        end
        if !isnothing(min_overlap)
            push!(cmd_args, "--min-overlap", string(min_overlap))
        end
        cmd = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n flye $(cmd_args)`
        )
        script = join([
            "set -euo pipefail",
            "mkdir -p \"$(outdir)\"",
            "if [ ! -f \"$(joinpath(outdir, "assembly.fasta"))\" ]; then",
            "  $(cmd)",
            "fi"
        ], "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, assembly = joinpath(outdir, "assembly.fasta"),
            graph = joinpath(outdir, "assembly_graph.gfa"))
    end

    if !isfile(joinpath(outdir, "assembly.fasta"))
        cmd_args = ["flye", "--$(read_type)", fastq, "--out-dir",
            outdir, "--threads", string(threads)]
        if !isnothing(genome_size)
            push!(cmd_args, "--genome-size", string(genome_size))
        end
        if !isnothing(min_overlap)
            push!(cmd_args, "--min-overlap", string(min_overlap))
        end
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n flye $(cmd_args)`)
    end
    return (; outdir, assembly = joinpath(outdir, "assembly.fasta"),
        graph = joinpath(outdir, "assembly_graph.gfa"))
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
function run_metaflye(;
        fastq::String,
        outdir::String = "metaflye_output",
        genome_size::Union{Nothing, Integer, AbstractString} = nothing,
        read_type::String = "pacbio-hifi",
        meta::Bool = true,
        min_overlap::Union{Nothing, Int} = nothing,
        iterations::Int = 1,
        keep_haplotypes::Bool = false,
        no_alt_contigs::Bool = false,
        threads::Int = min(get_default_threads(), FLYE_MAX_THREADS),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "metaflye",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("flye")
    mkpath(outdir)

    if executor !== nothing
        cmd_args = ["flye", "--$(read_type)", fastq, "--out-dir",
            outdir, "--threads", string(threads), "--iterations", string(iterations)]
        if !isnothing(genome_size)
            push!(cmd_args, "--genome-size", string(genome_size))
        end
        if meta
            push!(cmd_args, "--meta")
        end
        keep_haplotypes && push!(cmd_args, "--keep-haplotypes")
        no_alt_contigs && push!(cmd_args, "--no-alt-contigs")
        if !isnothing(min_overlap)
            push!(cmd_args, "--min-overlap", string(min_overlap))
        end
        cmd = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n flye $(cmd_args)`
        )
        script = join([
            "set -euo pipefail",
            "mkdir -p \"$(outdir)\"",
            "if [ ! -f \"$(joinpath(outdir, "assembly.fasta"))\" ]; then",
            "  $(cmd)",
            "fi"
        ], "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, assembly = joinpath(outdir, "assembly.fasta"),
            graph = joinpath(outdir, "assembly_graph.gfa"))
    end

    if !isfile(joinpath(outdir, "assembly.fasta"))
        cmd_args = ["flye", "--$(read_type)", fastq, "--out-dir",
            outdir, "--threads", string(threads), "--iterations", string(iterations)]
        if !isnothing(genome_size)
            push!(cmd_args, "--genome-size", string(genome_size))
        end

        if meta
            push!(cmd_args, "--meta")
        end
        keep_haplotypes && push!(cmd_args, "--keep-haplotypes")
        no_alt_contigs && push!(cmd_args, "--no-alt-contigs")

        if !isnothing(min_overlap)
            push!(cmd_args, "--min-overlap", string(min_overlap))
        end

        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n flye $(cmd_args)`)
    end
    return (; outdir, assembly = joinpath(outdir, "assembly.fasta"),
        graph = joinpath(outdir, "assembly_graph.gfa"))
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
function run_canu(;
        fastq,
        outdir = "canu_output",
        genome_size,
        read_type = "pacbio",
        stopOnLowCoverage = 10,
        threads::Int = get_default_threads(),
        cor_threads::Int = threads,
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "canu",
        time_limit::String = "3-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("canu")
    mkpath(outdir)
    threads = max(threads, 1)
    cor_threads = clamp(cor_threads, 1, threads)

    prefix = splitext(basename(fastq))[1]
    if executor !== nothing
        cmd = if read_type == "pacbio"
            Mycelia.command_string(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -pacbio $(fastq) maxThreads=$(threads) corThreads=$(cor_threads) stopOnLowCoverage=$(stopOnLowCoverage) saveReads=false useGrid=false`
            )
        elseif read_type == "nanopore"
            Mycelia.command_string(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -nanopore $(fastq) maxThreads=$(threads) corThreads=$(cor_threads) stopOnLowCoverage=$(stopOnLowCoverage) saveReads=false useGrid=false`
            )
        else
            error("Unsupported read type: $(read_type). Use 'pacbio' or 'nanopore'.")
        end
        script = join([
            "set -euo pipefail",
            "mkdir -p \"$(outdir)\"",
            "if [ ! -f \"$(joinpath(outdir, "$(prefix).contigs.fasta"))\" ]; then",
            "  $(cmd)",
            "fi"
        ], "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, assembly = joinpath(outdir, "$(prefix).contigs.fasta"),
            graph = joinpath(outdir, "$(prefix).contigs.gfa"))
    end

    if !isfile(joinpath(outdir, "$(prefix).contigs.fasta"))
        if read_type == "pacbio"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -pacbio $(fastq) maxThreads=$(threads) corThreads=$(cor_threads) stopOnLowCoverage=$(stopOnLowCoverage) saveReads=false useGrid=false`)
        elseif read_type == "nanopore"
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n canu canu -p $(prefix) -d $(outdir) genomeSize=$(genome_size) -nanopore $(fastq) maxThreads=$(threads) corThreads=$(cor_threads) stopOnLowCoverage=$(stopOnLowCoverage) saveReads=false useGrid=false`)
        else
            error("Unsupported read type: $(read_type). Use 'pacbio' or 'nanopore'.")
        end
    end
    return (; outdir, assembly = joinpath(outdir, "$(prefix).contigs.fasta"),
        graph = joinpath(outdir, "$(prefix).contigs.gfa"))
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
function run_hifiasm(;
        fastq,
        outdir = basename(fastq) * "_hifiasm",
        bloom_filter = -1,
        threads::Int = get_default_threads(),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "hifiasm",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("hifiasm")
    hifiasm_outprefix = joinpath(outdir, basename(fastq) * ".hifiasm")
    # Check if output directory exists before trying to read it
    hifiasm_outputs = isdir(outdir) ?
                      filter(x -> occursin(hifiasm_outprefix, x), readdir(outdir, join = true)) :
                      String[]
    if executor !== nothing
        cmd_args = [
            "hifiasm", "--primary", "-l0", "-o", hifiasm_outprefix, "-t", string(threads)]
        if bloom_filter >= 0
            push!(cmd_args, "-f$(bloom_filter)")
        end
        push!(cmd_args, fastq)
        cmd = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n hifiasm $(cmd_args)`
        )
        script = join([
            "set -euo pipefail",
            "if ! compgen -G \"$(hifiasm_outprefix)*\" > /dev/null; then",
            "  mkdir -p \"$(outdir)\"",
            "  $(cmd)",
            "fi"
        ], "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, hifiasm_outprefix)
    end

    # https://hifiasm.readthedocs.io/en/latest/faq.html#are-inbred-homozygous-genomes-supported
    if isempty(hifiasm_outputs)
        mkpath(outdir)
        # Build command with optional bloom filter flag
        cmd_args = [
            "hifiasm", "--primary", "-l0", "-o", hifiasm_outprefix, "-t", string(threads)]
        if bloom_filter >= 0
            push!(cmd_args, "-f$(bloom_filter)")
        end
        push!(cmd_args, fastq)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hifiasm $(cmd_args)`)
    end
    return (; outdir, hifiasm_outprefix)
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
function run_hifiasm_meta(;
        fastq,
        outdir = basename(fastq) * "_hifiasm_meta",
        bloom_filter = -1,
        read_selection = false,
        threads::Int = get_default_threads(),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "hifiasm_meta",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("hifiasm_meta")
    hifiasm_outprefix = joinpath(outdir, basename(fastq) * ".hifiasm_meta")
    # Check if output directory exists before trying to read it
    hifiasm_outputs = isdir(outdir) ?
                      filter(x -> occursin(hifiasm_outprefix, x), readdir(outdir, join = true)) :
                      String[]

    if executor !== nothing
        cmd_args = ["hifiasm_meta", "-t", string(threads), "-o", hifiasm_outprefix]
        if bloom_filter >= 0
            push!(cmd_args, "-f$(bloom_filter)")
        end
        if read_selection
            push!(cmd_args, "-S")
        end
        push!(cmd_args, fastq)
        cmd = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n hifiasm_meta $(cmd_args)`
        )
        script = join([
            "set -euo pipefail",
            "if ! compgen -G \"$(hifiasm_outprefix)*\" > /dev/null; then",
            "  mkdir -p \"$(outdir)\"",
            "  $(cmd)",
            "fi"
        ], "\n")
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, hifiasm_outprefix)
    end

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
    return (; outdir, hifiasm_outprefix)
end

const _UNICYCLER_ENVIRONMENT_LOCK_STALE_SECONDS = 7 * 24 * 60 * 60
const _UNICYCLER_ENVIRONMENT_LOCK_REFRESH_SECONDS = 60
const _UNICYCLER_CONTRACT_FILENAME = ".mycelia-unicycler-contract.json"
const _UNICYCLER_CONTRACT_SCHEMA = "mycelia-unicycler-run-contract-v2"
const _UNICYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING = 1_000_000_000_000

function _normalize_unicycler_package_inventory(
        package_records::Any,
)::Vector{NamedTuple}
    package_records isa AbstractVector || error(
        "Unicycler Conda package inventory was not a JSON array.",
    )
    packages = NamedTuple[]
    for (record_index, package_record) in enumerate(package_records)
        if package_record isa AbstractDict
            name = get(package_record, "name", nothing)
            version = get(package_record, "version", nothing)
            build = get(
                package_record,
                "build_string",
                get(package_record, "build", nothing),
            )
            channel = get(package_record, "channel", nothing)
        elseif package_record isa NamedTuple
            name = hasproperty(package_record, :name) ?
                   package_record.name : nothing
            version = hasproperty(package_record, :version) ?
                      package_record.version : nothing
            build = if hasproperty(package_record, :build_string)
                package_record.build_string
            elseif hasproperty(package_record, :build)
                package_record.build
            else
                nothing
            end
            channel = hasproperty(package_record, :channel) ?
                      package_record.channel : nothing
        else
            error(
                "Unicycler Conda package inventory record $(record_index) " *
                "is not an object.",
            )
        end
        all(
            value -> value isa AbstractString && !isempty(value),
            (name, version, build, channel),
        ) || error(
            "Unicycler Conda package inventory record $(record_index) must " *
            "report non-empty name, version, build, and channel fields.",
        )
        push!(packages, (
            name = String(name),
            version = String(version),
            build = String(build),
            channel = String(channel),
        ))
    end
    isempty(packages) && error("Unicycler Conda package inventory is empty.")
    sort!(
        packages;
        by = package -> (
            package.name,
            package.version,
            package.build,
            package.channel,
        ),
    )
    length(unique(package.name for package in packages)) == length(packages) ||
        error("Unicycler Conda package inventory contains duplicate package names.")
    return packages
end

function _unicycler_conda_package_inventory(;
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::AbstractString =
            _unicycler_environment_prefix(conda_runner),
        command_reader::Function = command -> read(command, String),
)::Vector{NamedTuple}
    resolved_runner = _canonical_unicycler_conda_runner(conda_runner)
    resolved_prefix = _canonical_unicycler_environment_prefix(
        environment_prefix,
    )
    command = Cmd(String[
        resolved_runner,
        "list",
        "-p",
        resolved_prefix,
        "--json",
    ])
    package_records = JSON.parse(command_reader(command))
    return _normalize_unicycler_package_inventory(package_records)
end

function _unicycler_package_inventory_sha256(package_records::Any)::String
    packages = _normalize_unicycler_package_inventory(package_records)
    canonical = IOBuffer()
    write(canonical, "mycelia-conda-package-inventory-v1")
    print(canonical, length(packages), ':')
    for package in packages
        for field in (
                package.name,
                package.version,
                package.build,
                package.channel,
        )
            print(canonical, ncodeunits(field), ':')
            write(canonical, field)
        end
    end
    return bytes2hex(SHA.sha256(take!(canonical)))
end

function _unicycler_toolchain_metadata(
        package_records::Any,
)::Dict{String, Any}
    packages = _normalize_unicycler_package_inventory(package_records)
    package_names = Set(package.name for package in packages)
    for required_package in ("unicycler", "spades")
        required_package in package_names || error(
            "Unicycler realized Conda package inventory is missing required " *
            "package $(repr(required_package)).",
        )
    end
    return Dict{String, Any}(
        "environment_name" => "unicycler",
        "inventory_schema" => "conda-name-version-build-channel-v1",
        "package_inventory_sha256" =>
            _unicycler_package_inventory_sha256(packages),
        "packages" => Dict{String, Any}[
            Dict{String, Any}(
                "name" => package.name,
                "version" => package.version,
                "build" => package.build,
                "channel" => package.channel,
            ) for package in packages
        ],
    )
end

function _unicycler_toolchain_provenance(;
        inventory_reader::Function = _unicycler_conda_package_inventory,
)::Dict{String, Any}
    return _unicycler_toolchain_metadata(inventory_reader())
end

function _require_unicycler_toolchain_provenance(
        toolchain::Any,
)::Dict{String, Any}
    toolchain isa AbstractDict || error(
        "Unicycler workflow did not report realized toolchain provenance.",
    )
    normalized = Dict{String, Any}(
        String(key) => value for (key, value) in pairs(toolchain)
    )
    get(normalized, "environment_name", nothing) == "unicycler" || error(
        "Unicycler workflow toolchain provenance has the wrong environment name.",
    )
    get(normalized, "inventory_schema", nothing) ==
        "conda-name-version-build-channel-v1" || error(
        "Unicycler workflow toolchain provenance has an unsupported inventory schema.",
    )
    digest = get(normalized, "package_inventory_sha256", nothing)
    packages = get(normalized, "packages", nothing)
    digest isa AbstractString && occursin(r"^[0-9a-f]{64}$", digest) || error(
        "Unicycler workflow toolchain provenance is missing a valid package " *
        "inventory SHA-256 digest.",
    )
    packages isa AbstractVector && !isempty(packages) || error(
        "Unicycler workflow toolchain provenance has no realized package inventory.",
    )
    normalized_packages = _normalize_unicycler_package_inventory(packages)
    expected = _unicycler_toolchain_metadata(normalized_packages)
    digest == expected["package_inventory_sha256"] || error(
        "Unicycler workflow package inventory digest does not match its " *
        "reported name/version/build/channel inventory.",
    )
    return expected
end

function _canonical_unicycler_conda_runner(
        conda_runner::AbstractString,
)::String
    executable = Sys.which(String(conda_runner))
    candidate = executable === nothing ?
                abspath(String(conda_runner)) : String(executable)
    return ispath(candidate) ? realpath(candidate) : normpath(candidate)
end

function _canonical_unicycler_environment_prefix(
        environment_prefix::AbstractString,
)::String
    normalized_prefix = normpath(abspath(String(environment_prefix)))
    existing_ancestor = normalized_prefix
    missing_components = String[]
    while !ispath(existing_ancestor) && !islink(existing_ancestor)
        parent = dirname(existing_ancestor)
        parent == existing_ancestor && break
        pushfirst!(missing_components, basename(existing_ancestor))
        existing_ancestor = parent
    end
    isdir(existing_ancestor) || throw(ArgumentError(
        "Unicycler environment prefix has a non-directory existing " *
        "ancestor: $(repr(existing_ancestor)).",
    ))
    canonical_ancestor = realpath(existing_ancestor)
    return isempty(missing_components) ? canonical_ancestor :
           normpath(joinpath(canonical_ancestor, missing_components...))
end

function _unicycler_environment_prefix(
        conda_runner::AbstractString,
)::String
    canonical_runner = _canonical_unicycler_conda_runner(conda_runner)
    existing_prefix = _conda_env_prefix_from_runner(
        "unicycler",
        canonical_runner,
    )
    candidate = if existing_prefix === nothing
        conda_root = normpath(joinpath(dirname(canonical_runner), ".."))
        joinpath(conda_root, "envs", "unicycler")
    else
        String(existing_prefix)
    end
    return _canonical_unicycler_environment_prefix(candidate)
end

function _unicycler_environment_lock_path(
        conda_runner::AbstractString = _conda_runner(),
)::String
    environment_prefix = _unicycler_environment_prefix(conda_runner)
    return _unicycler_environment_lock_path_from_prefix(environment_prefix)
end

function _unicycler_environment_lock_path_from_prefix(
        environment_prefix::AbstractString,
)::String
    normalized_prefix = _canonical_unicycler_environment_prefix(
        environment_prefix,
    )
    return joinpath(
        dirname(normalized_prefix),
        ".$(basename(normalized_prefix)).mycelia-unicycler.pid",
    )
end

function _with_unicycler_environment_lock(
        action::Function,
        lock_path::AbstractString = _unicycler_environment_lock_path();
        pidlock_runner::Function = FileWatching.Pidfile.mkpidlock,
)::Any
    normalized_lock_path = abspath(String(lock_path))
    mkpath(dirname(normalized_lock_path))
    return pidlock_runner(
        normalized_lock_path;
        stale_age = _UNICYCLER_ENVIRONMENT_LOCK_STALE_SECONDS,
        refresh = _UNICYCLER_ENVIRONMENT_LOCK_REFRESH_SECONDS,
    ) do
        action()
    end
end

function _nearest_existing_unicycler_ancestor(path::AbstractString)::String
    ancestor = abspath(String(path))
    while !ispath(ancestor) && !islink(ancestor)
        parent = dirname(ancestor)
        parent == ancestor && return ancestor
        ancestor = parent
    end
    return ancestor
end

function _canonical_planned_unicycler_output_path(
        outdir::AbstractString,
)::String
    normalized_outdir = abspath(String(outdir))
    islink(normalized_outdir) && throw(ArgumentError(
        "Unicycler outdir must not be a symbolic link: " *
        repr(normalized_outdir),
    ))
    ancestor = _nearest_existing_unicycler_ancestor(normalized_outdir)
    isdir(ancestor) || throw(ArgumentError(
        "Unicycler outdir has a non-directory existing ancestor: " *
        repr(ancestor),
    ))
    canonical_ancestor = realpath(ancestor)
    relative_path = relpath(normalized_outdir, ancestor)
    return relative_path == "." ? canonical_ancestor :
           normpath(joinpath(canonical_ancestor, relative_path))
end

function _unicycler_output_lock_path(outdir::AbstractString)::String
    canonical_outdir = _canonical_planned_unicycler_output_path(outdir)
    return _output_root_reservation_lock_path_from_canonical(canonical_outdir)
end

function _with_unicycler_output_lock(
        action::Function,
        outdir::AbstractString;
        pidlock_runner::Function = FileWatching.Pidfile.trymkpidlock,
)::Any
    canonical_outdir = _canonical_planned_unicycler_output_path(outdir)
    lock_path = _unicycler_output_lock_path(canonical_outdir)
    _require_exclusive_output_root_reservation(
        canonical_outdir,
        lock_path;
        subject = "Unicycler outdir",
        stale_age = _UNICYCLER_ENVIRONMENT_LOCK_STALE_SECONDS,
    )
    mkpath(dirname(lock_path))
    locked_action = function ()
        _require_exclusive_output_root_reservation(
            canonical_outdir,
            lock_path;
            subject = "Unicycler outdir",
            stale_age = _UNICYCLER_ENVIRONMENT_LOCK_STALE_SECONDS,
        )
        return action(canonical_outdir)
    end
    if pidlock_runner !== FileWatching.Pidfile.trymkpidlock
        return pidlock_runner(
            locked_action,
            lock_path;
            stale_age = _UNICYCLER_ENVIRONMENT_LOCK_STALE_SECONDS,
            refresh = _UNICYCLER_ENVIRONMENT_LOCK_REFRESH_SECONDS,
        )
    end
    lock_handle = pidlock_runner(
        lock_path;
        stale_age = _UNICYCLER_ENVIRONMENT_LOCK_STALE_SECONDS,
        refresh = _UNICYCLER_ENVIRONMENT_LOCK_REFRESH_SECONDS,
    )
    lock_handle === false && throw(ArgumentError(
        "Unicycler outdir is already reserved by another output-root " *
        "workflow: $(lock_path)",
    ))
    try
        return locked_action()
    finally
        Base.close(lock_handle)
    end
end

function _unicycler_output_adjacent_spool_parent(
        outdir::AbstractString,
)::String
    parent = dirname(normpath(abspath(String(outdir))))
    while !ispath(parent) && !islink(parent)
        next_parent = dirname(parent)
        next_parent == parent && break
        parent = next_parent
    end
    isdir(parent) && !islink(parent) || throw(ArgumentError(
        "Unicycler output-adjacent input spool parent must resolve to an " *
        "existing non-symlink directory: $(parent).",
    ))
    return realpath(parent)
end

function _require_unchanged_unicycler_output_plan(
        requested_outdir::AbstractString,
        planned_outdir::AbstractString,
)::String
    normalized_plan = normpath(abspath(String(planned_outdir)))
    observed_plan = _canonical_planned_unicycler_output_path(requested_outdir)
    observed_plan == normalized_plan || throw(ArgumentError(
        "Unicycler outdir changed physical identity after reservation: " *
        "planned $(repr(normalized_plan)), observed $(repr(observed_plan)).",
    ))
    return normalized_plan
end

function _require_unicycler_pre_environment_output_shape(
        outdir::AbstractString,
)::Nothing
    normalized_outdir = normpath(abspath(String(outdir)))
    if islink(normalized_outdir) ||
       (ispath(normalized_outdir) && !isdir(normalized_outdir))
        throw(ArgumentError(
            "Unicycler outdir must be absent or an empty directory: " *
            repr(normalized_outdir),
        ))
    end
    isdir(normalized_outdir) || return nothing
    entries = readdir(normalized_outdir)
    isempty(entries) && return nothing
    required = (
        joinpath(normalized_outdir, "assembly.fasta"),
        joinpath(normalized_outdir, "assembly.gfa"),
        joinpath(normalized_outdir, _UNICYCLER_CONTRACT_FILENAME),
    )
    all(path -> isfile(path) && !islink(path), required) ||
        throw(ArgumentError(
            "Refusing incomplete non-empty Unicycler outdir before " *
            "environment preparation: $(repr(normalized_outdir)).",
        ))
    return nothing
end

function _require_exact_unicycler_output_dir(
        outdir::AbstractString,
)::String
    normalized_outdir = normpath(abspath(String(outdir)))
    if islink(normalized_outdir) || !isdir(normalized_outdir)
        error(
            "Unicycler did not preserve its exact reserved output directory: " *
            repr(normalized_outdir),
        )
    end
    canonical_outdir = realpath(normalized_outdir)
    canonical_outdir == normalized_outdir || error(
        "Unicycler reserved output directory resolves through a symbolic-link " *
        "component: $(repr(normalized_outdir)) resolves to " *
        "$(repr(canonical_outdir)).",
    )
    return normalized_outdir
end

function _unicycler_command(;
        conda_runner::AbstractString,
        environment_prefix::AbstractString,
        short_1::AbstractString,
        short_2::Union{Nothing, AbstractString},
        long_reads::AbstractString,
        outdir::AbstractString,
        threads::Int,
        spades_options::Union{Nothing, String},
        kmers::Union{Nothing, String},
)::Cmd
    resolved_runner = _canonical_unicycler_conda_runner(conda_runner)
    resolved_prefix = _canonical_unicycler_environment_prefix(
        environment_prefix,
    )
    arguments = String[
        resolved_runner,
        "run",
        "--live-stream",
        "-p",
        resolved_prefix,
        "unicycler",
    ]
    if isnothing(short_2)
        append!(arguments, ["-s", String(short_1)])
    else
        append!(arguments, ["-1", String(short_1), "-2", String(short_2)])
    end
    append!(arguments, [
        "-l",
        String(long_reads),
        "-o",
        String(outdir),
        "-t",
        string(threads),
    ])
    isnothing(kmers) || push!(arguments, "--kmers=$(kmers)")
    isnothing(spades_options) ||
        push!(arguments, "--spades_options=$(spades_options)")
    return Cmd(arguments)
end

function _prepare_unicycler_environment(
        conda_runner::AbstractString,
        environment_prefix::AbstractString,
)::Nothing
    resolved_runner = _canonical_unicycler_conda_runner(conda_runner)
    resolved_prefix = _canonical_unicycler_environment_prefix(
        environment_prefix,
    )
    if !isdir(joinpath(resolved_prefix, "conda-meta"))
        run(Cmd([
            resolved_runner,
            "create",
            "-c",
            "conda-forge",
            "-c",
            "bioconda",
            "-c",
            "defaults",
            "--strict-channel-priority",
            "-p",
            resolved_prefix,
            "unicycler",
            "-y",
        ]))
        run(Cmd([resolved_runner, "clean", "--all", "-y"]))
    end
    return nothing
end

function _require_unicycler_artifact(
        path::AbstractString,
        label::AbstractString,
        outdir::AbstractString,
)::String
    normalized_outdir = _require_exact_unicycler_output_dir(outdir)
    normalized_path = normpath(abspath(String(path)))
    islink(normalized_path) && error(
        "Unicycler $(label) must be a regular non-symlink file: " *
        repr(normalized_path),
    )
    isfile(normalized_path) && filesize(normalized_path) > 0 || error(
        "Unicycler did not produce a non-empty $(label) at " *
        repr(normalized_path),
    )
    canonical_path = realpath(normalized_path)
    canonical_path == normalized_path || error(
        "Unicycler $(label) resolves through a symbolic-link component: " *
        "$(repr(normalized_path)) resolves to $(repr(canonical_path)).",
    )
    relative_path = relpath(canonical_path, normalized_outdir)
    relative_parts = splitpath(relative_path)
    if relative_path == "." ||
       isabspath(relative_path) ||
       (!isempty(relative_parts) && first(relative_parts) == "..")
        error(
            "Unicycler $(label) escapes its exact reserved output directory " *
            "$(repr(normalized_outdir)): $(repr(normalized_path)).",
        )
    end
    return normalized_path
end

function _unicycler_sha256(path::AbstractString)::String
    return open(path, "r") do input
        SHA.bytes2hex(SHA.sha256(input))
    end
end

function _unicycler_input_fingerprint(
        path::AbstractString,
        label::AbstractString,
)::Dict{String, Any}
    normalized_path = normpath(abspath(String(path)))
    isfile(normalized_path) || throw(ArgumentError(
        "Unicycler $(label) is not a regular file: $(repr(normalized_path)).",
    ))
    canonical_path = realpath(normalized_path)
    isfile(canonical_path) && filesize(canonical_path) > 0 || throw(
        ArgumentError(
            "Unicycler $(label) must be a non-empty regular file: " *
            repr(canonical_path),
        ),
    )
    return Dict{String, Any}(
        "canonical_path" => canonical_path,
        "size_bytes" => filesize(canonical_path),
        "sha256" => _unicycler_sha256(canonical_path),
    )
end

function _unicycler_input_fingerprints(
        short_1::AbstractString,
        short_2::Union{Nothing, AbstractString},
        long_reads::AbstractString;
        fingerprinter::Function = _unicycler_input_fingerprint,
)::Dict{String, Any}
    return Dict{String, Any}(
        "short_1" => fingerprinter(short_1, "short-read input R1"),
        "short_2" => short_2 === nothing ? nothing :
                     fingerprinter(short_2, "short-read input R2"),
        "long_reads" => fingerprinter(long_reads, "long-read input"),
    )
end

function _unicycler_pair_identifier(identifier::AbstractString)::String
    first_token = first(split(String(identifier)))
    return replace(first_token, r"/[12]$" => "")
end

function _unicycler_identifier_pair_role(
        identifier::AbstractString,
)::Union{Nothing, Int}
    first_token = first(split(String(identifier)))
    role_match = match(r"/([12])$", first_token)
    return role_match === nothing ? nothing :
           parse(Int, something(only(role_match.captures)))
end

function _unicycler_casava_pair_role(
        description::AbstractString,
)::Union{Nothing, Int}
    description_tokens = split(String(description))
    length(description_tokens) >= 2 || return nothing
    role_match = match(
        r"^([12]):[YN]:[0-9]+:[A-Za-z0-9+_-]+$",
        description_tokens[2],
    )
    return role_match === nothing ? nothing :
           parse(Int, something(only(role_match.captures)))
end

function _unicycler_pair_role(
        identifier::AbstractString,
        description::AbstractString,
)::Union{Nothing, Int}
    identifier_role = _unicycler_identifier_pair_role(identifier)
    casava_role = _unicycler_casava_pair_role(description)
    if identifier_role !== nothing && casava_role !== nothing &&
       identifier_role != casava_role
        throw(ArgumentError(
            "Unicycler FASTQ identifier and CASAVA description contain " *
            "conflicting explicit mate roles: " *
            "identifier=$(repr(String(identifier))), " *
            "description=$(repr(String(description))).",
        ))
    end
    return identifier_role === nothing ? casava_role : identifier_role
end

function _require_unicycler_fastq_path(
        path::AbstractString,
        label::AbstractString,
)::String
    normalized_path = normpath(abspath(String(path)))
    isfile(normalized_path) || throw(ArgumentError(
        "Unicycler $(label) is not a regular file: $(repr(normalized_path)).",
    ))
    filesize(normalized_path) > 0 || throw(ArgumentError(
        "Unicycler $(label) must be non-empty: $(repr(normalized_path)).",
    ))
    return normalized_path
end

function _open_unicycler_fastq(
        path::AbstractString,
        label::AbstractString,
)::Any
    normalized_path = _require_unicycler_fastq_path(path, label)
    return try
        Mycelia.open_fastx(normalized_path)
    catch caught
        caught isa InterruptException && rethrow()
        throw(ArgumentError(
            "Unicycler $(label) is not valid FASTQ: " *
            "$(repr(normalized_path)). Cause: $(sprint(showerror, caught))",
        ))
    end
end

function _validate_unicycler_fastq(
        path::AbstractString,
        label::AbstractString,
)::Int
    reader = _open_unicycler_fastq(path, label)
    record_count = 0
    try
        for record in reader
            record isa FASTX.FASTQ.Record || throw(ArgumentError(
                "Unicycler $(label) must be a FASTQ file: " *
                repr(normpath(abspath(String(path)))),
            ))
            record_count += 1
        end
    catch caught
        caught isa InterruptException && rethrow()
        caught isa ArgumentError && rethrow()
        throw(ArgumentError(
            "Unicycler $(label) is not valid FASTQ: " *
            "$(repr(normpath(abspath(String(path))))). " *
            "Cause: $(sprint(showerror, caught))",
        ))
    finally
        close(reader)
    end
    record_count > 0 || throw(ArgumentError(
        "Unicycler $(label) must contain at least one FASTQ record.",
    ))
    return record_count
end

function _validate_unicycler_paired_fastqs(
        short_1::AbstractString,
        short_2::AbstractString,
)::Int
    reader_1 = _open_unicycler_fastq(short_1, "short-read input R1")
    reader_2 = try
        _open_unicycler_fastq(short_2, "short-read input R2")
    catch
        close(reader_1)
        rethrow()
    end
    pair_count = 0
    try
        next_1 = iterate(reader_1)
        next_2 = iterate(reader_2)
        while next_1 !== nothing || next_2 !== nothing
            if next_1 === nothing || next_2 === nothing
                throw(ArgumentError(
                    "Unicycler paired short reads have different counts " *
                    "after $(pair_count) complete pairs.",
                ))
            end
            record_1, state_1 = next_1
            record_2, state_2 = next_2
            pair_count += 1
            if !(record_1 isa FASTX.FASTQ.Record) ||
               !(record_2 isa FASTX.FASTQ.Record)
                throw(ArgumentError(
                    "Unicycler paired short-read inputs must be FASTQ files.",
                ))
            end
            identifier_1 = String(FASTX.identifier(record_1))
            identifier_2 = String(FASTX.identifier(record_2))
            role_1 = _unicycler_pair_role(
                identifier_1,
                String(FASTX.description(record_1)),
            )
            role_2 = _unicycler_pair_role(
                identifier_2,
                String(FASTX.description(record_2)),
            )
            roles_valid = (role_1 === nothing && role_2 === nothing) ||
                          (role_1 == 1 && role_2 == 2)
            roles_valid || throw(ArgumentError(
                "Unicycler paired short reads have invalid explicit mate " *
                "roles at record $(pair_count): " *
                "R1=$(repr(identifier_1)), R2=$(repr(identifier_2)); " *
                "expected R1 role 1 followed by R2 role 2.",
            ))
            _unicycler_pair_identifier(identifier_1) ==
                _unicycler_pair_identifier(identifier_2) ||
                throw(ArgumentError(
                    "Unicycler paired short reads are out of sync at record " *
                    "$(pair_count): R1=$(repr(identifier_1)), " *
                    "R2=$(repr(identifier_2)).",
                ))
            next_1 = iterate(reader_1, state_1)
            next_2 = iterate(reader_2, state_2)
        end
    catch caught
        caught isa InterruptException && rethrow()
        caught isa ArgumentError && rethrow()
        throw(ArgumentError(
            "Unicycler paired short-read inputs are not valid FASTQ. " *
            "Cause: $(sprint(showerror, caught))",
        ))
    finally
        try
            close(reader_1)
        finally
            close(reader_2)
        end
    end
    pair_count > 0 || throw(ArgumentError(
        "Unicycler paired short reads must be non-empty.",
    ))
    return pair_count
end

function _validate_unicycler_inputs(
        short_1::AbstractString,
        short_2::Union{Nothing, AbstractString},
        long_reads::AbstractString,
)::NamedTuple
    normalized_short_1 = _require_unicycler_fastq_path(
        short_1,
        "short-read input R1",
    )
    normalized_long_reads = _require_unicycler_fastq_path(
        long_reads,
        "long-read input",
    )
    normalized_short_2 = short_2 === nothing ? nothing :
                         _require_unicycler_fastq_path(
        short_2,
        "short-read input R2",
    )
    physical_paths = Pair{String, String}[
        "short-read input R1" => normalized_short_1,
        "long-read input" => normalized_long_reads,
    ]
    normalized_short_2 === nothing || push!(
        physical_paths,
        "short-read input R2" => normalized_short_2,
    )
    for first_index in eachindex(physical_paths)
        for second_index in (first_index + 1):lastindex(physical_paths)
            first_label, first_path = physical_paths[first_index]
            second_label, second_path = physical_paths[second_index]
            Base.Filesystem.samefile(first_path, second_path) && throw(
                ArgumentError(
                    "Unicycler inputs must be physically distinct; " *
                    "$(first_label) and $(second_label) refer to the same " *
                    "file.",
                ),
            )
        end
    end
    short_count = if normalized_short_2 === nothing
        _validate_unicycler_fastq(
            normalized_short_1,
            "single short-read input",
        )
    else
        _validate_unicycler_paired_fastqs(
            normalized_short_1,
            normalized_short_2,
        )
    end
    long_count = _validate_unicycler_fastq(
        normalized_long_reads,
        "long-read input",
    )
    return (;
        short_1 = normalized_short_1,
        short_2 = normalized_short_2,
        long_reads = normalized_long_reads,
        short_count,
        long_count,
    )
end

function _unicycler_input_source_snapshot(
        path::AbstractString,
        label::AbstractString,
)::NamedTuple
    normalized_path = _require_unicycler_fastq_path(path, label)
    canonical_path = realpath(normalized_path)
    source_status = stat(canonical_path)
    isfile(source_status) || throw(ArgumentError(
        "Unicycler $(label) is not a regular file: $(normalized_path).",
    ))
    return (;
        path = normalized_path,
        canonical_path,
        size_bytes = filesize(canonical_path),
        device = source_status.device,
        inode = source_status.inode,
    )
end

function _require_unchanged_unicycler_input_source(
        source::NamedTuple,
        label::AbstractString,
)::Nothing
    observed = _unicycler_input_source_snapshot(source.path, label)
    observed == source || throw(ErrorException(
        "Unicycler $(label) changed physical identity or size after spool " *
        "preflight.",
    ))
    return nothing
end

function _unicycler_spool_root_identity(
        spool_root::AbstractString,
)::NamedTuple
    normalized_root = normpath(abspath(String(spool_root)))
    isdir(normalized_root) && !islink(normalized_root) || error(
        "Unicycler private input spool root is missing, replaced, or not a " *
        "directory: $(normalized_root).",
    )
    root_status = stat(normalized_root)
    return (;
        path = normalized_root,
        device = root_status.device,
        inode = root_status.inode,
    )
end

function _require_unchanged_unicycler_spool_root(
        expected::NamedTuple,
)::Nothing
    observed = _unicycler_spool_root_identity(expected.path)
    observed == expected || error(
        "Unicycler private input spool root changed physical identity.",
    )
    return nothing
end

function _remove_exact_private_input_spool_root!(
        expected::NamedTuple,
        label::AbstractString,
)::Nothing
    parent_path = dirname(expected.path)
    parent = Base.Filesystem.open(
        parent_path,
        Base.JL_O_RDONLY |
        Base.JL_O_DIRECTORY |
        Base.JL_O_NOFOLLOW |
        Base.JL_O_CLOEXEC,
    )
    root = nothing
    try
        root = _metamdbg_openat(
            parent,
            basename(expected.path),
            Base.JL_O_RDONLY |
            Base.JL_O_DIRECTORY |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
        )
        root_status = stat(root)
        root_status.device == expected.device &&
            root_status.inode == expected.inode || error(
                "$(label) private input spool root changed before cleanup.",
            )
        for entry in _metamdbg_descriptor_directory_entries(root)
            _metamdbg_unlinkat(root, entry)
        end
        isempty(_metamdbg_descriptor_directory_entries(root)) || error(
            "$(label) private input spool root retained unexpected entries.",
        )
        observed = _metamdbg_openat_entry_identity(
            parent,
            basename(expected.path),
            expected.path,
            true,
        )
        observed.device == expected.device &&
            observed.inode == expected.inode || error(
                "$(label) private input spool root changed during cleanup.",
            )
        close(root)
        root = nothing
        _metamdbg_unlinkat(
            parent,
            basename(expected.path);
            directory_entry = true,
        )
    finally
        root === nothing || close(root)
        close(parent)
    end
    return nothing
end

function _unicycler_stable_input_snapshot(
        source::NamedTuple,
        label::AbstractString,
        spool_root::AbstractString,
        spool_name::AbstractString;
        after_copy_hook::Function = binding -> nothing,
)::NamedTuple
    _require_unchanged_unicycler_input_source(source, label)
    extension = endswith(lowercase(source.path), ".gz") ? ".fastq.gz" :
                ".fastq"
    consumed_path = joinpath(spool_root, "$(spool_name)$(extension)")
    input = Base.Filesystem.open(
        source.canonical_path,
        Base.JL_O_RDONLY |
        Base.JL_O_NOFOLLOW |
        Base.JL_O_CLOEXEC,
    )
    output = nothing
    try
        input_status = stat(input)
        input_status.device == source.device &&
            input_status.inode == source.inode &&
            filesize(input) == source.size_bytes || throw(ErrorException(
                "Unicycler $(label) changed physical identity or size before " *
                "its stable snapshot was copied.",
            ))
        output = Base.Filesystem.open(
            consumed_path,
            Base.JL_O_WRONLY |
            Base.JL_O_CREAT |
            Base.JL_O_EXCL |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
            0o600,
        )
        remaining = source.size_bytes
        buffer = Vector{UInt8}(undef, 1024 * 1024)
        while remaining > 0
            requested = min(length(buffer), remaining)
            count = readbytes!(input, buffer, requested)
            count > 0 || throw(ErrorException(
                "Unicycler $(label) shrank while its stable input " *
                "snapshot was materialized.",
            ))
            write(output, @view buffer[1:count])
            remaining -= count
        end
        eof(input) || throw(ErrorException(
            "Unicycler $(label) grew while its stable input snapshot " *
            "was materialized.",
        ))
    finally
        output === nothing || close(output)
        close(input)
    end
    after_copy_hook((;
        source_path = source.path,
        consumed_path,
        label = String(label),
    ))
    _require_unchanged_unicycler_input_source(source, label)
    isfile(consumed_path) && !islink(consumed_path) || error(
        "Unicycler stable $(label) snapshot was replaced after copying.",
    )
    consumed_status = stat(consumed_path)
    consumed_status.size == source.size_bytes || error(
        "Unicycler stable $(label) snapshot has an unexpected size after " *
        "copying.",
    )
    chmod(consumed_path, 0o400)
    consumed_sha256 = _unicycler_sha256(consumed_path)
    provenance = Dict{String, Any}(
        "canonical_path" => source.canonical_path,
        "size_bytes" => source.size_bytes,
        "sha256" => consumed_sha256,
        "consumed_snapshot" => Dict{String, Any}(
            "size_bytes" => source.size_bytes,
            "sha256" => consumed_sha256,
        ),
    )
    return (; path = consumed_path, provenance)
end

function _with_unicycler_stable_input_snapshots(
        action::Function,
        short_1::AbstractString,
        short_2::Union{Nothing, AbstractString},
        long_reads::AbstractString;
        spool_parent::Union{Nothing, AbstractString} = nothing,
        byte_ceiling::Integer =
            _UNICYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
        available_bytes_reader::Function = parent -> diskstat(parent).available,
        after_copy_hook::Function = binding -> nothing,
)::Any
    byte_ceiling > 0 || throw(ArgumentError(
        "Unicycler input spool byte ceiling must be positive, got " *
        "$(byte_ceiling).",
    ))
    source_paths = short_2 === nothing ?
                   (short_1, long_reads) :
                   (short_1, short_2, long_reads)
    source_labels = short_2 === nothing ?
                    ("short-read input R1", "long-read input") :
                    (
        "short-read input R1",
        "short-read input R2",
        "long-read input",
    )
    sources = Tuple(
        _unicycler_input_source_snapshot(path, label) for
        (path, label) in zip(source_paths, source_labels)
    )
    total_bytes = sum(
        BigInt(source.size_bytes) for source in sources;
        init = BigInt(0),
    )
    total_bytes <= BigInt(byte_ceiling) || throw(ArgumentError(
        "Unicycler inputs require $(total_bytes) spool bytes, exceeding the " *
        "configured cumulative ceiling of $(byte_ceiling) bytes.",
    ))
    resolved_parent = spool_parent === nothing ? tempdir() :
                      normpath(abspath(String(spool_parent)))
    isdir(resolved_parent) && !islink(resolved_parent) || throw(ArgumentError(
        "Unicycler input spool parent must be an existing non-symlink " *
        "directory: $(resolved_parent).",
    ))
    available_bytes = BigInt(available_bytes_reader(resolved_parent))
    available_bytes >= total_bytes || throw(ArgumentError(
        "Unicycler input spool preflight requires $(total_bytes) bytes but " *
        "only $(available_bytes) bytes are available under " *
        "$(resolved_parent).",
    ))
    spool_root = mktempdir(
        resolved_parent;
        prefix = "mycelia-unicycler-input-",
    )
    chmod(spool_root, 0o700)
    spool_root_identity = _unicycler_spool_root_identity(spool_root)
    try
        _require_unchanged_unicycler_spool_root(spool_root_identity)
        short_1_binding = _unicycler_stable_input_snapshot(
            sources[1],
            "short-read input R1",
            spool_root,
            "short_1";
            after_copy_hook,
        )
        copied_bytes = BigInt(sources[1].size_bytes)
        copied_bytes <= BigInt(byte_ceiling) &&
            copied_bytes <= available_bytes || error(
                "Unicycler stable input copy exceeded its bound budget.",
            )
        _require_unchanged_unicycler_spool_root(spool_root_identity)
        short_2_binding = short_2 === nothing ? nothing :
                          _unicycler_stable_input_snapshot(
            sources[2],
            "short-read input R2",
            spool_root,
            "short_2";
            after_copy_hook,
        )
        short_2 === nothing ||
            (copied_bytes += BigInt(sources[2].size_bytes))
        copied_bytes <= BigInt(byte_ceiling) &&
            copied_bytes <= available_bytes || error(
                "Unicycler stable input copies exceeded their bound budget.",
            )
        _require_unchanged_unicycler_spool_root(spool_root_identity)
        long_binding = _unicycler_stable_input_snapshot(
            sources[end],
            "long-read input",
            spool_root,
            "long_reads";
            after_copy_hook,
        )
        copied_bytes += BigInt(sources[end].size_bytes)
        copied_bytes == total_bytes || error(
            "Unicycler stable input copies did not match their exact " *
            "preflight byte budget.",
        )
        _require_unchanged_unicycler_spool_root(spool_root_identity)
        input_fingerprints = Dict{String, Any}(
            "short_1" => short_1_binding.provenance,
            "short_2" => short_2_binding === nothing ? nothing :
                         short_2_binding.provenance,
            "long_reads" => long_binding.provenance,
        )
        return action((;
            short_1 = short_1_binding.path,
            short_2 = short_2_binding === nothing ? nothing :
                      short_2_binding.path,
            long_reads = long_binding.path,
            input_fingerprints,
            spool_root,
        ))
    catch caught
        caught isa InterruptException && rethrow()
        if caught isa Base.IOError || caught isa SystemError
            throw(ErrorException(
                "Unicycler input spooling failed after cleanup, possibly " *
                "because scratch space was exhausted: " *
                sprint(showerror, caught),
            ))
        end
        rethrow()
    finally
        if ispath(spool_root) || islink(spool_root)
            _remove_exact_private_input_spool_root!(
                spool_root_identity,
                "Unicycler",
            )
        end
    end
end

function _validate_unicycler_stable_input_snapshots(
        snapshots::NamedTuple,
        expected_short_count::Int,
        expected_long_count::Int,
)::NamedTuple
    short_count = if snapshots.short_2 === nothing
        _validate_unicycler_fastq(
            snapshots.short_1,
            "stable single short-read input snapshot",
        )
    else
        _validate_unicycler_paired_fastqs(
            snapshots.short_1,
            snapshots.short_2,
        )
    end
    short_count == expected_short_count || throw(ArgumentError(
        "Unicycler stable short-read input snapshots contain $(short_count) " *
        "records, but semantic preflight validated $(expected_short_count).",
    ))
    long_count = _validate_unicycler_fastq(
        snapshots.long_reads,
        "stable long-read input snapshot",
    )
    long_count == expected_long_count || throw(ArgumentError(
        "Unicycler stable long-read input snapshot contains $(long_count) " *
        "records, but semantic preflight validated $(expected_long_count).",
    ))
    snapshots.input_fingerprints["short_1"]["consumed_snapshot"][
        "record_count"
    ] = short_count
    if snapshots.short_2 !== nothing
        snapshots.input_fingerprints["short_2"]["consumed_snapshot"][
            "record_count"
        ] = short_count
    end
    snapshots.input_fingerprints["long_reads"]["consumed_snapshot"][
        "record_count"
    ] = long_count
    _require_unchanged_unicycler_stable_input_snapshots(snapshots)
    return snapshots
end

function _require_unchanged_unicycler_stable_input_snapshots(
        snapshots::NamedTuple,
)::Nothing
    paths = Dict{String, Union{Nothing, String}}(
        "short_1" => String(snapshots.short_1),
        "short_2" => snapshots.short_2 === nothing ? nothing :
                     String(snapshots.short_2),
        "long_reads" => String(snapshots.long_reads),
    )
    for label in ("short_1", "short_2", "long_reads")
        expected = snapshots.input_fingerprints[label]
        path = paths[label]
        if expected === nothing
            path === nothing || error(
                "Unicycler stable $(label) snapshot appeared unexpectedly.",
            )
            continue
        end
        path === nothing && error(
            "Unicycler stable $(label) snapshot is unexpectedly missing.",
        )
        isfile(path) && !islink(path) || error(
            "Unicycler stable $(label) snapshot is missing or is not a " *
            "regular non-symlink file: $(path).",
        )
        consumed = expected["consumed_snapshot"]
        observed_size = filesize(path)
        observed_sha256 = _unicycler_sha256(path)
        observed_size == consumed["size_bytes"] &&
            observed_sha256 == consumed["sha256"] || error(
                "Unicycler stable $(label) snapshot changed after " *
                "materialization and before environment preparation.",
            )
    end
    return nothing
end

function _require_unicycler_sources_match_stable_snapshots(
        snapshots::NamedTuple,
)::Nothing
    for label in ("short_1", "short_2", "long_reads")
        expected = snapshots.input_fingerprints[label]
        expected === nothing && continue
        observed = _unicycler_input_fingerprint(
            String(expected["canonical_path"]),
            "source $(label)",
        )
        observed["canonical_path"] == expected["canonical_path"] &&
            observed["size_bytes"] == expected["size_bytes"] &&
            observed["sha256"] == expected["sha256"] || error(
                "Unicycler source $(label) changed after its stable " *
                "consumed snapshot was materialized.",
            )
    end
    return nothing
end

function _unicycler_observed_input_fingerprints(
        short_1::AbstractString,
        short_2::Union{Nothing, AbstractString},
        long_reads::AbstractString,
        consumed_short_1::AbstractString,
        consumed_short_2::Union{Nothing, AbstractString},
        consumed_long_reads::AbstractString,
        consumed_input_fingerprints::Union{Nothing, AbstractDict},
        fingerprinter::Function,
)::Dict{String, Any}
    observed = _unicycler_input_fingerprints(
        short_1,
        short_2,
        long_reads;
        fingerprinter,
    )
    consumed_input_fingerprints === nothing && return observed
    consumed_paths = Dict{String, Union{Nothing, String}}(
        "short_1" => String(consumed_short_1),
        "short_2" => consumed_short_2 === nothing ? nothing :
                     String(consumed_short_2),
        "long_reads" => String(consumed_long_reads),
    )
    for label in ("short_1", "short_2", "long_reads")
        expected = consumed_input_fingerprints[label]
        if expected === nothing
            observed[label] === nothing || error(
                "Unicycler $(label) input contract changed unexpectedly.",
            )
            consumed_paths[label] === nothing || error(
                "Unicycler $(label) consumed snapshot appeared unexpectedly.",
            )
            continue
        end
        observed_identity = observed[label]
        observed_identity["size_bytes"] == expected["size_bytes"] &&
            observed_identity["sha256"] == expected["sha256"] || error(
                "Unicycler $(label) source changed after its stable consumed " *
                "snapshot was materialized.",
            )
        expected_consumed = expected["consumed_snapshot"]
        consumed_path = something(consumed_paths[label])
        observed_consumed = _unicycler_input_fingerprint(
            consumed_path,
            "stable consumed $(label) snapshot",
        )
        observed_consumed["size_bytes"] ==
            expected_consumed["size_bytes"] &&
            observed_consumed["sha256"] ==
            expected_consumed["sha256"] || error(
                "Unicycler stable consumed $(label) snapshot changed before " *
                "or during assembler consumption.",
            )
        observed_identity["consumed_snapshot"] = deepcopy(
            expected_consumed,
        )
    end
    return observed
end

function _unicycler_fasta_sequences_from_reader(
        reader::FASTX.FASTA.Reader,
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    records = Dict{String, String}()
    try
        for record in reader
            record isa FASTX.FASTA.Record || error(
                "Unicycler $(label) is not valid FASTA: $(repr(path)).",
            )
            identifier = strip(String(FASTX.identifier(record)))
            isempty(identifier) && error(
                "Unicycler $(label) contains an empty FASTA identifier.",
            )
            haskey(records, identifier) && error(
                "Unicycler $(label) contains duplicate FASTA identifier " *
                repr(identifier),
            )
            sequence = uppercase(String(FASTX.sequence(record)))
            isempty(sequence) && error(
                "Unicycler $(label) contains an empty FASTA sequence for " *
                repr(identifier),
            )
            occursin(r"^[ACGTRYSWKMBDHVN]+$", sequence) || error(
                "Unicycler $(label) contains invalid DNA for FASTA record " *
                repr(identifier),
            )
            records[identifier] = sequence
        end
    catch caught
        caught isa InterruptException && rethrow()
        if caught isa ErrorException && startswith(caught.msg, "Unicycler ")
            rethrow()
        end
        error(
            "Unicycler $(label) is not valid FASTA: $(repr(path)). Cause: " *
            sprint(showerror, caught),
        )
    finally
        close(reader)
    end
    isempty(records) && error(
        "Unicycler $(label) contains no FASTA records: $(repr(path)).",
    )
    return records
end

function _unicycler_fasta_sequences(
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    reader = try
        Mycelia.open_fastx(path)
    catch caught
        caught isa InterruptException && rethrow()
        error(
            "Unicycler $(label) is not valid FASTA: $(repr(path)). Cause: " *
            sprint(showerror, caught),
        )
    end
    return _unicycler_fasta_sequences_from_reader(reader, path, label)
end

function _unicycler_fasta_sequences(
        bytes::AbstractVector{UInt8},
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    reader = FASTX.FASTA.Reader(IOBuffer(bytes))
    return _unicycler_fasta_sequences_from_reader(reader, path, label)
end

function _unicycler_gfa_segment_sequences_from_input(
        input::IO,
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    _require_valid_metamdbg_gfa_input(input, "Unicycler $(label)", path)
    seekstart(input)
    segments = Dict{String, String}()
    for (line_number, line) in enumerate(eachline(input))
        isempty(line) && continue
        startswith(line, "#") && continue
        fields = split(line, '\t'; keepempty = true)
        isempty(fields) && continue
        fields[1] == "S" || continue
        length(fields) >= 3 || error(
            "Unicycler $(label) has a malformed segment at line " *
            "$(line_number): $(repr(path)).",
        )
        identifier = strip(fields[2])
        isempty(identifier) && error(
            "Unicycler $(label) has an empty segment identifier at line " *
            "$(line_number).",
        )
        haskey(segments, identifier) && error(
            "Unicycler $(label) contains duplicate segment identifier " *
            repr(identifier),
        )
        sequence = uppercase(strip(fields[3]))
        (isempty(sequence) || sequence == "*") && error(
            "Unicycler $(label) segment $(repr(identifier)) has no sequence.",
        )
        occursin(r"^[ACGTRYSWKMBDHVN]+$", sequence) || error(
            "Unicycler $(label) segment $(repr(identifier)) contains invalid DNA.",
        )
        segments[identifier] = sequence
    end
    isempty(segments) && error(
        "Unicycler $(label) contains no GFA segment records: $(repr(path)).",
    )
    return segments
end

function _unicycler_gfa_segment_sequences(
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    normalized_path = _require_metamdbg_regular_file(
        path,
        "Unicycler $(label)",
    )
    return open(normalized_path, "r") do input
        _unicycler_gfa_segment_sequences_from_input(
            input,
            normalized_path,
            label,
        )
    end
end

function _unicycler_gfa_segment_sequences(
        bytes::AbstractVector{UInt8},
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    return _unicycler_gfa_segment_sequences_from_input(
        IOBuffer(bytes),
        path,
        label,
    )
end

function _require_valid_unicycler_artifacts(
        assembly_path::AbstractString,
        graph_path::AbstractString,
        outdir::AbstractString,
)::NamedTuple
    assembly = _require_unicycler_artifact(
        assembly_path,
        "assembly FASTA",
        outdir,
    )
    graph = _require_unicycler_artifact(
        graph_path,
        "assembly GFA",
        outdir,
    )
    fasta_sequences = _unicycler_fasta_sequences(
        assembly,
        "assembly FASTA",
    )
    gfa_sequences = _unicycler_gfa_segment_sequences(
        graph,
        "assembly GFA",
    )
    _require_matching_unicycler_artifact_sequences(
        fasta_sequences,
        gfa_sequences,
    )
    return (; assembly, graph)
end

function _require_matching_unicycler_artifact_sequences(
        fasta_sequences::AbstractDict,
        gfa_sequences::AbstractDict,
)::Nothing
    fasta_identifiers = Set(keys(fasta_sequences))
    gfa_identifiers = Set(keys(gfa_sequences))
    fasta_identifiers == gfa_identifiers || error(
        "Unicycler assembly FASTA and GFA contain different contig " *
        "identifiers: FASTA-only=" *
        "$(repr(sort!(collect(setdiff(fasta_identifiers, gfa_identifiers))))), " *
        "GFA-only=" *
        "$(repr(sort!(collect(setdiff(gfa_identifiers, fasta_identifiers))))).",
    )
    mismatched = sort!(String[
        identifier for identifier in fasta_identifiers if
        fasta_sequences[identifier] != gfa_sequences[identifier]
    ])
    isempty(mismatched) || error(
        "Unicycler assembly FASTA and GFA contain different sequences for " *
        "contigs $(repr(mismatched)).",
    )
    return nothing
end

function _unicycler_artifact_fingerprint(path::AbstractString)::Dict{String, Any}
    normalized_path = normpath(abspath(String(path)))
    return Dict{String, Any}(
        "canonical_path" => realpath(normalized_path),
        "size_bytes" => filesize(normalized_path),
        "sha256" => _unicycler_sha256(normalized_path),
    )
end

function _unicycler_buffered_artifact_snapshot(
        path::AbstractString,
        label::AbstractString,
        outdir::AbstractString,
)::NamedTuple
    normalized_path = _require_unicycler_artifact(path, label, outdir)
    bytes = open(normalized_path, "r") do input
        read(input)
    end
    isempty(bytes) && error(
        "Unicycler $(label) became empty while its stable byte snapshot was read.",
    )
    fingerprint = Dict{String, Any}(
        "canonical_path" => normalized_path,
        "size_bytes" => length(bytes),
        "sha256" => SHA.bytes2hex(SHA.sha256(bytes)),
    )
    return (; path = normalized_path, bytes, fingerprint)
end

function _unicycler_semantic_artifact_snapshots(
        assembly_path::AbstractString,
        graph_path::AbstractString,
        outdir::AbstractString;
        after_semantic_validation_hook::Function = artifacts -> nothing,
)::NamedTuple
    assembly = _unicycler_buffered_artifact_snapshot(
        assembly_path,
        "assembly FASTA",
        outdir,
    )
    graph = _unicycler_buffered_artifact_snapshot(
        graph_path,
        "assembly GFA",
        outdir,
    )
    fasta_sequences = _unicycler_fasta_sequences(
        assembly.bytes,
        assembly.path,
        "assembly FASTA",
    )
    gfa_sequences = _unicycler_gfa_segment_sequences(
        graph.bytes,
        graph.path,
        "assembly GFA",
    )
    _require_matching_unicycler_artifact_sequences(
        fasta_sequences,
        gfa_sequences,
    )
    artifacts = (; assembly = assembly.path, graph = graph.path)
    snapshots = (;
        assembly = assembly.fingerprint,
        graph = graph.fingerprint,
    )
    after_semantic_validation_hook(artifacts)
    observed_snapshots = (;
        assembly = _unicycler_artifact_fingerprint(artifacts.assembly),
        graph = _unicycler_artifact_fingerprint(artifacts.graph),
    )
    observed_snapshots == snapshots || error(
        "Unicycler assembly FASTA or GFA changed while its stable semantic " *
        "artifact snapshot was being bound; refusing unvalidated bytes.",
    )
    return (; artifacts..., snapshots)
end

function _require_unchanged_unicycler_artifacts(
        expected_snapshots::NamedTuple,
        assembly_path::AbstractString,
        graph_path::AbstractString,
        outdir::AbstractString,
)::NamedTuple
    observed = _unicycler_semantic_artifact_snapshots(
        assembly_path,
        graph_path,
        outdir,
    )
    observed.snapshots == expected_snapshots || error(
        "Unicycler assembly FASTA or GFA changed after its validated semantic " *
        "artifact snapshot; refusing stale lifecycle provenance.",
    )
    return (; assembly = observed.assembly, graph = observed.graph)
end

function _unicycler_run_contract(
        input_fingerprints::AbstractDict,
        threads::Int,
        spades_options::Union{Nothing, String},
        kmers::Union{Nothing, String},
        environment_prefix::AbstractString,
        toolchain::AbstractDict,
        artifact_snapshots::NamedTuple,
)::Dict{String, Any}
    return Dict{String, Any}(
        "schema" => _UNICYCLER_CONTRACT_SCHEMA,
        "inputs" => deepcopy(input_fingerprints),
        "parameters" => Dict{String, Any}(
            "threads" => threads,
            "spades_options" => spades_options,
            "kmers" => kmers,
            "environment_prefix" => String(environment_prefix),
        ),
        "toolchain" => deepcopy(toolchain),
        "artifacts" => Dict{String, Any}(
            "assembly" => deepcopy(artifact_snapshots.assembly),
            "graph" => deepcopy(artifact_snapshots.graph),
        ),
    )
end

function _write_unicycler_run_contract(
        outdir::AbstractString,
        contract::AbstractDict,
)::String
    normalized_outdir = _require_exact_unicycler_output_dir(outdir)
    marker = joinpath(normalized_outdir, _UNICYCLER_CONTRACT_FILENAME)
    (ispath(marker) || islink(marker)) && error(
        "Refusing to overwrite an existing Unicycler run-contract marker: " *
        repr(marker),
    )
    temporary_path, temporary_io = mktemp(normalized_outdir)
    try
        JSON.print(temporary_io, contract, 2)
        write(temporary_io, '\n')
        flush(temporary_io)
        _fsync_metamdbg_file(temporary_io, temporary_path)
        close(temporary_io)
        mv(temporary_path, marker)
        _fsync_metamdbg_directory(normalized_outdir)
    finally
        isopen(temporary_io) && close(temporary_io)
        (ispath(temporary_path) || islink(temporary_path)) &&
            rm(temporary_path; force = true)
    end
    _require_unicycler_artifact(marker, "run-contract marker", normalized_outdir)
    return marker
end

function _require_matching_unicycler_run_contract(
        outdir::AbstractString,
        expected_contract::AbstractDict,
)::String
    marker = joinpath(
        _require_exact_unicycler_output_dir(outdir),
        _UNICYCLER_CONTRACT_FILENAME,
    )
    (ispath(marker) || islink(marker)) || error(
        "Refusing unbound legacy Unicycler output reuse: missing durable " *
        "run-contract marker $(repr(marker)).",
    )
    _require_unicycler_artifact(marker, "run-contract marker", outdir)
    observed_contract = try
        JSON.parsefile(marker)
    catch caught
        caught isa InterruptException && rethrow()
        error(
            "Unicycler run-contract marker is not valid JSON: " *
            "$(repr(marker)). Cause: $(sprint(showerror, caught))",
        )
    end
    observed_contract == expected_contract || error(
        "Refusing Unicycler output reuse because its durable input, parameter, " *
        "toolchain, or artifact contract does not match this request.",
    )
    return marker
end

function _invoke_unicycler_environment_callback(
        callback::Function,
        conda_runner::AbstractString,
        environment_prefix::AbstractString,
)::Any
    if applicable(callback, conda_runner, environment_prefix)
        return callback(conda_runner, environment_prefix)
    end
    return callback(conda_runner)
end

function _run_unicycler_with_contract(;
        short_1::AbstractString,
        short_2::Union{Nothing, AbstractString} = nothing,
        long_reads::AbstractString,
        consumed_short_1::AbstractString = short_1,
        consumed_short_2::Union{Nothing, AbstractString} = short_2,
        consumed_long_reads::AbstractString = long_reads,
        consumed_input_fingerprints::Union{Nothing, AbstractDict} = nothing,
        outdir::AbstractString = "unicycler_output",
        threads::Int = get_default_threads(),
        spades_options::Union{Nothing, String} = nothing,
        kmers::Union{Nothing, String} = nothing,
        executor::Any = nothing,
        site::Symbol = :local,
        job_name::String = "unicycler",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing,
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        environment_preparer::Function = _prepare_unicycler_environment,
        toolchain_inspector::Function = (runner, prefix) ->
            _unicycler_toolchain_provenance(
                inventory_reader = () ->
                    _unicycler_conda_package_inventory(
                        conda_runner = runner,
                        environment_prefix = prefix,
                    ),
            ),
        input_fingerprinter::Function = _unicycler_input_fingerprint,
        after_artifact_semantic_validation_hook::Function =
            artifacts -> nothing,
        pre_environment_validation::Function = () -> nothing,
        environment_lock_path::Union{Nothing, AbstractString} = nothing,
        environment_lock_runner::Function = _with_unicycler_environment_lock,
        output_lock_runner::Function = _with_unicycler_output_lock,
        command_runner::Function = run,
)::NamedTuple
    reported_outdir = String(outdir)
    isempty(strip(reported_outdir)) && throw(ArgumentError(
        "Unicycler outdir must be a non-empty path.",
    ))
    normalized_outdir = abspath(String(outdir))
    planned_outdir = _canonical_planned_unicycler_output_path(normalized_outdir)
    resolved_conda_runner = _canonical_unicycler_conda_runner(conda_runner)
    resolved_environment_prefix = environment_prefix === nothing ?
                                  _unicycler_environment_prefix(
        resolved_conda_runner,
    ) : _canonical_unicycler_environment_prefix(environment_prefix)
    threads > 0 || throw(ArgumentError(
        "Unicycler threads must be positive, got $(threads).",
    ))
    resolved_executor = resolve_executor(executor)
    assembly = joinpath(planned_outdir, "assembly.fasta")
    graph = joinpath(planned_outdir, "assembly.gfa")
    command = _unicycler_command(;
        conda_runner = resolved_conda_runner,
        environment_prefix = resolved_environment_prefix,
        short_1 = consumed_short_1,
        short_2 = consumed_short_2,
        long_reads = consumed_long_reads,
        outdir = planned_outdir,
        threads,
        spades_options,
        kmers,
    )
    if !(resolved_executor isa LocalExecutor)
        plannable = resolved_executor isa CollectExecutor ||
                    resolved_executor isa DryRunExecutor ||
                    (resolved_executor isa SlurmExecutor &&
                     resolved_executor.dry_run)
        plannable || throw(ArgumentError(
            "Unicycler real nonlocal execution is disabled until the remote " *
            "runtime can hold its output and environment locks. Use " *
            "CollectExecutor, DryRunExecutor, or SlurmExecutor(dry_run=true) " *
            "to obtain an explicit plan.",
        ))
        job = Mycelia.build_execution_job(
            cmd = Mycelia.command_string(command),
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user,
        )
        submission = Mycelia.execute(job, resolved_executor)
        return (;
            status = :planned,
            submission,
            outdir = planned_outdir,
            assembly,
            graph,
            contract = nothing,
            expected_artifacts = (;
                assembly,
                graph,
            ),
            requested_outdir = reported_outdir,
            toolchain = nothing,
            provenance_status = "planned-unrealized",
            conda_runner = resolved_conda_runner,
            environment_prefix = resolved_environment_prefix,
        )
    end
    lock_path = environment_lock_path === nothing ?
                _unicycler_environment_lock_path_from_prefix(
        resolved_environment_prefix,
    ) :
                String(environment_lock_path)
    return output_lock_runner(planned_outdir) do reserved_outdir
        reserved_outdir = _require_unchanged_unicycler_output_plan(
            normalized_outdir,
            reserved_outdir,
        )
        return environment_lock_runner(lock_path) do
            pre_environment_validation()
            _require_unicycler_pre_environment_output_shape(reserved_outdir)
            _invoke_unicycler_environment_callback(
                environment_preparer,
                resolved_conda_runner,
                resolved_environment_prefix,
            )
            toolchain_before = _require_unicycler_toolchain_provenance(
                _invoke_unicycler_environment_callback(
                    toolchain_inspector,
                    resolved_conda_runner,
                    resolved_environment_prefix,
                ),
            )
            assembly = joinpath(reserved_outdir, "assembly.fasta")
            graph = joinpath(reserved_outdir, "assembly.gfa")
            contract_marker = joinpath(
                reserved_outdir,
                _UNICYCLER_CONTRACT_FILENAME,
            )
            if ispath(assembly) || islink(assembly) ||
               ispath(graph) || islink(graph) ||
               ispath(contract_marker) || islink(contract_marker)
                bound_artifacts = _unicycler_semantic_artifact_snapshots(
                    assembly,
                    graph,
                    reserved_outdir;
                    after_semantic_validation_hook =
                        after_artifact_semantic_validation_hook,
                )
                (ispath(contract_marker) || islink(contract_marker)) || error(
                    "Refusing unbound legacy Unicycler output reuse: missing " *
                    "durable run-contract marker $(repr(contract_marker)).",
                )
                assembly = bound_artifacts.assembly
                graph = bound_artifacts.graph
                artifact_snapshots = bound_artifacts.snapshots
                input_fingerprints = _unicycler_observed_input_fingerprints(
                    short_1,
                    short_2,
                    long_reads,
                    consumed_short_1,
                    consumed_short_2,
                    consumed_long_reads,
                    consumed_input_fingerprints,
                    input_fingerprinter,
                )
                expected_contract = _unicycler_run_contract(
                    input_fingerprints,
                    threads,
                    spades_options,
                    kmers,
                    resolved_environment_prefix,
                    toolchain_before,
                    artifact_snapshots,
                )
                contract_marker = _require_matching_unicycler_run_contract(
                    reserved_outdir,
                    expected_contract,
                )
                final_input_fingerprints = _unicycler_observed_input_fingerprints(
                    short_1,
                    short_2,
                    long_reads,
                    consumed_short_1,
                    consumed_short_2,
                    consumed_long_reads,
                    consumed_input_fingerprints,
                    input_fingerprinter,
                )
                final_input_fingerprints == input_fingerprints || error(
                    "Unicycler input content changed while contract-verified " *
                    "output was being reused; refusing stale provenance.",
                )
                final_artifacts = _require_unchanged_unicycler_artifacts(
                    artifact_snapshots,
                    assembly,
                    graph,
                    reserved_outdir,
                )
                assembly = final_artifacts.assembly
                graph = final_artifacts.graph
                contract_marker = _require_matching_unicycler_run_contract(
                    reserved_outdir,
                    expected_contract,
                )
                _require_unchanged_unicycler_output_plan(
                    normalized_outdir,
                    reserved_outdir,
                )
                return (;
                    status = :reused,
                    outdir = reserved_outdir,
                    assembly,
                    graph,
                    contract = contract_marker,
                    requested_outdir = reported_outdir,
                    toolchain = toolchain_before,
                    provenance_status = "reused-verified-contract",
                    conda_runner = resolved_conda_runner,
                    environment_prefix = resolved_environment_prefix,
                )
            end

            command = _unicycler_command(;
                conda_runner = resolved_conda_runner,
                environment_prefix = resolved_environment_prefix,
                short_1 = consumed_short_1,
                short_2 = consumed_short_2,
                long_reads = consumed_long_reads,
                outdir = reserved_outdir,
                threads,
                spades_options,
                kmers,
            )
            if islink(reserved_outdir) ||
               (ispath(reserved_outdir) && !isdir(reserved_outdir))
                throw(ArgumentError(
                    "Unicycler outdir must be absent or an empty directory: " *
                    repr(reserved_outdir),
                ))
            end
            if isdir(reserved_outdir)
                _require_exact_unicycler_output_dir(reserved_outdir)
                isempty(readdir(reserved_outdir)) || throw(ArgumentError(
                    "Refusing to remove non-empty Unicycler outdir without a " *
                    "completed assembly: $(repr(reserved_outdir)).",
                ))
                try
                    rm(reserved_outdir)
                catch removal_error
                    removal_error isa InterruptException && rethrow()
                    throw(ArgumentError(
                        "Unicycler outdir changed or was not empty during " *
                        "nonrecursive removal: $(repr(reserved_outdir)).",
                    ))
                end
            end
            input_fingerprints = _unicycler_observed_input_fingerprints(
                short_1,
                short_2,
                long_reads,
                consumed_short_1,
                consumed_short_2,
                consumed_long_reads,
                consumed_input_fingerprints,
                input_fingerprinter,
            )
            command_runner(command)
            _require_unchanged_unicycler_output_plan(
                normalized_outdir,
                reserved_outdir,
            )
            bound_artifacts = _unicycler_semantic_artifact_snapshots(
                assembly,
                graph,
                reserved_outdir;
                after_semantic_validation_hook =
                    after_artifact_semantic_validation_hook,
            )
            assembly = bound_artifacts.assembly
            graph = bound_artifacts.graph
            artifact_snapshots = bound_artifacts.snapshots
            toolchain_after = _require_unicycler_toolchain_provenance(
                _invoke_unicycler_environment_callback(
                    toolchain_inspector,
                    resolved_conda_runner,
                    resolved_environment_prefix,
                ),
            )
            toolchain_before == toolchain_after || error(
                "Unicycler realized Conda package inventory changed while the " *
                "assembler ran; refusing stale toolchain provenance.",
            )
            _require_unchanged_unicycler_output_plan(
                normalized_outdir,
                reserved_outdir,
            )
            final_artifacts = _require_unchanged_unicycler_artifacts(
                artifact_snapshots,
                assembly,
                graph,
                reserved_outdir,
            )
            assembly = final_artifacts.assembly
            graph = final_artifacts.graph
            observed_inputs = _unicycler_observed_input_fingerprints(
                short_1,
                short_2,
                long_reads,
                consumed_short_1,
                consumed_short_2,
                consumed_long_reads,
                consumed_input_fingerprints,
                input_fingerprinter,
            )
            observed_inputs == input_fingerprints || error(
                "Unicycler input content changed while the assembler ran; " *
                "refusing a stale run contract.",
            )
            final_artifacts = _require_unchanged_unicycler_artifacts(
                artifact_snapshots,
                assembly,
                graph,
                reserved_outdir,
            )
            assembly = final_artifacts.assembly
            graph = final_artifacts.graph
            contract = _unicycler_run_contract(
                input_fingerprints,
                threads,
                spades_options,
                kmers,
                resolved_environment_prefix,
                toolchain_after,
                artifact_snapshots,
            )
            contract_marker = _write_unicycler_run_contract(
                reserved_outdir,
                contract,
            )
            final_input_fingerprints = _unicycler_observed_input_fingerprints(
                short_1,
                short_2,
                long_reads,
                consumed_short_1,
                consumed_short_2,
                consumed_long_reads,
                consumed_input_fingerprints,
                input_fingerprinter,
            )
            final_input_fingerprints == input_fingerprints || error(
                "Unicycler input content changed while its durable contract " *
                "was published; refusing stale lifecycle provenance.",
            )
            final_artifacts = _require_unchanged_unicycler_artifacts(
                artifact_snapshots,
                assembly,
                graph,
                reserved_outdir,
            )
            assembly = final_artifacts.assembly
            graph = final_artifacts.graph
            contract_marker = _require_matching_unicycler_run_contract(
                reserved_outdir,
                contract,
            )
            _require_unchanged_unicycler_output_plan(
                normalized_outdir,
                reserved_outdir,
            )
            return (;
                status = :completed,
                outdir = reserved_outdir,
                assembly,
                graph,
                contract = contract_marker,
                requested_outdir = reported_outdir,
                toolchain = toolchain_after,
                provenance_status = "realized-local-exact",
                conda_runner = resolved_conda_runner,
                environment_prefix = resolved_environment_prefix,
            )
        end
    end
end

function _run_unicycler(;
        short_1::AbstractString,
        short_2::Union{Nothing, AbstractString} = nothing,
        long_reads::AbstractString,
        outdir::AbstractString = "unicycler_output",
        threads::Int = get_default_threads(),
        spades_options::Union{Nothing, String} = nothing,
        kmers::Union{Nothing, String} = nothing,
        executor::Any = nothing,
        site::Symbol = :local,
        job_name::String = "unicycler",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing,
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        input_spool_parent::Union{Nothing, AbstractString} = nothing,
        input_spool_byte_ceiling::Integer =
            _UNICYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
        input_spool_available_bytes_reader::Function =
            parent -> diskstat(parent).available,
        after_input_spool_copy_hook::Function = binding -> nothing,
        after_input_spool_materialization_hook::Function = snapshots -> nothing,
        environment_preparer::Function = _prepare_unicycler_environment,
        toolchain_inspector::Function = (runner, prefix) ->
            _unicycler_toolchain_provenance(
                inventory_reader = () ->
                    _unicycler_conda_package_inventory(
                        conda_runner = runner,
                        environment_prefix = prefix,
                    ),
            ),
        input_fingerprinter::Function = _unicycler_input_fingerprint,
        after_artifact_semantic_validation_hook::Function =
            artifacts -> nothing,
        environment_lock_path::Union{Nothing, AbstractString} = nothing,
        environment_lock_runner::Function = _with_unicycler_environment_lock,
        output_lock_runner::Function = _with_unicycler_output_lock,
        command_runner::Function = run,
)::NamedTuple
    validated_inputs = _validate_unicycler_inputs(
        short_1,
        short_2,
        long_reads,
    )
    resolved_executor = resolve_executor(executor)
    function run_with_inputs(
            consumed_short_1::AbstractString,
            consumed_short_2::Union{Nothing, AbstractString},
            consumed_long_reads::AbstractString,
            consumed_input_fingerprints::Union{Nothing, AbstractDict},
            reserved_output_lock_runner::Function,
            pre_environment_validation::Function,
    )::NamedTuple
        return _run_unicycler_with_contract(;
            short_1 = validated_inputs.short_1,
            short_2 = validated_inputs.short_2,
            long_reads = validated_inputs.long_reads,
            consumed_short_1,
            consumed_short_2,
            consumed_long_reads,
            consumed_input_fingerprints,
            outdir,
            threads,
            spades_options,
            kmers,
            executor = resolved_executor,
            site,
            job_name,
            time_limit,
            partition,
            account,
            mem_gb,
            qos,
            mail_user,
            conda_runner,
            environment_prefix,
            environment_preparer,
            toolchain_inspector,
            input_fingerprinter,
            after_artifact_semantic_validation_hook,
            pre_environment_validation,
            environment_lock_path,
            environment_lock_runner,
            output_lock_runner = reserved_output_lock_runner,
            command_runner,
        )
    end
    if !(resolved_executor isa LocalExecutor)
        return run_with_inputs(
            validated_inputs.short_1,
            validated_inputs.short_2,
            validated_inputs.long_reads,
            nothing,
            output_lock_runner,
            () -> nothing,
        )
    end

    normalized_outdir = normpath(abspath(String(outdir)))
    planned_outdir = _canonical_planned_unicycler_output_path(
        normalized_outdir,
    )
    return output_lock_runner(planned_outdir) do reserved_outdir
        reserved_outdir = _require_unchanged_unicycler_output_plan(
            normalized_outdir,
            reserved_outdir,
        )
        resolved_spool_parent = input_spool_parent === nothing ?
                                _unicycler_output_adjacent_spool_parent(
            reserved_outdir,
        ) : String(input_spool_parent)
        function held_output_lock_runner(
                action::Function,
                requested_outdir::AbstractString,
        )::Any
            inner_reserved_outdir = _require_unchanged_unicycler_output_plan(
                requested_outdir,
                reserved_outdir,
            )
            inner_reserved_outdir == reserved_outdir || error(
                "Unicycler inner lifecycle escaped its held output " *
                "reservation.",
            )
            return action(reserved_outdir)
        end
        return _with_unicycler_stable_input_snapshots(
            validated_inputs.short_1,
            validated_inputs.short_2,
            validated_inputs.long_reads;
            spool_parent = resolved_spool_parent,
            byte_ceiling = input_spool_byte_ceiling,
            available_bytes_reader = input_spool_available_bytes_reader,
            after_copy_hook = after_input_spool_copy_hook,
        ) do snapshots
            after_input_spool_materialization_hook(snapshots)
            validated_snapshots =
                _validate_unicycler_stable_input_snapshots(
                snapshots,
                validated_inputs.short_count,
                validated_inputs.long_count,
            )
            _require_unicycler_sources_match_stable_snapshots(
                validated_snapshots,
            )
            function validate_held_inputs_before_environment()::Nothing
                _require_unchanged_unicycler_stable_input_snapshots(
                    validated_snapshots,
                )
                _require_unicycler_sources_match_stable_snapshots(
                    validated_snapshots,
                )
                return nothing
            end
            return run_with_inputs(
                validated_snapshots.short_1,
                validated_snapshots.short_2,
                validated_snapshots.long_reads,
                validated_snapshots.input_fingerprints,
                held_output_lock_runner,
                validate_held_inputs_before_environment,
            )
        end
    end
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run hybrid assembly combining short and long reads using Unicycler.

# Arguments
- `short_1::AbstractString`: Path to first short read FASTQ file
- `short_2::Union{Nothing,AbstractString}`: Optional second short read FASTQ file
- `long_reads::AbstractString`: Path to long read FASTQ file
- `outdir::AbstractString`: Output directory path (default: "unicycler_output")
- `threads::Int`: Number of threads to use (default: `get_default_threads()`)
- `spades_options::Union{Nothing,String}`: Extra SPAdes options (default: `nothing`)
- `kmers::Union{Nothing,String}`: Explicit SPAdes k-mers (e.g., "21,33,55")
- `conda_runner::AbstractString`: Conda-compatible executable used for every
  environment operation.
- `environment_prefix::Union{Nothing,AbstractString}`: Exact environment prefix
  used for provisioning, inventory, execution, and lifecycle locking.
- `input_spool_parent::Union{Nothing,AbstractString}`: Existing non-symlink
  scratch directory for local workflow-owned input snapshots. The default uses
  an output-adjacent parent selected only after the output reservation is held.
- `input_spool_byte_ceiling::Integer`: Maximum cumulative bytes permitted for
  the local input snapshots before any scratch directory or environment is
  created (default: 1 TB).

# Returns
Named tuple containing:
- `status::Symbol`: `:completed`, `:reused`, or `:planned`
- `outdir::String`: Canonical physical output path; validated and realized for
  completed or reused work and expected but unrealized for planned work
- `assembly::String`: Canonical final-assembly path, expected but unrealized for
  planned work
- `graph::String`: Canonical assembly-graph path, expected but unrealized for
  planned work
- `contract::Union{Nothing,String}`: Durable contract marker path for completed
  or reused local work; `nothing` for planned work because its raw command does
  not create Mycelia's reusable lifecycle marker
- `requested_outdir::String`: Caller-supplied output spelling retained for audit
- `toolchain::Union{Nothing,Dict}`: Exact realized package inventory for a
  completed or contract-verified reused run; `nothing` for a plan
- `provenance_status::String`: Why exact toolchain provenance is available or
  intentionally unavailable
- Planned results additionally contain `submission` and `expected_artifacts`.

# Details
- Automatically creates and uses one explicit prefix-addressed Conda
  environment with Unicycler.
- Combines short read accuracy with long read scaffolding
- Reuses output only when the semantic FASTA/GFA companion pair and an atomic
  durable contract exactly match canonical input paths, sizes, SHA-256 values,
  command-relevant parameters, realized toolchain, and artifact hashes. Legacy
  unbound output fails loudly. Fresh local runs report a locked before/after
  package inventory.
- `CollectExecutor`, `DryRunExecutor`, and `SlurmExecutor(dry_run=true)` return
  `status=:planned`, the submission result, and only the FASTA/GFA paths
  actually produced by the planned raw command. Such output has no reusable
  Mycelia contract marker and therefore cannot be reused by a later local call.
  Real nonlocal submission remains disabled until the remote runtime can hold
  the same lifecycle locks and finalize the durable contract.
- Rejects an incomplete non-empty output directory instead of deleting
  caller-owned contents. An empty output directory may be removed because
  Unicycler requires the output path to be absent.
- Streams every local FASTQ before any environment or output mutation, rejects
  malformed data, physical aliases, paired-count/order/identifier drift, and
  inverted explicit mate roles. After reserving the output, local execution
  creates bounded workflow-owned stable byte snapshots in configured scratch,
  revalidates their FASTQ semantics, record counts, and hashes, and only then
  prepares the environment. The command never consumes caller paths; the
  durable contract retains their canonical provenance and the exact
  consumed-snapshot hashes. Snapshot cleanup is guaranteed on success or
  failure.
- Returns `outdir`, `assembly`, and `graph` as canonical absolute safety fields,
  even when the caller supplied a relative path or a stable parent-directory
  alias. `requested_outdir` is audit-only caller spelling and must not be used
  for artifact access or cleanup authority.
- Utilizes all available CPU threads
"""
function run_unicycler(;
        short_1::AbstractString,
        short_2::Union{Nothing, AbstractString} = nothing,
        long_reads::AbstractString,
        outdir::AbstractString = "unicycler_output",
        threads::Int = get_default_threads(),
        spades_options::Union{Nothing, String} = nothing,
        kmers::Union{Nothing, String} = nothing,
        executor::Any = nothing,
        site::Symbol = :local,
        job_name::String = "unicycler",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing,
        conda_runner::AbstractString = _conda_runner(),
        environment_prefix::Union{Nothing, AbstractString} = nothing,
        input_spool_parent::Union{Nothing, AbstractString} = nothing,
        input_spool_byte_ceiling::Integer =
            _UNICYCLER_DEFAULT_INPUT_SPOOL_BYTE_CEILING,
)::NamedTuple
    return _run_unicycler(;
        short_1,
        short_2,
        long_reads,
        outdir,
        threads,
        spades_options,
        kmers,
        executor,
        site,
        job_name,
        time_limit,
        partition,
        account,
        mem_gb,
        qos,
        mail_user,
        conda_runner,
        environment_prefix,
        input_spool_parent,
        input_spool_byte_ceiling,
    )
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
function install_hybracter_db(db_dir::AbstractString; force::Bool = false)
    db_dir = abspath(String(db_dir))
    if isdir(db_dir) && !force && !isempty(readdir(db_dir))
        @info "Hybracter database already present; skipping download." db_dir
        return db_dir
    end
    if force && isdir(db_dir)
        rm(db_dir; recursive = true, force = true)
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
        db_dir::Union{Nothing, AbstractString} = nothing,
        threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[],
        force::Bool = false)
    isfile(long_reads) || error("Long read file not found: $(long_reads)")
    isfile(read1) || error("Read 1 file not found: $(read1)")
    isfile(read2) || error("Read 2 file not found: $(read2)")
    isempty(sample_name) && error("sample_name cannot be empty")

    outdir = abspath(String(outdir))
    sample_dir = joinpath(outdir, sample_name)
    complete_final = joinpath(outdir, "FINAL_OUTPUT", "complete", sample_name *
                                                                  "_final.fasta")
    incomplete_final = joinpath(outdir, "FINAL_OUTPUT", "incomplete", sample_name *
                                                                      "_final.fasta")
    legacy_final = joinpath(sample_dir, "final_assembly.fasta")
    complete_summary = joinpath(outdir, "FINAL_OUTPUT", "complete", sample_name *
                                                                    "_summary.tsv")
    incomplete_summary = joinpath(outdir, "FINAL_OUTPUT", "incomplete", sample_name *
                                                                        "_summary.tsv")
    legacy_summary = joinpath(sample_dir, "summary.tsv")

    if !force &&
       (isfile(complete_final) || isfile(incomplete_final) || isfile(legacy_final))
        final_fasta = isfile(complete_final) ? complete_final :
                      (isfile(incomplete_final) ? incomplete_final : legacy_final)
        summary_tsv = isfile(complete_summary) ? complete_summary :
                      (isfile(incomplete_summary) ? incomplete_summary : legacy_summary)
        return (; outdir, sample_dir, summary_tsv, final_fasta)
    end

    if force && isdir(sample_dir)
        rm(sample_dir; recursive = true, force = true)
    end

    if db_dir === nothing || isempty(db_dir)
        db_dir = joinpath(homedir(), "workspace", ".hybracter")
    end
    if !isdir(db_dir) || isempty(readdir(db_dir))
        install_hybracter_db(db_dir; force = force)
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
    final_fasta = isfile(complete_final) ? complete_final :
                  (isfile(incomplete_final) ? incomplete_final : legacy_final)
    summary_tsv = isfile(complete_summary) ? complete_summary :
                  (isfile(incomplete_summary) ? incomplete_summary : legacy_summary)
    return (; outdir, sample_dir, summary_tsv, final_fasta)
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
        db_dir::Union{Nothing, AbstractString} = nothing,
        threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[],
        force::Bool = false)
    isfile(long_reads) || error("Long read file not found: $(long_reads)")
    isempty(sample_name) && error("sample_name cannot be empty")

    outdir = abspath(String(outdir))
    sample_dir = joinpath(outdir, sample_name)
    complete_final = joinpath(outdir, "FINAL_OUTPUT", "complete", sample_name *
                                                                  "_final.fasta")
    incomplete_final = joinpath(outdir, "FINAL_OUTPUT", "incomplete", sample_name *
                                                                      "_final.fasta")
    legacy_final = joinpath(sample_dir, "final_assembly.fasta")
    complete_summary = joinpath(outdir, "FINAL_OUTPUT", "complete", sample_name *
                                                                    "_summary.tsv")
    incomplete_summary = joinpath(outdir, "FINAL_OUTPUT", "incomplete", sample_name *
                                                                        "_summary.tsv")
    legacy_summary = joinpath(sample_dir, "summary.tsv")

    if !force &&
       (isfile(complete_final) || isfile(incomplete_final) || isfile(legacy_final))
        final_fasta = isfile(complete_final) ? complete_final :
                      (isfile(incomplete_final) ? incomplete_final : legacy_final)
        summary_tsv = isfile(complete_summary) ? complete_summary :
                      (isfile(incomplete_summary) ? incomplete_summary : legacy_summary)
        return (; outdir, sample_dir, summary_tsv, final_fasta)
    end

    if force && isdir(sample_dir)
        rm(sample_dir; recursive = true, force = true)
    end

    if db_dir === nothing || isempty(db_dir)
        db_dir = joinpath(homedir(), "workspace", ".hybracter")
    end
    if !isdir(db_dir) || isempty(readdir(db_dir))
        install_hybracter_db(db_dir; force = force)
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
    final_fasta = isfile(complete_final) ? complete_final :
                  (isfile(incomplete_final) ? incomplete_final : legacy_final)
    summary_tsv = isfile(complete_summary) ? complete_summary :
                  (isfile(incomplete_summary) ? incomplete_summary : legacy_summary)
    return (; outdir, sample_dir, summary_tsv, final_fasta)
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
        prefix::AbstractString = "dnaapler",
        threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[],
        force::Bool = false)
    isfile(input_file) || error("Input file not found: $(input_file)")

    outdir = abspath(String(outdir))
    mkpath(outdir)

    lower_name = lowercase(input_file)
    ext = (endswith(lower_name, ".gfa") || endswith(lower_name, ".gfa.gz")) ? ".gfa" :
          ".fasta"
    reoriented = joinpath(outdir, "$(prefix)_reoriented$(ext)")

    if !force && isfile(reoriented)
        return (; outdir, reoriented)
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
    return (; outdir, reoriented)
end

function _plassembler_db_ready(db_dir::AbstractString)
    return isdir(db_dir) &&
           any(path -> endswith(path, ".msh"), readdir(db_dir; join = true))
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
function install_plassembler_db(db_dir::AbstractString; force::Bool = false)
    db_dir = abspath(String(db_dir))
    db_ready = _plassembler_db_ready(db_dir)
    if db_ready && !force
        @info "Plassembler database already present; skipping download." db_dir
        return db_dir
    end
    if isdir(db_dir) && (!db_ready || force)
        rm(db_dir; recursive = true, force = true)
    end
    mkpath(db_dir)
    Mycelia.add_bioconda_env("plassembler")
    try
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plassembler plassembler download -d $(db_dir) --force`)
    catch e
        tarballs = filter(path -> endswith(path, ".tar.gz"), readdir(db_dir; join = true))
        if isempty(tarballs)
            rethrow(e)
        end
        tarball = length(tarballs) == 1 ? tarballs[1] : sort(tarballs)[end]
        @warn "Plassembler download failed; attempting manual extraction." exception=(
            e, catch_backtrace())
        tmp_dir = mktempdir(dirname(db_dir))
        try
            Tar.extract(tarball, tmp_dir)
        catch extract_error
            rm(tmp_dir; recursive = true, force = true)
            rethrow(extract_error)
        end
        rm(db_dir; recursive = true, force = true)
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
        db_dir::Union{Nothing, AbstractString} = nothing,
        threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[],
        force::Bool = false)
    isfile(long_reads) || error("Long read file not found: $(long_reads)")
    isfile(read1) || error("Read 1 file not found: $(read1)")
    isfile(read2) || error("Read 2 file not found: $(read2)")

    outdir = abspath(String(outdir))
    plasmids_fasta = joinpath(outdir, "plassembler_plasmids.fasta")
    chromosome_fasta = joinpath(outdir, "chromosome.fasta")
    summary_tsv = joinpath(outdir, "plassembler_summary.tsv")

    if !force && isfile(plasmids_fasta)
        return (; outdir, plasmids_fasta, chromosome_fasta, summary_tsv)
    end

    if force && isdir(outdir)
        rm(outdir; recursive = true, force = true)
    end
    mkpath(dirname(outdir))

    if db_dir === nothing || isempty(db_dir)
        db_dir = joinpath(homedir(), "workspace", ".plassembler")
    end
    if !_plassembler_db_ready(db_dir)
        install_plassembler_db(db_dir; force = force)
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
    return (; outdir, plasmids_fasta, chromosome_fasta, summary_tsv)
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
        db_dir::Union{Nothing, AbstractString} = nothing,
        threads::Int = get_default_threads(),
        extra_args::Vector{String} = String[],
        force::Bool = false)
    isfile(long_reads) || error("Long read file not found: $(long_reads)")

    outdir = abspath(String(outdir))
    plasmids_fasta = joinpath(outdir, "plassembler_plasmids.fasta")
    chromosome_fasta = joinpath(outdir, "chromosome.fasta")
    summary_tsv = joinpath(outdir, "plassembler_summary.tsv")

    if !force && isfile(plasmids_fasta)
        return (; outdir, plasmids_fasta, chromosome_fasta, summary_tsv)
    end

    if force && isdir(outdir)
        rm(outdir; recursive = true, force = true)
    end
    mkpath(dirname(outdir))

    if db_dir === nothing || isempty(db_dir)
        db_dir = joinpath(homedir(), "workspace", ".plassembler")
    end
    if !_plassembler_db_ready(db_dir)
        install_plassembler_db(db_dir; force = force)
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
    return (; outdir, plasmids_fasta, chromosome_fasta, summary_tsv)
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
function run_plass_assemble(; reads1::String, reads2::Union{String, Nothing} = nothing,
        outdir::String = "plass_output",
        min_seq_id::Union{Real, Nothing} = nothing, min_length::Union{Int, Nothing} = nothing, evalue::Union{
            Real, Nothing} = nothing,
        num_iterations::Union{Int, Nothing} = nothing, filter_proteins::Bool = true,
        threads::Int = get_default_threads(),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "plass_assemble",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
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

    if executor !== nothing
        cmd_string = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`
        )
        job = Mycelia.build_execution_job(
            cmd = cmd_string,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, assembly = assembly_out, tmpdir)
    end

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`)

    return (; outdir, assembly = assembly_out, tmpdir)
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
function run_penguin_guided_nuclassemble(;
        reads1::String, reads2::Union{String, Nothing} = nothing,
        outdir::String = "penguin_guided_output",
        threads::Int = get_default_threads(),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "penguin_guided_nuclassemble",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
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

    if executor !== nothing
        cmd_string = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`
        )
        job = Mycelia.build_execution_job(
            cmd = cmd_string,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, assembly = assembly_out, tmpdir)
    end

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`)

    return (; outdir, assembly = assembly_out, tmpdir)
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
function run_penguin_nuclassemble(;
        reads1::String, reads2::Union{String, Nothing} = nothing,
        outdir::String = "penguin_output",
        threads::Int = get_default_threads(),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "penguin_nuclassemble",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
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

    if executor !== nothing
        cmd_string = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`
        )
        job = Mycelia.build_execution_job(
            cmd = cmd_string,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return (; outdir, assembly = assembly_out, tmpdir)
    end

    run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n plass $cmd`)

    return (; outdir, assembly = assembly_out, tmpdir)
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
function run_apollo(assembly_file::String, reads_file::String; outdir::String = assembly_file *
                                                                                "_apollo")
    Mycelia.add_bioconda_env("apollo")
    Mycelia.add_bioconda_env("minimap2")
    Mycelia.add_bioconda_env("samtools")
    mkpath(outdir)

    basename_assembly = splitext(basename(assembly_file))[1]
    polished_assembly = joinpath(outdir, basename_assembly * "_polished.fasta")

    apollo_prefix = _conda_env_prefix("apollo")
    apollo_bin = apollo_prefix === nothing ? nothing :
                 joinpath(apollo_prefix, "bin", "apollo")
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
            if unicycler_prefix !== nothing &&
               isfile(joinpath(unicycler_prefix, "bin", "racon"))
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
                run(pipeline(racon_cmd, stdout = io))
            end
        end
    end

    return (; outdir, polished_assembly)
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

function run_homopolish(assembly_file::String, local_db_path::String;
        outdir::String = assembly_file * "_homopolish", model_path::String = "",
        threads = get_default_threads(), sketch_path::String = "")
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
        cmd_args = ["homopolish", "polish", "-a", assembly_file, "-s", sketch_path,
            "-m", model_path, "-o", outdir, "-t", string(threads)]

        # Run Homopolish
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n homopolish $(cmd_args)`)

        # Move output to expected location if needed
        homopolish_output = joinpath(outdir, "homopolished_" * basename_assembly)
        if isfile(homopolish_output) && !isfile(polished_assembly)
            mv(homopolish_output, polished_assembly)
        end
    end

    return (; outdir, polished_assembly)
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
- `nsplit::Union{Nothing,Int}`: Number of split input chunks (default: nothing)
- `min_identity::Union{Nothing,Float64}`: Minimum overlap identity (default: nothing)
- `min_ovlp_len::Union{Nothing,Int}`: Minimum overlap length (default: nothing)
- `size::Union{Nothing,Int}`: Max cluster size for short reads (default: nothing)
- `max_tip_len::Union{Nothing,Int}`: Max tip length (default: nothing)
- `insert_size::Union{Nothing,Int}`: Short-read insert size (default: nothing)
- `average_read_len::Union{Nothing,Int}`: Short-read average length (default: nothing)
- `corrected::Bool`: Whether reads are already corrected (default: false)
- `low_q::Bool`: Whether short reads are low quality (default: false)

# Returns
Named tuple containing:
- `outdir::String`: Path to output directory
- `strain_assemblies::String`: Path to strain-resolved assemblies directory

# Details
- Uses hybrid overlap graphs combining short and long reads
- Implements "cross hybrid" mutual support strategy for strain resolution
- Separates closely related strains in metagenomic samples
- Automatically creates and uses a conda environment with hylight
- Concatenates paired short-read inputs to match HyLight's single-file interface
- Utilizes requested CPU threads for HyLight run
- Skips assembly if output directory already exists
"""
function run_hylight(short_reads_1::String, short_reads_2::String, long_reads::String;
        outdir::String = "hylight_output",
        threads = get_default_threads(),
        nsplit::Union{Nothing, Int} = nothing,
        min_identity::Union{Nothing, Float64} = nothing,
        min_ovlp_len::Union{Nothing, Int} = nothing,
        size::Union{Nothing, Int} = nothing,
        max_tip_len::Union{Nothing, Int} = nothing,
        insert_size::Union{Nothing, Int} = nothing,
        average_read_len::Union{Nothing, Int} = nothing,
        corrected::Bool = false,
        low_q::Bool = false
)
    isfile(short_reads_1) || error("Short read 1 file not found: $(short_reads_1)")
    isfile(short_reads_2) || error("Short read 2 file not found: $(short_reads_2)")
    isfile(long_reads) || error("Long read file not found: $(long_reads)")

    Mycelia.add_bioconda_env("hylight")
    mkpath(outdir)

    # HyLight expects a single short-read FASTQ; combine paired inputs.
    short_reads = joinpath(outdir, "short_reads.fq")
    open(short_reads, "w") do io
        for path in (short_reads_1, short_reads_2)
            open(path, "r") do input
                Base.copyto!(io, input)
            end
        end
    end

    strain_assemblies = joinpath(outdir, "strain_assemblies")

    if !isdir(strain_assemblies)
        # Run HyLight assembly
        cmd_args = [
            "hylight",
            "-s", short_reads,
            "-l", long_reads,
            "-o", outdir,
            "-t", string(threads)
        ]
        corrected && push!(cmd_args, "--corrected")
        low_q && push!(cmd_args, "--low_q")
        isnothing(nsplit) || push!(cmd_args, "--nsplit", string(nsplit))
        isnothing(min_identity) || push!(cmd_args, "--min_identity", string(min_identity))
        isnothing(min_ovlp_len) || push!(cmd_args, "--min_ovlp_len", string(min_ovlp_len))
        isnothing(size) || push!(cmd_args, "--size", string(size))
        isnothing(max_tip_len) || push!(cmd_args, "--max_tip_len", string(max_tip_len))
        isnothing(insert_size) || push!(cmd_args, "--insert_size", string(insert_size))
        isnothing(average_read_len) ||
            push!(cmd_args, "--average_read_len", string(average_read_len))

        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n hylight $(cmd_args)`)
    end

    return (; outdir, strain_assemblies)
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
function run_strong(assembly_graph::String, reads_file::String;
        outdir::String = "strong_output", nb_strains::Int = 0)
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

    return (; outdir, strain_unitigs)
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
function run_velvet(fastq1::String; fastq2::Union{String, Nothing} = nothing,
        outdir::String = "velvet_output", k::Int = 31,
        exp_cov::String = "auto", min_contig_lgth::Int = 200)
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

    return (; outdir, contigs = contigs_file)
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
function run_metavelvet(fastq1::String; fastq2::Union{String, Nothing} = nothing,
        outdir::String = "metavelvet_output", k::Int = 31,
        exp_cov::String = "auto", min_contig_lgth::Int = 200)
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
        velvetg_args = ["velvetg", outdir, "-read_trkg", "yes",
            "-min_contig_lgth", string(min_contig_lgth)]
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

    return (; outdir, contigs = contigs_file)
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
function run_velvetg(velvet_dir::String; outdir::String = velvet_dir,
        exp_cov::Union{String, Float64} = "auto",
        cov_cutoff::Union{String, Float64} = "auto",
        min_contig_lgth::Int = 200, scaffolding::Bool = true)
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
                cp(src_file, dst_file, force = true)
            end
        end
    end

    return (; outdir, contigs = contigs_file, stats = stats_file)
end

function _metamdbg_selected_input(
        hifi_reads::Union{String, Vector{String}, Nothing},
        ont_reads::Union{String, Vector{String}, Nothing},
        ont_r10_4_plus::Bool = false,
)::NamedTuple
    has_hifi_reads = hifi_reads !== nothing
    has_ont_reads = ont_reads !== nothing
    if has_hifi_reads == has_ont_reads
        throw(ArgumentError(
            "metaMDBG requires exactly one input technology: provide " *
            "either hifi_reads or ont_reads, but not both.",
        ))
    end
    if has_hifi_reads && ont_r10_4_plus
        throw(ArgumentError(
            "ont_r10_4_plus applies only when ont_reads is selected.",
        ))
    elseif has_ont_reads && !ont_r10_4_plus
        throw(ArgumentError(
            "metaMDBG ont_reads require ont_r10_4_plus=true as an explicit " *
            "attestation that every input was generated with supported " *
            "Nanopore R10.4 or later chemistry; generic, R9, and unknown ONT " *
            "inputs are not supported.",
        ))
    end

    selected_reads = has_hifi_reads ? hifi_reads : ont_reads
    selected_paths = selected_reads isa String ? [selected_reads] : selected_reads
    technology = has_hifi_reads ? "hifi_reads" : "ont_reads"
    flag = has_hifi_reads ? "--in-hifi" : "--in-ont"
    isempty(selected_paths) && throw(ArgumentError(
        "metaMDBG $(technology) must contain at least one non-empty path.",
    ))

    normalized_paths = String[]
    for path in selected_paths
        stripped_path = strip(path)
        isempty(stripped_path) && throw(ArgumentError(
            "metaMDBG $(technology) must contain at least one non-empty path.",
        ))
        normalized_path = normpath(abspath(stripped_path))
        isfile(normalized_path) || throw(ArgumentError(
            "metaMDBG $(technology) input file does not exist: $(normalized_path)",
        ))
        filesize(normalized_path) > 0 || throw(ArgumentError(
            "metaMDBG $(technology) input file is empty: $(normalized_path)",
        ))
        for existing_path in normalized_paths
            Base.Filesystem.samefile(normalized_path, existing_path) && throw(
                ArgumentError(
                    "metaMDBG $(technology) input files must refer to physically " *
                    "distinct files; $(normalized_path) aliases $(existing_path).",
                ),
            )
        end
        push!(normalized_paths, normalized_path)
    end
    platform_attestation = has_ont_reads ? "nanopore-r10.4-or-later" : nothing
    return (; flag, paths = normalized_paths, platform_attestation)
end

const METAMDBG_VERSION = "1.4"
const METAMDBG_ENVIRONMENT_SPEC_SHA256 =
    "3b51b282e8aa768da12e253af01dee43fa96a320baee96755ebb3c123723ff87"
const METAMDBG_ENV_NAME =
    "metamdbg-$(METAMDBG_VERSION)-$(first(METAMDBG_ENVIRONMENT_SPEC_SHA256, 16))"
const _METAMDBG_INSTALL_LOCK_STALE_SECONDS = 600
const _METAMDBG_INSTALL_LOCK_REFRESH_SECONDS = 60
const _METAMDBG_CONTRACT_SCHEMA_VERSION = 5
const _METAMDBG_CONTRACT_FILENAME = "mycelia_metamdbg_contract.json"
const _METAMDBG_COMPLETION_SCHEMA_VERSION = 1
const _METAMDBG_COMPLETION_FILENAME = "mycelia_metamdbg_completion.json"
const _METAMDBG_SUBMISSION_RESERVATION_SCHEMA_VERSION = 3
const _METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME = "contract.json"
const _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION = 1
const _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION = 2
const _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME = "job.json"
const _METAMDBG_PENDING_SUBMISSION_JOB_PREFIX =
    ".mycelia-metamdbg-job-pending."
const _METAMDBG_OUTPUT_ROOT_RESERVATION_SCHEMA_VERSION = 1
const _METAMDBG_OUTPUT_LOCK_RETRY_ATTEMPTS = 300
const _METAMDBG_OUTPUT_LOCK_RETRY_DELAY_SECONDS = 1.0
const _METAMDBG_ENUMERATION_MAX_DEPTH = 256
const _METAMDBG_ENUMERATION_MAX_ENTRIES = 1_000_000
const _METAMDBG_ENUMERATION_MAX_MANIFEST_BYTES = 268_435_456
const _METAMDBG_IUPAC_DNA_REGEX =
    r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$"
const _METAMDBG_GFA_IDENTIFIER_REGEX = r"^[!-)+-<>-~][!-~]*$"
const _METAMDBG_GFA_CIGAR_REGEX =
    r"^(?:\*|(?:[0-9]+[MIDNSHPX=])+)$"
const _METAMDBG_GFA_TAG_NAME_REGEX = r"^[A-Za-z][A-Za-z0-9]$"
const _METAMDBG_GFA_TAG_INTEGER_REGEX = r"^[-+]?[0-9]+$"
const _METAMDBG_GFA_TAG_FLOAT_REGEX =
    r"^[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?$"
const _METAMDBG_GFA_TAG_PRINTABLE_REGEX = r"^[ -~]+$"
const _METAMDBG_GFA_TAG_HEX_REGEX = r"^[0-9A-F]+$"
const _METAMDBG_GFA_B_INTEGER_RANGES = (
    c = (Int64(-128), Int64(127)),
    C = (Int64(0), Int64(255)),
    s = (Int64(-32_768), Int64(32_767)),
    S = (Int64(0), Int64(65_535)),
    i = (Int64(-2_147_483_648), Int64(2_147_483_647)),
    I = (Int64(0), Int64(4_294_967_295)),
)

function _metamdbg_paths()::Tuple{String, String}
    install_dir = joinpath(dirname(dirname(pathof(Mycelia))), "deps", "metamdbg")
    return install_dir, joinpath(install_dir, "environment.yml")
end

function _metamdbg_sha256(path::AbstractString)::String
    return open(path, "r") do input
        SHA.bytes2hex(SHA.sha256(input))
    end
end

function _metamdbg_string_sha256(value::AbstractString)::String
    return SHA.bytes2hex(SHA.sha256(String(value)))
end

function _require_verified_metamdbg_environment_spec(
        path::AbstractString,
)::String
    normalized_path = abspath(path)
    if !isfile(normalized_path) || filesize(normalized_path) == 0
        error(
            "Bundled metaMDBG environment specification is missing or empty: " *
            "$(normalized_path). Reinstall Mycelia before running metaMDBG.",
        )
    end
    actual_sha256 = _metamdbg_sha256(normalized_path)
    actual_sha256 == METAMDBG_ENVIRONMENT_SPEC_SHA256 || error(
        "metaMDBG environment specification checksum mismatch: expected " *
        "$(METAMDBG_ENVIRONMENT_SPEC_SHA256), got $(actual_sha256) for " *
        "$(normalized_path).",
    )
    return normalized_path
end

function _metamdbg_expected_toolchain()::NamedTuple
    return (
        metamdbg_version = METAMDBG_VERSION,
        environment_name = METAMDBG_ENV_NAME,
        environment_spec_sha256 = METAMDBG_ENVIRONMENT_SPEC_SHA256,
    )
end

function _metamdbg_package_record_field(
        package_record::Union{NamedTuple, AbstractDict},
        field::Symbol,
)::Any
    if package_record isa NamedTuple
        return hasproperty(package_record, field) ?
               getproperty(package_record, field) : nothing
    end
    return get(
        package_record,
        String(field),
        get(package_record, field, nothing),
    )
end

function _metamdbg_normalized_package_inventory(
        package_records::AbstractVector,
)::Vector{NamedTuple}
    isempty(package_records) && error(
        "metaMDBG conda package inventory must contain at least one package.",
    )
    inventory = NamedTuple[]
    for (record_index, package_record) in enumerate(package_records)
        package_record isa Union{NamedTuple, AbstractDict} || error(
            "metaMDBG conda package inventory record $(record_index) is not " *
            "an object.",
        )
        name = _metamdbg_package_record_field(package_record, :name)
        version = _metamdbg_package_record_field(package_record, :version)
        build = _metamdbg_package_record_field(package_record, :build_string)
        build === nothing &&
            (build = _metamdbg_package_record_field(package_record, :build))
        channel = _metamdbg_package_record_field(package_record, :channel)
        fields = (; name, version, build, channel)
        for (field_name, value) in pairs(fields)
            value isa AbstractString && !isempty(value) || error(
                "metaMDBG conda package inventory record $(record_index) has " *
                "a missing or empty $(field_name) field.",
            )
            occursin(r"[\t\r\n]", value) && error(
                "metaMDBG conda package inventory record $(record_index) has " *
                "a noncanonical $(field_name) field.",
            )
        end
        push!(inventory, (;
            name = String(name),
            version = String(version),
            build = String(build),
            channel = String(channel),
        ))
    end
    sort!(
        inventory;
        by = record -> (
            record.name,
            record.version,
            record.build,
            record.channel,
        ),
    )
    names = getproperty.(inventory, :name)
    allunique(names) || error(
        "metaMDBG conda package inventory contains duplicate package names.",
    )
    return inventory
end

function _metamdbg_package_inventory_contents(
        package_records::AbstractVector,
)::String
    inventory = _metamdbg_normalized_package_inventory(package_records)
    return join(
        map(inventory) do record
            join(
                (record.name, record.version, record.build, record.channel),
                '\t',
            )
        end,
        '\n',
    ) * "\n"
end

function _metamdbg_package_inventory_sha256(
        package_records::AbstractVector,
)::String
    contents = _metamdbg_package_inventory_contents(package_records)
    return SHA.bytes2hex(SHA.sha256(contents))
end

function _metamdbg_runtime_inventory_canonicalizer_python()::String
    return raw"""
import json
import pathlib
import sys


def fail(message: str) -> None:
    raise ValueError(message)


source = pathlib.Path(sys.argv[1])
destination = pathlib.Path(sys.argv[2])
with source.open("r", encoding="utf-8") as stream:
    records = json.load(stream)
if not isinstance(records, list) or not records:
    fail("package inventory must be a nonempty JSON array")
inventory = []
for index, record in enumerate(records, start=1):
    if not isinstance(record, dict):
        fail(f"package record {index} is not an object")
    build = record.get("build_string")
    if build is None:
        build = record.get("build")
    values = (
        record.get("name"),
        record.get("version"),
        build,
        record.get("channel"),
    )
    if any(
        not isinstance(value, str)
        or not value
        or any(character in value for character in "\t\r\n")
        for value in values
    ):
        fail(f"package record {index} has missing or noncanonical fields")
    inventory.append(values)
inventory.sort()
names = [record[0] for record in inventory]
if len(names) != len(set(names)):
    fail("package inventory contains duplicate package names")
with destination.open("w", encoding="utf-8", newline="\n") as stream:
    for record in inventory:
        stream.write("\t".join(record) + "\n")
"""
end

function _metamdbg_runtime_fsync_python()::String
    return raw"""
import os
import pathlib
import sys


path = pathlib.Path(sys.argv[1])
file_descriptor = os.open(path, os.O_RDONLY)
try:
    os.fsync(file_descriptor)
finally:
    os.close(file_descriptor)
directory_flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
directory_descriptor = os.open(path.parent, directory_flags)
try:
    os.fsync(directory_descriptor)
finally:
    os.close(directory_descriptor)
"""
end

function _metamdbg_runtime_directory_enumerator_python(
        ;
        post_open_hook::Union{Nothing, AbstractString} = nothing,
        post_enumeration_hook::Union{Nothing, AbstractString} = nothing,
        post_inventory_hook::Union{Nothing, AbstractString} = nothing,
        post_directory_open_hook::Union{Nothing, AbstractString} = nothing,
        post_directory_enumeration_hook::Union{Nothing, AbstractString} =
            nothing,
)::String
    function indent_hook(
            hook::Union{Nothing, AbstractString},
            spaces::Int = 4,
    )::String
        hook === nothing && return ""
        lines = split(String(hook), '\n'; keepempty = true)
        indentation = repeat(" ", spaces)
        return join((indentation * line for line in lines), '\n') * "\n"
    end
    open_hook = indent_hook(post_open_hook)
    enumeration_hook = indent_hook(post_enumeration_hook)
    inventory_hook = indent_hook(post_inventory_hook)
    directory_open_hook = indent_hook(post_directory_open_hook)
    directory_enumeration_hook = indent_hook(
        post_directory_enumeration_hook,
    )
    limit_declarations = """
MAX_DEPTH = $(_METAMDBG_ENUMERATION_MAX_DEPTH)
MAX_ENTRIES = $(_METAMDBG_ENUMERATION_MAX_ENTRIES)
MAX_MANIFEST_BYTES = $(_METAMDBG_ENUMERATION_MAX_MANIFEST_BYTES)


"""
    return raw"""
import os
import stat
import sys
import time

""" * limit_declarations * raw"""
def fail(message: str) -> None:
    raise RuntimeError(message)


def identity(status: os.stat_result) -> tuple[int, int, int]:
    return (status.st_dev, status.st_ino, stat.S_IFMT(status.st_mode))


def directory_names(directory_descriptor: int) -> tuple[bytes, ...]:
    names = []
    name_bytes = 0
    with os.scandir(directory_descriptor) as iterator:
        for entry in iterator:
            name = os.fsencode(entry.name)
            if not name or b"/" in name or b"\x00" in name:
                fail("directory enumeration returned an invalid entry name")
            if len(names) >= MAX_ENTRIES:
                fail(f"directory enumeration exceeds {MAX_ENTRIES} entries")
            name_bytes += len(name) + 1
            if name_bytes > MAX_MANIFEST_BYTES:
                fail(
                    "one directory name set exceeds the manifest byte bound "
                    f"of {MAX_MANIFEST_BYTES}"
                )
            names.append(name)
    names.sort()
    return tuple(names)


def post_directory_open(
    manifest_generation: int,
    relative_parts: tuple[bytes, ...],
) -> None:
""" * directory_open_hook * raw"""    return None


def post_directory_enumeration(
    manifest_generation: int,
    relative_parts: tuple[bytes, ...],
) -> None:
""" * directory_enumeration_hook * raw"""    return None


def account_manifest_bytes(state: dict[str, int], additional_bytes: int) -> None:
    state["bytes"] += additional_bytes
    if state["bytes"] > MAX_MANIFEST_BYTES:
        fail(f"directory manifest exceeds its byte bound of {MAX_MANIFEST_BYTES}")


def entry_status(
    directory_descriptor: int,
    name: bytes,
) -> os.stat_result:
    return os.stat(
        name,
        dir_fd=directory_descriptor,
        follow_symlinks=False,
    )


def build_manifest(
    root_descriptor: int,
    recursive: bool,
    manifest_generation: int,
    include_all_entries: bool,
    exact_names: tuple[bytes, ...],
    name_prefixes: tuple[bytes, ...],
) -> tuple[tuple[tuple, tuple], tuple[bytes, ...]]:
    entries = {}
    directories = {}
    state = {"entries": 0, "traversed": 0, "bytes": 0}

    def name_is_relevant(name: bytes) -> bool:
        return (
            include_all_entries
            or name in exact_names
            or any(name.startswith(prefix) for prefix in name_prefixes)
        )

    def account_traversed_path(relative_path: bytes) -> None:
        state["traversed"] += 1
        if state["traversed"] > MAX_ENTRIES:
            fail(f"directory traversal exceeds {MAX_ENTRIES} entries")
        account_manifest_bytes(state, len(relative_path) + 32)

    def record_entry(
        relative_path: bytes,
        status: os.stat_result,
    ) -> None:
        if relative_path in entries:
            fail("directory manifest contains a duplicate relative path")
        state["entries"] += 1
        if state["entries"] > MAX_ENTRIES:
            fail(f"directory manifest exceeds {MAX_ENTRIES} entries")
        account_manifest_bytes(state, len(relative_path) + 32)
        entries[relative_path] = identity(status)

    def frame(
        descriptor: int,
        relative_parts: tuple[bytes, ...],
        expected_status: os.stat_result,
        parent_descriptor: int | None,
        parent_name: bytes | None,
    ) -> dict[str, object]:
        if len(relative_parts) > MAX_DEPTH:
            fail(f"directory manifest exceeds depth bound {MAX_DEPTH}")
        descriptor_status = os.fstat(descriptor)
        if identity(descriptor_status) != identity(expected_status):
            fail("directory changed while its no-follow descriptor was opened")
        post_directory_open(manifest_generation, relative_parts)
        names = directory_names(descriptor)
        relative_path = b"/".join(relative_parts)
        if include_all_entries:
            if relative_path in directories:
                fail("directory manifest contains a duplicate directory path")
            directories[relative_path] = (identity(descriptor_status), names)
            account_manifest_bytes(
                state,
                len(relative_path) + sum(len(name) + 1 for name in names) + 32,
            )
        return {
            "descriptor": descriptor,
            "relative_parts": relative_parts,
            "status": descriptor_status,
            "names": names,
            "index": 0,
            "relevant_entries": {},
            "parent_descriptor": parent_descriptor,
            "parent_name": parent_name,
        }

    root_copy = os.dup(root_descriptor)
    stack = []
    try:
        try:
            root_frame = frame(
                root_copy,
                (),
                os.fstat(root_descriptor),
                None,
                None,
            )
        except BaseException:
            os.close(root_copy)
            raise
        stack.append(root_frame)
        while stack:
            current = stack[-1]
            names = current["names"]
            index = current["index"]
            if index < len(names):
                name = names[index]
                current["index"] = index + 1
                descriptor = current["descriptor"]
                relative_parts = current["relative_parts"] + (name,)
                if len(relative_parts) > MAX_DEPTH:
                    fail(f"directory manifest exceeds depth bound {MAX_DEPTH}")
                relative_path = b"/".join(relative_parts)
                if not include_all_entries:
                    account_traversed_path(relative_path)
                    if not recursive and not name_is_relevant(name):
                        continue
                status = entry_status(descriptor, name)
                if name_is_relevant(name):
                    record_entry(relative_path, status)
                    current["relevant_entries"][relative_path] = identity(status)
                if recursive and stat.S_ISDIR(status.st_mode):
                    child_descriptor = os.open(
                        name,
                        directory_flags,
                        dir_fd=descriptor,
                    )
                    try:
                        child_frame = frame(
                            child_descriptor,
                            relative_parts,
                            status,
                            descriptor,
                            name,
                        )
                    except BaseException:
                        os.close(child_descriptor)
                        raise
                    stack.append(child_frame)
                continue

            descriptor = current["descriptor"]
            relative_parts = current["relative_parts"]
            if include_all_entries:
                if directory_names(descriptor) != names:
                    fail("directory name set changed during manifest capture")
                for name in names:
                    relative_path = b"/".join(relative_parts + (name,))
                    if (
                        identity(entry_status(descriptor, name))
                        != entries[relative_path]
                    ):
                        fail("directory entry identity changed during manifest capture")
            else:
                relevant_entries = tuple(
                    sorted(
                        (
                            b"/".join(relative_parts + (name,)),
                            identity(entry_status(descriptor, name)),
                        )
                        for name in directory_names(descriptor)
                        if name_is_relevant(name)
                    )
                )
                captured_relevant_entries = tuple(
                    sorted(current["relevant_entries"].items())
                )
                if relevant_entries != captured_relevant_entries:
                    fail("reservation entry set changed during manifest capture")
            post_directory_enumeration(manifest_generation, relative_parts)
            if include_all_entries:
                if directory_names(descriptor) != names:
                    fail("directory name set changed during manifest hook")
                for name in names:
                    relative_path = b"/".join(relative_parts + (name,))
                    if (
                        identity(entry_status(descriptor, name))
                        != entries[relative_path]
                    ):
                        fail("directory entry identity changed during manifest hook")
            else:
                relevant_entries_after_hook = tuple(
                    sorted(
                        (
                            b"/".join(relative_parts + (name,)),
                            identity(entry_status(descriptor, name)),
                        )
                        for name in directory_names(descriptor)
                        if name_is_relevant(name)
                    )
                )
                if relevant_entries_after_hook != captured_relevant_entries:
                    fail("reservation entry set changed during manifest hook")
            parent_descriptor = current["parent_descriptor"]
            parent_name = current["parent_name"]
            if parent_descriptor is not None:
                observed_parent_identity = identity(
                    entry_status(parent_descriptor, parent_name)
                )
                if observed_parent_identity != identity(current["status"]):
                    fail("descendant directory path changed after traversal")
            os.close(descriptor)
            stack.pop()
    finally:
        for current in reversed(stack):
            try:
                os.close(current["descriptor"])
            except OSError:
                pass

    frozen_manifest = (
        tuple(sorted(entries.items())),
        tuple(sorted(directories.items())),
    )
    return frozen_manifest, tuple(sorted(entries))


def require_root_path_identity(
    root_path: bytes,
    expected_status: os.stat_result,
) -> None:
    path_status = os.stat(root_path, follow_symlinks=False)
    if identity(path_status) != identity(expected_status):
        fail("enumeration root path does not retain its opened directory inode")


root = os.fsencode(sys.argv[1])
inventory = os.fsencode(sys.argv[2])
mode = sys.argv[3]
for required_flag in ("O_DIRECTORY", "O_NOFOLLOW"):
    if not hasattr(os, required_flag):
        fail(f"platform lacks required directory flag {required_flag}")
directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
root_descriptor = os.open(root, directory_flags)
inventory_created = False
try:
    root_status = os.fstat(root_descriptor)
    if not stat.S_ISDIR(root_status.st_mode):
        fail("enumeration descriptor is not a directory")
""" * open_hook * raw"""    if mode == "one-level":
        recursive = False
        include_all_entries = True
        exact_names = ()
        reservation_prefixes = ()
    elif mode == "reservations-one-level":
        if len(sys.argv) != 7:
            fail(
                "one-level reservation enumeration requires an exact name "
                "and two prefixes"
            )
        recursive = False
        include_all_entries = False
        exact_names = (os.fsencode(sys.argv[4]),)
        reservation_prefixes = tuple(
            prefix
            for prefix in (
                os.fsencode(sys.argv[5]),
                os.fsencode(sys.argv[6]),
            )
            if prefix
        )
        if any(not name or b"/" in name for name in exact_names) or any(
            b"/" in prefix for prefix in reservation_prefixes
        ):
            fail("one-level reservation name or prefix is invalid")
    elif mode == "reservations-recursive":
        if len(sys.argv) not in (6, 7):
            fail("recursive reservation enumeration requires two prefixes")
        recursive = True
        include_all_entries = False
        exact_names = ()
        reservation_prefixes = (
            os.fsencode(sys.argv[4]),
            os.fsencode(sys.argv[5]),
        )
        if any(not prefix or b"/" in prefix for prefix in reservation_prefixes):
            fail("recursive reservation prefix is invalid")
    else:
        fail(f"unsupported directory enumeration mode: {mode}")
    first_manifest, first_entry_paths = build_manifest(
        root_descriptor,
        recursive,
        1,
        include_all_entries,
        exact_names,
        reservation_prefixes,
    )
""" * enumeration_hook * raw"""    require_root_path_identity(root, root_status)
    normalized_root = root.rstrip(b"/") or b"/"
    separator = b"" if normalized_root == b"/" else b"/"
    inventory_entries = first_entry_paths
    inventory_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW
    inventory_descriptor = os.open(inventory, inventory_flags, 0o600)
    inventory_created = True
    try:
        with os.fdopen(inventory_descriptor, "wb", closefd=True) as stream:
            for relative_path in inventory_entries:
                stream.write(normalized_root + separator + relative_path)
                stream.write(b"\x00")
    except BaseException:
        try:
            os.close(inventory_descriptor)
        except OSError:
            pass
        raise
""" * inventory_hook * raw"""    second_manifest, _second_entry_paths = build_manifest(
        root_descriptor,
        recursive,
        2,
        include_all_entries,
        exact_names,
        reservation_prefixes,
    )
    if second_manifest != first_manifest:
        fail("directory manifest changed before inventory publication")
    require_root_path_identity(root, root_status)
except BaseException:
    if inventory_created:
        try:
            os.unlink(inventory)
        except FileNotFoundError:
            pass
    raise
finally:
    os.close(root_descriptor)
"""
end

function _metamdbg_runtime_gfa_validator_python()::String
    return raw"""
import json
import math
import pathlib
import re
import struct
import sys


INTEGER_RE = re.compile(r"^[-+]?[0-9]+$")
FLOAT_RE = re.compile(r"^[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?$")
CIGAR_RE = re.compile(r"^(?:\*|(?:[0-9]+[MIDNSHPX=])+)$")
DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$")
B_RANGES = {
    "c": (-128, 127),
    "C": (0, 255),
    "s": (-32768, 32767),
    "S": (0, 65535),
    "i": (-2147483648, 2147483647),
    "I": (0, 4294967295),
}


def fail(message: str) -> None:
    raise ValueError(message)


def valid_name(name: str) -> bool:
    return (
        bool(name)
        and all(33 <= ord(character) <= 126 for character in name)
        and name[0] not in "*="
        and "+," not in name
        and "-," not in name
    )


def valid_float(value: str, *, single_precision: bool = False) -> bool:
    if not FLOAT_RE.fullmatch(value):
        return False
    try:
        parsed = float(value)
        if not math.isfinite(parsed):
            return False
        if single_precision:
            parsed = struct.unpack("!f", struct.pack("!f", parsed))[0]
            return math.isfinite(parsed)
        return True
    except (OverflowError, ValueError):
        return False


def valid_b_array(value: str) -> bool:
    components = value.split(",")
    if len(components) < 2 or len(components[0]) != 1:
        return False
    subtype = components[0]
    values = components[1:]
    if subtype == "f":
        return all(valid_float(item, single_precision=True) for item in values)
    if subtype not in B_RANGES:
        return False
    lower, upper = B_RANGES[subtype]
    for item in values:
        if not INTEGER_RE.fullmatch(item):
            return False
        try:
            parsed = int(item)
        except ValueError:
            return False
        if not lower <= parsed <= upper:
            return False
    return True


def valid_tag_value(tag_type: str, value: str) -> bool:
    if tag_type == "A":
        return len(value) == 1 and 33 <= ord(value) <= 126
    if tag_type == "i":
        if not INTEGER_RE.fullmatch(value):
            return False
        try:
            parsed = int(value)
        except ValueError:
            return False
        return -(2**63) <= parsed <= 2**63 - 1
    if tag_type == "f":
        return valid_float(value)
    if tag_type == "Z":
        return bool(value) and all(
            32 <= ord(character) <= 126 for character in value
        )
    if tag_type == "J":
        if not value or not all(32 <= ord(character) <= 126 for character in value):
            return False
        try:
            json.loads(
                value,
                parse_constant=lambda constant: fail(
                    f"nonstandard JSON constant: {constant}"
                ),
            )
        except (json.JSONDecodeError, ValueError):
            return False
        return True
    if tag_type == "H":
        return bool(value) and len(value) % 2 == 0 and bool(re.fullmatch(r"[0-9A-F]+", value))
    if tag_type == "B":
        return valid_b_array(value)
    return False


def validate_tags(fields: list[str], first_index: int) -> dict[str, tuple[str, str]]:
    tags = {}
    for tag in fields[first_index:]:
        components = tag.split(":", 2)
        if len(components) != 3:
            fail("invalid optional tag")
        name, tag_type, value = components
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]", name):
            fail("invalid optional tag name")
        if name in tags:
            fail("duplicate optional tag name")
        if not valid_tag_value(tag_type, value):
            fail("invalid optional tag value")
        tags[name] = (tag_type, value)
    return tags


def path_steps(field: str) -> list[str]:
    if not field:
        fail("empty path step list")
    identifiers = []
    step_start = 0
    position = 0
    field_length = len(field)
    while position < field_length:
        character = field[position]
        ends_step = character in "+-" and (
            position + 1 == field_length or field[position + 1] == ","
        )
        if not ends_step:
            position += 1
            continue
        if position == step_start:
            fail("empty path step identifier")
        identifier = field[step_start:position]
        if not valid_name(identifier):
            fail("invalid path step identifier")
        identifiers.append(identifier)
        if position + 1 == field_length:
            step_start = field_length
            break
        step_start = position + 2
        if step_start >= field_length:
            fail("trailing path step separator")
        position = step_start
    if step_start != field_length:
        fail("unterminated path step")
    return identifiers


def records(path: pathlib.Path):
    with path.open("r", encoding="utf-8", newline=None) as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            yield line_number, line.split("\t")


path = pathlib.Path(sys.argv[1])
segments = set()
record_names = set()
for line_number, fields in records(path):
    record_type = fields[0]
    if record_type == "H":
        tags = validate_tags(fields, 1)
        if "VN" in tags:
            tag_type, value = tags["VN"]
            if tag_type != "Z" or not re.fullmatch(r"1\.[0-9]+", value):
                fail(f"non-GFA1 VN header at line {line_number}")
    elif record_type == "S":
        if len(fields) < 3 or not valid_name(fields[1]):
            fail(f"invalid segment at line {line_number}")
        identifier = fields[1]
        if identifier in record_names:
            fail(f"duplicate segment/path name at line {line_number}")
        if not DNA_RE.fullmatch(fields[2]) or fields[2] == "*":
            fail(f"invalid segment sequence at line {line_number}")
        validate_tags(fields, 3)
        segments.add(identifier)
        record_names.add(identifier)
    elif record_type == "L":
        if len(fields) < 6:
            fail(f"malformed link at line {line_number}")
        if not valid_name(fields[1]) or not valid_name(fields[3]):
            fail(f"invalid link identifier at line {line_number}")
        if fields[2] not in ("+", "-") or fields[4] not in ("+", "-"):
            fail(f"invalid link orientation at line {line_number}")
        if not CIGAR_RE.fullmatch(fields[5]):
            fail(f"invalid link CIGAR at line {line_number}")
        validate_tags(fields, 6)
    elif record_type == "P":
        if len(fields) < 4 or not valid_name(fields[1]):
            fail(f"invalid path at line {line_number}")
        identifier = fields[1]
        if identifier in record_names:
            fail(f"duplicate segment/path name at line {line_number}")
        steps = path_steps(fields[2])
        if fields[3] != "*":
            overlaps = fields[3].split(",")
            if len(overlaps) != len(steps) - 1:
                fail(f"wrong path overlap count at line {line_number}")
            if any(not CIGAR_RE.fullmatch(overlap) for overlap in overlaps):
                fail(f"invalid path overlap at line {line_number}")
        validate_tags(fields, 4)
        record_names.add(identifier)
    else:
        fail(f"unknown GFA record type at line {line_number}")
if not segments:
    fail("no sequence-bearing segments")
for line_number, fields in records(path):
    if fields[0] == "L":
        references = (fields[1], fields[3])
    elif fields[0] == "P":
        references = path_steps(fields[2])
    else:
        continue
    if any(identifier not in segments for identifier in references):
        fail(f"dangling segment reference at line {line_number}")
"""
end

function _metamdbg_environment_packages(;
        conda_runner::AbstractString = _conda_runner(),
        command_reader::Function = command -> read(command, String),
)::Vector{NamedTuple}
    command = Cmd(String[
        String(conda_runner),
        "list",
        "-n",
        METAMDBG_ENV_NAME,
        "--json",
    ])
    package_records = JSON.parse(command_reader(command))
    package_records isa AbstractVector || error(
        "metaMDBG conda package inventory was not a JSON array.",
    )
    return _metamdbg_normalized_package_inventory(package_records)
end

function _require_metamdbg_package_version(
        package_records::AbstractVector,
)::NamedTuple
    inventory = _metamdbg_normalized_package_inventory(package_records)
    metamdbg_records = filter(
        record -> record.name == "metamdbg",
        inventory,
    )
    actual_version = isempty(metamdbg_records) ?
                     nothing : only(metamdbg_records).version
    displayed_version =
        actual_version === nothing ? "missing" : repr(actual_version)
    actual_version == METAMDBG_VERSION || error(
        "metaMDBG spec-addressed environment $(repr(METAMDBG_ENV_NAME)) must " *
        "contain metamdbg exactly $(METAMDBG_VERSION), got " *
        "$(displayed_version). " *
        "Refusing to repair it in place; remove it manually before reinstalling.",
    )
    expected = _metamdbg_expected_toolchain()
    return (;
        expected...,
        package_inventory = inventory,
        package_inventory_sha256 =
            _metamdbg_package_inventory_sha256(inventory),
        package_count = length(inventory),
    )
end

function _require_expected_metamdbg_toolchain(toolchain::Any)::NamedTuple
    expected = _metamdbg_expected_toolchain()
    toolchain isa NamedTuple || error(
        "metaMDBG dependency validation did not return toolchain provenance.",
    )
    for field in propertynames(expected)
        actual_value = hasproperty(toolchain, field) ?
                       getproperty(toolchain, field) : nothing
        if actual_value != getproperty(expected, field)
            error(
                "metaMDBG dependency validation returned incompatible " *
                "$(field) provenance: expected " *
                "$(repr(getproperty(expected, field))), got " *
                "$(repr(actual_value)).",
            )
        end
    end
    hasproperty(toolchain, :package_inventory) || error(
        "metaMDBG dependency validation did not return its resolved package " *
        "inventory.",
    )
    realized = _require_metamdbg_package_version(toolchain.package_inventory)
    if !hasproperty(toolchain, :package_inventory_sha256) ||
       toolchain.package_inventory_sha256 !=
       realized.package_inventory_sha256
        error(
            "metaMDBG dependency validation returned an incompatible " *
            "normalized package inventory digest.",
        )
    end
    if !hasproperty(toolchain, :package_count) ||
       toolchain.package_count != realized.package_count
        error(
            "metaMDBG dependency validation returned an incompatible package " *
            "count.",
        )
    end
    propertynames(toolchain) == propertynames(realized) || error(
        "metaMDBG dependency validation returned unexpected toolchain fields.",
    )
    return realized
end

function _require_unchanged_metamdbg_toolchain(
        before::NamedTuple,
        after::NamedTuple,
)::NamedTuple
    normalized_before = _require_expected_metamdbg_toolchain(before)
    normalized_after = _require_expected_metamdbg_toolchain(after)
    normalized_after == normalized_before || error(
        "metaMDBG resolved package inventory changed while assembly tools " *
        "were running. Refusing to publish mixed-toolchain provenance.",
    )
    return normalized_after
end

function _canonical_metamdbg_conda_runner(
        conda_runner::AbstractString,
)::String
    executable = Sys.which(String(conda_runner))
    candidate = executable === nothing ?
                abspath(String(conda_runner)) : String(executable)
    return ispath(candidate) ? realpath(candidate) : normpath(candidate)
end

function _metamdbg_environment_prefix(
        conda_runner::AbstractString,
)::String
    canonical_runner = _canonical_metamdbg_conda_runner(conda_runner)
    conda_root = normpath(joinpath(dirname(canonical_runner), ".."))
    return joinpath(conda_root, "envs", METAMDBG_ENV_NAME)
end

function _metamdbg_install_lock_path(
        conda_runner::AbstractString = _conda_runner(),
)::String
    environment_prefix = _metamdbg_environment_prefix(conda_runner)
    conda_root = dirname(dirname(environment_prefix))
    return joinpath(
        conda_root,
        ".mycelia-locks",
        "$(METAMDBG_ENV_NAME).pid",
    )
end

function _metamdbg_environment_is_installed(
        conda_runner::AbstractString,
)::Bool
    environment_prefix = _metamdbg_environment_prefix(conda_runner)
    environment_names = _conda_environment_names(
        _canonical_metamdbg_conda_runner(conda_runner),
        dirname(environment_prefix),
    )
    return METAMDBG_ENV_NAME in environment_names
end

function _create_metamdbg_environment(
        environment_path::AbstractString,
        environment_name::AbstractString,
        conda_runner::AbstractString;
        force::Bool = false,
)::String
    force && error(
        "metaMDBG immutable environments cannot be recreated in place.",
    )
    resolved_runner = _canonical_metamdbg_conda_runner(conda_runner)
    verified_environment_path =
        _require_verified_metamdbg_environment_spec(environment_path)
    Base.run(Cmd(String[
        resolved_runner,
        "env",
        "create",
        "-f",
        verified_environment_path,
        "-n",
        String(environment_name),
    ]))
    Base.run(Cmd(String[resolved_runner, "clean", "--all", "-y"]))
    return String(environment_name)
end

function _with_metamdbg_install_lock(
        action::Function,
        lock_path::AbstractString,
)::Any
    normalized_lock_path = abspath(lock_path)
    mkpath(dirname(normalized_lock_path))
    return FileWatching.Pidfile.mkpidlock(
        normalized_lock_path;
        stale_age = _METAMDBG_INSTALL_LOCK_STALE_SECONDS,
        refresh = _METAMDBG_INSTALL_LOCK_REFRESH_SECONDS,
    ) do
        action()
    end
end

function _ensure_metamdbg_installed_locked(;
        conda_runner::AbstractString = _conda_runner(),
        paths::Tuple{String, String} = _metamdbg_paths(),
        environment_checker::Function = _metamdbg_environment_is_installed,
        environment_creator::Function = _create_metamdbg_environment,
        package_inspector::Function = runner ->
            _metamdbg_environment_packages(conda_runner = runner),
)::NamedTuple
    resolved_runner = _canonical_metamdbg_conda_runner(conda_runner)
    install_dir, environment_path = paths
    mkpath(install_dir)
    verified_environment_path =
        _require_verified_metamdbg_environment_spec(environment_path)
    if !Bool(environment_checker(resolved_runner))
        environment_creator(
            verified_environment_path,
            METAMDBG_ENV_NAME,
            resolved_runner;
            force = false,
        )
    end
    return _require_metamdbg_package_version(package_inspector(resolved_runner))
end

function _ensure_metamdbg_installed(;
        conda_runner::AbstractString = _conda_runner(),
        paths::Tuple{String, String} = _metamdbg_paths(),
        environment_checker::Function = _metamdbg_environment_is_installed,
        environment_creator::Function = _create_metamdbg_environment,
        package_inspector::Function = runner ->
            _metamdbg_environment_packages(conda_runner = runner),
        lock_path::Union{Nothing, AbstractString} = nothing,
        lock_runner::Function = _with_metamdbg_install_lock,
)::NamedTuple
    resolved_runner = _canonical_metamdbg_conda_runner(conda_runner)
    resolved_lock_path = lock_path === nothing ?
                         _metamdbg_install_lock_path(resolved_runner) :
                         String(lock_path)
    return lock_runner(resolved_lock_path) do
        _ensure_metamdbg_installed_locked(;
            conda_runner = resolved_runner,
            paths,
            environment_checker,
            environment_creator,
            package_inspector,
        )
    end
end

function _metamdbg_output_paths(
        outdir::AbstractString,
        graph_k::Int,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    return (
        outdir = normalized_outdir,
        contigs_plain = joinpath(normalized_outdir, "contigs.fasta"),
        contigs_gz = joinpath(normalized_outdir, "contigs.fasta.gz"),
        graph_alias = joinpath(
            normalized_outdir,
            "assemblyGraph_k$(graph_k).gfa",
        ),
        contract_marker = joinpath(
            normalized_outdir,
            _METAMDBG_CONTRACT_FILENAME,
        ),
        completion_marker = joinpath(
            normalized_outdir,
            _METAMDBG_COMPLETION_FILENAME,
        ),
    )
end

function _require_positive_metamdbg_graph_k(graph_k::Int)::Int
    graph_k > 0 || throw(ArgumentError("graph_k must be positive."))
    return graph_k
end

function _metamdbg_canonical_output_path(outdir::AbstractString)::String
    stripped_outdir = strip(String(outdir))
    isempty(stripped_outdir) && throw(ArgumentError("outdir must be nonempty."))
    requested_path = normpath(abspath(stripped_outdir))
    islink(requested_path) && throw(ArgumentError(
        "metaMDBG outdir must not be a symbolic link: $(requested_path)",
    ))

    missing_components = String[]
    existing_ancestor = requested_path
    while !ispath(existing_ancestor)
        islink(existing_ancestor) && throw(ArgumentError(
            "metaMDBG outdir has a dangling symbolic-link ancestor: " *
            "$(existing_ancestor)",
        ))
        parent = dirname(existing_ancestor)
        parent == existing_ancestor && error(
            "Unable to find an existing ancestor for metaMDBG outdir: " *
            "$(requested_path)",
        )
        push!(missing_components, basename(existing_ancestor))
        existing_ancestor = parent
    end
    isdir(existing_ancestor) || throw(ArgumentError(
        "metaMDBG outdir has a non-directory ancestor: $(existing_ancestor)",
    ))

    canonical_path = realpath(existing_ancestor)
    for component in Iterators.reverse(missing_components)
        canonical_path = joinpath(canonical_path, component)
    end
    canonical_path = normpath(canonical_path)
    if ispath(canonical_path) && !isdir(canonical_path)
        throw(ArgumentError(
            "metaMDBG outdir exists but is not a directory: $(canonical_path)",
        ))
    end
    islink(canonical_path) && throw(ArgumentError(
        "metaMDBG outdir must not be a symbolic link: $(canonical_path)",
    ))
    return canonical_path
end

function _metamdbg_output_root_identity(
        outdir::AbstractString,
)::NamedTuple
    normalized_outdir = normpath(abspath(String(outdir)))
    canonical_outdir = _metamdbg_canonical_output_path(normalized_outdir)
    canonical_outdir == normalized_outdir || error(
        "metaMDBG output root is not bound to its canonical physical " *
        "location: $(normalized_outdir) resolves to $(canonical_outdir).",
    )
    isdir(canonical_outdir) && !islink(canonical_outdir) || error(
        "metaMDBG output root must be a regular directory before binding its " *
        "physical identity: $(canonical_outdir).",
    )
    metadata = stat(canonical_outdir)
    return (;
        canonical_outdir,
        device = metadata.device,
        inode = metadata.inode,
    )
end

function _require_unchanged_metamdbg_output_root(
        identity::NamedTuple,
        outdir::AbstractString,
)::String
    expected_outdir = String(identity.canonical_outdir)
    normalized_outdir = normpath(abspath(String(outdir)))
    normalized_outdir == expected_outdir || error(
        "metaMDBG output root path changed after its physical identity was " *
        "bound: expected $(expected_outdir), observed $(normalized_outdir).",
    )
    observed_outdir = try
        _metamdbg_canonical_output_path(normalized_outdir)
    catch caught
        caught isa InterruptException && rethrow()
        error(
            "metaMDBG output root physical identity changed after it was " *
            "bound: $(normalized_outdir). Cause: $(sprint(showerror, caught))",
        )
    end
    observed_outdir == expected_outdir || error(
        "metaMDBG output root canonical physical location changed after it " *
        "was bound: expected $(expected_outdir), observed " *
        "$(observed_outdir).",
    )
    isdir(observed_outdir) && !islink(observed_outdir) || error(
        "metaMDBG output root is no longer a regular directory: " *
        "$(observed_outdir).",
    )
    metadata = stat(observed_outdir)
    if metadata.device != identity.device || metadata.inode != identity.inode
        error(
            "metaMDBG output root physical identity changed after it was " *
            "bound: $(observed_outdir).",
        )
    end
    return observed_outdir
end

function _with_metamdbg_output_root_guard(
        action::Function,
        output_root_identity::NamedTuple,
        operation::AbstractString,
)::Any
    _require_unchanged_metamdbg_output_root(
        output_root_identity,
        output_root_identity.canonical_outdir,
    )
    result = try
        action()
    catch primary_error
        try
            _require_unchanged_metamdbg_output_root(
                output_root_identity,
                output_root_identity.canonical_outdir,
            )
        catch identity_error
            identity_error isa InterruptException && rethrow()
            error(
                "metaMDBG output root changed while $(operation) failed. " *
                "Operation cause: $(sprint(showerror, primary_error)). " *
                "Identity cause: $(sprint(showerror, identity_error))",
            )
        end
        Base.rethrow()
    end
    _require_unchanged_metamdbg_output_root(
        output_root_identity,
        output_root_identity.canonical_outdir,
    )
    return result
end

function _require_metamdbg_artifact_containment(
        path::AbstractString,
        label::AbstractString,
        output_root_identity::NamedTuple,
)::String
    canonical_outdir = _require_unchanged_metamdbg_output_root(
        output_root_identity,
        output_root_identity.canonical_outdir,
    )
    normalized_path = normpath(abspath(String(path)))
    canonical_path = try
        realpath(normalized_path)
    catch caught
        caught isa InterruptException && rethrow()
        error(
            "$(label) cannot be resolved within the bound metaMDBG output " *
            "root: $(normalized_path). Cause: $(sprint(showerror, caught))",
        )
    end
    dirname(canonical_path) == canonical_outdir || error(
        "$(label) escaped the bound metaMDBG output root: " *
        "$(canonical_path) is not directly contained by " *
        "$(canonical_outdir).",
    )
    return canonical_path
end

function _metamdbg_output_lock_path(outdir::AbstractString)::String
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    output_identity = _output_root_reservation_identity(normalized_outdir)
    return joinpath(
        dirname(normalized_outdir),
        ".mycelia-metamdbg.lock.$(output_identity)",
    )
end

function _metamdbg_output_lock_identity(
        lock_path::AbstractString,
)::NamedTuple
    normalized_lock_path = normpath(abspath(String(lock_path)))
    isdir(normalized_lock_path) && !islink(normalized_lock_path) || error(
        "metaMDBG lifecycle lock is missing, replaced, or not a regular " *
        "directory: $(normalized_lock_path).",
    )
    metadata = stat(normalized_lock_path)
    return (; device = metadata.device, inode = metadata.inode)
end

function _require_unchanged_metamdbg_output_lock(
        lock_path::AbstractString,
        expected_identity::NamedTuple,
)::String
    normalized_lock_path = normpath(abspath(String(lock_path)))
    observed_identity = _metamdbg_output_lock_identity(normalized_lock_path)
    observed_identity == expected_identity || error(
        "metaMDBG lifecycle lock was replaced before cleanup: " *
        "$(normalized_lock_path).",
    )
    return normalized_lock_path
end

function _with_metamdbg_output_lock(
        action::Function,
        outdir::AbstractString,
        ;
        lock_remover::Union{Nothing, Function} = nothing,
        lock_cleanup_failure_handler::Function =
            (_path::AbstractString, _cleanup_error::Any) -> nothing,
)::Any
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    lock_path = _metamdbg_output_lock_path(canonical_outdir)
    mkpath(dirname(lock_path))
    try
        Base.Filesystem.mkdir(lock_path; mode = 0o700)
    catch caught
        caught isa InterruptException && rethrow()
        error(
            "metaMDBG output is locked by another lifecycle: $(lock_path). " *
            "If no workflow is active, remove this stale lock directory. Cause: " *
            sprint(showerror, caught),
        )
    end
    lock_identity = _metamdbg_output_lock_identity(lock_path)
    cleanup_lock = function ()
        _require_unchanged_metamdbg_output_lock(lock_path, lock_identity)
        # Retained only as an internal fault-injection hook. Production
        # cleanup never delegates pathname removal to this callback.
        lock_remover === nothing || lock_remover(lock_path)
        _remove_exact_metamdbg_durable_directory!(
            lock_path,
            merge(lock_identity, (; path = lock_path)),
        )
        return nothing
    end
    result = try
        _metamdbg_canonical_output_path(canonical_outdir) == canonical_outdir ||
            error("metaMDBG output path changed while acquiring its lock.")
        action()
    catch primary_error
        try
            cleanup_lock()
        catch cleanup_error
            try
                lock_cleanup_failure_handler(lock_path, cleanup_error)
            catch preservation_error
                @error "metaMDBG could not preserve a shared fail-closed " *
                       "reservation after private lifecycle cleanup failed" lock_path primary_error cleanup_error preservation_error
            end
            @warn "metaMDBG failed to release its output lock while preserving " *
                  "the primary workflow failure" lock_path primary_error cleanup_error
        end
        Base.rethrow()
    end
    try
        cleanup_lock()
    catch cleanup_error
        try
            lock_cleanup_failure_handler(lock_path, cleanup_error)
        catch preservation_error
            error(
                "metaMDBG private lifecycle cleanup failed and its shared " *
                "fail-closed reservation could not be published. Private " *
                "cause: $(sprint(showerror, cleanup_error)). Preservation " *
                "cause: $(sprint(showerror, preservation_error))",
            )
        end
        throw(cleanup_error)
    end
    return result
end

function _metamdbg_lifecycle_cleanup_reservation_path(
        outdir::AbstractString,
)::String
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    return _output_root_durable_reservation_path_from_canonical(
        canonical_outdir,
        "metamdbg-private-lifecycle-cleanup-failed",
    )
end

function _publish_metamdbg_lifecycle_cleanup_reservation!(
        outdir::AbstractString,
)::String
    reservation_path =
        _metamdbg_lifecycle_cleanup_reservation_path(outdir)
    !_output_root_path_entry_exists(reservation_path) || error(
        "metaMDBG refuses to replace an existing lifecycle cleanup " *
        "reservation: $(reservation_path).",
    )
    Base.Filesystem.mkdir(reservation_path; mode = 0o700)
    chmod(reservation_path, 0o700)
    reservation_status = stat(reservation_path)
    reservation_status.uid == Base.Libc.getuid() || error(
        "metaMDBG lifecycle cleanup reservation is not owned by the current " *
        "user: $(reservation_path).",
    )
    (reservation_status.mode & 0o777) == 0o700 || error(
        "metaMDBG lifecycle cleanup reservation must have mode 0700: " *
        "$(reservation_path).",
    )
    _fsync_metamdbg_directory(dirname(reservation_path))
    return reservation_path
end

function _require_unchanged_metamdbg_lifecycle_cleanup_reservation(
        reservation_path::AbstractString,
        expected_identity::NamedTuple,
)::String
    normalized_path = normpath(abspath(String(reservation_path)))
    isdir(normalized_path) && !islink(normalized_path) || error(
        "metaMDBG lifecycle cleanup reservation is missing, replaced, or not " *
        "a regular directory: $(normalized_path).",
    )
    reservation_status = stat(normalized_path)
    observed_identity = (;
        device = reservation_status.device,
        inode = reservation_status.inode,
    )
    observed_identity == expected_identity || error(
        "metaMDBG lifecycle cleanup reservation was replaced: " *
        "$(normalized_path).",
    )
    reservation_status.uid == Base.Libc.getuid() || error(
        "metaMDBG lifecycle cleanup reservation changed owner: " *
        "$(normalized_path).",
    )
    (reservation_status.mode & 0o777) == 0o700 || error(
        "metaMDBG lifecycle cleanup reservation changed mode: " *
        "$(normalized_path).",
    )
    return normalized_path
end

function _remove_metamdbg_lifecycle_cleanup_reservation!(
        reservation_path::AbstractString,
        expected_identity::NamedTuple,
)::Nothing
    normalized_path =
        _require_unchanged_metamdbg_lifecycle_cleanup_reservation(
            reservation_path,
            expected_identity,
        )
    isempty(readdir(normalized_path)) || error(
        "metaMDBG lifecycle cleanup reservation contains unexpected entries: " *
        "$(normalized_path).",
    )
    _remove_exact_metamdbg_durable_directory!(
        normalized_path,
        merge(expected_identity, (; path = normalized_path)),
    )
    return nothing
end

function _metamdbg_output_root_pid_lock_identity(
        lock_path::AbstractString,
)::NamedTuple
    normalized_path = normpath(abspath(String(lock_path)))
    return _metamdbg_output_root_pid_lock_snapshot(normalized_path).identity
end

function _metamdbg_output_root_pid_lock_snapshot(
        lock_path::AbstractString;
        post_parent_descriptor_open_hook::Function =
            (_path::AbstractString, _identity::NamedTuple) -> nothing,
)::NamedTuple
    normalized_path = normpath(abspath(String(lock_path)))
    parent_identity =
        _metamdbg_directory_path_identity(dirname(normalized_path))
    parent = Base.Filesystem.open(
        parent_identity.path,
        Base.JL_O_RDONLY |
        Base.JL_O_DIRECTORY |
        Base.JL_O_NOFOLLOW |
        Base.JL_O_CLOEXEC,
    )
    lock = nothing
    try
        _require_unchanged_metamdbg_directory_descriptor(
            parent,
            parent_identity,
            parent_identity.path,
        )
        post_parent_descriptor_open_hook(
            parent_identity.path,
            parent_identity,
        )
        _metamdbg_directory_path_identity(parent_identity.path) ==
            parent_identity || error(
            "metaMDBG output-root PID lock parent changed after its descriptor " *
            "was opened: $(parent_identity.path).",
        )
        lock = _metamdbg_openat(
            parent,
            basename(normalized_path),
            Base.JL_O_RDONLY |
            Base.JL_O_NONBLOCK |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
        )
        lock_status = stat(lock)
        isfile(lock_status) || error(
            "metaMDBG output-root PID lock is missing, replaced, or not a " *
            "regular file: $(normalized_path).",
        )
        lock_status.uid == Base.Libc.getuid() || error(
            "metaMDBG output-root PID lock is not owned by the current user: " *
            "$(normalized_path).",
        )
        identity = (;
            device = lock_status.device,
            inode = lock_status.inode,
        )
        contents = read(lock, String)
        observed_status = stat(lock)
        observed_status.device == identity.device &&
            observed_status.inode == identity.inode || error(
            "metaMDBG output-root PID lock descriptor changed while being " *
            "read: $(normalized_path).",
        )
        observed_identity = _metamdbg_openat_entry_identity(
            parent,
            basename(normalized_path),
            normalized_path,
            false,
        )
        _metamdbg_path_identity_matches(
            observed_identity,
            merge(identity, (; path = normalized_path)),
        ) || error(
            "metaMDBG output-root PID lock was replaced while being read: " *
            "$(normalized_path).",
        )
        _require_unchanged_metamdbg_directory_descriptor(
            parent,
            parent_identity,
            parent_identity.path,
        )
        _metamdbg_directory_path_identity(parent_identity.path) ==
            parent_identity || error(
            "metaMDBG output-root PID lock parent changed while its snapshot " *
            "was read: $(parent_identity.path).",
        )
        return (; identity, contents)
    finally
        if lock !== nothing && isopen(lock)
            close(lock)
        end
        close(parent)
    end
end

function _remove_dead_metamdbg_output_root_pid_lock!(
        lock_path::AbstractString;
        pre_quarantine_unlink_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
)::Nothing
    normalized_path = normpath(abspath(String(lock_path)))
    expected_snapshot =
        _metamdbg_output_root_pid_lock_snapshot(normalized_path)
    expected_identity = expected_snapshot.identity
    lock_contents = expected_snapshot.contents
    lock_match = match(r"^([1-9][0-9]*) ([^[:space:]]+)$", lock_contents)
    lock_match === nothing && error(
        "metaMDBG refuses dead-process recovery because its output-root PID " *
        "lock is empty or malformed: $(normalized_path).",
    )
    parsed_pid = Base.tryparse(UInt64, lock_match.captures[1])
    parsed_pid !== nothing && parsed_pid <= UInt64(Base.typemax(Cint)) || error(
        "metaMDBG refuses dead-process recovery because its output-root PID " *
        "lock does not contain a locally verifiable process ID: " *
        "$(normalized_path).",
    )
    hostname = lock_match.captures[2]
    hostname == Base.gethostname() || error(
        "metaMDBG refuses dead-process recovery because its output-root PID " *
        "lock names a remote or unverifiable host: $(normalized_path).",
    )
    pid = Cuint(parsed_pid)
    expected_contents = "$(pid) $(hostname)"
    lock_contents == expected_contents || error(
        "metaMDBG refuses dead-process recovery because its output-root PID " *
        "lock is not canonical: $(normalized_path).",
    )
    !FileWatching.Pidfile.isvalidpid(hostname, pid) || error(
        "metaMDBG refuses dead-process recovery because its output-root PID " *
        "lock still names a live or remotely unverifiable process: " *
        "$(normalized_path).",
    )
    observed_snapshot =
        _metamdbg_output_root_pid_lock_snapshot(normalized_path)
    observed_snapshot.identity == expected_identity || error(
        "metaMDBG output-root PID lock was replaced during dead-process " *
        "recovery: $(normalized_path).",
    )
    observed_snapshot.contents == expected_contents || error(
        "metaMDBG output-root PID lock contents changed during dead-process " *
        "recovery: $(normalized_path).",
    )
    !FileWatching.Pidfile.isvalidpid(hostname, pid) || error(
        "metaMDBG output-root PID became live again during dead-process " *
        "recovery: $(normalized_path).",
    )
    function quarantined_payload_validator(
            payload::Base.Filesystem.File,
            payload_path::AbstractString,
    )::Nothing
        seekstart(payload)
        first_contents = read(payload, String)
        seekstart(payload)
        second_contents = read(payload, String)
        first_contents == expected_contents &&
            second_contents == expected_contents || error(
            "metaMDBG quarantined output-root PID lock contents changed " *
            "during dead-process recovery: $(payload_path).",
        )
        !FileWatching.Pidfile.isvalidpid(hostname, pid) || error(
            "metaMDBG quarantined output-root PID became live during " *
            "dead-process recovery: $(payload_path).",
        )
        return nothing
    end
    _remove_exact_metamdbg_durable_file!(
        normalized_path,
        merge(expected_identity, (; path = normalized_path));
        pre_quarantine_unlink_hook,
        quarantined_payload_validator,
    )
    return nothing
end

function _recover_dead_metamdbg_private_output_lock!(
        lock_path::AbstractString,
)::Nothing
    normalized_path = normpath(abspath(String(lock_path)))
    expected_identity = _metamdbg_output_lock_identity(normalized_path)
    lock_status = stat(normalized_path)
    lock_status.uid == Base.Libc.getuid() || error(
        "metaMDBG private lifecycle lock is not owned by the current user: " *
        "$(normalized_path).",
    )
    (lock_status.mode & 0o777) == 0o700 || error(
        "metaMDBG private lifecycle lock must have mode 0700 for explicit " *
        "dead-process recovery: $(normalized_path).",
    )
    isempty(readdir(normalized_path)) || error(
        "metaMDBG private lifecycle lock contains unexpected entries: " *
        "$(normalized_path).",
    )
    _require_unchanged_metamdbg_output_lock(
        normalized_path,
        expected_identity,
    )
    _remove_exact_metamdbg_durable_directory!(
        normalized_path,
        merge(expected_identity, (; path = normalized_path)),
    )
    return nothing
end

function _recover_dead_metamdbg_lifecycle_locks!(
        outdir::AbstractString;
        pending_recovery_function::Function =
            _recover_metamdbg_pending_submission_job_records!,
        post_recovery_pid_acquisition_hook::Function =
            (_outdir::AbstractString) -> nothing,
)::Nothing
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    reservation_lock_path =
        _output_root_reservation_lock_path_from_canonical(canonical_outdir)
    private_lock_path = _metamdbg_output_lock_path(canonical_outdir)
    cleanup_reservation_path =
        _metamdbg_lifecycle_cleanup_reservation_path(canonical_outdir)
    runtime_owner_exists = any(
        path -> _metamdbg_submission_reservation_path_state(
            path,
            canonical_outdir,
        ) == :runtime,
        _metamdbg_submission_reservation_paths(canonical_outdir),
    )
    runtime_marker_exists =
        _metamdbg_has_runtime_scheduler_evidence(canonical_outdir)
    (runtime_owner_exists || runtime_marker_exists) && error(
        "metaMDBG refuses dead-submitter lock recovery while durable runtime " *
        "scheduler ownership exists. Use exact job-ID-bound terminal or " *
        "cancelled recovery instead.",
    )
    had_pid_lock = _output_root_path_entry_exists(reservation_lock_path)
    has_private_state =
        _output_root_path_entry_exists(private_lock_path) ||
        _output_root_path_entry_exists(cleanup_reservation_path)
    if !had_pid_lock
        has_private_state && error(
            "metaMDBG refuses dead-process recovery because private " *
            "lifecycle state exists without its pre-existing local PID " *
            "ownership record. The state may belong to a scheduler runtime; " *
            "use explicit job-ID recovery evidence.",
        )
        error(
            "metaMDBG confirm_process_dead requires a pre-existing canonical " *
            "local PID ownership record; refusing a no-op recovery.",
        )
    end
    _remove_dead_metamdbg_output_root_pid_lock!(reservation_lock_path)
    stale_age = _OUTPUT_ROOT_RESERVATION_STALE_AGE_SECONDS
    lock_handle = FileWatching.Pidfile.trymkpidlock(
        reservation_lock_path;
        stale_age,
        refresh = stale_age / 2,
    )
    lock_handle === false && error(
        "metaMDBG could not acquire its output-root PID lock after explicit " *
        "dead-process recovery: $(reservation_lock_path).",
    )
    try
        post_recovery_pid_acquisition_hook(canonical_outdir)
        post_acquisition_paths = _metamdbg_submission_reservation_paths(
            canonical_outdir;
            include_consumed = true,
        )
        post_acquisition_runtime_owner = any(post_acquisition_paths) do path
            return _metamdbg_submission_reservation_path_state(
                path,
                canonical_outdir,
            ) == :runtime
        end
        post_acquisition_runtime_shared =
            _metamdbg_has_runtime_scheduler_evidence(canonical_outdir)
        (post_acquisition_runtime_owner || post_acquisition_runtime_shared) &&
            error(
                "metaMDBG runtime scheduler ownership appeared after the " *
                "replacement recovery PID was acquired; refusing destructive " *
                "private, sentinel, or pending cleanup.",
            )
        if _output_root_path_entry_exists(private_lock_path)
            _recover_dead_metamdbg_private_output_lock!(private_lock_path)
        end
        if _output_root_path_entry_exists(cleanup_reservation_path)
            cleanup_identity =
                _metamdbg_output_lock_identity(cleanup_reservation_path)
            _remove_metamdbg_lifecycle_cleanup_reservation!(
                cleanup_reservation_path,
                cleanup_identity,
            )
        end
        pending_recovery_function(canonical_outdir)
    finally
        Base.close(lock_handle)
    end
    return nothing
end

function _metamdbg_has_runtime_scheduler_evidence(
        outdir::AbstractString,
)::Bool
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    reservation_paths = _metamdbg_submission_reservation_paths(
        canonical_outdir;
        include_consumed = true,
    )
    known_queued_markers = Set{String}()
    for path in reservation_paths
        path_state = _metamdbg_submission_reservation_path_state(
            path,
            canonical_outdir,
        )
        path_state == :runtime && return true
        reservation = _metamdbg_submission_reservation_from_path(
            path,
            canonical_outdir;
            allow_provisional = true,
        )
        push!(
            known_queued_markers,
            reservation.output_root_reservation_marker,
        )
        _output_root_path_entry_exists(
            reservation.runtime_output_root_reservation_marker,
        ) && return true
    end
    cleanup_path =
        _metamdbg_lifecycle_cleanup_reservation_path(canonical_outdir)
    reservation_lock_path =
        _output_root_reservation_lock_path_from_canonical(canonical_outdir)
    for path in _same_output_root_reservation_paths(canonical_outdir)
        path in (cleanup_path, reservation_lock_path) && continue
        path in known_queued_markers && continue
        _output_root_path_entry_exists(path) && return true
    end
    return false
end

function _with_metamdbg_output_domain_lock(
        action::Function,
        outdir::AbstractString,
        ;
        allowed_same_root_locks::Tuple = (),
)::Any
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    _require_no_metamdbg_removal_quarantine_evidence!(canonical_outdir)
    reservation_lock_path =
        _output_root_reservation_lock_path_from_canonical(canonical_outdir)
    stale_age = _OUTPUT_ROOT_RESERVATION_STALE_AGE_SECONDS
    _require_exclusive_output_root_reservation(
        canonical_outdir,
        reservation_lock_path;
        subject = "metaMDBG outdir",
        stale_age,
        allowed_same_root_locks,
    )
    mkpath(dirname(reservation_lock_path))
    lock_handle = FileWatching.Pidfile.trymkpidlock(
        reservation_lock_path;
        stale_age,
        refresh = stale_age / 2,
    )
    lock_handle === false && throw(ArgumentError(
        "metaMDBG output is already reserved by another output-root " *
        "workflow: $(reservation_lock_path)",
    ))
    cleanup_reservation_path =
        _metamdbg_lifecycle_cleanup_reservation_path(canonical_outdir)
    cleanup_reservation_identity = nothing
    private_cleanup_failed = Ref(false)
    cleanup_failure_handler = function (
            _lock_path::AbstractString,
            _cleanup_error::Any,
    )
        private_cleanup_failed[] = true
        return nothing
    end
    try
        _require_exclusive_output_root_reservation(
            canonical_outdir,
            reservation_lock_path;
            subject = "metaMDBG outdir",
            stale_age,
            allowed_same_root_locks,
        )
        _publish_metamdbg_lifecycle_cleanup_reservation!(canonical_outdir)
        cleanup_reservation_identity =
            _metamdbg_output_lock_identity(cleanup_reservation_path)
        return _with_metamdbg_output_lock(
            action,
            canonical_outdir;
            lock_cleanup_failure_handler = cleanup_failure_handler,
        )
    finally
        try
            if cleanup_reservation_identity !== nothing
                if private_cleanup_failed[]
                    _require_unchanged_metamdbg_lifecycle_cleanup_reservation(
                        cleanup_reservation_path,
                        cleanup_reservation_identity,
                    )
                else
                    _remove_metamdbg_lifecycle_cleanup_reservation!(
                        cleanup_reservation_path,
                        cleanup_reservation_identity,
                    )
                end
            end
        finally
            Base.close(lock_handle)
        end
    end
end

function _metamdbg_input_snapshot(
        selected_input::NamedTuple,
)::Vector{NamedTuple}
    return map(selected_input.paths) do path
        isfile(path) || error(
            "metaMDBG input disappeared while capturing invocation metadata: " *
            "$(path).",
        )
        file_status = stat(path)
        modification_time_ns = round(
            Int64,
            file_status.mtime * 1_000_000_000,
        )
        return (;
            path,
            size_bytes = Int64(file_status.size),
            modification_time_ns,
        )
    end
end

function _require_current_metamdbg_input_snapshot!(
        selected_input::NamedTuple,
        expected_snapshot::AbstractVector{<:NamedTuple},
)::Vector{NamedTuple}
    current_snapshot = _metamdbg_input_snapshot(selected_input)
    current_snapshot == expected_snapshot || error(
        "metaMDBG input content contract changed after this invocation " *
        "captured its size and modification-time snapshot. Refusing stale " *
        "artifact reuse or execution.",
    )
    return current_snapshot
end

function _metamdbg_input_contract(
        selected_input::NamedTuple,
        abundance_min::Int,
        toolchain::NamedTuple = _metamdbg_expected_toolchain(),
        ;
        digest_function::Function = _metamdbg_sha256,
        input_snapshot::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
)::NamedTuple
    input_technology = if selected_input.flag == "--in-hifi"
        "hifi"
    elseif selected_input.flag == "--in-ont"
        "ont"
    else
        error("Unsupported metaMDBG input flag: $(selected_input.flag)")
    end
    captured_snapshot = input_snapshot === nothing ?
                        _metamdbg_input_snapshot(selected_input) :
                        collect(input_snapshot)
    snapshot_paths = getproperty.(captured_snapshot, :path)
    snapshot_paths == selected_input.paths || error(
        "metaMDBG invocation snapshot paths do not match the selected inputs.",
    )
    inputs = map(captured_snapshot) do input
        return (;
            path = input.path,
            size_bytes = input.size_bytes,
            sha256 = String(digest_function(input.path)),
        )
    end
    _require_current_metamdbg_input_snapshot!(
        selected_input,
        captured_snapshot,
    )
    contract = (
        schema_version = _METAMDBG_CONTRACT_SCHEMA_VERSION,
        input_technology,
        input_flag = selected_input.flag,
        platform_attestation = selected_input.platform_attestation,
        inputs,
        abundance_min,
        toolchain,
    )
    serialized_contract = JSON.json(contract)
    signature = SHA.bytes2hex(SHA.sha256(serialized_contract))
    contents = JSON.json((;
        schema_version = _METAMDBG_CONTRACT_SCHEMA_VERSION,
        signature_algorithm = "sha256",
        signature,
        contract,
    )) * "\n"
    return (;
        contract,
        signature,
        contents,
        input_snapshot = captured_snapshot,
    )
end

function _require_current_metamdbg_input_contract!(
        selected_input::NamedTuple,
        abundance_min::Int,
        expected_contract::NamedTuple,
        toolchain::NamedTuple = _metamdbg_expected_toolchain(),
        ;
        digest_function::Function = _metamdbg_sha256,
)::NamedTuple
    _require_current_metamdbg_input_snapshot!(
        selected_input,
        expected_contract.input_snapshot,
    )
    current_contract = _metamdbg_input_contract(
        selected_input,
        abundance_min,
        toolchain;
        digest_function,
        input_snapshot = expected_contract.input_snapshot,
    )
    current_contract.contents == expected_contract.contents || error(
        "metaMDBG input content contract changed after this invocation " *
        "captured its schema-v$(_METAMDBG_CONTRACT_SCHEMA_VERSION) " *
        "provenance. Refusing stale artifact reuse or execution.",
    )
    return current_contract
end

function _copy_metamdbg_input!(
        source::AbstractString,
        destination::AbstractString,
)::String
    open(source, "r") do input
        open(destination, "w") do output
            buffer = Vector{UInt8}(undef, 1024 * 1024)
            while !eof(input)
                bytes_read = readbytes!(input, buffer)
                bytes_read == 0 && break
                write(output, view(buffer, 1:bytes_read))
            end
            flush(output)
        end
    end
    chmod(destination, 0o400)
    return String(destination)
end

function _metamdbg_staged_input_name(
        input_index::Int,
        input_path::AbstractString,
)::String
    input_index > 0 || throw(ArgumentError(
        "metaMDBG staged input indexes must be positive.",
    ))
    filename = basename(String(input_path))
    normalized_filename = lowercase(filename)
    compression_suffix = ""
    for suffix in (".gz", ".bgz", ".bz2", ".xz", ".zst")
        if endswith(normalized_filename, suffix)
            compression_suffix = last(filename, length(suffix))
            break
        end
    end
    path_identity = first(_metamdbg_string_sha256(String(input_path)), 20)
    return "input-$(lpad(input_index, 6, '0'))-$(path_identity)" *
           compression_suffix
end

function _stage_metamdbg_inputs!(
        selected_input::NamedTuple,
        input_contract::NamedTuple,
        staging_parent::AbstractString,
)::NamedTuple
    staging_root = mktempdir(
        staging_parent;
        prefix = ".mycelia-metamdbg-inputs.",
    )
    chmod(staging_root, 0o700)
    staged_paths = String[]
    try
        for (input_index, input) in enumerate(input_contract.contract.inputs)
            input.path == selected_input.paths[input_index] || error(
                "metaMDBG staged-input contract path mismatch.",
            )
            staged_path = joinpath(
                staging_root,
                _metamdbg_staged_input_name(input_index, input.path),
            )
            _copy_metamdbg_input!(input.path, staged_path)
            if !isfile(staged_path) || islink(staged_path) ||
               filesize(staged_path) != input.size_bytes ||
               _metamdbg_sha256(staged_path) != input.sha256
                error(
                    "metaMDBG staged input does not match its captured size " *
                    "and SHA-256 contract: $(input.path).",
                )
            end
            push!(staged_paths, staged_path)
        end
    catch
        rm(staging_root; recursive = true, force = true)
        rethrow()
    end
    return (;
        root = staging_root,
        selected_input = (;
            flag = selected_input.flag,
            paths = staged_paths,
        ),
    )
end

function _require_unchanged_staged_metamdbg_inputs!(
        staged::NamedTuple,
        input_contract::NamedTuple,
)::Nothing
    for (staged_path, input) in
        zip(staged.selected_input.paths, input_contract.contract.inputs)
        if !isfile(staged_path) || islink(staged_path) ||
           filesize(staged_path) != input.size_bytes ||
           _metamdbg_sha256(staged_path) != input.sha256
            error(
                "metaMDBG staged input changed while assembly tools were " *
                "running: $(staged_path).",
            )
        end
    end
    return nothing
end

function _metamdbg_submission_reservation_prefix(
        outdir::AbstractString,
)::String
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    output_identity = _output_root_reservation_identity(normalized_outdir)
    return ".mycelia-metamdbg-submission.$(output_identity)."
end

function _metamdbg_pending_submission_job_prefix(
        outdir::AbstractString,
)::String
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    output_identity = _output_root_reservation_identity(normalized_outdir)
    return "$(_METAMDBG_PENDING_SUBMISSION_JOB_PREFIX)$(output_identity)."
end

function _metamdbg_pending_submission_job_path(
        outdir::AbstractString,
        capability::AbstractString,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString} = nothing,
)::String
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    normalized_capability = String(capability)
    occursin(r"^[0-9a-f]{64}$", normalized_capability) || error(
        "metaMDBG pending submission job record requires a lowercase " *
        "SHA-256 owner capability.",
    )
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    normalized_job_cluster = _normalize_metamdbg_job_cluster(job_cluster)
    cluster_component = normalized_job_cluster === nothing ?
                        "" :
                        ".$(normalized_job_cluster)"
    return joinpath(
        dirname(normalized_outdir),
        _metamdbg_pending_submission_job_prefix(normalized_outdir) *
        "$(normalized_capability).$(normalized_job_id)$(cluster_component).json",
    )
end

function _metamdbg_pending_submission_job_paths(
        outdir::AbstractString,
)::Vector{String}
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    parent = dirname(normalized_outdir)
    isdir(parent) || return String[]
    prefix = _metamdbg_pending_submission_job_prefix(normalized_outdir)
    return sort!(filter(
        path -> startswith(basename(path), prefix),
        readdir(parent; join = true),
    ))
end

function _metamdbg_pending_submission_job_path_parts(
        path::AbstractString,
        outdir::AbstractString,
)::NamedTuple
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    normalized_path = normpath(abspath(String(path)))
    dirname(normalized_path) == dirname(normalized_outdir) || error(
        "metaMDBG pending submission job record is outside its output-root " *
        "parent: $(normalized_path).",
    )
    prefix = _metamdbg_pending_submission_job_prefix(normalized_outdir)
    filename = basename(normalized_path)
    startswith(filename, prefix) || error(
        "metaMDBG pending submission job record has an invalid prefix: " *
        "$(normalized_path).",
    )
    ncodeunits(filename) > ncodeunits(prefix) || error(
        "metaMDBG found a malformed pending submission job record: " *
        "$(normalized_path).",
    )
    suffix = filename[(ncodeunits(prefix) + 1):end]
    parsed = match(
        r"^([0-9a-f]{64})\.([0-9]+)(?:\.([A-Za-z0-9_.-]+))?\.json$",
        suffix,
    )
    parsed === nothing && error(
        "metaMDBG found a malformed pending submission job record: " *
        "$(normalized_path). Remove it only after confirming no scheduler " *
        "job owns it.",
    )
    capability = parsed.captures[1]
    job_id = _normalize_metamdbg_job_id(parsed.captures[2])
    job_cluster = _normalize_metamdbg_job_cluster(parsed.captures[3])
    expected_path = _metamdbg_pending_submission_job_path(
        normalized_outdir,
        capability,
        job_id,
        job_cluster,
    )
    expected_path == normalized_path || error(
        "metaMDBG pending submission job record path is not canonical: " *
        "$(normalized_path).",
    )
    return (; path = normalized_path, capability, job_id, job_cluster)
end

function _metamdbg_submission_reservation_path(
        outdir::AbstractString,
        workflow_signature::AbstractString,
)::String
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    signature = String(workflow_signature)
    occursin(r"^[0-9a-f]{64}$", signature) || error(
        "metaMDBG submission reservation requires a lowercase SHA-256 " *
        "workflow signature, got $(repr(signature)).",
    )
    return joinpath(
        dirname(normalized_outdir),
        _metamdbg_submission_reservation_prefix(normalized_outdir) * signature,
    )
end

function _metamdbg_runtime_submission_reservation_path(
        outdir::AbstractString,
        capability::AbstractString,
)::String
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    normalized_capability = String(capability)
    occursin(r"^[0-9a-f]{64}$", normalized_capability) || error(
        "metaMDBG runtime reservation requires a lowercase SHA-256 owner " *
        "capability, got $(repr(normalized_capability)).",
    )
    return joinpath(
        dirname(normalized_outdir),
        _metamdbg_submission_reservation_prefix(normalized_outdir) *
        "runtime.$(normalized_capability)",
    )
end

function _metamdbg_consumed_submission_reservation_path(
        outdir::AbstractString,
        capability::AbstractString,
)::String
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    normalized_capability = String(capability)
    occursin(r"^[0-9a-f]{64}$", normalized_capability) || error(
        "metaMDBG consumed reservation requires a lowercase SHA-256 owner " *
        "capability, got $(repr(normalized_capability)).",
    )
    return joinpath(
        dirname(normalized_outdir),
        _metamdbg_submission_reservation_prefix(normalized_outdir) *
        "consumed.$(normalized_capability)",
    )
end

function _metamdbg_reclaiming_submission_reservation_path(
        outdir::AbstractString,
        capability::AbstractString,
)::String
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    normalized_capability = String(capability)
    occursin(r"^[0-9a-f]{64}$", normalized_capability) || error(
        "metaMDBG reclaiming reservation requires a lowercase SHA-256 owner " *
        "capability, got $(repr(normalized_capability)).",
    )
    return joinpath(
        dirname(normalized_outdir),
        _metamdbg_submission_reservation_prefix(normalized_outdir) *
        "reclaiming.$(normalized_capability)",
    )
end

function _metamdbg_scheduler_job_name(
        requested_job_name::AbstractString,
        workflow_signature::AbstractString,
)::String
    signature = String(workflow_signature)
    occursin(r"^[0-9a-f]{64}$", signature) || error(
        "metaMDBG scheduler job name requires a lowercase SHA-256 workflow " *
        "signature, got $(repr(signature)).",
    )
    requested_prefix = replace(
        strip(String(requested_job_name)),
        r"[^A-Za-z0-9_.-]+" => "-",
    )
    requested_prefix = replace(requested_prefix, r"-+" => "-")
    requested_prefix = replace(requested_prefix, r"^-|-$" => "")
    isempty(requested_prefix) && (requested_prefix = "metamdbg")
    suffix = "-mycelia-metamdbg-$(signature)"
    maximum_prefix_length = 128 - ncodeunits(suffix)
    prefix = first(
        requested_prefix,
        min(length(requested_prefix), maximum_prefix_length),
    )
    return "$(prefix)$(suffix)"
end

function _metamdbg_scheduler_job_name_prefix(
        scheduler_job_name::AbstractString,
        workflow_signature::AbstractString,
)::String
    normalized_job_name = String(scheduler_job_name)
    signature = String(workflow_signature)
    suffix = "-mycelia-metamdbg-$(signature)"
    occursin(r"^[A-Za-z0-9_.-]+$", normalized_job_name) &&
        ncodeunits(normalized_job_name) <= 128 &&
        endswith(normalized_job_name, suffix) || error(
        "metaMDBG scheduler job name is not bound to its workflow signature.",
    )
    prefix_length = ncodeunits(normalized_job_name) - ncodeunits(suffix)
    prefix_length > 0 || error(
        "metaMDBG scheduler job name has no requested-name prefix.",
    )
    prefix = normalized_job_name[1:prefix_length]
    _metamdbg_scheduler_job_name(prefix, signature) == normalized_job_name ||
        error("metaMDBG scheduler job name is not canonical.")
    return prefix
end

function _metamdbg_workflow_contract(
        outputs::NamedTuple,
        input_contract::NamedTuple,
        graph_k::Int,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    contract = (;
        schema_version = _METAMDBG_SUBMISSION_RESERVATION_SCHEMA_VERSION,
        canonical_outdir = outputs.outdir,
        input_contract_signature = input_contract.signature,
        graph_k,
    )
    signature = _metamdbg_string_sha256(JSON.json(contract))
    return (; contract, signature)
end

function _metamdbg_output_root_reservation_capability(
        workflow_signature::AbstractString,
        owner_token::AbstractString,
)::String
    normalized_signature = String(workflow_signature)
    occursin(r"^[0-9a-f]{64}$", normalized_signature) || error(
        "metaMDBG output-root reservation requires a lowercase SHA-256 " *
        "workflow signature.",
    )
    normalized_owner_token = String(owner_token)
    isempty(normalized_owner_token) && error(
        "metaMDBG output-root reservation owner token must be nonempty.",
    )
    return _metamdbg_string_sha256(JSON.json((;
        workflow_signature = normalized_signature,
        owner_token = normalized_owner_token,
    )))
end

function _metamdbg_output_root_reservation_metadata(
        outputs::NamedTuple,
        workflow_signature::AbstractString,
        owner_token::AbstractString,
)::NamedTuple
    capability = _metamdbg_output_root_reservation_capability(
        workflow_signature,
        owner_token,
    )
    queued_marker = _output_root_durable_reservation_path_from_canonical(
        outputs.outdir,
        "metamdbg-queued-$(capability)",
    )
    runtime_marker = _output_root_durable_reservation_path_from_canonical(
        outputs.outdir,
        "metamdbg-runtime-$(capability)",
    )
    contents = JSON.json((;
        schema_version = _METAMDBG_OUTPUT_ROOT_RESERVATION_SCHEMA_VERSION,
        canonical_outdir = outputs.outdir,
        workflow_signature = String(workflow_signature),
        capability,
    )) * "\n"
    return (;
        output_root_reservation_capability = capability,
        output_root_reservation_marker = queued_marker,
        runtime_output_root_reservation_marker = runtime_marker,
        output_root_reservation_contents = contents,
    )
end

function _metamdbg_submission_reservation(
        outputs::NamedTuple,
        input_contract::NamedTuple,
        graph_k::Int,
        ;
        owner_token::AbstractString = string(UUIDs.uuid4()),
        job_name::AbstractString = "metamdbg",
)::NamedTuple
    workflow = _metamdbg_workflow_contract(
        outputs,
        input_contract,
        graph_k,
    )
    reservation_contract = workflow.contract
    workflow_signature = workflow.signature
    scheduler_job_name = _metamdbg_scheduler_job_name(
        job_name,
        workflow_signature,
    )
    path = _metamdbg_submission_reservation_path(
        outputs.outdir,
        workflow_signature,
    )
    contents = JSON.json((;
        reservation_contract...,
        workflow_signature,
        scheduler_job_name,
        owner_token = String(owner_token),
    )) * "\n"
    output_root_reservation =
        _metamdbg_output_root_reservation_metadata(
            outputs,
            workflow_signature,
            owner_token,
        )
    return (;
        path,
        runtime_path = _metamdbg_runtime_submission_reservation_path(
            outputs.outdir,
            output_root_reservation.output_root_reservation_capability,
        ),
        consumed_path = _metamdbg_consumed_submission_reservation_path(
            outputs.outdir,
            output_root_reservation.output_root_reservation_capability,
        ),
        reclaiming_path = _metamdbg_reclaiming_submission_reservation_path(
            outputs.outdir,
            output_root_reservation.output_root_reservation_capability,
        ),
        contract_marker = joinpath(
            path,
            _METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
        ),
        job_marker = joinpath(
            path,
            _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME,
        ),
        canonical_outdir = outputs.outdir,
        input_contract_signature = input_contract.signature,
        graph_k,
        workflow_signature,
        scheduler_job_name,
        owner_token = String(owner_token),
        output_root_reservation...,
        contents,
        job_id = nothing,
        job_cluster = nothing,
        job_schema_version = nothing,
        job_contents = nothing,
        pending_job_marker = nothing,
        submission_state = :reserved,
    )
end

function _normalize_metamdbg_job_id(job_id::AbstractString)::String
    normalized_job_id = strip(String(job_id))
    occursin(r"^[0-9]+$", normalized_job_id) || throw(
        ArgumentError(
            "metaMDBG scheduler job_id must contain only decimal digits.",
        ),
    )
    return normalized_job_id
end

function _normalize_metamdbg_job_cluster(
        job_cluster::Union{Nothing, AbstractString},
)::Union{Nothing, String}
    job_cluster === nothing && return nothing
    normalized_job_cluster = strip(String(job_cluster))
    occursin(r"^[A-Za-z0-9_.-]+$", normalized_job_cluster) || throw(
        ArgumentError(
            "metaMDBG scheduler job_cluster must contain only letters, " *
            "decimal digits, underscore, dot, or hyphen.",
        ),
    )
    return normalized_job_cluster
end

function _metamdbg_submission_job_parts(
        reservation::NamedTuple,
        job_cluster::Union{Nothing, AbstractString} = nothing,
        job_schema_version::Int =
            _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
)::Tuple{String, String}
    job_schema_version in (
        _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
        _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
    ) || throw(ArgumentError(
        "metaMDBG submission job schema version is unsupported: " *
        "$(job_schema_version).",
    ))
    job_schema_version ==
        _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION &&
        job_cluster !== nothing && throw(ArgumentError(
        "Legacy metaMDBG submission job records cannot name a cluster.",
    ))
    prefix_object = JSON.json((;
        schema_version = job_schema_version,
        workflow_signature = reservation.workflow_signature,
        owner_token = reservation.owner_token,
    ))
    endswith(prefix_object, "}") || error(
        "metaMDBG failed to construct its submission job record prefix.",
    )
    normalized_job_cluster = _normalize_metamdbg_job_cluster(job_cluster)
    prefix = chop(prefix_object; tail = 1) * ",\"job_id\":\""
    suffix = if job_schema_version ==
                _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION
        "\"}\n"
    else
        "\",\"job_cluster\":" * JSON.json(normalized_job_cluster) * "}\n"
    end
    return prefix, suffix
end

function _metamdbg_submission_job_contents(
        reservation::NamedTuple,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString} = nothing,
        job_schema_version::Int =
            _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
)::String
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    prefix, suffix = _metamdbg_submission_job_parts(
        reservation,
        job_cluster,
        job_schema_version,
    )
    return prefix * normalized_job_id * suffix
end

function _metamdbg_bound_submission_reservation(
        reservation::NamedTuple,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString} = nothing,
        job_schema_version::Int =
            _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
)::NamedTuple
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    normalized_job_cluster = _normalize_metamdbg_job_cluster(job_cluster)
    job_schema_version ==
        _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION &&
        normalized_job_cluster !== nothing && throw(ArgumentError(
        "Legacy metaMDBG submission job records cannot name a cluster.",
    ))
    return merge(reservation, (;
        job_id = normalized_job_id,
        job_cluster = normalized_job_cluster,
        job_schema_version,
        job_contents = _metamdbg_submission_job_contents(
            reservation,
            normalized_job_id,
            normalized_job_cluster,
            job_schema_version,
        ),
        pending_job_marker = _metamdbg_pending_submission_job_path(
            reservation.canonical_outdir,
            reservation.output_root_reservation_capability,
            normalized_job_id,
            normalized_job_cluster,
        ),
        submission_state = :submitted,
    ))
end

function _metamdbg_pending_submission_job_identity(
        path::AbstractString,
)::NamedTuple
    normalized_path = normpath(abspath(String(path)))
    isfile(normalized_path) && !islink(normalized_path) || error(
        "metaMDBG pending submission job record must be a regular, " *
        "non-symlink file: $(normalized_path).",
    )
    pending_status = stat(normalized_path)
    pending_status.uid == Base.Libc.getuid() || error(
        "metaMDBG pending submission job record is not owned by the current " *
        "user: $(normalized_path).",
    )
    (pending_status.mode & 0o777) == 0o600 || error(
        "metaMDBG pending submission job record must have mode 0600: " *
        "$(normalized_path).",
    )
    return (; device = pending_status.device, inode = pending_status.inode)
end

function _metamdbg_pending_submission_job_descriptor_identity(
        pending_io::Base.Filesystem.File,
        path::AbstractString,
)::NamedTuple
    pending_status = stat(pending_io)
    pending_status.uid == Base.Libc.getuid() || error(
        "metaMDBG pending submission job descriptor is not owned by the " *
        "current user: $(path).",
    )
    (pending_status.mode & 0o777) == 0o600 || error(
        "metaMDBG pending submission job descriptor must have mode 0600: " *
        "$(path).",
    )
    return (; device = pending_status.device, inode = pending_status.inode)
end

function _metamdbg_submission_job_content_candidates(
        reservation::NamedTuple,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString},
)::Vector{NamedTuple}
    normalized_job_cluster = _normalize_metamdbg_job_cluster(job_cluster)
    candidates = NamedTuple[(;
        job_schema_version =
            _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
        contents = _metamdbg_submission_job_contents(
            reservation,
            job_id,
            normalized_job_cluster,
            _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
        ),
    )]
    if normalized_job_cluster === nothing
        push!(candidates, (;
            job_schema_version =
                _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
            contents = _metamdbg_submission_job_contents(
                reservation,
                job_id,
                nothing,
                _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
            ),
        ))
    end
    return candidates
end

function _metamdbg_pending_submission_job_record(
        reservation::NamedTuple,
        path::AbstractString;
        allow_empty::Bool = false,
        allow_incomplete::Bool = false,
)::NamedTuple
    parts = _metamdbg_pending_submission_job_path_parts(
        path,
        reservation.canonical_outdir,
    )
    parts.capability == reservation.output_root_reservation_capability || error(
        "metaMDBG pending submission job record does not match its exact " *
        "owner capability: $(parts.path).",
    )
    identity = _metamdbg_pending_submission_job_identity(parts.path)
    pending_io = Base.Filesystem.open(
        parts.path,
        Base.JL_O_RDONLY | Base.JL_O_CLOEXEC,
    )
    bytes = try
        _metamdbg_pending_submission_job_descriptor_identity(
            pending_io,
            parts.path,
        ) == identity || error(
            "metaMDBG pending submission job record was replaced before its " *
            "contents were read: $(parts.path).",
        )
        captured_bytes = read(pending_io)
        _metamdbg_pending_submission_job_descriptor_identity(
            pending_io,
            parts.path,
        ) == identity || error(
            "metaMDBG pending submission job descriptor changed while its " *
            "contents were read: $(parts.path).",
        )
        captured_bytes
    finally
        close(pending_io)
    end
    _metamdbg_pending_submission_job_identity(parts.path) == identity || error(
        "metaMDBG pending submission job record was replaced while being " *
        "validated: $(parts.path).",
    )
    content_candidates = _metamdbg_submission_job_content_candidates(
        reservation,
        parts.job_id,
        parts.job_cluster,
    )
    candidate_bytes = map(content_candidates) do candidate
        return collect(codeunits(candidate.contents))
    end
    complete_indices = findall(expected -> bytes == expected, candidate_bytes)
    prefix_indices = findall(candidate_bytes) do expected
        return length(bytes) <= length(expected) &&
               bytes == expected[1:length(bytes)]
    end
    length(complete_indices) <= 1 || error(
        "metaMDBG pending submission job record matches multiple complete " *
        "schemas: $(parts.path).",
    )
    is_complete = length(complete_indices) == 1
    is_recoverable_prefix = allow_incomplete && !isempty(prefix_indices)
    if !is_complete && !is_recoverable_prefix
        if isempty(bytes) && allow_empty
            nothing
        else
            error(
                "metaMDBG pending submission job record is not canonical or " *
                "does not match its workflow owner: $(parts.path).",
            )
        end
    end
    selected_index = is_complete ? only(complete_indices) :
                     (isempty(prefix_indices) ? 1 : first(prefix_indices))
    selected_candidate = content_candidates[selected_index]
    return merge(parts, (;
        bytes,
        identity,
        is_complete,
        job_schema_version = selected_candidate.job_schema_version,
    ))
end

function _require_unchanged_metamdbg_pending_submission_job_record(
        reservation::NamedTuple,
        expected::NamedTuple;
        allow_empty::Bool = false,
        allow_incomplete::Bool = false,
)::NamedTuple
    observed = _metamdbg_pending_submission_job_record(
        reservation,
        expected.path;
        allow_empty,
        allow_incomplete,
    )
    for field in (
            :path,
            :capability,
            :job_id,
            :job_cluster,
            :job_schema_version,
            :bytes,
            :identity,
            :is_complete,
        )
        getproperty(observed, field) == getproperty(expected, field) || error(
            "metaMDBG pending submission job record $(field) changed during " *
            "recovery inspection: $(expected.path).",
        )
    end
    return observed
end

function _metamdbg_pending_submission_job_record_for_reservation(
        reservation::NamedTuple,
        pending_paths::AbstractVector{<:AbstractString};
        allow_empty::Bool = false,
        allow_incomplete::Bool = false,
)::Union{Nothing, NamedTuple}
    candidates = filter(pending_paths) do path
        parts = _metamdbg_pending_submission_job_path_parts(
            path,
            reservation.canonical_outdir,
        )
        return parts.capability ==
               reservation.output_root_reservation_capability
    end
    length(candidates) <= 1 || error(
        "metaMDBG found multiple pending scheduler job records for one exact " *
        "owner capability: $(join(candidates, ", ")).",
    )
    isempty(candidates) && return nothing
    return _metamdbg_pending_submission_job_record(
        reservation,
        only(candidates);
        allow_empty,
        allow_incomplete,
    )
end

function _complete_metamdbg_pending_submission_job_record!(
        reservation::NamedTuple,
        pending::NamedTuple;
        pre_descriptor_open_hook::Function =
            (_pending::NamedTuple) -> nothing,
)::NamedTuple
    pending.is_complete && return pending
    _require_unchanged_metamdbg_pending_submission_job_record(
        reservation,
        pending;
        allow_incomplete = true,
    )
    expected_bytes = collect(codeunits(_metamdbg_submission_job_contents(
        reservation,
        pending.job_id,
        pending.job_cluster,
        pending.job_schema_version,
    )))
    length(pending.bytes) < length(expected_bytes) || error(
        "metaMDBG incomplete pending scheduler job record is not a strict " *
        "canonical prefix: $(pending.path).",
    )
    remaining_bytes = expected_bytes[(length(pending.bytes) + 1):end]
    pre_descriptor_open_hook(pending)
    pending_io = Base.Filesystem.open(
        pending.path,
        Base.JL_O_RDWR | Base.JL_O_CLOEXEC,
    )
    try
        _metamdbg_pending_submission_job_descriptor_identity(
            pending_io,
            pending.path,
        ) == pending.identity || error(
            "metaMDBG pending scheduler job descriptor was replaced before " *
            "prefix completion: $(pending.path).",
        )
        descriptor_bytes = read(pending_io)
        descriptor_bytes == pending.bytes || error(
            "metaMDBG pending scheduler job contents changed before prefix " *
            "completion: $(pending.path).",
        )
        _metamdbg_pending_submission_job_identity(pending.path) ==
            pending.identity || error(
            "metaMDBG pending scheduler job path was replaced before prefix " *
            "completion: $(pending.path).",
        )
        seekend(pending_io)
        write(pending_io, remaining_bytes)
        _fsync_metamdbg_descriptor(Base.fd(pending_io), pending.path)
        _metamdbg_pending_submission_job_descriptor_identity(
            pending_io,
            pending.path,
        ) == pending.identity || error(
            "metaMDBG pending scheduler job descriptor changed during prefix " *
            "completion: $(pending.path).",
        )
    finally
        close(pending_io)
    end
    _fsync_metamdbg_directory(dirname(pending.path))
    completed = _metamdbg_pending_submission_job_record(
        reservation,
        pending.path,
    )
    completed.identity == pending.identity || error(
        "metaMDBG pending scheduler job record was replaced while completing " *
        "its canonical prefix: $(pending.path).",
    )
    return completed
end

function _remove_exact_metamdbg_job_marker!(
        marker::AbstractString,
        expected_identity::NamedTuple,
)::Nothing
    _output_root_path_entry_exists(marker) || return nothing
    _metamdbg_pending_submission_job_identity(marker) == expected_identity ||
        error(
            "metaMDBG refuses to remove a replacement scheduler job marker: " *
            "$(marker).",
        )
    _remove_exact_metamdbg_durable_file!(
        marker,
        merge(expected_identity, (; path = normpath(abspath(String(marker))))),
    )
    return nothing
end

function _publish_metamdbg_pending_submission_job_record!(
        reservation::NamedTuple,
        pending::NamedTuple,
        reservation_identity::NamedTuple,
        shared_reservation_identity::NamedTuple,
)::NamedTuple
    pending.is_complete || error(
        "metaMDBG refuses to publish an incomplete pending scheduler job " *
        "record: $(pending.path).",
    )
    _require_unchanged_metamdbg_recovery_identities(
        reservation,
        reservation_identity,
        shared_reservation_identity,
    )
    pending = _require_unchanged_metamdbg_pending_submission_job_record(
        reservation,
        pending,
    )
    _metamdbg_pending_submission_job_paths(reservation.canonical_outdir) ==
        [pending.path] || error(
        "metaMDBG found ambiguous pending scheduler job evidence before exact " *
        "publication.",
    )
    marker = reservation.job_marker
    bound_reservation = _metamdbg_bound_submission_reservation(
        reservation,
        pending.job_id,
        pending.job_cluster,
        pending.job_schema_version,
    )
    !_output_root_path_entry_exists(marker) || error(
        "metaMDBG refuses to overwrite an existing durable submission job " *
        "record: $(marker).",
    )
    linked = false
    committed = false
    try
        Base.Filesystem.hardlink(pending.path, marker)
        linked = true
        _metamdbg_pending_submission_job_identity(marker) ==
            pending.identity || error(
            "metaMDBG newly linked scheduler job marker does not retain the " *
            "validated pending inode.",
        )
        _metamdbg_pending_submission_job_identity(pending.path) ==
            pending.identity || error(
            "metaMDBG pending scheduler job path was replaced during " *
            "publication.",
        )
        Base.Filesystem.samefile(marker, pending.path) || error(
            "metaMDBG scheduler job marker and pending record are not the " *
            "same publication inode.",
        )
        _metamdbg_pending_submission_job_paths(reservation.canonical_outdir) ==
            [pending.path] || error(
            "metaMDBG pending scheduler job inventory changed during exact " *
            "publication.",
        )
        _require_unchanged_metamdbg_recovery_identities(
            bound_reservation,
            reservation_identity,
            shared_reservation_identity,
        )
        _fsync_metamdbg_directory(reservation.path)
        committed = true
        _metamdbg_pending_submission_job_identity(marker) ==
            pending.identity || error(
            "metaMDBG committed scheduler job marker changed after its owner " *
            "directory fsync.",
        )
        _require_unchanged_metamdbg_recovery_identities(
            bound_reservation,
            reservation_identity,
            shared_reservation_identity,
        )
    catch primary_error
        if linked && !committed
            try
                _remove_exact_metamdbg_job_marker!(
                    marker,
                    pending.identity,
                )
            catch cleanup_error
                @warn "metaMDBG failed to roll back an uncommitted exact " *
                      "scheduler job marker while preserving the primary " *
                      "publication failure" marker primary_error cleanup_error
            end
        end
        Base.rethrow()
    end
    return _metamdbg_submission_reservation_from_path(
        reservation.path,
        reservation.canonical_outdir,
    )
end

function _remove_metamdbg_pending_submission_job_record!(
        reservation::NamedTuple,
        pending::NamedTuple;
        pending_remover::Union{Nothing, Function} = nothing,
)::Nothing
    _require_unchanged_metamdbg_pending_submission_job_record(
        reservation,
        pending,
    )
    # Retained only as an internal fault-injection hook. Production cleanup
    # never delegates pathname removal to this callback.
    pending_remover === nothing || pending_remover(pending.path)
    _remove_exact_metamdbg_durable_file!(
        pending.path,
        merge(pending.identity, (; path = pending.path)),
    )
    return nothing
end

function _resume_metamdbg_pending_submission_job_record!(
        reservation::NamedTuple,
        pending::NamedTuple,
        reservation_identity::NamedTuple,
        shared_reservation_identity::NamedTuple;
        pre_promotion_hook::Function =
            (_reservation::NamedTuple, _pending::NamedTuple) -> nothing,
        pending_remover::Union{Nothing, Function} = nothing,
)::NamedTuple
    _require_unchanged_metamdbg_recovery_identities(
        reservation,
        reservation_identity,
        shared_reservation_identity,
    )
    pending = _complete_metamdbg_pending_submission_job_record!(
        reservation,
        pending,
    )
    _require_unchanged_metamdbg_recovery_identities(
        reservation,
        reservation_identity,
        shared_reservation_identity,
    )
    pre_promotion_hook(reservation, pending)
    _require_unchanged_metamdbg_recovery_identities(
        reservation,
        reservation_identity,
        shared_reservation_identity,
    )
    pending = _require_unchanged_metamdbg_pending_submission_job_record(
        reservation,
        pending,
    )
    bound = if reservation.job_id === nothing
        _publish_metamdbg_pending_submission_job_record!(
            reservation,
            pending,
            reservation_identity,
            shared_reservation_identity,
        )
    else
        reservation.job_id == pending.job_id || error(
            "metaMDBG committed scheduler job ID conflicts with its pending " *
            "record.",
        )
        reservation.job_cluster == pending.job_cluster || error(
            "metaMDBG committed scheduler cluster conflicts with its pending " *
            "record.",
        )
        reservation.job_schema_version == pending.job_schema_version || error(
            "metaMDBG committed scheduler job schema conflicts with its " *
            "pending record.",
        )
        pending.bytes == collect(codeunits(reservation.job_contents)) || error(
            "metaMDBG committed scheduler job content conflicts with its " *
            "pending record.",
        )
        _metamdbg_pending_submission_job_identity(reservation.job_marker) ==
            pending.identity || error(
            "metaMDBG committed scheduler job marker does not retain the " *
            "validated pending inode.",
        )
        Base.Filesystem.samefile(
            reservation.job_marker,
            pending.path,
        ) || error(
            "metaMDBG committed scheduler job marker is not the pending " *
            "publication inode.",
        )
        reservation
    end
    _remove_metamdbg_pending_submission_job_record!(
        bound,
        pending;
        pending_remover,
    )
    _require_unchanged_metamdbg_recovery_identities(
        bound,
        reservation_identity,
        shared_reservation_identity,
    )
    return _metamdbg_submission_reservation_from_path(
        bound.path,
        bound.canonical_outdir,
    )
end

function _recover_metamdbg_pending_submission_job_records!(
        outdir::AbstractString;
        pre_promotion_hook::Function =
            (_reservation::NamedTuple, _pending::NamedTuple) -> nothing,
        post_initial_snapshot_hook::Function =
            (_records::Vector{NamedTuple}) -> nothing,
)::Nothing
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    pending_paths =
        _metamdbg_pending_submission_job_paths(canonical_outdir)
    isempty(pending_paths) && return nothing
    owner_paths = _metamdbg_submission_reservation_paths(
        canonical_outdir;
        include_consumed = true,
    )
    owners_by_capability = Dict{String, NamedTuple}()
    for owner_path in owner_paths
        owner = _metamdbg_submission_reservation_from_path(
            owner_path,
            canonical_outdir;
            allow_provisional = true,
        )
        capability = owner.output_root_reservation_capability
        haskey(owners_by_capability, capability) && error(
            "metaMDBG found multiple owner records for one pending scheduler " *
            "job capability: $(capability).",
        )
        owners_by_capability[capability] = owner
    end
    recovery_records = NamedTuple[]
    seen_capabilities = Set{String}()
    for pending_path in pending_paths
        parts = _metamdbg_pending_submission_job_path_parts(
            pending_path,
            canonical_outdir,
        )
        owner = get(owners_by_capability, parts.capability, nothing)
        owner isa NamedTuple || error(
            "metaMDBG pending scheduler job record has no exact durable owner " *
            "record: $(pending_path).",
        )
        owner.reservation_state == :queued || error(
            "metaMDBG pending scheduler job record is paired with unsupported " *
            "owner state $(repr(owner.reservation_state)): $(pending_path).",
        )
        parts.capability in seen_capabilities && error(
            "metaMDBG found multiple pending scheduler job records for one " *
            "exact owner capability: $(parts.capability).",
        )
        push!(seen_capabilities, parts.capability)
        pending = _metamdbg_pending_submission_job_record(
            owner,
            pending_path;
            allow_incomplete = true,
        )
        reservation_identity =
            _metamdbg_submission_reservation_identity(owner)
        shared_reservation_identity =
            _metamdbg_shared_reservation_identity(owner)
        push!(recovery_records, (;
            owner,
            pending,
            reservation_identity,
            shared_reservation_identity,
        ))
    end
    post_initial_snapshot_hook(recovery_records)
    for recovery in recovery_records
        recovered_owner = _resume_metamdbg_pending_submission_job_record!(
            recovery.owner,
            recovery.pending,
            recovery.reservation_identity,
            recovery.shared_reservation_identity;
            pre_promotion_hook,
        )
        _require_metamdbg_submission_reservation!(recovered_owner)
    end
    isempty(_metamdbg_pending_submission_job_paths(canonical_outdir)) || error(
        "metaMDBG pending scheduler job records remain after confirmed-dead " *
        "recovery.",
    )
    return nothing
end

function _metamdbg_submission_reservation_path_state(
        path::AbstractString,
        canonical_outdir::AbstractString,
)::Symbol
    normalized_path = normpath(abspath(String(path)))
    normalized_outdir = _metamdbg_canonical_output_path(canonical_outdir)
    dirname(normalized_path) == dirname(normalized_outdir) || error(
        "metaMDBG submission reservation is outside its output-root parent: " *
        "$(normalized_path).",
    )
    prefix = _metamdbg_submission_reservation_prefix(normalized_outdir)
    filename = basename(normalized_path)
    startswith(filename, prefix) || error(
        "metaMDBG submission reservation path has an invalid prefix: " *
        "$(normalized_path).",
    )
    suffix = filename[(ncodeunits(prefix) + 1):end]
    if occursin(r"^[0-9a-f]{64}$", suffix)
        return :published
    elseif occursin(r"^runtime\.[0-9a-f]{64}$", suffix)
        return :runtime
    elseif occursin(r"^reclaiming\.[0-9a-f]{64}$", suffix)
        return :reclaiming
    elseif startswith(suffix, "tmp.") && ncodeunits(suffix) > 4
        return :provisional
    elseif occursin(r"^consumed\.[0-9a-f]{64}$", suffix)
        return :consumed
    elseif startswith(suffix, "consumed.") && ncodeunits(suffix) > 9
        return :legacy_consumed
    end
    error(
        "metaMDBG found a malformed submission reservation entry: " *
        "$(normalized_path). Remove it only after confirming no job owns it.",
    )
end

function _metamdbg_has_potential_queued_output_root_reservation(
        outdir::AbstractString,
)::Bool
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    reservation_lock_path =
        _output_root_reservation_lock_path_from_canonical(canonical_outdir)
    cleanup_reservation_path =
        _metamdbg_lifecycle_cleanup_reservation_path(canonical_outdir)
    for path in _same_output_root_reservation_paths(canonical_outdir)
        path == reservation_lock_path && continue
        path == cleanup_reservation_path && continue
        _is_output_root_durable_reservation_path(path) || continue
        if _output_root_path_entry_exists(path)
            return true
        end
    end
    return false
end

function _metamdbg_submission_reservation_paths(
        outdir::AbstractString,
        ;
        include_consumed::Bool = false,
)::Vector{String}
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    parent = dirname(normalized_outdir)
    isdir(parent) || return String[]
    prefix = _metamdbg_submission_reservation_prefix(normalized_outdir)
    matching_paths = sort!(filter(
        path -> startswith(basename(path), prefix),
        readdir(parent; join = true),
    ))
    published_paths = String[]
    provisional_candidates = String[]
    for path in matching_paths
        path_state = _metamdbg_submission_reservation_path_state(
            path,
            normalized_outdir,
        )
        if path_state == :published
            push!(published_paths, path)
        elseif path_state == :provisional
            push!(provisional_candidates, path)
        elseif path_state == :runtime
            push!(published_paths, path)
        elseif path_state == :reclaiming
            push!(published_paths, path)
        elseif path_state == :consumed
            include_consumed && push!(published_paths, path)
        elseif path_state == :legacy_consumed
            continue
        end
    end
    active_provisional_paths = String[]
    for path in provisional_candidates
        provisional = try
            _metamdbg_submission_reservation_from_path(
                path,
                normalized_outdir;
                allow_provisional = true,
                require_output_root_reservation = false,
            )
        catch caught
            caught isa InterruptException && rethrow()
            if _metamdbg_has_potential_queued_output_root_reservation(
                    normalized_outdir,
            )
                error(
                    "metaMDBG could not validate a provisional submission " *
                    "reservation while a potentially paired shared " *
                    "output-root reservation remains: $(path). Refusing " *
                    "capability-blind recovery. Cause: " *
                    sprint(showerror, caught),
                )
            end
            # A partially written temporary directory is nonblocking only
            # when no queued shared marker could pair with it.
            continue
        end
        if _output_root_path_entry_exists(
                provisional.output_root_reservation_marker,
        )
            _require_metamdbg_submission_reservation!(provisional)
            push!(active_provisional_paths, path)
        end
    end
    return sort!(vcat(published_paths, active_provisional_paths))
end

function _fsync_metamdbg_descriptor(
        descriptor::Union{Integer, Base.RawFD},
        label::AbstractString,
)::Nothing
    raw_descriptor = descriptor isa Base.RawFD ?
                     reinterpret(Cint, descriptor) : Cint(descriptor)
    result = ccall(:fsync, Cint, (Cint,), raw_descriptor)
    if result != 0
        saved_errno = Base.Libc.errno()
        throw(SystemError("fsync $(label)", saved_errno))
    end
    return nothing
end

function _fsync_metamdbg_file(
        output::IOStream,
        label::AbstractString,
)::Nothing
    flush(output)
    _fsync_metamdbg_descriptor(Base.fd(output), label)
    return nothing
end

function _metamdbg_directory_path_identity(
        path::AbstractString,
)::NamedTuple
    normalized_path = normpath(abspath(String(path)))
    isdir(normalized_path) && !islink(normalized_path) || error(
        "metaMDBG durable directory must be a regular, non-symlink " *
        "directory: $(normalized_path).",
    )
    directory_status = stat(normalized_path)
    return (;
        path = normalized_path,
        device = directory_status.device,
        inode = directory_status.inode,
    )
end

function _metamdbg_directory_descriptor_identity(
        directory::Base.Filesystem.File,
        path::AbstractString,
)::NamedTuple
    directory_status = stat(directory)
    return (;
        path = normpath(abspath(String(path))),
        device = directory_status.device,
        inode = directory_status.inode,
    )
end

function _require_unchanged_metamdbg_directory_descriptor(
        directory::Base.Filesystem.File,
        expected_identity::NamedTuple,
        label::AbstractString,
)::Nothing
    observed_identity = _metamdbg_directory_descriptor_identity(
        directory,
        expected_identity.path,
    )
    observed_identity == expected_identity || error(
        "metaMDBG durable directory descriptor changed for $(label): " *
        "$(expected_identity.path).",
    )
    return nothing
end

function _fsync_bound_metamdbg_directory(
        directory::Base.Filesystem.File,
        expected_identity::NamedTuple,
        label::AbstractString,
)::Nothing
    _require_unchanged_metamdbg_directory_descriptor(
        directory,
        expected_identity,
        label,
    )
    _fsync_metamdbg_descriptor(Base.fd(directory), label)
    _require_unchanged_metamdbg_directory_descriptor(
        directory,
        expected_identity,
        label,
    )
    return nothing
end

function _fsync_metamdbg_directory(
        path::AbstractString;
        expected_identity::Union{Nothing, NamedTuple} = nothing,
        post_descriptor_open_hook::Function =
            (_path::AbstractString, _identity::NamedTuple) -> nothing,
        post_fsync_hook::Function =
            (_path::AbstractString, _identity::NamedTuple) -> nothing,
)::Nothing
    bound_identity = _metamdbg_directory_path_identity(path)
    if expected_identity !== nothing
        bound_identity == expected_identity || error(
            "metaMDBG durable directory path changed before descriptor " *
            "acquisition: $(bound_identity.path).",
        )
    end
    directory = Base.Filesystem.open(
        bound_identity.path,
        Base.JL_O_RDONLY |
        Base.JL_O_DIRECTORY |
        Base.JL_O_NOFOLLOW |
        Base.JL_O_CLOEXEC,
    )
    try
        descriptor_identity = _metamdbg_directory_descriptor_identity(
            directory,
            bound_identity.path,
        )
        descriptor_identity == bound_identity || error(
            "metaMDBG durable directory was replaced before its descriptor " *
            "could be bound: $(bound_identity.path).",
        )
        post_descriptor_open_hook(bound_identity.path, bound_identity)
        _metamdbg_directory_path_identity(bound_identity.path) ==
            bound_identity || error(
            "metaMDBG durable directory path changed before fsync: " *
            "$(bound_identity.path).",
        )
        _fsync_bound_metamdbg_directory(
            directory,
            bound_identity,
            bound_identity.path,
        )
        post_fsync_hook(bound_identity.path, bound_identity)
        _metamdbg_directory_path_identity(bound_identity.path) ==
            bound_identity || error(
            "metaMDBG durable directory path changed during fsync: " *
            "$(bound_identity.path).",
        )
    finally
        close(directory)
    end
    return nothing
end

function _metamdbg_regular_file_identity(
        path::AbstractString,
)::NamedTuple
    normalized_path = normpath(abspath(String(path)))
    isfile(normalized_path) && !islink(normalized_path) || error(
        "metaMDBG durable file must be a regular, non-symlink file: " *
        "$(normalized_path).",
    )
    file_status = stat(normalized_path)
    return (;
        path = normalized_path,
        device = file_status.device,
        inode = file_status.inode,
    )
end

function _metamdbg_path_identity_matches(
        observed::NamedTuple,
        expected::NamedTuple,
)::Bool
    required_fields = (:path, :device, :inode)
    all(field -> hasproperty(observed, field), required_fields) || return false
    all(field -> hasproperty(expected, field), required_fields) || return false
    return all(
        getproperty(observed, field) == getproperty(expected, field)
        for field in required_fields
    )
end

function _metamdbg_file_descriptor_number(
        file::Base.Filesystem.File,
)::Cint
    return reinterpret(Cint, Base.fd(file))
end

function _metamdbg_require_single_path_component(
        component::AbstractString,
)::String
    normalized_component = String(component)
    normalized_component == basename(normalized_component) &&
        normalized_component ∉ ("", ".", "..") || throw(ArgumentError(
        "metaMDBG descriptor-relative path must be one non-special basename.",
    ))
    return normalized_component
end

function _metamdbg_openat(
        directory::Base.Filesystem.File,
        component::AbstractString,
        flags::Integer;
        mode::Union{Nothing, Integer} = nothing,
)::Base.Filesystem.File
    normalized_component =
        _metamdbg_require_single_path_component(component)
    descriptor = if mode === nothing
        ccall(
            :openat,
            Cint,
            (Cint, Cstring, Cint),
            _metamdbg_file_descriptor_number(directory),
            normalized_component,
            Cint(flags),
        )
    else
        ccall(
            :openat,
            Cint,
            (Cint, Cstring, Cint, Base.Cmode_t),
            _metamdbg_file_descriptor_number(directory),
            normalized_component,
            Cint(flags),
            Base.Cmode_t(mode),
        )
    end
    if descriptor < 0
        saved_errno = Base.Libc.errno()
        throw(SystemError(
            "openat $(repr(normalized_component))",
            saved_errno,
        ))
    end
    return Base.Filesystem.File(Base.RawFD(descriptor))
end

function _metamdbg_openat_entry_identity(
        directory::Base.Filesystem.File,
        component::AbstractString,
        path::AbstractString,
        expect_directory::Bool,
)::NamedTuple
    entry = _metamdbg_openat(
        directory,
        component,
        Base.JL_O_RDONLY |
        Base.JL_O_NONBLOCK |
        Base.JL_O_NOFOLLOW |
        Base.JL_O_CLOEXEC,
    )
    try
        entry_status = stat(entry)
        type_matches = expect_directory ?
                       isdir(entry_status) : isfile(entry_status)
        type_matches || error(
            "metaMDBG descriptor-relative durable entry changed type: " *
            "$(path).",
        )
        return (;
            path = normpath(abspath(String(path))),
            device = entry_status.device,
            inode = entry_status.inode,
        )
    finally
        close(entry)
    end
end

function _metamdbg_openat_entry_exists(
        directory::Base.Filesystem.File,
        component::AbstractString,
)::Bool
    normalized_component =
        _metamdbg_require_single_path_component(component)
    descriptor = ccall(
        :openat,
        Cint,
        (Cint, Cstring, Cint),
        _metamdbg_file_descriptor_number(directory),
        normalized_component,
        Cint(
            Base.JL_O_RDONLY |
            Base.JL_O_NONBLOCK |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
        ),
    )
    saved_errno = descriptor < 0 ? Base.Libc.errno() : Cint(0)
    if descriptor >= 0
        close(Base.Filesystem.File(Base.RawFD(descriptor)))
        return true
    end
    saved_errno == Base.Libc.ENOENT && return false
    return true
end

function _metamdbg_mkdirat(
        directory::Base.Filesystem.File,
        component::AbstractString,
        mode::Integer,
)::Nothing
    normalized_component =
        _metamdbg_require_single_path_component(component)
    result = ccall(
        :mkdirat,
        Cint,
        (Cint, Cstring, Base.Cmode_t),
        _metamdbg_file_descriptor_number(directory),
        normalized_component,
        Base.Cmode_t(mode),
    )
    if result != 0
        saved_errno = Base.Libc.errno()
        throw(SystemError(
            "mkdirat $(repr(normalized_component))",
            saved_errno,
        ))
    end
    return nothing
end

function _metamdbg_linux_renameat2_syscall_number(
        architecture::Symbol = Sys.ARCH,
)::Clong
    number = if architecture == :x86_64
        316
    elseif architecture in (:aarch64, :riscv64)
        276
    elseif architecture == :i686
        353
    elseif architecture == :armv7l
        382
    elseif architecture in (:powerpc64le, :ppc64le)
        357
    elseif architecture == :s390x
        347
    else
        error(
            "metaMDBG has no verified Linux renameat2 syscall number for " *
            "architecture $(architecture).",
        )
    end
    return Clong(number)
end

function _metamdbg_linux_renameat2_libc(
        source_descriptor::Cint,
        source::AbstractString,
        destination_descriptor::Cint,
        destination::AbstractString,
        flags::Cuint,
)::Cint
    return ccall(
        :renameat2,
        Cint,
        (Cint, Cstring, Cint, Cstring, Cuint),
        source_descriptor,
        String(source),
        destination_descriptor,
        String(destination),
        flags,
    )
end

function _metamdbg_linux_renameat2_syscall(
        syscall_number::Clong,
        source_descriptor::Cint,
        source::AbstractString,
        destination_descriptor::Cint,
        destination::AbstractString,
        flags::Cuint,
)::Clong
    return ccall(
        :syscall,
        Clong,
        (Clong, Cint, Cstring, Cint, Cstring, Cuint),
        syscall_number,
        source_descriptor,
        String(source),
        destination_descriptor,
        String(destination),
        flags,
    )
end

function _metamdbg_linux_renameat2(
        source_descriptor::Cint,
        source::AbstractString,
        destination_descriptor::Cint,
        destination::AbstractString,
        flags::Cuint;
        libc_runner::Function = _metamdbg_linux_renameat2_libc,
        syscall_runner::Function = _metamdbg_linux_renameat2_syscall,
)::Cint
    try
        return libc_runner(
            source_descriptor,
            String(source),
            destination_descriptor,
            String(destination),
            flags,
        )
    catch caught
        caught isa InterruptException && rethrow()
        occursin("renameat2", sprint(showerror, caught)) || rethrow()
    end
    result = syscall_runner(
        _metamdbg_linux_renameat2_syscall_number(),
        source_descriptor,
        String(source),
        destination_descriptor,
        String(destination),
        flags,
    )
    return Cint(result)
end

function _metamdbg_renameat_noreplace(
        source_directory::Base.Filesystem.File,
        source_component::AbstractString,
        destination_directory::Base.Filesystem.File,
        destination_component::AbstractString,
)::Nothing
    normalized_source =
        _metamdbg_require_single_path_component(source_component)
    normalized_destination =
        _metamdbg_require_single_path_component(destination_component)
    source_descriptor =
        _metamdbg_file_descriptor_number(source_directory)
    destination_descriptor =
        _metamdbg_file_descriptor_number(destination_directory)
    result = if Sys.isapple()
        rename_exclusive = Cuint(0x00000004)
        ccall(
            :renameatx_np,
            Cint,
            (Cint, Cstring, Cint, Cstring, Cuint),
            source_descriptor,
            normalized_source,
            destination_descriptor,
            normalized_destination,
            rename_exclusive,
        )
    elseif Sys.islinux()
        rename_noreplace = Cuint(0x00000001)
        _metamdbg_linux_renameat2(
            source_descriptor,
            normalized_source,
            destination_descriptor,
            normalized_destination,
            rename_noreplace,
        )
    else
        error(
            "metaMDBG exact durable removal requires an operating-system " *
            "descriptor-relative no-replace rename primitive; refusing " *
            "cleanup on $(Sys.KERNEL).",
        )
    end
    if result != 0
        saved_errno = Base.Libc.errno()
        throw(SystemError(
            "no-replace renameat $(repr(normalized_source)) to " *
            repr(normalized_destination),
            saved_errno,
        ))
    end
    return nothing
end

function _metamdbg_unlinkat(
        directory::Base.Filesystem.File,
        component::AbstractString;
        directory_entry::Bool = false,
)::Nothing
    normalized_component =
        _metamdbg_require_single_path_component(component)
    remove_directory_flag = directory_entry ?
                            (Sys.isapple() ? 0x0080 : 0x0200) : 0
    result = ccall(
        :unlinkat,
        Cint,
        (Cint, Cstring, Cint),
        _metamdbg_file_descriptor_number(directory),
        normalized_component,
        Cint(remove_directory_flag),
    )
    if result != 0
        saved_errno = Base.Libc.errno()
        throw(SystemError(
            "unlinkat $(repr(normalized_component))",
            saved_errno,
        ))
    end
    return nothing
end

function _metamdbg_descriptor_directory_entries(
        directory::Base.Filesystem.File,
)::Vector{String}
    descriptor = _metamdbg_file_descriptor_number(directory)
    independent_descriptor = ccall(
        :openat,
        Cint,
        (Cint, Cstring, Cint),
        descriptor,
        ".",
        Cint(
            Base.JL_O_RDONLY |
            Base.JL_O_DIRECTORY |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
        ),
    )
    if independent_descriptor < 0
        saved_errno = Base.Libc.errno()
        throw(SystemError(
            "openat independent metaMDBG directory scan descriptor",
            saved_errno,
        ))
    end
    directory_stream = ccall(
        :fdopendir,
        Ptr{Cvoid},
        (Cint,),
        independent_descriptor,
    )
    if directory_stream == C_NULL
        saved_errno = Base.Libc.errno()
        ccall(:close, Cint, (Cint,), independent_descriptor)
        throw(SystemError(
            "fdopendir metaMDBG quarantine directory descriptor",
            saved_errno,
        ))
    end
    ccall(:rewinddir, Cvoid, (Ptr{Cvoid},), directory_stream)
    name_offset = if Sys.isapple()
        21
    elseif Sys.WORD_SIZE == 64
        19
    else
        11
    end
    entries = String[]
    try
        while true
            Base.Libc.errno(0)
            entry = ccall(
                :readdir,
                Ptr{UInt8},
                (Ptr{Cvoid},),
                directory_stream,
            )
            if entry == C_NULL
                saved_errno = Base.Libc.errno()
                saved_errno == 0 || throw(SystemError(
                    "readdir metaMDBG quarantine directory descriptor",
                    saved_errno,
                ))
                break
            end
            name = unsafe_string(entry + name_offset)
            name in (".", "..") && continue
            push!(
                entries,
                _metamdbg_require_single_path_component(name),
            )
        end
    finally
        close_result =
            ccall(:closedir, Cint, (Ptr{Cvoid},), directory_stream)
        if close_result != 0
            saved_errno = Base.Libc.errno()
            throw(SystemError(
                "closedir metaMDBG quarantine directory descriptor",
                saved_errno,
            ))
        end
    end
    sort!(entries)
    return entries
end

function _remove_metamdbg_quarantined_directory_contents!(
        directory::Base.Filesystem.File,
        label::AbstractString,
)::Nothing
    for _attempt in 1:3
        entries = _metamdbg_descriptor_directory_entries(directory)
        isempty(entries) && return nothing
        for entry_name in entries
            result = ccall(
                :unlinkat,
                Cint,
                (Cint, Cstring, Cint),
                _metamdbg_file_descriptor_number(directory),
                entry_name,
                Cint(0),
            )
            result == 0 && continue
            child = try
                _metamdbg_openat(
                    directory,
                    entry_name,
                    Base.JL_O_RDONLY |
                    Base.JL_O_DIRECTORY |
                    Base.JL_O_NOFOLLOW |
                    Base.JL_O_CLOEXEC,
                )
            catch caught
                error(
                    "metaMDBG could not remove quarantined entry " *
                    "$(repr(entry_name)) from $(label): " *
                    sprint(showerror, caught),
                )
            end
            try
                _remove_metamdbg_quarantined_directory_contents!(
                    child,
                    "$(label)/$(entry_name)",
                )
                _fsync_metamdbg_descriptor(
                    Base.fd(child),
                    "$(label)/$(entry_name)",
                )
            finally
                close(child)
            end
            _metamdbg_unlinkat(
                directory,
                entry_name;
                directory_entry = true,
            )
        end
        _fsync_metamdbg_descriptor(Base.fd(directory), label)
    end
    isempty(_metamdbg_descriptor_directory_entries(directory)) || error(
        "metaMDBG quarantined directory remained nonempty after three " *
        "descriptor-relative cleanup passes: $(label).",
    )
    return nothing
end

function _metamdbg_removal_quarantine_name(
        path::AbstractString,
)::String
    return basename(String(path)) * ".removing." *
           bytes2hex(Random.rand(UInt8, 16))
end

function _is_metamdbg_removal_quarantine_name(
        name::AbstractString,
)::Bool
    return occursin(r"\.removing\.[0-9a-f]{32}$", String(name))
end

function _is_metamdbg_parent_removal_quarantine_for_output_root(
        name::AbstractString,
        canonical_outdir::AbstractString,
)::Bool
    normalized_name = String(name)
    _is_metamdbg_removal_quarantine_name(normalized_name) || return false
    source_name = replace(
        normalized_name,
        r"\.removing\.[0-9a-f]{32}$" => "",
    )
    exact_names = Set(String[
        basename(String(canonical_outdir)),
        basename(_metamdbg_output_lock_path(canonical_outdir)),
        basename(_output_root_reservation_lock_path_from_canonical(
            canonical_outdir,
        )),
        basename(_metamdbg_lifecycle_cleanup_reservation_path(
            canonical_outdir,
        )),
    ])
    source_name in exact_names && return true
    prefixes = String[
        _metamdbg_submission_reservation_prefix(canonical_outdir),
        _metamdbg_pending_submission_job_prefix(canonical_outdir),
        _output_root_durable_reservation_prefix_from_canonical(
            canonical_outdir,
        ),
    ]
    return any(prefix -> startswith(source_name, prefix), prefixes)
end

function _metamdbg_removal_quarantine_evidence(
        outdir::AbstractString;
        post_reservation_scan_hook::Function =
            (_path::AbstractString, _identity::NamedTuple) -> nothing,
)::Vector{String}
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    parent_path = dirname(canonical_outdir)
    if !ispath(parent_path)
        islink(parent_path) && error(
            "metaMDBG removal-evidence parent is a dangling symbolic link: " *
            "$(parent_path).",
        )
        return String[]
    end
    parent_identity = _metamdbg_directory_path_identity(parent_path)
    parent = Base.Filesystem.open(
        parent_identity.path,
        Base.JL_O_RDONLY |
        Base.JL_O_DIRECTORY |
        Base.JL_O_NOFOLLOW |
        Base.JL_O_CLOEXEC,
    )
    evidence = String[]
    outdir_identity = nothing
    outdir_removal_entries = String[]
    try
        _require_unchanged_metamdbg_directory_descriptor(
            parent,
            parent_identity,
            parent_identity.path,
        )
        entries = _metamdbg_descriptor_directory_entries(parent)
        for entry_name in entries
            if _is_metamdbg_parent_removal_quarantine_for_output_root(
                    entry_name,
                    canonical_outdir,
            )
                push!(evidence, joinpath(parent_identity.path, entry_name))
            end
        end

        if ispath(canonical_outdir) || islink(canonical_outdir)
            outdir_identity =
                _metamdbg_directory_path_identity(canonical_outdir)
            output_directory = Base.Filesystem.open(
                canonical_outdir,
                Base.JL_O_RDONLY |
                Base.JL_O_DIRECTORY |
                Base.JL_O_NONBLOCK |
                Base.JL_O_NOFOLLOW |
                Base.JL_O_CLOEXEC,
            )
            try
                _require_unchanged_metamdbg_directory_descriptor(
                    output_directory,
                    outdir_identity,
                    canonical_outdir,
                )
                outdir_entries =
                    _metamdbg_descriptor_directory_entries(output_directory)
                outdir_removal_entries = sort!(String[
                    nested_name for nested_name in outdir_entries if
                    _is_metamdbg_removal_quarantine_name(nested_name)
                ])
                for nested_name in outdir_removal_entries
                    push!(
                        evidence,
                        joinpath(canonical_outdir, nested_name),
                    )
                end
                _require_unchanged_metamdbg_directory_descriptor(
                    output_directory,
                    outdir_identity,
                    canonical_outdir,
                )
                _metamdbg_directory_path_identity(canonical_outdir) ==
                    outdir_identity || error(
                    "metaMDBG canonical output directory changed during " *
                    "removal-evidence scan: $(canonical_outdir).",
                )
            finally
                close(output_directory)
            end
        end

        reservation_prefix =
            _metamdbg_submission_reservation_prefix(canonical_outdir)
        relevant_parent_entries = sort!(String[
            entry_name for entry_name in entries if
            _is_metamdbg_parent_removal_quarantine_for_output_root(
                entry_name,
                canonical_outdir,
            ) ||
            (
                startswith(entry_name, reservation_prefix) &&
                !_is_metamdbg_removal_quarantine_name(entry_name)
            )
        ])
        for entry_name in entries
            startswith(entry_name, reservation_prefix) || continue
            _is_metamdbg_removal_quarantine_name(entry_name) && continue
            reservation = _metamdbg_openat(
                parent,
                entry_name,
                Base.JL_O_RDONLY |
                Base.JL_O_DIRECTORY |
                Base.JL_O_NONBLOCK |
                Base.JL_O_NOFOLLOW |
                Base.JL_O_CLOEXEC,
            )
            reservation_path = joinpath(parent_identity.path, entry_name)
            reservation_status = stat(reservation)
            isdir(reservation_status) || error(
                "metaMDBG submission reservation changed type during " *
                "removal-evidence scan: $(reservation_path).",
            )
            reservation_identity = (;
                path = reservation_path,
                device = reservation_status.device,
                inode = reservation_status.inode,
            )
            try
                nested_entries =
                    _metamdbg_descriptor_directory_entries(reservation)
                for nested_name in nested_entries
                    _is_metamdbg_removal_quarantine_name(nested_name) ||
                        continue
                    push!(
                        evidence,
                        joinpath(
                            parent_identity.path,
                            entry_name,
                            nested_name,
                        ),
                    )
                end
                post_reservation_scan_hook(
                    reservation_path,
                    reservation_identity,
                )
                _require_unchanged_metamdbg_directory_descriptor(
                    reservation,
                    reservation_identity,
                    reservation_path,
                )
                observed_reservation_identity =
                    _metamdbg_openat_entry_identity(
                        parent,
                        entry_name,
                        reservation_path,
                        true,
                    )
                _metamdbg_path_identity_matches(
                    observed_reservation_identity,
                    reservation_identity,
                ) || error(
                    "metaMDBG submission reservation changed during " *
                    "removal-evidence scan: $(reservation_path).",
                )
                _metamdbg_descriptor_directory_entries(reservation) ==
                    nested_entries || error(
                    "metaMDBG submission reservation entry set changed " *
                    "during removal-evidence scan: $(reservation_path).",
                )
            finally
                close(reservation)
            end
        end
        final_parent_entries = _metamdbg_descriptor_directory_entries(parent)
        final_relevant_parent_entries = sort!(String[
            entry_name for entry_name in final_parent_entries if
            _is_metamdbg_parent_removal_quarantine_for_output_root(
                entry_name,
                canonical_outdir,
            ) ||
            (
                startswith(entry_name, reservation_prefix) &&
                !_is_metamdbg_removal_quarantine_name(entry_name)
            )
        ])
        final_relevant_parent_entries == relevant_parent_entries || error(
            "metaMDBG output-root removal-evidence entry set changed " *
            "during restart scan: $(parent_identity.path).",
        )
        if outdir_identity === nothing
            if ispath(canonical_outdir) || islink(canonical_outdir)
                error(
                    "metaMDBG canonical output directory appeared during " *
                    "removal-evidence scan: $(canonical_outdir).",
                )
            end
        else
            output_directory = _metamdbg_openat(
                parent,
                basename(canonical_outdir),
                Base.JL_O_RDONLY |
                Base.JL_O_DIRECTORY |
                Base.JL_O_NONBLOCK |
                Base.JL_O_NOFOLLOW |
                Base.JL_O_CLOEXEC,
            )
            try
                _require_unchanged_metamdbg_directory_descriptor(
                    output_directory,
                    outdir_identity,
                    canonical_outdir,
                )
                final_outdir_removal_entries = sort!(String[
                    nested_name for nested_name in
                    _metamdbg_descriptor_directory_entries(output_directory) if
                    _is_metamdbg_removal_quarantine_name(nested_name)
                ])
                final_outdir_removal_entries ==
                    outdir_removal_entries || error(
                    "metaMDBG canonical output directory removal-evidence " *
                    "entry set changed during scan: $(canonical_outdir).",
                )
            finally
                close(output_directory)
            end
        end
        _require_unchanged_metamdbg_directory_descriptor(
            parent,
            parent_identity,
            parent_identity.path,
        )
        _metamdbg_directory_path_identity(parent_identity.path) ==
            parent_identity || error(
            "metaMDBG removal-evidence parent changed during restart scan: " *
            "$(parent_identity.path).",
        )
    finally
        close(parent)
    end
    return sort!(unique!(evidence))
end

function _require_no_metamdbg_removal_quarantine_evidence!(
        outdir::AbstractString,
)::Nothing
    evidence = _metamdbg_removal_quarantine_evidence(outdir)
    isempty(evidence) || error(
        "metaMDBG found fail-closed .removing evidence from an interrupted " *
        "exact durable cleanup: $(join(evidence, ", ")). Inspect the bound " *
        "manifest and payload before manually resolving it; refusing " *
        "automatic restart or mutation.",
    )
    return nothing
end

function _create_metamdbg_removal_quarantine!(
        parent::Base.Filesystem.File,
        parent_identity::NamedTuple,
        path::AbstractString,
        expected_identity::NamedTuple,
        entry_kind::Symbol,
)::NamedTuple
    quarantine_name = _metamdbg_removal_quarantine_name(path)
    _metamdbg_mkdirat(parent, quarantine_name, 0o700)
    quarantine_path = joinpath(parent_identity.path, quarantine_name)
    quarantine = nothing
    try
        quarantine = _metamdbg_openat(
            parent,
            quarantine_name,
            Base.JL_O_RDONLY |
            Base.JL_O_DIRECTORY |
            Base.JL_O_NONBLOCK |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
        )
        quarantine_status = stat(quarantine)
        quarantine_status.uid == Base.Libc.getuid() || error(
            "metaMDBG removal quarantine is not owned by the current user: " *
            "$(quarantine_path).",
        )
        (quarantine_status.mode & 0o777) == 0o700 || error(
            "metaMDBG removal quarantine must have mode 0700: " *
            "$(quarantine_path).",
        )
        quarantine_identity = (;
            path = quarantine_path,
            device = quarantine_status.device,
            inode = quarantine_status.inode,
        )
        manifest_name = "manifest.json"
        manifest = _metamdbg_openat(
            quarantine,
            manifest_name,
            Base.JL_O_WRONLY |
            Base.JL_O_CREAT |
            Base.JL_O_EXCL |
            Base.JL_O_CLOEXEC;
            mode = 0o600,
        )
        try
            write(manifest, JSON.json((;
                schema_version = 1,
                operation = "remove-exact-metamdbg-durable-entry",
                target_path = normpath(abspath(String(path))),
                target_kind = String(entry_kind),
                expected_device = expected_identity.device,
                expected_inode = expected_identity.inode,
            )) * "\n")
            _fsync_metamdbg_descriptor(
                Base.fd(manifest),
                joinpath(quarantine_path, manifest_name),
            )
        finally
            close(manifest)
        end
        _fsync_bound_metamdbg_directory(
            quarantine,
            quarantine_identity,
            quarantine_path,
        )
        _fsync_bound_metamdbg_directory(
            parent,
            parent_identity,
            parent_identity.path,
        )
        return (;
            quarantine,
            quarantine_name,
            quarantine_path,
            quarantine_identity,
            manifest_name,
        )
    catch caught
        if quarantine !== nothing && isopen(quarantine)
            close(quarantine)
        end
        try
            _fsync_bound_metamdbg_directory(
                parent,
                parent_identity,
                parent_identity.path,
            )
        catch sync_error
            caught isa InterruptException && rethrow()
            error(
                "metaMDBG removal-quarantine creation failed and its " *
                "fail-closed evidence could not be durably synchronized at " *
                "$(quarantine_path). Creation cause: " *
                "$(sprint(showerror, caught)). Sync cause: " *
                "$(sprint(showerror, sync_error))",
            )
        end
        caught isa InterruptException && rethrow()
        evidence_exists =
            _metamdbg_openat_entry_exists(parent, quarantine_name)
        evidence_state = evidence_exists ? "retained" : "not observable"
        error(
            "metaMDBG removal-quarantine creation failed; fail-closed " *
            "evidence is $(evidence_state) at $(quarantine_path). Cause: " *
            sprint(showerror, caught),
        )
    end
end

function _require_metamdbg_quarantined_payload_descriptor(
        payload::Base.Filesystem.File,
        payload_path::AbstractString,
        expected_identity::NamedTuple,
        entry_kind::Symbol,
)::Nothing
    payload_status = stat(payload)
    type_matches = entry_kind == :directory ?
                   isdir(payload_status) : isfile(payload_status)
    type_matches || error(
        "metaMDBG quarantined durable $(entry_kind) descriptor changed type: " *
        "$(payload_path).",
    )
    observed_identity = (;
        path = normpath(abspath(String(payload_path))),
        device = payload_status.device,
        inode = payload_status.inode,
    )
    _metamdbg_path_identity_matches(
        observed_identity,
        merge(expected_identity, (; path = observed_identity.path)),
    ) || error(
        "metaMDBG quarantined durable $(entry_kind) descriptor does not bind " *
        "the expected payload: $(payload_path).",
    )
    return nothing
end

function _remove_exact_metamdbg_durable_entry!(
        path::AbstractString,
        expected_identity::NamedTuple,
        entry_kind::Symbol;
        recursive::Bool = false,
        remover::Union{Nothing, Function} = nothing,
        pre_quarantine_rename_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        pre_quarantine_unlink_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        post_quarantine_validation_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        pre_quarantine_final_unlink_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        post_quarantine_destination_fsync_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        pre_source_parent_fsync_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        quarantined_payload_validator::Function =
            (_payload::Base.Filesystem.File, _path::AbstractString) -> nothing,
        pre_quarantine_directory_unlink_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
)::Nothing
    entry_kind in (:file, :directory) || throw(ArgumentError(
        "metaMDBG durable entry kind must be :file or :directory.",
    ))
    normalized_path = normpath(abspath(String(path)))
    expected_identity.path == normalized_path || error(
        "metaMDBG durable-entry removal identity targets another path: " *
        "$(normalized_path).",
    )
    parent_identity = _metamdbg_directory_path_identity(dirname(normalized_path))
    parent = Base.Filesystem.open(
        parent_identity.path,
        Base.JL_O_RDONLY |
        Base.JL_O_DIRECTORY |
        Base.JL_O_NOFOLLOW |
        Base.JL_O_CLOEXEC,
    )
    quarantine = nothing
    quarantine_path = nothing
    quarantine_removed = false
    payload = nothing
    source_name = basename(normalized_path)
    payload_name = "payload"
    try
        _require_unchanged_metamdbg_directory_descriptor(
            parent,
            parent_identity,
            parent_identity.path,
        )
        observed_identity = _metamdbg_openat_entry_identity(
            parent,
            source_name,
            normalized_path,
            entry_kind == :directory,
        )
        _metamdbg_path_identity_matches(observed_identity, expected_identity) ||
            error(
                "metaMDBG refuses to remove a replacement durable " *
                "$(entry_kind): $(normalized_path).",
            )
        if entry_kind == :directory && !recursive && remover === nothing
            source_directory = _metamdbg_openat(
                parent,
                source_name,
                Base.JL_O_RDONLY |
                Base.JL_O_DIRECTORY |
                Base.JL_O_NONBLOCK |
                Base.JL_O_NOFOLLOW |
                Base.JL_O_CLOEXEC,
            )
            try
                _require_metamdbg_quarantined_payload_descriptor(
                    source_directory,
                    normalized_path,
                    expected_identity,
                    entry_kind,
                )
                isempty(_metamdbg_descriptor_directory_entries(
                    source_directory,
                )) || error(
                    "metaMDBG refuses nonrecursive removal of a nonempty " *
                    "durable directory: $(normalized_path).",
                )
            finally
                close(source_directory)
            end
        end
        quarantine = _create_metamdbg_removal_quarantine!(
            parent,
            parent_identity,
            normalized_path,
            expected_identity,
            entry_kind,
        )
        quarantine_path = quarantine.quarantine_path
        payload_path = joinpath(quarantine_path, payload_name)
        pre_quarantine_rename_hook(normalized_path, quarantine_path)
        _metamdbg_renameat_noreplace(
            parent,
            source_name,
            quarantine.quarantine,
            payload_name,
        )
        moved_identity = _metamdbg_openat_entry_identity(
            quarantine.quarantine,
            payload_name,
            payload_path,
            entry_kind == :directory,
        )
        _metamdbg_path_identity_matches(
            moved_identity,
            merge(expected_identity, (; path = payload_path)),
        ) || error(
            "metaMDBG moved a replacement into its fail-closed removal " *
            "quarantine; the replacement was retained at $(payload_path).",
        )
        payload_flags = Base.JL_O_RDONLY |
                        Base.JL_O_NONBLOCK |
                        Base.JL_O_NOFOLLOW |
                        Base.JL_O_CLOEXEC
        entry_kind == :directory &&
            (payload_flags |= Base.JL_O_DIRECTORY)
        payload = _metamdbg_openat(
            quarantine.quarantine,
            payload_name,
            payload_flags,
        )
        _require_metamdbg_quarantined_payload_descriptor(
            payload,
            payload_path,
            expected_identity,
            entry_kind,
        )
        # A cross-directory rename becomes fail-closed only when its destination
        # name is durable before the source-name deletion is made durable.
        _fsync_bound_metamdbg_directory(
            quarantine.quarantine,
            quarantine.quarantine_identity,
            quarantine_path,
        )
        post_quarantine_destination_fsync_hook(
            normalized_path,
            quarantine_path,
        )
        pre_source_parent_fsync_hook(normalized_path, quarantine_path)
        _fsync_bound_metamdbg_directory(
            parent,
            parent_identity,
            parent_identity.path,
        )
        !_metamdbg_openat_entry_exists(parent, source_name) || error(
            "metaMDBG durable source name reappeared after quarantine; " *
            "retaining exact removal evidence at $(quarantine_path).",
        )
        pre_quarantine_unlink_hook(normalized_path, payload_path)
        moved_identity = _metamdbg_openat_entry_identity(
            quarantine.quarantine,
            payload_name,
            payload_path,
            entry_kind == :directory,
        )
        _metamdbg_path_identity_matches(
            moved_identity,
            merge(expected_identity, (; path = payload_path)),
        ) || error(
            "metaMDBG quarantined durable $(entry_kind) was replaced before " *
            "unlink; replacement retained at $(payload_path).",
        )
        _require_metamdbg_quarantined_payload_descriptor(
            payload,
            payload_path,
            expected_identity,
            entry_kind,
        )
        post_quarantine_validation_hook(normalized_path, payload_path)
        _require_metamdbg_quarantined_payload_descriptor(
            payload,
            payload_path,
            expected_identity,
            entry_kind,
        )
        post_validation_identity = _metamdbg_openat_entry_identity(
            quarantine.quarantine,
            payload_name,
            payload_path,
            entry_kind == :directory,
        )
        _metamdbg_path_identity_matches(
            post_validation_identity,
            merge(expected_identity, (; path = payload_path)),
        ) || error(
            "metaMDBG quarantined durable $(entry_kind) was replaced after " *
            "descriptor validation; replacement retained at $(payload_path).",
        )
        quarantined_payload_validator(payload, payload_path)
        if remover !== nothing
            if entry_kind == :directory
                remover(payload_path; recursive)
            else
                remover(payload_path)
            end
        elseif entry_kind == :directory && recursive
            _remove_metamdbg_quarantined_directory_contents!(
                payload,
                payload_path,
            )
            pre_quarantine_final_unlink_hook(
                normalized_path,
                payload_path,
            )
            _require_metamdbg_quarantined_payload_descriptor(
                payload,
                payload_path,
                expected_identity,
                entry_kind,
            )
            final_identity = _metamdbg_openat_entry_identity(
                quarantine.quarantine,
                payload_name,
                payload_path,
                true,
            )
            _metamdbg_path_identity_matches(
                final_identity,
                merge(expected_identity, (; path = payload_path)),
            ) || error(
                "metaMDBG quarantined durable directory was replaced before " *
                "final unlink; replacement retained at $(payload_path).",
            )
            quarantined_payload_validator(payload, payload_path)
            _metamdbg_unlinkat(
                quarantine.quarantine,
                payload_name;
                directory_entry = true,
            )
        else
            pre_quarantine_final_unlink_hook(
                normalized_path,
                payload_path,
            )
            _require_metamdbg_quarantined_payload_descriptor(
                payload,
                payload_path,
                expected_identity,
                entry_kind,
            )
            final_identity = _metamdbg_openat_entry_identity(
                quarantine.quarantine,
                payload_name,
                payload_path,
                entry_kind == :directory,
            )
            _metamdbg_path_identity_matches(
                final_identity,
                merge(expected_identity, (; path = payload_path)),
            ) || error(
                "metaMDBG quarantined durable $(entry_kind) was replaced " *
                "before final unlink; replacement retained at " *
                "$(payload_path).",
            )
            quarantined_payload_validator(payload, payload_path)
            _metamdbg_unlinkat(
                quarantine.quarantine,
                payload_name;
                directory_entry = entry_kind == :directory,
            )
        end
        !_metamdbg_openat_entry_exists(
            quarantine.quarantine,
            payload_name,
        ) || error(
            "metaMDBG quarantined durable $(entry_kind) remains after removal; " *
            "evidence retained at $(quarantine_path).",
        )
        !_metamdbg_openat_entry_exists(parent, source_name) || error(
            "metaMDBG durable source name reappeared after exact removal; " *
            "evidence retained at $(quarantine_path).",
        )
        _fsync_bound_metamdbg_directory(
            quarantine.quarantine,
            quarantine.quarantine_identity,
            quarantine_path,
        )
        _metamdbg_unlinkat(
            quarantine.quarantine,
            quarantine.manifest_name,
        )
        !_metamdbg_openat_entry_exists(
            quarantine.quarantine,
            quarantine.manifest_name,
        ) || error(
            "metaMDBG removal-quarantine manifest remains after unlink: " *
            "$(quarantine_path).",
        )
        _fsync_bound_metamdbg_directory(
            quarantine.quarantine,
            quarantine.quarantine_identity,
            quarantine_path,
        )
        pre_quarantine_directory_unlink_hook(
            normalized_path,
            quarantine_path,
        )
        _require_unchanged_metamdbg_directory_descriptor(
            quarantine.quarantine,
            quarantine.quarantine_identity,
            quarantine_path,
        )
        observed_quarantine_identity = _metamdbg_openat_entry_identity(
            parent,
            quarantine.quarantine_name,
            quarantine_path,
            true,
        )
        _metamdbg_path_identity_matches(
            observed_quarantine_identity,
            quarantine.quarantine_identity,
        ) || error(
            "metaMDBG removal quarantine was replaced before final unlink: " *
            "$(quarantine_path).",
        )
        close(quarantine.quarantine)
        _metamdbg_unlinkat(
            parent,
            quarantine.quarantine_name;
            directory_entry = true,
        )
        quarantine_removed = true
        !_metamdbg_openat_entry_exists(
            parent,
            quarantine.quarantine_name,
        ) || error(
            "metaMDBG removal quarantine reappeared after deletion: " *
            "$(quarantine_path).",
        )
        !_metamdbg_openat_entry_exists(parent, source_name) || error(
            "metaMDBG durable source name reappeared after quarantine cleanup: " *
            "$(normalized_path).",
        )
        _fsync_bound_metamdbg_directory(
            parent,
            parent_identity,
            parent_identity.path,
        )
        _metamdbg_directory_path_identity(parent_identity.path) ==
            parent_identity || error(
            "metaMDBG durable removal parent path changed during cleanup: " *
            "$(parent_identity.path).",
        )
    catch caught
        caught isa InterruptException && rethrow()
        if quarantine_path === nothing || quarantine_removed
            rethrow()
        end
        error(
            "metaMDBG retained fail-closed durable removal evidence at " *
            "$(quarantine_path). Cause: $(sprint(showerror, caught))",
        )
    finally
        if payload !== nothing && isopen(payload)
            close(payload)
        end
        if quarantine !== nothing && isopen(quarantine.quarantine)
            close(quarantine.quarantine)
        end
        close(parent)
    end
    return nothing
end

function _remove_exact_metamdbg_durable_file!(
        path::AbstractString,
        expected_identity::NamedTuple;
        remover::Union{Nothing, Function} = nothing,
        pre_quarantine_rename_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        pre_quarantine_unlink_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        post_quarantine_validation_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        pre_quarantine_final_unlink_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        post_quarantine_destination_fsync_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        pre_source_parent_fsync_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        quarantined_payload_validator::Function =
            (_payload::Base.Filesystem.File, _path::AbstractString) -> nothing,
        pre_quarantine_directory_unlink_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
)::Nothing
    return _remove_exact_metamdbg_durable_entry!(
        path,
        expected_identity,
        :file;
        remover,
        pre_quarantine_rename_hook,
        pre_quarantine_unlink_hook,
        post_quarantine_validation_hook,
        pre_quarantine_final_unlink_hook,
        post_quarantine_destination_fsync_hook,
        pre_source_parent_fsync_hook,
        quarantined_payload_validator,
        pre_quarantine_directory_unlink_hook,
    )
end

function _remove_exact_metamdbg_durable_directory!(
        path::AbstractString,
        expected_identity::NamedTuple;
        recursive::Bool = false,
        remover::Union{Nothing, Function} = nothing,
        pre_quarantine_rename_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        pre_quarantine_unlink_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        post_quarantine_validation_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        pre_quarantine_final_unlink_hook::Function =
            (_path::AbstractString, _payload::AbstractString) -> nothing,
        post_quarantine_destination_fsync_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        pre_source_parent_fsync_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
        quarantined_payload_validator::Function =
            (_payload::Base.Filesystem.File, _path::AbstractString) -> nothing,
        pre_quarantine_directory_unlink_hook::Function =
            (_path::AbstractString, _quarantine::AbstractString) -> nothing,
)::Nothing
    return _remove_exact_metamdbg_durable_entry!(
        path,
        expected_identity,
        :directory;
        recursive,
        remover,
        pre_quarantine_rename_hook,
        pre_quarantine_unlink_hook,
        post_quarantine_validation_hook,
        pre_quarantine_final_unlink_hook,
        post_quarantine_destination_fsync_hook,
        pre_source_parent_fsync_hook,
        quarantined_payload_validator,
        pre_quarantine_directory_unlink_hook,
    )
end

function _require_metamdbg_output_root_reservation_marker!(
        reservation::NamedTuple,
)::String
    marker = reservation.output_root_reservation_marker
    if !isfile(marker) || islink(marker) || filesize(marker) == 0
        error(
            "metaMDBG shared output-root reservation must be a nonempty, " *
            "regular, non-symlink file: $(marker).",
        )
    end
    marker_status = stat(marker)
    marker_status.uid == Base.Libc.getuid() || error(
        "metaMDBG shared output-root reservation is not owned by the " *
        "current user: $(marker).",
    )
    (marker_status.mode & 0o777) == 0o600 || error(
        "metaMDBG shared output-root reservation must have mode 0600: " *
        "$(marker).",
    )
    read(marker, String) == reservation.output_root_reservation_contents ||
        error(
            "metaMDBG shared output-root reservation does not match this " *
            "exact owner capability: $(marker).",
        )
    return marker
end

function _require_metamdbg_runtime_output_root_reservation_marker!(
        reservation::NamedTuple,
)::String
    marker = reservation.runtime_output_root_reservation_marker
    if !isdir(marker) || islink(marker)
        error(
            "metaMDBG runtime output-root reservation must be a regular, " *
            "non-symlink directory: $(marker).",
        )
    end
    marker_status = stat(marker)
    marker_status.uid == Base.Libc.getuid() || error(
        "metaMDBG runtime output-root reservation is not owned by the " *
        "current user: $(marker).",
    )
    (marker_status.mode & 0o777) == 0o700 || error(
        "metaMDBG runtime output-root reservation must have mode 0700: " *
        "$(marker).",
    )
    isempty(readdir(marker)) || error(
        "metaMDBG runtime output-root reservation contains unexpected " *
        "entries: $(marker).",
    )
    return marker
end

function _publish_metamdbg_output_root_reservation_marker!(
        reservation::NamedTuple,
)::String
    marker = reservation.output_root_reservation_marker
    if ispath(marker) || islink(marker)
        error(
            "metaMDBG refuses to overwrite an existing shared output-root " *
            "reservation: $(marker).",
        )
    end
    temporary_path, temporary_io = mktemp(dirname(marker))
    temporary_identity = _metamdbg_regular_file_identity(temporary_path)
    published = false
    try
        write(temporary_io, reservation.output_root_reservation_contents)
        chmod(temporary_path, 0o600)
        _fsync_metamdbg_file(temporary_io, temporary_path)
        close(temporary_io)
        Base.Filesystem.hardlink(temporary_path, marker)
        published = true
        _fsync_metamdbg_directory(dirname(marker))
        _remove_exact_metamdbg_durable_file!(
            temporary_path,
            temporary_identity,
        )
        _require_metamdbg_output_root_reservation_marker!(reservation)
    catch primary_error
        isopen(temporary_io) && close(temporary_io)
        if _output_root_path_entry_exists(temporary_path)
            try
                _remove_exact_metamdbg_durable_file!(
                    temporary_path,
                    temporary_identity,
                )
            catch cleanup_error
                @warn "metaMDBG failed to clean an unpublished shared " *
                      "output-root marker while preserving the primary " *
                      "publication failure" temporary_path primary_error cleanup_error
            end
        end
        if published && (ispath(marker) || islink(marker))
            try
                _remove_exact_metamdbg_durable_file!(
                    marker,
                    merge(temporary_identity, (; path = marker)),
                )
            catch cleanup_error
                @warn "metaMDBG failed to remove an invalid newly published " *
                      "shared output-root marker while preserving the " *
                      "primary publication failure" marker primary_error cleanup_error
            end
        end
        Base.rethrow()
    end
    return marker
end

function _require_no_active_metamdbg_submission_reservation!(
        outdir::AbstractString,
)::Nothing
    _require_no_metamdbg_removal_quarantine_evidence!(outdir)
    active_reservations = _metamdbg_submission_reservation_paths(outdir)
    isempty(active_reservations) || error(
        "metaMDBG output has an active nonlocal submission reservation: " *
        "$(join(active_reservations, ", ")). Refusing competing execution.",
    )
    pending_job_records = _metamdbg_pending_submission_job_paths(outdir)
    isempty(pending_job_records) || error(
        "metaMDBG output has pending scheduler job evidence: " *
        "$(join(pending_job_records, ", ")). Recover the exact accepted job " *
        "before competing execution.",
    )
    return nothing
end

function _require_metamdbg_private_submission_reservation!(
        reservation::NamedTuple,
)::NamedTuple
    reservation_path = reservation.path
    if !isdir(reservation_path) || islink(reservation_path)
        error(
            "metaMDBG submission reservation must be a regular, non-symlink " *
            "directory: $(reservation_path).",
        )
    end
    reservation_status = stat(reservation_path)
    reservation_status.uid == Base.Libc.getuid() || error(
        "metaMDBG submission reservation is not owned by the current user: " *
        "$(reservation_path).",
    )
    (reservation_status.mode & 0o777) == 0o700 || error(
        "metaMDBG submission reservation directory must have mode 0700: " *
        "$(reservation_path).",
    )
    has_bound_job = reservation.job_id isa AbstractString
    expected_entries = has_bound_job ?
                       String[
        _METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
        _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME,
    ] :
                       String[
        _METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
    ]
    readdir(reservation_path) == expected_entries || error(
        "metaMDBG submission reservation contains unexpected entries: " *
        "$(reservation_path).",
    )
    marker = reservation.contract_marker
    if !isfile(marker) || islink(marker) || filesize(marker) == 0
        error(
            "metaMDBG submission reservation contract must be a nonempty, " *
            "regular, non-symlink file: $(marker).",
        )
    end
    marker_status = stat(marker)
    marker_status.uid == Base.Libc.getuid() || error(
        "metaMDBG submission reservation contract is not owned by the " *
        "current user: $(marker).",
    )
    (marker_status.mode & 0o777) == 0o600 || error(
        "metaMDBG submission reservation contract must have mode 0600: " *
        "$(marker).",
    )
    read(marker, String) == reservation.contents || error(
        "metaMDBG submission reservation contract does not match this " *
        "invocation: $(marker).",
    )
    if has_bound_job
        job_marker = reservation.job_marker
        if !isfile(job_marker) || islink(job_marker) || filesize(job_marker) == 0
            error(
                "metaMDBG submission job record must be a nonempty, regular, " *
                "non-symlink file: $(job_marker).",
            )
        end
        job_status = stat(job_marker)
        job_status.uid == Base.Libc.getuid() || error(
            "metaMDBG submission job record is not owned by the current user: " *
            "$(job_marker).",
        )
        (job_status.mode & 0o777) == 0o600 || error(
            "metaMDBG submission job record must have mode 0600: " *
            "$(job_marker).",
        )
        read(job_marker, String) == reservation.job_contents || error(
            "metaMDBG submission job record does not match this reservation: " *
            "$(job_marker).",
        )
    end
    return reservation
end

function _require_metamdbg_submission_reservation!(
        reservation::NamedTuple,
)::NamedTuple
    _require_metamdbg_output_root_reservation_marker!(reservation)
    return _require_metamdbg_private_submission_reservation!(reservation)
end

function _metamdbg_submission_reservation_identity(
        reservation::NamedTuple,
)::NamedTuple
    _require_metamdbg_private_submission_reservation!(reservation)
    metadata = stat(reservation.path)
    return (; device = metadata.device, inode = metadata.inode)
end

function _metamdbg_shared_reservation_identity(
        reservation::NamedTuple,
)::NamedTuple
    marker = _require_metamdbg_output_root_reservation_marker!(reservation)
    metadata = stat(marker)
    return (; device = metadata.device, inode = metadata.inode)
end

function _metamdbg_runtime_shared_reservation_identity(
        reservation::NamedTuple,
)::NamedTuple
    marker =
        _require_metamdbg_runtime_output_root_reservation_marker!(reservation)
    metadata = stat(marker)
    return (; device = metadata.device, inode = metadata.inode)
end

function _require_unchanged_metamdbg_recovery_identities(
        reservation::NamedTuple,
        reservation_identity::NamedTuple,
        shared_reservation_identity::NamedTuple,
)::Nothing
    observed_reservation_identity =
        _metamdbg_submission_reservation_identity(reservation)
    observed_reservation_identity == reservation_identity || error(
        "metaMDBG submission reservation was replaced after recovery " *
        "inspection: $(reservation.path).",
    )
    observed_shared_identity = _metamdbg_shared_reservation_identity(reservation)
    observed_shared_identity == shared_reservation_identity || error(
        "metaMDBG shared output-root reservation was replaced after recovery " *
        "inspection: $(reservation.output_root_reservation_marker).",
    )
    return nothing
end

function _metamdbg_submission_reservation_from_path(
        path::AbstractString,
        canonical_outdir::AbstractString,
        ;
        allow_provisional::Bool = false,
        require_output_root_reservation::Bool = true,
)::NamedTuple
    normalized_path = normpath(abspath(path))
    path_state = _metamdbg_submission_reservation_path_state(
        normalized_path,
        canonical_outdir,
    )
    path_state == :legacy_consumed && error(
        "metaMDBG legacy consumed reservation remnants are not recovery " *
        "capabilities: $(normalized_path).",
    )
    if path_state == :provisional && !allow_provisional
        error(
            "metaMDBG provisional submission reservations require explicit " *
            "recovery inspection: $(normalized_path).",
        )
    end
    marker = joinpath(
        normalized_path,
        _METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
    )
    if !isdir(normalized_path) || islink(normalized_path)
        error(
            "metaMDBG submission reservation must be a regular, non-symlink " *
            "directory: $(normalized_path).",
        )
    end
    entries = readdir(normalized_path)
    allowed_entries = (
        String[_METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME],
        String[
            _METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
            _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME,
        ],
    )
    entries in allowed_entries || error(
        "metaMDBG submission reservation contains unexpected entries: " *
        "$(normalized_path).",
    )
    if !isfile(marker) || islink(marker) || filesize(marker) == 0
        error(
            "metaMDBG submission reservation contract must be a nonempty, " *
            "regular, non-symlink file: $(marker).",
        )
    end
    contents = read(marker, String)
    parsed = try
        JSON.parse(contents)
    catch caught
        error(
            "metaMDBG submission reservation contract is not valid JSON: " *
            "$(marker). Cause: $(sprint(showerror, caught))",
        )
    end
    parsed isa AbstractDict || error(
        "metaMDBG submission reservation contract is not a JSON object: " *
        "$(marker).",
    )
    expected_keys = Set(String[
        "schema_version",
        "canonical_outdir",
        "input_contract_signature",
        "graph_k",
        "workflow_signature",
        "scheduler_job_name",
        "owner_token",
    ])
    Set(String.(keys(parsed))) == expected_keys || error(
        "metaMDBG submission reservation contract has unexpected fields: " *
        "$(marker).",
    )
    schema_version = get(parsed, "schema_version", nothing)
    schema_version == _METAMDBG_SUBMISSION_RESERVATION_SCHEMA_VERSION || error(
        "metaMDBG submission reservation has an unsupported schema version: " *
        "$(repr(schema_version)).",
    )
    recorded_outdir = get(parsed, "canonical_outdir", nothing)
    recorded_outdir isa AbstractString || error(
        "metaMDBG submission reservation canonical_outdir must be a string.",
    )
    expected_outdir = _metamdbg_canonical_output_path(canonical_outdir)
    String(recorded_outdir) == expected_outdir || error(
        "metaMDBG submission reservation targets a different output root.",
    )
    input_contract_signature = get(
        parsed,
        "input_contract_signature",
        nothing,
    )
    input_contract_signature isa AbstractString &&
        occursin(r"^[0-9a-f]{64}$", input_contract_signature) || error(
        "metaMDBG submission reservation input contract signature is invalid.",
    )
    graph_k_value = get(parsed, "graph_k", nothing)
    graph_k_value isa Integer || error(
        "metaMDBG submission reservation graph_k must be an integer.",
    )
    graph_k = Int(graph_k_value)
    graph_k > 0 || error(
        "metaMDBG submission reservation graph_k must be positive.",
    )
    workflow_signature = get(parsed, "workflow_signature", nothing)
    workflow_signature isa AbstractString &&
        occursin(r"^[0-9a-f]{64}$", workflow_signature) || error(
        "metaMDBG submission reservation workflow signature is invalid.",
    )
    scheduler_job_name = get(parsed, "scheduler_job_name", nothing)
    scheduler_job_name isa AbstractString || error(
        "metaMDBG submission reservation scheduler_job_name must be a string.",
    )
    requested_job_name = _metamdbg_scheduler_job_name_prefix(
        scheduler_job_name,
        workflow_signature,
    )
    owner_token = get(parsed, "owner_token", nothing)
    owner_token isa AbstractString && !isempty(owner_token) || error(
        "metaMDBG submission reservation owner token is invalid.",
    )
    outputs = _metamdbg_output_paths(expected_outdir, graph_k)
    expected = _metamdbg_submission_reservation(
        outputs,
        (; signature = String(input_contract_signature)),
        graph_k;
        owner_token = String(owner_token),
        job_name = requested_job_name,
    )
    if path_state == :published
        expected.path == normalized_path || error(
            "metaMDBG submission reservation path does not match its workflow " *
            "signature.",
        )
        expected = merge(expected, (;
            publication_state = :published,
            reservation_state = :queued,
            lifecycle_owner = :submitter,
        ))
    elseif path_state == :provisional
        _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME in entries && error(
            "metaMDBG provisional pre-submit reservation must not contain a " *
            "scheduler job record: $(normalized_path).",
        )
        expected = merge(expected, (;
            path = normalized_path,
            contract_marker = marker,
            job_marker = joinpath(
                normalized_path,
                _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME,
            ),
            publication_state = :provisional,
            reservation_state = :provisional,
            lifecycle_owner = :submitter,
        ))
    elseif path_state in (:runtime, :consumed, :reclaiming)
        expected_path = if path_state == :runtime
            expected.runtime_path
        elseif path_state == :consumed
            expected.consumed_path
        else
            expected.reclaiming_path
        end
        expected_path == normalized_path || error(
            "metaMDBG $(path_state) reservation path does not match its " *
            "owner capability.",
        )
        if path_state in (:runtime, :consumed)
            _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME in entries || error(
                "metaMDBG $(path_state) reservation must retain its exact " *
                "scheduler job record: $(normalized_path).",
            )
        end
        expected = merge(expected, (;
            path = normalized_path,
            contract_marker = marker,
            job_marker = joinpath(
                normalized_path,
                _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME,
            ),
            publication_state = :published,
            reservation_state = path_state,
            lifecycle_owner = path_state == :reclaiming ?
                              :recovery : :runtime,
        ))
    else
        error(
            "metaMDBG cannot reconstruct reservation path state " *
            "$(repr(path_state)).",
        )
    end
    expected.workflow_signature == workflow_signature || error(
        "metaMDBG submission reservation workflow signature is inconsistent.",
    )
    expected.contents == contents || error(
        "metaMDBG submission reservation contract is not canonical.",
    )
    if _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME in entries
        if !isfile(expected.job_marker) || islink(expected.job_marker) ||
           filesize(expected.job_marker) == 0
            error(
                "metaMDBG submission job record must be a nonempty, regular, " *
                "non-symlink file: $(expected.job_marker).",
            )
        end
        job_contents = read(expected.job_marker, String)
        parsed_job = try
            JSON.parse(job_contents)
        catch caught
            error(
                "metaMDBG submission job record is not valid JSON: " *
                "$(expected.job_marker). Cause: $(sprint(showerror, caught))",
            )
        end
        parsed_job isa AbstractDict || error(
            "metaMDBG submission job record is not a JSON object: " *
            "$(expected.job_marker).",
        )
        job_schema_version = get(parsed_job, "schema_version", nothing)
        job_schema_version isa Integer && job_schema_version in (
            _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
            _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
        ) || error(
            "metaMDBG submission job record schema version is unsupported: " *
            "$(expected.job_marker).",
        )
        common_job_fields = String[
            "schema_version",
            "workflow_signature",
            "owner_token",
            "job_id",
        ]
        expected_job_fields = job_schema_version ==
                              _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION ?
                              Set(common_job_fields) :
                              Set(vcat(common_job_fields, ["job_cluster"]))
        Set(String.(keys(parsed_job))) == expected_job_fields || error(
            "metaMDBG submission job record has unexpected fields: " *
            "$(expected.job_marker).",
        )
        job_id = get(parsed_job, "job_id", nothing)
        job_id isa AbstractString || error(
            "metaMDBG submission job record job_id must be a string.",
        )
        job_cluster = job_schema_version ==
                      _METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION ?
                      nothing : get(parsed_job, "job_cluster", nothing)
        (job_cluster === nothing || job_cluster isa AbstractString) || error(
            "metaMDBG submission job record job_cluster must be a string or " *
            "null.",
        )
        expected = _metamdbg_bound_submission_reservation(
            expected,
            job_id,
            job_cluster,
            Int(job_schema_version),
        )
        expected.job_contents == job_contents || error(
            "metaMDBG submission job record is not canonical or does not " *
            "match its workflow owner.",
        )
    end
    if path_state in (:published, :provisional)
        if require_output_root_reservation
            _require_metamdbg_submission_reservation!(expected)
        else
            _require_metamdbg_private_submission_reservation!(expected)
        end
        runtime_exists = _output_root_path_entry_exists(
            expected.runtime_output_root_reservation_marker,
        )
        if runtime_exists
            _require_metamdbg_runtime_output_root_reservation_marker!(expected)
            return merge(expected, (;
                reservation_state = :runtime_claiming,
                lifecycle_owner = :runtime,
                submission_state = :runtime,
            ))
        end
        return expected
    elseif path_state == :runtime
        _require_metamdbg_private_submission_reservation!(expected)
        queued_exists = _output_root_path_entry_exists(
            expected.output_root_reservation_marker,
        )
        runtime_exists = _output_root_path_entry_exists(
            expected.runtime_output_root_reservation_marker,
        )
        queued_exists &&
            _require_metamdbg_output_root_reservation_marker!(expected)
        runtime_exists &&
            _require_metamdbg_runtime_output_root_reservation_marker!(expected)
        reservation_state = if queued_exists && runtime_exists
            :runtime_transition
        elseif runtime_exists
            :runtime
        elseif queued_exists
            :runtime_transition_ambiguous
        else
            :runtime_release_pending
        end
        return merge(expected, (;
            reservation_state,
            submission_state = :runtime,
        ))
    elseif path_state == :reclaiming
        _require_metamdbg_private_submission_reservation!(expected)
        queued_exists = _output_root_path_entry_exists(
            expected.output_root_reservation_marker,
        )
        runtime_exists = _output_root_path_entry_exists(
            expected.runtime_output_root_reservation_marker,
        )
        runtime_exists && error(
            "metaMDBG explicit reclaim transition unexpectedly has a runtime " *
            "shared marker.",
        )
        queued_exists &&
            _require_metamdbg_output_root_reservation_marker!(expected)
        return merge(expected, (;
            reservation_state = queued_exists ?
                                :reclaiming : :reclaim_release_pending,
        ))
    end
    _require_metamdbg_private_submission_reservation!(expected)
    queued_exists = _output_root_path_entry_exists(
        expected.output_root_reservation_marker,
    )
    runtime_exists = _output_root_path_entry_exists(
        expected.runtime_output_root_reservation_marker,
    )
    !queued_exists && !runtime_exists || error(
        "metaMDBG consumed reservation still has a queued or runtime shared " *
        "marker; refusing to treat it as durably consumed.",
    )
    return merge(expected, (;
        reservation_state = :consumed,
        submission_state = :consumed,
    ))
end

function _create_metamdbg_submission_reservation!(
        reservation::NamedTuple,
        outdir::AbstractString;
        pre_rename_hook::Function =
            (_reservation::NamedTuple) -> nothing,
        post_rename_hook::Function =
            (_reservation::NamedTuple) -> nothing,
)::NamedTuple
    existing_reservations =
        _metamdbg_submission_reservation_paths(outdir)
    if !isempty(existing_reservations)
        if existing_reservations == String[reservation.path]
            has_regular_marker = isfile(reservation.contract_marker) &&
                                 !islink(reservation.contract_marker)
            existing_contents = if has_regular_marker
                read(reservation.contract_marker, String)
            else
                nothing
            end
            if existing_contents == reservation.contents
                error(
                    "metaMDBG submission is already reserved by this " *
                    "invocation: $(reservation.path).",
                )
            end
            error(
                "metaMDBG submission is already reserved for this exact " *
                "workflow contract: $(reservation.path).",
            )
        end
        error(
            "metaMDBG output has a conflicting active submission " *
            "reservation: $(join(existing_reservations, ", ")).",
        )
    end

    reservation_parent = dirname(reservation.path)
    temporary_path = mktempdir(
        reservation_parent;
        prefix = _metamdbg_submission_reservation_prefix(outdir) * "tmp.",
    )
    chmod(temporary_path, 0o700)
    temporary_identity = _metamdbg_directory_path_identity(temporary_path)
    temporary_reservation = merge(reservation, (;
        path = temporary_path,
        contract_marker = joinpath(
            temporary_path,
            _METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
        ),
        job_marker = joinpath(
            temporary_path,
            _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME,
        ),
    ))
    output_root_marker_published = false
    output_root_marker_identity = nothing
    try
        open(temporary_reservation.contract_marker, "w") do output
            write(output, reservation.contents)
            chmod(temporary_reservation.contract_marker, 0o600)
            _fsync_metamdbg_file(
                output,
                temporary_reservation.contract_marker,
            )
        end
        _fsync_metamdbg_directory(temporary_path)
        _fsync_metamdbg_directory(reservation_parent)
        _publish_metamdbg_output_root_reservation_marker!(reservation)
        output_root_marker_published = true
        output_root_marker_identity = _metamdbg_regular_file_identity(
            reservation.output_root_reservation_marker,
        )
        _require_metamdbg_submission_reservation!(temporary_reservation)
        pre_rename_hook(temporary_reservation)
        mv(temporary_path, reservation.path)
        post_rename_hook(reservation)
        _fsync_metamdbg_directory(reservation_parent)
        _require_metamdbg_submission_reservation!(reservation)
    catch primary_error
        try
            if _output_root_path_entry_exists(reservation.path)
                _remove_metamdbg_submission_reservation!(reservation)
            else
                if _output_root_path_entry_exists(temporary_path)
                    _remove_exact_metamdbg_durable_directory!(
                        temporary_path,
                        temporary_identity;
                        recursive = true,
                    )
                end
                if output_root_marker_published &&
                   (ispath(reservation.output_root_reservation_marker) ||
                    islink(reservation.output_root_reservation_marker))
                    _remove_exact_metamdbg_durable_file!(
                        reservation.output_root_reservation_marker,
                        output_root_marker_identity,
                    )
                end
            end
        catch cleanup_error
            @warn "metaMDBG failed to clean an incomplete submission " *
                  "reservation while preserving the primary reservation " *
                  "failure" reservation primary_error cleanup_error
        end
        Base.rethrow()
    end
    return reservation
end

function _bind_metamdbg_submission_job!(
        reservation::NamedTuple,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString} = nothing,
        ;
        post_pending_job_record_creation_hook::Function =
            (_reservation::NamedTuple, _path::AbstractString) -> nothing,
        pre_job_record_publication_hook::Function =
            (_reservation::NamedTuple, _path::AbstractString) -> nothing,
        post_job_record_publication_hook::Function =
            (_reservation::NamedTuple) -> nothing,
        pending_job_record_remover::Union{Nothing, Function} = nothing,
)::NamedTuple
    _require_metamdbg_submission_reservation!(reservation)
    bound_reservation = _metamdbg_bound_submission_reservation(
        reservation,
        job_id,
        job_cluster,
    )
    marker = bound_reservation.job_marker
    if ispath(marker) || islink(marker)
        error(
            "metaMDBG refuses to overwrite an existing durable submission " *
            "job record: $(marker).",
        )
    end
    existing_pending_paths = _metamdbg_pending_submission_job_paths(
        reservation.canonical_outdir,
    )
    isempty(existing_pending_paths) || error(
        "metaMDBG refuses a new scheduler job binding while pending accepted-ID " *
        "evidence exists: $(join(existing_pending_paths, ", ")).",
    )
    pending_path = bound_reservation.pending_job_marker
    _output_root_path_entry_exists(pending_path) && error(
        "metaMDBG refuses to overwrite an existing pending scheduler job " *
        "record: $(pending_path).",
    )
    pending_io = nothing
    pending_record = nothing
    bound_owner = nothing
    committed = false
    try
        pending_io = Base.Filesystem.open(
            pending_path,
            Base.JL_O_WRONLY |
            Base.JL_O_CREAT |
            Base.JL_O_EXCL |
            Base.JL_O_CLOEXEC,
            0o600,
        )
        _fsync_metamdbg_descriptor(Base.fd(pending_io), pending_path)
        _fsync_metamdbg_directory(dirname(pending_path))
        post_pending_job_record_creation_hook(
            bound_reservation,
            pending_path,
        )
        write(pending_io, bound_reservation.job_contents)
        _fsync_metamdbg_descriptor(Base.fd(pending_io), pending_path)
        close(pending_io)
        pending_io = nothing
        _fsync_metamdbg_directory(dirname(pending_path))
        pending_record = _metamdbg_pending_submission_job_record(
            reservation,
            pending_path,
        )
        pending_record.job_id == bound_reservation.job_id || error(
            "metaMDBG pending scheduler job record lost its exact job ID.",
        )
        pending_record.job_cluster == bound_reservation.job_cluster || error(
            "metaMDBG pending scheduler job record lost its exact cluster.",
        )
        reservation_identity =
            _metamdbg_submission_reservation_identity(reservation)
        shared_reservation_identity =
            _metamdbg_shared_reservation_identity(reservation)
        pre_job_record_publication_hook(
            bound_reservation,
            pending_path,
        )
        _require_unchanged_metamdbg_recovery_identities(
            reservation,
            reservation_identity,
            shared_reservation_identity,
        )
        pending_record =
            _require_unchanged_metamdbg_pending_submission_job_record(
                reservation,
                pending_record,
            )
        bound_owner = _publish_metamdbg_pending_submission_job_record!(
            reservation,
            pending_record,
            reservation_identity,
            shared_reservation_identity,
        )
        committed = true
        post_job_record_publication_hook(bound_owner)
        _remove_metamdbg_pending_submission_job_record!(
            bound_owner,
            pending_record;
            pending_remover = pending_job_record_remover,
        )
        _require_metamdbg_submission_reservation!(bound_owner)
    catch primary_error
        if pending_io !== nothing && isopen(pending_io)
            close(pending_io)
        end
        if committed && pending_record isa NamedTuple &&
           bound_owner isa NamedTuple &&
           _output_root_path_entry_exists(pending_path)
            try
                _remove_metamdbg_pending_submission_job_record!(
                    bound_owner,
                    pending_record;
                    pending_remover = pending_job_record_remover,
                )
            catch cleanup_error
                @warn "metaMDBG failed to clean its committed pending " *
                      "scheduler job record while preserving the primary " *
                      "bind failure" pending_path primary_error cleanup_error
            end
        end
        Base.rethrow()
    end
    return bound_owner
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Bind a scheduler job ID to an inspected metaMDBG submission reservation.

This recovery operation is intended for the narrow crash window after the
scheduler accepted a job but before `run_metamdbg` durably recorded the returned
job identity. Set `confirm_submitted = true` only after independently confirming
that the exact scheduler job belongs to this reservation. For a federated Slurm
job, `job_cluster` is mandatory and the durable identity is the exact
`(job_id, job_cluster)` pair; leave it `nothing` only for an unscoped job.
Binding requires the fully published record and both filesystem identities
returned by inspection, so a
same-content private or shared replacement fails before `job.json` publication.
The job record is published atomically under the output-domain lock and becomes
available through `inspect_metamdbg_submission_reservations`. Binding holds the
full canonical local PID record, cleanup sentinel, and private lifecycle lock,
and first publishes an adjacent capability-and-job-ID-bound pending record.
Hard termination before or after durable `job.json` publication therefore
retains the exact accepted scheduler ID and cluster. Confirmed-dead inspection validates
and promotes pending-only evidence or verifies the exact hardlink before
durably removing the pending name. An ordinary precommit error retains pending
evidence for idempotent exact-job resume through this function. Once `job.json`
and its owner directory are fsynced, that commit is irreversible: a later
ordinary error preserves the exact bound owner and removes only the unchanged
pending inode when possible. A committed pending-cleanup remnant can be resumed
idempotently through the same API.
Complete either pending-recovery path before manually releasing the held Slurm
job; generated runtimes reject exact bound-job pending evidence before runtime
ownership publication and any other same-output pending evidence before owner
consumption.
"""
function bind_metamdbg_submission_reservation_job!(
        metadata::NamedTuple;
        owner_token::AbstractString,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString} = nothing,
        confirm_submitted::Bool = false,
)::NamedTuple
    confirm_submitted || throw(ArgumentError(
        "Set confirm_submitted=true only after independently confirming the " *
        "scheduler job belongs to this metaMDBG reservation.",
    ))
    required_fields = (
        :canonical_outdir,
        :path,
        :workflow_signature,
        :scheduler_job_name,
        :input_contract_signature,
        :graph_k,
        :owner_token,
        :job_id,
    )
    all(field -> hasproperty(metadata, field), required_fields) || throw(
        ArgumentError(
            "metaMDBG reservation metadata is incomplete; inspect the exact " *
            "durable reservation before binding its scheduler job.",
        ),
    )
    hasproperty(metadata, :publication_state) &&
        metadata.publication_state == :published || throw(ArgumentError(
        "metaMDBG scheduler job binding requires an inspected, fully " *
        "published submission reservation.",
    ))
    all(
        field -> hasproperty(metadata, field),
        (:reservation_identity, :shared_reservation_identity),
    ) || throw(ArgumentError(
        "metaMDBG scheduler job binding requires the exact private and " *
        "shared filesystem identities returned by inspection.",
    ))
    submission_state = hasproperty(metadata, :submission_state) ?
                       metadata.submission_state : :reserved
    has_pending_evidence = submission_state in (
        :submission_pending,
        :submission_commit_cleanup_pending,
    )
    if has_pending_evidence
        all(
            field -> hasproperty(metadata, field),
            (
                :pending_job_id,
                :pending_job_cluster,
                :pending_job_schema_version,
                :pending_job_path,
                :pending_job_identity,
                :pending_job_complete,
            ),
        ) || throw(ArgumentError(
            "metaMDBG pending-job recovery requires every exact pending " *
            "identity returned by inspection.",
        ))
        metadata.pending_job_id isa AbstractString || error(
            "metaMDBG pending-job recovery metadata has no exact scheduler " *
            "job ID.",
        )
    else
        metadata.job_id === nothing || error(
            "metaMDBG submission reservation already has a scheduler job id.",
        )
    end
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    normalized_job_cluster = _normalize_metamdbg_job_cluster(job_cluster)
    metadata_job_cluster = hasproperty(metadata, :job_cluster) ?
                           _normalize_metamdbg_job_cluster(
        metadata.job_cluster,
    ) : nothing
    if has_pending_evidence
        normalized_job_id == metadata.pending_job_id || error(
            "metaMDBG requested scheduler job ID does not match the exact " *
            "pending job record.",
        )
        normalized_job_cluster == metadata.pending_job_cluster || error(
            "metaMDBG requested scheduler cluster does not match the exact " *
            "pending job record.",
        )
        if metadata.job_id isa AbstractString
            normalized_job_id == metadata.job_id || error(
                "metaMDBG committed and pending scheduler job IDs conflict.",
            )
            normalized_job_cluster == metadata_job_cluster || error(
                "metaMDBG committed and pending scheduler clusters conflict.",
            )
        end
    end
    normalized_owner_token = String(owner_token)
    normalized_owner_token == metadata.owner_token || error(
        "metaMDBG reservation owner token does not match the durable " *
        "reservation capability.",
    )
    graph_k = metadata.graph_k
    graph_k isa Int || throw(ArgumentError(
        "metaMDBG reservation graph_k metadata must be an Int.",
    ))
    outputs = _metamdbg_output_paths(String(metadata.canonical_outdir), graph_k)
    requested_job_name = _metamdbg_scheduler_job_name_prefix(
        String(metadata.scheduler_job_name),
        String(metadata.workflow_signature),
    )
    expected = _metamdbg_submission_reservation(
        outputs,
        (; signature = String(metadata.input_contract_signature)),
        graph_k;
        owner_token = normalized_owner_token,
        job_name = requested_job_name,
    )
    expected.path == metadata.path || error(
        "metaMDBG reservation path does not match its recomputed workflow path.",
    )
    expected.workflow_signature == metadata.workflow_signature || error(
        "metaMDBG reservation workflow signature does not match metadata.",
    )
    bound = _with_metamdbg_output_domain_lock(
        outputs.outdir;
        allowed_same_root_locks = (
            expected.output_root_reservation_marker,
        ),
    ) do
        current = _metamdbg_submission_reservation_from_path(
            expected.path,
            outputs.outdir,
        )
        for field in (
                :canonical_outdir,
                :path,
                :workflow_signature,
                :scheduler_job_name,
                :input_contract_signature,
                :graph_k,
                :owner_token,
        )
            getproperty(current, field) == getproperty(expected, field) || error(
                "metaMDBG submission reservation $(field) changed before its " *
                "scheduler job could be bound.",
            )
        end
        if current.job_id != metadata.job_id ||
           current.job_cluster != metadata_job_cluster
            metadata.job_id === nothing && current.job_id isa AbstractString &&
                error(
                    "metaMDBG submission reservation already has a scheduler " *
                    "job id.",
                )
            error(
                "metaMDBG submission scheduler job state changed after " *
                "inspection.",
            )
        end
        _require_unchanged_metamdbg_recovery_identities(
            current,
            metadata.reservation_identity,
            metadata.shared_reservation_identity,
        )
        if has_pending_evidence
            pending_paths = _metamdbg_pending_submission_job_paths(
                outputs.outdir,
            )
            pending =
                _metamdbg_pending_submission_job_record_for_reservation(
                    current,
                    pending_paths;
                    allow_incomplete = true,
                )
            pending isa NamedTuple || error(
                "metaMDBG inspected pending scheduler job record disappeared " *
                "before exact recovery.",
            )
            pending.path == metadata.pending_job_path || error(
                "metaMDBG pending scheduler job path changed after inspection.",
            )
            pending.job_id == normalized_job_id || error(
                "metaMDBG pending scheduler job ID changed after inspection.",
            )
            pending.job_cluster == normalized_job_cluster || error(
                "metaMDBG pending scheduler cluster changed after inspection.",
            )
            pending.job_schema_version ==
                metadata.pending_job_schema_version || error(
                "metaMDBG pending scheduler job schema changed after " *
                "inspection.",
            )
            pending.identity == metadata.pending_job_identity || error(
                "metaMDBG pending scheduler job inode changed after inspection.",
            )
            pending.is_complete == metadata.pending_job_complete || error(
                "metaMDBG pending scheduler job completion state changed after " *
                "inspection.",
            )
            pending_paths == [pending.path] || error(
                "metaMDBG found additional pending scheduler job evidence after " *
                "inspection; refusing ambiguous exact-job recovery.",
            )
            _resume_metamdbg_pending_submission_job_record!(
                current,
                pending,
                metadata.reservation_identity,
                metadata.shared_reservation_identity,
            )
        else
            pending_paths = _metamdbg_pending_submission_job_paths(
                outputs.outdir,
            )
            isempty(pending_paths) || error(
                "metaMDBG pending scheduler job evidence appeared after " *
                "inspection; reinspect and recover the exact accepted job ID.",
            )
            current.job_id === nothing || error(
                "metaMDBG submission reservation already has a scheduler job " *
                "id.",
            )
            _bind_metamdbg_submission_job!(
                current,
                normalized_job_id,
                normalized_job_cluster,
            )
        end
    end
    return (;
        canonical_outdir = bound.canonical_outdir,
        path = bound.path,
        workflow_signature = bound.workflow_signature,
        scheduler_job_name = bound.scheduler_job_name,
        input_contract_signature = bound.input_contract_signature,
        graph_k = bound.graph_k,
        owner_token = bound.owner_token,
        job_id = bound.job_id,
        job_cluster = bound.job_cluster,
        job_schema_version = bound.job_schema_version,
        submission_state = bound.submission_state,
        publication_state = :published,
        pending_job_id = nothing,
        pending_job_cluster = nothing,
        pending_job_schema_version = nothing,
        pending_job_path = nothing,
        pending_job_identity = nothing,
        pending_job_complete = nothing,
        reservation_identity = metadata.reservation_identity,
        shared_reservation_identity = metadata.shared_reservation_identity,
    )
end

function _remove_metamdbg_submission_reservation!(
        reservation::NamedTuple,
        ;
        post_private_rename_hook::Function =
            (_reservation::NamedTuple) -> nothing,
        post_shared_release_hook::Function =
            (_reservation::NamedTuple) -> nothing,
)::Nothing
    output_root_marker = reservation.output_root_reservation_marker
    pending_job_records = _metamdbg_pending_submission_job_paths(
        reservation.canonical_outdir,
    )
    isempty(pending_job_records) || error(
        "metaMDBG refuses reservation removal while pending scheduler job " *
        "evidence exists: $(join(pending_job_records, ", ")). Recover and " *
        "verify the exact accepted job first.",
    )
    if !ispath(reservation.path) && !islink(reservation.path)
        if ispath(output_root_marker) || islink(output_root_marker)
            error(
                "metaMDBG submission reservation directory is missing while " *
                "its shared output-root reservation remains: " *
                "$(output_root_marker). Refusing capability-blind cleanup.",
            )
        end
        return nothing
    end
    path_state = _metamdbg_submission_reservation_path_state(
        reservation.path,
        reservation.canonical_outdir,
    )
    transition = if path_state in (:published, :provisional)
        _require_metamdbg_submission_reservation!(reservation)
        reclaiming_path = reservation.reclaiming_path
        if _output_root_path_entry_exists(reclaiming_path)
            error(
                "metaMDBG refuses to overwrite an existing identity-bound " *
                "reclaim transition: $(reclaiming_path).",
            )
        end
        mv(reservation.path, reclaiming_path)
        _fsync_metamdbg_directory(dirname(reclaiming_path))
        _metamdbg_submission_reservation_from_path(
            reclaiming_path,
            reservation.canonical_outdir,
        )
    elseif path_state == :reclaiming
        _metamdbg_submission_reservation_from_path(
            reservation.path,
            reservation.canonical_outdir,
        )
    else
        error(
            "metaMDBG cannot explicitly reclaim reservation state " *
            "$(repr(path_state)).",
        )
    end
    post_private_rename_hook(transition)
    if _output_root_path_entry_exists(output_root_marker)
        _require_metamdbg_output_root_reservation_marker!(transition)
        output_root_marker_identity = _metamdbg_regular_file_identity(
            output_root_marker,
        )
        _remove_exact_metamdbg_durable_file!(
            output_root_marker,
            output_root_marker_identity,
        )
    end
    !_output_root_path_entry_exists(output_root_marker) || error(
        "metaMDBG shared output-root reservation reappeared during explicit " *
        "reclaim: $(output_root_marker).",
    )
    post_shared_release_hook(transition)
    current = _metamdbg_submission_reservation_from_path(
        transition.path,
        transition.canonical_outdir,
    )
    current.reservation_state == :reclaim_release_pending || error(
        "metaMDBG explicit reclaim did not reach its durable shared-release " *
        "state.",
    )
    current_identity = _metamdbg_directory_path_identity(current.path)
    _remove_exact_metamdbg_durable_directory!(
        current.path,
        current_identity;
        recursive = true,
    )
    return nothing
end

function _require_metamdbg_optional_recovery_identity(
        path::AbstractString,
        expected_identity::Union{Nothing, NamedTuple},
        identity_function::Function,
        label::AbstractString,
)::Nothing
    exists = _output_root_path_entry_exists(path)
    (expected_identity === nothing) == !exists || error(
        "metaMDBG $(label) presence changed after recovery inspection: $(path).",
    )
    if exists
        identity_function(path) == expected_identity || error(
            "metaMDBG $(label) identity changed after recovery inspection: " *
            "$(path).",
        )
    end
    return nothing
end

function _remove_metamdbg_runtime_recovery_state!(
        reservation::NamedTuple,
        metadata::NamedTuple,
)::Nothing
    _metamdbg_submission_reservation_identity(reservation) ==
        metadata.reservation_identity || error(
        "metaMDBG runtime owner record changed after recovery inspection.",
    )
    _require_metamdbg_optional_recovery_identity(
        reservation.output_root_reservation_marker,
        metadata.queued_reservation_identity,
        path -> begin
            marker = _require_metamdbg_output_root_reservation_marker!(
                reservation,
            )
            marker == path || error("metaMDBG queued marker path changed.")
            status = stat(marker)
            return (; device = status.device, inode = status.inode)
        end,
        "queued shared reservation",
    )
    _require_metamdbg_optional_recovery_identity(
        reservation.runtime_output_root_reservation_marker,
        metadata.runtime_reservation_identity,
        path -> begin
            marker = _require_metamdbg_runtime_output_root_reservation_marker!(
                reservation,
            )
            marker == path || error("metaMDBG runtime marker path changed.")
            status = stat(marker)
            return (; device = status.device, inode = status.inode)
        end,
        "runtime shared reservation",
    )
    private_lock_path = _metamdbg_output_lock_path(reservation.canonical_outdir)
    _require_metamdbg_optional_recovery_identity(
        private_lock_path,
        metadata.private_lock_identity,
        _metamdbg_output_lock_identity,
        "runtime private lifecycle lock",
    )
    cleanup_path = _metamdbg_lifecycle_cleanup_reservation_path(
        reservation.canonical_outdir,
    )
    _require_metamdbg_optional_recovery_identity(
        cleanup_path,
        metadata.cleanup_reservation_identity,
        _metamdbg_output_lock_identity,
        "runtime lifecycle cleanup reservation",
    )
    if metadata.private_lock_identity !== nothing
        _recover_dead_metamdbg_private_output_lock!(private_lock_path)
    end
    if metadata.cleanup_reservation_identity !== nothing
        _remove_metamdbg_lifecycle_cleanup_reservation!(
            cleanup_path,
            metadata.cleanup_reservation_identity,
        )
    end
    if metadata.queued_reservation_identity !== nothing
        _remove_exact_metamdbg_durable_file!(
            reservation.output_root_reservation_marker,
            merge(
                metadata.queued_reservation_identity,
                (; path = reservation.output_root_reservation_marker),
            ),
        )
    end
    if metadata.runtime_reservation_identity !== nothing
        _remove_exact_metamdbg_durable_directory!(
            reservation.runtime_output_root_reservation_marker,
            merge(
                metadata.runtime_reservation_identity,
                (; path = reservation.runtime_output_root_reservation_marker),
            ),
        )
    end
    _metamdbg_submission_reservation_identity(reservation) ==
        metadata.reservation_identity || error(
        "metaMDBG runtime owner record changed during recovery.",
    )
    _remove_exact_metamdbg_durable_directory!(
        reservation.path,
        merge(metadata.reservation_identity, (; path = reservation.path));
        recursive = true,
    )
    return nothing
end

function _require_metamdbg_contract!(
        outputs::NamedTuple,
        expected_contract::NamedTuple,
)::String
    marker = outputs.contract_marker
    if !isfile(marker) || islink(marker)
        error(
            "metaMDBG provenance contract marker must be a regular, non-symlink " *
            "file: $(marker). Remove the stale output directory and rerun.",
        )
    end
    filesize(marker) > 0 || error(
        "metaMDBG provenance contract marker is empty: $(marker).",
    )
    actual_contents = read(marker, String)
    actual_contents == expected_contract.contents || error(
        "metaMDBG existing output contract does not match the normalized " *
        "input paths, input technology, input file metadata or SHA-256 content " *
        "digests, abundance_min, or exact toolchain for this invocation: " *
        "$(marker). Refusing stale output reuse.",
    )
    return marker
end

function _prepare_metamdbg_output_root!(
        outputs::NamedTuple,
        expected_contract::NamedTuple,
)::Bool
    canonical_outdir = _metamdbg_canonical_output_path(outputs.outdir)
    canonical_outdir == outputs.outdir || error(
        "metaMDBG output path is not canonical: $(outputs.outdir)",
    )
    if !ispath(canonical_outdir)
        mkpath(canonical_outdir)
    end
    isdir(canonical_outdir) && !islink(canonical_outdir) || error(
        "metaMDBG outdir must be a regular directory, not a symbolic link: " *
        "$(canonical_outdir)",
    )

    isempty(readdir(canonical_outdir)) && return false
    marker = outputs.contract_marker
    if !isfile(marker) || islink(marker)
        error(
            "metaMDBG refuses to adopt a nonempty output root without its " *
            "regular provenance contract marker: $(canonical_outdir). " *
            "Remove the partial or legacy output directory and rerun.",
        )
    end
    _require_metamdbg_contract!(outputs, expected_contract)
    return true
end

function _write_metamdbg_contract!(
        outputs::NamedTuple,
        input_contract::NamedTuple,
)::String
    marker = outputs.contract_marker
    if ispath(marker) || islink(marker)
        error(
            "metaMDBG refuses to overwrite an existing provenance contract " *
            "marker: $(marker).",
        )
    end
    temporary_path, temporary_io = mktemp(dirname(marker))
    temporary_identity = _metamdbg_regular_file_identity(temporary_path)
    published = false
    try
        write(temporary_io, input_contract.contents)
        chmod(temporary_path, 0o600)
        _fsync_metamdbg_file(temporary_io, temporary_path)
        close(temporary_io)
        filesize(temporary_path) > 0 || error(
            "Failed to write metaMDBG provenance contract marker: $(marker).",
        )
        Base.Filesystem.hardlink(temporary_path, marker)
        published = true
        _fsync_metamdbg_directory(dirname(marker))
        _remove_exact_metamdbg_durable_file!(
            temporary_path,
            temporary_identity,
        )
        read(marker, String) == input_contract.contents || error(
            "Failed to verify metaMDBG provenance contract marker: $(marker).",
        )
    catch primary_error
        if published && (ispath(marker) || islink(marker))
            try
                _remove_exact_metamdbg_durable_file!(
                    marker,
                    merge(temporary_identity, (; path = marker)),
                )
            catch cleanup_error
                @warn "metaMDBG failed to remove an invalid newly published " *
                      "contract marker while preserving the primary " *
                      "publication failure" marker primary_error cleanup_error
            end
        end
        Base.rethrow()
    finally
        isopen(temporary_io) && close(temporary_io)
        if _output_root_path_entry_exists(temporary_path)
            _remove_exact_metamdbg_durable_file!(
                temporary_path,
                temporary_identity,
            )
        end
    end
    return marker
end

function _metamdbg_nonempty_file(path::AbstractString)::Bool
    return isfile(path) && !islink(path) && filesize(path) > 0
end

function _require_metamdbg_regular_file(
        path::AbstractString,
        label::AbstractString,
)::String
    normalized_path = abspath(path)
    if !isfile(normalized_path) || islink(normalized_path)
        error("$(label) is not a regular, non-symlink file: $(normalized_path).")
    end
    filesize(normalized_path) > 0 || error(
        "$(label) is empty: $(normalized_path).",
    )
    return normalized_path
end

function _require_valid_metamdbg_fasta(
        path::AbstractString,
        label::AbstractString,
)::String
    normalized_path = _require_metamdbg_regular_file(path, label)
    reader = try
        Mycelia.open_fastx(normalized_path)
    catch caught
        caught isa InterruptException && rethrow()
        error(
            "$(label) is not valid FASTA: $(normalized_path). Cause: " *
            sprint(showerror, caught),
        )
    end
    record_count = 0
    identifiers = Set{String}()
    try
        for record in reader
            record isa FASTX.FASTA.Record || error(
                "$(label) is not valid FASTA: $(normalized_path).",
            )
            identifier = strip(String(FASTX.identifier(record)))
            isempty(identifier) && error(
                "$(label) contains an empty FASTA identifier at record " *
                "$(record_count + 1): $(normalized_path).",
            )
            identifier in identifiers && error(
                "$(label) contains duplicate FASTA identifier " *
                "$(repr(identifier)) at record $(record_count + 1): " *
                "$(normalized_path).",
            )
            sequence = FASTX.sequence(String, record)
            isempty(sequence) && error(
                "$(label) contains an empty FASTA sequence at record " *
                "$(record_count + 1): $(normalized_path).",
            )
            occursin(_METAMDBG_IUPAC_DNA_REGEX, sequence) || error(
                "$(label) contains invalid DNA at FASTA record " *
                "$(record_count + 1): $(normalized_path).",
            )
            try
                BioSequences.LongDNA{4}(sequence)
            catch caught
                caught isa InterruptException && rethrow()
                error(
                    "$(label) contains invalid DNA at FASTA record " *
                    "$(record_count + 1): $(normalized_path). Cause: " *
                    sprint(showerror, caught),
                )
            end
            push!(identifiers, identifier)
            record_count += 1
        end
    catch caught
        caught isa InterruptException && rethrow()
        if caught isa ErrorException && startswith(caught.msg, String(label))
            rethrow()
        end
        error(
            "$(label) is not valid FASTA: $(normalized_path). Cause: " *
            sprint(showerror, caught),
        )
    finally
        close(reader)
    end
    record_count > 0 || error(
        "$(label) contains no FASTA records: $(normalized_path).",
    )
    return normalized_path
end

function _require_valid_metamdbg_gfa_identifier(
        identifier::AbstractString,
        record_type::AbstractString,
        line_number::Int,
        label::AbstractString,
        path::AbstractString,
)::String
    normalized_identifier = String(identifier)
    is_valid_identifier =
        occursin(_METAMDBG_GFA_IDENTIFIER_REGEX, normalized_identifier) &&
        !occursin("+,", normalized_identifier) &&
        !occursin("-,", normalized_identifier)
    is_valid_identifier || error(
        "$(label) has an invalid $(record_type) identifier at line " *
        "$(line_number): $(path).",
    )
    return normalized_identifier
end

function _require_valid_metamdbg_gfa_orientation(
        orientation::AbstractString,
        record_type::AbstractString,
        line_number::Int,
        label::AbstractString,
        path::AbstractString,
)::String
    normalized_orientation = String(orientation)
    normalized_orientation in ("+", "-") || error(
        "$(label) has an invalid $(record_type) orientation at line " *
        "$(line_number): $(path).",
    )
    return normalized_orientation
end

function _require_valid_metamdbg_gfa_cigar(
        cigar::AbstractString,
        record_type::AbstractString,
        line_number::Int,
        label::AbstractString,
        path::AbstractString,
)::String
    normalized_cigar = String(cigar)
    occursin(_METAMDBG_GFA_CIGAR_REGEX, normalized_cigar) || error(
        "$(label) has an invalid $(record_type) CIGAR at line " *
        "$(line_number): $(path).",
    )
    return normalized_cigar
end

function _metamdbg_gfa_json_value_is_finite(value::Any)::Bool
    if value isa AbstractFloat
        return isfinite(value)
    elseif value isa AbstractVector
        return all(_metamdbg_gfa_json_value_is_finite, value)
    elseif value isa AbstractDict
        return all(_metamdbg_gfa_json_value_is_finite, values(value))
    end
    return true
end

function _is_valid_metamdbg_gfa_json(value::AbstractString)::Bool
    return try
        parsed = JSON.parse(String(value))
        _metamdbg_gfa_json_value_is_finite(parsed)
    catch caught
        caught isa InterruptException && rethrow()
        false
    end
end

function _require_valid_metamdbg_gfa_tag_value(
        tag_type::AbstractString,
        value::AbstractString,
        record_type::AbstractString,
        line_number::Int,
        label::AbstractString,
        path::AbstractString,
)::Nothing
    normalized_type = String(tag_type)
    normalized_value = String(value)
    is_valid = if normalized_type == "A"
        ncodeunits(normalized_value) == 1 &&
            occursin(r"^[!-~]$", normalized_value)
    elseif normalized_type == "i"
        occursin(_METAMDBG_GFA_TAG_INTEGER_REGEX, normalized_value) &&
            tryparse(Int64, normalized_value) !== nothing
    elseif normalized_type == "f"
        parsed = tryparse(Float64, normalized_value)
        occursin(_METAMDBG_GFA_TAG_FLOAT_REGEX, normalized_value) &&
            parsed !== nothing && isfinite(parsed)
    elseif normalized_type == "Z"
        occursin(_METAMDBG_GFA_TAG_PRINTABLE_REGEX, normalized_value)
    elseif normalized_type == "J"
        occursin(_METAMDBG_GFA_TAG_PRINTABLE_REGEX, normalized_value) &&
            _is_valid_metamdbg_gfa_json(normalized_value)
    elseif normalized_type == "H"
        iseven(ncodeunits(normalized_value)) &&
            occursin(_METAMDBG_GFA_TAG_HEX_REGEX, normalized_value)
    elseif normalized_type == "B"
        _is_valid_metamdbg_gfa_b_array(normalized_value)
    else
        false
    end
    is_valid || error(
        "$(label) has an invalid or unsupported $(record_type) optional tag " *
        "value at line $(line_number): $(path).",
    )
    return nothing
end

function _is_valid_metamdbg_gfa_b_array(value::AbstractString)::Bool
    components = split(String(value), ','; keepempty = true)
    length(components) >= 2 || return false
    subtype = first(components)
    ncodeunits(subtype) == 1 || return false
    values = @view components[2:end]
    if subtype == "f"
        return all(values) do item
            parsed = tryparse(Float32, item)
            return occursin(_METAMDBG_GFA_TAG_FLOAT_REGEX, item) &&
                   parsed !== nothing && isfinite(parsed)
        end
    end
    subtype_symbol = Symbol(subtype)
    hasproperty(_METAMDBG_GFA_B_INTEGER_RANGES, subtype_symbol) || return false
    lower, upper = getproperty(
        _METAMDBG_GFA_B_INTEGER_RANGES,
        subtype_symbol,
    )
    return all(values) do item
        occursin(_METAMDBG_GFA_TAG_INTEGER_REGEX, item) || return false
        parsed = tryparse(Int64, item)
        return parsed !== nothing && lower <= parsed <= upper
    end
end

function _require_valid_metamdbg_gfa_tags(
        fields::AbstractVector{<:AbstractString},
        first_tag_index::Int,
        record_type::AbstractString,
        line_number::Int,
        label::AbstractString,
        path::AbstractString,
)::Nothing
    tag_names = Set{String}()
    for tag_index in first_tag_index:length(fields)
        tag = fields[tag_index]
        components = split(tag, ':'; limit = 3, keepempty = true)
        length(components) == 3 || error(
            "$(label) has an invalid $(record_type) optional tag at line " *
            "$(line_number): $(path).",
        )
        tag_name, tag_type, value = String.(components)
        occursin(_METAMDBG_GFA_TAG_NAME_REGEX, tag_name) || error(
            "$(label) has an invalid $(record_type) optional tag name at line " *
            "$(line_number): $(path).",
        )
        tag_name in tag_names && error(
            "$(label) has duplicate $(record_type) optional tag name " *
            "$(repr(tag_name)) at line $(line_number): $(path).",
        )
        _require_valid_metamdbg_gfa_tag_value(
            tag_type,
            value,
            record_type,
            line_number,
            label,
            path,
        )
        push!(tag_names, tag_name)
    end
    return nothing
end

function _metamdbg_gfa_path_step_identifiers(
        field::AbstractString,
        line_number::Int,
        label::AbstractString,
        path::AbstractString,
)::Vector{String}
    normalized_field = String(field)
    isempty(normalized_field) && error(
        "$(label) has an empty GFA path step list at line $(line_number): " *
        "$(path).",
    )
    bytes = codeunits(normalized_field)
    byte_count = length(bytes)
    identifiers = String[]
    step_start = 1
    position = 1
    while position <= byte_count
        byte = bytes[position]
        is_orientation = byte == UInt8('+') || byte == UInt8('-')
        ends_step = is_orientation &&
                    (position == byte_count ||
                     bytes[position + 1] == UInt8(','))
        if !ends_step
            position += 1
            continue
        end
        position > step_start || error(
            "$(label) has a malformed GFA path step at line $(line_number): " *
            "$(path).",
        )
        identifier = String(Vector{UInt8}(bytes[step_start:(position - 1)]))
        push!(identifiers, _require_valid_metamdbg_gfa_identifier(
            identifier,
            "GFA path step",
            line_number,
            label,
            path,
        ))
        if position == byte_count
            step_start = byte_count + 1
            break
        end
        step_start = position + 2
        step_start <= byte_count || error(
            "$(label) has a trailing GFA path-step separator at line " *
            "$(line_number): $(path).",
        )
        position = step_start
    end
    step_start == byte_count + 1 || error(
        "$(label) has a malformed GFA path step at line $(line_number): " *
        "$(path).",
    )
    return identifiers
end

function _require_valid_metamdbg_gfa_input(
        input::IO,
        label::AbstractString,
        normalized_path::AbstractString,
)::Nothing
    seekstart(input)
    segment_identifiers = Set{String}()
    record_identifiers = Set{String}()
    for (line_number, raw_line) in enumerate(eachline(input))
        line = chomp(raw_line)
        isempty(strip(line)) && continue
        startswith(line, "#") && continue
        fields = split(line, '\t'; keepempty = true)
        record_type = String(first(fields))
        if record_type == "H"
            _require_valid_metamdbg_gfa_tags(
                fields,
                2,
                "GFA header",
                line_number,
                label,
                normalized_path,
            )
            for tag in @view fields[2:end]
                components = split(tag, ':'; limit = 3, keepempty = true)
                if first(components) == "VN"
                    components[2] == "Z" &&
                        occursin(r"^1\.[0-9]+$", components[3]) || error(
                        "$(label) declares a non-GFA1 VN header at line " *
                        "$(line_number): $(normalized_path).",
                    )
                end
            end
        elseif record_type == "S"
            length(fields) >= 3 || error(
                "$(label) has a malformed GFA segment at line " *
                "$(line_number): $(normalized_path).",
            )
            identifier = _require_valid_metamdbg_gfa_identifier(
                fields[2],
                "GFA segment",
                line_number,
                label,
                normalized_path,
            )
            identifier in record_identifiers && error(
                "$(label) has duplicate GFA segment/path name " *
                "$(repr(identifier)): $(normalized_path).",
            )
            sequence = String(fields[3])
            if isempty(sequence) || sequence == "*"
                error(
                    "$(label) has no sequence for GFA segment " *
                    "$(repr(identifier)): $(normalized_path).",
                )
            end
            occursin(_METAMDBG_IUPAC_DNA_REGEX, sequence) || error(
                "$(label) has invalid DNA for GFA segment " *
                "$(repr(identifier)): $(normalized_path).",
            )
            try
                BioSequences.LongDNA{4}(sequence)
            catch caught
                caught isa InterruptException && rethrow()
                error(
                    "$(label) has invalid DNA for GFA segment " *
                    "$(repr(identifier)): $(normalized_path). Cause: " *
                    sprint(showerror, caught),
                )
            end
            _require_valid_metamdbg_gfa_tags(
                fields,
                4,
                "GFA segment",
                line_number,
                label,
                normalized_path,
            )
            push!(segment_identifiers, identifier)
            push!(record_identifiers, identifier)
        elseif record_type == "L"
            length(fields) >= 6 || error(
                "$(label) has a malformed GFA link at line " *
                "$(line_number): $(normalized_path).",
            )
            _require_valid_metamdbg_gfa_identifier(
                fields[2],
                "GFA link source",
                line_number,
                label,
                normalized_path,
            )
            _require_valid_metamdbg_gfa_orientation(
                fields[3],
                "GFA link source",
                line_number,
                label,
                normalized_path,
            )
            _require_valid_metamdbg_gfa_identifier(
                fields[4],
                "GFA link destination",
                line_number,
                label,
                normalized_path,
            )
            _require_valid_metamdbg_gfa_orientation(
                fields[5],
                "GFA link destination",
                line_number,
                label,
                normalized_path,
            )
            _require_valid_metamdbg_gfa_cigar(
                fields[6],
                "GFA link",
                line_number,
                label,
                normalized_path,
            )
            _require_valid_metamdbg_gfa_tags(
                fields,
                7,
                "GFA link",
                line_number,
                label,
                normalized_path,
            )
        elseif record_type == "P"
            length(fields) >= 4 || error(
                "$(label) has a malformed GFA path at line " *
                "$(line_number): $(normalized_path).",
            )
            path_identifier = _require_valid_metamdbg_gfa_identifier(
                fields[2],
                "GFA path",
                line_number,
                label,
                normalized_path,
            )
            path_identifier in record_identifiers && error(
                "$(label) has duplicate GFA segment/path name " *
                "$(repr(path_identifier)): $(normalized_path).",
            )
            steps = _metamdbg_gfa_path_step_identifiers(
                fields[3],
                line_number,
                label,
                normalized_path,
            )
            overlap_field = String(fields[4])
            if overlap_field != "*"
                overlaps = split(overlap_field, ','; keepempty = true)
                length(overlaps) == length(steps) - 1 || error(
                    "$(label) has the wrong number of GFA path overlaps " *
                    "at line $(line_number): $(normalized_path).",
                )
                for overlap in overlaps
                    _require_valid_metamdbg_gfa_cigar(
                        overlap,
                        "GFA path",
                        line_number,
                        label,
                        normalized_path,
                    )
                end
            end
            _require_valid_metamdbg_gfa_tags(
                fields,
                5,
                "GFA path",
                line_number,
                label,
                normalized_path,
            )
            push!(record_identifiers, path_identifier)
        else
            error(
                "$(label) has unknown GFA record type " *
                "$(repr(record_type)) at line $(line_number): " *
                "$(normalized_path).",
            )
        end
    end
    isempty(segment_identifiers) && error(
        "$(label) contains no sequence-bearing GFA segments: " *
        "$(normalized_path).",
    )
    seekstart(input)
    for (line_number, raw_line) in enumerate(eachline(input))
        line = chomp(raw_line)
        isempty(strip(line)) && continue
        startswith(line, "#") && continue
        fields = split(line, '\t'; keepempty = true)
        record_type = String(first(fields))
        referenced_identifiers = if record_type == "L"
            String[fields[2], fields[4]]
        elseif record_type == "P"
            _metamdbg_gfa_path_step_identifiers(
                fields[3],
                line_number,
                label,
                normalized_path,
            )
        else
            continue
        end
        reference_type = record_type == "L" ? "link" : "path"
        for identifier in referenced_identifiers
            identifier in segment_identifiers || error(
                "$(label) has a dangling GFA $(reference_type) segment " *
                "reference $(repr(identifier)) at line " *
                "$(line_number): $(normalized_path).",
            )
        end
    end
    return nothing
end

function _require_valid_metamdbg_gfa(
        path::AbstractString,
        label::AbstractString,
)::String
    normalized_path = _require_metamdbg_regular_file(path, label)
    open(normalized_path, "r") do input
        _require_valid_metamdbg_gfa_input(
            input,
            label,
            normalized_path,
        )
    end
    return normalized_path
end

function _existing_metamdbg_contigs(
        outputs::NamedTuple,
)::Union{Nothing, String}
    if ispath(outputs.contigs_gz) || islink(outputs.contigs_gz)
        _require_valid_metamdbg_fasta(
            outputs.contigs_gz,
            "metaMDBG compressed contigs artifact",
        )
        return outputs.contigs_gz
    elseif ispath(outputs.contigs_plain) || islink(outputs.contigs_plain)
        _require_valid_metamdbg_fasta(
            outputs.contigs_plain,
            "metaMDBG plain contigs artifact",
        )
        return outputs.contigs_plain
    end
    return nothing
end

function _gzip_metamdbg_contigs!(
        plain_path::AbstractString,
        gzip_path::AbstractString,
)::String
    temporary_path, temporary_io = mktemp(dirname(gzip_path))
    close(temporary_io)
    try
        open(plain_path, "r") do input
            open(temporary_path, "w") do raw_output
                gzip_output = CodecZlib.GzipCompressorStream(raw_output)
                try
                    buffer = Vector{UInt8}(undef, 1024 * 1024)
                    while !eof(input)
                        bytes_read = readbytes!(input, buffer)
                        bytes_read == 0 && break
                        write(gzip_output, view(buffer, 1:bytes_read))
                    end
                finally
                    close(gzip_output)
                end
            end
        end
        _metamdbg_nonempty_file(temporary_path) || error(
            "Failed to gzip metaMDBG contigs from $(plain_path).",
        )
        mv(temporary_path, gzip_path; force = true)
    finally
        rm(temporary_path; force = true)
    end
    return String(gzip_path)
end

function _normalize_metamdbg_contigs!(
        outputs::NamedTuple,
)::Union{Nothing, String}
    existing_contigs = _existing_metamdbg_contigs(outputs)
    existing_contigs === nothing && return nothing
    existing_contigs == outputs.contigs_gz && return outputs.contigs_gz
    compressed_contigs = _gzip_metamdbg_contigs!(
        outputs.contigs_plain,
        outputs.contigs_gz,
    )
    return _require_valid_metamdbg_fasta(
        compressed_contigs,
        "metaMDBG compressed contigs artifact",
    )
end

function _metamdbg_graph_candidates(
        outputs::NamedTuple,
        graph_k::Int,
)::Vector{String}
    _require_positive_metamdbg_graph_k(graph_k)
    isdir(outputs.outdir) || return String[]
    prefix = "assemblyGraph_k$(graph_k)_"
    return sort!(filter(
        path -> begin
            filename = basename(path)
            startswith(filename, prefix) && endswith(filename, "bps.gfa")
        end,
        readdir(outputs.outdir; join = true),
    ))
end

function _metamdbg_all_graph_candidates(
        outputs::NamedTuple,
)::Vector{String}
    isdir(outputs.outdir) || return String[]
    return sort!(filter(
        path -> occursin(
            r"^assemblyGraph_k[0-9]+_.*bps\.gfa$",
            basename(path),
        ),
        readdir(outputs.outdir; join = true),
    ))
end

function _normalize_metamdbg_graph!(
        outputs::NamedTuple,
        graph_k::Int,
)::Union{Nothing, String}
    _require_positive_metamdbg_graph_k(graph_k)
    candidates = _metamdbg_graph_candidates(outputs, graph_k)
    isempty(candidates) && return nothing
    length(candidates) == 1 || error(
        "metaMDBG produced multiple graph artifacts for k=$(graph_k): " *
        join(candidates, ", "),
    )
    all_candidates = _metamdbg_all_graph_candidates(outputs)
    length(all_candidates) == 1 || error(
        "metaMDBG requires exactly one total graph artifact per output root; " *
        "found $(length(all_candidates)): $(join(all_candidates, ", ")).",
    )
    only(all_candidates) == only(candidates) || error(
        "metaMDBG graph inventory does not match requested k=$(graph_k).",
    )
    graph_source = _require_valid_metamdbg_gfa(
        only(candidates),
        "metaMDBG graph artifact",
    )

    if ispath(outputs.graph_alias) || islink(outputs.graph_alias)
        rm(outputs.graph_alias; force = true)
    end
    symlink(basename(graph_source), outputs.graph_alias)
    if !islink(outputs.graph_alias) || !isfile(outputs.graph_alias) ||
       filesize(outputs.graph_alias) == 0
        error(
            "metaMDBG graph alias does not resolve to a nonempty artifact: " *
            outputs.graph_alias,
        )
    end
    return outputs.graph_alias
end

function _require_metamdbg_artifacts!(
        outputs::NamedTuple,
        graph_k::Int,
        ;
        output_root_identity::Union{Nothing, NamedTuple} = nothing,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    if output_root_identity !== nothing
        _require_unchanged_metamdbg_output_root(
            output_root_identity,
            outputs.outdir,
        )
    end
    contigs = _normalize_metamdbg_contigs!(outputs)
    contigs === nothing && error(
        "metaMDBG produced neither contigs.fasta.gz nor contigs.fasta in " *
        outputs.outdir,
    )
    graph = _normalize_metamdbg_graph!(outputs, graph_k)
    graph === nothing && error(
        "metaMDBG produced no nonempty assemblyGraph_k$(graph_k)_*bps.gfa " *
        "artifact in $(outputs.outdir).",
    )
    if output_root_identity !== nothing
        _require_unchanged_metamdbg_output_root(
            output_root_identity,
            outputs.outdir,
        )
        _require_metamdbg_artifact_containment(
            contigs,
            "metaMDBG contigs artifact",
            output_root_identity,
        )
        _require_metamdbg_artifact_containment(
            graph,
            "metaMDBG graph artifact",
            output_root_identity,
        )
    end
    return (; outdir = outputs.outdir, contigs, graph)
end

function _metamdbg_completion_toolchain_summary(
        toolchain::NamedTuple,
)::NamedTuple
    expected = _metamdbg_expected_toolchain()
    for field in propertynames(expected)
        if !hasproperty(toolchain, field) ||
           getproperty(toolchain, field) != getproperty(expected, field)
            error(
                "metaMDBG completion toolchain has incompatible $(field) " *
                "provenance.",
            )
        end
    end
    hasproperty(toolchain, :package_inventory_sha256) || error(
        "metaMDBG completion toolchain is missing its normalized package " *
        "inventory digest.",
    )
    inventory_sha256 = toolchain.package_inventory_sha256
    if !(inventory_sha256 isa AbstractString) ||
       !occursin(r"^[0-9a-f]{64}$", inventory_sha256)
        error(
            "metaMDBG completion toolchain package inventory digest is " *
            "invalid.",
        )
    end
    hasproperty(toolchain, :package_count) || error(
        "metaMDBG completion toolchain is missing its package count.",
    )
    package_count = toolchain.package_count
    if !(package_count isa Integer) || package_count <= 0
        error("metaMDBG completion toolchain package count must be positive.")
    end
    return (;
        expected...,
        package_inventory_sha256 = String(inventory_sha256),
        package_count = Int(package_count),
    )
end

function _metamdbg_completion_manifest(
        outputs::NamedTuple,
        artifacts::NamedTuple,
        input_contract::NamedTuple,
        graph_k::Int,
        toolchain::NamedTuple,
        ;
        digest_function::Function = _metamdbg_sha256,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    toolchain_summary = _metamdbg_completion_toolchain_summary(toolchain)
    workflow = _metamdbg_workflow_contract(outputs, input_contract, graph_k)
    contigs_path = realpath(artifacts.contigs)
    graph_alias_path = normpath(abspath(artifacts.graph))
    graph_source_path = realpath(artifacts.graph)
    manifest = (;
        schema_version = _METAMDBG_COMPLETION_SCHEMA_VERSION,
        workflow = (;
            canonical_outdir_sha256 =
                _metamdbg_string_sha256(outputs.outdir),
            input_contract_signature = input_contract.signature,
            graph_k,
            workflow_signature = workflow.signature,
        ),
        toolchain = toolchain_summary,
        artifacts = (;
            contigs = (;
                canonical_path_sha256 =
                    _metamdbg_string_sha256(contigs_path),
                size_bytes = Int64(filesize(contigs_path)),
                sha256 = String(digest_function(contigs_path)),
            ),
            graph = (;
                alias_path_sha256 =
                    _metamdbg_string_sha256(graph_alias_path),
                source_path_sha256 =
                    _metamdbg_string_sha256(graph_source_path),
                size_bytes = Int64(filesize(graph_source_path)),
                sha256 = String(digest_function(graph_source_path)),
            ),
        ),
    )
    serialized_manifest = JSON.json(manifest)
    signature = _metamdbg_string_sha256(serialized_manifest)
    contents = JSON.json((;
        schema_version = _METAMDBG_COMPLETION_SCHEMA_VERSION,
        signature_algorithm = "sha256",
        signature,
        manifest,
    )) * "\n"
    return (; manifest, signature, contents)
end

function _write_metamdbg_completion_manifest!(
        outputs::NamedTuple,
        completion::NamedTuple,
)::String
    marker = outputs.completion_marker
    if ispath(marker) || islink(marker)
        error(
            "metaMDBG refuses to overwrite an existing completion manifest: " *
            "$(marker).",
        )
    end
    temporary_path, temporary_io = mktemp(dirname(marker))
    temporary_identity = _metamdbg_regular_file_identity(temporary_path)
    published = false
    try
        write(temporary_io, completion.contents)
        chmod(temporary_path, 0o600)
        _fsync_metamdbg_file(temporary_io, temporary_path)
        close(temporary_io)
        Base.Filesystem.hardlink(temporary_path, marker)
        published = true
        _fsync_metamdbg_directory(dirname(marker))
        _remove_exact_metamdbg_durable_file!(
            temporary_path,
            temporary_identity,
        )
        read(marker, String) == completion.contents || error(
            "Failed to write metaMDBG completion manifest: $(marker).",
        )
    catch primary_error
        isopen(temporary_io) && close(temporary_io)
        if _output_root_path_entry_exists(temporary_path)
            try
                _remove_exact_metamdbg_durable_file!(
                    temporary_path,
                    temporary_identity,
                )
            catch cleanup_error
                @warn "metaMDBG failed to clean an unpublished completion " *
                      "manifest while preserving the primary publication " *
                      "failure" temporary_path primary_error cleanup_error
            end
        end
        if published && _output_root_path_entry_exists(marker)
            try
                _remove_exact_metamdbg_durable_file!(
                    marker,
                    merge(temporary_identity, (; path = marker)),
                )
            catch cleanup_error
                @warn "metaMDBG failed to remove an invalid newly published " *
                      "completion manifest while preserving the primary " *
                      "publication failure" marker primary_error cleanup_error
            end
        end
        Base.rethrow()
    end
    return marker
end

function _metamdbg_completion_toolchain_from_json(
        manifest::AbstractDict,
)::NamedTuple
    toolchain = get(manifest, "toolchain", nothing)
    toolchain isa AbstractDict || error(
        "metaMDBG completion manifest toolchain must be a JSON object.",
    )
    expected_keys = Set(String[
        "metamdbg_version",
        "environment_name",
        "environment_spec_sha256",
        "package_inventory_sha256",
        "package_count",
    ])
    Set(String.(keys(toolchain))) == expected_keys || error(
        "metaMDBG completion manifest toolchain has unexpected fields.",
    )
    return _metamdbg_completion_toolchain_summary((;
        metamdbg_version = get(toolchain, "metamdbg_version", nothing),
        environment_name = get(toolchain, "environment_name", nothing),
        environment_spec_sha256 =
            get(toolchain, "environment_spec_sha256", nothing),
        package_inventory_sha256 =
            get(toolchain, "package_inventory_sha256", nothing),
        package_count = get(toolchain, "package_count", nothing),
    ))
end

function _require_metamdbg_completion_manifest!(
        outputs::NamedTuple,
        artifacts::NamedTuple,
        input_contract::NamedTuple,
        graph_k::Int,
        ;
        digest_function::Function = _metamdbg_sha256,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    marker = outputs.completion_marker
    if !isfile(marker) || islink(marker) || filesize(marker) == 0
        error(
            "metaMDBG complete output is missing its regular, nonempty " *
            "completion manifest: $(marker). Refusing unbound artifact reuse.",
        )
    end
    actual_contents = read(marker, String)
    parsed = try
        JSON.parse(actual_contents)
    catch caught
        error(
            "metaMDBG completion manifest is not valid JSON: $(marker). " *
            "Cause: $(sprint(showerror, caught))",
        )
    end
    parsed isa AbstractDict || error(
        "metaMDBG completion manifest must be a JSON object: $(marker).",
    )
    Set(String.(keys(parsed))) == Set(String[
        "schema_version",
        "signature_algorithm",
        "signature",
        "manifest",
    ]) || error(
        "metaMDBG completion manifest has unexpected top-level fields: " *
        "$(marker).",
    )
    get(parsed, "schema_version", nothing) ==
        _METAMDBG_COMPLETION_SCHEMA_VERSION || error(
        "metaMDBG completion manifest has an unsupported schema version.",
    )
    get(parsed, "signature_algorithm", nothing) == "sha256" || error(
        "metaMDBG completion manifest has an unsupported signature algorithm.",
    )
    manifest = get(parsed, "manifest", nothing)
    manifest isa AbstractDict || error(
        "metaMDBG completion manifest payload must be a JSON object.",
    )
    toolchain = _metamdbg_completion_toolchain_from_json(manifest)
    expected = _metamdbg_completion_manifest(
        outputs,
        artifacts,
        input_contract,
        graph_k,
        toolchain;
        digest_function,
    )
    actual_contents == expected.contents || error(
        "metaMDBG completion manifest does not match graph_k, workflow, exact " *
        "resolved toolchain inventory, or current artifact identities, sizes, " *
        "and SHA-256 digests: $(marker). Refusing stale output reuse.",
    )
    return expected
end

function _metamdbg_shell_literal(value::AbstractString)::String
    occursin('\0', value) && throw(ArgumentError(
        "metaMDBG cannot embed a NUL byte in generated shell source.",
    ))
    return "'$(replace(String(value), "'" => "'\"'\"'"))'"
end

function _metamdbg_executor_script(
        asm_cmd::String,
        gfa_cmd::String,
        outputs::NamedTuple,
        graph_k::Int,
        input_contract::NamedTuple,
        ;
        conda_runner::AbstractString = _conda_runner(),
        environment_path::AbstractString = last(_metamdbg_paths()),
        submission_reservation::Union{Nothing, NamedTuple} = nothing,
        threads::Union{Nothing, Int} = nothing,
        post_runtime_marker_publication_hook::Union{Nothing, String} = nothing,
        post_runtime_private_rename_hook::Union{Nothing, String} = nothing,
        post_queued_marker_release_hook::Union{Nothing, String} = nothing,
        post_runtime_marker_release_hook::Union{Nothing, String} = nothing,
        post_consumed_owner_rename_hook::Union{Nothing, String} = nothing,
        post_completion_publication_hook::Union{Nothing, String} = nothing,
        directory_enumerator_post_open_hook::Union{Nothing, String} = nothing,
        directory_enumerator_post_enumeration_hook::Union{Nothing, String} =
            nothing,
        directory_enumerator_post_inventory_hook::Union{Nothing, String} =
            nothing,
        directory_enumerator_post_directory_open_hook::Union{Nothing, String} =
            nothing,
        directory_enumerator_post_directory_enumeration_hook::Union{
            Nothing,
            String,
        } = nothing,
        lock_retry_attempts::Int = _METAMDBG_OUTPUT_LOCK_RETRY_ATTEMPTS,
        lock_retry_delay_seconds::Real =
            _METAMDBG_OUTPUT_LOCK_RETRY_DELAY_SECONDS,
)::String
    _require_positive_metamdbg_graph_k(graph_k)
    lock_retry_attempts > 0 || throw(ArgumentError(
        "metaMDBG executor lock_retry_attempts must be positive.",
    ))
    isfinite(lock_retry_delay_seconds) && lock_retry_delay_seconds >= 0 ||
        throw(ArgumentError(
            "metaMDBG executor lock_retry_delay_seconds must be finite and " *
            "nonnegative.",
        ))
    threads === nothing || threads > 0 || throw(ArgumentError(
        "metaMDBG executor threads must be positive when provided.",
    ))
    verified_environment_path =
        _require_verified_metamdbg_environment_spec(environment_path)
    outdir_parent = dirname(outputs.outdir)
    outdir_identity = _output_root_reservation_identity(outputs.outdir)
    pending_submission_job_prefix = joinpath(
        outdir_parent,
        _metamdbg_pending_submission_job_prefix(outputs.outdir),
    )
    lock_path = _metamdbg_output_lock_path(outputs.outdir)
    effective_submission_reservation = if submission_reservation === nothing
        _metamdbg_submission_reservation(outputs, input_contract, graph_k)
    else
        submission_reservation
    end
    pending_submission_job_owner_prefix =
        pending_submission_job_prefix *
        effective_submission_reservation.output_root_reservation_capability *
        "."
    directory_enumerator = _metamdbg_shell_literal(
        _metamdbg_runtime_directory_enumerator_python(
            ;
            post_open_hook = directory_enumerator_post_open_hook,
            post_enumeration_hook =
                directory_enumerator_post_enumeration_hook,
            post_inventory_hook = directory_enumerator_post_inventory_hook,
            post_directory_open_hook =
                directory_enumerator_post_directory_open_hook,
            post_directory_enumeration_hook =
                directory_enumerator_post_directory_enumeration_hook,
        ),
    )
    input_declaration_lines = String[]
    input_validation_function_lines = String[
        "validate_metamdbg_inputs() {",
        "  require_unchanged_metamdbg_output_root",
    ]
    for (input_index, input) in enumerate(input_contract.contract.inputs)
        input_path_variable = "input_path_$(input_index)"
        expected_digest_variable = "expected_input_sha256_$(input_index)"
        actual_digest_variable = "actual_input_sha256_$(input_index)"
        append!(input_declaration_lines, String[
            "$(input_path_variable)=$(_metamdbg_shell_literal(input.path))",
            "$(expected_digest_variable)=$(_metamdbg_shell_literal(input.sha256))",
        ])
        append!(input_validation_function_lines, String[
            "if [ ! -f \"\$$(input_path_variable)\" ]; then",
            "  echo \"metaMDBG input is missing during execution: \$$(input_path_variable)\" >&2",
            "  return 1",
            "fi",
            "$(actual_digest_variable)=\$(sha256_file \"\$$(input_path_variable)\")",
            "if [ \"\$$(actual_digest_variable)\" != \"\$$(expected_digest_variable)\" ]; then",
            "  echo \"metaMDBG input content changed during execution: \$$(input_path_variable)\" >&2",
            "  return 1",
            "fi",
        ])
    end
    append!(input_validation_function_lines, String[
        "  require_unchanged_metamdbg_output_root",
        "}",
    ])
    staged_input_declaration_lines = String[]
    staged_input_creation_lines = String[
        "staged_inputs_dir=\"\$secure_tmpdir/inputs\"",
        "mkdir -m 700 -- \"\$staged_inputs_dir\"",
    ]
    staged_input_validation_lines = String[
        "validate_staged_metamdbg_inputs() {",
        "  require_unchanged_metamdbg_output_root",
    ]
    staged_input_references = String[]
    for (input_index, input) in enumerate(input_contract.contract.inputs)
        input_path_variable = "input_path_$(input_index)"
        expected_digest_variable = "expected_input_sha256_$(input_index)"
        staged_path_variable = "staged_input_path_$(input_index)"
        staged_digest_variable = "staged_input_sha256_$(input_index)"
        staged_name = _metamdbg_staged_input_name(input_index, input.path)
        push!(staged_input_declaration_lines,
            "$(staged_path_variable)=\"\$staged_inputs_dir/\"$(_metamdbg_shell_literal(staged_name))",
        )
        append!(staged_input_creation_lines, String[
            "cp -- \"\$$(input_path_variable)\" \"\$$(staged_path_variable)\"",
            "chmod 400 \"\$$(staged_path_variable)\"",
            "$(staged_digest_variable)=\$(sha256_file \"\$$(staged_path_variable)\")",
            "if [ \"\$$(staged_digest_variable)\" != \"\$$(expected_digest_variable)\" ]; then",
            "  echo \"metaMDBG staged input does not match its SHA-256 contract\" >&2",
            "  exit 1",
            "fi",
        ])
        append!(staged_input_validation_lines, String[
            "if [ ! -f \"\$$(staged_path_variable)\" ] || [ -L \"\$$(staged_path_variable)\" ]; then",
            "  echo \"metaMDBG staged input is missing or not regular\" >&2",
            "  return 1",
            "fi",
            "$(staged_digest_variable)=\$(sha256_file \"\$$(staged_path_variable)\")",
            "if [ \"\$$(staged_digest_variable)\" != \"\$$(expected_digest_variable)\" ]; then",
            "  echo \"metaMDBG staged input changed during execution\" >&2",
            "  return 1",
            "fi",
        ])
        push!(staged_input_references, "\"\$$(staged_path_variable)\"")
    end
    append!(staged_input_validation_lines, String[
        "  require_unchanged_metamdbg_output_root",
        "}",
    ])
    runtime_asm_cmd = if threads === nothing
        asm_cmd
    else
        join(String[
            "\"\$conda_runner\" run --live-stream -n \"\$environment_name\"",
            "metaMDBG asm --out-dir \"\$outdir\"",
            input_contract.contract.input_flag,
            staged_input_references...,
            "--min-abundance $(input_contract.contract.abundance_min)",
            "--threads $(threads)",
        ], " ")
    end
    runtime_gfa_cmd = if threads === nothing
        gfa_cmd
    else
        join(String[
            "\"\$conda_runner\" run --live-stream -n \"\$environment_name\"",
            "metaMDBG gfa --assembly-dir \"\$outdir\"",
            "--k $(graph_k)",
            "--threads $(threads)",
        ], " ")
    end
    job_record_parts = _metamdbg_submission_job_parts(
        effective_submission_reservation,
    )
    job_record_cluster_prefix = "\",\"job_cluster\":\""
    job_record_cluster_suffix = "\"}\n"
    reservation_declaration_lines = String[
        "submission_reservation_dir=$(_metamdbg_shell_literal(effective_submission_reservation.path))",
        "runtime_submission_reservation_dir=$(_metamdbg_shell_literal(effective_submission_reservation.runtime_path))",
        "consumed_submission_reservation_dir=$(_metamdbg_shell_literal(effective_submission_reservation.consumed_path))",
        "submission_reservation_contract=$(_metamdbg_shell_literal(effective_submission_reservation.contract_marker))",
        "submission_reservation_job=$(_metamdbg_shell_literal(effective_submission_reservation.job_marker))",
        "submission_output_root_reservation=$(_metamdbg_shell_literal(effective_submission_reservation.output_root_reservation_marker))",
        "runtime_output_root_reservation=$(_metamdbg_shell_literal(effective_submission_reservation.runtime_output_root_reservation_marker))",
        "runtime_cleanup_failure_reservation=$(_metamdbg_shell_literal(_metamdbg_lifecycle_cleanup_reservation_path(outputs.outdir)))",
        "pending_submission_job_prefix=$(_metamdbg_shell_literal(pending_submission_job_prefix))",
        "pending_submission_job_owner_prefix=$(_metamdbg_shell_literal(pending_submission_job_owner_prefix))",
    ]
    reservation_validation_function_lines = String[
        "metamdbg_file_mode() {",
        "if stat -c '%a' -- \"\$1\" >/dev/null 2>&1; then",
        "  stat -c '%a' -- \"\$1\"",
        "else",
        "  stat -f '%Lp' -- \"\$1\"",
        "fi",
        "}",
        "metamdbg_file_uid() {",
        "if stat -c '%u' -- \"\$1\" >/dev/null 2>&1; then",
        "  stat -c '%u' -- \"\$1\"",
        "else",
        "  stat -f '%u' -- \"\$1\"",
        "fi",
        "}",
        "validate_submission_reservation() {",
        "if [ ! -d \"\$submission_reservation_dir\" ] || [ -L \"\$submission_reservation_dir\" ]; then",
        "  echo \"metaMDBG submission reservation is missing or not a regular directory\" >&2",
        "  return 1",
        "fi",
        "if [ \"\$(metamdbg_file_uid \"\$submission_reservation_dir\")\" != \"\$(id -u)\" ] || [ \"\$(metamdbg_file_mode \"\$submission_reservation_dir\")\" != \"700\" ]; then",
        "  echo \"metaMDBG submission reservation owner or mode is invalid\" >&2",
        "  return 1",
        "fi",
        "if [ ! -f \"\$submission_reservation_contract\" ] || [ -L \"\$submission_reservation_contract\" ] || [ ! -s \"\$submission_reservation_contract\" ]; then",
        "  echo \"metaMDBG submission reservation contract is missing, empty, or not regular\" >&2",
        "  return 1",
        "fi",
        "if [ \"\$(metamdbg_file_uid \"\$submission_reservation_contract\")\" != \"\$(id -u)\" ] || [ \"\$(metamdbg_file_mode \"\$submission_reservation_contract\")\" != \"600\" ]; then",
        "  echo \"metaMDBG submission reservation contract owner or mode is invalid\" >&2",
        "  return 1",
        "fi",
        "cmp -s -- \"\$expected_reservation\" \"\$submission_reservation_contract\" || {",
        "  echo \"metaMDBG submission reservation contract does not match this job\" >&2",
        "  return 1",
        "}",
        "if [ ! -f \"\$submission_output_root_reservation\" ] || [ -L \"\$submission_output_root_reservation\" ] || [ ! -s \"\$submission_output_root_reservation\" ]; then",
        "  echo \"metaMDBG shared output-root reservation is missing, empty, or not regular\" >&2",
        "  return 1",
        "fi",
        "if [ \"\$(metamdbg_file_uid \"\$submission_output_root_reservation\")\" != \"\$(id -u)\" ] || [ \"\$(metamdbg_file_mode \"\$submission_output_root_reservation\")\" != \"600\" ]; then",
        "  echo \"metaMDBG shared output-root reservation owner or mode is invalid\" >&2",
        "  return 1",
        "fi",
        "cmp -s -- \"\$expected_output_root_reservation\" \"\$submission_output_root_reservation\" || {",
        "  echo \"metaMDBG shared output-root reservation does not match this owner capability\" >&2",
        "  return 1",
        "}",
        "if [ -z \"\${SLURM_JOB_ID:-}\" ]; then",
        "  echo \"metaMDBG bound submission reservation requires SLURM_JOB_ID\" >&2",
        "  return 1",
        "fi",
        "pending_submission_job=\"\${pending_submission_job_owner_prefix}\${SLURM_JOB_ID}.json\"",
        "pending_cluster_submission_job=\"\${pending_submission_job_owner_prefix}\${SLURM_JOB_ID}.\${SLURM_CLUSTER_NAME:-invalid/cluster}.json\"",
        "if [ -e \"\$pending_submission_job\" ] || [ -L \"\$pending_submission_job\" ] || [ -e \"\$pending_cluster_submission_job\" ] || [ -L \"\$pending_cluster_submission_job\" ]; then",
        "  echo \"metaMDBG runtime refuses pending scheduler job evidence; complete exact binding recovery before scheduler release\" >&2",
        "  return 1",
        "fi",
        "collect_metamdbg_directory_entries \"\$submission_reservation_dir\" \"\$reservation_entries_file\" \"submission reservation\"",
        "reservation_entry_count=0",
        "while IFS= read -r -d '' reservation_entry; do",
        "  reservation_entry_count=\$((reservation_entry_count + 1))",
        "done < \"\$reservation_entries_file\"",
        "if [ \"\$reservation_entry_count\" -ne 2 ] || [ ! -f \"\$submission_reservation_job\" ] || [ -L \"\$submission_reservation_job\" ] || [ ! -s \"\$submission_reservation_job\" ]; then",
        "  echo \"metaMDBG runtime refuses an unbound or invalid submission reservation\" >&2",
        "  return 1",
        "fi",
        "if [ \"\$(metamdbg_file_uid \"\$submission_reservation_job\")\" != \"\$(id -u)\" ] || [ \"\$(metamdbg_file_mode \"\$submission_reservation_job\")\" != \"600\" ]; then",
        "  echo \"metaMDBG submission reservation job record owner or mode is invalid\" >&2",
        "  return 1",
        "fi",
        "printf '%s' $(_metamdbg_shell_literal(first(job_record_parts))) > \"\$expected_reservation_job\"",
        "printf '%s' \"\$SLURM_JOB_ID\" >> \"\$expected_reservation_job\"",
        "printf '%s' $(_metamdbg_shell_literal(last(job_record_parts))) >> \"\$expected_reservation_job\"",
        "if ! cmp -s -- \"\$expected_reservation_job\" \"\$submission_reservation_job\"; then",
        "  if [[ ! \"\${SLURM_CLUSTER_NAME:-}\" =~ ^[A-Za-z0-9_.-]+\$ ]]; then",
        "    echo \"metaMDBG federated submission reservation requires canonical SLURM_CLUSTER_NAME\" >&2",
        "    return 1",
        "  fi",
        "  printf '%s' $(_metamdbg_shell_literal(first(job_record_parts))) > \"\$expected_reservation_job\"",
        "  printf '%s' \"\$SLURM_JOB_ID\" >> \"\$expected_reservation_job\"",
        "  printf '%s' $(_metamdbg_shell_literal(job_record_cluster_prefix)) >> \"\$expected_reservation_job\"",
        "  printf '%s' \"\$SLURM_CLUSTER_NAME\" >> \"\$expected_reservation_job\"",
        "  printf '%s' $(_metamdbg_shell_literal(job_record_cluster_suffix)) >> \"\$expected_reservation_job\"",
        "  cmp -s -- \"\$expected_reservation_job\" \"\$submission_reservation_job\" || {",
        "    echo \"metaMDBG submission reservation job record does not match exact SLURM job and cluster\" >&2",
        "    return 1",
        "  }",
        "fi",
        "}",
    ]
    output_root_identity_function_lines = String[
        "metamdbg_output_root_identity() {",
        "  local requested_root=\"\$1\"",
    ]
    reservation_identity_root = outputs.outdir
    while true
        push!(output_root_identity_function_lines,
            "  if [ \"\$requested_root\" = " *
            "$(_metamdbg_shell_literal(reservation_identity_root)) ]; then",
            "    printf '%s' " *
            _metamdbg_shell_literal(_output_root_reservation_identity(
                reservation_identity_root,
            )),
            "    return 0",
            "  fi",
        )
        reservation_identity_parent = dirname(reservation_identity_root)
        reservation_identity_parent == reservation_identity_root && break
        reservation_identity_root = reservation_identity_parent
    end
    append!(output_root_identity_function_lines, String[
        "  echo \"metaMDBG cannot identify an unexpected output-root domain: \$requested_root\" >&2",
        "  return 1",
        "}",
    ])
    output_root_domain_function_lines = String[
        output_root_identity_function_lines...,
        "require_owned_runtime_output_root_reservation() {",
        "  local observed_runtime_reservation_identity",
        "  if [ ! -d \"\$runtime_output_root_reservation\" ] || [ -L \"\$runtime_output_root_reservation\" ]; then",
        "    echo \"metaMDBG runtime output-root reservation is missing, replaced, or not a regular directory\" >&2",
        "    return 1",
        "  fi",
        "  observed_runtime_reservation_identity=\$(metamdbg_directory_identity \"\$runtime_output_root_reservation\") || observed_runtime_reservation_identity=",
        "  if [ -z \"\$observed_runtime_reservation_identity\" ] || [ \"\$observed_runtime_reservation_identity\" != \"\$bound_runtime_output_root_reservation_identity\" ]; then",
        "    echo \"metaMDBG runtime output-root reservation identity changed after acquisition\" >&2",
        "    return 1",
        "  fi",
        "}",
        "publish_metamdbg_runtime_cleanup_reservation() {",
        "  if [ -e \"\$runtime_cleanup_failure_reservation\" ] || [ -L \"\$runtime_cleanup_failure_reservation\" ]; then",
        "    return 0",
        "  fi",
        "  if mkdir -m 700 -- \"\$runtime_cleanup_failure_reservation\" 2>/dev/null; then",
        "    return 0",
        "  fi",
        "  if [ -e \"\$runtime_cleanup_failure_reservation\" ] || [ -L \"\$runtime_cleanup_failure_reservation\" ]; then",
        "    return 0",
        "  fi",
        "  echo \"metaMDBG could not publish its runtime cleanup failure reservation\" >&2",
        "  return 1",
        "}",
        "collect_metamdbg_directory_entries() {",
        "  local enumeration_root=\"\$1\"",
        "  local inventory_path=\"\$2\"",
        "  local inventory_label=\"\$3\"",
        "  local enumeration_mode=\"\${4:-one-level}\"",
        "  local reservation_lock_prefix=\"\${5:-}\"",
        "  local reservation_durable_prefix=\"\${6:-}\"",
        "  local reservation_pending_prefix=\"\${7:-}\"",
        "  local attempt_inventory=\"\${inventory_path}.attempt\"",
        "  local attempt_errors=\"\${inventory_path}.errors\"",
        "  local enumeration_attempt=1",
        "  local last_inventory_byte",
        "  while [ \"\$enumeration_attempt\" -le 3 ]; do",
        "    rm -f -- \"\$attempt_inventory\" \"\$attempt_errors\"",
        "    if [ ! -d \"\$enumeration_root\" ] || [ -L \"\$enumeration_root\" ]; then",
        "      rm -f -- \"\$inventory_path\"",
        "      echo \"metaMDBG cannot enumerate \$inventory_label because its root is missing, replaced, or not a regular directory: \$enumeration_root\" >&2",
        "      return 1",
        "    fi",
        "    if \"\$conda_runner\" run -n \"\$environment_name\" python -c $(directory_enumerator) \"\$enumeration_root\" \"\$attempt_inventory\" \"\$enumeration_mode\" \"\$reservation_lock_prefix\" \"\$reservation_durable_prefix\" \"\$reservation_pending_prefix\" 2> \"\$attempt_errors\"; then",
        "      if [ ! -f \"\$attempt_inventory\" ] || [ -L \"\$attempt_inventory\" ]; then",
        "        printf '%s\\n' \"metaMDBG enumerator exited successfully without a regular inventory\" > \"\$attempt_errors\"",
        "      elif [ -s \"\$attempt_inventory\" ]; then",
        "        last_inventory_byte=\$(tail -c 1 \"\$attempt_inventory\" | od -An -t u1 | tr -d '[:space:]')",
        "        if [ \"\$last_inventory_byte\" = 0 ]; then",
        "          mv -- \"\$attempt_inventory\" \"\$inventory_path\"",
        "          rm -f -- \"\$attempt_errors\"",
        "          return 0",
        "        fi",
        "        printf '%s\\n' \"metaMDBG successful inventory contains an unterminated NUL record\" > \"\$attempt_errors\"",
        "      else",
        "        mv -- \"\$attempt_inventory\" \"\$inventory_path\"",
        "        rm -f -- \"\$attempt_errors\"",
        "        return 0",
        "      fi",
        "    fi",
        "    if [ \"\$enumeration_attempt\" -eq 3 ] && [ -s \"\$attempt_errors\" ]; then",
        "      cat \"\$attempt_errors\" >&2",
        "    fi",
        "    rm -f -- \"\$attempt_inventory\" \"\$attempt_errors\"",
        "    enumeration_attempt=\$((enumeration_attempt + 1))",
        "  done",
        "  rm -f -- \"\$inventory_path\"",
        "  echo \"metaMDBG could not completely enumerate \$inventory_label after 3 attempts; refusing execution\" >&2",
        "  return 1",
        "}",
        "require_exclusive_output_root_domain() {",
        "  local reservation_root=\"\$outdir\"",
        "  local reservation_parent",
        "  local root_identity",
        "  local pid_reservation",
        "  local durable_prefix",
        "  local pending_prefix",
        "  local candidate",
        "  local inventory_path",
        "  while :; do",
        "    reservation_parent=\$(dirname -- \"\$reservation_root\")",
        "    root_identity=\$(metamdbg_output_root_identity \"\$reservation_root\")",
        "    if [ \"\$reservation_parent\" = / ]; then",
        "      pid_reservation=\"/$(_OUTPUT_ROOT_RESERVATION_LOCK_PREFIX)\$root_identity\"",
        "      durable_prefix=\"/$(_OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX)\$root_identity.\"",
        "    else",
        "      pid_reservation=\"\$reservation_parent/$(_OUTPUT_ROOT_RESERVATION_LOCK_PREFIX)\$root_identity\"",
        "      durable_prefix=\"\$reservation_parent/$(_OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX)\$root_identity.\"",
        "    fi",
        "    output_root_scan_index=\$((output_root_scan_index + 1))",
        "    inventory_path=\"\$secure_tmpdir/output-root-domain-\$output_root_scan_index.paths\"",
        "    pending_prefix=",
        "    if [ \"\$reservation_parent\" = \"\$outdir_parent\" ]; then",
        "      pending_prefix=\"\${pending_submission_job_prefix##*/}\"",
        "    fi",
        "    collect_metamdbg_directory_entries \"\$reservation_parent\" \"\$inventory_path\" \"same-root or ancestor output reservations\" \"reservations-one-level\" \"\${pid_reservation##*/}\" \"\${durable_prefix##*/}\" \"\$pending_prefix\"",
        "    while IFS= read -r -d '' candidate; do",
        "      if [[ \"\$candidate\" == \"\$pending_submission_job_prefix\"* ]]; then",
        "        echo \"metaMDBG runtime refuses ambiguous pending scheduler job evidence before owner consumption\" >&2",
        "        return 1",
        "      fi",
        "      if [ \"\$candidate\" = \"\$pid_reservation\" ]; then",
        "        echo \"metaMDBG overlaps an active shared output-root reservation: \$candidate\" >&2",
        "        return 1",
        "      fi",
        "      if [[ \"\$candidate\" == \"\$durable_prefix\"* ]]; then",
        "        if [ \"\$candidate\" = \"\$submission_output_root_reservation\" ]; then",
        "          continue",
        "        fi",
        "        if [ \"\$candidate\" = \"\$runtime_output_root_reservation\" ]; then",
        "          require_owned_runtime_output_root_reservation",
        "          continue",
        "        fi",
        "        echo \"metaMDBG overlaps an active shared output-root reservation: \$candidate\" >&2",
        "        return 1",
        "      fi",
        "    done < \"\$inventory_path\"",
        "    if [ \"\$reservation_parent\" = \"\$reservation_root\" ]; then",
        "      break",
        "    fi",
        "    reservation_root=\"\$reservation_parent\"",
        "  done",
        "  if [ -L \"\$outdir\" ] || { [ -e \"\$outdir\" ] && [ ! -d \"\$outdir\" ]; }; then",
        "    echo \"metaMDBG cannot enumerate descendant output reservations because outdir is replaced or not a regular directory\" >&2",
        "    return 1",
        "  fi",
        "  if [ -d \"\$outdir\" ]; then",
        "    output_root_scan_index=\$((output_root_scan_index + 1))",
        "    inventory_path=\"\$secure_tmpdir/output-root-domain-\$output_root_scan_index.paths\"",
        "    collect_metamdbg_directory_entries \"\$outdir\" \"\$inventory_path\" \"descendant output reservations\" \"reservations-recursive\" \"$(_OUTPUT_ROOT_RESERVATION_LOCK_PREFIX)\" \"$(_OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX)\"",
        "    while IFS= read -r -d '' candidate; do",
        "      if [ \"\$candidate\" = \"\$submission_output_root_reservation\" ]; then",
        "        continue",
        "      fi",
        "      if [ \"\$candidate\" = \"\$runtime_output_root_reservation\" ]; then",
        "        require_owned_runtime_output_root_reservation",
        "        continue",
        "      fi",
        "      echo \"metaMDBG overlaps an active descendant shared output-root reservation: \$candidate\" >&2",
        "      return 1",
        "    done < \"\$inventory_path\"",
        "  fi",
        "}",
    ]
    reservation_consumption_lines = String[
        "validate_submission_reservation",
        "if [ -e \"\$runtime_submission_reservation_dir\" ] || [ -L \"\$runtime_submission_reservation_dir\" ] || [ -e \"\$consumed_submission_reservation_dir\" ] || [ -L \"\$consumed_submission_reservation_dir\" ]; then",
        "  echo \"metaMDBG runtime or consumed owner state already exists for this capability\" >&2",
        "  exit 1",
        "fi",
        "mv -- \"\$submission_reservation_dir\" \"\$runtime_submission_reservation_dir\"",
        "runtime_owner_record_acquired=1",
        "bound_runtime_owner_record_identity=\$(metamdbg_directory_identity \"\$runtime_submission_reservation_dir\") || bound_runtime_owner_record_identity=",
        "if [ -z \"\$bound_runtime_owner_record_identity\" ]; then",
        "  echo \"metaMDBG could not bind its durable runtime owner record identity\" >&2",
        "  exit 1",
        "fi",
        "fsync_file_and_parent \"\$runtime_submission_reservation_dir\"",
        (post_runtime_private_rename_hook === nothing ? "true" :
         post_runtime_private_rename_hook),
        "if ! rm -- \"\$submission_output_root_reservation\"; then",
        "  echo \"metaMDBG failed to consume its queued shared output-root reservation\" >&2",
        "  cleanup_safe_to_release_domain=0",
        "  exit 1",
        "fi",
        "fsync_file_and_parent \"\$runtime_submission_reservation_dir\"",
        "if [ -e \"\$submission_output_root_reservation\" ] || [ -L \"\$submission_output_root_reservation\" ]; then",
        "  echo \"metaMDBG queued shared output-root reservation reappeared after durable release\" >&2",
        "  cleanup_safe_to_release_domain=0",
        "  exit 1",
        "fi",
        "queued_output_root_reservation_released=1",
        (post_queued_marker_release_hook === nothing ? "true" :
         post_queued_marker_release_hook),
    ]
    inventory_canonicalizer = _metamdbg_shell_literal(
        _metamdbg_runtime_inventory_canonicalizer_python(),
    )
    durability_fsyncer =
        _metamdbg_shell_literal(_metamdbg_runtime_fsync_python())
    gfa_validator =
        _metamdbg_shell_literal(_metamdbg_runtime_gfa_validator_python())
    completion_publication_hook = post_completion_publication_hook === nothing ?
                                  String[] :
                                  String[post_completion_publication_hook]
    lines = String[
        "set -euo pipefail",
        "umask 077",
        "outdir=$(_metamdbg_shell_literal(outputs.outdir))",
        "outdir_parent=$(_metamdbg_shell_literal(outdir_parent))",
        "outdir_identity=$(_metamdbg_shell_literal(outdir_identity))",
        "lock_dir=$(_metamdbg_shell_literal(lock_path))",
        "contigs_plain=\"\$outdir/contigs.fasta\"",
        "contigs_gz=\"\$outdir/contigs.fasta.gz\"",
        "graph_alias=\"\$outdir/assemblyGraph_k$(graph_k).gfa\"",
        "contract_marker=\"\$outdir/$(_METAMDBG_CONTRACT_FILENAME)\"",
        "completion_marker=\"\$outdir/$(_METAMDBG_COMPLETION_FILENAME)\"",
        "environment_spec=$(_metamdbg_shell_literal(verified_environment_path))",
        "expected_spec_sha256=$(_metamdbg_shell_literal(METAMDBG_ENVIRONMENT_SPEC_SHA256))",
        "expected_metamdbg_version=$(_metamdbg_shell_literal(METAMDBG_VERSION))",
        "environment_name=$(_metamdbg_shell_literal(METAMDBG_ENV_NAME))",
        "conda_runner=$(_metamdbg_shell_literal(conda_runner))",
        "expected_input_contract_signature=$(_metamdbg_shell_literal(input_contract.signature))",
        "expected_workflow_signature=$(_metamdbg_shell_literal(effective_submission_reservation.workflow_signature))",
        reservation_declaration_lines...,
        input_declaration_lines...,
        "secure_tmpdir=",
        "bound_outdir_identity=",
        "bound_lock_identity=",
        "bound_runtime_output_root_reservation_identity=",
        "bound_runtime_owner_record_identity=",
        "lock_acquired=0",
        "output_root_runtime_reservation_acquired=0",
        "runtime_owner_record_acquired=0",
        "queued_output_root_reservation_released=0",
        "output_root_scan_index=0",
        "preserve_secure_tmpdir=0",
        "cleanup_safe_to_release_domain=1",
        "lock_retry_attempts=$(lock_retry_attempts)",
        "lock_retry_delay_seconds=$(Float64(lock_retry_delay_seconds))",
        "sha256_stream() {",
        "  if command -v sha256sum >/dev/null 2>&1; then",
        "    sha256sum | awk '{print \$1}'",
        "  elif command -v shasum >/dev/null 2>&1; then",
        "    shasum -a 256 | awk '{print \$1}'",
        "  else",
        "    echo \"metaMDBG requires sha256sum or shasum to validate content\" >&2",
        "    return 1",
        "  fi",
        "}",
        "sha256_file() {",
        "  local path=\"\$1\"",
        "  sha256_stream < \"\$path\"",
        "}",
        "sha256_text() {",
        "  printf '%s' \"\$1\" | sha256_stream",
        "}",
        "fsync_file_and_parent() {",
        "  if ! \"\$conda_runner\" run -n \"\$environment_name\" python -c $(durability_fsyncer) \"\$1\"; then",
        "    echo \"metaMDBG could not durably fsync lifecycle path and containing directory: \$1\" >&2",
        "    return 1",
        "  fi",
        "}",
        "metamdbg_directory_identity() {",
        "  local path=\"\$1\"",
        "  if stat -Lc '%d:%i' -- \"\$path\" >/dev/null 2>&1; then",
        "    stat -Lc '%d:%i' -- \"\$path\"",
        "  else",
        "    stat -f '%d:%i' -- \"\$path\"",
        "  fi",
        "}",
        "bind_metamdbg_output_root() {",
        "  if [ ! -d \"\$outdir\" ] || [ -L \"\$outdir\" ]; then",
        "    echo \"metaMDBG output root must be a regular directory before binding its physical identity: \$outdir\" >&2",
        "    return 1",
        "  fi",
        "  local observed_outdir",
        "  observed_outdir=\$(cd -- \"\$outdir\" && pwd -P)",
        "  if [ \"\$observed_outdir\" != \"\$outdir\" ]; then",
        "    echo \"metaMDBG output root is not at its canonical physical location: \$outdir resolves to \$observed_outdir\" >&2",
        "    return 1",
        "  fi",
        "  bound_outdir_identity=\$(metamdbg_directory_identity \"\$outdir\")",
        "  if [ -z \"\$bound_outdir_identity\" ]; then",
        "    echo \"metaMDBG could not bind the output root physical identity\" >&2",
        "    return 1",
        "  fi",
        "}",
        "require_unchanged_metamdbg_output_root() {",
        "  if [ ! -d \"\$outdir\" ] || [ -L \"\$outdir\" ]; then",
        "    echo \"metaMDBG output root physical identity changed after binding: \$outdir\" >&2",
        "    return 1",
        "  fi",
        "  local observed_outdir",
        "  local observed_identity",
        "  observed_outdir=\$(cd -- \"\$outdir\" && pwd -P)",
        "  observed_identity=\$(metamdbg_directory_identity \"\$outdir\")",
        "  if [ \"\$observed_outdir\" != \"\$outdir\" ] || [ \"\$observed_identity\" != \"\$bound_outdir_identity\" ]; then",
        "    echo \"metaMDBG output root physical identity changed after binding: \$outdir\" >&2",
        "    return 1",
        "  fi",
        "}",
        "require_metamdbg_artifact_containment() {",
        "  local path=\"\$1\"",
        "  local label=\"\$2\"",
        "  local physical_parent",
        "  require_unchanged_metamdbg_output_root",
        "  if [ ! -e \"\$path\" ] && [ ! -L \"\$path\" ]; then",
        "    echo \"\$label is missing from the bound metaMDBG output root: \$path\" >&2",
        "    return 1",
        "  fi",
        "  physical_parent=\$(cd -- \"\$(dirname -- \"\$path\")\" && pwd -P)",
        "  if [ \"\$physical_parent\" != \"\$outdir\" ]; then",
        "    echo \"\$label escaped the bound metaMDBG output root: \$path\" >&2",
        "    return 1",
        "  fi",
        "}",
        input_validation_function_lines...,
        reservation_validation_function_lines...,
        output_root_domain_function_lines...,
        "mkdir -p -- \"\$outdir_parent\"",
        "runtime_outdir_parent=\$(cd -- \"\$outdir_parent\" && pwd -P)",
        "if [ \"\$runtime_outdir_parent\" != \"\$outdir_parent\" ]; then",
        "  echo \"metaMDBG canonical output parent changed before execution\" >&2",
        "  exit 1",
        "fi",
        "cleanup_metamdbg_lifecycle() {",
        "  local status=\$?",
        "  local observed_lock_identity",
        "  trap - EXIT INT TERM HUP",
        "  if [ -n \"\$secure_tmpdir\" ]; then",
        "    if [ \"\$preserve_secure_tmpdir\" -eq 1 ]; then",
        "      echo \"metaMDBG retained secure lifecycle recovery state: \$secure_tmpdir\" >&2",
        "      cleanup_safe_to_release_domain=0",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    elif ! rm -rf -- \"\$secure_tmpdir\"; then",
        "      echo \"metaMDBG failed to remove secure lifecycle temporary directory\" >&2",
        "      cleanup_safe_to_release_domain=0",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    fi",
        "  fi",
        "  if [ \"\$output_root_runtime_reservation_acquired\" -eq 1 ] && ! require_owned_runtime_output_root_reservation; then",
        "    publish_metamdbg_runtime_cleanup_reservation || true",
        "    cleanup_safe_to_release_domain=0",
        "    [ \"\$status\" -ne 0 ] || status=1",
        "  fi",
        "  if [ \"\$lock_acquired\" -eq 1 ]; then",
        "    if [ \"\$cleanup_safe_to_release_domain\" -eq 1 ] && [ -d \"\$lock_dir\" ] && [ ! -L \"\$lock_dir\" ]; then",
        "      observed_lock_identity=\$(metamdbg_directory_identity \"\$lock_dir\") || observed_lock_identity=",
        "    else",
        "      observed_lock_identity=",
        "    fi",
        "    if [ \"\$cleanup_safe_to_release_domain\" -ne 1 ] || [ -z \"\$observed_lock_identity\" ] || [ \"\$observed_lock_identity\" != \"\$bound_lock_identity\" ] || ! rmdir -- \"\$lock_dir\"; then",
        "      echo \"metaMDBG failed to release its exact private lifecycle lock; retaining the shared runtime reservation fail-closed: \$lock_dir\" >&2",
        "      cleanup_safe_to_release_domain=0",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    elif [ -e \"\$lock_dir\" ] || [ -L \"\$lock_dir\" ]; then",
        "      echo \"metaMDBG private lifecycle lock reappeared after cleanup; retaining the shared runtime reservation fail-closed: \$lock_dir\" >&2",
        "      cleanup_safe_to_release_domain=0",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    else",
        "      lock_acquired=0",
        "    fi",
        "  fi",
        "  if [ \"\$output_root_runtime_reservation_acquired\" -eq 1 ]; then",
        "    if [ \"\$cleanup_safe_to_release_domain\" -eq 1 ] && [ \"\$lock_acquired\" -eq 0 ]; then",
        "      if ! require_owned_runtime_output_root_reservation; then",
        "        publish_metamdbg_runtime_cleanup_reservation || true",
        "        cleanup_safe_to_release_domain=0",
        "        [ \"\$status\" -ne 0 ] || status=1",
        "      elif rmdir -- \"\$runtime_output_root_reservation\"; then",
        "        if [ \"\$runtime_owner_record_acquired\" -eq 1 ]; then",
        "          fsync_file_and_parent \"\$runtime_submission_reservation_dir\"",
        "        else",
        "          fsync_file_and_parent \"\$submission_reservation_dir\"",
        "        fi",
        (post_runtime_marker_release_hook === nothing ? "        true" :
         "        $(post_runtime_marker_release_hook)"),
        "        output_root_runtime_reservation_acquired=0",
        "      else",
        "        publish_metamdbg_runtime_cleanup_reservation || true",
        "        echo \"metaMDBG failed to release shared output-root reservation: \$runtime_output_root_reservation\" >&2",
        "        [ \"\$status\" -ne 0 ] || status=1",
        "      fi",
        "    else",
        "      echo \"metaMDBG retained shared output-root reservation after inner cleanup failure: \$runtime_output_root_reservation\" >&2",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    fi",
        "  fi",
        "  if [ \"\$runtime_owner_record_acquired\" -eq 1 ]; then",
        "    if [ \"\$cleanup_safe_to_release_domain\" -eq 1 ] && [ \"\$queued_output_root_reservation_released\" -eq 1 ] && [ \"\$lock_acquired\" -eq 0 ] && [ \"\$output_root_runtime_reservation_acquired\" -eq 0 ]; then",
        "      observed_runtime_owner_identity=",
        "      if [ -d \"\$runtime_submission_reservation_dir\" ] && [ ! -L \"\$runtime_submission_reservation_dir\" ]; then",
        "        observed_runtime_owner_identity=\$(metamdbg_directory_identity \"\$runtime_submission_reservation_dir\") || observed_runtime_owner_identity=",
        "      fi",
        "      if [ -z \"\$observed_runtime_owner_identity\" ] || [ \"\$observed_runtime_owner_identity\" != \"\$bound_runtime_owner_record_identity\" ]; then",
        "        echo \"metaMDBG durable runtime owner record identity changed during cleanup\" >&2",
        "        [ \"\$status\" -ne 0 ] || status=1",
        "      elif [ -e \"\$consumed_submission_reservation_dir\" ] || [ -L \"\$consumed_submission_reservation_dir\" ]; then",
        "        echo \"metaMDBG refuses to overwrite an existing consumed owner record\" >&2",
        "        [ \"\$status\" -ne 0 ] || status=1",
        "      elif mv -- \"\$runtime_submission_reservation_dir\" \"\$consumed_submission_reservation_dir\"; then",
        (post_consumed_owner_rename_hook === nothing ? "        true" :
         "        $(post_consumed_owner_rename_hook)"),
        "        fsync_file_and_parent \"\$consumed_submission_reservation_dir\"",
        "        runtime_owner_record_acquired=0",
        "      else",
        "        echo \"metaMDBG failed to publish its durable consumed owner record\" >&2",
        "        [ \"\$status\" -ne 0 ] || status=1",
        "      fi",
        "    else",
        "      echo \"metaMDBG retained its durable runtime owner record after lifecycle cleanup failure: \$runtime_submission_reservation_dir\" >&2",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    fi",
        "  fi",
        "  exit \"\$status\"",
        "}",
        "trap cleanup_metamdbg_lifecycle EXIT",
        "trap 'exit 129' HUP",
        "trap 'exit 130' INT",
        "trap 'exit 143' TERM",
        "secure_tmpdir=\$(mktemp -d \"\$outdir_parent/.mycelia-metamdbg-tmp.\${outdir_identity}.XXXXXXXXXX\")",
        "chmod 700 \"\$secure_tmpdir\"",
        "expected_contract=\"\$secure_tmpdir/expected-contract.json\"",
        "contract_new=\"\$secure_tmpdir/contract.new\"",
        "contigs_new=\"\$secure_tmpdir/contigs.fasta.new.gz\"",
        "package_inventory_before=\"\$secure_tmpdir/conda-packages.before.json\"",
        "package_inventory_normalized=\"\$secure_tmpdir/conda-packages.before.tsv\"",
        "package_inventory_after=\"\$secure_tmpdir/conda-packages.after.json\"",
        "package_inventory_after_normalized=\"\$secure_tmpdir/conda-packages.after.tsv\"",
        "completion_payload=\"\$secure_tmpdir/completion-payload.json\"",
        "completion_new=\"\$secure_tmpdir/completion.new\"",
        "printf '%s' $(_metamdbg_shell_literal(input_contract.contents)) > \"\$expected_contract\"",
        "chmod 600 \"\$expected_contract\"",
        "expected_reservation=\"\$secure_tmpdir/expected-reservation.json\"",
        "expected_reservation_job=\"\$secure_tmpdir/expected-reservation-job.json\"",
        "expected_output_root_reservation=\"\$secure_tmpdir/expected-output-root-reservation.json\"",
        "reservation_entries_file=\"\$secure_tmpdir/submission-reservation.paths\"",
        "printf '%s' $(_metamdbg_shell_literal(effective_submission_reservation.contents)) > \"\$expected_reservation\"",
        "chmod 600 \"\$expected_reservation\"",
        "printf '%s' $(_metamdbg_shell_literal(effective_submission_reservation.output_root_reservation_contents)) > \"\$expected_output_root_reservation\"",
        "chmod 600 \"\$expected_output_root_reservation\"",
        "validate_submission_reservation",
        "if ! mkdir -m 700 -- \"\$runtime_output_root_reservation\" 2>/dev/null; then",
        "  echo \"metaMDBG could not acquire its exact runtime output-root reservation: \$runtime_output_root_reservation\" >&2",
        "  exit 1",
        "fi",
        "output_root_runtime_reservation_acquired=1",
        "bound_runtime_output_root_reservation_identity=\$(metamdbg_directory_identity \"\$runtime_output_root_reservation\") || bound_runtime_output_root_reservation_identity=",
        "if [ -z \"\$bound_runtime_output_root_reservation_identity\" ]; then",
        "  echo \"metaMDBG could not bind its runtime output-root reservation identity\" >&2",
        "  exit 1",
        "fi",
        "fsync_file_and_parent \"\$runtime_output_root_reservation\"",
        (post_runtime_marker_publication_hook === nothing ? "true" :
         post_runtime_marker_publication_hook),
        "require_exclusive_output_root_domain",
        "lock_attempt=1",
        "while ! mkdir -m 700 -- \"\$lock_dir\" 2>/dev/null; do",
        "  validate_submission_reservation",
        "  require_owned_runtime_output_root_reservation",
        "  if [ \"\$lock_attempt\" -ge \"\$lock_retry_attempts\" ]; then",
        "    echo \"metaMDBG output remained locked after \$lock_retry_attempts attempts: \$lock_dir\" >&2",
        "    exit 1",
        "  fi",
        "  sleep \"\$lock_retry_delay_seconds\"",
        "  lock_attempt=\$((lock_attempt + 1))",
        "done",
        "lock_acquired=1",
        "bound_lock_identity=\$(metamdbg_directory_identity \"\$lock_dir\")",
        "if [ -z \"\$bound_lock_identity\" ]; then",
        "  echo \"metaMDBG could not bind its private lifecycle lock identity\" >&2",
        "  exit 1",
        "fi",
        "require_exclusive_output_root_domain",
        reservation_consumption_lines...,
        "if [ -L \"\$outdir\" ]; then",
        "  echo \"metaMDBG outdir must not be a symbolic link: \$outdir\" >&2",
        "  exit 1",
        "fi",
        "if [ -e \"\$outdir\" ] && [ ! -d \"\$outdir\" ]; then",
        "  echo \"metaMDBG outdir exists but is not a directory: \$outdir\" >&2",
        "  exit 1",
        "fi",
        "mkdir -p -- \"\$outdir\"",
        "bind_metamdbg_output_root",
        "validate_metamdbg_inputs",
        "outdir_entries_file=\"\$secure_tmpdir/outdir.paths\"",
        "collect_metamdbg_directory_entries \"\$outdir\" \"\$outdir_entries_file\" \"metaMDBG output root\"",
        "outdir_entry_count=0",
        "partial_entry=",
        "while IFS= read -r -d '' outdir_entry; do",
        "  outdir_entry_count=\$((outdir_entry_count + 1))",
        "  if [ \"\$outdir_entry\" != \"\$contract_marker\" ] && [ -z \"\$partial_entry\" ]; then",
        "    partial_entry=\"\$outdir_entry\"",
        "  fi",
        "done < \"\$outdir_entries_file\"",
        "contract_exists=0",
        "if [ \"\$outdir_entry_count\" -gt 0 ]; then",
        "  if [ ! -f \"\$contract_marker\" ] || [ -L \"\$contract_marker\" ]; then",
        "    echo \"metaMDBG refuses to adopt a nonempty output root without its regular provenance contract marker\" >&2",
        "    exit 1",
        "  fi",
        "  test -s \"\$contract_marker\" || {",
        "    echo \"metaMDBG provenance contract marker is empty\" >&2",
        "    exit 1",
        "  }",
        "  cmp -s -- \"\$expected_contract\" \"\$contract_marker\" || {",
        "    echo \"metaMDBG existing output contract does not match this invocation\" >&2",
        "    exit 1",
        "  }",
        "  contract_exists=1",
        "fi",
        "if [ \"\$contract_exists\" -eq 1 ] && [ ! -e \"\$completion_marker\" ]; then",
        "  if [ -n \"\$partial_entry\" ]; then",
        "    echo \"metaMDBG refuses partial contracted output without realized-stage completion provenance; use a fresh output root\" >&2",
        "    exit 1",
        "  fi",
        "fi",
        "if [ ! -f \"\$environment_spec\" ] || [ -L \"\$environment_spec\" ]; then",
        "  echo \"metaMDBG environment specification is missing or not regular\" >&2",
        "  exit 1",
        "fi",
        "actual_spec_sha256=\$(sha256_file \"\$environment_spec\")",
        "if [ \"\$actual_spec_sha256\" != \"\$expected_spec_sha256\" ]; then",
        "  echo \"metaMDBG environment specification checksum mismatch\" >&2",
        "  exit 1",
        "fi",
        "capture_package_inventory() {",
        "  local json_path=\"\$1\"",
        "  local normalized_path=\"\$2\"",
        "  require_unchanged_metamdbg_output_root",
        "  \"\$conda_runner\" list -n \"\$environment_name\" --json > \"\$json_path\"",
        "  require_unchanged_metamdbg_output_root",
        "  \"\$conda_runner\" run -n \"\$environment_name\" python -c $(inventory_canonicalizer) \"\$json_path\" \"\$normalized_path\" || {",
        "    echo \"metaMDBG environment package inventory is incomplete or malformed\" >&2",
        "    return 1",
        "  }",
        "  require_unchanged_metamdbg_output_root",
        "}",
        "capture_package_inventory \"\$package_inventory_before\" \"\$package_inventory_normalized\"",
        "package_count=\$(awk 'END { print NR }' \"\$package_inventory_normalized\")",
        "if [ \"\$package_count\" -le 0 ]; then",
        "  echo \"metaMDBG normalized package inventory is empty\" >&2",
        "  exit 1",
        "fi",
        "if cut -f 1 \"\$package_inventory_normalized\" | uniq -d | grep -q .; then",
        "  echo \"metaMDBG normalized package inventory has duplicate package names\" >&2",
        "  exit 1",
        "fi",
        "awk -F '\\t' -v expected=\"\$expected_metamdbg_version\" '",
        "  \$1 == \"metamdbg\" && \$2 == expected { count += 1 }",
        "  END { exit count == 1 ? 0 : 1 }",
        "' \"\$package_inventory_normalized\" || {",
        "  echo \"metaMDBG environment must contain one metamdbg record at exactly \$expected_metamdbg_version\" >&2",
        "  exit 1",
        "}",
        "package_inventory_sha256=\$(sha256_file \"\$package_inventory_normalized\")",
        "require_unchanged_metamdbg_output_root",
        staged_input_creation_lines[1:2]...,
        staged_input_declaration_lines...,
        staged_input_creation_lines[3:end]...,
        "require_unchanged_metamdbg_output_root",
        staged_input_validation_lines...,
        "validate_fasta_stream() {",
        "  awk '",
        "    BEGIN { records = 0; sequence_bases = 0 }",
        "    { sub(/\\r\$/, \"\", \$0) }",
        "    /^>/ {",
        "      if (records > 0 && sequence_bases == 0) exit 10",
        "      identifier = substr(\$0, 2)",
        "      if (identifier == \"\" || identifier ~ /^[[:space:]]/) exit 14",
        "      sub(/[[:space:]].*\$/, \"\", identifier)",
        "      if (identifier_seen[identifier]++) exit 15",
        "      records += 1",
        "      sequence_bases = 0",
        "      next",
        "    }",
        "    /^[[:space:]]*\$/ { next }",
        "    {",
        "      if (records == 0) exit 11",
        "      if (\$0 !~ /^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+\$/) exit 12",
        "      sequence_bases += length(\$0)",
        "    }",
        "    END { if (records == 0 || sequence_bases == 0) exit 13 }",
        "  ' || return 1",
        "}",
        "validate_contigs() {",
        "  local path=\$1",
        "  require_unchanged_metamdbg_output_root",
        "  if [ ! -f \"\$path\" ] || [ -L \"\$path\" ] || [ ! -s \"\$path\" ]; then",
        "    echo \"metaMDBG contigs are missing, empty, or not regular: \$path\" >&2",
        "    return 1",
        "  fi",
        "  if [[ \"\$path\" == *.gz ]]; then",
        "    gzip -cd -- \"\$path\" | validate_fasta_stream || {",
        "      echo \"metaMDBG contigs are not sequence-bearing FASTA: \$path\" >&2",
        "      return 1",
        "    }",
        "  else",
        "    validate_fasta_stream < \"\$path\" || {",
        "      echo \"metaMDBG contigs are not sequence-bearing FASTA: \$path\" >&2",
        "      return 1",
        "    }",
        "  fi",
        "  require_unchanged_metamdbg_output_root",
        "}",
        "validate_gfa() {",
        "  local path=\$1",
        "  require_unchanged_metamdbg_output_root",
        "  if [ ! -f \"\$path\" ] || [ -L \"\$path\" ] || [ ! -s \"\$path\" ]; then",
        "    echo \"metaMDBG graph is missing, empty, or not regular: \$path\" >&2",
        "    return 1",
        "  fi",
        "  \"\$conda_runner\" run -n \"\$environment_name\" python -c $(gfa_validator) \"\$path\" || {",
        "    echo \"metaMDBG graph has malformed, unknown, duplicate, typed-tag, version, or dangling GFA1 records: \$path\" >&2",
        "    return 1",
        "  }",
        "  require_unchanged_metamdbg_output_root",
        "}",
        "find_metamdbg_graph() {",
        "  local all_candidates=()",
        "  local requested_candidates=()",
        "  local candidate",
        "  local filename",
        "  shopt -s nullglob",
        "  for candidate in \"\$outdir\"/assemblyGraph_k*_*bps.gfa; do",
        "    filename=\${candidate##*/}",
        "    if [[ ! \"\$filename\" =~ ^assemblyGraph_k[0-9]+_.*bps\\.gfa\$ ]]; then",
        "      continue",
        "    fi",
        "    if [ ! -f \"\$candidate\" ] || [ -L \"\$candidate\" ] || [ ! -s \"\$candidate\" ]; then",
        "      echo \"metaMDBG graph candidate is empty or not regular: \$candidate\" >&2",
        "      shopt -u nullglob",
        "      return 2",
        "    fi",
        "    all_candidates+=(\"\$candidate\")",
        "    if [[ \"\$filename\" == assemblyGraph_k$(graph_k)_*bps.gfa ]]; then",
        "      requested_candidates+=(\"\$candidate\")",
        "    fi",
        "  done",
        "  shopt -u nullglob",
        "  if [ \"\${#all_candidates[@]}\" -eq 0 ]; then",
        "    return 1",
        "  elif [ \"\${#all_candidates[@]}\" -ne 1 ]; then",
        "    echo \"metaMDBG requires exactly one total graph artifact per output root\" >&2",
        "    return 2",
        "  elif [ \"\${#requested_candidates[@]}\" -ne 1 ]; then",
        "    echo \"metaMDBG graph inventory does not match requested k=$(graph_k)\" >&2",
        "    return 2",
        "  fi",
        "  printf '%s\\n' \"\${requested_candidates[0]}\"",
        "}",
        "if [ -e \"\$contigs_gz\" ] || [ -L \"\$contigs_gz\" ]; then",
        "  validate_contigs \"\$contigs_gz\"",
        "elif [ -e \"\$contigs_plain\" ] || [ -L \"\$contigs_plain\" ]; then",
        "  validate_contigs \"\$contigs_plain\"",
        "else",
        "  require_unchanged_metamdbg_output_root",
        "  $(runtime_asm_cmd)",
        "  require_unchanged_metamdbg_output_root",
        "fi",
        "if [ -e \"\$contigs_gz\" ] || [ -L \"\$contigs_gz\" ]; then",
        "  validate_contigs \"\$contigs_gz\"",
        "else",
        "  validate_contigs \"\$contigs_plain\"",
        "  gzip -c -- \"\$contigs_plain\" > \"\$contigs_new\"",
        "  chmod 600 \"\$contigs_new\"",
        "  validate_contigs \"\$contigs_new\"",
        "  require_unchanged_metamdbg_output_root",
        "  mv -f -- \"\$contigs_new\" \"\$contigs_gz\"",
        "  require_unchanged_metamdbg_output_root",
        "fi",
        "validate_contigs \"\$contigs_gz\"",
        "graph_status=0",
        "graph_source=\"\"",
        "graph_source=\$(find_metamdbg_graph) || graph_status=\$?",
        "if [ \"\$graph_status\" -eq 1 ]; then",
        "  require_unchanged_metamdbg_output_root",
        "  $(runtime_gfa_cmd)",
        "  require_unchanged_metamdbg_output_root",
        "  graph_status=0",
        "  graph_source=\$(find_metamdbg_graph) || graph_status=\$?",
        "fi",
        "if [ \"\$graph_status\" -ne 0 ]; then",
        "  echo \"metaMDBG produced no unique nonempty graph for k=$(graph_k)\" >&2",
        "  exit \"\$graph_status\"",
        "fi",
        "validate_gfa \"\$graph_source\"",
        "require_unchanged_metamdbg_output_root",
        "require_metamdbg_artifact_containment \"\$contigs_gz\" \"metaMDBG contigs artifact\"",
        "require_metamdbg_artifact_containment \"\$graph_source\" \"metaMDBG graph artifact\"",
        "rm -f -- \"\$graph_alias\"",
        "ln -s -- \"\$(basename \"\$graph_source\")\" \"\$graph_alias\"",
        "require_unchanged_metamdbg_output_root",
        "test -L \"\$graph_alias\" && test -f \"\$graph_alias\" && test -s \"\$graph_alias\"",
        "validate_contigs \"\$contigs_gz\"",
        "validate_gfa \"\$graph_source\"",
        "validate_staged_metamdbg_inputs",
        "validate_metamdbg_inputs",
        "capture_package_inventory \"\$package_inventory_after\" \"\$package_inventory_after_normalized\"",
        "cmp -s -- \"\$package_inventory_normalized\" \"\$package_inventory_after_normalized\" || {",
        "  echo \"metaMDBG resolved package inventory changed while assembly tools were running\" >&2",
        "  exit 1",
        "}",
        "if [ -e \"\$completion_marker\" ] || [ -L \"\$completion_marker\" ]; then",
        "  echo \"metaMDBG refuses to overwrite an existing completion manifest\" >&2",
        "  exit 1",
        "fi",
        "require_unchanged_metamdbg_output_root",
        "require_metamdbg_artifact_containment \"\$contigs_gz\" \"metaMDBG contigs artifact before provenance publication\"",
        "require_metamdbg_artifact_containment \"\$graph_source\" \"metaMDBG graph artifact before provenance publication\"",
        "contigs_canonical=\"\$outdir/contigs.fasta.gz\"",
        "graph_source_parent=\$(cd -- \"\$(dirname -- \"\$graph_source\")\" && pwd -P)",
        "graph_source_canonical=\"\$graph_source_parent/\$(basename -- \"\$graph_source\")\"",
        "outdir_sha256=\$(sha256_text \"\$outdir\")",
        "contigs_path_sha256=\$(sha256_text \"\$contigs_canonical\")",
        "graph_alias_path_sha256=\$(sha256_text \"\$graph_alias\")",
        "graph_source_path_sha256=\$(sha256_text \"\$graph_source_canonical\")",
        "contigs_size=\$(wc -c < \"\$contigs_canonical\" | tr -d '[:space:]')",
        "graph_size=\$(wc -c < \"\$graph_source_canonical\" | tr -d '[:space:]')",
        "contigs_sha256=\$(sha256_file \"\$contigs_canonical\")",
        "graph_sha256=\$(sha256_file \"\$graph_source_canonical\")",
        "printf '%s' '{\"schema_version\":$(_METAMDBG_COMPLETION_SCHEMA_VERSION),\"workflow\":{' > \"\$completion_payload\"",
        "printf '\"canonical_outdir_sha256\":\"%s\",' \"\$outdir_sha256\" >> \"\$completion_payload\"",
        "printf '\"input_contract_signature\":\"%s\",' \"\$expected_input_contract_signature\" >> \"\$completion_payload\"",
        "printf '\"graph_k\":%s,' $(graph_k) >> \"\$completion_payload\"",
        "printf '\"workflow_signature\":\"%s\"},\"toolchain\":{' \"\$expected_workflow_signature\" >> \"\$completion_payload\"",
        "printf '\"metamdbg_version\":\"%s\",' \"\$expected_metamdbg_version\" >> \"\$completion_payload\"",
        "printf '\"environment_name\":\"%s\",' \"\$environment_name\" >> \"\$completion_payload\"",
        "printf '\"environment_spec_sha256\":\"%s\",' \"\$expected_spec_sha256\" >> \"\$completion_payload\"",
        "printf '\"package_inventory_sha256\":\"%s\",' \"\$package_inventory_sha256\" >> \"\$completion_payload\"",
        "printf '\"package_count\":%s},\"artifacts\":{\"contigs\":{' \"\$package_count\" >> \"\$completion_payload\"",
        "printf '\"canonical_path_sha256\":\"%s\",' \"\$contigs_path_sha256\" >> \"\$completion_payload\"",
        "printf '\"size_bytes\":%s,' \"\$contigs_size\" >> \"\$completion_payload\"",
        "printf '\"sha256\":\"%s\"},\"graph\":{' \"\$contigs_sha256\" >> \"\$completion_payload\"",
        "printf '\"alias_path_sha256\":\"%s\",' \"\$graph_alias_path_sha256\" >> \"\$completion_payload\"",
        "printf '\"source_path_sha256\":\"%s\",' \"\$graph_source_path_sha256\" >> \"\$completion_payload\"",
        "printf '\"size_bytes\":%s,' \"\$graph_size\" >> \"\$completion_payload\"",
        "printf '\"sha256\":\"%s\"}}}' \"\$graph_sha256\" >> \"\$completion_payload\"",
        "completion_signature=\$(sha256_file \"\$completion_payload\")",
        "printf '%s' '{\"schema_version\":$(_METAMDBG_COMPLETION_SCHEMA_VERSION),\"signature_algorithm\":\"sha256\",' > \"\$completion_new\"",
        "printf '\"signature\":\"%s\",\"manifest\":' \"\$completion_signature\" >> \"\$completion_new\"",
        "cat \"\$completion_payload\" >> \"\$completion_new\"",
        "printf '}\\n' >> \"\$completion_new\"",
        "chmod 600 \"\$completion_new\"",
        "expected_completion_sha256=\$(sha256_file \"\$completion_new\")",
        "validate_metamdbg_inputs",
        "if [ \"\$contract_exists\" -eq 1 ]; then",
        "  require_unchanged_metamdbg_output_root",
        "  cmp -s -- \"\$expected_contract\" \"\$contract_marker\"",
        "else",
        "  if [ -e \"\$contract_marker\" ] || [ -L \"\$contract_marker\" ]; then",
        "    echo \"metaMDBG provenance marker appeared during locked execution\" >&2",
        "    exit 1",
        "  fi",
        "  cp -- \"\$expected_contract\" \"\$contract_new\"",
        "  chmod 600 \"\$contract_new\"",
        "  fsync_file_and_parent \"\$contract_new\"",
        "  require_unchanged_metamdbg_output_root",
        "  mv -n -- \"\$contract_new\" \"\$contract_marker\"",
        "  fsync_file_and_parent \"\$contract_marker\"",
        "  require_unchanged_metamdbg_output_root",
        "  if [ -e \"\$contract_new\" ]; then",
        "    echo \"metaMDBG refused to overwrite a concurrent provenance marker\" >&2",
        "    exit 1",
        "  fi",
        "  cmp -s -- \"\$expected_contract\" \"\$contract_marker\"",
        "fi",
        "require_unchanged_metamdbg_output_root",
        "require_metamdbg_artifact_containment \"\$contigs_gz\" \"metaMDBG contigs artifact before completion publication\"",
        "require_metamdbg_artifact_containment \"\$graph_source\" \"metaMDBG graph artifact before completion publication\"",
        "fsync_file_and_parent \"\$completion_new\"",
        "mv -n -- \"\$completion_new\" \"\$completion_marker\"",
        "fsync_file_and_parent \"\$completion_marker\"",
        "require_unchanged_metamdbg_output_root",
        "if [ -e \"\$completion_new\" ]; then",
        "  echo \"metaMDBG refused to overwrite a concurrent completion manifest\" >&2",
        "  exit 1",
        "fi",
        "if [ ! -f \"\$completion_marker\" ] || [ -L \"\$completion_marker\" ] || [ ! -s \"\$completion_marker\" ]; then",
        "  echo \"metaMDBG failed to publish a regular completion manifest\" >&2",
        "  exit 1",
        "fi",
        completion_publication_hook...,
        "require_unchanged_metamdbg_output_root",
        "actual_completion_sha256=\$(sha256_file \"\$completion_marker\")",
        "if [ \"\$actual_completion_sha256\" != \"\$expected_completion_sha256\" ]; then",
        "  rm -f -- \"\$completion_marker\"",
        "  echo \"metaMDBG published completion manifest changed unexpectedly\" >&2",
        "  exit 1",
        "fi",
        "published_graph_source_parent=\$(cd -- \"\$(dirname -- \"\$graph_alias\")\" && pwd -P)",
        "published_graph_source_canonical=\"\$published_graph_source_parent/\$(readlink -- \"\$graph_alias\")\"",
        "published_contigs_size=\$(wc -c < \"\$contigs_canonical\" | tr -d '[:space:]')",
        "published_graph_size=\$(wc -c < \"\$published_graph_source_canonical\" | tr -d '[:space:]')",
        "published_contigs_sha256=\$(sha256_file \"\$contigs_canonical\")",
        "published_graph_sha256=\$(sha256_file \"\$published_graph_source_canonical\")",
        "if [ \"\$published_graph_source_canonical\" != \"\$graph_source_canonical\" ] || [ \"\$published_contigs_size\" != \"\$contigs_size\" ] || [ \"\$published_graph_size\" != \"\$graph_size\" ] || [ \"\$published_contigs_sha256\" != \"\$contigs_sha256\" ] || [ \"\$published_graph_sha256\" != \"\$graph_sha256\" ]; then",
        "  rm -f -- \"\$completion_marker\"",
        "  echo \"metaMDBG artifacts changed after completion publication; removed stale completion manifest\" >&2",
        "  exit 1",
        "fi",
    ]
    return join(lines, "\n")
end

function _metamdbg_provenance(
        outputs::NamedTuple,
        input_contract::NamedTuple,
        graph_k::Int,
        ;
        completion::Union{Nothing, NamedTuple} = nothing,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    toolchain = _metamdbg_expected_toolchain()
    workflow = _metamdbg_workflow_contract(outputs, input_contract, graph_k)
    provenance = (;
        contract_signature = input_contract.signature,
        input_platform_attestation =
            input_contract.contract.platform_attestation,
        workflow_signature = workflow.signature,
        graph_k,
        metamdbg_version = toolchain.metamdbg_version,
        environment_name = toolchain.environment_name,
        environment_spec_sha256 = toolchain.environment_spec_sha256,
    )
    completion === nothing && return provenance
    manifest = completion.manifest
    return (;
        provenance...,
        completion_signature = completion.signature,
        package_inventory_sha256 =
            manifest.toolchain.package_inventory_sha256,
        package_count = manifest.toolchain.package_count,
        canonical_outdir_sha256 =
            manifest.workflow.canonical_outdir_sha256,
        contigs_path_sha256 =
            manifest.artifacts.contigs.canonical_path_sha256,
        contigs_size_bytes = manifest.artifacts.contigs.size_bytes,
        contigs_sha256 = manifest.artifacts.contigs.sha256,
        graph_alias_path_sha256 =
            manifest.artifacts.graph.alias_path_sha256,
        graph_source_path_sha256 =
            manifest.artifacts.graph.source_path_sha256,
        graph_size_bytes = manifest.artifacts.graph.size_bytes,
        graph_sha256 = manifest.artifacts.graph.sha256,
    )
end

function _metamdbg_complete_result(
        outputs::NamedTuple,
        artifacts::NamedTuple,
        input_contract::NamedTuple,
        graph_k::Int,
        completion::NamedTuple,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    return (;
        status = :complete,
        outdir = outputs.outdir,
        contigs = artifacts.contigs,
        graph = artifacts.graph,
        contract_marker = outputs.contract_marker,
        completion_marker = outputs.completion_marker,
        provenance = _metamdbg_provenance(
            outputs,
            input_contract,
            graph_k;
            completion,
        ),
    )
end

function _metamdbg_planned_result(
        status::Symbol,
        submission::Any,
        outputs::NamedTuple,
        input_contract::NamedTuple,
        graph_k::Int,
        ;
        submission_reservation::Union{Nothing, NamedTuple} = nothing,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    status in (:planned, :submitted) || throw(ArgumentError(
        "metaMDBG asynchronous status must be :planned or :submitted.",
    ))
    result = (;
        status,
        submission,
        outdir = outputs.outdir,
        expected_artifacts = (;
            contigs = outputs.contigs_gz,
            graph = outputs.graph_alias,
            contract_marker = outputs.contract_marker,
            completion_marker = outputs.completion_marker,
        ),
        provenance = _metamdbg_provenance(outputs, input_contract, graph_k),
    )
    status == :planned && return result
    submission_reservation === nothing && error(
        "Submitted metaMDBG result is missing its reclamation capability.",
    )
    job_id = submission_reservation.job_id
    job_cluster = submission_reservation.job_cluster
    return merge(result, (;
        submission_reservation = (;
            canonical_outdir = submission_reservation.canonical_outdir,
            path = submission_reservation.path,
            workflow_signature = submission_reservation.workflow_signature,
            scheduler_job_name =
                submission_reservation.scheduler_job_name,
            input_contract_signature =
                submission_reservation.input_contract_signature,
            graph_k = submission_reservation.graph_k,
            owner_token = submission_reservation.owner_token,
            job_id,
            job_cluster,
            job_schema_version =
                submission_reservation.job_schema_version,
            submission_state = submission_reservation.submission_state,
            publication_state = hasproperty(
                submission_reservation,
                :publication_state,
            ) ? submission_reservation.publication_state : :published,
            reservation_state = hasproperty(
                submission_reservation,
                :reservation_state,
            ) ? submission_reservation.reservation_state : :queued,
            lifecycle_owner = hasproperty(
                submission_reservation,
                :lifecycle_owner,
            ) ? submission_reservation.lifecycle_owner : :submitter,
        ),
    ))
end

function _metamdbg_existing_artifacts(
        outputs::NamedTuple,
        graph_k::Int,
        ;
        output_root_identity::Union{Nothing, NamedTuple} = nothing,
)::Union{Nothing, NamedTuple}
    _require_positive_metamdbg_graph_k(graph_k)
    if output_root_identity !== nothing
        _require_unchanged_metamdbg_output_root(
            output_root_identity,
            outputs.outdir,
        )
    end
    existing_contigs = _normalize_metamdbg_contigs!(outputs)
    requested_graphs = _metamdbg_graph_candidates(outputs, graph_k)
    all_graphs = _metamdbg_all_graph_candidates(outputs)
    existing_graph = if isempty(all_graphs)
        nothing
    elseif isempty(requested_graphs)
        error(
            "metaMDBG supports exactly one graph_k lifecycle per output " *
            "root. Existing graph artifacts do not match the requested " *
            "graph_k; choose a fresh output root.",
        )
    else
        _normalize_metamdbg_graph!(outputs, graph_k)
    end
    if existing_contigs !== nothing && existing_graph !== nothing
        if output_root_identity !== nothing
            _require_unchanged_metamdbg_output_root(
                output_root_identity,
                outputs.outdir,
            )
            _require_metamdbg_artifact_containment(
                existing_contigs,
                "metaMDBG existing contigs artifact",
                output_root_identity,
            )
            _require_metamdbg_artifact_containment(
                existing_graph,
                "metaMDBG existing graph artifact",
                output_root_identity,
            )
        end
        return (;
            outdir = outputs.outdir,
            contigs = existing_contigs,
            graph = existing_graph,
        )
    elseif existing_contigs === nothing && existing_graph !== nothing
        error(
            "metaMDBG output is inconsistent: graph exists without contigs in " *
            outputs.outdir,
        )
    elseif existing_contigs !== nothing
        error(
            "metaMDBG refuses to adopt partial contracted contigs without a " *
            "realized-stage completion manifest. Remove the partial output " *
            "root and rerun from clean staged inputs.",
        )
    end
    if ispath(outputs.completion_marker) || islink(outputs.completion_marker)
        error(
            "metaMDBG output has a completion manifest without both validated " *
            "artifacts in $(outputs.outdir).",
        )
    end
    return nothing
end

function _require_successful_metamdbg_submission(
        submission::Any,
        ;
        require_sbatch::Bool = false,
)::Any
    if require_sbatch && !(submission isa Mycelia.SubmitResult)
        error(
            "metaMDBG real Slurm submission did not return a verifiable " *
            "SubmitResult.",
        )
    end
    if submission isa Mycelia.SubmitResult
        if !submission.ok
            error_details = isempty(submission.errors) ?
                            "the backend reported no error details" :
                            join(submission.errors, "; ")
            error("metaMDBG submission failed: $(error_details)")
        end
        submission.dry_run && error(
            "metaMDBG real submission unexpectedly returned a dry-run result.",
        )
        submission.backend == :sbatch || error(
            "metaMDBG real nonlocal submission did not use the sbatch " *
            "backend: $(submission.backend).",
        )
        submission.scheduler_acceptance == :accepted || error(
            "metaMDBG real Slurm submission has unverified scheduler " *
            "acceptance state $(repr(submission.scheduler_acceptance)).",
        )
        submission.held || error(
            "metaMDBG real Slurm submission was not verified as held; " *
            "refusing to bind or release an immediate-start job.",
        )
        job_id = submission.job_id
        (job_id isa AbstractString && !isempty(strip(job_id))) || error(
            "metaMDBG sbatch submission returned no job id.",
        )
        normalized_job_id = _normalize_metamdbg_job_id(job_id)
        normalized_job_cluster = _normalize_metamdbg_job_cluster(
            submission.job_cluster,
        )
        String(job_id) == normalized_job_id || error(
            "metaMDBG sbatch submission returned a non-exact job id: " *
            "$(repr(job_id)).",
        )
        if submission.job_cluster !== nothing
            String(submission.job_cluster) == normalized_job_cluster || error(
                "metaMDBG sbatch submission returned a non-exact job cluster: " *
                "$(repr(submission.job_cluster)).",
            )
        end
        if submission.stdout !== nothing
            stdout_job_reference = Mycelia._extract_sbatch_job_reference(
                submission.stdout,
            )
            stdout_job_reference == (;
                job_id = normalized_job_id,
                job_cluster = normalized_job_cluster,
            ) || error(
                "metaMDBG sbatch submission stdout did not contain the exact " *
                "returned job reference.",
            )
        end
    end
    return submission
end

function _metamdbg_submission_recovery_guidance(
        reservation::NamedTuple,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString} = nothing,
)::String
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    normalized_job_cluster = _normalize_metamdbg_job_cluster(job_cluster)
    display_reference = normalized_job_cluster === nothing ?
                        normalized_job_id :
                        "$(normalized_job_id);$(normalized_job_cluster)"
    scontrol_cluster = normalized_job_cluster === nothing ?
                       "" : "-M $(normalized_job_cluster) "
    return "SLURM job $(display_reference) remains held or has an unknown " *
           "release state. Its durable reservation is $(reservation.path) " *
           "and scheduler name is $(reservation.scheduler_job_name). Inspect " *
           "that exact job with `scontrol $(scontrol_cluster)show job " *
           "$(normalized_job_id)`, then " *
           "inspect the reservation with " *
           "`inspect_metamdbg_submission_reservations`. If pending evidence " *
           "remains, first complete confirmed-dead recovery when required and " *
           "resume the exact ID with " *
           "`bind_metamdbg_submission_reservation_job!`; verify a fresh " *
           "inspection has no pending evidence before running `scontrol " *
           "$(scontrol_cluster)release $(normalized_job_id)`. To cancel instead, confirm exact-job " *
           "cancellation, complete the same pending cleanup, and only then " *
           "reclaim the reservation."
end

function _release_metamdbg_submission_job!(
        reservation::NamedTuple,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString},
        release_runner::Function,
)::NamedTuple
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    normalized_job_cluster = _normalize_metamdbg_job_cluster(job_cluster)
    release_result = try
        normalized_job_cluster === nothing ?
        release_runner(normalized_job_id) :
        release_runner(
            normalized_job_id;
            job_cluster = normalized_job_cluster,
        )
    catch caught
        error(
            "metaMDBG failed to release its durably bound scheduler job. " *
            _metamdbg_submission_recovery_guidance(
                reservation,
                normalized_job_id,
                normalized_job_cluster,
            ) * " Cause: $(sprint(showerror, caught))",
        )
    end
    if !(release_result isa AbstractString) ||
       String(release_result) != normalized_job_id
        error(
            "metaMDBG scheduler release was absent or ambiguous for exact job " *
            "$(normalized_job_id); got $(repr(release_result)). " *
            _metamdbg_submission_recovery_guidance(
                reservation,
                normalized_job_id,
                normalized_job_cluster,
            ),
        )
    end
    queued_exists = _output_root_path_entry_exists(reservation.path)
    runtime_exists = _output_root_path_entry_exists(reservation.runtime_path)
    consumed_exists = _output_root_path_entry_exists(reservation.consumed_path)
    state_count = count(identity, Bool[
        queued_exists,
        runtime_exists,
        consumed_exists,
    ])
    state_count <= 1 || error(
        "metaMDBG scheduler release observed multiple private reservation " *
        "states for exact job $(normalized_job_id); refusing an ambiguous " *
        "lifecycle transition.",
    )
    if queued_exists
        current = _metamdbg_submission_reservation_from_path(
            reservation.path,
            reservation.canonical_outdir,
        )
        current.job_id == normalized_job_id || error(
            "metaMDBG queued reservation job identity changed after scheduler " *
            "release.",
        )
        current.job_cluster == normalized_job_cluster || error(
            "metaMDBG queued reservation scheduler cluster changed after " *
            "release.",
        )
        return current
    elseif runtime_exists
        current = _metamdbg_submission_reservation_from_path(
            reservation.runtime_path,
            reservation.canonical_outdir,
        )
        current.job_id == normalized_job_id || error(
            "metaMDBG runtime reservation job identity changed after scheduler " *
            "release.",
        )
        current.job_cluster == normalized_job_cluster || error(
            "metaMDBG runtime reservation scheduler cluster changed after " *
            "release.",
        )
        return current
    elseif consumed_exists
        current = _metamdbg_submission_reservation_from_path(
            reservation.consumed_path,
            reservation.canonical_outdir,
        )
        current.job_id == normalized_job_id || error(
            "metaMDBG consumed reservation job identity changed after " *
            "scheduler release.",
        )
        current.job_cluster == normalized_job_cluster || error(
            "metaMDBG consumed reservation scheduler cluster changed after " *
            "release.",
        )
        return current
    end
    error(
        "metaMDBG scheduler release could not discover a queued, runtime, or " *
        "durably consumed owner record for exact job $(normalized_job_id). " *
        "Refusing to infer :consumed from path absence.",
    )
end

function _cleanup_metamdbg_submission_reservation_after_failure!(
        outputs::NamedTuple,
        reservation::NamedTuple,
)::Nothing
    if !ispath(reservation.path) && !islink(reservation.path)
        return nothing
    end
    _with_metamdbg_output_domain_lock(
        outputs.outdir;
        allowed_same_root_locks = (
            reservation.output_root_reservation_marker,
        ),
    ) do
        _remove_metamdbg_submission_reservation!(reservation)
    end
    return nothing
end

function _bind_metamdbg_submission_job_after_submit!(
        reservation::NamedTuple,
        job_id::AbstractString,
        job_cluster::Union{Nothing, AbstractString} = nothing,
)::NamedTuple
    return _with_metamdbg_output_domain_lock(
        reservation.canonical_outdir;
        allowed_same_root_locks = (
            reservation.output_root_reservation_marker,
        ),
    ) do
        current = _metamdbg_submission_reservation_from_path(
            reservation.path,
            reservation.canonical_outdir,
        )
        current.workflow_signature == reservation.workflow_signature || error(
            "metaMDBG submission reservation workflow changed before its " *
            "held scheduler job could be bound.",
        )
        current.owner_token == reservation.owner_token || error(
            "metaMDBG submission reservation owner changed before its held " *
            "scheduler job could be bound.",
        )
        current.job_id === nothing || error(
            "metaMDBG submission reservation already has a scheduler job id.",
        )
        return _bind_metamdbg_submission_job!(current, job_id, job_cluster)
    end
end

function _metamdbg_dependency_toolchain(
        dependency_checker::Function,
        conda_runner::AbstractString,
)::Any
    if applicable(dependency_checker, conda_runner)
        return dependency_checker(conda_runner)
    end
    return dependency_checker()
end

function _metamdbg_submission_acceptance_is_ambiguous(
        submission::Any,
        attempt_started::Bool,
)::Bool
    attempt_started || return false
    submission isa Mycelia.SubmitResult || return true
    acceptance = submission.scheduler_acceptance
    return acceptance != :not_attempted
end

function _run_metamdbg(;
        hifi_reads::Union{String, Vector{String}, Nothing} = nothing,
        ont_reads::Union{String, Vector{String}, Nothing} = nothing,
        ont_r10_4_plus::Bool = false,
        outdir::String = "metamdbg_output",
        abundance_min::Int = 3,
        threads::Int = get_default_threads(),
        graph_k::Int = 21,
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "metamdbg",
        time_limit::String = "3-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing,
        conda_runner::AbstractString = _conda_runner(),
        dependency_checker::Function = runner ->
            _ensure_metamdbg_installed(conda_runner = runner),
        local_runner::Function = Base.run,
        submission_runner::Function = Mycelia.execute,
        submission_job_binder::Function =
            _bind_metamdbg_submission_job_after_submit!,
        submission_release_runner::Function = Mycelia.release_slurm_job,
        input_digest_function::Function = _metamdbg_sha256,
)::NamedTuple
    _require_positive_metamdbg_graph_k(graph_k)
    selected_input = _metamdbg_selected_input(
        hifi_reads,
        ont_reads,
        ont_r10_4_plus,
    )
    abundance_min > 0 || throw(ArgumentError("abundance_min must be positive."))
    threads > 0 || throw(ArgumentError("threads must be positive."))
    outputs = _metamdbg_output_paths(outdir, graph_k)
    resolved_conda_runner =
        _canonical_metamdbg_conda_runner(conda_runner)
    toolchain = _metamdbg_expected_toolchain()
    resolved_executor = executor === nothing ?
                        Mycelia.LocalExecutor() :
                        Mycelia.resolve_executor(executor)
    supported_executor =
        resolved_executor isa Mycelia.LocalExecutor ||
        resolved_executor isa Mycelia.CollectExecutor ||
        resolved_executor isa Mycelia.DryRunExecutor ||
        resolved_executor isa Mycelia.SlurmExecutor
    supported_executor || throw(ArgumentError(
        "metaMDBG supports only LocalExecutor, CollectExecutor, " *
        "DryRunExecutor, or SlurmExecutor; refusing an unverifiable custom " *
        "nonlocal executor before creating a submission reservation.",
    ))
    _require_no_active_metamdbg_submission_reservation!(outputs.outdir)
    input_contract = _metamdbg_input_contract(
        selected_input,
        abundance_min,
        toolchain;
        digest_function = input_digest_function,
    )
    asm_cmd_args = String[
        "metaMDBG",
        "asm",
        "--out-dir",
        outputs.outdir,
        selected_input.flag,
        selected_input.paths...,
        "--min-abundance",
        string(abundance_min),
        "--threads",
        string(threads),
    ]
    gfa_cmd_args = String[
        "metaMDBG",
        "gfa",
        "--assembly-dir",
        outputs.outdir,
        "--k",
        string(graph_k),
        "--threads",
        string(threads),
    ]
    asm_command = `$(resolved_conda_runner) run --live-stream -n $(METAMDBG_ENV_NAME) $(asm_cmd_args)`
    gfa_command = `$(resolved_conda_runner) run --live-stream -n $(METAMDBG_ENV_NAME) $(gfa_cmd_args)`
    if resolved_executor isa Mycelia.LocalExecutor
        return _with_metamdbg_output_domain_lock(outputs.outdir) do
            _require_no_active_metamdbg_submission_reservation!(
                outputs.outdir,
            )
            _require_current_metamdbg_input_snapshot!(
                selected_input,
                input_contract.input_snapshot,
            )
            has_contract = _prepare_metamdbg_output_root!(
                outputs,
                input_contract,
            )
            output_root_identity =
                _metamdbg_output_root_identity(outputs.outdir)
            existing_artifacts = _metamdbg_existing_artifacts(
                outputs,
                graph_k;
                output_root_identity,
            )
            if existing_artifacts !== nothing
                has_contract || error(
                    "metaMDBG complete output is missing its provenance contract.",
                )
                _require_metamdbg_contract!(outputs, input_contract)
                _require_current_metamdbg_input_contract!(
                    selected_input,
                    abundance_min,
                    input_contract,
                    toolchain,
                    ;
                    digest_function = input_digest_function,
                )
                completion = _require_metamdbg_completion_manifest!(
                    outputs,
                    existing_artifacts,
                    input_contract,
                    graph_k,
                )
                _require_unchanged_metamdbg_output_root(
                    output_root_identity,
                    outputs.outdir,
                )
                return _metamdbg_complete_result(
                    outputs,
                    existing_artifacts,
                    input_contract,
                    graph_k,
                    completion,
                )
            end

            realized_toolchain_before = _require_expected_metamdbg_toolchain(
                _with_metamdbg_output_root_guard(
                    output_root_identity,
                    "checking the pre-execution dependency toolchain",
                ) do
                    _metamdbg_dependency_toolchain(
                        dependency_checker,
                        resolved_conda_runner,
                    )
                end,
            )
            staged = _stage_metamdbg_inputs!(
                selected_input,
                input_contract,
                dirname(outputs.outdir),
            )
            artifacts, realized_toolchain = try
                staged_asm_args = String[
                    "metaMDBG",
                    "asm",
                    "--out-dir",
                    outputs.outdir,
                    staged.selected_input.flag,
                    staged.selected_input.paths...,
                    "--min-abundance",
                    string(abundance_min),
                    "--threads",
                    string(threads),
                ]
                staged_asm_command =
                    `$(resolved_conda_runner) run --live-stream -n $(METAMDBG_ENV_NAME) $(staged_asm_args)`
                _require_current_metamdbg_input_snapshot!(
                    selected_input,
                    input_contract.input_snapshot,
                )
                _with_metamdbg_output_root_guard(
                    output_root_identity,
                    "running metaMDBG assembly",
                ) do
                    local_runner(staged_asm_command)
                end
                existing_contigs = _normalize_metamdbg_contigs!(outputs)
                existing_contigs === nothing && error(
                    "metaMDBG assembly produced no sequence-bearing contigs " *
                    "artifact in $(outputs.outdir).",
                )
                _require_current_metamdbg_input_snapshot!(
                    selected_input,
                    input_contract.input_snapshot,
                )
                _with_metamdbg_output_root_guard(
                    output_root_identity,
                    "running metaMDBG graph generation",
                ) do
                    local_runner(gfa_command)
                end
                _require_unchanged_staged_metamdbg_inputs!(
                    staged,
                    input_contract,
                )
                completed_artifacts = _require_metamdbg_artifacts!(
                    outputs,
                    graph_k;
                    output_root_identity,
                )
                toolchain_after = _require_unchanged_metamdbg_toolchain(
                    realized_toolchain_before,
                    _with_metamdbg_output_root_guard(
                        output_root_identity,
                        "checking the post-execution dependency toolchain",
                    ) do
                        _metamdbg_dependency_toolchain(
                            dependency_checker,
                            resolved_conda_runner,
                        )
                    end,
                )
                completed_artifacts, toolchain_after
            finally
                rm(staged.root; recursive = true, force = true)
            end
            completion = _metamdbg_completion_manifest(
                outputs,
                artifacts,
                input_contract,
                graph_k,
                realized_toolchain,
            )
            _require_current_metamdbg_input_contract!(
                selected_input,
                abundance_min,
                input_contract,
                toolchain,
                ;
                digest_function = input_digest_function,
            )
            _require_unchanged_metamdbg_output_root(
                output_root_identity,
                outputs.outdir,
            )
            _require_metamdbg_artifact_containment(
                artifacts.contigs,
                "metaMDBG contigs artifact before provenance publication",
                output_root_identity,
            )
            _require_metamdbg_artifact_containment(
                artifacts.graph,
                "metaMDBG graph artifact before provenance publication",
                output_root_identity,
            )
            wrote_contract = false
            wrote_completion = false
            contract_identity = nothing
            completion_identity = nothing
            try
                if has_contract
                    _require_metamdbg_contract!(outputs, input_contract)
                else
                    _require_unchanged_metamdbg_output_root(
                        output_root_identity,
                        outputs.outdir,
                    )
                    _write_metamdbg_contract!(outputs, input_contract)
                    wrote_contract = true
                    contract_identity = _metamdbg_regular_file_identity(
                        outputs.contract_marker,
                    )
                end
                _require_metamdbg_contract!(outputs, input_contract)
                _require_unchanged_metamdbg_output_root(
                    output_root_identity,
                    outputs.outdir,
                )
                _write_metamdbg_completion_manifest!(outputs, completion)
                wrote_completion = true
                completion_identity = _metamdbg_regular_file_identity(
                    outputs.completion_marker,
                )
                completion = _require_metamdbg_completion_manifest!(
                    outputs,
                    artifacts,
                    input_contract,
                    graph_k,
                )
                _require_unchanged_metamdbg_output_root(
                    output_root_identity,
                    outputs.outdir,
                )
                return _metamdbg_complete_result(
                    outputs,
                    artifacts,
                    input_contract,
                    graph_k,
                    completion,
                )
            catch primary_error
                if wrote_completion
                    try
                        _remove_exact_metamdbg_durable_file!(
                            outputs.completion_marker,
                            completion_identity,
                        )
                    catch cleanup_error
                        @warn "metaMDBG failed to remove a newly written stale " *
                              "completion manifest while preserving the " *
                              "primary completion failure" outputs primary_error cleanup_error
                    end
                end
                if wrote_contract
                    try
                        _remove_exact_metamdbg_durable_file!(
                            outputs.contract_marker,
                            contract_identity,
                        )
                    catch cleanup_error
                        @warn "metaMDBG failed to remove a newly written stale " *
                              "contract while preserving the primary input " *
                              "contract failure" outputs primary_error cleanup_error
                    end
                end
                Base.rethrow()
            end
        end
    end

    complete_result = _with_metamdbg_output_domain_lock(outputs.outdir) do
        _require_no_active_metamdbg_submission_reservation!(outputs.outdir)
        _require_current_metamdbg_input_snapshot!(
            selected_input,
            input_contract.input_snapshot,
        )
        has_contract = _prepare_metamdbg_output_root!(outputs, input_contract)
        existing_artifacts = _metamdbg_existing_artifacts(outputs, graph_k)
        if existing_artifacts === nothing
            return nothing
        end
        has_contract || error(
            "metaMDBG complete output is missing its provenance contract.",
        )
        _require_metamdbg_contract!(outputs, input_contract)
        _require_current_metamdbg_input_contract!(
            selected_input,
            abundance_min,
            input_contract,
            toolchain,
            ;
            digest_function = input_digest_function,
        )
        completion = _require_metamdbg_completion_manifest!(
            outputs,
            existing_artifacts,
            input_contract,
            graph_k,
        )
        return _metamdbg_complete_result(
            outputs,
            existing_artifacts,
            input_contract,
            graph_k,
            completion,
        )
    end
    complete_result !== nothing && return complete_result

    is_planned =
        resolved_executor isa Mycelia.CollectExecutor ||
        resolved_executor isa Mycelia.DryRunExecutor ||
        (resolved_executor isa Mycelia.SlurmExecutor &&
         resolved_executor.dry_run)
    submission_reservation = _metamdbg_submission_reservation(
        outputs,
        input_contract,
        graph_k;
        job_name,
    )

    script = _metamdbg_executor_script(
        Mycelia.command_string(asm_command),
        Mycelia.command_string(gfa_command),
        outputs,
        graph_k,
        input_contract;
        conda_runner = resolved_conda_runner,
        submission_reservation,
        threads,
    )
    scheduler_executor = if resolved_executor isa Mycelia.SlurmExecutor
        Mycelia.SlurmExecutor(
            dry_run = resolved_executor.dry_run,
            hold = true,
        )
    else
        resolved_executor
    end
    job = Mycelia.build_execution_job(
        cmd = script,
        job_name = submission_reservation.scheduler_job_name,
        site = site,
        time_limit = time_limit,
        cpus_per_task = threads,
        mem_gb = mem_gb,
        partition = partition,
        qos = qos,
        account = account,
        mail_user = mail_user,
    )
    if is_planned
        submission = submission_runner(job, scheduler_executor)
        return _metamdbg_planned_result(
            :planned,
            submission,
            outputs,
            input_contract,
            graph_k,
        )
    end
    if resolved_executor isa Mycelia.SlurmExecutor
        job.site in (:nersc, :lawrencium, :scg) || throw(ArgumentError(
            "metaMDBG real Slurm submission requires a supported batch site.",
        ))
        Mycelia._is_interactive_job(job) && throw(ArgumentError(
            "metaMDBG real Slurm submission does not support interactive jobs.",
        ))
    end

    newly_complete = _with_metamdbg_output_domain_lock(outputs.outdir) do
        _require_no_active_metamdbg_submission_reservation!(outputs.outdir)
        _require_current_metamdbg_input_snapshot!(
            selected_input,
            input_contract.input_snapshot,
        )
        has_contract = _prepare_metamdbg_output_root!(outputs, input_contract)
        existing_artifacts = _metamdbg_existing_artifacts(outputs, graph_k)
        if existing_artifacts !== nothing
            has_contract || error(
                "metaMDBG complete output is missing its provenance contract.",
            )
            _require_metamdbg_contract!(outputs, input_contract)
            _require_current_metamdbg_input_contract!(
                selected_input,
                abundance_min,
                input_contract,
                toolchain,
                ;
                digest_function = input_digest_function,
            )
            completion = _require_metamdbg_completion_manifest!(
                outputs,
                existing_artifacts,
                input_contract,
                graph_k,
            )
            return _metamdbg_complete_result(
                outputs,
                existing_artifacts,
                input_contract,
                graph_k,
                completion,
            )
        end
        _require_current_metamdbg_input_snapshot!(
            selected_input,
            input_contract.input_snapshot,
        )
        if resolved_executor isa Mycelia.SlurmExecutor
            _require_expected_metamdbg_toolchain(
                _metamdbg_dependency_toolchain(
                    dependency_checker,
                    resolved_conda_runner,
                ),
            )
        end
        _create_metamdbg_submission_reservation!(
            submission_reservation,
            outputs.outdir,
        )
        return nothing
    end
    newly_complete !== nothing && return newly_complete

    submission_attempt = nothing
    submission_attempt_started = false
    submission = try
        submission_attempt_started = true
        submission_attempt = submission_runner(job, scheduler_executor)
        _require_successful_metamdbg_submission(
            submission_attempt;
            require_sbatch = resolved_executor isa Mycelia.SlurmExecutor,
        )
    catch primary_error
        if _metamdbg_submission_acceptance_is_ambiguous(
            submission_attempt,
            submission_attempt_started,
        )
            error(
                "metaMDBG received an ambiguous response after attempting a " *
                "held SLURM submission. The unbound reservation was preserved " *
                "at $(submission_reservation.path); scheduler name " *
                "$(submission_reservation.scheduler_job_name). Inspect that " *
                "exact name with `squeue --name " *
                "$(submission_reservation.scheduler_job_name)` before binding " *
                "or reclaiming the reservation. Cause: " *
                sprint(showerror, primary_error),
            )
        end
        try
            _cleanup_metamdbg_submission_reservation_after_failure!(
                outputs,
                submission_reservation,
            )
        catch cleanup_error
            @warn "metaMDBG failed to clean its submission reservation " *
                  "while preserving the primary submission failure" outputs primary_error cleanup_error
        end
        Base.rethrow()
    end
    submission_job_id = String(submission.job_id)
    submission_job_cluster = submission.job_cluster
    bound_submission_reservation = try
        submission_job_cluster === nothing ?
        submission_job_binder(
            submission_reservation,
            submission_job_id,
        ) :
        submission_job_binder(
            submission_reservation,
            submission_job_id,
            submission_job_cluster,
        )
    catch caught
        error(
            "metaMDBG submitted SLURM job $(submission_job_id) held but could " *
            "not durably bind job.json; the job was not released. " *
            _metamdbg_submission_recovery_guidance(
                submission_reservation,
                submission_job_id,
                submission_job_cluster,
            ) * " Cause: $(sprint(showerror, caught))",
        )
    end
    released_submission_reservation = _release_metamdbg_submission_job!(
        bound_submission_reservation,
        submission_job_id,
        submission_job_cluster,
        submission_release_runner,
    )
    return _metamdbg_planned_result(
        :submitted,
        submission,
        outputs,
        input_contract,
        graph_k,
        ;
        submission_reservation = released_submission_reservation,
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Inspect durable metaMDBG submission reservations for an output root.

Each returned record is reconstructed from and verified against its mode-0600
on-disk owner contract. Normal inspection is fully read-only and lockless for
queued, scheduler-owned runtime, reclaiming, and immutable consumed-audit
records. It verifies complete before-and-after owner and pending-path
inventories, double-parsed owner state, and exact owner, queued-marker,
runtime-marker, pending-job, private-lock, and cleanup-sentinel identities. A
concurrent lifecycle transition therefore produces either one coherent stable
state or a fail-loud state-change error; inspection never publishes a PID,
sentinel, or private lock that could interfere with a released runtime. A
shared-marker-paired temporary record is reported with
`publication_state = :provisional`. A durable pending scheduler-job record is
reported as `submission_state = :submission_pending` with its exact
`pending_job_id` and `pending_job_cluster`; if `job.json` is already committed,
the cleanup window is reported as `:submission_commit_cleanup_pending`.
Published, runtime, and consumed records likewise expose `job_cluster`, so a
federated job is always identified by its exact `(job_id, job_cluster)` pair.

After independently proving a pre-submit caller dead, set
`confirm_process_dead = true` to remove its exact same-user submitter locks. That
path requires a pre-existing canonical local dead PID record. Runtime evidence
or private state without that PID fails closed and requires exact scheduler
job-ID cancellation or terminal-state recovery instead.
"""
function inspect_metamdbg_submission_reservations(
        outdir::AbstractString,
        ;
        confirm_process_dead::Bool = false,
)::Vector{NamedTuple}
    return _inspect_metamdbg_submission_reservations(
        outdir;
        confirm_process_dead,
    )
end

function _inspect_metamdbg_submission_reservations(
        outdir::AbstractString,
        ;
        confirm_process_dead::Bool = false,
        post_initial_snapshot_hook::Function =
            (_snapshot::NamedTuple) -> nothing,
        pending_recovery_function::Function =
            _recover_metamdbg_pending_submission_job_records!,
        post_recovery_pid_acquisition_hook::Function =
            (_outdir::AbstractString) -> nothing,
)::Vector{NamedTuple}
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    _require_no_metamdbg_removal_quarantine_evidence!(canonical_outdir)
    all_paths = _metamdbg_submission_reservation_paths(
        canonical_outdir;
        include_consumed = true,
    )
    pending_paths =
        _metamdbg_pending_submission_job_paths(canonical_outdir)
    has_runtime_owner = any(all_paths) do path
        state = _metamdbg_submission_reservation_path_state(
            path,
            canonical_outdir,
        )
        return state == :runtime
    end
    has_runtime_marker =
        _metamdbg_has_runtime_scheduler_evidence(canonical_outdir)
    if confirm_process_dead
        (has_runtime_owner || has_runtime_marker) && error(
            "metaMDBG confirm_process_dead cannot recover scheduler-owned " *
            "runtime state. Confirm the exact job terminal or cancelled and " *
            "use job-ID-bound reservation recovery.",
        )
        _recover_dead_metamdbg_lifecycle_locks!(
            canonical_outdir;
            pending_recovery_function,
            post_recovery_pid_acquisition_hook,
        )
        all_paths = _metamdbg_submission_reservation_paths(
            canonical_outdir;
            include_consumed = true,
        )
        pending_paths =
            _metamdbg_pending_submission_job_paths(canonical_outdir)
    end
    private_lock_path = _metamdbg_output_lock_path(canonical_outdir)
    private_lock_identity = if _output_root_path_entry_exists(private_lock_path)
        _metamdbg_output_lock_identity(private_lock_path)
    else
        nothing
    end
    cleanup_reservation_path =
        _metamdbg_lifecycle_cleanup_reservation_path(canonical_outdir)
    cleanup_reservation_identity = if _output_root_path_entry_exists(
            cleanup_reservation_path,
        )
        _metamdbg_output_lock_identity(cleanup_reservation_path)
    else
        nothing
    end
    function pending_submission_state(
            reservation::NamedTuple,
            pending::Union{Nothing, NamedTuple},
    )::Symbol
        pending === nothing && return reservation.submission_state
        reservation.job_id === nothing && return :submission_pending
        pending.is_complete || error(
            "metaMDBG committed scheduler job record is paired with an " *
            "incomplete pending record.",
        )
        pending.job_id == reservation.job_id || error(
            "metaMDBG committed scheduler job ID conflicts with its pending " *
            "record.",
        )
        pending.job_cluster == reservation.job_cluster || error(
            "metaMDBG committed scheduler cluster conflicts with its pending " *
            "record.",
        )
        pending.job_schema_version == reservation.job_schema_version || error(
            "metaMDBG committed scheduler job schema conflicts with its " *
            "pending record.",
        )
        pending.bytes == collect(codeunits(reservation.job_contents)) || error(
            "metaMDBG committed scheduler job content conflicts with its " *
            "pending record.",
        )
        _metamdbg_pending_submission_job_identity(reservation.job_marker) ==
            pending.identity || error(
            "metaMDBG committed scheduler job marker does not retain its " *
            "pending publication inode.",
        )
        Base.Filesystem.samefile(
            reservation.job_marker,
            pending.path,
        ) || error(
            "metaMDBG committed and pending scheduler job records are not the " *
            "same publication inode.",
        )
        return :submission_commit_cleanup_pending
    end
    phase_one = map(all_paths) do path
        reservation = _metamdbg_submission_reservation_from_path(
            path,
            canonical_outdir;
            allow_provisional = true,
        )
        reservation_identity =
            _metamdbg_submission_reservation_identity(reservation)
        queued_reservation_identity = if _output_root_path_entry_exists(
                reservation.output_root_reservation_marker,
            )
            _metamdbg_shared_reservation_identity(reservation)
        else
            nothing
        end
        runtime_reservation_identity = if _output_root_path_entry_exists(
                reservation.runtime_output_root_reservation_marker,
            )
            _metamdbg_runtime_shared_reservation_identity(reservation)
        else
            nothing
        end
        pending = _metamdbg_pending_submission_job_record_for_reservation(
            reservation,
            pending_paths;
            allow_incomplete = true,
        )
        submission_state = pending_submission_state(reservation, pending)
        return (;
            reservation,
            reservation_identity,
            queued_reservation_identity,
            runtime_reservation_identity,
            pending,
            submission_state,
        )
    end
    accounted_pending_paths = Set(String[
        snapshot.pending.path for snapshot in phase_one
        if snapshot.pending !== nothing
    ])
    accounted_pending_paths == Set(pending_paths) || error(
        "metaMDBG found pending scheduler job evidence without an exact " *
        "durable owner record.",
    )
    post_initial_snapshot_hook((;
        reservations = getproperty.(phase_one, :reservation),
        pending_paths = copy(pending_paths),
    ))
    phase_two_paths = _metamdbg_submission_reservation_paths(
        canonical_outdir;
        include_consumed = true,
    )
    phase_two_paths == all_paths || error(
        "metaMDBG reservation path inventory changed during recovery " *
        "inspection.",
    )
    phase_two_pending_paths =
        _metamdbg_pending_submission_job_paths(canonical_outdir)
    phase_two_pending_paths == pending_paths || error(
        "metaMDBG pending scheduler job path inventory changed during " *
        "recovery inspection.",
    )
    _require_metamdbg_optional_recovery_identity(
        private_lock_path,
        private_lock_identity,
        _metamdbg_output_lock_identity,
        "private lifecycle lock",
    )
    _require_metamdbg_optional_recovery_identity(
        cleanup_reservation_path,
        cleanup_reservation_identity,
        _metamdbg_output_lock_identity,
        "lifecycle cleanup reservation",
    )
    records = map(phase_one) do snapshot
        reservation = snapshot.reservation
        final_reservation = _metamdbg_submission_reservation_from_path(
            reservation.path,
            canonical_outdir;
            allow_provisional = true,
        )
        for field in (
                :canonical_outdir,
                :path,
                :workflow_signature,
                :scheduler_job_name,
                :input_contract_signature,
                :graph_k,
                :owner_token,
                :job_id,
                :job_cluster,
                :job_schema_version,
                :submission_state,
                :publication_state,
                :reservation_state,
                :lifecycle_owner,
            )
            getproperty(final_reservation, field) ==
                getproperty(reservation, field) || error(
                    "metaMDBG reservation $(field) changed during recovery " *
                    "inspection: $(reservation.path).",
                )
        end
        _metamdbg_submission_reservation_identity(final_reservation) ==
            snapshot.reservation_identity || error(
                "metaMDBG reservation owner identity changed during recovery " *
                "inspection: $(reservation.path).",
            )
        _require_metamdbg_optional_recovery_identity(
            reservation.output_root_reservation_marker,
            snapshot.queued_reservation_identity,
            _path -> _metamdbg_shared_reservation_identity(final_reservation),
            "queued shared reservation",
        )
        _require_metamdbg_optional_recovery_identity(
            reservation.runtime_output_root_reservation_marker,
            snapshot.runtime_reservation_identity,
            _path ->
                _metamdbg_runtime_shared_reservation_identity(final_reservation),
            "runtime shared reservation",
        )
        final_pending =
            _metamdbg_pending_submission_job_record_for_reservation(
                final_reservation,
                phase_two_pending_paths;
                allow_incomplete = true,
            )
        (snapshot.pending === nothing) == (final_pending === nothing) || error(
            "metaMDBG pending scheduler job presence changed during recovery " *
            "inspection: $(reservation.path).",
        )
        if snapshot.pending !== nothing
            final_pending =
                _require_unchanged_metamdbg_pending_submission_job_record(
                    final_reservation,
                    snapshot.pending;
                    allow_incomplete = true,
                )
        end
        pending_submission_state(final_reservation, final_pending) ==
            snapshot.submission_state || error(
                "metaMDBG pending scheduler submission state changed during " *
                "recovery inspection: $(reservation.path).",
            )
        pending = snapshot.pending
        return (;
            canonical_outdir = reservation.canonical_outdir,
            path = reservation.path,
            workflow_signature = reservation.workflow_signature,
            scheduler_job_name = reservation.scheduler_job_name,
            input_contract_signature = reservation.input_contract_signature,
            graph_k = reservation.graph_k,
            owner_token = reservation.owner_token,
            job_id = reservation.job_id,
            job_cluster = reservation.job_cluster,
            job_schema_version = reservation.job_schema_version,
            pending_job_id = pending === nothing ? nothing : pending.job_id,
            pending_job_cluster = pending === nothing ?
                                  nothing : pending.job_cluster,
            pending_job_schema_version = pending === nothing ?
                                         nothing : pending.job_schema_version,
            pending_job_path = pending === nothing ? nothing : pending.path,
            pending_job_identity = pending === nothing ?
                                   nothing : pending.identity,
            pending_job_complete = pending === nothing ?
                                   nothing : pending.is_complete,
            submission_state = snapshot.submission_state,
            publication_state = reservation.publication_state,
            reservation_state = reservation.reservation_state,
            lifecycle_owner = reservation.lifecycle_owner,
            reservation_identity = snapshot.reservation_identity,
            shared_reservation_identity =
                snapshot.queued_reservation_identity,
            queued_reservation_identity =
                snapshot.queued_reservation_identity,
            runtime_reservation_identity =
                snapshot.runtime_reservation_identity,
            private_lock_identity,
            cleanup_reservation_identity,
        )
    end
    final_paths = _metamdbg_submission_reservation_paths(
        canonical_outdir;
        include_consumed = true,
    )
    final_paths == all_paths || error(
        "metaMDBG reservation path inventory changed during final recovery " *
        "inspection validation.",
    )
    _metamdbg_pending_submission_job_paths(canonical_outdir) == pending_paths ||
        error(
            "metaMDBG pending scheduler job path inventory changed during " *
            "final recovery inspection validation.",
        )
    _require_metamdbg_optional_recovery_identity(
        private_lock_path,
        private_lock_identity,
        _metamdbg_output_lock_identity,
        "private lifecycle lock",
    )
    _require_metamdbg_optional_recovery_identity(
        cleanup_reservation_path,
        cleanup_reservation_identity,
        _metamdbg_output_lock_identity,
        "lifecycle cleanup reservation",
    )
    return records
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Reclaim a durable metaMDBG submission reservation after explicit confirmation.

This is an explicit cancellation capability, not an expiry mechanism. Pass the
`submission_reservation` metadata returned by `run_metamdbg` with
`status = :submitted`, the exact random `owner_token`, the exact scheduler
`job_id`, and `confirm_cancelled = true` only after the scheduler has confirmed
that the job cannot start. For a federated job, also pass its exact nonempty
`job_cluster`; the cancellation or terminal-state capability is bound to the
full `(job_id, job_cluster)` pair. Leave `job_cluster = nothing` only for an
unscoped job. For a process death before submission, first call
`inspect_metamdbg_submission_reservations`, independently confirm that no job
was submitted, then pass that inspected record, its exact owner token, and
`confirm_not_submitted = true`. The confirmation modes are mutually exclusive.
If the exact submitted job instead reaches a terminal failed or completed
state, pass its exact `job_id` and `confirm_terminal = :failed` or `:completed`
after independently confirming that scheduler state. Runtime and consumed
recovery additionally requires all identities from a fresh inspection. The
reservation is removed only when the complete on-disk owner state still matches
exactly. Pre-submit reclaim transitions run under a canonical local PID record,
cleanup sentinel, and private lifecycle lock. A hard-killed transition remains
discoverable as `:reclaiming` or `:reclaim_release_pending`; while its PID is
live, another reclaim cannot mutate any owner state. After independently
confirming process death, reinspect with `confirm_process_dead = true` before
takeover. Missing or replacement-owner reservations fail loudly and no state
expires automatically.
"""
function reclaim_metamdbg_submission_reservation!(
        metadata::NamedTuple,
        ;
        owner_token::AbstractString,
        job_id::Union{Nothing, AbstractString} = nothing,
        job_cluster::Union{Nothing, AbstractString} = nothing,
        confirm_cancelled::Bool = false,
        confirm_not_submitted::Bool = false,
        confirm_terminal::Union{Nothing, Symbol} = nothing,
)::NamedTuple
    confirmation_count = count(identity, Bool[
        confirm_cancelled,
        confirm_not_submitted,
        confirm_terminal !== nothing,
    ])
    confirmation_count == 1 || throw(ArgumentError(
        "Set exactly one of confirm_cancelled=true, " *
        "confirm_not_submitted=true, or confirm_terminal=:failed/:completed " *
        "after " *
        "independently verifying the corresponding scheduler state.",
    ))
    confirm_terminal === nothing ||
        confirm_terminal in (:failed, :completed) || throw(
        ArgumentError(
            "metaMDBG confirm_terminal accepts only :failed or :completed.",
        ),
    )
    required_fields = (
        :canonical_outdir,
        :path,
        :workflow_signature,
        :input_contract_signature,
        :graph_k,
        :owner_token,
        :job_id,
    )
    all(field -> hasproperty(metadata, field), required_fields) || throw(
        ArgumentError(
            "metaMDBG reservation metadata is incomplete; use the exact " *
            "submission_reservation value returned by run_metamdbg.",
        ),
    )
    publication_state = hasproperty(metadata, :publication_state) ?
                        metadata.publication_state : :published
    publication_state in (:published, :provisional) || throw(ArgumentError(
        "metaMDBG reservation publication_state must be :published or " *
        ":provisional.",
    ))
    submission_state = hasproperty(metadata, :submission_state) ?
                       metadata.submission_state : :reserved
    if submission_state in (
            :submission_pending,
            :submission_commit_cleanup_pending,
        )
        throw(ArgumentError(
            "metaMDBG pending scheduler job evidence may represent an " *
            "accepted job and cannot be reclaimed as unsubmitted or stable " *
            "submitted state. Resume the exact job ID with " *
            "bind_metamdbg_submission_reservation_job! and " *
            "confirm_submitted=true. If a hard-killed binder left lifecycle " *
            "locks, independently confirm it is dead and reinspect with " *
            "confirm_process_dead=true first.",
        ))
    end
    if publication_state == :provisional
        confirm_not_submitted || throw(ArgumentError(
            "A provisional metaMDBG reservation can be reclaimed only after " *
            "confirm_not_submitted=true.",
        ))
    end
    reservation_state = if hasproperty(metadata, :reservation_state)
        metadata.reservation_state
    elseif publication_state == :provisional
        :provisional
    else
        :queued
    end
    runtime_states = (
        :runtime_claiming,
        :runtime_transition,
        :runtime_transition_ambiguous,
        :runtime,
        :runtime_release_pending,
        :consumed,
    )
    reclaim_states = (:reclaiming, :reclaim_release_pending)
    identity_bound_states = (runtime_states..., reclaim_states...)
    reservation_state in (
        :queued,
        :provisional,
        identity_bound_states...,
    ) || throw(
        ArgumentError(
            "metaMDBG reservation metadata has an unsupported durable state: " *
            "$(repr(reservation_state)).",
        ),
    )
    if reservation_state in runtime_states
        confirm_not_submitted && throw(ArgumentError(
            "Scheduler-owned metaMDBG runtime or consumed state cannot be " *
            "reclaimed as not submitted; provide exact job-ID cancellation " *
            "or terminal-state evidence.",
        ))
    end
    if reservation_state in identity_bound_states
        all(
            field -> hasproperty(metadata, field),
            (
                :reservation_identity,
                :queued_reservation_identity,
                :runtime_reservation_identity,
                :private_lock_identity,
                :cleanup_reservation_identity,
            ),
        ) || throw(ArgumentError(
            "Identity-bound metaMDBG recovery requires every exact filesystem " *
            "identity returned by inspection.",
        ))
    end
    if reservation_state in reclaim_states &&
       (metadata.private_lock_identity !== nothing ||
        metadata.cleanup_reservation_identity !== nothing)
        throw(ArgumentError(
            "A metaMDBG reclaim transition still has process-owned lifecycle " *
            "state. Independently confirm that its local process is dead, " *
            "then reinspect with confirm_process_dead=true before takeover.",
        ))
    end
    if confirm_not_submitted
        all(
            field -> hasproperty(metadata, field),
            (:reservation_identity, :shared_reservation_identity),
        ) || throw(ArgumentError(
            "Pre-submit metaMDBG recovery requires the exact private and " *
            "shared filesystem identities returned by inspection.",
        ))
    end
    normalized_owner_token = String(owner_token)
    isempty(normalized_owner_token) && throw(ArgumentError(
        "metaMDBG reservation owner_token must be nonempty.",
    ))
    normalized_owner_token == metadata.owner_token || error(
        "metaMDBG reservation owner token does not match the durable " *
        "reservation capability.",
    )
    metadata_job_cluster = hasproperty(metadata, :job_cluster) ?
                           _normalize_metamdbg_job_cluster(
        metadata.job_cluster,
    ) : nothing
    normalized_job_id = if confirm_cancelled || confirm_terminal !== nothing
        job_id isa AbstractString || throw(ArgumentError(
            "metaMDBG cancelled scheduler job_id must be provided.",
        ))
        normalized = strip(String(job_id))
        isempty(normalized) && throw(ArgumentError(
            "metaMDBG cancelled scheduler job_id must be nonempty.",
        ))
        metadata.job_id isa AbstractString || error(
            "metaMDBG submitted reservation metadata has no scheduler job id.",
        )
        normalized == strip(String(metadata.job_id)) || error(
            "metaMDBG cancelled scheduler job id does not match the submitted " *
            "job.",
        )
        _normalize_metamdbg_job_cluster(job_cluster) == metadata_job_cluster ||
            error(
                "metaMDBG scheduler cluster does not match the submitted job.",
            )
        normalized
    else
        job_id === nothing || throw(ArgumentError(
            "Do not provide a scheduler job_id when confirming that submission " *
            "never occurred.",
        ))
        metadata.job_id === nothing || error(
            "metaMDBG reservation metadata contains a scheduler job id and " *
            "cannot be reclaimed as not submitted.",
        )
        job_cluster === nothing || throw(ArgumentError(
            "Do not provide a scheduler job_cluster when confirming that " *
            "submission never occurred.",
        ))
        metadata_job_cluster === nothing || error(
            "metaMDBG reservation metadata contains a scheduler cluster and " *
            "cannot be reclaimed as not submitted.",
        )
        nothing
    end
    graph_k = metadata.graph_k
    graph_k isa Int || throw(ArgumentError(
        "metaMDBG reservation graph_k metadata must be an Int.",
    ))
    outputs = _metamdbg_output_paths(String(metadata.canonical_outdir), graph_k)
    requested_job_name = hasproperty(metadata, :scheduler_job_name) ?
                         _metamdbg_scheduler_job_name_prefix(
        String(metadata.scheduler_job_name),
        String(metadata.workflow_signature),
    ) : "metamdbg"
    expected_reservation = _metamdbg_submission_reservation(
        outputs,
        (; signature = String(metadata.input_contract_signature)),
        graph_k;
        owner_token = normalized_owner_token,
        job_name = requested_job_name,
    )
    if metadata.job_id isa AbstractString
        expected_reservation = _metamdbg_bound_submission_reservation(
            expected_reservation,
            String(metadata.job_id),
            metadata_job_cluster,
        )
    end
    expected_reservation.workflow_signature == metadata.workflow_signature ||
        error("metaMDBG reservation workflow signature does not match metadata.")
    metadata_path = normpath(abspath(String(metadata.path)))
    if reservation_state == :queued
        expected_reservation.path == metadata_path || error(
            "metaMDBG queued reservation path does not match its recomputed " *
            "workflow path.",
        )
    elseif reservation_state in identity_bound_states
        expected_path = if reservation_state == :consumed
            expected_reservation.consumed_path
        elseif reservation_state == :runtime_claiming
            expected_reservation.path
        elseif reservation_state in reclaim_states
            expected_reservation.reclaiming_path
        else
            expected_reservation.runtime_path
        end
        expected_path == metadata_path || error(
            "metaMDBG runtime reservation path does not match its recomputed " *
            "owner capability.",
        )
    else
        _metamdbg_submission_reservation_path_state(
            metadata_path,
            outputs.outdir,
        ) == :provisional || error(
            "metaMDBG provisional recovery path is not a valid temporary owner " *
            "record: $(metadata_path).",
        )
    end
    ispath(metadata_path) || error(
        "metaMDBG submission reservation is missing or was already consumed: " *
        "$(metadata_path).",
    )
    recover_current = function ()
        ispath(metadata_path) || error(
            "metaMDBG submission reservation is missing or was already consumed: " *
            "$(metadata_path).",
        )
        current = _metamdbg_submission_reservation_from_path(
            metadata_path,
            outputs.outdir;
            allow_provisional = publication_state == :provisional,
        )
        current.publication_state == publication_state || error(
            "metaMDBG reservation publication state changed after inspection.",
        )
        for field in (
                :canonical_outdir,
                :workflow_signature,
                :scheduler_job_name,
                :input_contract_signature,
                :graph_k,
                :owner_token,
                :job_id,
                :job_cluster,
        )
            getproperty(current, field) == getproperty(expected_reservation, field) ||
                error(
                    "metaMDBG submission reservation $(field) changed before " *
                    "explicit recovery.",
                )
        end
        hasproperty(current, :reservation_state) &&
            current.reservation_state == reservation_state || error(
            "metaMDBG durable reservation state changed after inspection.",
        )
        return current
    end
    if reservation_state in runtime_states
        reservation_lock_path =
            _output_root_reservation_lock_path_from_canonical(outputs.outdir)
        stale_age = _OUTPUT_ROOT_RESERVATION_STALE_AGE_SECONDS
        lock_handle = FileWatching.Pidfile.trymkpidlock(
            reservation_lock_path;
            stale_age,
            refresh = stale_age / 2,
        )
        lock_handle === false && error(
            "metaMDBG runtime recovery could not acquire its exact output-root " *
            "recovery lock; another lifecycle may still be active.",
        )
        try
            current = recover_current()
            _remove_metamdbg_runtime_recovery_state!(current, metadata)
        finally
            Base.close(lock_handle)
        end
    elseif reservation_state in reclaim_states
        _with_metamdbg_output_domain_lock(
            outputs.outdir;
            allowed_same_root_locks = (
                expected_reservation.output_root_reservation_marker,
            ),
        ) do
            current = recover_current()
            _metamdbg_submission_reservation_identity(current) ==
                metadata.reservation_identity || error(
                "metaMDBG reclaim owner changed after recovery inspection.",
            )
            _require_metamdbg_optional_recovery_identity(
                current.output_root_reservation_marker,
                metadata.queued_reservation_identity,
                path -> begin
                    marker =
                        _require_metamdbg_output_root_reservation_marker!(
                            current,
                        )
                    marker == path || error(
                        "metaMDBG queued reclaim marker path changed.",
                    )
                    marker_status = stat(marker)
                    return (;
                        device = marker_status.device,
                        inode = marker_status.inode,
                    )
                end,
                "queued reclaim reservation",
            )
            _require_metamdbg_optional_recovery_identity(
                current.runtime_output_root_reservation_marker,
                metadata.runtime_reservation_identity,
                _path ->
                    _metamdbg_runtime_shared_reservation_identity(current),
                "runtime reclaim reservation",
            )
            _remove_metamdbg_submission_reservation!(current)
        end
    else
        _with_metamdbg_output_domain_lock(
            outputs.outdir;
            allowed_same_root_locks = (
                expected_reservation.output_root_reservation_marker,
            ),
        ) do
            current = recover_current()
            has_reservation_identity =
                hasproperty(metadata, :reservation_identity)
            has_shared_identity =
                hasproperty(metadata, :shared_reservation_identity)
            has_reservation_identity == has_shared_identity || throw(
                ArgumentError(
                    "metaMDBG recovery metadata must provide both private " *
                    "and shared filesystem identities together.",
                ),
            )
            if has_reservation_identity
                _require_unchanged_metamdbg_recovery_identities(
                    current,
                    metadata.reservation_identity,
                    metadata.shared_reservation_identity,
                )
            end
            _remove_metamdbg_submission_reservation!(current)
        end
    end
    return (;
        status = :reclaimed,
        job_id = normalized_job_id,
        job_cluster = metadata_job_cluster,
        recovery_reason = if confirm_cancelled
            :cancelled
        elseif confirm_terminal !== nothing
            confirm_terminal
        else
            :not_submitted
        end,
        path = metadata_path,
        publication_state,
        reservation_state,
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run metaMDBG assembler for metagenomic long-read assembly.

# Arguments
- `hifi_reads::Union{String,Vector{String},Nothing}`: Existing, nonempty HiFi
  (PacBio) read file(s) (default: nothing).
- `ont_reads::Union{String,Vector{String},Nothing}`: Existing, nonempty ONT
  (Nanopore) read file(s) (default: nothing).
- `ont_r10_4_plus::Bool`: Required explicit attestation for `ont_reads` that
  every input was generated with Nanopore R10.4 or later chemistry. Generic,
  R9, and unknown ONT inputs are rejected (default: false). It must remain
  false for HiFi input.
- `outdir::String`: Output directory path (default: "metamdbg_output").
- `abundance_min::Int`: Minimum abundance threshold (default: 3).
- `threads::Int`: Number of threads to use (default: get_default_threads()).
- `graph_k::Int`: Strictly positive graph resolution requested from
  `metaMDBG gfa` (default: 21).

Exactly one input technology is required. metaMDBG v1.4 rejects simultaneous
`--in-hifi` and `--in-ont`; mixed HiFi-plus-ONT assembly is therefore excluded
from this wrapper rather than advertised as a false combined-input contract.

# Returns
Synchronous execution and complete reuse return `status = :complete`, `outdir`,
the sequence-bearing `contigs.fasta.gz`, a stable `assemblyGraph_k<k>.gfa`
alias, the input-contract marker, an atomic completion manifest, and exact
resolved toolchain-inventory provenance. Nonlocal execution
returns `status = :planned` for collected/dry-run jobs or `:submitted` for a
real submission, together with the backend `submission` result and
`expected_artifacts`; submitted results also include the exact random-owner
`submission_reservation` cancellation capability. Planned paths are not
reported as completed artifacts.
The graph alias resolves to metaMDBG's validated dynamic
`assemblyGraph_k<k>_<length>bps.gfa` artifact.

# Details
- Uses metaMDBG's minimizer-space de Bruijn graphs for metagenomic assembly.
- Accepts both upstream `contigs.fasta.gz` and legacy plain `contigs.fasta`; a
  plain artifact is retained and normalized to `contigs.fasta.gz`.
- Reuses a complete contigs-plus-requested-graph result without provisioning or
  rerunning metaMDBG. Partial contigs or graphs are never resumed because they
  lack realized-stage provenance. Complete reuse requires an exact durable
  contract marker for the normalized input paths, technology flag, input file
  sizes and SHA-256 content digests, `abundance_min`, metaMDBG 1.4, the
  spec-addressed environment name, and the bundled environment-spec checksum.
  Modification times exist only in the invocation snapshot and are never
  serialized into the durable schema-v5 contract. For ONT input, that contract
  and its signature bind the required Nanopore R10.4-or-later attestation. Tool
  execution uses private,
  mode-0400 staged input copies verified against that contract, so transient
  source mutation cannot affect assembled bytes. Local and runtime lifecycles
  compare the complete normalized Conda inventory before and after all tool
  commands. Nonlocal callers hash once, while runtime jobs hash source inputs at
  queued start and immediately before finalization. Any
  nonempty uncontracted output root fails before provisioning. Complete reuse
  additionally verifies the atomic completion manifest against `graph_k`, the
  workflow signature, normalized resolved package-inventory digest, and current
  canonical artifact identities, sizes, and SHA-256 digests.
  One output root owns exactly one `graph_k`; request another graph in a fresh
  output root.
- Every nonlocal job embeds an adjacent, contract-addressed submission
  reservation. Real submissions create it atomically before submission and the
  runtime consumes it under the output lock; collected and dry-run jobs persist
  no reservation and therefore fail closed if executed directly. Active
  federated submissions preserve the exact `(job_id, job_cluster)` identity in
  pending, committed, runtime, inspection, release, and reclaim state; generated
  runtimes require matching `SLURM_JOB_ID` and `SLURM_CLUSTER_NAME`, and held
  jobs are released with cluster-scoped `scontrol -M`.
  Active
  reservations block competing execution before input hashing. Runtime jobs
  validate their exact reservation, publish and fsync their exact runtime
  marker, and check output-domain exclusivity before bounded private-lock
  retries. The durable marker makes every subsequent hard-kill window
  scheduler-inspectable. After acquiring the private lock and rechecking the
  domain, they atomically rename the owner to a same-parent, capability-keyed
  runtime record.
  Reservation creation makes the complete temporary owner record and its parent
  entry durable, then publishes the shared marker before its private atomic
  rename. An incomplete temporary remnant is nonblocking only when no durable
  same-root marker could pair with it; ambiguous corruption fails loudly. A
  temporary owner record paired to its exact shared marker is an active,
  inspectable `publication_state = :provisional` recovery capability. This
  covers hard termination before rename and rename durability. Runtime
  transition records remain inspectable across hard termination. Successful
  cleanup durably renames the owner record to a persistent, nonblocking
  `consumed.<capability>` audit state only after releasing the queued and runtime
  shared markers.
  Reservations never auto-expire:
  after scheduler-confirmed cancellation, reclaim one explicitly with
  `reclaim_metamdbg_submission_reservation!`, the returned owner token and job
  id plus cluster when federated, and `confirm_cancelled = true`. If a caller dies before submission,
  recover the durable capability with
  `inspect_metamdbg_submission_reservations`; after proving the caller process
  dead, `confirm_process_dead = true` removes its exact identity-checked private
  lock, cleanup sentinel, and canonical local PID file. Live, remote, malformed,
  empty, absent, and replacement PID files fail closed; the flag refuses a
  no-op when no pre-existing PID ownership record exists. Pre-submit reclaim and recovery
  job binding require the exact private and shared identities returned by
  inspection and refuse same-content replacements. Pre-submit cleanup first
  renames its owner to the same-parent, capability-keyed `reclaiming` state and
  fsyncs that transition before releasing the queued marker. It holds the local
  PID, sentinel, and private lock throughout; a hard-killed transition therefore
  remains inspectable, blocks a second live reclaim without inode changes, and
  requires confirmed-dead inspection before takeover. Reclaim only after
  explicit independent confirmation that no job was submitted.
  A submitted runtime or consumed record confirmed terminal-failed or
  terminal-completed can be reclaimed from freshly inspected, identity-bound
  metadata with its exact job id and cluster and `confirm_terminal = :failed` or
  `:completed`.
- Installs metaMDBG exactly 1.4 from a checksum-verified, spec-hash-addressed
  environment and records a normalized digest over every resolved package's
  name, version, build, and channel before execution.
- Local and executor lifecycles validate nonempty unique FASTA identifiers,
  gap-free IUPAC DNA FASTA, and strict H/S/L/P GFA1 structure, shared segment/
  path naming, and typed A/i/f/Z/J/H/B optional tags. GFA validation is two-pass
  and retains segment/path identifiers rather than every edge or path reference.
  Completion is stamped only after semantic validation, exact inventory
  capture, and a final input-content digest pass.
"""
function run_metamdbg(;
        hifi_reads::Union{String, Vector{String}, Nothing} = nothing,
        ont_reads::Union{String, Vector{String}, Nothing} = nothing,
        ont_r10_4_plus::Bool = false,
        outdir::String = "metamdbg_output",
        abundance_min::Int = 3,
        threads::Int = get_default_threads(),
        graph_k::Int = 21,
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "metamdbg",
        time_limit::String = "3-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing,
        conda_runner::AbstractString = _conda_runner(),
)::NamedTuple
    return _run_metamdbg(;
        hifi_reads,
        ont_reads,
        ont_r10_4_plus,
        outdir,
        abundance_min,
        threads,
        graph_k,
        executor,
        site,
        job_name,
        time_limit,
        partition,
        account,
        mem_gb,
        qos,
        mail_user,
        conda_runner,
    )
end

function _strainy_output_paths(outdir::String, stage::Symbol)::NamedTuple
    stage in (:phase, :transform, :e2e) || throw(
        ArgumentError("stage must be :phase, :transform, or :e2e"))
    has_phase_outputs = stage in (:phase, :e2e)
    has_transform_outputs = stage in (:transform, :e2e)
    alignment_phased_bam = joinpath(outdir, "alignment_phased.bam")
    consensus_cache = joinpath(outdir, "intermediate", "consensus_dict.pkl")
    strain_unitigs_gfa = joinpath(outdir, "strain_unitigs.gfa")
    strain_contigs_gfa = joinpath(outdir, "strain_contigs.gfa")
    strain_variants_vcf = joinpath(outdir, "strain_variants.vcf")
    phased_unitig_info = joinpath(outdir, "phased_unitig_info_table.csv")
    multiplicity_stats = joinpath(outdir, "multiplicity_stats.txt")
    experimental_segment_candidates_fasta =
        joinpath(outdir, "experimental_segment_candidates.fasta")
    return (;
        outdir,
        stage,
        alignment_phased_bam =
            has_phase_outputs ? alignment_phased_bam : nothing,
        alignment_phased_bai =
            has_phase_outputs ? alignment_phased_bam * ".bai" : nothing,
        # Transform consumes the cache created by a prior phase run.
        consensus_cache,
        strain_unitigs_gfa = has_transform_outputs ? strain_unitigs_gfa : nothing,
        strain_contigs_gfa = has_transform_outputs ? strain_contigs_gfa : nothing,
        strain_variants_vcf = has_transform_outputs ? strain_variants_vcf : nothing,
        phased_unitig_info = has_transform_outputs ? phased_unitig_info : nothing,
        multiplicity_stats = has_transform_outputs ? multiplicity_stats : nothing,
        experimental_segment_candidates_fasta =
            has_transform_outputs ? experimental_segment_candidates_fasta : nothing,
        # Deprecated compatibility alias. This FASTA contains experimental
        # Strainy S-record segments, not complete strain assemblies.
        strain_assemblies =
            has_transform_outputs ? experimental_segment_candidates_fasta : nothing,
    )
end

function _strainy_required_outputs(outputs::NamedTuple)::Vector{String}
    required = String[]
    if outputs.stage in (:phase, :e2e)
        append!(required, [
            outputs.alignment_phased_bam,
            outputs.alignment_phased_bai,
            outputs.consensus_cache,
        ])
    end
    if outputs.stage === :transform
        push!(required, outputs.consensus_cache)
    end
    if outputs.stage in (:transform, :e2e)
        append!(required, [
            outputs.strain_unitigs_gfa,
            outputs.strain_contigs_gfa,
            outputs.strain_variants_vcf,
            outputs.phased_unitig_info,
            outputs.multiplicity_stats,
        ])
    end
    return required
end

function _write_strainy_segment_candidates_fasta(
        gfa_path::String,
        fasta_path::String
)::String
    isfile(gfa_path) || error("Strainy GFA does not exist: $(gfa_path)")
    mkpath(dirname(fasta_path))
    mktemp(dirname(fasta_path)) do temporary_path, output
        identifiers = Set{String}()
        n_segments = 0
        open(gfa_path, "r") do input
            for (line_number, line) in enumerate(eachline(input))
                fields = split(line, '\t'; keepempty = true)
                isempty(fields) && continue
                fields[1] == "S" || continue
                length(fields) >= 3 || error(
                    "malformed Strainy GFA S record at line $(line_number)")
                identifier = String(fields[2])
                isempty(identifier) && error(
                    "empty Strainy GFA segment identifier at line $(line_number)")
                identifier in identifiers && error(
                    "duplicate Strainy GFA segment identifier: $(identifier)")
                push!(identifiers, identifier)
                sequence = String(fields[3])
                sequence == "*" && continue
                isempty(sequence) && error(
                    "empty Strainy GFA segment sequence at line $(line_number)")
                try
                    BioSequences.LongDNA{4}(sequence)
                catch exception
                    error(
                        "invalid DNA in Strainy GFA segment $(identifier): " *
                        sprint(showerror, exception),
                    )
                end
                println(
                    output,
                    ">experimental_segment_candidate|id=$(identifier)|" *
                    "source=strainy_S_record",
                )
                println(output, sequence)
                n_segments += 1
            end
        end
        n_segments > 0 || error(
            "Strainy GFA contains no sequence-bearing S records: $(gfa_path)")
        close(output)
        mv(temporary_path, fasta_path; force = true)
    end
    return fasta_path
end

function _strainy_executor_script(
        strainy_cmd::String,
        outputs::NamedTuple
)::String
    required = _strainy_required_outputs(outputs)
    tests = ["[ -f $(Base.shell_escape(path)) ]" for path in required]
    lines = [
        "set -euo pipefail",
        "mkdir -p $(Base.shell_escape(outputs.outdir))",
        "if ! { $(join(tests, " && ")); }; then",
        "  $(strainy_cmd)",
        "fi",
    ]
    append!(lines, ["test -f $(Base.shell_escape(path))" for path in required])
    if outputs.stage in (:transform, :e2e)
        gfa = Base.shell_escape(outputs.strain_contigs_gfa)
        fasta = Base.shell_escape(outputs.experimental_segment_candidates_fasta)
        temporary = Base.shell_escape(
            outputs.experimental_segment_candidates_fasta * ".tmp")
        awk_program = "'BEGIN { n = 0 } " *
                      "\$1 == \"S\" { " *
                      "if (NF < 3 || \$2 == \"\" || seen[\$2]++) exit 2; " *
                      "if (\$3 == \"*\") next; " *
                      "if (\$3 !~ /^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+\$/) " *
                      "exit 3; " *
                      "print \">experimental_segment_candidate|id=\" \$2 " *
                      "\"|source=strainy_S_record\"; print \$3; n += 1 } " *
                      "END { if (n == 0) exit 1 }'"
        push!(lines, "awk -F '\\t' $(awk_program) $(gfa) > $(temporary)")
        push!(lines, "mv $(temporary) $(fasta)")
        push!(lines, "test -s $(fasta)")
    end
    return join(lines, "\n")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Run a Strainy phase, transform, or end-to-end stage from long reads.

The returned named tuple has stable fields for every stage. Fields unavailable
for the selected stage are `nothing`. `:phase` returns the phased BAM, BAM
index, and consensus cache. `:transform` requires a prior phase consensus cache
and returns Strainy's GFA, VCF, coverage, and multiplicity outputs plus a
deterministic FASTA conversion of sequence-bearing GFA `S` records. `:e2e`
returns both sets.

Stable fields are `outdir`, `stage`, `alignment_phased_bam`,
`alignment_phased_bai`, `consensus_cache`, `strain_unitigs_gfa`,
`strain_contigs_gfa`, `strain_variants_vcf`, `phased_unitig_info`,
`multiplicity_stats`, and `experimental_segment_candidates_fasta`.

The converted FASTA contains **experimental segment candidates**, not complete
strain paths, consensus genomes, or inferred haplotypes. `strain_assemblies` is
a deprecated compatibility alias for that FASTA.
"""
function run_strainy(
        ; gfa_ref::Union{Nothing, String} = nothing,
        fasta_ref::Union{Nothing, String} = nothing,
        fastq::String,
        outdir::String = "strainy_output",
        read_mode::Symbol = :nano,
        stage::Symbol = :e2e,
        min_unitig_coverage::Int = 3,
        unitig_split_length::Int = 0,
        threads::Int = get_default_threads(),
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "strainy",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing
)::NamedTuple
    (gfa_ref === nothing) == (fasta_ref === nothing) &&
        throw(ArgumentError("exactly one of gfa_ref or fasta_ref must be provided"))
    read_mode in (:nano, :hifi) ||
        throw(ArgumentError("read_mode must be :nano or :hifi, got :$(read_mode)"))
    stage in (:phase, :transform, :e2e) ||
        throw(ArgumentError("stage must be :phase, :transform, or :e2e, got :$(stage)"))
    min_unitig_coverage > 0 ||
        throw(ArgumentError("min_unitig_coverage must be positive"))
    unitig_split_length >= 0 ||
        throw(ArgumentError("unitig_split_length must be nonnegative"))
    threads > 0 || throw(ArgumentError("threads must be positive"))

    outputs = _strainy_output_paths(outdir, stage)
    if stage === :transform && !isfile(outputs.consensus_cache)
        error(
            "Strainy stage=:transform requires a prior phase consensus cache: " *
            outputs.consensus_cache,
        )
    end
    Mycelia.add_bioconda_env("strainy")
    mkpath(outdir)

    required = _strainy_required_outputs(outputs)
    ref_flag = gfa_ref === nothing ? "--fasta_ref" : "--gfa_ref"
    ref_path = something(gfa_ref, fasta_ref)
    cmd_args = ["strainy", ref_flag, ref_path, "--fastq", fastq,
        "--output", outdir, "--mode", String(read_mode), "--stage", String(stage),
        "--threads", string(threads), "--min-unitig-coverage",
        string(min_unitig_coverage), "--unitig-split-length", string(unitig_split_length)]
    if executor !== nothing
        strainy_cmd = Mycelia.command_string(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy $(cmd_args)`
        )
        script = _strainy_executor_script(strainy_cmd, outputs)
        job = Mycelia.build_execution_job(
            cmd = script,
            job_name = job_name,
            site = site,
            time_limit = time_limit,
            cpus_per_task = threads,
            mem_gb = mem_gb,
            partition = partition,
            qos = qos,
            account = account,
            mail_user = mail_user
        )
        Mycelia.execute(job, Mycelia.resolve_executor(executor))
        return outputs
    end

    if !all(isfile, required)
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n strainy $(cmd_args)`)
    end

    missing_outputs = filter(path -> !isfile(path), required)
    isempty(missing_outputs) ||
        error("Strainy did not produce required outputs: $(join(missing_outputs, ", "))")
    if stage in (:transform, :e2e)
        _write_strainy_segment_candidates_fasta(
            outputs.strain_contigs_gfa,
            outputs.experimental_segment_candidates_fasta,
        )
    end
    return outputs
end

function _legacy_strainy_stage(mode::String)::Symbol
    return mode == "phase" ? :e2e : Symbol(mode)
end

function run_strainy(
        assembly_file::String,
        long_reads::String;
        outdir::String = "strainy_output",
        mode::String = "phase",
        read_mode::Symbol = :nano,
        kwargs...
)::NamedTuple
    Base.depwarn(
        "run_strainy(assembly_file, reads; mode=...) is deprecated; use " *
        "run_strainy(; fasta_ref=..., fastq=..., read_mode=..., stage=...)",
        :run_strainy)
    legacy_stage = _legacy_strainy_stage(mode)
    return run_strainy(;
        fasta_ref = assembly_file,
        fastq = long_reads,
        outdir = outdir,
        read_mode = read_mode,
        stage = legacy_stage,
        kwargs...)
end
