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
const _UNICYCLER_CONTRACT_SCHEMA = "mycelia-unicycler-run-contract-v1"

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

function _unicycler_fasta_sequences(
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    records = Dict{String, String}()
    reader = try
        Mycelia.open_fastx(path)
    catch caught
        caught isa InterruptException && rethrow()
        error(
            "Unicycler $(label) is not valid FASTA: $(repr(path)). Cause: " *
            sprint(showerror, caught),
        )
    end
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

function _unicycler_gfa_segment_sequences(
        path::AbstractString,
        label::AbstractString,
)::Dict{String, String}
    _require_valid_metamdbg_gfa(path, "Unicycler $(label)")
    segments = Dict{String, String}()
    for (line_number, line) in enumerate(eachline(path))
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
    return (; assembly, graph)
end

function _unicycler_artifact_fingerprint(path::AbstractString)::Dict{String, Any}
    normalized_path = normpath(abspath(String(path)))
    return Dict{String, Any}(
        "canonical_path" => realpath(normalized_path),
        "size_bytes" => filesize(normalized_path),
        "sha256" => _unicycler_sha256(normalized_path),
    )
end

function _unicycler_run_contract(
        input_fingerprints::AbstractDict,
        threads::Int,
        spades_options::Union{Nothing, String},
        kmers::Union{Nothing, String},
        environment_prefix::AbstractString,
        toolchain::AbstractDict,
        assembly::AbstractString,
        graph::AbstractString,
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
            "assembly" => _unicycler_artifact_fingerprint(assembly),
            "graph" => _unicycler_artifact_fingerprint(graph),
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
        close(temporary_io)
        mv(temporary_path, marker)
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
    contract_marker = joinpath(planned_outdir, _UNICYCLER_CONTRACT_FILENAME)
    command = _unicycler_command(;
        conda_runner = resolved_conda_runner,
        environment_prefix = resolved_environment_prefix,
        short_1,
        short_2,
        long_reads,
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
                validated_artifacts = _require_valid_unicycler_artifacts(
                    assembly,
                    graph,
                    reserved_outdir,
                )
                (ispath(contract_marker) || islink(contract_marker)) || error(
                    "Refusing unbound legacy Unicycler output reuse: missing " *
                    "durable run-contract marker $(repr(contract_marker)).",
                )
                assembly = validated_artifacts.assembly
                graph = validated_artifacts.graph
                input_fingerprints = _unicycler_input_fingerprints(
                    short_1,
                    short_2,
                    long_reads;
                    fingerprinter = input_fingerprinter,
                )
                expected_contract = _unicycler_run_contract(
                    input_fingerprints,
                    threads,
                    spades_options,
                    kmers,
                    resolved_environment_prefix,
                    toolchain_before,
                    assembly,
                    graph,
                )
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
                short_1,
                short_2,
                long_reads,
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
            input_fingerprints = _unicycler_input_fingerprints(
                short_1,
                short_2,
                long_reads;
                fingerprinter = input_fingerprinter,
            )
            command_runner(command)
            _require_unchanged_unicycler_output_plan(
                normalized_outdir,
                reserved_outdir,
            )
            validated_artifacts = _require_valid_unicycler_artifacts(
                assembly,
                graph,
                reserved_outdir,
            )
            assembly = validated_artifacts.assembly
            graph = validated_artifacts.graph
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
            validated_artifacts = _require_valid_unicycler_artifacts(
                assembly,
                graph,
                reserved_outdir,
            )
            assembly = validated_artifacts.assembly
            graph = validated_artifacts.graph
            observed_inputs = _unicycler_input_fingerprints(
                short_1,
                short_2,
                long_reads;
                fingerprinter = input_fingerprinter,
            )
            observed_inputs == input_fingerprints || error(
                "Unicycler input content changed while the assembler ran; " *
                "refusing a stale run contract.",
            )
            contract = _unicycler_run_contract(
                input_fingerprints,
                threads,
                spades_options,
                kmers,
                resolved_environment_prefix,
                toolchain_after,
                assembly,
                graph,
            )
            contract_marker = _write_unicycler_run_contract(
                reserved_outdir,
                contract,
            )
            _require_unchanged_unicycler_output_plan(
                normalized_outdir,
                reserved_outdir,
            )
            _require_matching_unicycler_run_contract(
                reserved_outdir,
                contract,
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
)::NamedTuple
    return _run_unicycler_with_contract(;
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
)::NamedTuple
    has_hifi_reads = hifi_reads !== nothing
    has_ont_reads = ont_reads !== nothing
    if has_hifi_reads == has_ont_reads
        throw(ArgumentError(
            "metaMDBG requires exactly one input technology: provide " *
            "either hifi_reads or ont_reads, but not both.",
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
    return (; flag, paths = normalized_paths)
end

const METAMDBG_VERSION = "1.4"
const METAMDBG_ENVIRONMENT_SPEC_SHA256 =
    "3b51b282e8aa768da12e253af01dee43fa96a320baee96755ebb3c123723ff87"
const METAMDBG_ENV_NAME =
    "metamdbg-$(METAMDBG_VERSION)-$(first(METAMDBG_ENVIRONMENT_SPEC_SHA256, 16))"
const _METAMDBG_INSTALL_LOCK_STALE_SECONDS = 600
const _METAMDBG_INSTALL_LOCK_REFRESH_SECONDS = 60
const _METAMDBG_CONTRACT_SCHEMA_VERSION = 4
const _METAMDBG_CONTRACT_FILENAME = "mycelia_metamdbg_contract.json"
const _METAMDBG_COMPLETION_SCHEMA_VERSION = 1
const _METAMDBG_COMPLETION_FILENAME = "mycelia_metamdbg_completion.json"
const _METAMDBG_SUBMISSION_RESERVATION_SCHEMA_VERSION = 3
const _METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME = "contract.json"
const _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION = 1
const _METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME = "job.json"
const _METAMDBG_OUTPUT_ROOT_RESERVATION_SCHEMA_VERSION = 1
const _METAMDBG_OUTPUT_LOCK_RETRY_ATTEMPTS = 300
const _METAMDBG_OUTPUT_LOCK_RETRY_DELAY_SECONDS = 1.0
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
    return joinpath(
        dirname(normalized_outdir),
        ".$(basename(normalized_outdir)).mycelia-metamdbg.lock",
    )
end

function _with_metamdbg_output_lock(
        action::Function,
        outdir::AbstractString,
        ;
        lock_remover::Function = path -> rm(path),
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
    result = try
        _metamdbg_canonical_output_path(canonical_outdir) == canonical_outdir ||
            error("metaMDBG output path changed while acquiring its lock.")
        action()
    catch primary_error
        try
            lock_remover(lock_path)
        catch cleanup_error
            @warn "metaMDBG failed to release its output lock while preserving " *
                  "the primary workflow failure" lock_path primary_error cleanup_error
        end
        Base.rethrow()
    end
    lock_remover(lock_path)
    return result
end

function _with_metamdbg_output_domain_lock(
        action::Function,
        outdir::AbstractString,
)::Any
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    reservation_lock_path =
        _output_root_reservation_lock_path_from_canonical(canonical_outdir)
    mkpath(dirname(reservation_lock_path))
    stale_age = _OUTPUT_ROOT_RESERVATION_STALE_AGE_SECONDS
    lock_handle = FileWatching.Pidfile.trymkpidlock(
        reservation_lock_path;
        stale_age,
        refresh = stale_age / 2,
    )
    lock_handle === false && throw(ArgumentError(
        "metaMDBG output is already reserved by another output-root " *
        "workflow: $(reservation_lock_path)",
    ))
    try
        _require_exclusive_output_root_reservation(
            canonical_outdir,
            reservation_lock_path;
            subject = "metaMDBG outdir",
            stale_age,
        )
        return _with_metamdbg_output_lock(action, canonical_outdir)
    finally
        Base.close(lock_handle)
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
                "$(lpad(input_index, 6, '0'))--$(basename(input.path))",
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
    return ".$(basename(normalized_outdir)).mycelia-metamdbg-submission."
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
        job_contents = nothing,
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

function _metamdbg_submission_job_parts(
        reservation::NamedTuple,
)::Tuple{String, String}
    prefix_object = JSON.json((;
        schema_version = _METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION,
        workflow_signature = reservation.workflow_signature,
        owner_token = reservation.owner_token,
    ))
    endswith(prefix_object, "}") || error(
        "metaMDBG failed to construct its submission job record prefix.",
    )
    prefix = chop(prefix_object; tail = 1) * ",\"job_id\":\""
    return prefix, "\"}\n"
end

function _metamdbg_submission_job_contents(
        reservation::NamedTuple,
        job_id::AbstractString,
)::String
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    prefix, suffix = _metamdbg_submission_job_parts(reservation)
    return prefix * normalized_job_id * suffix
end

function _metamdbg_bound_submission_reservation(
        reservation::NamedTuple,
        job_id::AbstractString,
)::NamedTuple
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    return merge(reservation, (;
        job_id = normalized_job_id,
        job_contents =
            _metamdbg_submission_job_contents(reservation, normalized_job_id),
        submission_state = :submitted,
    ))
end

function _metamdbg_submission_reservation_paths(
        outdir::AbstractString,
)::Vector{String}
    normalized_outdir = _metamdbg_canonical_output_path(outdir)
    parent = dirname(normalized_outdir)
    isdir(parent) || return String[]
    prefix = _metamdbg_submission_reservation_prefix(normalized_outdir)
    matching_paths = sort!(filter(
        path -> startswith(basename(path), prefix),
        readdir(parent; join = true),
    ))
    active_paths = String[]
    for path in matching_paths
        filename = basename(path)
        suffix = filename[(ncodeunits(prefix) + 1):end]
        if occursin(r"^[0-9a-f]{64}$", suffix)
            push!(active_paths, path)
        elseif startswith(suffix, "tmp.") || startswith(suffix, "consumed.")
            continue
        else
            error(
                "metaMDBG found a malformed submission reservation entry: " *
                "$(path). Remove it only after confirming no job owns it.",
            )
        end
    end
    return active_paths
end

function _fsync_metamdbg_descriptor(
        descriptor::Union{Integer, Base.RawFD},
        label::AbstractString,
)::Nothing
    raw_descriptor = descriptor isa Base.RawFD ?
                     reinterpret(Cint, descriptor) : Cint(descriptor)
    result = ccall(:fsync, Cint, (Cint,), raw_descriptor)
    result == 0 || throw(SystemError("fsync $(label)", true))
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

function _fsync_metamdbg_directory(path::AbstractString)::Nothing
    directory = Base.Filesystem.open(String(path), Base.JL_O_RDONLY)
    try
        _fsync_metamdbg_descriptor(Base.fd(directory), path)
    finally
        close(directory)
    end
    return nothing
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
    published = false
    try
        write(temporary_io, reservation.output_root_reservation_contents)
        chmod(temporary_path, 0o600)
        _fsync_metamdbg_file(temporary_io, temporary_path)
        close(temporary_io)
        Base.Filesystem.hardlink(temporary_path, marker)
        published = true
        _fsync_metamdbg_directory(dirname(marker))
        rm(temporary_path)
        _require_metamdbg_output_root_reservation_marker!(reservation)
    catch primary_error
        isopen(temporary_io) && close(temporary_io)
        if ispath(temporary_path)
            try
                rm(temporary_path)
            catch cleanup_error
                @warn "metaMDBG failed to clean an unpublished shared " *
                      "output-root marker while preserving the primary " *
                      "publication failure" temporary_path primary_error cleanup_error
            end
        end
        if published && (ispath(marker) || islink(marker))
            try
                rm(marker)
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
    active_reservations = _metamdbg_submission_reservation_paths(outdir)
    isempty(active_reservations) || error(
        "metaMDBG output has an active nonlocal submission reservation: " *
        "$(join(active_reservations, ", ")). Refusing competing execution.",
    )
    return nothing
end

function _require_metamdbg_submission_reservation!(
        reservation::NamedTuple,
)::NamedTuple
    _require_metamdbg_output_root_reservation_marker!(reservation)
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

function _metamdbg_submission_reservation_from_path(
        path::AbstractString,
        canonical_outdir::AbstractString,
)::NamedTuple
    normalized_path = normpath(abspath(path))
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
    graph_k >= 0 || error(
        "metaMDBG submission reservation graph_k must be nonnegative.",
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
    expected.path == normalized_path || error(
        "metaMDBG submission reservation path does not match its workflow " *
        "signature.",
    )
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
        Set(String.(keys(parsed_job))) == Set(String[
            "schema_version",
            "workflow_signature",
            "owner_token",
            "job_id",
        ]) || error(
            "metaMDBG submission job record has unexpected fields: " *
            "$(expected.job_marker).",
        )
        job_id = get(parsed_job, "job_id", nothing)
        job_id isa AbstractString || error(
            "metaMDBG submission job record job_id must be a string.",
        )
        expected = _metamdbg_bound_submission_reservation(expected, job_id)
        expected.job_contents == job_contents || error(
            "metaMDBG submission job record is not canonical or does not " *
            "match its workflow owner.",
        )
    end
    return _require_metamdbg_submission_reservation!(expected)
end

function _create_metamdbg_submission_reservation!(
        reservation::NamedTuple,
        outdir::AbstractString,
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
    published = false
    output_root_marker_published = false
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
        _publish_metamdbg_output_root_reservation_marker!(reservation)
        output_root_marker_published = true
        _require_metamdbg_submission_reservation!(temporary_reservation)
        mv(temporary_path, reservation.path)
        published = true
        _fsync_metamdbg_directory(reservation_parent)
        _require_metamdbg_submission_reservation!(reservation)
    catch primary_error
        try
            if published && ispath(reservation.path)
                _remove_metamdbg_submission_reservation!(reservation)
            else
                if ispath(temporary_path)
                    rm(temporary_path; recursive = true)
                end
                if output_root_marker_published &&
                   (ispath(reservation.output_root_reservation_marker) ||
                    islink(reservation.output_root_reservation_marker))
                    rm(reservation.output_root_reservation_marker)
                    _fsync_metamdbg_directory(reservation_parent)
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
)::NamedTuple
    _require_metamdbg_submission_reservation!(reservation)
    bound_reservation =
        _metamdbg_bound_submission_reservation(reservation, job_id)
    marker = bound_reservation.job_marker
    if ispath(marker) || islink(marker)
        error(
            "metaMDBG refuses to overwrite an existing durable submission " *
            "job record: $(marker).",
        )
    end
    temporary_path, temporary_io = mktemp(dirname(reservation.path))
    published = false
    try
        write(temporary_io, bound_reservation.job_contents)
        chmod(temporary_path, 0o600)
        _fsync_metamdbg_file(temporary_io, temporary_path)
        close(temporary_io)
        Base.Filesystem.hardlink(temporary_path, marker)
        published = true
        _fsync_metamdbg_directory(reservation.path)
        rm(temporary_path)
        _fsync_metamdbg_directory(dirname(reservation.path))
        _require_metamdbg_submission_reservation!(bound_reservation)
    catch primary_error
        isopen(temporary_io) && close(temporary_io)
        if ispath(temporary_path)
            try
                rm(temporary_path)
            catch cleanup_error
                @warn "metaMDBG failed to clean an unpublished submission " *
                      "job record while preserving the primary bind failure" temporary_path primary_error cleanup_error
            end
        end
        if published && ispath(marker)
            try
                rm(marker)
            catch cleanup_error
                @warn "metaMDBG failed to remove an invalid newly published " *
                      "submission job record while preserving the primary " *
                      "bind failure" marker primary_error cleanup_error
            end
        end
        Base.rethrow()
    end
    return bound_reservation
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Bind a scheduler job ID to an inspected metaMDBG submission reservation.

This recovery operation is intended for the narrow crash window after the
scheduler accepted a job but before `run_metamdbg` durably recorded the returned
job ID. Set `confirm_submitted = true` only after independently confirming that
the exact scheduler job belongs to this reservation. The job record is published
atomically under the output lock and becomes available through
`inspect_metamdbg_submission_reservations`.
"""
function bind_metamdbg_submission_reservation_job!(
        metadata::NamedTuple;
        owner_token::AbstractString,
        job_id::AbstractString,
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
    metadata.job_id === nothing || error(
        "metaMDBG submission reservation already has a scheduler job id.",
    )
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
    bound = _with_metamdbg_output_lock(outputs.outdir) do
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
        current.job_id === nothing || error(
            "metaMDBG submission reservation already has a scheduler job id.",
        )
        _bind_metamdbg_submission_job!(current, job_id)
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
        submission_state = bound.submission_state,
    )
end

function _remove_metamdbg_submission_reservation!(
        reservation::NamedTuple,
)::Nothing
    output_root_marker = reservation.output_root_reservation_marker
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
    _require_metamdbg_submission_reservation!(reservation)
    tombstone_root = mktempdir(
        dirname(reservation.path);
        prefix = _metamdbg_submission_reservation_prefix(
            reservation.canonical_outdir,
        ) * "consumed.",
    )
    chmod(tombstone_root, 0o700)
    tombstone = joinpath(tombstone_root, "reservation")
    try
        mv(reservation.path, tombstone)
    catch primary_error
        try
            rm(tombstone_root; recursive = true)
        catch cleanup_error
            @warn "metaMDBG failed to clean an unused reservation tombstone " *
                  "while preserving the primary atomic-consumption failure" tombstone_root primary_error cleanup_error
        end
        Base.rethrow()
    end
    try
        rm(output_root_marker)
        _fsync_metamdbg_directory(dirname(output_root_marker))
    catch cleanup_error
        error(
            "metaMDBG atomically consumed its private submission " *
            "reservation but could not durably release its shared " *
            "output-root reservation. The output domain remains fail-closed " *
            "or its release durability is unknown: $(output_root_marker). " *
            "Cause: $(sprint(showerror, cleanup_error))",
        )
    end
    try
        rm(tombstone; recursive = true)
        rm(tombstone_root)
    catch cleanup_error
        @warn "metaMDBG atomically consumed its submission reservation but " *
              "left a nonblocking tombstone remnant" tombstone_root cleanup_error
    end
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
    try
        write(temporary_io, input_contract.contents)
        flush(temporary_io)
        close(temporary_io)
        chmod(temporary_path, 0o600)
        filesize(temporary_path) > 0 || error(
            "Failed to write metaMDBG provenance contract marker: $(marker).",
        )
        mv(temporary_path, marker)
    finally
        isopen(temporary_io) && close(temporary_io)
        rm(temporary_path; force = true)
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

function _require_valid_metamdbg_gfa(
        path::AbstractString,
        label::AbstractString,
)::String
    normalized_path = _require_metamdbg_regular_file(path, label)
    segment_identifiers = Set{String}()
    record_identifiers = Set{String}()
    open(normalized_path, "r") do input
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
    end
    isempty(segment_identifiers) && error(
        "$(label) contains no sequence-bearing GFA segments: " *
        "$(normalized_path).",
    )
    open(normalized_path, "r") do input
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
    published = false
    try
        write(temporary_io, completion.contents)
        flush(temporary_io)
        close(temporary_io)
        chmod(temporary_path, 0o600)
        Base.Filesystem.hardlink(temporary_path, marker)
        published = true
        rm(temporary_path)
        read(marker, String) == completion.contents || error(
            "Failed to write metaMDBG completion manifest: $(marker).",
        )
    catch primary_error
        isopen(temporary_io) && close(temporary_io)
        if ispath(temporary_path)
            try
                rm(temporary_path)
            catch cleanup_error
                @warn "metaMDBG failed to clean an unpublished completion " *
                      "manifest while preserving the primary publication " *
                      "failure" temporary_path primary_error cleanup_error
            end
        end
        if published && ispath(marker)
            try
                rm(marker)
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
        post_completion_publication_hook::Union{Nothing, String} = nothing,
        lock_retry_attempts::Int = _METAMDBG_OUTPUT_LOCK_RETRY_ATTEMPTS,
        lock_retry_delay_seconds::Real =
            _METAMDBG_OUTPUT_LOCK_RETRY_DELAY_SECONDS,
)::String
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
    outdir_base = basename(outputs.outdir)
    lock_path = _metamdbg_output_lock_path(outputs.outdir)
    effective_submission_reservation = if submission_reservation === nothing
        _metamdbg_submission_reservation(outputs, input_contract, graph_k)
    else
        submission_reservation
    end
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
        staged_name =
            "$(lpad(input_index, 6, '0'))--$(basename(input.path))"
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
    reservation_declaration_lines = String[
        "submission_reservation_dir=$(_metamdbg_shell_literal(effective_submission_reservation.path))",
        "submission_reservation_contract=$(_metamdbg_shell_literal(effective_submission_reservation.contract_marker))",
        "submission_reservation_job=$(_metamdbg_shell_literal(effective_submission_reservation.job_marker))",
        "submission_output_root_reservation=$(_metamdbg_shell_literal(effective_submission_reservation.output_root_reservation_marker))",
        "runtime_output_root_reservation=$(_metamdbg_shell_literal(effective_submission_reservation.runtime_output_root_reservation_marker))",
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
        "reservation_entry_count=\$(find \"\$submission_reservation_dir\" -mindepth 1 -maxdepth 1 -print | awk 'END { print NR }')",
        "if [ \"\$reservation_entry_count\" -ne 2 ] || [ ! -f \"\$submission_reservation_job\" ] || [ -L \"\$submission_reservation_job\" ] || [ ! -s \"\$submission_reservation_job\" ]; then",
        "  echo \"metaMDBG runtime refuses an unbound or invalid submission reservation\" >&2",
        "  return 1",
        "fi",
        "if [ \"\$(metamdbg_file_uid \"\$submission_reservation_job\")\" != \"\$(id -u)\" ] || [ \"\$(metamdbg_file_mode \"\$submission_reservation_job\")\" != \"600\" ]; then",
        "  echo \"metaMDBG submission reservation job record owner or mode is invalid\" >&2",
        "  return 1",
        "fi",
        "if [ -z \"\${SLURM_JOB_ID:-}\" ]; then",
        "  echo \"metaMDBG bound submission reservation requires SLURM_JOB_ID\" >&2",
        "  return 1",
        "fi",
        "printf '%s' $(_metamdbg_shell_literal(first(job_record_parts))) > \"\$expected_reservation_job\"",
        "printf '%s' \"\$SLURM_JOB_ID\" >> \"\$expected_reservation_job\"",
        "printf '%s' $(_metamdbg_shell_literal(last(job_record_parts))) >> \"\$expected_reservation_job\"",
        "cmp -s -- \"\$expected_reservation_job\" \"\$submission_reservation_job\" || {",
        "  echo \"metaMDBG submission reservation job record does not match SLURM_JOB_ID\" >&2",
        "  return 1",
        "}",
        "}",
    ]
    output_root_domain_function_lines = String[
        "require_exclusive_output_root_domain() {",
        "  local reservation_root=\"\$outdir\"",
        "  local reservation_parent",
        "  local reservation_base",
        "  local pid_reservation",
        "  local durable_prefix",
        "  local candidate",
        "  while [ \"\$(dirname -- \"\$reservation_root\")\" != \"\$reservation_root\" ]; do",
        "    reservation_parent=\$(dirname -- \"\$reservation_root\")",
        "    reservation_base=\$(basename -- \"\$reservation_root\")",
        "    pid_reservation=\"\$reservation_parent/.\$reservation_base$(_OUTPUT_ROOT_RESERVATION_LOCK_SUFFIX)\"",
        "    if [ -e \"\$pid_reservation\" ] || [ -L \"\$pid_reservation\" ]; then",
        "      echo \"metaMDBG overlaps an active shared output-root reservation: \$pid_reservation\" >&2",
        "      return 1",
        "    fi",
        "    durable_prefix=\"\$reservation_parent/.\$reservation_base$(_OUTPUT_ROOT_DURABLE_RESERVATION_SEPARATOR)\"",
        "    for candidate in \"\$durable_prefix\"*; do",
        "      if [ ! -e \"\$candidate\" ] && [ ! -L \"\$candidate\" ]; then",
        "        continue",
        "      fi",
        "      if [ \"\$candidate\" = \"\$submission_output_root_reservation\" ] || [ \"\$candidate\" = \"\$runtime_output_root_reservation\" ]; then",
        "        continue",
        "      fi",
        "      echo \"metaMDBG overlaps an active shared output-root reservation: \$candidate\" >&2",
        "      return 1",
        "    done",
        "    reservation_root=\"\$reservation_parent\"",
        "  done",
        "  if [ -d \"\$outdir\" ] && [ ! -L \"\$outdir\" ]; then",
        "    while IFS= read -r -d '' candidate; do",
        "      if [ \"\$candidate\" = \"\$submission_output_root_reservation\" ] || [ \"\$candidate\" = \"\$runtime_output_root_reservation\" ]; then",
        "        continue",
        "      fi",
        "      echo \"metaMDBG overlaps an active descendant shared output-root reservation: \$candidate\" >&2",
        "      return 1",
        "    done < <(find -P \"\$outdir\" \\( -name '.*$(_OUTPUT_ROOT_RESERVATION_LOCK_SUFFIX)' -o -name '.*$(_OUTPUT_ROOT_DURABLE_RESERVATION_SEPARATOR)*' \\) -print0)",
        "  fi",
        "}",
    ]
    reservation_consumption_lines = String[
        "validate_submission_reservation",
        "reservation_tombstone=\"\$secure_tmpdir/consumed-reservation\"",
        "mv -- \"\$submission_reservation_dir\" \"\$reservation_tombstone\"",
        "rm -- \"\$submission_output_root_reservation\"",
        "rm -rf -- \"\$reservation_tombstone\"",
    ]
    inventory_canonicalizer = _metamdbg_shell_literal(
        _metamdbg_runtime_inventory_canonicalizer_python(),
    )
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
        "outdir_base=$(_metamdbg_shell_literal(outdir_base))",
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
        "lock_acquired=0",
        "output_root_runtime_reservation_acquired=0",
        "lock_retry_attempts=$(lock_retry_attempts)",
        "lock_retry_delay_seconds=$(Float64(lock_retry_delay_seconds))",
        "sha256_file() {",
        "  local path=\"\$1\"",
        "  if command -v sha256sum >/dev/null 2>&1; then",
        "    sha256sum -- \"\$path\" | awk '{print \$1}'",
        "  elif command -v shasum >/dev/null 2>&1; then",
        "    shasum -a 256 -- \"\$path\" | awk '{print \$1}'",
        "  else",
        "    echo \"metaMDBG requires sha256sum or shasum to validate content\" >&2",
        "    return 1",
        "  fi",
        "}",
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
        "sha256_text() {",
        "  printf '%s' \"\$1\" | sha256_stream",
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
        "  trap - EXIT INT TERM HUP",
        "  if [ -n \"\$secure_tmpdir\" ]; then",
        "    if ! rm -rf -- \"\$secure_tmpdir\"; then",
        "      echo \"metaMDBG failed to remove secure lifecycle temporary directory\" >&2",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    fi",
        "  fi",
        "  if [ \"\$output_root_runtime_reservation_acquired\" -eq 1 ]; then",
        "    if ! rmdir -- \"\$runtime_output_root_reservation\"; then",
        "      echo \"metaMDBG failed to release shared output-root reservation: \$runtime_output_root_reservation\" >&2",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    fi",
        "  fi",
        "  if [ \"\$lock_acquired\" -eq 1 ]; then",
        "    if ! rmdir -- \"\$lock_dir\"; then",
        "      echo \"metaMDBG failed to release lifecycle lock: \$lock_dir\" >&2",
        "      [ \"\$status\" -ne 0 ] || status=1",
        "    fi",
        "  fi",
        "  exit \"\$status\"",
        "}",
        "trap cleanup_metamdbg_lifecycle EXIT",
        "trap 'exit 129' HUP",
        "trap 'exit 130' INT",
        "trap 'exit 143' TERM",
        "secure_tmpdir=\$(mktemp -d \"\$outdir_parent/.\${outdir_base}.mycelia.XXXXXXXXXX\")",
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
        "printf '%s' $(_metamdbg_shell_literal(effective_submission_reservation.contents)) > \"\$expected_reservation\"",
        "chmod 600 \"\$expected_reservation\"",
        "printf '%s' $(_metamdbg_shell_literal(effective_submission_reservation.output_root_reservation_contents)) > \"\$expected_output_root_reservation\"",
        "chmod 600 \"\$expected_output_root_reservation\"",
        "validate_submission_reservation",
        "lock_attempt=1",
        "while ! mkdir -m 700 -- \"\$lock_dir\" 2>/dev/null; do",
        "  validate_submission_reservation",
        "  if [ \"\$lock_attempt\" -ge \"\$lock_retry_attempts\" ]; then",
        "    echo \"metaMDBG output remained locked after \$lock_retry_attempts attempts: \$lock_dir\" >&2",
        "    exit 1",
        "  fi",
        "  sleep \"\$lock_retry_delay_seconds\"",
        "  lock_attempt=\$((lock_attempt + 1))",
        "done",
        "lock_acquired=1",
        "if ! mkdir -m 700 -- \"\$runtime_output_root_reservation\" 2>/dev/null; then",
        "  echo \"metaMDBG could not acquire its exact runtime output-root reservation: \$runtime_output_root_reservation\" >&2",
        "  exit 1",
        "fi",
        "output_root_runtime_reservation_acquired=1",
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
        "contract_exists=0",
        "if [ -n \"\$(find \"\$outdir\" -mindepth 1 -maxdepth 1 -print -quit)\" ]; then",
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
        "  partial_entry=\$(find \"\$outdir\" -mindepth 1 -maxdepth 1 ! -name '$(_METAMDBG_CONTRACT_FILENAME)' -print -quit)",
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
        "  require_unchanged_metamdbg_output_root",
        "  mv -n -- \"\$contract_new\" \"\$contract_marker\"",
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
        "mv -n -- \"\$completion_new\" \"\$completion_marker\"",
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
    toolchain = _metamdbg_expected_toolchain()
    workflow = _metamdbg_workflow_contract(outputs, input_contract, graph_k)
    provenance = (;
        contract_signature = input_contract.signature,
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
    job_id = submission isa Mycelia.SubmitResult ? submission.job_id : nothing
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
            submission_state = submission_reservation.submission_state,
        ),
    ))
end

function _metamdbg_existing_artifacts(
        outputs::NamedTuple,
        graph_k::Int,
        ;
        output_root_identity::Union{Nothing, NamedTuple} = nothing,
)::Union{Nothing, NamedTuple}
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
        String(job_id) == normalized_job_id || error(
            "metaMDBG sbatch submission returned a non-exact job id: " *
            "$(repr(job_id)).",
        )
        if submission.stdout !== nothing
            stdout_job_id = Mycelia._extract_sbatch_job_id(submission.stdout)
            stdout_job_id == normalized_job_id || error(
                "metaMDBG sbatch submission stdout did not contain the exact " *
                "returned job id $(normalized_job_id).",
            )
        end
    end
    return submission
end

function _metamdbg_submission_recovery_guidance(
        reservation::NamedTuple,
        job_id::AbstractString,
)::String
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    return "SLURM job $(normalized_job_id) remains held or has an unknown " *
           "release state. Its durable reservation is $(reservation.path) " *
           "and scheduler name is $(reservation.scheduler_job_name). Inspect " *
           "that exact job with `scontrol show job $(normalized_job_id)` before " *
           "running `scontrol release $(normalized_job_id)` or cancelling and " *
           "reclaiming the reservation."
end

function _release_metamdbg_submission_job!(
        reservation::NamedTuple,
        job_id::AbstractString,
        release_runner::Function,
)::NamedTuple
    normalized_job_id = _normalize_metamdbg_job_id(job_id)
    release_result = try
        release_runner(normalized_job_id)
    catch caught
        error(
            "metaMDBG failed to release its durably bound scheduler job. " *
            _metamdbg_submission_recovery_guidance(
                reservation,
                normalized_job_id,
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
            ),
        )
    end
    return ispath(reservation.path) ? reservation :
           merge(reservation, (; submission_state = :consumed))
end

function _cleanup_metamdbg_submission_reservation_after_failure!(
        outputs::NamedTuple,
        reservation::NamedTuple,
)::Nothing
    if !ispath(reservation.path) && !islink(reservation.path)
        return nothing
    end
    _with_metamdbg_output_lock(outputs.outdir) do
        _remove_metamdbg_submission_reservation!(reservation)
    end
    return nothing
end

function _bind_metamdbg_submission_job_after_submit!(
        reservation::NamedTuple,
        job_id::AbstractString,
)::NamedTuple
    return _with_metamdbg_output_lock(reservation.canonical_outdir) do
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
        return _bind_metamdbg_submission_job!(current, job_id)
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
    selected_input = _metamdbg_selected_input(hifi_reads, ont_reads)
    abundance_min > 0 || throw(ArgumentError("abundance_min must be positive."))
    threads > 0 || throw(ArgumentError("threads must be positive."))
    graph_k >= 0 || throw(ArgumentError("graph_k must be nonnegative."))
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
                end
                _require_metamdbg_contract!(outputs, input_contract)
                _require_unchanged_metamdbg_output_root(
                    output_root_identity,
                    outputs.outdir,
                )
                _write_metamdbg_completion_manifest!(outputs, completion)
                wrote_completion = true
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
                        rm(outputs.completion_marker)
                    catch cleanup_error
                        @warn "metaMDBG failed to remove a newly written stale " *
                              "completion manifest while preserving the " *
                              "primary completion failure" outputs primary_error cleanup_error
                    end
                end
                if wrote_contract
                    try
                        rm(outputs.contract_marker)
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
    bound_submission_reservation = try
        submission_job_binder(
            submission_reservation,
            submission_job_id,
        )
    catch caught
        error(
            "metaMDBG submitted SLURM job $(submission_job_id) held but could " *
            "not durably bind job.json; the job was not released. " *
            _metamdbg_submission_recovery_guidance(
                submission_reservation,
                submission_job_id,
            ) * " Cause: $(sprint(showerror, caught))",
        )
    end
    released_submission_reservation = _release_metamdbg_submission_job!(
        bound_submission_reservation,
        submission_job_id,
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
on-disk owner contract under the output lock. This is the recovery path when a
caller dies after publishing the reservation but before receiving the returned
capability. Inspection never expires or removes a reservation.
"""
function inspect_metamdbg_submission_reservations(
        outdir::AbstractString,
)::Vector{NamedTuple}
    canonical_outdir = _metamdbg_canonical_output_path(outdir)
    return _with_metamdbg_output_lock(canonical_outdir) do
        map(_metamdbg_submission_reservation_paths(canonical_outdir)) do path
            reservation = _metamdbg_submission_reservation_from_path(
                path,
                canonical_outdir,
            )
            return (;
                canonical_outdir = reservation.canonical_outdir,
                path = reservation.path,
                workflow_signature = reservation.workflow_signature,
                scheduler_job_name = reservation.scheduler_job_name,
                input_contract_signature =
                    reservation.input_contract_signature,
                graph_k = reservation.graph_k,
                owner_token = reservation.owner_token,
                job_id = reservation.job_id,
                submission_state = reservation.submission_state,
            )
        end
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Reclaim a durable metaMDBG submission reservation after explicit confirmation.

This is an explicit cancellation capability, not an expiry mechanism. Pass the
`submission_reservation` metadata returned by `run_metamdbg` with
`status = :submitted`, the exact random `owner_token`, the exact scheduler
`job_id`, and `confirm_cancelled = true` only after the scheduler has confirmed
that the job cannot start. For a process death before submission, first call
`inspect_metamdbg_submission_reservations`, independently confirm that no job
was submitted, then pass that inspected record, its exact owner token, and
`confirm_not_submitted = true`. The confirmation modes are mutually exclusive.
If the exact submitted job instead reaches a terminal failed state,
pass its exact `job_id` and `confirm_terminal = :failed` after independently
confirming that scheduler state. The reservation is removed under the output
lock only when the
on-disk owner record still matches exactly. Missing, consumed, or
replacement-owner reservations fail loudly and are never removed automatically.
"""
function reclaim_metamdbg_submission_reservation!(
        metadata::NamedTuple,
        ;
        owner_token::AbstractString,
        job_id::Union{Nothing, AbstractString} = nothing,
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
        "confirm_not_submitted=true, or confirm_terminal=:failed after " *
        "independently verifying the corresponding scheduler state.",
    ))
    confirm_terminal === nothing || confirm_terminal == :failed || throw(
        ArgumentError(
            "metaMDBG confirm_terminal currently accepts only :failed.",
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
    normalized_owner_token = String(owner_token)
    isempty(normalized_owner_token) && throw(ArgumentError(
        "metaMDBG reservation owner_token must be nonempty.",
    ))
    normalized_owner_token == metadata.owner_token || error(
        "metaMDBG reservation owner token does not match the durable " *
        "reservation capability.",
    )
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
        )
    end
    expected_reservation.workflow_signature == metadata.workflow_signature ||
        error("metaMDBG reservation workflow signature does not match metadata.")
    expected_reservation.path == metadata.path || error(
        "metaMDBG reservation path does not match its recomputed workflow path.",
    )
    ispath(expected_reservation.path) || error(
        "metaMDBG submission reservation is missing or was already consumed: " *
        "$(expected_reservation.path).",
    )
    _with_metamdbg_output_lock(outputs.outdir) do
        ispath(expected_reservation.path) || error(
            "metaMDBG submission reservation is missing or was already consumed: " *
            "$(expected_reservation.path).",
        )
        _remove_metamdbg_submission_reservation!(expected_reservation)
    end
    return (;
        status = :reclaimed,
        job_id = normalized_job_id,
        recovery_reason = if confirm_cancelled
            :cancelled
        elseif confirm_terminal !== nothing
            confirm_terminal
        else
            :not_submitted
        end,
        path = expected_reservation.path,
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
- `outdir::String`: Output directory path (default: "metamdbg_output").
- `abundance_min::Int`: Minimum abundance threshold (default: 3).
- `threads::Int`: Number of threads to use (default: get_default_threads()).
- `graph_k::Int`: Graph resolution requested from `metaMDBG gfa` (default: 21).

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
  serialized into the durable schema-v4 contract. Tool execution uses private,
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
  reservations block competing execution before input hashing. Runtime jobs
  validate their exact reservation before bounded output-lock retries and
  atomically rename it to a private tombstone only after acquiring the lock.
  Reservation creation likewise publishes a complete marker by atomic rename;
  incomplete temporary or consumed-tombstone remnants are nonblocking.
  Reservations never auto-expire:
  after scheduler-confirmed cancellation, reclaim one explicitly with
  `reclaim_metamdbg_submission_reservation!`, the returned owner token and job
  id, and `confirm_cancelled = true`. If a caller dies before submission,
  recover the durable capability with
  `inspect_metamdbg_submission_reservations` and reclaim it only after explicit
  independent confirmation that no job was submitted.
  A submitted job confirmed terminal-failed can be reclaimed with its exact job
  id and `confirm_terminal = :failed`.
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
