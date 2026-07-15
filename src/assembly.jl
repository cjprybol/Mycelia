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
function run_unicycler(;
        short_1,
        short_2 = nothing,
        long_reads,
        outdir = "unicycler_output",
        threads::Int = get_default_threads(),
        spades_options::Union{Nothing, String} = nothing,
        kmers::Union{Nothing, String} = nothing,
        executor = nothing,
        site::Symbol = :local,
        job_name::String = "unicycler",
        time_limit::String = "2-00:00:00",
        partition::Union{Nothing, String} = nothing,
        account::Union{Nothing, String} = nothing,
        mem_gb::Union{Nothing, Real} = nothing,
        qos::Union{Nothing, String} = nothing,
        mail_user::Union{Nothing, String} = nothing)
    Mycelia.add_bioconda_env("unicycler")

    if executor !== nothing
        spades_args = isnothing(spades_options) ? String[] :
                      ["--spades_options=$(spades_options)"]
        kmer_args = isnothing(kmers) ? String[] : ["--kmers=$(kmers)"]
        cmd = if isnothing(short_2)
            Mycelia.command_string(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n unicycler unicycler -s $(short_1) -l $(long_reads) -o $(outdir) -t $(threads) $(kmer_args...) $(spades_args...)`
            )
        else
            Mycelia.command_string(
                `$(Mycelia.CONDA_RUNNER) run --live-stream -n unicycler unicycler -1 $(short_1) -2 $(short_2) -l $(long_reads) -o $(outdir) -t $(threads) $(kmer_args...) $(spades_args...)`
            )
        end
        script = join([
            "set -euo pipefail",
            "if [ ! -f \"$(joinpath(outdir, "assembly.fasta"))\" ]; then",
            "  rm -rf \"$(outdir)\" || true",
            "  $(cmd)",
            "fi",
            "mkdir -p \"$(outdir)\""
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
            graph = joinpath(outdir, "assembly.gfa"))
    end

    # Unicycler requires the output directory to not exist, so check output file first
    if !isfile(joinpath(outdir, "assembly.fasta"))
        # Remove output directory if it exists (Unicycler will create it)
        if isdir(outdir)
            rm(outdir, recursive = true)
        end

        spades_args = isnothing(spades_options) ? String[] :
                      ["--spades_options=$(spades_options)"]
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
    return (; outdir, assembly = joinpath(outdir, "assembly.fasta"),
        graph = joinpath(outdir, "assembly.gfa"))
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
    "c6fbbb3c2e85ffae22764424333cda0937bfcf13187f9a40c20282c40736db3f"
const METAMDBG_ENV_NAME =
    "metamdbg-$(METAMDBG_VERSION)-$(first(METAMDBG_ENVIRONMENT_SPEC_SHA256, 16))"
const _METAMDBG_INSTALL_LOCK_STALE_SECONDS = 600
const _METAMDBG_INSTALL_LOCK_REFRESH_SECONDS = 60
const _METAMDBG_CONTRACT_SCHEMA_VERSION = 3
const _METAMDBG_CONTRACT_FILENAME = "mycelia_metamdbg_contract.json"

function _metamdbg_paths()::Tuple{String, String}
    install_dir = joinpath(dirname(dirname(pathof(Mycelia))), "deps", "metamdbg")
    return install_dir, joinpath(install_dir, "environment.yml")
end

function _metamdbg_sha256(path::AbstractString)::String
    return open(path, "r") do input
        SHA.bytes2hex(SHA.sha256(input))
    end
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

function _metamdbg_environment_packages(;
        conda_runner::AbstractString = CONDA_RUNNER,
        command_reader::Function = command -> read(command, String),
)::Dict{String, String}
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
    versions = Dict{String, String}()
    for package_record in package_records
        package_record isa AbstractDict || continue
        name = get(package_record, "name", nothing)
        version = get(package_record, "version", nothing)
        if name isa AbstractString && version isa AbstractString
            versions[String(name)] = String(version)
        end
    end
    return versions
end

function _require_metamdbg_package_version(
        versions::AbstractDict{<:AbstractString, <:AbstractString},
)::NamedTuple
    actual_version = get(versions, "metamdbg", nothing)
    displayed_version =
        actual_version === nothing ? "missing" : repr(actual_version)
    actual_version == METAMDBG_VERSION || error(
        "metaMDBG spec-addressed environment $(repr(METAMDBG_ENV_NAME)) must " *
        "contain metamdbg exactly $(METAMDBG_VERSION), got " *
        "$(displayed_version). " *
        "Refusing to repair it in place; remove it manually before reinstalling.",
    )
    return _metamdbg_expected_toolchain()
end

function _require_expected_metamdbg_toolchain(toolchain::Any)::NamedTuple
    expected = _metamdbg_expected_toolchain()
    toolchain isa NamedTuple || error(
        "metaMDBG dependency validation did not return toolchain provenance.",
    )
    toolchain == expected || error(
        "metaMDBG dependency validation returned incompatible toolchain " *
        "provenance: expected $(repr(expected)), got $(repr(toolchain)).",
    )
    return expected
end

function _metamdbg_install_lock_path()::String
    return joinpath(
        first(Base.DEPOT_PATH),
        "locks",
        "mycelia-$(METAMDBG_ENV_NAME)-install.pid",
    )
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
        paths::Tuple{String, String} = _metamdbg_paths(),
        environment_checker::Function = () ->
            check_bioconda_env_is_installed(METAMDBG_ENV_NAME),
        environment_creator::Function = create_conda_env_from_yaml,
        package_inspector::Function = _metamdbg_environment_packages,
)::NamedTuple
    install_dir, environment_path = paths
    mkpath(install_dir)
    verified_environment_path =
        _require_verified_metamdbg_environment_spec(environment_path)
    if !Bool(environment_checker())
        environment_creator(
            verified_environment_path,
            METAMDBG_ENV_NAME;
            force = false,
        )
    end
    return _require_metamdbg_package_version(package_inspector())
end

function _ensure_metamdbg_installed(;
        paths::Tuple{String, String} = _metamdbg_paths(),
        environment_checker::Function = () ->
            check_bioconda_env_is_installed(METAMDBG_ENV_NAME),
        environment_creator::Function = create_conda_env_from_yaml,
        package_inspector::Function = _metamdbg_environment_packages,
        lock_path::AbstractString = _metamdbg_install_lock_path(),
        lock_runner::Function = _with_metamdbg_install_lock,
)::NamedTuple
    return lock_runner(String(lock_path)) do
        _ensure_metamdbg_installed_locked(;
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

function _metamdbg_input_contract(
        selected_input::NamedTuple,
        abundance_min::Int,
        toolchain::NamedTuple = _metamdbg_expected_toolchain(),
)::NamedTuple
    input_technology = if selected_input.flag == "--in-hifi"
        "hifi"
    elseif selected_input.flag == "--in-ont"
        "ont"
    else
        error("Unsupported metaMDBG input flag: $(selected_input.flag)")
    end
    inputs = map(selected_input.paths) do path
        file_status = stat(path)
        modification_time_ns = round(
            Int64,
            file_status.mtime * 1_000_000_000,
        )
        return (
            path = path,
            size_bytes = Int64(file_status.size),
            modification_time_ns,
            sha256 = _metamdbg_sha256(path),
        )
    end
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
    return (; contract, signature, contents)
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
    try
        for record in reader
            record isa FASTX.FASTA.Record || error(
                "$(label) is not valid FASTA: $(normalized_path).",
            )
            sequence = FASTX.sequence(String, record)
            isempty(sequence) && error(
                "$(label) contains an empty FASTA sequence at record " *
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

function _require_valid_metamdbg_gfa(
        path::AbstractString,
        label::AbstractString,
)::String
    normalized_path = _require_metamdbg_regular_file(path, label)
    segment_identifiers = Set{String}()
    open(normalized_path, "r") do input
        for (line_number, line) in enumerate(eachline(input))
            isempty(line) && continue
            fields = split(line, '\t'; keepempty = true)
            first(fields) == "S" || continue
            length(fields) >= 3 || error(
                "$(label) has a malformed GFA segment at line " *
                "$(line_number): $(normalized_path).",
            )
            identifier = String(fields[2])
            sequence = String(fields[3])
            isempty(identifier) && error(
                "$(label) has an empty GFA segment identifier at line " *
                "$(line_number): $(normalized_path).",
            )
            identifier in segment_identifiers && error(
                "$(label) has duplicate GFA segment identifier " *
                "$(repr(identifier)): $(normalized_path).",
            )
            if isempty(sequence) || sequence == "*"
                error(
                    "$(label) has no sequence for GFA segment " *
                    "$(repr(identifier)): $(normalized_path).",
                )
            end
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
            push!(segment_identifiers, identifier)
        end
    end
    isempty(segment_identifiers) && error(
        "$(label) contains no sequence-bearing GFA segments: " *
        "$(normalized_path).",
    )
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
)::NamedTuple
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
    return (; outdir = outputs.outdir, contigs, graph)
end

function _metamdbg_executor_script(
        asm_cmd::String,
        gfa_cmd::String,
        outputs::NamedTuple,
        graph_k::Int,
        input_contract::NamedTuple,
        ;
        conda_runner::AbstractString = CONDA_RUNNER,
        environment_path::AbstractString = last(_metamdbg_paths()),
)::String
    verified_environment_path =
        _require_verified_metamdbg_environment_spec(environment_path)
    outdir_parent = dirname(outputs.outdir)
    outdir_base = basename(outputs.outdir)
    lock_path = _metamdbg_output_lock_path(outputs.outdir)
    input_validation_lines = String[]
    for (input_index, input) in enumerate(input_contract.contract.inputs)
        input_path_variable = "input_path_$(input_index)"
        expected_digest_variable = "expected_input_sha256_$(input_index)"
        actual_digest_variable = "actual_input_sha256_$(input_index)"
        append!(input_validation_lines, String[
            "$(input_path_variable)=$(Base.shell_escape(input.path))",
            "$(expected_digest_variable)=$(Base.shell_escape(input.sha256))",
            "if [ ! -f \"\$$(input_path_variable)\" ]; then",
            "  echo \"metaMDBG input is missing before execution: \$$(input_path_variable)\" >&2",
            "  exit 1",
            "fi",
            "$(actual_digest_variable)=\$(sha256_file \"\$$(input_path_variable)\")",
            "if [ \"\$$(actual_digest_variable)\" != \"\$$(expected_digest_variable)\" ]; then",
            "  echo \"metaMDBG input content changed before execution: \$$(input_path_variable)\" >&2",
            "  exit 1",
            "fi",
        ])
    end
    lines = String[
        "set -euo pipefail",
        "umask 077",
        "outdir=$(Base.shell_escape(outputs.outdir))",
        "outdir_parent=$(Base.shell_escape(outdir_parent))",
        "outdir_base=$(Base.shell_escape(outdir_base))",
        "lock_dir=$(Base.shell_escape(lock_path))",
        "contigs_plain=\"\$outdir/contigs.fasta\"",
        "contigs_gz=\"\$outdir/contigs.fasta.gz\"",
        "graph_alias=\"\$outdir/assemblyGraph_k$(graph_k).gfa\"",
        "contract_marker=\"\$outdir/$(_METAMDBG_CONTRACT_FILENAME)\"",
        "environment_spec=$(Base.shell_escape(verified_environment_path))",
        "expected_spec_sha256=$(Base.shell_escape(METAMDBG_ENVIRONMENT_SPEC_SHA256))",
        "expected_metamdbg_version=$(Base.shell_escape(METAMDBG_VERSION))",
        "environment_name=$(Base.shell_escape(METAMDBG_ENV_NAME))",
        "conda_runner=$(Base.shell_escape(conda_runner))",
        "secure_tmpdir=",
        "lock_acquired=0",
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
        "mkdir -p -- \"\$outdir_parent\"",
        "runtime_outdir_parent=\$(cd -- \"\$outdir_parent\" && pwd -P)",
        "if [ \"\$runtime_outdir_parent\" != \"\$outdir_parent\" ]; then",
        "  echo \"metaMDBG canonical output parent changed before execution\" >&2",
        "  exit 1",
        "fi",
        "if ! mkdir -m 700 -- \"\$lock_dir\"; then",
        "  echo \"metaMDBG output is locked by another lifecycle: \$lock_dir\" >&2",
        "  exit 1",
        "fi",
        "lock_acquired=1",
        "cleanup_metamdbg_lifecycle() {",
        "  local status=\$?",
        "  trap - EXIT INT TERM HUP",
        "  if [ -n \"\$secure_tmpdir\" ]; then",
        "    if ! rm -rf -- \"\$secure_tmpdir\"; then",
        "      echo \"metaMDBG failed to remove secure lifecycle temporary directory\" >&2",
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
        input_validation_lines...,
        "if [ -L \"\$outdir\" ]; then",
        "  echo \"metaMDBG outdir must not be a symbolic link: \$outdir\" >&2",
        "  exit 1",
        "fi",
        "if [ -e \"\$outdir\" ] && [ ! -d \"\$outdir\" ]; then",
        "  echo \"metaMDBG outdir exists but is not a directory: \$outdir\" >&2",
        "  exit 1",
        "fi",
        "mkdir -p -- \"\$outdir\"",
        "secure_tmpdir=\$(mktemp -d \"\$outdir_parent/.\${outdir_base}.mycelia.XXXXXXXXXX\")",
        "chmod 700 \"\$secure_tmpdir\"",
        "expected_contract=\"\$secure_tmpdir/expected-contract.json\"",
        "contract_new=\"\$secure_tmpdir/contract.new\"",
        "contigs_new=\"\$secure_tmpdir/contigs.fasta.new.gz\"",
        "package_inventory=\"\$secure_tmpdir/conda-packages.json\"",
        "printf '%s' $(Base.shell_escape(input_contract.contents)) > \"\$expected_contract\"",
        "chmod 600 \"\$expected_contract\"",
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
        "if [ ! -f \"\$environment_spec\" ] || [ -L \"\$environment_spec\" ]; then",
        "  echo \"metaMDBG environment specification is missing or not regular\" >&2",
        "  exit 1",
        "fi",
        "actual_spec_sha256=\$(sha256_file \"\$environment_spec\")",
        "if [ \"\$actual_spec_sha256\" != \"\$expected_spec_sha256\" ]; then",
        "  echo \"metaMDBG environment specification checksum mismatch\" >&2",
        "  exit 1",
        "fi",
        "\"\$conda_runner\" list -n \"\$environment_name\" '^metamdbg\$' --json > \"\$package_inventory\"",
        "grep -Eq '\"name\"[[:space:]]*:[[:space:]]*\"metamdbg\"' \"\$package_inventory\" || {",
        "  echo \"metaMDBG package is missing from environment \$environment_name\" >&2",
        "  exit 1",
        "}",
        "grep -Eq '\"version\"[[:space:]]*:[[:space:]]*\"1\\.4\"' \"\$package_inventory\" || {",
        "  echo \"metaMDBG environment must contain metamdbg exactly \$expected_metamdbg_version\" >&2",
        "  exit 1",
        "}",
        "validate_fasta_stream() {",
        "  awk '",
        "    BEGIN { records = 0; sequence_bases = 0 }",
        "    { sub(/\\r\$/, \"\", \$0) }",
        "    /^>/ {",
        "      if (records > 0 && sequence_bases == 0) exit 10",
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
        "}",
        "validate_gfa() {",
        "  local path=\$1",
        "  if [ ! -f \"\$path\" ] || [ -L \"\$path\" ] || [ ! -s \"\$path\" ]; then",
        "    echo \"metaMDBG graph is missing, empty, or not regular: \$path\" >&2",
        "    return 1",
        "  fi",
        "  awk -F '\\t' '",
        "    \$1 == \"S\" {",
        "      sub(/\\r\$/, \"\", \$3)",
        "      if (NF < 3 || \$2 == \"\" || \$3 == \"\" || \$3 == \"*\") exit 20",
        "      if (\$3 !~ /^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+\$/) exit 21",
        "      if (seen[\$2]++) exit 22",
        "      segments += 1",
        "    }",
        "    END { if (segments == 0) exit 23 }",
        "  ' \"\$path\" || {",
        "    echo \"metaMDBG graph has no unique sequence-bearing GFA segments: \$path\" >&2",
        "    return 1",
        "  }",
        "}",
        "find_metamdbg_graph() {",
        "  local candidates=()",
        "  local candidate",
        "  shopt -s nullglob",
        "  for candidate in \"\$outdir\"/assemblyGraph_k$(graph_k)_*bps.gfa; do",
        "    if [ ! -f \"\$candidate\" ] || [ -L \"\$candidate\" ] || [ ! -s \"\$candidate\" ]; then",
        "      echo \"metaMDBG graph candidate is empty or not regular: \$candidate\" >&2",
        "      shopt -u nullglob",
        "      return 2",
        "    fi",
        "    candidates+=(\"\$candidate\")",
        "  done",
        "  shopt -u nullglob",
        "  if [ \"\${#candidates[@]}\" -eq 0 ]; then",
        "    return 1",
        "  elif [ \"\${#candidates[@]}\" -ne 1 ]; then",
        "    echo \"metaMDBG produced multiple graphs for k=$(graph_k)\" >&2",
        "    return 2",
        "  fi",
        "  printf '%s\\n' \"\${candidates[0]}\"",
        "}",
        "if [ -e \"\$contigs_gz\" ] || [ -L \"\$contigs_gz\" ]; then",
        "  validate_contigs \"\$contigs_gz\"",
        "elif [ -e \"\$contigs_plain\" ] || [ -L \"\$contigs_plain\" ]; then",
        "  validate_contigs \"\$contigs_plain\"",
        "else",
        "  $(asm_cmd)",
        "fi",
        "if [ -e \"\$contigs_gz\" ] || [ -L \"\$contigs_gz\" ]; then",
        "  validate_contigs \"\$contigs_gz\"",
        "else",
        "  validate_contigs \"\$contigs_plain\"",
        "  gzip -c -- \"\$contigs_plain\" > \"\$contigs_new\"",
        "  chmod 600 \"\$contigs_new\"",
        "  validate_contigs \"\$contigs_new\"",
        "  mv -f -- \"\$contigs_new\" \"\$contigs_gz\"",
        "fi",
        "validate_contigs \"\$contigs_gz\"",
        "graph_status=0",
        "graph_source=\"\"",
        "graph_source=\$(find_metamdbg_graph) || graph_status=\$?",
        "if [ \"\$graph_status\" -eq 1 ]; then",
        "  $(gfa_cmd)",
        "  graph_status=0",
        "  graph_source=\$(find_metamdbg_graph) || graph_status=\$?",
        "fi",
        "if [ \"\$graph_status\" -ne 0 ]; then",
        "  echo \"metaMDBG produced no unique nonempty graph for k=$(graph_k)\" >&2",
        "  exit \"\$graph_status\"",
        "fi",
        "validate_gfa \"\$graph_source\"",
        "rm -f -- \"\$graph_alias\"",
        "ln -s -- \"\$(basename \"\$graph_source\")\" \"\$graph_alias\"",
        "test -L \"\$graph_alias\" && test -f \"\$graph_alias\" && test -s \"\$graph_alias\"",
        "validate_contigs \"\$contigs_gz\"",
        "validate_gfa \"\$graph_source\"",
        "if [ \"\$contract_exists\" -eq 1 ]; then",
        "  cmp -s -- \"\$expected_contract\" \"\$contract_marker\"",
        "else",
        "  if [ -e \"\$contract_marker\" ] || [ -L \"\$contract_marker\" ]; then",
        "    echo \"metaMDBG provenance marker appeared during locked execution\" >&2",
        "    exit 1",
        "  fi",
        "  cp -- \"\$expected_contract\" \"\$contract_new\"",
        "  chmod 600 \"\$contract_new\"",
        "  mv -n -- \"\$contract_new\" \"\$contract_marker\"",
        "  if [ -e \"\$contract_new\" ]; then",
        "    echo \"metaMDBG refused to overwrite a concurrent provenance marker\" >&2",
        "    exit 1",
        "  fi",
        "  cmp -s -- \"\$expected_contract\" \"\$contract_marker\"",
        "fi",
    ]
    return join(lines, "\n")
end

function _metamdbg_provenance(input_contract::NamedTuple)::NamedTuple
    toolchain = _metamdbg_expected_toolchain()
    return (;
        contract_signature = input_contract.signature,
        metamdbg_version = toolchain.metamdbg_version,
        environment_name = toolchain.environment_name,
        environment_spec_sha256 = toolchain.environment_spec_sha256,
    )
end

function _metamdbg_complete_result(
        outputs::NamedTuple,
        artifacts::NamedTuple,
        input_contract::NamedTuple,
)::NamedTuple
    return (;
        status = :complete,
        outdir = outputs.outdir,
        contigs = artifacts.contigs,
        graph = artifacts.graph,
        contract_marker = outputs.contract_marker,
        provenance = _metamdbg_provenance(input_contract),
    )
end

function _metamdbg_planned_result(
        status::Symbol,
        submission::Any,
        outputs::NamedTuple,
        input_contract::NamedTuple,
)::NamedTuple
    status in (:planned, :submitted) || throw(ArgumentError(
        "metaMDBG asynchronous status must be :planned or :submitted.",
    ))
    return (;
        status,
        submission,
        outdir = outputs.outdir,
        expected_artifacts = (;
            contigs = outputs.contigs_gz,
            graph = outputs.graph_alias,
            contract_marker = outputs.contract_marker,
        ),
        provenance = _metamdbg_provenance(input_contract),
    )
end

function _metamdbg_existing_artifacts(
        outputs::NamedTuple,
        graph_k::Int,
)::Union{Nothing, NamedTuple}
    existing_contigs = _normalize_metamdbg_contigs!(outputs)
    existing_graph = _normalize_metamdbg_graph!(outputs, graph_k)
    if existing_contigs !== nothing && existing_graph !== nothing
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
    end
    return nothing
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
        dependency_checker::Function = _ensure_metamdbg_installed,
        local_runner::Function = Base.run,
)::NamedTuple
    selected_input = _metamdbg_selected_input(hifi_reads, ont_reads)
    abundance_min > 0 || throw(ArgumentError("abundance_min must be positive."))
    threads > 0 || throw(ArgumentError("threads must be positive."))
    graph_k >= 0 || throw(ArgumentError("graph_k must be nonnegative."))
    outputs = _metamdbg_output_paths(outdir, graph_k)
    toolchain = _metamdbg_expected_toolchain()
    input_contract =
        _metamdbg_input_contract(selected_input, abundance_min, toolchain)
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
    asm_command = `$(Mycelia.CONDA_RUNNER) run --live-stream -n $(METAMDBG_ENV_NAME) $(asm_cmd_args)`
    gfa_command = `$(Mycelia.CONDA_RUNNER) run --live-stream -n $(METAMDBG_ENV_NAME) $(gfa_cmd_args)`
    resolved_executor = executor === nothing ?
                        Mycelia.LocalExecutor() :
                        Mycelia.resolve_executor(executor)

    if resolved_executor isa Mycelia.LocalExecutor
        return _with_metamdbg_output_lock(outputs.outdir) do
            has_contract = _prepare_metamdbg_output_root!(
                outputs,
                input_contract,
            )
            existing_artifacts = _metamdbg_existing_artifacts(outputs, graph_k)
            if existing_artifacts !== nothing
                has_contract || error(
                    "metaMDBG complete output is missing its provenance contract.",
                )
                _require_metamdbg_contract!(outputs, input_contract)
                return _metamdbg_complete_result(
                    outputs,
                    existing_artifacts,
                    input_contract,
                )
            end

            _require_expected_metamdbg_toolchain(dependency_checker())
            existing_contigs = _normalize_metamdbg_contigs!(outputs)
            if existing_contigs === nothing
                local_runner(asm_command)
                existing_contigs = _normalize_metamdbg_contigs!(outputs)
                existing_contigs === nothing && error(
                    "metaMDBG assembly produced no sequence-bearing contigs " *
                    "artifact in $(outputs.outdir).",
                )
            end
            existing_graph = _normalize_metamdbg_graph!(outputs, graph_k)
            if existing_graph === nothing
                local_runner(gfa_command)
            end
            artifacts = _require_metamdbg_artifacts!(outputs, graph_k)
            if has_contract
                _require_metamdbg_contract!(outputs, input_contract)
            else
                _write_metamdbg_contract!(outputs, input_contract)
            end
            _require_metamdbg_contract!(outputs, input_contract)
            return _metamdbg_complete_result(
                outputs,
                artifacts,
                input_contract,
            )
        end
    end

    complete_result = _with_metamdbg_output_lock(outputs.outdir) do
        has_contract = _prepare_metamdbg_output_root!(outputs, input_contract)
        existing_artifacts = _metamdbg_existing_artifacts(outputs, graph_k)
        if existing_artifacts === nothing
            return nothing
        end
        has_contract || error(
            "metaMDBG complete output is missing its provenance contract.",
        )
        _require_metamdbg_contract!(outputs, input_contract)
        return _metamdbg_complete_result(
            outputs,
            existing_artifacts,
            input_contract,
        )
    end
    complete_result !== nothing && return complete_result

    script = _metamdbg_executor_script(
        Mycelia.command_string(asm_command),
        Mycelia.command_string(gfa_command),
        outputs,
        graph_k,
        input_contract,
    )
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
        mail_user = mail_user,
    )
    if resolved_executor isa Mycelia.SlurmExecutor &&
       !resolved_executor.dry_run
        _require_expected_metamdbg_toolchain(dependency_checker())
    end
    submission = Mycelia.execute(job, resolved_executor)
    status = if resolved_executor isa Mycelia.CollectExecutor ||
                resolved_executor isa Mycelia.DryRunExecutor ||
                (resolved_executor isa Mycelia.SlurmExecutor &&
                 resolved_executor.dry_run)
        :planned
    else
        :submitted
    end
    return _metamdbg_planned_result(
        status,
        submission,
        outputs,
        input_contract,
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
alias, the contract marker, and exact toolchain provenance. Nonlocal execution
returns `status = :planned` for collected/dry-run jobs or `:submitted` for a
real submission, together with the backend `submission` result and
`expected_artifacts`; planned paths are not reported as completed artifacts.
The graph alias resolves to metaMDBG's validated dynamic
`assemblyGraph_k<k>_<length>bps.gfa` artifact.

# Details
- Uses metaMDBG's minimizer-space de Bruijn graphs for metagenomic assembly.
- Accepts both upstream `contigs.fasta.gz` and legacy plain `contigs.fasta`; a
  plain artifact is retained and normalized to `contigs.fasta.gz`.
- Reuses a complete contigs-plus-requested-graph result without provisioning or
  rerunning metaMDBG. If contigs exist but the requested graph does not, only
  graph generation runs. Complete and partial reuse require an exact durable
  contract marker for the normalized input paths, technology flag, input file
  size/modification metadata and SHA-256 content digests, `abundance_min`,
  metaMDBG 1.4, the spec-addressed environment name, and the bundled
  environment-spec checksum. Any nonempty uncontracted output root fails before
  provisioning. Nonlocal jobs revalidate every input digest after acquiring the
  runtime output lock and before artifact reuse or assembly.
- Installs metaMDBG exactly 1.4 from a checksum-verified, spec-hash-addressed
  environment and validates the installed package inventory before execution.
- Local and executor lifecycles validate sequence-bearing DNA FASTA and at
  least one unique sequence-bearing DNA GFA segment under an adjacent atomic
  output lock. The contract is stamped only after both semantic validations
  pass.
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
