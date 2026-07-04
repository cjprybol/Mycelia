# Iterative-Correction vs Naive Fixed-k Assembly Benchmark (DECISIVE harness)
#
# Purpose
# -------
# Head-to-head comparison of Mycelia's FOUNDATIONAL iterative read-correction
# pipeline (`mycelia_iterative_assemble`) against the NAIVE fixed-k assembly path
# (`Mycelia.Rhizomorph.assemble_genome`) on the SAME simulated cells. For each
# cell we simulate Illumina reads ONCE, then run each arm in a FRESH julia
# subprocess (for an honest whole-process peak RSS per arm) and score every
# arm's contigs with QUAST against the reference. The only variable across arms
# is whether the reads fed to the SHARED final assembler were (a) raw (naive) or
# (b) iteratively error-corrected first (iterative).
#
# The two arms (what actually differs)
# ------------------------------------
#   naive     : assemble_genome(RAW reads, SHARED_CONFIG)
#   iterative : reads = mycelia_iterative_assemble(fastq; max_k=...)  # Viterbi
#               read-polishing across a prime-k progression, returns POLISHED
#               reads (not contigs); then assemble_genome(POLISHED reads,
#               SHARED_CONFIG). The re-assemble uses the IDENTICAL config as the
#               naive arm, so any QUAST delta is attributable purely to the
#               iterative read-correction, not to a different assembler setting.
#
# Why `mycelia_iterative_assemble` is a read-polisher, not a contig-producer
# --------------------------------------------------------------------------
# `mycelia_iterative_assemble` (src/iterative-assembly.jl:153) runs, per prime k:
# build_qualmer_graph -> improve_read_set_likelihood (Viterbi correct_observations)
# -> write corrected reads -> next_prime_k. `finalize_iterative_assembly`
# (src/iterative-assembly.jl:1075) returns a Dict:
#     :final_assembly => Vector{String}  # the FINAL corrected READS (one string
#                                        #   per read) — NOT assembled contigs
#     :k_progression  => Vector{Int}     # the prime-k ladder visited
#     :metadata       => Dict with :final_k, :total_improvements,
#                        :final_fastq_file (path to the last corrected FASTQ),
#                        :iteration_history (per-k improvement counts)
# So to obtain CONTIGS for QUAST we assemble `metadata[:final_fastq_file]` with
# the shared config. We extract, from the returned metadata:
#   final_k          = res[:metadata][:final_k]
#   total_corrections= res[:metadata][:total_improvements]
#   per_k_corrections= sum(iter[:improvements_made]) over iteration_history[k]
#
# Shared assembler config (identical for BOTH arms)
# -------------------------------------------------
# k=MYCELIA_ITERVN_K (default 31, must be 1..64 per AssemblyConfig validation),
# use_quality_scores=false (the k-mer path — the ONLY path that honors
# memory_profile, and the path that makes bacterial-scale E. coli feasible),
# graph_mode=DoubleStrand, memory_profile=MYCELIA_ITERVN_PROFILE (default
# :lightweight — required so E. coli 10x does not OOM the way memory_profile=:full
# did; see benchmarking/ecoli_lightweight_crashtest.jl). Holding this config
# constant across arms is what makes the comparison honest.
#
# Momentum: NOT a third arm (documented, not faked)
# -------------------------------------------------
# The momentum fork resolver (`Mycelia.Rhizomorph.MomentumForkResolver`) is NOT
# wired into `assemble_genome` — grep of src/rhizomorph/assembly.jl finds zero
# references to it, and AssemblyConfig exposes no path-finding/strategy field for
# it. Its only entry point is benchmarking/08_momentum_fork_resolution_benchmark.jl,
# a standalone SYNTHETIC repeat-fork calibration that operates on ActiveReadState
# directly and produces NO contigs — it cannot be scored by QUAST and is
# incommensurable with an assemble->QUAST row. Adding it here would be a fake
# third arm, so per the harness spec it is OMITTED, and its proper harness is
# referenced instead.
#
# Cells
# -----
#   Lambda (NC_001416, ~48.5 kb) @ 30x  — always run
#   E. coli K-12 MG1655 (U00096.3, ~4.64 Mb) @ 10x — run iff MYCELIA_ITERVN_ECOLI
#     is truthy (default "true"). The iterative E. coli arm is the heaviest cell
#     (qualmer graphs at bacterial scale); set MYCELIA_ITERVN_ECOLI=false for a
#     Lambda-only smoke run.
# Both: seed 42, Illumina 150bp paired.
#
# External tools
# --------------
# art_illumina (simulation) and QUAST (scoring) are Conda-backed and set up on
# Lovelace, NOT the laptop — this script CANNOT run locally and is static-verified
# (parse + symbol checks) only. QUAST is gated behind MYCELIA_RUN_QUAST
# (default "true"). LD_LIBRARY_PATH / thread caps / JULIA_DEPOT_PATH are the
# caller's job (the run command sets them); the driver propagates the environment
# to each worker subprocess via `addenv`.
#
# Env knobs (all optional):
#   MYCELIA_RUN_QUAST        default "true"        — run QUAST scoring
#   MYCELIA_ITERVN_K         default "31"          — SHARED assemble k (1..64)
#   MYCELIA_ITERVN_MAX_K     default "31"          — iterative polish k-ceiling
#   MYCELIA_ITERVN_PROFILE   default "lightweight" — shared memory_profile
#   MYCELIA_ITERVN_ECOLI     default "true"        — include the E. coli cell
#   MYCELIA_QUAST_MIN_CONTIG default "100"         — QUAST --min-contig
#   MYCELIA_ITERVN_SEED      default "42"          — ART seed (both cells)
#
# Internal env knobs (set by the driver when it re-invokes this script as a
# worker — do NOT set these by hand):
#   MYCELIA_ITERVN_WORKER, IVN_ARM, IVN_CELL, IVN_FWD, IVN_REV, IVN_FASTQ,
#   IVN_K, IVN_MAX_K, IVN_PROFILE, IVN_CONTIGS_OUT, IVN_METRICS_OUT
#
# Usage (on Lovelace — invoke ONCE; the script self-manages a fresh subprocess
# per arm; run Lambda first (E. coli off) to smoke it, then enable E. coli):
#   julia --project=. benchmarking/iterative_vs_naive.jl

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import DataFrames
import CSV
import Dates

# ===========================================================================
# Shared helpers (used by both driver and worker)
# ===========================================================================

"""
    parse_quast_metric(report_tsv, metric) -> Union{Missing, Float64}

Return the numeric value of `metric` from a QUAST `report.tsv`, or `missing` if
the file/metric is absent or the value is non-numeric (QUAST prints "-" when a
metric could not be computed, e.g. NGA50 when nothing aligns).
"""
function parse_quast_metric(report_tsv::String, metric::String)::Union{Missing, Float64}
    isfile(report_tsv) || return missing
    result::Union{Missing, Float64} = missing
    open(report_tsv, "r") do io
        for line in eachline(io)
            fields = split(line, '\t')
            if length(fields) >= 2 && strip(fields[1]) == metric
                parsed = tryparse(Float64, strip(fields[2]))
                result = parsed === nothing ? missing : parsed
                break
            end
        end
    end
    return result
end

"""
    load_reads_from_paths(fwd, rev) -> Vector{FASTX.FASTQ.Record}

Load forward and (optional) reverse simulated reads into a single flat record
vector, matching how the k-mer assembler consumes reads (no pairing metadata).
"""
function load_reads_from_paths(fwd::String, rev::Union{String, Nothing})::Vector{FASTX.FASTQ.Record}
    records = FASTX.FASTQ.Record[]
    read_paths = String[fwd]
    if rev !== nothing && !isempty(rev)
        push!(read_paths, rev)
    end
    for p in read_paths
        reader = Mycelia.open_fastx(p)
        for record in reader
            push!(records, record)
        end
        close(reader)
    end
    return records
end

"""
    build_shared_config(records, k, profile) -> Mycelia.Rhizomorph.AssemblyConfig

Construct the SINGLE assembler configuration used identically by both arms. The
k-mer path (use_quality_scores=false) is chosen because it is the only path that
honors `memory_profile` and the path that makes bacterial-scale E. coli feasible.
"""
function build_shared_config(records::Vector{FASTX.FASTQ.Record}, k::Int,
        profile::Symbol)::Mycelia.Rhizomorph.AssemblyConfig
    seqtype = Mycelia.Rhizomorph._detect_sequence_type(records)
    return Mycelia.Rhizomorph.AssemblyConfig(;
        k = k,
        sequence_type = seqtype,
        graph_mode = Mycelia.Rhizomorph.DoubleStrand,
        use_quality_scores = false,
        verbose = false,
        memory_profile = profile
    )
end

"""
    write_contigs(contigs, path)

Write a `Vector{String}` of contigs to a simple FASTA at `path`.
"""
function write_contigs(contigs::Vector{String}, path::String)
    open(path, "w") do io
        for (i, contig) in enumerate(contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)
        end
    end
    return nothing
end

# ===========================================================================
# WORKER MODE — run ONE arm in this fresh process and emit metrics
# ===========================================================================

"""
    per_k_correction_string(metadata) -> String

Encode per-k read-correction counts from the iterative pipeline's
`iteration_history` as "k1:count1;k2:count2;..." for the metrics file.
"""
function per_k_correction_string(metadata::Dict)::String
    history = metadata[:iteration_history]
    parts = String[]
    for k in sort(collect(keys(history)))
        c = sum(iter[:improvements_made] for iter in history[k])
        push!(parts, "$(k):$(c)")
    end
    return join(parts, ";")
end

"""
    run_worker()

Executed when this script is re-invoked with MYCELIA_ITERVN_WORKER set. Runs the
requested arm in isolation (fresh process => honest whole-process peak RSS),
writes the contigs FASTA, and writes a `key=value` metrics file. All try/catch
result vars are function-local here (NOT top-level), so the Julia soft-scope trap
that bit the earlier crashtest does not apply and no `global` is needed.
"""
function run_worker()
    arm = ENV["IVN_ARM"]
    cell = ENV["IVN_CELL"]
    fwd = ENV["IVN_FWD"]
    rev = get(ENV, "IVN_REV", "")
    fastq = get(ENV, "IVN_FASTQ", "")
    k = parse(Int, ENV["IVN_K"])
    max_k = parse(Int, ENV["IVN_MAX_K"])
    profile = Symbol(ENV["IVN_PROFILE"])
    contigs_out = ENV["IVN_CONTIGS_OUT"]
    metrics_out = ENV["IVN_METRICS_OUT"]

    # Metrics accumulators (function-local — no soft-scope trap).
    completed = false
    n_contigs = -1
    assemble_runtime_s = -1.0
    polish_runtime_s = -1.0
    final_k = -1
    total_corrections = -1
    per_k = ""

    GC.gc(true)

    if arm == "naive"
        println("  [worker:$(cell):naive] loading raw reads ...")
        records = load_reads_from_paths(fwd, isempty(rev) ? nothing : rev)
        isempty(records) && error("[worker:$(cell):naive] no reads loaded from $(fwd)")
        config = build_shared_config(records, k, profile)
        t0 = time()
        result = Mycelia.Rhizomorph.assemble_genome(records, config)
        assemble_runtime_s = round(time() - t0; digits = 3)
        n_contigs = length(result.contigs)
        write_contigs(result.contigs, contigs_out)
        completed = true
        println("  [worker:$(cell):naive] -> $(n_contigs) contigs in $(assemble_runtime_s)s")
    elseif arm == "iterative"
        # Phase A: iterative Viterbi read-polishing across the prime-k ladder.
        isempty(fastq) && error("[worker:$(cell):iterative] IVN_FASTQ not set")
        polish_dir = joinpath(dirname(contigs_out), "iterative_out")
        println("  [worker:$(cell):iterative] polishing reads (max_k=$(max_k)) ...")
        tp = time()
        # td-ve02 speed levers (zero accuracy change): (1) parallelize read
        # correction — embarrassingly parallel, ~Nthreads x; (2) cap iterations
        # per k — convergence is usually < 3, the default 10 is wasteful.
        max_iter = parse(Int, get(ENV, "MYCELIA_ITERVN_MAX_ITER", "3"))
        iter_verbose = get(ENV, "MYCELIA_ITERVN_VERBOSE", "false") == "true"
        res = Mycelia.mycelia_iterative_assemble(fastq;
            max_k = max_k,
            max_iterations_per_k = max_iter,
            enable_parallel = true,
            output_dir = polish_dir,
            verbose = iter_verbose,
            enable_checkpointing = false)
        polish_runtime_s = round(time() - tp; digits = 3)
        final_k = res[:metadata][:final_k]
        total_corrections = res[:metadata][:total_improvements]
        per_k = per_k_correction_string(res[:metadata])
        polished_fastq = res[:metadata][:final_fastq_file]
        println("  [worker:$(cell):iterative] polish done: final_k=$(final_k), " *
                "corrections=$(total_corrections), per_k=[$(per_k)] in $(polish_runtime_s)s")
        # Surface the per-k, per-iteration wall-clock already recorded in
        # iteration_history so the profile shows WHERE iterative time goes.
        ihist = get(res[:metadata], :iteration_history, Dict())
        for kk in sort(collect(keys(ihist)))
            for it in ihist[kk]
                println("    [profile] k=$(kk) iter=$(get(it, :iteration, 0)) " *
                        "improvements=$(get(it, :improvements_made, 0)) " *
                        "runtime_s=$(round(Float64(get(it, :runtime_seconds, 0.0)); digits=2)) " *
                        "kmers=$(get(it, :memory_kmers, 0))")
            end
        end
        flush(stdout)

        # Phase B: assemble the POLISHED reads with the SHARED config.
        records = load_reads_from_paths(polished_fastq, nothing)
        isempty(records) && error("[worker:$(cell):iterative] no polished reads at $(polished_fastq)")
        config = build_shared_config(records, k, profile)
        t0 = time()
        result = Mycelia.Rhizomorph.assemble_genome(records, config)
        assemble_runtime_s = round(time() - t0; digits = 3)
        n_contigs = length(result.contigs)
        write_contigs(result.contigs, contigs_out)
        completed = true
        println("  [worker:$(cell):iterative] -> $(n_contigs) contigs in $(assemble_runtime_s)s " *
                "(polish $(polish_runtime_s)s)")
    else
        error("Unknown IVN_ARM: $(arm)")
    end

    peak_rss_kb = Int(Sys.maxrss() ÷ 1024)   # honest whole-process peak for this arm

    open(metrics_out, "w") do io
        println(io, "arm=$(arm)")
        println(io, "cell=$(cell)")
        println(io, "completed=$(completed)")
        println(io, "n_contigs=$(n_contigs)")
        println(io, "assemble_runtime_s=$(assemble_runtime_s)")
        println(io, "polish_runtime_s=$(polish_runtime_s)")
        println(io, "peak_rss_kb=$(peak_rss_kb)")
        println(io, "assemble_k=$(k)")
        println(io, "final_k=$(final_k)")
        println(io, "total_corrections=$(total_corrections)")
        println(io, "per_k_corrections=$(per_k)")
    end
    return nothing
end

# ===========================================================================
# DRIVER MODE — simulate once per cell, spawn a worker per arm, QUAST + aggregate
# ===========================================================================

"""
    parse_metrics_file(path) -> Dict{String, String}

Parse a worker `key=value` metrics file into a Dict.
"""
function parse_metrics_file(path::String)::Dict{String, String}
    d = Dict{String, String}()
    isfile(path) || return d
    for line in eachline(path)
        isempty(strip(line)) && continue
        kv = split(line, '='; limit = 2)
        length(kv) == 2 && (d[strip(kv[1])] = strip(kv[2]))
    end
    return d
end

"""
    Cell — one benchmark organism/coverage cell.
"""
struct Cell
    name::String
    accession::String
    coverage::Int
end

function run_driver()
    seed = parse(Int, get(ENV, "MYCELIA_ITERVN_SEED", "42"))
    k = parse(Int, get(ENV, "MYCELIA_ITERVN_K", "31"))
    max_k = parse(Int, get(ENV, "MYCELIA_ITERVN_MAX_K", "31"))
    profile = get(ENV, "MYCELIA_ITERVN_PROFILE", "lightweight")
    run_quast = get(ENV, "MYCELIA_RUN_QUAST", "true") in ["1", "true", "yes"]
    quast_min_contig = parse(Int, get(ENV, "MYCELIA_QUAST_MIN_CONTIG", "100"))
    include_ecoli = get(ENV, "MYCELIA_ITERVN_ECOLI", "true") in ["1", "true", "yes"]

    project_dir = dirname(@__DIR__)   # worktree root (holds Project.toml)
    script_path = @__FILE__

    cells = Cell[Cell("Lambda", "NC_001416", 30)]
    if include_ecoli
        push!(cells, Cell("Ecoli_K12_MG1655", "U00096.3", 10))
    end
    arms = ["naive", "iterative"]

    println("=== Iterative-vs-Naive Assembly Benchmark ===")
    println("Start time (UTC): $(Dates.now(Dates.UTC))")
    println("Cells: $(join(["$(c.name)@$(c.coverage)x" for c in cells], ", "))")
    println("Arms: $(join(arms, ", ")) (momentum OMITTED — not an assemble_genome strategy; see header)")
    println("Shared assemble: k=$(k), use_quality_scores=false, DoubleStrand, memory_profile=:$(profile)")
    println("Iterative polish ceiling: max_k=$(max_k), seed=$(seed)")
    println("QUAST scoring: $(run_quast ? "ENABLED (min_contig=$(quast_min_contig))" : "DISABLED")")

    benchmark_dir = mktempdir(prefix = "iterative_vs_naive_")
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    println("Working directory: $(benchmark_dir)")

    results = DataFrames.DataFrame(
        cell = String[],
        accession = String[],
        coverage = Int[],
        arm = String[],
        completed = Bool[],
        n_contigs = Union{Missing, Int}[],
        NGA50 = Union{Missing, Float64}[],
        num_misassemblies = Union{Missing, Float64}[],
        dup_ratio = Union{Missing, Float64}[],
        genome_fraction_quast = Union{Missing, Float64}[],
        assemble_k = Union{Missing, Int}[],
        assemble_runtime_s = Union{Missing, Float64}[],
        polish_runtime_s = Union{Missing, Float64}[],
        peak_rss_kb = Union{Missing, Int}[],
        final_k = Union{Missing, Int}[],
        total_corrections = Union{Missing, Int}[],
        per_k_corrections = Union{Missing, String}[],
        quast_runtime_s = Union{Missing, Float64}[],
        worker_ok = Bool[]
    )

    for cell in cells
        println("\n" * "="^70)
        println("CELL: $(cell.name) ($(cell.accession)) @ $(cell.coverage)x")
        println("="^70)

        # --- Phase 1: download reference ---
        println("\n--- Phase 1: download $(cell.name) reference ---")
        ref_dir = joinpath(benchmark_dir, "$(cell.name)_reference")
        mkpath(ref_dir)
        reference_path = Mycelia.download_genome_by_accession(
            accession = cell.accession, outdir = ref_dir, compressed = false)
        if !isfile(reference_path) || filesize(reference_path) == 0
            error("Failed to download $(cell.name) reference ($(cell.accession)) to $(ref_dir)")
        end
        println("Reference: $(reference_path) ($(filesize(reference_path)) bytes)")

        # --- Phase 2: simulate reads ONCE (shared across both arms) ---
        println("\n--- Phase 2: simulate reads once (coverage=$(cell.coverage)x seed=$(seed)) ---")
        reads_dir = joinpath(benchmark_dir, "$(cell.name)_reads")
        mkpath(reads_dir)
        outbase = joinpath(reads_dir, "$(cell.name)_cov$(cell.coverage)x_seed$(seed)")
        sim = Mycelia.simulate_illumina_reads(
            fasta = reference_path,
            coverage = cell.coverage,
            rndSeed = seed,
            read_length = 150,
            paired = true,
            quiet = true,
            outbase = outbase)
        fwd = sim.forward_reads
        rev = sim.reverse_reads === nothing ? "" : sim.reverse_reads
        println("Forward reads: $(fwd)")
        println("Reverse reads: $(isempty(rev) ? "(none)" : rev)")

        # The iterative arm needs a SINGLE FASTQ path. Combine fwd+rev into one
        # file once, so both mates feed the polishing pipeline.
        combined_fastq = joinpath(reads_dir, "$(cell.name)_combined.fastq")
        combined_records = load_reads_from_paths(fwd, isempty(rev) ? nothing : rev)
        Mycelia.write_fastq(records = combined_records, filename = combined_fastq)
        println("Combined FASTQ (for iterative arm): $(combined_fastq) " *
                "($(length(combined_records)) reads)")

        # --- Phase 3: per-arm fresh-subprocess run + QUAST ---
        for arm in arms
            println("\n  --- arm: $(arm) ($(cell.name)) ---")
            arm_dir = joinpath(benchmark_dir, "$(cell.name)_arm_$(arm)")
            mkpath(arm_dir)
            contigs_out = joinpath(arm_dir, "contigs.fasta")
            metrics_out = joinpath(arm_dir, "metrics.txt")

            worker_ok = true
            worker_cmd = `$(Base.julia_cmd()) --project=$(project_dir) --threads=$(Threads.nthreads()) $(script_path)`
            worker_cmd = addenv(worker_cmd,
                "MYCELIA_ITERVN_WORKER" => "1",
                "IVN_ARM" => arm,
                "IVN_CELL" => cell.name,
                "IVN_FWD" => fwd,
                "IVN_REV" => rev,
                "IVN_FASTQ" => combined_fastq,
                "IVN_K" => string(k),
                "IVN_MAX_K" => string(max_k),
                "IVN_PROFILE" => profile,
                "IVN_CONTIGS_OUT" => contigs_out,
                "IVN_METRICS_OUT" => metrics_out)
            try
                run(worker_cmd)
            catch e
                worker_ok = false
                @warn "Worker process failed for cell=$(cell.name) arm=$(arm) — recording as failure" exception = (e,)
            end

            metrics = parse_metrics_file(metrics_out)
            getint(key) = haskey(metrics, key) ? tryparse(Int, metrics[key]) : missing
            getfloat(key) = haskey(metrics, key) ? tryparse(Float64, metrics[key]) : missing
            # Sentinels (-1 / "") from the worker mean "not applicable for this arm".
            nz(x) = (x === nothing || x === missing || x == -1) ? missing : x
            completed = get(metrics, "completed", "false") == "true"
            n_contigs = nz(getint("n_contigs"))
            assemble_k = nz(getint("assemble_k"))
            assemble_runtime_s = nz(getfloat("assemble_runtime_s"))
            polish_runtime_s = nz(getfloat("polish_runtime_s"))
            peak_rss_kb = nz(getint("peak_rss_kb"))
            final_k = nz(getint("final_k"))
            total_corrections = nz(getint("total_corrections"))
            per_k_raw = get(metrics, "per_k_corrections", "")
            per_k = isempty(per_k_raw) ? missing : per_k_raw

            # -- QUAST aligned-truth scoring on the produced contigs --
            nga50 = missing
            num_misassemblies = missing
            dup_ratio = missing
            genome_fraction = missing
            quast_runtime_s = missing
            have_contigs = (n_contigs !== missing && n_contigs > 0 &&
                            isfile(contigs_out) && filesize(contigs_out) > 0)
            if run_quast && have_contigs
                quast_dir = joinpath(arm_dir, "quast")
                tq = time()
                try
                    # Redirect QUAST's ~700 lines of verbose stdout/stderr to a
                    # per-arm tool logfile so only the parsed metrics reach the
                    # main log — the decisive numbers were getting buried (td-cluf).
                    quast_tool_log = joinpath(arm_dir, "quast_tool.log")
                    open(quast_tool_log, "w") do qio
                        redirect_stdout(qio) do
                            redirect_stderr(qio) do
                                Mycelia.run_quast(
                                    contigs_out;
                                    outdir = quast_dir,
                                    reference = reference_path,
                                    min_contig = quast_min_contig)
                            end
                        end
                    end
                    quast_runtime_s = round(time() - tq; digits = 3)
                    report_tsv = joinpath(quast_dir, "report.tsv")
                    nga50 = parse_quast_metric(report_tsv, "NGA50")
                    num_misassemblies = parse_quast_metric(report_tsv, "# misassemblies")
                    dup_ratio = parse_quast_metric(report_tsv, "Duplication ratio")
                    genome_fraction = parse_quast_metric(report_tsv, "Genome fraction (%)")
                catch e
                    @warn "QUAST failed for cell=$(cell.name) arm=$(arm) — leaving accuracy metrics missing" exception = (e,)
                end
            elseif run_quast && !have_contigs
                println("    (no contigs to score — arm likely FAILED; leaving QUAST metrics missing)")
            end

            push!(results,
                (
                    cell = cell.name,
                    accession = cell.accession,
                    coverage = cell.coverage,
                    arm = arm,
                    completed = completed,
                    n_contigs = n_contigs,
                    NGA50 = nga50,
                    num_misassemblies = num_misassemblies,
                    dup_ratio = dup_ratio,
                    genome_fraction_quast = genome_fraction,
                    assemble_k = assemble_k,
                    assemble_runtime_s = assemble_runtime_s,
                    polish_runtime_s = polish_runtime_s,
                    peak_rss_kb = peak_rss_kb,
                    final_k = final_k,
                    total_corrections = total_corrections,
                    per_k_corrections = per_k,
                    quast_runtime_s = quast_runtime_s,
                    worker_ok = worker_ok
                ))
        end
    end

    # --- Phase 4: write CSV (UTC, yyyymmdd_HHMMSS — no 'Z' directive) ---
    timestamp = Dates.format(Dates.now(Dates.UTC), "yyyymmdd_HHMMSS")
    csv_path = joinpath(results_dir, "iterative_vs_naive_$(timestamp).csv")
    CSV.write(csv_path, results)
    println("\nResults written to: $(csv_path)")

    # --- Phase 5: comparison table + verdict ---
    print_comparison_table(results)
    print_verdict(results)

    println("\n=== Iterative-vs-Naive complete (UTC): $(Dates.now(Dates.UTC)) ===")
    return nothing
end

"""
    fmt(x) -> String

Compact formatter for possibly-missing metric values.
"""
function fmt(x)::String
    ismissing(x) && return "NA"
    x isa AbstractFloat && return string(round(x; digits = 3))
    return string(x)
end

"""
    print_comparison_table(results)

Print a per-cell arm x {contigs, GF, NGA50, misasm, dup, runtime, peak_rss,
final_k, corrections} comparison table.
"""
function print_comparison_table(results::DataFrames.DataFrame)
    println("\n--- Comparison table (cell x arm) ---")
    header = rpad("cell", 20) * rpad("arm", 11) * rpad("contigs", 9) *
             rpad("GF%", 9) * rpad("NGA50", 10) * rpad("misasm", 8) *
             rpad("dup", 7) * rpad("asm_s", 9) * rpad("polish_s", 10) *
             rpad("peak_rss_kb", 13) * rpad("final_k", 9) * rpad("corrections", 12)
    println(header)
    println(repeat("-", length(header)))
    for row in eachrow(results)
        line = rpad(row.cell, 20) *
               rpad(row.arm, 11) *
               rpad(fmt(row.n_contigs), 9) *
               rpad(fmt(row.genome_fraction_quast), 9) *
               rpad(fmt(row.NGA50), 10) *
               rpad(fmt(row.num_misassemblies), 8) *
               rpad(fmt(row.dup_ratio), 7) *
               rpad(fmt(row.assemble_runtime_s), 9) *
               rpad(fmt(row.polish_runtime_s), 10) *
               rpad(fmt(row.peak_rss_kb), 13) *
               rpad(fmt(row.final_k), 9) *
               rpad(fmt(row.total_corrections), 12)
        println(line)
    end
end

"""
    print_verdict(results)

For each cell, compare the `iterative` arm against the `naive` arm and print the
three DECISIVE tests the harness exists to answer:
  (a) lower peak RSS (esp. at E. coli scale),
  (b) higher NGA50 / fewer contigs,
  (c) preserved-or-better genome fraction / misassemblies.
"""
function print_verdict(results::DataFrames.DataFrame)
    println("\n--- Verdict (iterative vs naive, per cell) ---")
    for cell in unique(results.cell)
        sub = results[results.cell .== cell, :]
        nrows = sub[sub.arm .== "naive", :]
        irows = sub[sub.arm .== "iterative", :]
        if DataFrames.nrow(nrows) == 0 || DataFrames.nrow(irows) == 0
            println("  [$(cell)] missing an arm — cannot compare.")
            continue
        end
        n = nrows[1, :]
        it = irows[1, :]

        println("  [$(cell)]")
        if !n.completed || !it.completed
            println("    completion: naive=$(n.completed), iterative=$(it.completed) " *
                    "(a FAILED arm blocks a clean comparison below)")
        end

        # (a) peak RSS
        if !ismissing(n.peak_rss_kb) && !ismissing(it.peak_rss_kb)
            delta = it.peak_rss_kb - n.peak_rss_kb
            verdict = delta < 0 ? "LOWER (iterative wins)" : delta == 0 ? "equal" : "HIGHER (naive wins)"
            println("    (a) peak RSS: iterative $(it.peak_rss_kb) vs naive $(n.peak_rss_kb) kb " *
                    "(delta $(delta) kb) -> $(verdict)")
        else
            println("    (a) peak RSS: NA (an arm did not report)")
        end

        # (b) NGA50 (higher better) + contigs (fewer better)
        if !ismissing(it.NGA50) && !ismissing(n.NGA50)
            v = it.NGA50 > n.NGA50 ? "HIGHER (iterative wins)" : it.NGA50 == n.NGA50 ? "equal" : "LOWER (naive wins)"
            println("    (b) NGA50: iterative $(fmt(it.NGA50)) vs naive $(fmt(n.NGA50)) -> $(v)")
        else
            println("    (b) NGA50: NA (an arm did not align/score)")
        end
        if !ismissing(it.n_contigs) && !ismissing(n.n_contigs)
            v = it.n_contigs < n.n_contigs ? "FEWER (iterative wins)" : it.n_contigs == n.n_contigs ? "equal" : "MORE (naive wins)"
            println("        contigs: iterative $(it.n_contigs) vs naive $(n.n_contigs) -> $(v)")
        end

        # (c) genome fraction (higher better) + misassemblies (fewer better)
        if !ismissing(it.genome_fraction_quast) && !ismissing(n.genome_fraction_quast)
            v = it.genome_fraction_quast >= n.genome_fraction_quast ? "PRESERVED/BETTER (iterative >= naive)" : "WORSE (iterative < naive)"
            println("    (c) genome fraction: iterative $(fmt(it.genome_fraction_quast))% vs naive " *
                    "$(fmt(n.genome_fraction_quast))% -> $(v)")
        else
            println("    (c) genome fraction: NA (an arm did not align/score)")
        end
        if !ismissing(it.num_misassemblies) && !ismissing(n.num_misassemblies)
            v = it.num_misassemblies <= n.num_misassemblies ? "PRESERVED/BETTER (iterative <= naive)" : "WORSE (iterative > naive)"
            println("        misassemblies: iterative $(fmt(it.num_misassemblies)) vs naive " *
                    "$(fmt(n.num_misassemblies)) -> $(v)")
        end

        if !ismissing(it.final_k) || !ismissing(it.total_corrections)
            println("    iterative internals: final_k=$(fmt(it.final_k)), " *
                    "total_corrections=$(fmt(it.total_corrections)), per_k=[$(fmt(it.per_k_corrections))]")
        end
    end

    println("\nNote: accuracy uses QUAST aligned-truth metrics; peak RSS is the honest")
    println("whole-process Sys.maxrss() of each arm's fresh subprocess (iterative includes")
    println("both the read-polishing and the re-assembly, as that is the arm's true cost).")
    return nothing
end

# ===========================================================================
# Entry point: worker vs driver
# ===========================================================================

if get(ENV, "MYCELIA_ITERVN_WORKER", "") == "1"
    run_worker()
else
    run_driver()
end
