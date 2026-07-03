# Rhizomorph Assembly-Efficiency Mode-Comparison Benchmark
#
# Purpose
# -------
# Empirically measures the assembly-efficiency modes introduced in PR #340
# (cp/rhizomorph-efficiency-modes) so we can build a data-informed
# accuracy-vs-efficiency recommendation table. For ONE fixed Lambda cell
# (30x, seed 42) we simulate Illumina reads ONCE, then assemble those SAME
# reads under each efficiency-mode configuration and score every assembly with
# QUAST. The only variable across rows is the AssemblyConfig mode; reads,
# reference, k, and QUAST settings are held constant.
#
# Why the NON-quality k-mer path (use_quality_scores = false)
# -----------------------------------------------------------
# The three efficiency modes under test — memory_profile, dedup_revcomp,
# compact_unitigs — are threaded ONLY into `_assemble_kmer_graph`
# (src/rhizomorph/assembly.jl, the use_quality_scores = false arm). The
# quality-aware qualmer path (`_assemble_qualmer_graph`, auto-selected for FASTQ
# input) does NOT honor any of them. To make the modes actually take effect —
# and to keep "the only variable is the mode" honest — every row here builds an
# AssemblyConfig with use_quality_scores = false, routing all modes through the
# same k-mer path. This deliberately differs from track_a_lambda_pilot.jl, which
# uses the qualmer path; absolute accuracy numbers are therefore not directly
# comparable to the pilot, but they ARE directly comparable to each other, which
# is the point of a mode comparison.
#
# Memory measurement (chosen approach + caveat)
# ---------------------------------------------
# `Sys.maxrss()` is process-cumulative (monotonic peak), so running several
# assembles sequentially in one process cannot reveal per-mode peak RSS. We
# therefore run each mode's assemble in a FRESH julia subprocess (the driver
# re-invokes THIS script in worker mode, once per mode) and report that
# process's whole-process peak `Sys.maxrss()` as `assemble_peak_rss_kb`. This is
# the most honest peak-memory figure available.
#   Caveat: the reported peak includes the CONSTANT cost of loading Julia +
#   Mycelia + the reads into the worker (hundreds of MB), which is the same for
#   every mode, so cross-mode DIFFERENCES are meaningful even though the absolute
#   value is inflated by fixed load. For transparency each row also records:
#     - baseline_rss_kb        : peak maxrss immediately BEFORE the assemble call
#     - assemble_rss_delta_kb  : peak_rss - baseline_rss (assembly-attributable
#                                peak growth)
#     - gc_live_delta_mb       : Base.gc_live_bytes() delta around the assemble
#                                (retained live-heap proxy, e.g. the graph)
#
# External tools
# --------------
# Requires Conda-backed external tools (art_illumina for simulation, QUAST for
# scoring) — set up on Lovelace, NOT the laptop. This script CANNOT be run
# locally; it is static-verified (parse + field-name checks) only. QUAST is
# gated behind MYCELIA_RUN_QUAST (default "true"). LD_LIBRARY_PATH / thread caps
# / JULIA_DEPOT_PATH are the caller's job (the run command sets them); the driver
# propagates the current environment to each worker subprocess via `addenv`.
#
# Env knobs (all optional):
#   MYCELIA_RUN_QUAST        default "true"  — run QUAST scoring
#   MYCELIA_MODE_K           default "31"    — k-mer size for assembly
#   MYCELIA_QUAST_MIN_CONTIG default "100"   — QUAST --min-contig (fixed across
#                                              rows so it adds no cross-row noise)
#   MYCELIA_MODE_COVERAGE    default "30"    — single coverage tier
#   MYCELIA_MODE_SEED        default "42"    — single ART seed
#   MYCELIA_MODE_COMPACT     default "true"  — include the optional compact_unitigs
#                                              no-op-confirmation row
#
# Internal env knobs (set by the driver when it re-invokes this script as a
# worker — do NOT set these by hand):
#   MYCELIA_MODE_WORKER, MMC_MODE_LABEL, MMC_GRAPH_MODE, MMC_DEDUP,
#   MMC_MEMPROFILE, MMC_COMPACT, MMC_K, MMC_FWD, MMC_REV, MMC_CONTIGS_OUT,
#   MMC_METRICS_OUT
#
# Usage (on Lovelace — invoke ONCE; the script self-manages a fresh subprocess
# per mode):
#   julia --project=. benchmarking/mode_comparison.jl

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import DataFrames
import CSV
import Dates
import BioSequences

# ===========================================================================
# Shared helpers (used by both driver and worker)
# ===========================================================================

"""
    parse_quast_metric(report_tsv, metric) -> Union{Missing, Float64}

Return the numeric value of `metric` from a QUAST `report.tsv`, or `missing` if
the file/metric is absent or the value is non-numeric (QUAST prints "-" when a
metric could not be computed, e.g. NGA50 when nothing aligns — the empirical
signature we expect for the BROKEN canonical mode).
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
    graph_mode_from_string(s) -> Mycelia.Rhizomorph.GraphMode

Map a lowercase mode string to the Rhizomorph GraphMode enum value.
"""
function graph_mode_from_string(s::String)::Mycelia.Rhizomorph.GraphMode
    if s == "doublestrand"
        return Mycelia.Rhizomorph.DoubleStrand
    elseif s == "canonical"
        return Mycelia.Rhizomorph.Canonical
    elseif s == "singlestrand"
        return Mycelia.Rhizomorph.SingleStrand
    else
        error("Unknown graph_mode string: $(s)")
    end
end

# ===========================================================================
# WORKER MODE — assemble ONE mode in this fresh process and emit metrics
# ===========================================================================

"""
    run_worker()

Executed when this script is re-invoked with MYCELIA_MODE_WORKER set. Loads the
pre-simulated reads, builds the AssemblyConfig for the requested mode, runs the
assemble in isolation, writes the contigs FASTA, and writes a `key=value`
metrics file. Because this is a fresh process, `Sys.maxrss()` here is this mode's
whole-process peak RSS.
"""
function run_worker()
    mode_label = ENV["MMC_MODE_LABEL"]
    graph_mode = graph_mode_from_string(ENV["MMC_GRAPH_MODE"])
    dedup = ENV["MMC_DEDUP"] == "true"
    memprofile = Symbol(ENV["MMC_MEMPROFILE"])
    compact = ENV["MMC_COMPACT"] == "true"
    k = parse(Int, ENV["MMC_K"])
    fwd = ENV["MMC_FWD"]
    rev = get(ENV, "MMC_REV", "")
    contigs_out = ENV["MMC_CONTIGS_OUT"]
    metrics_out = ENV["MMC_METRICS_OUT"]

    println("  [worker:$(mode_label)] loading reads ...")
    records = load_reads_from_paths(fwd, isempty(rev) ? nothing : rev)
    if isempty(records)
        error("[worker:$(mode_label)] no reads loaded from $(fwd)")
    end

    # Build the exact AssemblyConfig for this mode. use_quality_scores = false
    # forces the k-mer path where the efficiency modes are wired (see header).
    seqtype = Mycelia.Rhizomorph._detect_sequence_type(records)
    config = Mycelia.Rhizomorph.AssemblyConfig(;
        k = k,
        sequence_type = seqtype,
        graph_mode = graph_mode,
        use_quality_scores = false,
        verbose = false,
        dedup_revcomp = dedup,
        compact_unitigs = compact,
        memory_profile = memprofile
    )

    # Establish a clean pre-assemble memory baseline, then measure around the call.
    GC.gc(true)
    baseline_rss_kb = Int(Sys.maxrss() ÷ 1024)
    live_before = Base.gc_live_bytes()

    t0 = time()
    result = Mycelia.Rhizomorph.assemble_genome(records, config)
    runtime = time() - t0

    live_after = Base.gc_live_bytes()
    peak_rss_kb = Int(Sys.maxrss() ÷ 1024)   # whole-process peak for this mode

    n_contigs = length(result.contigs)
    open(contigs_out, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)  # contigs are String
        end
    end

    gc_live_delta_mb = round((live_after - live_before) / (1024 * 1024); digits = 3)
    assemble_rss_delta_kb = peak_rss_kb - baseline_rss_kb

    open(metrics_out, "w") do io
        println(io, "mode_label=$(mode_label)")
        println(io, "n_contigs=$(n_contigs)")
        println(io, "assemble_runtime_s=$(round(runtime; digits = 3))")
        println(io, "assemble_peak_rss_kb=$(peak_rss_kb)")
        println(io, "baseline_rss_kb=$(baseline_rss_kb)")
        println(io, "assemble_rss_delta_kb=$(assemble_rss_delta_kb)")
        println(io, "gc_live_delta_mb=$(gc_live_delta_mb)")
    end
    println("  [worker:$(mode_label)] -> $(n_contigs) contigs, " *
            "$(round(runtime; digits = 2))s, peak_rss=$(peak_rss_kb)kb " *
            "(delta $(assemble_rss_delta_kb)kb, live +$(gc_live_delta_mb)MB)")
    return nothing
end

# ===========================================================================
# DRIVER MODE — simulate once, spawn a worker per mode, QUAST + aggregate
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

function run_driver()
    organism = "Lambda"
    accession = "NC_001416"
    coverage = parse(Int, get(ENV, "MYCELIA_MODE_COVERAGE", "30"))
    seed = parse(Int, get(ENV, "MYCELIA_MODE_SEED", "42"))
    k = parse(Int, get(ENV, "MYCELIA_MODE_K", "31"))
    run_quast = get(ENV, "MYCELIA_RUN_QUAST", "true") in ["1", "true", "yes"]
    quast_min_contig = parse(Int, get(ENV, "MYCELIA_QUAST_MIN_CONTIG", "100"))
    include_compact = get(ENV, "MYCELIA_MODE_COMPACT", "true") in ["1", "true", "yes"]

    project_dir = dirname(@__DIR__)   # worktree root (holds Project.toml)
    script_path = @__FILE__

    println("=== Rhizomorph Mode-Comparison Benchmark ===")
    println("Start time (UTC): $(Dates.now(Dates.UTC))")
    println("Fixed cell: $(organism) ($(accession)) illumina, coverage=$(coverage)x seed=$(seed), k=$(k)")
    println("Assembly path: NON-quality k-mer path (use_quality_scores=false) — required for modes to take effect")
    println("Memory: fresh subprocess per mode, whole-process peak Sys.maxrss()")
    println("QUAST scoring: $(run_quast ? "ENABLED (min_contig=$(quast_min_contig))" : "DISABLED")")

    # --- Mode table: only the AssemblyConfig fields vary across rows ---
    modes = [
        (label = "baseline", graph_mode = "doublestrand", dedup = false, memprofile = "full", compact = false),
        (label = "dedup", graph_mode = "doublestrand", dedup = true, memprofile = "full", compact = false),
        (label = "lightweight", graph_mode = "doublestrand", dedup = false, memprofile = "lightweight", compact = false),
        (label = "dedup+lightweight", graph_mode = "doublestrand", dedup = true, memprofile = "lightweight", compact = false),
        (label = "canonical", graph_mode = "canonical", dedup = false, memprofile = "full", compact = false)
    ]
    if include_compact
        push!(modes,
            (label = "compact_unitigs", graph_mode = "doublestrand", dedup = false, memprofile = "full", compact = true))
    end

    benchmark_dir = mktempdir(prefix = "mode_comparison_")
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    println("Working directory: $(benchmark_dir)")

    # --- Phase 1: download reference (once) ---
    println("\n--- Phase 1: download Lambda reference ---")
    ref_dir = joinpath(benchmark_dir, "reference")
    mkpath(ref_dir)
    reference_path = Mycelia.download_genome_by_accession(
        accession = accession, outdir = ref_dir, compressed = false)
    if !isfile(reference_path) || filesize(reference_path) == 0
        error("Failed to download Lambda reference ($(accession)) to $(ref_dir)")
    end
    println("Reference: $(reference_path) ($(filesize(reference_path)) bytes)")

    # --- Phase 2: simulate reads ONCE (shared across all modes) ---
    println("\n--- Phase 2: simulate reads once (coverage=$(coverage)x seed=$(seed)) ---")
    reads_dir = joinpath(benchmark_dir, "reads")
    mkpath(reads_dir)
    outbase = joinpath(reads_dir, "$(organism)_cov$(coverage)x_seed$(seed)")
    sim = Mycelia.simulate_illumina_reads(
        fasta = reference_path,
        coverage = coverage,
        rndSeed = seed,
        read_length = 150,
        paired = true,
        quiet = true,
        outbase = outbase)
    fwd = sim.forward_reads
    rev = sim.reverse_reads === nothing ? "" : sim.reverse_reads
    println("Forward reads: $(fwd)")
    println("Reverse reads: $(isempty(rev) ? "(none)" : rev)")

    # --- Phase 3: per-mode fresh-subprocess assemble + QUAST ---
    results = DataFrames.DataFrame(
        organism = String[],
        mode_label = String[],
        graph_mode = String[],
        dedup = Bool[],
        memory_profile = String[],
        compact = Bool[],
        n_contigs = Union{Missing, Int}[],
        NGA50 = Union{Missing, Float64}[],
        num_misassemblies = Union{Missing, Float64}[],
        dup_ratio = Union{Missing, Float64}[],
        genome_fraction_quast = Union{Missing, Float64}[],
        assemble_runtime_s = Union{Missing, Float64}[],
        assemble_peak_rss_kb = Union{Missing, Int}[],
        baseline_rss_kb = Union{Missing, Int}[],
        assemble_rss_delta_kb = Union{Missing, Int}[],
        gc_live_delta_mb = Union{Missing, Float64}[],
        quast_runtime_s = Union{Missing, Float64}[],
        worker_ok = Bool[]
    )

    println("\n--- Phase 3: assemble each mode in a fresh subprocess, then QUAST ---")
    for m in modes
        println("\n  mode: $(m.label) (graph_mode=$(m.graph_mode), dedup=$(m.dedup), " *
                "memory_profile=$(m.memprofile), compact=$(m.compact))")
        mode_dir = joinpath(benchmark_dir, "mode_$(replace(m.label, r"[^A-Za-z0-9]" => "_"))")
        mkpath(mode_dir)
        contigs_out = joinpath(mode_dir, "contigs.fasta")
        metrics_out = joinpath(mode_dir, "metrics.txt")

        worker_ok = true
        # Re-invoke THIS script in a fresh worker process. Base.julia_cmd()
        # preserves the sysimage/flags; --project pins the environment; addenv
        # augments (never replaces) the inherited env so LD_LIBRARY_PATH,
        # JULIA_DEPOT_PATH, JULIA_NUM_THREADS set by the run command carry over.
        worker_cmd = `$(Base.julia_cmd()) --project=$(project_dir) --threads=$(Threads.nthreads()) $(script_path)`
        worker_cmd = addenv(worker_cmd,
            "MYCELIA_MODE_WORKER" => "1",
            "MMC_MODE_LABEL" => m.label,
            "MMC_GRAPH_MODE" => m.graph_mode,
            "MMC_DEDUP" => string(m.dedup),
            "MMC_MEMPROFILE" => m.memprofile,
            "MMC_COMPACT" => string(m.compact),
            "MMC_K" => string(k),
            "MMC_FWD" => fwd,
            "MMC_REV" => rev,
            "MMC_CONTIGS_OUT" => contigs_out,
            "MMC_METRICS_OUT" => metrics_out)
        try
            run(worker_cmd)
        catch e
            worker_ok = false
            @warn "Worker process failed for mode=$(m.label) — recording as failure" exception = (e,)
        end

        metrics = parse_metrics_file(metrics_out)
        n_contigs = haskey(metrics, "n_contigs") ? tryparse(Int, metrics["n_contigs"]) : missing
        assemble_runtime_s = haskey(metrics, "assemble_runtime_s") ? tryparse(Float64, metrics["assemble_runtime_s"]) : missing
        assemble_peak_rss_kb = haskey(metrics, "assemble_peak_rss_kb") ? tryparse(Int, metrics["assemble_peak_rss_kb"]) : missing
        baseline_rss_kb = haskey(metrics, "baseline_rss_kb") ? tryparse(Int, metrics["baseline_rss_kb"]) : missing
        assemble_rss_delta_kb = haskey(metrics, "assemble_rss_delta_kb") ? tryparse(Int, metrics["assemble_rss_delta_kb"]) : missing
        gc_live_delta_mb = haskey(metrics, "gc_live_delta_mb") ? tryparse(Float64, metrics["gc_live_delta_mb"]) : missing

        # -- QUAST aligned-truth scoring on the produced contigs --
        nga50 = missing
        num_misassemblies = missing
        dup_ratio = missing
        genome_fraction = missing
        quast_runtime_s = missing
        have_contigs = (n_contigs !== missing && n_contigs !== nothing && n_contigs > 0 &&
                        isfile(contigs_out) && filesize(contigs_out) > 0)
        if run_quast && have_contigs
            quast_dir = joinpath(mode_dir, "quast")
            tq = time()
            try
                Mycelia.run_quast(
                    contigs_out;
                    outdir = quast_dir,
                    reference = reference_path,
                    min_contig = quast_min_contig)
                quast_runtime_s = round(time() - tq; digits = 3)
                report_tsv = joinpath(quast_dir, "report.tsv")
                nga50 = parse_quast_metric(report_tsv, "NGA50")
                num_misassemblies = parse_quast_metric(report_tsv, "# misassemblies")
                dup_ratio = parse_quast_metric(report_tsv, "Duplication ratio")
                genome_fraction = parse_quast_metric(report_tsv, "Genome fraction (%)")
            catch e
                @warn "QUAST failed for mode=$(m.label) — leaving accuracy metrics missing" exception = (e,)
            end
        elseif run_quast && !have_contigs
            println("    (no contigs to score — mode likely BROKEN; leaving QUAST metrics missing)")
        end

        push!(results,
            (
                organism = organism,
                mode_label = m.label,
                graph_mode = m.graph_mode,
                dedup = m.dedup,
                memory_profile = m.memprofile,
                compact = m.compact,
                n_contigs = (n_contigs === nothing ? missing : n_contigs),
                NGA50 = nga50,
                num_misassemblies = num_misassemblies,
                dup_ratio = dup_ratio,
                genome_fraction_quast = genome_fraction,
                assemble_runtime_s = (assemble_runtime_s === nothing ? missing : assemble_runtime_s),
                assemble_peak_rss_kb = (assemble_peak_rss_kb === nothing ? missing : assemble_peak_rss_kb),
                baseline_rss_kb = (baseline_rss_kb === nothing ? missing : baseline_rss_kb),
                assemble_rss_delta_kb = (assemble_rss_delta_kb === nothing ? missing : assemble_rss_delta_kb),
                gc_live_delta_mb = (gc_live_delta_mb === nothing ? missing : gc_live_delta_mb),
                quast_runtime_s = quast_runtime_s,
                worker_ok = worker_ok
            ))
    end

    # --- Phase 4: write CSV ---
    # UTC timestamp, yyyymmdd_HHMMSS (no 'Z' directive — 'Z' is a TimeZones token
    # that errors on a plain Dates.DateTime).
    timestamp = Dates.format(Dates.now(Dates.UTC), "yyyymmdd_HHMMSS")
    csv_path = joinpath(results_dir, "mode_comparison_$(timestamp).csv")
    CSV.write(csv_path, results)
    println("\nResults written to: $(csv_path)")

    # --- Phase 5: comparison table + verdict ---
    print_comparison_table(results)
    print_verdict(results)

    println("\n=== Mode-Comparison complete (UTC): $(Dates.now(Dates.UTC)) ===")
    return nothing
end

"""
    fmt(x) -> String

Compact fixed-width-friendly formatter for possibly-missing metric values.
"""
function fmt(x)::String
    ismissing(x) && return "NA"
    x isa AbstractFloat && return string(round(x; digits = 3))
    return string(x)
end

"""
    print_comparison_table(results)

Print a mode x {dup, GF, NGA50, misasm, runtime, memory} comparison table.
"""
function print_comparison_table(results::DataFrames.DataFrame)
    println("\n--- Comparison table (mode x accuracy/efficiency) ---")
    header = rpad("mode", 18) * rpad("contigs", 9) * rpad("dup", 8) *
             rpad("GF%", 9) * rpad("NGA50", 10) * rpad("misasm", 8) *
             rpad("run_s", 9) * rpad("peak_rss_kb", 13) * rpad("rss_delta_kb", 14) *
             rpad("live_MB", 9)
    println(header)
    println(repeat("-", length(header)))
    for row in eachrow(results)
        line = rpad(row.mode_label, 18) *
               rpad(fmt(row.n_contigs), 9) *
               rpad(fmt(row.dup_ratio), 8) *
               rpad(fmt(row.genome_fraction_quast), 9) *
               rpad(fmt(row.NGA50), 10) *
               rpad(fmt(row.num_misassemblies), 8) *
               rpad(fmt(row.assemble_runtime_s), 9) *
               rpad(fmt(row.assemble_peak_rss_kb), 13) *
               rpad(fmt(row.assemble_rss_delta_kb), 14) *
               rpad(fmt(row.gc_live_delta_mb), 9)
        println(line)
    end
end

"""
    print_verdict(results)

Compare every mode against the `baseline` row and print a short verdict block:
which modes preserve accuracy (GF / misasm / NGA50 unchanged), which reduce the
duplication ratio, which reduce memory, and flag `canonical` as BROKEN if its
genome fraction collapsed (or its worker/QUAST produced no aligned metrics).
"""
function print_verdict(results::DataFrames.DataFrame)
    println("\n--- Verdict ---")
    base_rows = results[results.mode_label .== "baseline", :]
    if DataFrames.nrow(base_rows) == 0
        println("No baseline row present — cannot compute relative verdict.")
        return nothing
    end
    base = base_rows[1, :]
    b_gf = base.genome_fraction_quast
    b_misasm = base.num_misassemblies
    b_nga50 = base.NGA50
    b_dup = base.dup_ratio
    b_mem = base.assemble_rss_delta_kb   # assembly-attributable peak growth

    approx(a, b; rtol = 0.02) = (!ismissing(a) && !ismissing(b) && b != 0) ? (abs(a - b) / abs(b) <= rtol) : (a === b)

    for row in eachrow(results)
        row.mode_label == "baseline" && continue
        notes = String[]

        # Canonical brokenness flag: collapsed GF, no contigs, or failed worker.
        broken = row.mode_label == "canonical" &&
                 (!row.worker_ok ||
                  ismissing(row.genome_fraction_quast) ||
                  (!ismissing(row.genome_fraction_quast) && !ismissing(b_gf) && row.genome_fraction_quast < 0.5 * b_gf) ||
                  (!ismissing(row.n_contigs) && row.n_contigs == 0))
        if broken
            gf_str = ismissing(row.genome_fraction_quast) ? "NA (nothing aligned)" : "$(fmt(row.genome_fraction_quast))%"
            println("  [$(row.mode_label)] BROKEN — genome fraction $(gf_str) vs baseline $(fmt(b_gf))% " *
                    "(empirical proof the canonical assembly is incorrect).")
            continue
        end

        # Accuracy preservation.
        acc_preserved = approx(row.genome_fraction_quast, b_gf) &&
                        approx(row.NGA50, b_nga50) &&
                        (ismissing(row.num_misassemblies) == ismissing(b_misasm)) &&
                        (ismissing(row.num_misassemblies) || row.num_misassemblies == b_misasm)
        push!(notes, acc_preserved ? "accuracy PRESERVED (GF/NGA50/misasm ~ baseline)" :
              "accuracy CHANGED (GF=$(fmt(row.genome_fraction_quast)) vs $(fmt(b_gf)), " *
              "NGA50=$(fmt(row.NGA50)) vs $(fmt(b_nga50)), misasm=$(fmt(row.num_misassemblies)) vs $(fmt(b_misasm)))")

        # Duplication reduction.
        if !ismissing(row.dup_ratio) && !ismissing(b_dup)
            if row.dup_ratio < b_dup - 1.0e-6
                push!(notes, "dup REDUCED $(fmt(b_dup)) -> $(fmt(row.dup_ratio))")
            elseif approx(row.dup_ratio, b_dup; rtol = 0.001)
                push!(notes, "dup unchanged ($(fmt(row.dup_ratio)))")
            else
                push!(notes, "dup INCREASED $(fmt(b_dup)) -> $(fmt(row.dup_ratio))")
            end
        end

        # Memory reduction (assembly-attributable peak growth).
        if !ismissing(row.assemble_rss_delta_kb) && !ismissing(b_mem)
            if row.assemble_rss_delta_kb < b_mem
                push!(notes, "memory REDUCED $(b_mem) -> $(row.assemble_rss_delta_kb) kb (assemble delta)")
            elseif row.assemble_rss_delta_kb == b_mem
                push!(notes, "memory unchanged")
            else
                push!(notes, "memory HIGHER $(b_mem) -> $(row.assemble_rss_delta_kb) kb")
            end
        end

        println("  [$(row.mode_label)] " * join(notes, "; "))
    end

    println("\nNote: accuracy comparisons use QUAST aligned-truth metrics; memory uses the")
    println("fresh-subprocess assemble-attributable peak-RSS delta (baseline-subtracted).")
    return nothing
end

# ===========================================================================
# Entry point: worker vs driver
# ===========================================================================

if get(ENV, "MYCELIA_MODE_WORKER", "") == "1"
    run_worker()
else
    run_driver()
end
