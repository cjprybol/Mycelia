# Track A baseline benchmark — greedy Viterbi on the viral tier (rhizomorph-paper td-bblmi)
#
# Runs the CURRENT default Rhizomorph assembler (greedy path-finding) across the
# viral tier to establish the baseline for all H1-H7 comparisons and to confirm
# the pre-registration power analysis (assumed NGA50 CV ~ 0.15) before the
# registration is locked.
#
# Matrix (288 cells): 4 organisms x 3 technologies x 4 coverages x 3 seeds x 2 decoder arms.
#   - decoder arm "qualmer": simulated FASTQ passed as-is -> assemble_genome auto-enables
#     quality scores -> qualmer graph (the assembler's real default on real reads).
#   - decoder arm "kmer":    quality stripped (FASTQ -> FASTA records) -> plain k-mer graph,
#     matching the existing FASTA-based benchmark and the future DP arm's graph type.
#
# Each cell: simulate reads -> assemble (k=31) -> QUAST vs reference -> parse NGA50 /
# misassemblies / genome fraction / duplication ratio. Per-cell JSON checkpoint enables
# crash-safe resume. A final step computes NGA50 CV per (organism x tech x coverage x arm)
# and writes a pass/fail power-analysis summary.
#
# Usage:
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl              # full 288-cell run
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --smoke      # 1 cell (Lambda/illumina/30x/seed42/qualmer)
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --organisms Lambda,T4 --arms kmer
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --coverages 30,100 --seeds 42 --technologies illumina,ont
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --output-dir /scratch/track_a
#
# Shard flags (--organisms/--technologies/--coverages/--seeds/--arms) take comma-separated
# values and compose, so an HPC array job can split the matrix and share one results tree.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import DataFrames
import CSV
import JSON
import Dates
import Statistics

# Optional provenance writer (git SHA + tool versions). Degrade gracefully if unavailable.
const HAVE_ARTIFACT_WRITER = try
    include(joinpath(@__DIR__, "benchmark_artifacts.jl"))
    true
catch e
    @warn "benchmark_artifacts.jl unavailable; skipping provenance bundle" exception = e
    false
end

# === Configuration ===

# (display_name, NCBI accession, expected_size_bp). Lambda/T4/SARS-CoV-2 accessions match
# benchmarking/real_genome_benchmark.jl; phi29 = Bacillus phage phi29 (NC_011048, ~19.3 kb).
const ORGANISMS = [
    ("Lambda", "NC_001416", 48_502),
    ("T4", "NC_000866", 168_903),
    ("phi29", "NC_011048", 19_282),
    ("SARS-CoV-2", "NC_045512", 29_903)
]
const TECHNOLOGIES = ["illumina", "pacbio", "ont"]
const COVERAGES = [10, 30, 50, 100]
const SEEDS = [42, 123, 456]
const DECODER_ARMS = ["qualmer", "kmer"]
const K = 31
const CV_THRESHOLD = 0.15  # assumed NGA50 coefficient of variation in the power analysis

# Canonical row schema (fixed order so in-memory and JSON-reloaded rows align in the DataFrame).
#
# `peak_rss_bytes` carries its own provenance: `rss_baseline_bytes` (process RSS on entry
# to the cell) and `peak_rss_method` (HOW the peak was obtained). Rows measured by
# different methods are different quantities, so the method has to travel with the number
# rather than be inferred from the run that produced it. New columns are appended ahead of
# `:status` only; the preceding column order is unchanged so existing readers of
# `track_a_results.tsv` keep working.
const ROW_KEYS = (
    :organism, :accession, :technology, :coverage, :seed, :decoder_arm, :k,
    :n_reads, :n_contigs, :NGA50, :misassemblies, :genome_fraction,
    :duplication_ratio, :largest_contig, :wall_seconds, :peak_rss_bytes,
    :rss_baseline_bytes, :peak_rss_method, :status
)
const INT_KEYS = (
    :coverage, :seed, :k, :n_reads, :n_contigs, :largest_contig, :peak_rss_bytes,
    :rss_baseline_bytes)
const FLOAT_KEYS = (
    :NGA50, :misassemblies, :genome_fraction, :duplication_ratio, :wall_seconds)
const STR_KEYS = (
    :organism, :accession, :technology, :decoder_arm, :peak_rss_method, :status)

# === Argument parsing ===

arg_value(flag) =
    let i = findfirst(==(flag), ARGS)
        (i !== nothing && i < length(ARGS)) ? ARGS[i + 1] : nothing
    end
arg_list(flag) =
    let v = arg_value(flag)
        v === nothing ? nothing : String.(split(v, ","))
    end

const SMOKE = "--smoke" in ARGS

organisms = ORGANISMS
technologies = TECHNOLOGIES
coverages = COVERAGES
seeds = SEEDS
arms = DECODER_ARMS

if SMOKE
    organisms = ORGANISMS[1:1]      # Lambda (accession verified in-repo)
    technologies = ["illumina"]
    coverages = [30]
    seeds = [42]
    arms = ["qualmer"]
else
    # NOTE: assign directly (no `let`). A `let` block introduces a new scope, so
    # `coverages = ...` inside it would create a local shadow and silently leave
    # the global matrix at its default — i.e. the shard/filter flags would be
    # ignored. `if/else` bodies do not introduce scope, so these reach the globals.
    _f = arg_list("--organisms")
    _f !== nothing && (organisms = filter(o -> o[1] in _f, ORGANISMS))
    _f = arg_list("--technologies")
    _f !== nothing && (technologies = _f)
    _f = arg_list("--coverages")
    _f !== nothing && (coverages = parse.(Int, _f))
    _f = arg_list("--seeds")
    _f !== nothing && (seeds = parse.(Int, _f))
    _f = arg_list("--arms")
    _f !== nothing && (arms = _f)
end

const OUTPUT_DIR = let v = arg_value("--output-dir")
    v === nothing ? joinpath(@__DIR__, "results", "track_a_baseline") : v
end

const N_CELLS = length(organisms) * length(technologies) * length(coverages) *
                length(seeds) * length(arms)

# === Metrics + row helpers ===

function empty_metrics()
    (NGA50 = 0.0, misassemblies = 0.0, genome_fraction = 0.0,
        duplication_ratio = 0.0, largest_contig = 0.0)
end

function cell_row(org, acc, tech, cov, seed, arm; n_reads, n_contigs,
        wall_seconds, peak_rss_bytes, rss_baseline_bytes, peak_rss_method,
        metrics, status)
    return (
        organism = String(org), accession = String(acc), technology = String(tech),
        coverage = Int(cov), seed = Int(seed), decoder_arm = String(arm), k = K,
        n_reads = Int(n_reads), n_contigs = Int(n_contigs),
        NGA50 = Float64(metrics.NGA50), misassemblies = Float64(metrics.misassemblies),
        genome_fraction = Float64(metrics.genome_fraction),
        duplication_ratio = Float64(metrics.duplication_ratio),
        largest_contig = Int(round(Float64(metrics.largest_contig))),
        wall_seconds = round(Float64(wall_seconds); digits = 3),
        peak_rss_bytes = Int(peak_rss_bytes),
        rss_baseline_bytes = Int(rss_baseline_bytes),
        peak_rss_method = String(peak_rss_method), status = String(status)
    )
end

# Provenance columns a checkpoint written before they existed is allowed to lack, and the
# sentinel each takes. Every OTHER column is a measurement or a grouping key, and a missing
# one must fail loudly: `NGA50 = 0.0` and `genome_fraction = 0.0` are the GENUINE values for
# every ONT cell at 10x/30x, so a zero-filled truncated checkpoint is byte-identical to a
# real result and would be pooled straight into the CV that gates the pre-registration. A
# truncated checkpoint must not be able to move a verdict.
#
# Driving the defaulting from this table (rather than hardcoding the two cases inside
# `canonical`) is what makes the error message below true: adding a key here really is all
# it takes to make that key defaultable.
const OPTIONAL_KEY_DEFAULTS = Dict{Symbol, Any}(
    :peak_rss_method => "unknown", :rss_baseline_bytes => -1)
const OPTIONAL_KEYS = Tuple(sort(collect(keys(OPTIONAL_KEY_DEFAULTS))))

# Rebuild a canonical, type-coerced NamedTuple from a parsed JSON dict (resumed cells).
#
# An absent OPTIONAL_KEY takes an explicit "unknown"/-1 sentinel that no real measurement
# can produce, so a resumed pre-schema tree is reloadable but can never be mistaken for a
# measured one. Any other absent key raises: indexing `d[String(key)]` directly raised a
# bare KeyError naming only the key, which killed the run on the first cached cell with no
# indication that the tree simply predated a schema change.
function canonical(d::AbstractDict)
    vals = map(ROW_KEYS) do key
        name = String(key)
        if !haskey(d, name)
            haskey(OPTIONAL_KEY_DEFAULTS, key) && return OPTIONAL_KEY_DEFAULTS[key]
            error("checkpoint is missing required key $(name). This is a measurement or " *
                  "grouping column, so it cannot be defaulted — a zero here is " *
                  "indistinguishable from a real result and would silently enter the " *
                  "power analysis. Delete the checkpoint to recompute the cell, or add " *
                  "$(name) to OPTIONAL_KEY_DEFAULTS if it is genuinely provenance-only.")
        end
        v = d[name]
        if key in INT_KEYS
            v isa Integer ? Int(v) : Int(round(Float64(v)))
        elseif key in FLOAT_KEYS
            Float64(v)
        else
            String(v)
        end
    end
    return (; (ROW_KEYS .=> vals)...)
end

function save_cell_json(path, row)
    open(path, "w") do io
        JSON.print(io, Dict(string(k) => v for (k, v) in pairs(row)), 2)
    end
    return path
end

# === QUAST report parsing ===

# report.tsv: column 1 = metric label, column 2 = this assembly's value. NGA50 and
# misassembly rows appear only when QUAST is run with a reference.
function parse_quast_report(report_tsv)
    isfile(report_tsv) || return empty_metrics()
    df = CSV.read(report_tsv, DataFrames.DataFrame; delim = '\t', header = true)
    DataFrames.ncol(df) >= 2 || return empty_metrics()
    labelcol, valcol = names(df)[1], names(df)[2]
    getval(metric) =
        let i = findfirst(==(metric), df[!, labelcol])
            i === nothing ? missing : df[i, valcol]
        end
    tonum(x) = x === missing ? 0.0 :
               x isa Number ? Float64(x) :
               something(tryparse(Float64, strip(string(x))), 0.0)
    return (
        NGA50 = tonum(getval("NGA50")),
        misassemblies = tonum(getval("# misassemblies")),
        genome_fraction = tonum(getval("Genome fraction (%)")),
        duplication_ratio = tonum(getval("Duplication ratio")),
        largest_contig = tonum(getval("Largest contig"))
    )
end

# === Read simulation (per-technology adapter; APIs are non-uniform) ===

function simulate_reads(tech, ref_fasta, cov, seed, reads_dir)
    mkpath(reads_dir)
    if tech == "illumina"
        outbase = joinpath(reads_dir, "illumina_$(cov)x")
        res = Mycelia.simulate_illumina_reads(
            fasta = ref_fasta, coverage = cov, outbase = outbase,
            rndSeed = seed, paired = true, quiet = true)
        recs = FASTX.FASTQ.Record[]
        for p in (res.forward_reads, res.reverse_reads)
            p === nothing && continue
            reader = Mycelia.open_fastx(p)
            append!(recs, collect(reader))
            close(reader)
        end
        return recs
    elseif tech == "pacbio"
        fq = Mycelia.simulate_pacbio_reads(
            fasta = ref_fasta, quantity = "$(cov)x",
            outfile = joinpath(reads_dir, "pacbio_$(cov)x.fq.gz"),
            seed = seed, quiet = true)
        reader = Mycelia.open_fastx(fq)
        recs = collect(reader)
        close(reader)
        return recs
    elseif tech == "ont"
        fq = Mycelia.simulate_nanopore_reads(
            fasta = ref_fasta, quantity = "$(cov)x",
            outfile = joinpath(reads_dir, "ont_$(cov)x.fq.gz"),
            seed = seed, quiet = true)
        reader = Mycelia.open_fastx(fq)
        recs = collect(reader)
        close(reader)
        return recs
    else
        error("unknown technology: $tech")
    end
end

# Quality-on arm keeps FASTQ records; quality-off arm strips quality to FASTA records.
function shape_for_arm(records, arm)
    if arm == "qualmer"
        return records
    elseif arm == "kmer"
        return [FASTX.FASTA.Record(FASTX.identifier(r), FASTX.sequence(r)) for r in records]
    else
        error("unknown decoder arm: $arm")
    end
end

# === Memory measurement ===

# `Sys.maxrss()` is the process-lifetime high-water mark and is monotonically
# non-decreasing, so the previous `max(0, Sys.maxrss() - rss0)` could only be nonzero when
# a cell beat the ALL-TIME process peak. It measured "excess over historical peak", not
# per-cell peak: on the 2026-07-25 Lovelace run 271/288 rows (94%) were exactly 0,
# including 144/144 kmer cells, and the few survivors understate real usage by orders of
# magnitude because each is a delta over the previous peak rather than a measurement.
#
# Read VmRSS from /proc/self/status instead: it reports explicit kB, so unlike
# /proc/self/statm it needs no `Sys.PAGESIZE` (which does not exist in Julia 1.10).
# Returns `nothing` wherever /proc is unavailable (macOS), which selects the fallback.
function current_rss_bytes()
    status_path = "/proc/self/status"
    isfile(status_path) || return nothing
    # `open(...) do` rather than a bare `eachline(path)`: eachline closes its handle only
    # at end-of-iteration, and VmRSS sits near the top of the file so the early return
    # below always exits first — leaking one descriptor per poll, ~72k/hour at 20 Hz.
    return open(status_path) do io
        for line in eachline(io)
            startswith(line, "VmRSS:") || continue
            fields = split(line)
            length(fields) >= 2 || return nothing
            kilobytes = tryparse(Int, fields[2])
            return kilobytes === nothing ? nothing : kilobytes * 1024
        end
        return nothing
    end
end

"""
Run `f`, returning `(; value, wall_seconds, peak_rss_bytes, rss_baseline_bytes, method)`.

`peak_rss_bytes` is the highest WHOLE-PROCESS resident set observed while `f` ran. It is
deliberately not called a per-cell figure: it includes the Julia runtime, the loaded
package images, every other thread, and memory retained (not returned to the OS) by
earlier cells. `rss_baseline_bytes` is the process RSS immediately before `f`, so a caller
wanting a cell-attributable increment can take the difference — but the absolute level is
what an HPC memory request actually needs, so that is what the primary column holds.

The sampler runs via `Threads.@spawn`, which requires a second thread: the assembly is
CPU-bound Julia that never yields, so a same-thread `@async` sampler would never be
scheduled. Both committed sbatch wrappers export JULIA_NUM_THREADS=\$SLURM_CPUS_PER_TASK —
run_track_a_baseline_lrc.sbatch (16) and run_track_a_baseline_nersc.sbatch (4) — so
sampling is the production path under SLURM. Julia defaults to ONE thread, so a bare
`julia track_a_baseline_benchmark.jl` silently reverts to the broken high-water semantics;
launch with `julia -t N` or the column is meaningless.

`method` records how the number was obtained, so a reader never has to infer it and rows
measured different ways are never averaged together:
  "sampled"          - process-RSS high-water, polled while the cell ran
  "sampled-degraded" - the sampler never completed a single read (starved, or every /proc
                       read failed); the number is the entry baseline, NOT a measurement.
                       A flat cell whose RSS never exceeds its baseline is still "sampled"
                       — that is a real observation.
  "highwater-delta"  - single-threaded or no /proc: the old, known-broken semantics
  "unknown"          - checkpoint predates this column (see OPTIONAL_KEY_DEFAULTS)
Values under different methods are different quantities. Always filter on
`peak_rss_method` before aggregating.
"""
function timed_with_peak_rss(f; interval_seconds::Float64 = 0.05)
    # Guard the baseline read too: everything here is careful never to let telemetry fail
    # science, and an unguarded /proc read error would abort the cell before f() even ran.
    baseline = try
        current_rss_bytes()
    catch e
        @warn "baseline RSS read failed; falling back to high-water semantics" exception = e
        nothing
    end

    can_sample = baseline !== nothing && Threads.nthreads() > 1
    if !can_sample
        rss0 = Sys.maxrss()
        timed = @timed f()
        return (value = timed.value, wall_seconds = timed.time,
            peak_rss_bytes = max(0, Sys.maxrss() - rss0),
            rss_baseline_bytes = baseline === nothing ? -1 : baseline,
            method = "highwater-delta")
    end

    peak = Threads.Atomic{Int}(baseline)
    keep_sampling = Threads.Atomic{Bool}(true)
    # Count SUCCESSFUL READS, not increases over the baseline. The failure worth flagging
    # is a sampler that never observed anything — starved, or every /proc read failed —
    # because it then returns the baseline as though it were a measurement: a large,
    # stable, plausible byte count, strictly harder to spot than the exactly-0 signature
    # of the bug being replaced. A cell whose RSS never rises above its baseline is NOT
    # degraded; that happens whenever the allocation fits inside already-resident heap.
    reads = Threads.Atomic{Int}(0)
    sampler = Threads.@spawn begin
        # The sampler must never be able to fail a good cell: a transient /proc read error
        # is telemetry, not science. Swallow it and let the read count report the
        # degradation instead.
        try
            while keep_sampling[]
                observed = current_rss_bytes()
                if observed !== nothing
                    Threads.atomic_add!(reads, 1)
                    observed > peak[] && Threads.atomic_max!(peak, observed)
                end
                sleep(interval_seconds)
            end
        catch
            # deliberately swallowed; `reads == 0` reports the degradation
        end
    end

    timed = try
        @timed f()
    finally
        # Stop AND reap the sampler here, not after the block. With `wait` outside, a
        # sampler exception on the success path would raise TaskFailedException and
        # DISCARD a completed assembly — telemetry destroying science — while on the throw
        # path `wait` would never be reached at all and the sampler would outlive the cell.
        keep_sampling[] = false
        try
            wait(sampler)
        catch
            # sampler already reported via `reads`; never let it mask f()'s own outcome
        end
    end
    return (value = timed.value, wall_seconds = timed.time,
        peak_rss_bytes = peak[], rss_baseline_bytes = baseline,
        method = reads[] > 0 ? "sampled" : "sampled-degraded")
end

# === Per-cell execution ===

function run_cell(org, acc, ref, tech, cov, seed, arm, cell_dir)
    records = simulate_reads(tech, ref, cov, seed, joinpath(cell_dir, "reads"))
    n_reads = length(records)
    asm_input = shape_for_arm(records, arm)

    GC.gc()
    measured = timed_with_peak_rss(
        () -> Mycelia.Rhizomorph.assemble_genome(asm_input; k = K, verbose = false))
    result = measured.value
    wall_seconds = measured.wall_seconds
    peak_rss_bytes = measured.peak_rss_bytes
    rss_baseline_bytes = measured.rss_baseline_bytes
    peak_rss_method = measured.method

    contigs_path = joinpath(cell_dir, "contigs.fasta")
    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)  # Rhizomorph contigs are String
        end
    end
    n_contigs = length(result.contigs)

    metrics = if n_contigs > 0 && filesize(contigs_path) > 0
        try
            quast_dir = joinpath(cell_dir, "quast")
            Mycelia.run_quast([contigs_path]; outdir = quast_dir, reference = ref, min_contig = 500)
            parse_quast_report(joinpath(quast_dir, "report.tsv"))
        catch e
            @warn "QUAST failed" cell=basename(cell_dir) exception=(e, catch_backtrace())
            empty_metrics()
        end
    else
        empty_metrics()
    end

    status = n_contigs == 0 ? "empty_assembly" : "ok"
    return cell_row(org, acc, tech, cov, seed, arm;
        n_reads, n_contigs, wall_seconds, peak_rss_bytes, rss_baseline_bytes,
        peak_rss_method, metrics, status)
end

# === Aggregation ===

function write_aggregate(root, rows)
    df = DataFrames.DataFrame(rows)
    CSV.write(joinpath(root, "track_a_results.tsv"), df; delim = '\t')
    return df
end

function write_power_analysis(root, df)
    cv_rows = NamedTuple[]
    for g in DataFrames.groupby(df, [:organism, :technology, :coverage, :decoder_arm])
        nga = Float64.(g.NGA50)
        m = Statistics.mean(nga)
        s = length(nga) > 1 ? Statistics.std(nga; corrected = true) : NaN
        cv = (m == 0 || isnan(s)) ? NaN : s / m
        push!(cv_rows,
            (
                organism = g.organism[1], technology = g.technology[1],
                coverage = g.coverage[1], decoder_arm = g.decoder_arm[1],
                n = length(nga), mean_nga50 = round(m; digits = 1),
                sd_nga50 = round(s; digits = 1), cv_nga50 = round(cv; digits = 4),
                passes = (isfinite(cv) && cv <= CV_THRESHOLD)))
    end
    cv_df = DataFrames.DataFrame(cv_rows)
    CSV.write(joinpath(root, "power_analysis_cv.tsv"), cv_df; delim = '\t')

    evaluable = filter(r -> isfinite(r.cv_nga50), cv_rows)
    n_eval = length(evaluable)
    n_pass = count(r -> r.passes, evaluable)
    max_cv = isempty(evaluable) ? NaN : maximum(r.cv_nga50 for r in evaluable)
    verdict = (n_eval > 0 && n_pass == n_eval) ? "supported" :
              n_eval == 0 ? "indeterminate (no evaluable cells)" : "NOT fully supported"

    open(joinpath(root, "power_analysis_summary.md"), "w") do io
        println(io, "# Track A power-analysis check — NGA50 CV vs assumed $(CV_THRESHOLD)\n")
        println(io, "Generated: $(Dates.now())\n")
        println(io, "- Evaluable cells (organism×tech×coverage×arm, ≥2 seeds, nonzero NGA50): $n_eval")
        println(io, "- Cells with CV ≤ $(CV_THRESHOLD): $n_pass / $n_eval")
        println(io, "- Max CV observed: $(round(max_cv; digits = 4))")
        println(io, "- **Verdict: assumed CV ≈ $(CV_THRESHOLD) is $verdict.**\n")
        fails = filter(r -> isfinite(r.cv_nga50) && !r.passes, cv_rows)
        if !isempty(fails)
            println(io, "## Cells exceeding CV $(CV_THRESHOLD)\n")
            for r in fails
                println(io,
                    "- $(r.organism) / $(r.technology) / $(r.coverage)x / $(r.decoder_arm): " *
                    "CV=$(r.cv_nga50) (mean NGA50 $(r.mean_nga50), n=$(r.n))")
            end
            println(io)
        end
        println(io,
            "_Caveat: n=3 seeds makes each CV estimate noisy; treat this as directional " *
            "evidence for the power analysis, not a definitive variance estimate._")
    end
    return cv_df
end

# === Cell identity + error rows (extracted so both are unit-testable) ===

# k belongs in the cell id: without it, two runs at different k write the same
# cells/<id>/cell_result.json and the second silently resumes the first's result.
#
# But the suffix is applied ONLY for a non-default k. Interpolating it unconditionally
# renames the DEFAULT k=31 cells too, which orphans every checkpoint written before k was
# selectable — 288 on Lovelace and 144 on LRC — turning the sbatch wrappers' documented
# "just re-submit, completed cells are skipped" recovery into a silent full recompute.
# Legacy trees are k=31 by definition (k was a hardcoded const), so keying k=31 to the
# historical name is unambiguous and needs no migration.
#
# The unsuffixed branch is pinned to the literal 31, not to `K`: the reservation belongs to
# the names already on disk, so changing the default k must never re-point the unsuffixed
# name at a different k and republish a completed tree as if it were the new k. The
# `cell_id = ...` line is also the exact string track_a_merge_hosts_test.jl's drift guard
# asserts against this source, keeping the merge tool's mirrored id format verifiable.
function cell_id_for(org, tech, cov, seed, arm; k = K)
    cell_id = "$(org)__$(tech)__$(cov)x__seed$(seed)__$(arm)"
    return k == 31 ? cell_id : "$(cell_id)__k$(k)"
end

# The row recorded when a cell throws. Extracted from the catch block because it was
# unreachable from any test while it lived there: adding a required keyword to `cell_row`
# without updating this call makes the ERROR HANDLER ITSELF throw UndefKeywordError, so
# the first failing cell aborts the whole matrix from inside its own recovery path — and
# `--smoke` runs a single deliberately-successful cell, so no test could reach it.
error_row(org, acc, tech, cov, seed, arm) =
    cell_row(org, acc, tech, cov, seed, arm;
        n_reads = 0, n_contigs = 0, wall_seconds = 0.0, peak_rss_bytes = 0,
        rss_baseline_bytes = -1, peak_rss_method = "unknown",
        metrics = empty_metrics(), status = "error")

# === Main ===

# Guard the driver so a test can `include()` this file for its pure helpers without
# downloading genomes or running assemblies — which is what makes
# test/4_assembly/track_a_baseline_benchmark_test.jl possible, and matches the convention
# already used by 18 other harnesses in this directory. `if` does not introduce scope at
# top level, so the globals below stay global.
if abspath(PROGRAM_FILE) == @__FILE__

println("=== Track A baseline benchmark ===")
println("Start: $(Dates.now())")
println("Smoke mode: $SMOKE")
println("Organisms: $(join((o[1] for o in organisms), ", "))")
println("Technologies: $(join(technologies, ", "))")
println("Coverages: $(join(coverages, ", "))x")
println("Seeds: $(join(seeds, ", "))")
println("Decoder arms: $(join(arms, ", "))")
println("Cells to run: $N_CELLS")
println("Output dir: $OUTPUT_DIR")

mkpath(OUTPUT_DIR)
refs_dir = joinpath(OUTPUT_DIR, "refs")
cells_dir = joinpath(OUTPUT_DIR, "cells")
mkpath(refs_dir)
mkpath(cells_dir)

println("\n--- Phase 1: download reference genomes ---")
ref_paths = Dict{String, String}()
for (org, acc, expected) in organisms
    haskey(ref_paths, org) && continue
    print("  $org ($acc) ... ")
    rp = Mycelia.download_genome_by_accession(accession = acc, outdir = refs_dir, compressed = false)
    if !isfile(rp) || filesize(rp) == 0
        error("download failed for $org ($acc): $rp")
    end
    actual = try
        Mycelia.total_fasta_size(rp)
    catch
        -1
    end
    if actual > 0 && abs(actual - expected) > 0.2 * expected
        @warn "reference size mismatch" organism=org expected=expected actual=actual
    end
    ref_paths[org] = rp
    println("$(actual > 0 ? actual : "?") bp")
end

println("\n--- Phase 2: matrix ($N_CELLS cells) ---")
rows = NamedTuple[]
cell_index = 0
for (org, acc, _expected) in organisms,
    tech in technologies,
    cov in coverages,
    seed in seeds,
    arm in arms
    global cell_index += 1
    cell_id = cell_id_for(org, tech, cov, seed, arm)
    cell_dir = joinpath(cells_dir, cell_id)
    ckpt = joinpath(cell_dir, "cell_result.json")

    if isfile(ckpt)
        println("  [$cell_index/$N_CELLS] $cell_id — cached, skipping")
        push!(rows, canonical(JSON.parsefile(ckpt)))
        continue
    end

    mkpath(cell_dir)
    print("  [$cell_index/$N_CELLS] $cell_id ... ")
    row = try
        run_cell(org, acc, ref_paths[org], tech, cov, seed, arm, cell_dir)
    catch e
        @warn "cell failed" cell_id exception = (e, catch_backtrace())
        error_row(org, acc, tech, cov, seed, arm)
    end
    save_cell_json(ckpt, row)
    push!(rows, row)
    write_aggregate(OUTPUT_DIR, rows)  # rewrite aggregate each cell (cheap at this scale)
    println("$(row.status): $(row.n_contigs) contigs, NGA50=$(row.NGA50), " *
            "GF=$(row.genome_fraction)%, $(round(row.wall_seconds; digits = 1))s")
end

println("\n--- Phase 3: aggregate + power analysis ---")
results_df = write_aggregate(OUTPUT_DIR, rows)
cv_df = write_power_analysis(OUTPUT_DIR, results_df)

if HAVE_ARTIFACT_WRITER
    try
        write_benchmark_artifacts(
            ["track_a_results" => results_df, "track_a_power_analysis_cv" => cv_df];
            output_dir = joinpath(OUTPUT_DIR, "artifacts"),
            run_id = "track_a_baseline_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))",
            scale = SMOKE ? "smoke" : "full",
            dataset_ids = [o[2] for o in organisms],
            command_args = ARGS,
            metadata = Dict("benchmark" => "track_a_baseline", "k" => K,
                "cv_threshold" => CV_THRESHOLD)
        )
    catch e
        @warn "provenance artifact bundle failed (results TSVs still written)" exception = e
    end
end

n_ok = count(r -> r.status == "ok", rows)
println("\nDone: $(length(rows)) cells, $n_ok ok. Results in $OUTPUT_DIR")
println("End: $(Dates.now())")

end  # PROGRAM_FILE guard
