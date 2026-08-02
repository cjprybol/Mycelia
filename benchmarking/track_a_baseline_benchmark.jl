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
# Each cell: simulate reads -> assemble (k = --k, default 31) -> QUAST vs reference -> parse NGA50 /
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
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --k 19       # k-mer size; positive ODD integer, default 31
#
# --k changes the per-cell checkpoint namespace for any k != 31 (see cell_id_for), so a
# k-sweep and the k=31 baseline can share one --output-dir without colliding.
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
# k is settable so an ONT k-sweep can probe whether ANY k assembles high-error long
# reads usefully (at k=31 every ONT cell at 10x/30x yields NGA50 = 0). Default unchanged.
#
# A bare trailing `--k` must ERROR, not fall back to the default. Because k is part of
# the cell id, a silent fallback to 31 would make a `--k` run resume and republish the
# COMPLETED k=31 tree as if it were the requested k — the exact wrong-answer path the
# k-in-cell-id change exists to close, reintroduced through the argument parser.
const K = let v = findfirst(==("--k"), ARGS)
    if v === nothing
        31
    elseif v == length(ARGS)
        error("--k requires a value (got a bare trailing --k). Refusing to default to 31: " *
              "k is part of the cell id, so defaulting would silently resume the k=31 tree.")
    else
        raw = ARGS[v + 1]
        parsed = tryparse(Int, raw)
        parsed === nothing && error("--k expects an integer, got $(repr(raw))")
        parsed < 1 && error("--k must be >= 1, got $(parsed)")
        iseven(parsed) && error("--k must be odd (canonical k-mers require it), got $(parsed)")
        parsed
    end
end
const CV_THRESHOLD = 0.15  # assumed NGA50 coefficient of variation in the power analysis

# Canonical row schema (fixed order so in-memory and JSON-reloaded rows align in the DataFrame).
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
const STR_KEYS = (:organism, :accession, :technology, :decoder_arm, :peak_rss_method, :status)

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
        wall_seconds, peak_rss_bytes, rss_baseline_bytes, peak_rss_method, metrics, status)
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

# Keys whose absence from an on-disk checkpoint is genuinely benign: they are
# PROVENANCE about how a cell was measured, added after some checkpoints were
# written, and no analysis groups or averages on them.
#
# Everything NOT in this set is a scientific measurement or a grouping key, and a
# missing one must fail loudly. Zero-filling those is not a safe default here: in
# this dataset `NGA50 = 0.0` and `genome_fraction = 0.0` are the GENUINE values for
# every ONT cell at 10x/30x, and `peak_rss_bytes = 0` was genuine for 271/288 cells
# in the pre-fix run — so a zero-filled corrupt checkpoint is byte-identical to a
# real result and would be pooled straight into the CV that gates the
# pre-registration. A truncated checkpoint must not be able to move a verdict.
# Key => sentinel. Driving the defaulting from this table (rather than hardcoding the
# two cases in canonical) is what makes the error message below TRUE: adding a key here
# really does make it defaultable. Previously OPTIONAL_KEYS was declared but never read,
# so a maintainer following the message's advice would still have hit the hard error.
const OPTIONAL_KEY_DEFAULTS = Dict{Symbol, Any}(
    :peak_rss_method => "unknown", :rss_baseline_bytes => -1)
const OPTIONAL_KEYS = Tuple(sort(collect(keys(OPTIONAL_KEY_DEFAULTS))))

# Rebuild a canonical, type-coerced NamedTuple from a parsed JSON dict (resumed cells).
#
# Absent OPTIONAL_KEYS take an explicit "unknown"/-1 sentinel that no real
# measurement can produce. Any other absent key raises: indexing `d[String(key)]`
# directly raised a bare KeyError naming only the key, which crashed the run on the
# first cached cell with no indication that the tree simply predated a schema change.
function canonical(d::AbstractDict)
    vals = map(ROW_KEYS) do key
        name = String(key)
        if !haskey(d, name)
            haskey(OPTIONAL_KEY_DEFAULTS, key) && return OPTIONAL_KEY_DEFAULTS[key]
            error("checkpoint is missing required key $(name). This is a measurement or " *
                  "grouping column, so it cannot be defaulted — a zero here is " *
                  "indistinguishable from a real result and would silently enter the " *
                  "power analysis. Delete the checkpoint to recompute the cell, or add " *
                  "$(name) to OPTIONAL_KEYS if it is genuinely provenance-only.")
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

# Write via temp + rename. A SIGKILL (SLURM timeout, OOM-killer) mid-write otherwise
# leaves a truncated checkpoint that `isfile` accepts as a valid cache on the next
# run. Same-filesystem rename is atomic, so a checkpoint is either absent or complete.
function save_cell_json(path, row)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        JSON.print(io, Dict(string(k) => v for (k, v) in pairs(row)), 2)
    end
    mv(tmp, path; force = true)
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

# === Per-cell memory measurement ===

# `Sys.maxrss()` is the process-lifetime high-water mark (ru_maxrss) and is
# monotonically non-decreasing, so the old `max(0, Sys.maxrss() - rss0)` could only
# be nonzero when a cell beat the all-time PROCESS peak. It measured "excess over
# historical peak", not per-cell peak.
#
# Measured on the 2026-07-25 Lovelace run (288 cells): 271/288 rows (94%) were
# exactly 0 — 144/144 kmer AND 127/144 qualmer. The effect is a global ratchet, not
# a per-pair alternation: all 17 nonzero rows are qualmer only because qualmer is
# the inner loop's first arm, they sit in the first 48% of the run, and the ratchet
# closes for good at T4/ont/100x/seed42 = 51.3 GiB — after which nothing beats it.
# The 6% that survived understate true usage by orders of magnitude (e.g.
# T4/illumina/100x/seed456 = 11 MiB for a multi-GiB cell), because each is a
# delta over the previous all-time peak rather than a measurement. LRC's 144-cell
# run agrees: 131/144 zero (72/72 kmer, 59/72 qualmer).
#
# Read VmRSS from /proc/self/status instead: it reports explicit kB, so unlike
# /proc/self/statm it needs no Sys.PAGESIZE (which does not exist in Julia 1.10).
function current_rss_bytes()
    status_path = "/proc/self/status"
    isfile(status_path) || return nothing
    # `open(...) do` rather than bare `eachline(path)`: eachline closes its handle at
    # end-of-iteration, and VmRSS sits near the top of the file so the early return
    # below ALWAYS exits first. Measured: 2000 polls leaked 2000 descriptors, released
    # only by GC. At 20 Hz that is ~72k/hour, and EMFILE would surface inside the
    # sampler task where it is hardest to diagnose.
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
`julia track_a_baseline_benchmark.jl` reverts to the broken high-water semantics; launch
with `julia -t N` or the column is meaningless. That revert is no longer silent — a
one-time `@warn` fires from the fallback branch below whenever it happens.

`method` records how the number was obtained, so a reader never has to infer it and rows
measured different ways are never averaged together:
  "sampled"          - process-RSS high-water, polled while the cell ran, with EVERY poll
                       succeeding. A flat cell whose RSS never exceeds its baseline is
                       still "sampled" — that is a real observation.
  "sampled-partial"  - at least one poll succeeded and at least one THREW. Sampling
                       continued, but the peak may have been missed while the probe was
                       failing, so the number is a lower bound, not a measurement.
  "sampled-degraded" - the sampler never completed a single read (starved, or every /proc
                       read failed); the number is the entry baseline, NOT a measurement.
  "highwater-delta"  - single-threaded or no /proc: the old, known-broken semantics
  "unknown"          - checkpoint predates this column (see OPTIONAL_KEY_DEFAULTS)
Values under different methods are different quantities. Always filter on
`peak_rss_method` before aggregating.

`probe` is the RSS reader, injectable ONLY so the failure labelling above can be tested:
`current_rss_bytes` returns `nothing` on macOS and there is no portable way to make a
real /proc read fail on the third poll.
"""
function timed_with_peak_rss(f; interval_seconds::Float64 = 0.05,
        probe = current_rss_bytes)
    # Guard the baseline read too: everything here is careful never to let telemetry fail
    # science, and an unguarded /proc read error would abort the cell before f() even ran.
    baseline = try
        probe()
    catch e
        @warn "baseline RSS read failed; falling back to high-water semantics" exception = e
        nothing
    end

    can_sample = baseline !== nothing && Threads.nthreads() > 1
    if !can_sample
        # The docstring says a bare `julia` invocation "silently reverts to the broken
        # high-water semantics". `can_sample == false` knows that at runtime, so say it
        # instead of leaving the operator to infer it from a column of near-zeros after a
        # multi-day run. `maxlog` keeps it to one line, not one per cell.
        if Threads.nthreads() == 1
            @warn "peak_rss_bytes is falling back to the known-broken high-water " *
                  "semantics: Julia has ONE thread, so the RSS sampler cannot be " *
                  "scheduled against a CPU-bound assembly. Relaunch with `julia -t N` " *
                  "(both sbatch wrappers export JULIA_NUM_THREADS) or treat the column " *
                  "as unusable." maxlog = 1
        elseif baseline === nothing
            @warn "peak_rss_bytes is falling back to the known-broken high-water " *
                  "semantics: /proc/self/status is unreadable on this host." maxlog = 1
        end
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
    # Count FAILED reads separately, because "never observed anything" and "stopped
    # observing partway" are different lies. With the try OUTSIDE the loop, the first
    # throwing poll ended sampling for the rest of the cell, yet `reads > 0` still
    # labelled the frozen partial maximum "sampled" — a plausible number presented as a
    # measurement, which is the exact failure mode the read count was added to catch.
    failed_reads = Threads.Atomic{Int}(0)
    sampler = Threads.@spawn begin
        # The sampler must never be able to fail a good cell: a transient /proc read error
        # is telemetry, not science. The try is INSIDE the loop, so one bad poll costs one
        # sample rather than the rest of the cell, and the counts below report it.
        try
            while keep_sampling[]
                observed = try
                    probe()
                catch
                    Threads.atomic_add!(failed_reads, 1)
                    nothing
                end
                if observed !== nothing
                    Threads.atomic_add!(reads, 1)
                    observed > peak[] && Threads.atomic_max!(peak, observed)
                end
                sleep(interval_seconds)
            end
        catch
            # Last-resort net for a failure OUTSIDE the probe (an interrupted sleep, say).
            # Deliberately swallowed; the two read counts report the degradation.
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
    # Read the counters only after the `finally` above has reaped the sampler, so the
    # label describes the whole cell rather than a racing snapshot of it.
    method = if reads[] == 0
        "sampled-degraded"
    elseif failed_reads[] > 0
        "sampled-partial"
    else
        "sampled"
    end
    return (value = timed.value, wall_seconds = timed.time,
        peak_rss_bytes = peak[], rss_baseline_bytes = baseline, method = method)
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
    peak_rss_method = measured.method
    rss_baseline_bytes = measured.rss_baseline_bytes

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

# Rebuild the aggregate by re-reading EVERY checkpoint under root/cells, not from the
# in-memory rows of this invocation.
#
# `rows` holds only the cells in the current matrix. Since the file is truncate-written,
# any narrower run — a shard (`--organisms Lambda`), a re-run, or another `--k` — used to
# replace a complete table with its own slice. That is not hypothetical: six `--k`
# invocations against one --output-dir left a 6-row TSV over 36 completed cells (83%
# of the sweep), and regenerated power_analysis_summary.md to describe the surviving
# sixth. The per-cell store was correct the whole time; only the published layer on
# top of it lost data.
#
# Re-reading the store makes the aggregate monotone and idempotent by construction:
# shards, sweeps and resumes all converge on the union rather than the last writer's
# slice, and a no-op re-run repairs a previously truncated table.
# Atomic publish with a PER-PROCESS temp name.
#
# A FIXED temp path (`target * ".tmp"`) is safe for a single writer and unsafe here:
# write_aggregate runs once per cell, and sharding is the documented workflow — the
# harness header advertises that "an HPC array job can split the matrix and share one
# results tree", and run_track_a_baseline_lrc.sbatch shows a shard invocation with no
# --output-dir, so every shard defaults to the same directory. Two shards then open
# the SAME temp inode and both rename it, and the rename publishes interleaved content
# as a complete table. That is worse than the truncation bug the store-rebuild fixed:
# a short table is detectable, a corrupt one is not.
#
# tempname(dir) is unique per call, so concurrent shards can no longer collide, and the
# rename stays atomic because the temp is on the same filesystem as the target.
function _publish_atomically(write!, target::AbstractString)
    dir = dirname(target)
    mkpath(dir)
    tmp = tempname(dir)
    try
        write!(tmp)
        mv(tmp, target; force = true)
    catch
        isfile(tmp) && rm(tmp; force = true)
        rethrow()
    end
    return target
end

function write_aggregate(root, rows)
    cells_dir = joinpath(root, "cells")
    all_rows = NamedTuple[]
    if isdir(cells_dir)
        skipped = String[]
        for entry in sort(readdir(cells_dir))
            ckpt = joinpath(cells_dir, entry, "cell_result.json")
            isfile(ckpt) || continue
            # Guard per checkpoint. canonical() now raises on a missing measurement
            # (deliberately — a zero-filled corrupt cell is indistinguishable from a
            # real ONT failure), and this loop touches every checkpoint in the tree, so
            # ONE stale or truncated file would otherwise abort the aggregate and with
            # it the whole run. Skip it loudly and keep the other cells.
            try
                push!(all_rows, canonical(JSON.parsefile(ckpt)))
            catch e
                push!(skipped, entry)
                @warn "skipping unreadable checkpoint in aggregate" entry exception = e
            end
        end
        isempty(skipped) ||
            @warn "aggregate omits $(length(skipped)) unreadable checkpoint(s); " *
                  "delete them to force recompute" skipped
    end
    # Fall back to the in-memory rows only when the store is unreadable, so an aggregate
    # is still produced rather than silently emptied.
    isempty(all_rows) && (all_rows = collect(rows))
    df = DataFrames.DataFrame(all_rows)
    target = joinpath(root, "track_a_results.tsv")
    _publish_atomically(tmp -> CSV.write(tmp, df; delim = '\t'), target)
    return df
end

function write_power_analysis(root, df)
    cv_rows = NamedTuple[]
    # Group by :k as well. With --k in play a single tree can hold several k, and
    # pooling them turns a between-SEED CV into a between-K CV — inflating it by up to
    # ~39x in testing and flipping the pre-registration verdict with no warning.
    #
    # Exclude rows that are not measurements. An error row is not a measurement and
    # must not enter a variance estimate.
    #
    # `status` alone is insufficient, and this is the second door into the same bug:
    # on a QUAST exception run_cell substitutes empty_metrics() (NGA50 and
    # largest_contig both 0.0) and then derives status from n_contigs ALONE, so a
    # cell QUAST never scored is recorded as "ok" carrying a full row of zeros.
    # Those zeros would be pooled into the CV as if measured.
    #
    # `n_contigs > 0 && largest_contig == 0` is reachable ONLY through
    # empty_metrics(), so every row it matches is a non-measurement — never a real
    # assembly that merely scored badly. n_contigs is the UNFILTERED assembly count
    # (bench-side, before QUAST), while QUAST runs with min_contig = 500 and refuses
    # to emit a report at all when nothing clears it. Verified against the 2026-07-25
    # Lovelace run: both matching cells (phi29/ont/10x/seed456, 6331 contigs each)
    # have NO quast/report.tsv, and their quast.log ends
    #     WARNING: Skipping contigs because it doesn't contain contigs >= 500 bp.
    #     ERROR! None of the assembly files contains correct contigs.
    # i.e. QUAST exited non-zero, run_cell caught it, and empty_metrics() supplied the
    # zeros. The min-contig threshold is the CAUSE of that error, not evidence that a
    # measurement happened. track_a_merge_hosts.jl states the same implication.
    #
    # So the exclusion needs no judgement call about a borderline class: every row it
    # matches is provably a non-measurement.
    #
    # It IS one-directional on the verdict, though, and not through the CV magnitude —
    # through the DENOMINATOR. Dropping a row shrinks its group's n; at n == 1 the
    # sd is NaN, the cv is NaN, and the group leaves `evaluable` entirely, so a group
    # that would have FAILED the threshold is removed from `n_pass == n_eval` rather
    # than counted against it. Measured on a two-group fixture, varying only the
    # exclusion: with it, 1 evaluable cell and "supported"; without it, 2 evaluable
    # cells, max CV 1.41 and "NOT fully supported". The verdict flips toward
    # "supported" — the favourable direction for the pre-registration this gates.
    #
    # Excluding non-measurements is still right; a reader just has to be able to see
    # how much of the matrix went unscored, which is why n_excluded is printed in the
    # summary artifact next to the counts rather than left in a @warn.
    #
    # Same implication track_a_merge_hosts.jl's `quast_evidence` uses, and the same
    # predicate as ok_cells() in track_a_harvest_figures.jl — the figures and this
    # power analysis must not disagree about which rows are measurements.
    measured = (df.status .== "ok") .&
               .!((df.n_contigs .> 0) .& (df.largest_contig .== 0))
    n_excluded = DataFrames.nrow(df) - count(measured)
    n_excluded > 0 &&
        @warn "power analysis excludes non-ok and QUAST-unscored cells" n_excluded
    df = df[measured, :]
    for g in DataFrames.groupby(df, [:organism, :technology, :coverage, :decoder_arm, :k])
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
    # Same publish discipline as the aggregate: these two derived files were bare
    # truncating writes, so a concurrent shard (or a crash mid-write) could leave a
    # half-written verdict beside a complete results table.
    _publish_atomically(tmp -> CSV.write(tmp, cv_df; delim = '\t'),
        joinpath(root, "power_analysis_cv.tsv"))

    evaluable = filter(r -> isfinite(r.cv_nga50), cv_rows)
    n_eval = length(evaluable)
    n_pass = count(r -> r.passes, evaluable)
    max_cv = isempty(evaluable) ? NaN : maximum(r.cv_nga50 for r in evaluable)
    verdict = (n_eval > 0 && n_pass == n_eval) ? "supported" :
              n_eval == 0 ? "indeterminate (no evaluable cells)" : "NOT fully supported"

    _publish_atomically(joinpath(root, "power_analysis_summary.md")) do path
        open(path, "w") do io
            println(io, "# Track A power-analysis check — NGA50 CV vs assumed $(CV_THRESHOLD)\n")
            println(io, "Generated: $(Dates.now())\n")
            println(io, "- Evaluable cells (organism×tech×coverage×arm, ≥2 seeds, nonzero NGA50): $n_eval")
            println(io, "- Cells with CV ≤ $(CV_THRESHOLD): $n_pass / $n_eval")
            println(io, "- Max CV observed: $(round(max_cv; digits = 4))")
            println(io,
                "- Rows excluded as non-measurements (status != ok, or QUAST " *
                "unscored): $n_excluded")
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
    end
    return cv_df
end

# === Cell identity + error rows (extracted so both are unit-testable) ===

# k is part of the cell id: without it, two runs at different --k write the same
# cells/<id>/cell_result.json and the second silently resumes the first's result.
#
# But the suffix is applied ONLY for non-default k. Interpolating it unconditionally
# renamed the DEFAULT k=31 cells too, which orphaned every checkpoint written before
# the flag existed — 288 on Lovelace and 144 on LRC — turning the sbatch wrappers'
# documented "just re-submit, completed cells are skipped" recovery into a silent
# full recompute. Legacy trees are k=31 by definition (k was a hardcoded const), so
# keying k=31 to the historical name is unambiguous and needs no migration.
cell_id_for(org, tech, cov, seed, arm; k = K) =
    k == 31 ? "$(org)__$(tech)__$(cov)x__seed$(seed)__$(arm)" :
    "$(org)__$(tech)__$(cov)x__seed$(seed)__$(arm)__k$(k)"

# The row recorded when a cell throws. Extracted from the catch block because it was
# unreachable from any test there: adding a required `peak_rss_method` keyword to
# cell_row without updating this call made the ERROR HANDLER ITSELF throw
# UndefKeywordError, so the first failing cell aborted the whole matrix from inside
# its own recovery path — and `--smoke` runs a single deliberately-successful cell,
# so no test could reach it.
error_row(org, acc, tech, cov, seed, arm) =
    cell_row(org, acc, tech, cov, seed, arm;
        n_reads = 0, n_contigs = 0, wall_seconds = 0.0, peak_rss_bytes = 0,
        rss_baseline_bytes = -1, peak_rss_method = "unknown",
        metrics = empty_metrics(), status = "error")

# === Main ===

# Guard the driver so a test can `include()` this file for its pure helpers without
# downloading genomes or running assemblies. Matches the convention already used by
# rhizomorph_benchmark_harness.jl, viterbi_accuracy_benchmark.jl and others, and is
# what makes test/4_assembly/track_a_baseline_benchmark_test.jl possible. `if` does
# not introduce scope at top level, so the globals below stay global.
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
    # k belongs in the cell id: without it, two runs at different --k write the same
    # cells/<id>/cell_result.json and the second silently resumes the first's result.
    # Existing k=31 trees keep their old names, so point a k-sweep at its own
    # --output-dir rather than reusing a completed one.
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

const ARTIFACT_RUN_ID = "track_a_baseline_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))"
if HAVE_ARTIFACT_WRITER
    try
        write_benchmark_artifacts(
            ["track_a_results" => results_df, "track_a_power_analysis_cv" => cv_df];
            # Per-run subdirectory: a fixed "artifacts" path meant six --k invocations
            # left ONE bundle whose provenance (git SHA, tool versions, command_args)
            # described only the last run. Unlike the cell JSONs there is no second
            # copy, so that loss is unrecoverable.
            output_dir = joinpath(OUTPUT_DIR, "artifacts", ARTIFACT_RUN_ID),
            run_id = ARTIFACT_RUN_ID,
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
