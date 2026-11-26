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
# READ THIS BEFORE INTERPRETING THE TWO ARMS AS A COMPARISON (td-4e19d.2).
# The two arms are NOT two decoders as far as assembly quality is concerned. They take
# genuinely different code paths, but under the DEFAULT traversal_weighting=:evidence the
# qualmer arm's per-base quality reaches contig EMISSION only and never a traversal
# decision, while contig extraction itself (find_eulerian_paths_next / find_contigs_next)
# is purely topological and shared with the k-mer arm. The arms therefore produce
# IDENTICAL n_contigs / NGA50 / misassemblies / genome_fraction / duplication_ratio /
# largest_contig by construction; they differ only in wall_seconds and peak_rss_bytes.
#
# This is not a conjecture about the code. In the committed 2026-07-24 pilot table
# (track_a_pilot_results_20260724.tsv) all 66 paired runs agree byte-for-byte on every
# assembly metric, and `benchmarking/qualmer_quality_channel_probe.jl` reproduces the
# same invariance from first principles in 18/18 cells across all three chemistries by
# re-assembling identical reads under four different quality vectors.
#
# Consequences for anyone using this harness:
#   - Do NOT report a qualmer-vs-kmer difference in assembly quality; there is none to
#     find, and an apparent one indicates a harness bug rather than a decoder effect.
#   - The arms ARE a valid comparison of COST (runtime, memory) for carrying quality.
#   - To benchmark a genuinely quality-dependent decoder, run this harness with
#     `--traversal-weighting quality` (opt-in; it changes assembly output, so its rows
#     are NOT comparable to any :evidence row), or use the corrector/Viterbi route.
# See benchmarking/results/qualmer_quality_channel_probe/README.md for the full scoping.
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
#   julia --project=. benchmarking/track_a_baseline_benchmark.jl --traversal-weighting quality
#
# --k changes the per-cell checkpoint namespace for any k != 31 (see cell_id_for), so a
# k-sweep and the k=31 baseline can share one --output-dir without colliding.
#
# --traversal-weighting behaves the same way: `quality` namespaces the checkpoint AND is
# recorded as a row column, so an :evidence baseline and a :quality run can share one
# results tree without either resuming the other's cells or being pooled into the other's
# NGA50 CV. The two are different measurements of different assemblers, not replicates.
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
        iseven(parsed) &&
            error("--k must be odd (canonical k-mers require it), got $(parsed)")
        parsed
    end
end
const CV_THRESHOLD = 0.15  # assumed NGA50 coefficient of variation in the power analysis

# Traversal weighting, propagated to `Mycelia.Rhizomorph.assemble_genome`. Without this
# flag the header's instruction to "pass traversal_weighting=:quality" named something no
# invocation of this script could reach — the assemble call took no such keyword and the
# CLI exposed no option for it. Flagged in PR #453 review.
#
# Parsed with the same strictness as --k, for the same reason: the value is part of the
# cell id, so a bare trailing `--traversal-weighting` silently falling back to `evidence`
# would make a requested :quality run resume and republish the :evidence tree.
const TRAVERSAL_WEIGHTINGS = ("evidence", "quality")
const TRAVERSAL_WEIGHTING = let v = findfirst(==("--traversal-weighting"), ARGS)
    if v === nothing
        "evidence"
    elseif v == length(ARGS)
        error("--traversal-weighting requires a value (got a bare trailing flag). " *
              "Refusing to default to evidence: the weighting is part of the cell id, " *
              "so defaulting would silently resume the :evidence tree.")
    else
        raw = ARGS[v + 1]
        if !(raw in TRAVERSAL_WEIGHTINGS)
            allowed = join(TRAVERSAL_WEIGHTINGS, ", ")
            error("--traversal-weighting expects one of $(allowed), got $(repr(raw))")
        end
        raw
    end
end

# Resolve Mycelia checkout provenance stamped into every emitted row.
# Same source as `benchmark_artifacts.jl` (`git rev-parse HEAD` + dirty via
# `git status --porcelain` at the repo root). Empty string when git is
# unavailable — provenance-only, never a measurement. stderr is redirected so a
# non-repo cwd does not spam `fatal: not a git repository` into harness logs.
function resolve_benchmark_git_commit(repo_root::AbstractString = dirname(@__DIR__))
    try
        return chomp(read(pipeline(`git -C $(repo_root) rev-parse HEAD`;
            stderr = Base.devnull), String))
    catch
        return ""
    end
end

function resolve_benchmark_git_dirty(repo_root::AbstractString = dirname(@__DIR__))
    try
        # Exclude `benchmarking/results`: default `--output-dir` is
        # `benchmarking/results/track_a_baseline`, and knowledge-tier files there
        # (track_a_results.tsv, power_analysis_*, artifacts/) are tracked. Without
        # the exclude, a prior run stamps dirty=true on an unchanged source tree and
        # causes spurious cross-host merge collisions on `benchmark_git_dirty`.
        status = chomp(read(pipeline(
            `git -C $(repo_root) status --porcelain -- . ':(exclude)benchmarking/results'`;
            stderr = Base.devnull), String))
        return isempty(status) ? "false" : "true"
    catch
        return ""
    end
end

const BENCHMARK_GIT_COMMIT = resolve_benchmark_git_commit()
const BENCHMARK_GIT_DIRTY = resolve_benchmark_git_dirty()

# Canonical row schema (fixed order so in-memory and JSON-reloaded rows align in the DataFrame).
const ROW_KEYS = (
    :organism, :accession, :technology, :coverage, :seed, :decoder_arm, :k,
    :traversal_weighting,
    :n_reads, :n_contigs, :NGA50, :misassemblies, :genome_fraction,
    :duplication_ratio, :largest_contig, :wall_seconds, :peak_rss_bytes,
    :rss_baseline_bytes, :peak_rss_method, :benchmark_git_commit,
    :benchmark_git_dirty, :status
)
const INT_KEYS = (
    :coverage, :seed, :k, :n_reads, :n_contigs, :largest_contig, :peak_rss_bytes,
    :rss_baseline_bytes)
const FLOAT_KEYS = (
    :NGA50, :misassemblies, :genome_fraction, :duplication_ratio, :wall_seconds)
# `benchmark_git_dirty` is a String ("true"/"false"/""), not Bool: STR_KEYS + legacy
# absent sentinel "" need one type. Differs from benchmark_artifacts.jl's Bool/nothing
# `is_dirty` on purpose.
const STR_KEYS = (
    :organism, :accession, :technology, :decoder_arm, :traversal_weighting,
    :peak_rss_method, :benchmark_git_commit, :benchmark_git_dirty, :status)

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
        traversal_weighting = TRAVERSAL_WEIGHTING,
        n_reads = Int(n_reads), n_contigs = Int(n_contigs),
        NGA50 = Float64(metrics.NGA50), misassemblies = Float64(metrics.misassemblies),
        genome_fraction = Float64(metrics.genome_fraction),
        duplication_ratio = Float64(metrics.duplication_ratio),
        largest_contig = Int(round(Float64(metrics.largest_contig))),
        wall_seconds = round(Float64(wall_seconds); digits = 3),
        peak_rss_bytes = Int(peak_rss_bytes),
        rss_baseline_bytes = Int(rss_baseline_bytes),
        peak_rss_method = String(peak_rss_method),
        benchmark_git_commit = BENCHMARK_GIT_COMMIT,
        benchmark_git_dirty = BENCHMARK_GIT_DIRTY,
        status = String(status)
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
#
# `traversal_weighting` is the one entry here that IS grouped on, and it earns the
# exemption on the same ground the k=31 cell-id name does: the option did not exist when
# any legacy checkpoint was written, so "evidence" is a FACT about how those cells were
# produced rather than a guess standing in for a lost measurement. Every checkpoint
# written from here on carries the key explicitly, so the default can only ever apply to
# a pre-option tree. That is the bar for adding a grouping key to this table — a default
# that is provably the value, not merely a plausible one.
#
# `benchmark_git_commit` / `benchmark_git_dirty` are provenance-only: legacy
# checkpoints predate the columns, so an empty string means "not recorded" rather
# than inventing a SHA or dirty bit for a past run. Dirty is "true"/"false" when
# known — a SHA alone overstates provenance on an uncommitted tree.
const OPTIONAL_KEY_DEFAULTS = Dict{Symbol, Any}(
    :peak_rss_method => "unknown", :rss_baseline_bytes => -1,
    :traversal_weighting => "evidence", :benchmark_git_commit => "",
    :benchmark_git_dirty => "")
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
# current_rss_bytes / timed_with_peak_rss live in a shared include (the D5 grid
# driver reuses them); see that file for the measured rationale.
include(joinpath(@__DIR__, "peak_rss_sampler.jl"))

# === Per-cell execution ===

function run_cell(org, acc, ref, tech, cov, seed, arm, cell_dir)
    records = simulate_reads(tech, ref, cov, seed, joinpath(cell_dir, "reads"))
    n_reads = length(records)
    asm_input = shape_for_arm(records, arm)

    GC.gc()
    measured = timed_with_peak_rss(
        () -> Mycelia.Rhizomorph.assemble_genome(asm_input; k = K, verbose = false,
        traversal_weighting = Symbol(TRAVERSAL_WEIGHTING)))
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
# Unique-temp rename lives in atomic_publish.jl so the merge script calls the same
# function. This file cannot be included from there (it imports Mycelia and reads ARGS).
include(joinpath(@__DIR__, "atomic_publish.jl"))

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
    # :traversal_weighting joins :k as a grouping key for the same reason — an :evidence
    # cell and a :quality cell are different assemblers on the same inputs, not two draws
    # from one distribution, so pooling them reports a between-ASSEMBLER spread as the
    # between-seed CV the pre-registration is about.
    #
    # A table assembled outside `canonical` (a hand-built frame, or a results TSV written
    # before the option existed) has no such column. Fill it with the same value
    # OPTIONAL_KEY_DEFAULTS uses rather than erroring: pre-option rows are :evidence by
    # definition, because the keyword had no CLI surface here. `df` is already a copy
    # from the `measured` filter above, so this does not mutate the caller's frame.
    hasproperty(df, :traversal_weighting) ||
        (df[!, :traversal_weighting] = fill("evidence", DataFrames.nrow(df)))
    for g in DataFrames.groupby(df,
        [:organism, :technology, :coverage, :decoder_arm, :k, :traversal_weighting])
        nga = Float64.(g.NGA50)
        m = Statistics.mean(nga)
        s = length(nga) > 1 ? Statistics.std(nga; corrected = true) : NaN
        cv = (m == 0 || isnan(s)) ? NaN : s / m
        push!(cv_rows,
            (
                organism = g.organism[1], technology = g.technology[1],
                coverage = g.coverage[1], decoder_arm = g.decoder_arm[1],
                traversal_weighting = g.traversal_weighting[1],
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
#
# `traversal_weighting` is namespaced on exactly the same terms and for exactly the same
# reason. A :quality run assembles a DIFFERENT graph traversal from the :evidence
# baseline, so sharing a cell id would let one resume and republish the other's contigs.
# The suffix is again applied only for the non-default value, so every checkpoint written
# before the option existed keeps its historical name — legacy trees are :evidence by
# definition, because the keyword had no CLI surface here at all.
function cell_id_for(org, tech, cov, seed, arm; k = K,
        traversal_weighting = TRAVERSAL_WEIGHTING)
    base = k == 31 ? "$(org)__$(tech)__$(cov)x__seed$(seed)__$(arm)" :
           "$(org)__$(tech)__$(cov)x__seed$(seed)__$(arm)__k$(k)"
    return traversal_weighting == "evidence" ? base : "$(base)__w$(traversal_weighting)"
end

# The row recorded when a cell throws. Extracted from the catch block because it was
# unreachable from any test there: adding a required `peak_rss_method` keyword to
# cell_row without updating this call made the ERROR HANDLER ITSELF throw
# UndefKeywordError, so the first failing cell aborted the whole matrix from inside
# its own recovery path — and `--smoke` runs a single deliberately-successful cell,
# so no test could reach it.
function error_row(org, acc, tech, cov, seed, arm)
    cell_row(org, acc, tech, cov, seed, arm;
        n_reads = 0, n_contigs = 0, wall_seconds = 0.0, peak_rss_bytes = 0,
        rss_baseline_bytes = -1, peak_rss_method = "unknown",
        metrics = empty_metrics(), status = "error")
end

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
    println("Traversal weighting: $TRAVERSAL_WEIGHTING")
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
