# ONT k-selection sweep — is the pilot's ONT degeneracy a k artifact? (td-4e19d.28)
#
# WHY THIS EXISTS
# ---------------
# The Track-A greedy-baseline pilot (rhizomorph-paper, 2026-07-24) ran 132 cells
# and its ONT cells are largely degenerate: Lambda and T4 both return NGA50 = 0
# in every replicate at 10x and 30x. Two facts about that pilot make the obvious
# explanation ("ONT is too noisy for naked k-mer assembly") non-obvious:
#
#   1. All 132 pilot rows used k = 31 — the Illumina-tuned primary k — including
#      all 36 ONT rows. k was a hardcoded `const` in that harness.
#   2. ONT reads came from `Mycelia.simulate_nanopore_reads`, which invokes
#      Badread with NO --error_model / --qscore_model / --identity flags
#      (src/simulation.jl:1091). It therefore inherits the INSTALLED Badread's
#      compiled-in defaults, which for v0.4.1 are nanopore2023 / identity
#      95,99,2.5 — a modern, relatively low-error chemistry, not R9.
#
# At a per-base error rate e = 0.05, P(error-free k-mer) = (1-e)^k = 0.204 for
# k = 31, so 30x raw coverage still carries ~6.1x error-free 31-mer coverage.
# Degeneracy is therefore NOT arithmetically forced at 30x, which is what makes
# the pilot's result worth resolving rather than explaining away.
#
# HYPOTHESES UNDER TEST
#   H-a  k-selection artifact — k = 31 is simply too long for this error rate,
#        and the degeneracy dissolves into a per-chemistry k choice.
#   H-b  degeneracy persists at every appropriate k — a genuine finding about
#        naked k-mer assembly on long reads.
#   H-c  an assembler defect that manifests at low error rates — a bug, and the
#        most consequential outcome.
#
# DESIGN
# ------
# Lambda (NC_001416) x {ont, illumina} x k in {11,15,21,31} x coverage
# {10,30,50,100} x seeds {42,123,456} = 96 cells.
#
#   * SINGLE k PER CELL. The pilot's raw baseline used a single k while the
#     iterative corrector walks a k-ladder; running a ladder here would confound
#     k with correction. `corrector` is left at its default `:none`, which
#     `assemble_genome` documents as the byte-identical single-k pipeline
#     (src/rhizomorph/assembly.jl:1083).
#   * ILLUMINA CONTROL over the SAME k ladder and coverages, so a k effect can
#     be separated from a chemistry effect. Illumina and ONT are different error
#     processes and are NEVER pooled into an aggregate — every summary in this
#     harness is stratified by technology.
#   * ONE DECODER ARM (`kmer`, quality stripped). The pilot ran both `qualmer`
#     and `kmer` arms and all 66 complete arm-pairs produced identical
#     n_contigs / NGA50 / genome_fraction, so the arm is not a live factor and
#     doubling the grid would buy nothing. The column is still recorded so the
#     choice is visible in the data rather than only in this comment.
#
# CENSORING
# ---------
# QUAST prints "-" for NGA50 when NOTHING aligns. The pilot's parser coerced
# that to 0.0 (`tryparse(Float64, "-") === nothing` -> 0.0), which is why its
# ONT rows read NGA50 = 0. A censored floor and a measured zero are different
# facts, so this harness keeps them apart: `NGA50` holds `NA` when QUAST could
# not compute it, and `nga50_status` records WHY. Every consumer of this table
# must read genome_fraction and contig counts beside NGA50.
#
# Relatedly, QUAST reports "# misassemblies" = 0 on UNALIGNED contigs. A 0 in
# that column on a degenerate cell means non-alignment, not correctness, so
# `misassemblies` is only interpretable where `nga50_status == "measured"`.
#
# Usage:
#   julia --project=. benchmarking/ont_k_sweep.jl                     # full 96-cell grid
#   julia --project=. benchmarking/ont_k_sweep.jl --smoke             # 1 cheap cell
#   julia --project=. benchmarking/ont_k_sweep.jl --technologies illumina
#   julia --project=. benchmarking/ont_k_sweep.jl --ks 11,21 --coverages 10,30
#   julia --project=. benchmarking/ont_k_sweep.jl --output-dir /scratch/ont_k_sweep
#
# Per-cell JSON checkpoints make the run crash-safe and resumable: re-invoking
# with the same --output-dir skips completed cells.

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

# === Configuration ===

# (display name, NCBI accession, expected size in bp). Accessions match the
# Track-A pilot so cells are directly comparable to it.
#
# Lambda is the default: it is the smaller, faster genome and spans both the
# degenerate and non-degenerate regimes. T4 is the generality check — it is 3.5x
# larger with different repeat structure, and repeat structure is what sets BOTH
# ends of the usable k range (the short-k floor where unitigs cannot grow, and
# the saturation ceiling where a short k can no longer resolve repeats). A k
# recommendation derived from one genome does not automatically transfer, so the
# organism has to be a variable rather than a constant.
const ORGANISMS = [
    ("Lambda", "NC_001416", 48_502),
    ("T4", "NC_000866", 168_903)
]

const TECHNOLOGIES = ["illumina", "ont"]
const KS = [11, 15, 21, 31]
const COVERAGES = [10, 30, 50, 100]
const SEEDS = [42, 123, 456]
const DECODER_ARM = "kmer"
# Matches the Track-A pilot's QUAST setting, so cells are directly comparable.
const MIN_CONTIG = 500
# QUAST threads PER CELL. Deliberately small and fixed rather than
# `get_default_threads()`, because this harness is designed to be run as many
# concurrent single-cell shards: on a 112-core host, 32 shards each taking
# QUAST's default of 16 threads would ask for 512. Assembly is single-threaded
# Julia and dominates the runtime, so QUAST's thread count is a scheduling
# detail — it changes no metric, only how much the shards fight each other.
const QUAST_THREADS = 2

# --- Outcome classification -------------------------------------------------
#
# THE ONE JUDGMENT CALL IN THIS HARNESS. Step 5 of the task asks whether there
# is a k at which ONT "stops being degenerate", which requires an operational
# definition of degenerate. Several are defensible; this one is deliberately
# stated in exactly one place so it can be changed and the whole table
# reclassified by re-running `classify_outcome` over the committed TSV.
#
# The ladder is anchored on genome fraction (how much of the reference is
# recovered at all) and on NGA50 as a fraction of the genome (whether that
# recovery is contiguous), because those are the two quantities that survive the
# censoring problem above. Contig count is deliberately NOT part of the rule:
# this assembler emits thousands of contigs even on a fully successful cell
# (Illumina/30x/k=31 recovered NGA50 48058 with 5357 contigs), so contig count
# does not separate success from failure here.
const GF_DEGENERATE_MAX = 25.0    # below this, essentially nothing was recovered
const GF_SUBSTANTIAL_MIN = 90.0   # at/above this, the genome is mostly present
const GF_NEAR_COMPLETE_MIN = 95.0
const NGA50_SUBSTANTIAL_FRAC = 0.10   # NGA50 >= 10% of genome
const NGA50_NEAR_COMPLETE_FRAC = 0.50 # NGA50 >= 50% of genome

"""
    classify_outcome(nga50, genome_fraction, nga50_status, genome_size) -> String

Map a scored cell onto one of four ordered outcome tiers.

`genome_size` is an explicit argument rather than a constant because the NGA50
tiers are defined as a FRACTION of the genome: "NGA50 >= 10% of the genome" is
4,850 bp on Lambda and 16,890 bp on T4. Reusing one genome's size to classify
another would silently rescale the tiers and make cross-organism rows
incomparable — which is the whole point of running a second organism.

The ordered tiers are:

  * `degenerate`    — nothing aligned (censored NGA50) or genome fraction below
                      `GF_DEGENERATE_MAX`. This is the pilot's ONT signature.
  * `partial`       — some of the genome is recovered, but either fragmented or
                      incomplete.
  * `substantial`   — genome fraction >= 90% AND NGA50 >= 10% of the genome.
  * `near_complete` — genome fraction >= 95% AND NGA50 >= 50% of the genome.

A censored NGA50 always classifies as `degenerate` regardless of genome
fraction: if QUAST could not compute NGA50, no contig aligned, so any nonzero
genome fraction would be internally inconsistent rather than informative.
"""
function classify_outcome(nga50, genome_fraction, nga50_status::AbstractString,
        genome_size::Real)
    nga50_status == "measured" || return "degenerate"
    gf = Float64(genome_fraction)
    n = Float64(nga50)
    gf < GF_DEGENERATE_MAX && return "degenerate"
    if gf >= GF_NEAR_COMPLETE_MIN && n >= NGA50_NEAR_COMPLETE_FRAC * genome_size
        return "near_complete"
    elseif gf >= GF_SUBSTANTIAL_MIN && n >= NGA50_SUBSTANTIAL_FRAC * genome_size
        return "substantial"
    else
        return "partial"
    end
end

"""
    genome_size_for(organism) -> Int

Expected reference length for `organism`, from the ORGANISMS table.

Errors on an unknown organism rather than defaulting: a wrong genome size does
not fail loudly, it silently rescales the NGA50 outcome tiers and produces a
well-formed table of misclassified cells.
"""
function genome_size_for(organism::AbstractString)
    for (name, _accession, size) in ORGANISMS
        name == organism && return size
    end
    error("unknown organism $(organism); known: " *
          join((o[1] for o in ORGANISMS), ", "))
end

# Canonical row schema. Fixed order so in-memory rows and JSON-reloaded
# checkpoints align in the DataFrame.
#
# The `asm_*` columns are computed FROM THE CONTIGS IN JULIA, never from QUAST.
# They exist because QUAST refuses to score an assembly in which no contig
# reaches --min-contig, and that refusal is the empirical signature of the
# low-k regime: at ONT/30x/k=11 this assembler emits 149,115 contigs whose
# LONGEST is 338 bp. If the only contig statistics came from QUAST, every cell
# in that regime would be blank and the most informative stratum of the sweep
# would read as missing data rather than as the result it is.
const ROW_KEYS = (
    :organism, :accession, :technology, :k, :coverage, :seed, :decoder_arm,
    :n_reads, :n_contigs, :asm_contigs_ge_min, :asm_max_contig, :asm_total_bp,
    :quast_contigs, :total_length, :N50, :largest_contig,
    :NGA50, :nga50_status, :NA50, :largest_alignment, :genome_fraction,
    :duplication_ratio, :misassemblies, :unaligned_contigs, :unaligned_length,
    :outcome, :wall_seconds, :status
)
const INT_KEYS = (:k, :coverage, :seed, :n_reads, :n_contigs,
    :asm_contigs_ge_min, :asm_max_contig, :asm_total_bp)
const STR_KEYS = (:organism, :accession, :technology, :decoder_arm,
    :nga50_status, :outcome, :status)
# Everything else is a Float64-or-missing measurement column: QUAST legitimately
# emits "-" for several of these on a degenerate assembly, and `missing` is the
# honest representation of "QUAST could not compute this", never 0.

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

# Default to Lambda only. T4 is opt-in via --organisms because it is 3.5x larger
# and its ONT high-coverage cells are the most expensive in the grid; a bare
# invocation should not silently commit to them.
organisms = ORGANISMS[1:1]
technologies = TECHNOLOGIES
ks = KS
coverages = COVERAGES
seeds = SEEDS

if SMOKE
    technologies = ["illumina"]
    ks = [31]
    coverages = [10]
    seeds = [42]
else
    # Assign directly, not inside a `let`: a `let` body introduces a new scope,
    # so these assignments would create locals and silently leave the globals at
    # their defaults — i.e. the shard flags would be accepted and then ignored.
    _f = arg_list("--organisms")
    if _f !== nothing
        organisms = filter(o -> o[1] in _f, ORGANISMS)
        # An unrecognised name would otherwise silently yield an empty grid that
        # completes instantly and reports success over zero cells.
        isempty(organisms) && error(
            "no known organism matched --organisms $(join(_f, ",")); known: " *
            join((o[1] for o in ORGANISMS), ", "))
    end
    _f = arg_list("--technologies")
    _f !== nothing && (technologies = _f)
    _f = arg_list("--ks")
    _f !== nothing && (ks = parse.(Int, _f))
    _f = arg_list("--coverages")
    _f !== nothing && (coverages = parse.(Int, _f))
    _f = arg_list("--seeds")
    _f !== nothing && (seeds = parse.(Int, _f))
end

const OUTPUT_DIR = let v = arg_value("--output-dir")
    v === nothing ? joinpath(@__DIR__, "results", "ont_k_sweep") : v
end

const N_CELLS = length(organisms) * length(technologies) * length(ks) *
                length(coverages) * length(seeds)

# === QUAST report parsing ===
#
# Returns `missing` — never 0.0 — for any metric QUAST could not compute. This
# is the whole point of the harness: the pilot's coercion of "-" to 0.0 is what
# made a censored floor indistinguishable from a measured zero.

function parse_quast_metrics(report_tsv)
    isfile(report_tsv) && return _parse_quast_metrics(report_tsv)
    return empty_metrics()
end

function _parse_quast_metrics(report_tsv)
    values = Dict{String, Union{Missing, Float64}}()
    open(report_tsv, "r") do io
        for line in eachline(io)
            fields = split(line, '\t')
            length(fields) >= 2 || continue
            label = strip(fields[1])
            parsed = tryparse(Float64, strip(fields[2]))
            values[label] = parsed === nothing ? missing : parsed
        end
    end
    get_metric(name) = get(values, name, missing)
    return (
        quast_contigs = get_metric("# contigs"),
        total_length = get_metric("Total length"),
        N50 = get_metric("N50"),
        largest_contig = get_metric("Largest contig"),
        NGA50 = get_metric("NGA50"),
        NA50 = get_metric("NA50"),
        largest_alignment = get_metric("Largest alignment"),
        genome_fraction = get_metric("Genome fraction (%)"),
        duplication_ratio = get_metric("Duplication ratio"),
        misassemblies = get_metric("# misassemblies"),
        unaligned_contigs = get_metric("# unaligned contigs"),
        unaligned_length = get_metric("Unaligned length")
    )
end

function empty_metrics()
    return (
        quast_contigs = missing, total_length = missing, N50 = missing,
        largest_contig = missing, NGA50 = missing, NA50 = missing,
        largest_alignment = missing, genome_fraction = missing,
        duplication_ratio = missing, misassemblies = missing,
        unaligned_contigs = missing, unaligned_length = missing
    )
end

"""
    nga50_status_for(metrics, n_contigs, contigs_ge_min, quast_ran) -> String

Why NGA50 is absent, when it is. Every one of these is a CENSORED FLOOR, not a
measured zero, and they are kept distinct because they have different causes and
different implications:

  * `measured`            — QUAST computed NGA50; the number is real.
  * `no_contigs`          — the assembler emitted nothing at all.
  * `no_contigs_ge_min`   — contigs exist but NONE reaches `--min-contig`
                            (500 bp), so QUAST declines to score the assembly.
                            This is a RESULT, not a tool failure: it is the
                            signature of the low-k regime, where the graph is so
                            tangled that unitigs never grow past a few tens of
                            bases. The `asm_*` columns quantify it.
  * `censored_no_alignment`  — contigs long enough to score exist, but QUAST
                            reported NO genome fraction at all, meaning nothing
                            aligned to the reference. `# misassemblies` does NOT
                            qualify this — QUAST reports 0 misassemblies on
                            unaligned contigs, so a 0 there means non-alignment,
                            not correctness.
  * `censored_gf_below_50`  — contigs ALIGNED (genome fraction is a real,
                            nonzero measurement) but NGA50 is still absent,
                            because NGA50 is UNDEFINED below 50% genome
                            fraction. NGA50 is the aligned-block length at which
                            blocks of that size or larger cover half the
                            REFERENCE; if aligned blocks total under half the
                            genome, no such length exists and QUAST prints "-".
                            This is a definitional boundary, NOT an assembly
                            that failed to align, and conflating the two
                            understates how much of the genome was recovered.
  * `quast_failed`        — QUAST genuinely errored on an assembly it should
                            have been able to score. This is the ONLY value
                            here that indicates an infrastructure problem, which
                            is why `no_contigs_ge_min` was split out of it.

Because NGA50 is a step function of genome fraction at the 50% boundary, it is
an UNSTABLE endpoint wherever cells sit near that line: two replicate seeds of
the same condition can land on opposite sides and report "measured" versus
"absent" for what is nearly the same assembly. Genome fraction is continuous
across that boundary and does not have the pathology, which is why every
summary in this harness reports it beside NGA50.

`contigs_ge_min` is computed in Julia from the contigs themselves, so this
classification never depends on QUAST having succeeded.
"""
function nga50_status_for(metrics, n_contigs::Int, contigs_ge_min::Int,
        quast_ran::Bool)
    n_contigs == 0 && return "no_contigs"
    contigs_ge_min == 0 && return "no_contigs_ge_min"
    quast_ran || return "quast_failed"
    ismissing(metrics.NGA50) || return "measured"
    # NGA50 absent: distinguish "nothing aligned" from "aligned, but under the
    # 50% genome-fraction floor below which NGA50 has no definition".
    ismissing(metrics.genome_fraction) && return "censored_no_alignment"
    return metrics.genome_fraction > 0 ? "censored_gf_below_50" :
           "censored_no_alignment"
end

"""
    contig_stats(contigs, min_contig) -> (; n, n_ge_min, max_length, total_bp)

Assembler-side contig statistics, independent of QUAST.

These are the numbers that survive when QUAST cannot score a cell, and they are
what makes the low-k regime legible: `n` large with `max_length` small is a
tangled graph, whereas `n` large with `max_length` large is a graph that
assembled long but erroneous paths. Those are opposite failure modes and the
sweep has to be able to tell them apart.
"""
function contig_stats(contigs, min_contig::Int)
    lengths = [length(c) for c in contigs]
    return (
        n = length(lengths),
        n_ge_min = count(>=(min_contig), lengths),
        max_length = isempty(lengths) ? 0 : maximum(lengths),
        total_bp = isempty(lengths) ? 0 : sum(lengths)
    )
end

# === Read simulation ===
#
# Deliberately the same call path the pilot used, so the ONT reads here are
# drawn from the same generator whose degeneracy is being explained. Both
# wrappers take an explicit seed, so cells are reproducible.

function simulate_reads(tech, ref_fasta, cov, seed, reads_dir)
    mkpath(reads_dir)
    if tech == "illumina"
        outbase = joinpath(reads_dir, "illumina_$(cov)x")
        res = Mycelia.simulate_illumina_reads(
            fasta = ref_fasta, coverage = cov, outbase = outbase,
            rndSeed = seed, paired = true, quiet = true)
        records = FASTX.FASTQ.Record[]
        for path in (res.forward_reads, res.reverse_reads)
            path === nothing && continue
            reader = Mycelia.open_fastx(path)
            append!(records, collect(reader))
            close(reader)
        end
        return records
    elseif tech == "ont"
        # No --error_model / --identity: inherits the installed Badread's
        # defaults, exactly as the pilot did. See the header note.
        fq = Mycelia.simulate_nanopore_reads(
            fasta = ref_fasta, quantity = "$(cov)x",
            outfile = joinpath(reads_dir, "ont_$(cov)x.fq.gz"),
            seed = seed, quiet = true)
        reader = Mycelia.open_fastx(fq)
        records = collect(reader)
        close(reader)
        return records
    else
        error("unknown technology: $(tech)")
    end
end

# Strip quality: FASTQ -> FASTA gives the plain k-mer graph (the `kmer` arm).
function strip_quality(records)
    [FASTX.FASTA.Record(FASTX.identifier(r), FASTX.sequence(r)) for r in records]
end

# === Per-cell execution ===

function cell_row(org, acc, tech, k, cov, seed; n_reads, asm, metrics,
        nga50_status, outcome, wall_seconds, status)
    return (
        organism = String(org), accession = String(acc), technology = String(tech),
        k = Int(k), coverage = Int(cov), seed = Int(seed),
        decoder_arm = DECODER_ARM,
        n_reads = Int(n_reads), n_contigs = Int(asm.n),
        asm_contigs_ge_min = Int(asm.n_ge_min),
        asm_max_contig = Int(asm.max_length),
        asm_total_bp = Int(asm.total_bp),
        quast_contigs = metrics.quast_contigs,
        total_length = metrics.total_length,
        N50 = metrics.N50,
        largest_contig = metrics.largest_contig,
        NGA50 = metrics.NGA50,
        nga50_status = String(nga50_status),
        NA50 = metrics.NA50,
        largest_alignment = metrics.largest_alignment,
        genome_fraction = metrics.genome_fraction,
        duplication_ratio = metrics.duplication_ratio,
        misassemblies = metrics.misassemblies,
        unaligned_contigs = metrics.unaligned_contigs,
        unaligned_length = metrics.unaligned_length,
        outcome = String(outcome),
        wall_seconds = round(Float64(wall_seconds); digits = 3),
        status = String(status)
    )
end

function run_cell(org, acc, ref, tech, k, cov, seed, cell_dir)
    records = simulate_reads(tech, ref, cov, seed, joinpath(cell_dir, "reads"))
    n_reads = length(records)
    asm_input = strip_quality(records)

    GC.gc()
    # `corrector` is left at its default `:none` — the single-k pipeline. A
    # k-ladder here would confound k with correction, which is the specific
    # confound this sweep exists to avoid.
    timed = @timed Mycelia.Rhizomorph.assemble_genome(asm_input; k = k, verbose = false)
    result = timed.value
    wall_seconds = timed.time

    contigs_path = joinpath(cell_dir, "contigs.fasta")
    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)  # Rhizomorph contigs are String
        end
    end
    asm = contig_stats(result.contigs, MIN_CONTIG)

    # QUAST is invoked only when at least one contig can actually be scored.
    # Calling it with nothing >= --min-contig makes it exit non-zero ("None of
    # the assembly files contains correct contigs"), which is a correct refusal
    # rather than an error — skipping the call keeps that refusal from being
    # logged as a tool failure and keeps the log readable in the low-k regime.
    quast_ran = false
    metrics = empty_metrics()
    if asm.n_ge_min > 0 && filesize(contigs_path) > 0
        try
            quast_dir = joinpath(cell_dir, "quast")
            Mycelia.run_quast([contigs_path]; outdir = quast_dir, reference = ref,
                min_contig = MIN_CONTIG, threads = QUAST_THREADS)
            report = joinpath(quast_dir, "report.tsv")
            if isfile(report)
                metrics = parse_quast_metrics(report)
                quast_ran = true
            end
        catch e
            @warn "QUAST failed" cell=basename(cell_dir) exception=(e, catch_backtrace())
        end
    end

    nga50_status = nga50_status_for(metrics, asm.n, asm.n_ge_min, quast_ran)
    outcome = classify_outcome(
        ismissing(metrics.NGA50) ? 0.0 : metrics.NGA50,
        ismissing(metrics.genome_fraction) ? 0.0 : metrics.genome_fraction,
        nga50_status, genome_size_for(org))

    status = asm.n == 0 ? "empty_assembly" : "ok"
    # Route the fresh row through `reclassify` as well, so a cell computed in
    # this run and the same cell reloaded from its checkpoint are classified by
    # exactly one code path. Without this, the two paths derive `quast_ran`
    # differently (here from whether report.tsv appeared, there from whether a
    # contig count survived into the row) and could disagree on an edge case.
    # `reclassify` is idempotent, so applying it twice is harmless.
    return reclassify(cell_row(org, acc, tech, k, cov, seed; n_reads, asm,
        metrics, nga50_status, outcome, wall_seconds, status))
end

function error_row(org, acc, tech, k, cov, seed)
    cell_row(org, acc, tech, k, cov, seed; n_reads = 0,
        asm = (n = 0, n_ge_min = 0, max_length = 0, total_bp = 0),
        metrics = empty_metrics(), nga50_status = "quast_failed",
        outcome = "degenerate", wall_seconds = 0.0, status = "error")
end

# === Checkpointing ===

cell_id_for(org, tech, k, cov, seed) = "$(org)__$(tech)__k$(k)__$(cov)x__seed$(seed)"

function save_cell_json(path, row)
    open(path, "w") do io
        JSON.print(io, Dict(string(key) => value for (key, value) in pairs(row)), 2)
    end
    return path
end

# Rebuild a canonical, type-coerced NamedTuple from a parsed checkpoint. A
# missing key raises rather than defaulting: every column here is a measurement
# or a grouping key, and `missing`/0 silently substituted for an absent
# measurement is exactly the failure this harness was written to avoid.
function canonical(d::AbstractDict)
    values = map(ROW_KEYS) do key
        name = String(key)
        haskey(d, name) || error(
            "checkpoint is missing required key $(name). Every column in this " *
            "schema is a measurement or a grouping key, so none can be " *
            "defaulted — delete the checkpoint to recompute the cell.")
        v = d[name]
        if key in INT_KEYS
            v isa Integer ? Int(v) : Int(round(Float64(v)))
        elseif key in STR_KEYS
            String(v)
        elseif v === nothing
            missing
        else
            Float64(v)
        end
    end
    return reclassify((; (ROW_KEYS .=> values)...))
end

"""
    reclassify(row) -> NamedTuple

Recompute `nga50_status` and `outcome` from the row's own recorded measurements
rather than trusting whatever the checkpoint stored.

Both are DERIVED columns — pure functions of (NGA50, genome_fraction,
n_contigs, asm_contigs_ge_min), all of which are stored. Deriving them on read
means a classification bug can be fixed by editing this file and re-running the
aggregation, instead of by recomputing 96 assemblies. That mattered in practice:
the first version of `nga50_status_for` labelled every absent NGA50
`censored_unaligned`, which was wrong for cells that aligned but fell under the
50% genome-fraction floor where NGA50 is undefined — and those cells had
genome fractions as high as 31.7%, so the label was asserting "nothing aligned"
about assemblies that had recovered a third of the genome.

The raw measurements are never touched, so this is lossless relabelling.
"""
function reclassify(row)
    metrics = (NGA50 = row.NGA50, genome_fraction = row.genome_fraction)
    quast_ran = !ismissing(row.quast_contigs)
    status = nga50_status_for(metrics, row.n_contigs, row.asm_contigs_ge_min,
        quast_ran)
    # Genome size comes from the ROW's own organism, not from a constant, so a
    # multi-organism table is reclassified against the right tier scale for each
    # row rather than against whichever organism happened to be first.
    outcome = classify_outcome(
        ismissing(row.NGA50) ? 0.0 : row.NGA50,
        ismissing(row.genome_fraction) ? 0.0 : row.genome_fraction,
        status, genome_size_for(row.organism))
    return merge(row, (; nga50_status = status, outcome = outcome))
end

# === Aggregation ===

function write_aggregate(root, rows)
    df = DataFrames.DataFrame(rows)
    sort!(df, [:organism, :technology, :k, :coverage, :seed])
    CSV.write(joinpath(root, "ont_k_sweep_results.tsv"), df; delim = '\t', missingstring = "NA")
    return df
end

"""
Per-(technology, k, coverage) summary across seeds.

STRATIFIED BY TECHNOLOGY BY CONSTRUCTION — the grouping key leads with
`technology`, so no row of this table ever mixes ONT with Illumina. Those are
different error processes and an aggregate over both would describe neither.

`n_degenerate` is reported beside the medians because a median NGA50 over a
group containing censored cells is not a measurement of that group, and
`nga50_status_mix` names WHICH censoring applies so a reader never has to guess
whether a blank median means "nothing aligned" or "nothing was long enough to
score". `median_max_contig` is carried for the same reason: it is the one
contiguity number that is defined in every regime, including the ones QUAST
declines to score.
"""
function write_summary(root, df)
    # Cells that threw are EXCLUDED. An error row carries outcome="degenerate"
    # and genome_fraction=missing, which is shaped exactly like a genuinely
    # degenerate cell — so leaving it in would let a transient infrastructure
    # failure count as evidence for degeneracy, the very thing the sweep is
    # measuring. `n_seeds` then reports how many cells actually contributed,
    # so a stratum thinned by a failure is visible rather than silently averaged.
    df = df[df.status .== "ok", :]
    summary_rows = NamedTuple[]
    for g in DataFrames.groupby(df, [:organism, :technology, :k, :coverage])
        measured = g[g.nga50_status .== "measured", :]
        gf = collect(skipmissing(g.genome_fraction))
        push!(summary_rows,
            (
                organism = g.organism[1], technology = g.technology[1],
                k = g.k[1], coverage = g.coverage[1],
                n_seeds = DataFrames.nrow(g),
                n_measured_nga50 = DataFrames.nrow(measured),
                n_degenerate = count(==("degenerate"), g.outcome),
                median_nga50 = DataFrames.nrow(measured) == 0 ? missing :
                               Statistics.median(collect(skipmissing(measured.NGA50))),
                median_genome_fraction = isempty(gf) ? missing : Statistics.median(gf),
                median_contigs = Statistics.median(Float64.(g.n_contigs)),
                median_contigs_ge_min = Statistics.median(Float64.(g.asm_contigs_ge_min)),
                median_max_contig = Statistics.median(Float64.(g.asm_max_contig)),
                nga50_status_mix = join(sort(unique(g.nga50_status)), "|"),
                outcome = join(sort(unique(g.outcome)), "|"),
                median_wall_seconds = Statistics.median(Float64.(g.wall_seconds))
            ))
    end
    summary_df = DataFrames.DataFrame(summary_rows)
    sort!(summary_df, [:organism, :technology, :k, :coverage])
    CSV.write(joinpath(root, "ont_k_sweep_summary.tsv"), summary_df;
        delim = '\t', missingstring = "NA")
    return summary_df
end

# === Main ===
#
# Guarded so a test can `include()` this file for the pure helpers
# (`classify_outcome`, `nga50_status_for`, `parse_quast_metrics`, `canonical`)
# without downloading a genome or running an assembly. Matches the convention
# used by the other harnesses in this directory.
if abspath(PROGRAM_FILE) == @__FILE__
    println("=== ONT k-selection sweep (td-4e19d.28) ===")
    println("Start: $(Dates.now())")
    println("Smoke mode: $(SMOKE)")
    println("Organisms: $(join((o[1] for o in organisms), ", "))")
    println("Technologies: $(join(technologies, ", "))")
    println("k values: $(join(ks, ", "))")
    println("Coverages: $(join(coverages, ", "))x")
    println("Seeds: $(join(seeds, ", "))")
    println("Decoder arm: $(DECODER_ARM)")
    println("Cells to run: $(N_CELLS)")
    println("Output dir: $(OUTPUT_DIR)")

    mkpath(OUTPUT_DIR)
    refs_dir = joinpath(OUTPUT_DIR, "refs")
    cells_dir = joinpath(OUTPUT_DIR, "cells")
    mkpath(refs_dir)
    mkpath(cells_dir)

    println("\n--- Phase 1: reference genomes ---")
    ref_paths = Dict{String, String}()
    for (org, acc, expected) in organisms
        path = Mycelia.download_genome_by_accession(
            accession = acc, outdir = refs_dir, compressed = false)
        if !isfile(path) || filesize(path) == 0
            error("reference download failed for $(org) ($(acc)): $(path)")
        end
        actual = Mycelia.total_fasta_size(path)
        # A silently-wrong reference would rescale every NGA50 tier for this
        # organism, so check the size against the table rather than trusting the
        # accession to have resolved to what we meant.
        if abs(actual - expected) > 0.2 * expected
            error("reference size mismatch for $(org) ($(acc)): expected " *
                  "~$(expected) bp, got $(actual) bp")
        end
        ref_paths[org] = path
        println("  $(org) ($(acc)): $(actual) bp -> $(path)")
    end

    println("\n--- Phase 2: sweep ($(N_CELLS) cells) ---")
    # Loop order runs the cheap strata first (Illumina before ONT, low coverage
    # before high) so the control arm is banked early and a run that is interrupted
    # still yields complete low-coverage strata rather than a ragged fringe.
    rows = NamedTuple[]
    cell_index = 0
    for (org, acc, _expected) in organisms,
        tech in technologies,
        k in ks,
        cov in coverages,
        seed in seeds

        global cell_index += 1
        cell_id = cell_id_for(org, tech, k, cov, seed)
        cell_dir = joinpath(cells_dir, cell_id)
        ckpt = joinpath(cell_dir, "cell_result.json")

        # A checkpoint is reused ONLY if it recorded a real attempt. An `error`
        # row is a cell that threw — under 32 concurrent shards the observed
        # cause was a transient collision in read simulation, with n_reads = 0
        # and wall_seconds = 0 — and caching that would freeze a transient
        # failure into a permanent hole in the grid that no amount of
        # re-running could fill. Errors are retried; successes are never
        # recomputed.
        if isfile(ckpt)
            cached = canonical(JSON.parsefile(ckpt))
            if cached.status == "error"
                println("  [$(cell_index)/$(N_CELLS)] $(cell_id) — cached row is an " *
                        "error; retrying")
            else
                println("  [$(cell_index)/$(N_CELLS)] $(cell_id) — cached, skipping")
                push!(rows, cached)
                continue
            end
        end

        mkpath(cell_dir)
        print("  [$(cell_index)/$(N_CELLS)] $(cell_id) ... ")
        flush(stdout)
        row = try
            run_cell(org, acc, ref_paths[org], tech, k, cov, seed, cell_dir)
        catch e
            @warn "cell failed" cell_id exception = (e, catch_backtrace())
            error_row(org, acc, tech, k, cov, seed)
        end
        save_cell_json(ckpt, row)
        push!(rows, row)
        write_aggregate(OUTPUT_DIR, rows)  # rewrite each cell; cheap at this scale
        println("$(row.outcome) [$(row.status)]: $(row.n_contigs) contigs " *
                "($(row.asm_contigs_ge_min) >= $(MIN_CONTIG)bp, max $(row.asm_max_contig)bp), " *
                "NGA50=$(ismissing(row.NGA50) ? row.nga50_status : row.NGA50), " *
                "GF=$(ismissing(row.genome_fraction) ? "NA" : row.genome_fraction)%, " *
                "$(round(row.wall_seconds; digits = 1))s")
        flush(stdout)
    end

    println("\n--- Phase 3: aggregate ---")
    results_df = write_aggregate(OUTPUT_DIR, rows)
    summary_df = write_summary(OUTPUT_DIR, results_df)

    println("\nPer-(technology, k, coverage) summary:")
    show(stdout, summary_df; allrows = true, allcols = true)
    println()

    n_ok = count(r -> r.status == "ok", rows)
    println("\nDone: $(length(rows)) cells, $(n_ok) ok. Results in $(OUTPUT_DIR)")

    # Say plainly whether the grid is complete. An errored cell carries
    # `degenerate` in its outcome column and is otherwise shaped exactly like a
    # real degenerate result, so without this it would be counted as evidence
    # for the very conclusion the sweep is testing. Name the cells so a reader
    # of the log never has to reconstruct which ones are missing.
    failed = filter(r -> r.status != "ok", rows)
    if isempty(failed)
        println("All $(length(rows)) cells completed. Grid is COMPLETE.")
    else
        println("\n*** WARNING: $(length(failed)) cell(s) did NOT complete. Their " *
                "`outcome` is a placeholder, NOT a measurement, and must be " *
                "excluded from any summary. Re-run to retry them. ***")
        for r in failed
            println("  - $(r.technology) k=$(r.k) $(r.coverage)x seed$(r.seed) " *
                    "[$(r.status)]")
        end
    end
    println("End: $(Dates.now())")
end  # PROGRAM_FILE guard
