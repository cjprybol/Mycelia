# Track-A λ Pilot Benchmark — greedy-only (current-algorithm) real-data pilot
#
# Purpose
# -------
# This is the ungated FIRST step of the Track-A rhizomorph benchmark program.
# It runs the CURRENT greedy assembly algorithm on real Lambda-phage Illumina
# reads across a coverage x seed grid, scores each assembly against the Lambda
# reference with QUAST, and reports the NGA50 coefficient of variation (CV)
# across cells. That observed NGA50 CV is the quantity the rhizomorph
# pre-registration power analysis needs; producing it lifts the pre-reg lock.
#
# Design
# ------
# - Organism: Lambda phage (accession NC_001416, ~48.5 kb), Illumina reads.
# - Grid: coverage {10, 30, 50, 100}x  x  seed {42, 123, 456} = 12 cells.
# - Metrics per cell come from QUAST's ALIGNED-truth report (report.tsv):
#   NGA50, # misassemblies, Duplication ratio, Largest alignment,
#   Genome fraction (%). NGA50 here is the HONEST aligned N50 — NOT the
#   size-ratio "genome fraction" proxy used in t4_ksweep.jl /
#   real_genome_benchmark.jl (that Phase-4 proxy divides assembled length by
#   reference size and is a known over-count; this driver does not use it).
#
# Read-simulation path: DIRECT-PER-COVERAGE (see note below).
#   `Mycelia.simulate_illumina_reads` accepts `coverage` (ART --fcov) AND
#   `rndSeed` (ART --rndSeed) directly, so each (coverage, seed) cell is
#   simulated independently with full seed control. We deliberately do NOT
#   simulate-once-then-subsample, because `Mycelia.subsample_reads_seqkit`
#   exposes no rand-seed argument (seqkit's --rand-seed is hard-defaulted),
#   so a subsample path could not deliver seed-controlled draws per cell.
#   Direct simulation is the honest seed-controlled path the real API supports.
#
# External tools
# --------------
# The whole pilot requires Conda-backed external tools (art_illumina for
# simulation, QUAST for scoring); these are set up on Lovelace, not the laptop.
# QUAST specifically is gated behind MYCELIA_RUN_QUAST (default "true" — QUAST
# is the point of the pilot). Set MYCELIA_RUN_QUAST=false to run assembly-only
# (no aligned-truth metrics). LD_LIBRARY_PATH / thread caps are the caller's
# job (the run command sets them).
#
# Env knobs (all optional):
#   MYCELIA_RUN_QUAST        default "true"  — run QUAST scoring
#   MYCELIA_PILOT_K          default "31"    — k-mer size for greedy assembly
#   MYCELIA_QUAST_MIN_CONTIG default "100"   — QUAST --min-contig (fixed across
#                                              cells so it does not add CV noise)
#
# Usage (on Lovelace, see PR for the exact launch command):
#   julia --project=. benchmarking/track_a_lambda_pilot.jl

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import DataFrames
import CSV
import Dates
import Statistics

# === Configuration ===

const ORGANISM::String = "Lambda"
const ACCESSION::String = "NC_001416"
const TECHNOLOGY::String = "illumina"
const COVERAGES::Vector{Int} = [10, 30, 50, 100]
const SEEDS::Vector{Int} = [42, 123, 456]

const RUN_QUAST::Bool = get(ENV, "MYCELIA_RUN_QUAST", "true") in ["1", "true", "yes"]
const K::Int = parse(Int, get(ENV, "MYCELIA_PILOT_K", "31"))
const QUAST_MIN_CONTIG::Int = parse(Int, get(ENV, "MYCELIA_QUAST_MIN_CONTIG", "100"))

println("=== Track-A Lambda Pilot ===")
println("Start time (UTC): $(Dates.now(Dates.UTC))")
println("Organism: $ORGANISM ($ACCESSION), technology: $TECHNOLOGY")
println("Coverage tiers: $COVERAGES  x  seeds: $SEEDS  = $(length(COVERAGES) * length(SEEDS)) cells")
println("Greedy assembly k: $K")
println("QUAST scoring: $(RUN_QUAST ? "ENABLED (min_contig=$QUAST_MIN_CONTIG)" : "DISABLED (MYCELIA_RUN_QUAST=false)")")

# Working + results directories
benchmark_dir = mktempdir(prefix = "track_a_lambda_pilot_")
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
println("Working directory: $benchmark_dir")

# === QUAST report.tsv parsing ===
# QUAST writes a tab-separated report.tsv whose first column is the metric name
# and second column is the value for the (single) assembly we pass. Exact metric
# names below are QUAST's canonical labels for reference-based (aligned) metrics.

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

# === Phase 1: Download Lambda reference (once) ===

println("\n--- Phase 1: Download Lambda reference ---")
ref_dir = joinpath(benchmark_dir, "reference")
mkpath(ref_dir)
reference_path = Mycelia.download_genome_by_accession(
    accession = ACCESSION,
    outdir = ref_dir,
    compressed = false
)
if !isfile(reference_path) || filesize(reference_path) == 0
    error("Failed to download Lambda reference ($ACCESSION) to $ref_dir")
end
println("Reference: $reference_path ($(filesize(reference_path)) bytes)")

# === Phase 2: coverage x seed sweep ===

results = DataFrames.DataFrame(
    organism = String[],
    technology = String[],
    coverage = Int[],
    seed = Int[],
    contigs = Union{Missing, Int}[],
    NGA50 = Union{Missing, Float64}[],
    num_misassemblies = Union{Missing, Float64}[],
    dup_ratio = Union{Missing, Float64}[],
    largest_alignment = Union{Missing, Float64}[],
    genome_fraction_quast = Union{Missing, Float64}[],
    runtime_s = Union{Missing, Float64}[],
    peak_rss_kb = Union{Missing, Int}[]
)

"""
    load_reads(sim) -> Vector{FASTX.FASTQ.Record}

Load the forward and (paired) reverse simulated reads into a single flat record
vector, matching how the greedy assembler consumes reads (no pairing metadata).
"""
function load_reads(sim)::Vector{FASTX.FASTQ.Record}
    records = FASTX.FASTQ.Record[]
    read_paths = String[sim.forward_reads]
    if sim.reverse_reads !== nothing
        push!(read_paths, sim.reverse_reads)
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
    run_cell(coverage, seed, cell_dir) -> NamedTuple

Simulate reads, run the current greedy assembly, and (optionally) score with
QUAST for one (coverage, seed) cell. Returns the per-cell metric row.
"""
function run_cell(coverage::Int, seed::Int, cell_dir::String)
    mkpath(cell_dir)

    # -- Simulate Illumina reads directly at this coverage + seed --
    outbase = joinpath(cell_dir, "$(ORGANISM)_cov$(coverage)x_seed$(seed)")
    sim = Mycelia.simulate_illumina_reads(
        fasta = reference_path,
        coverage = coverage,
        rndSeed = seed,
        read_length = 150,
        paired = true,
        quiet = true,
        outbase = outbase
    )
    records = load_reads(sim)
    if isempty(records)
        @warn "No reads simulated for coverage=$coverage seed=$seed"
        return nothing
    end

    # -- Current greedy assembly (quality-aware qualmer path auto-detected for FASTQ) --
    t0 = time()
    result = Mycelia.Rhizomorph.assemble_genome(records; k = K, verbose = false)
    runtime = time() - t0
    peak_rss_kb = Int(Sys.maxrss() ÷ 1024)  # process-cumulative peak RSS (not per-cell isolated)

    n_contigs = length(result.contigs)
    contigs_path = joinpath(cell_dir, "$(ORGANISM)_cov$(coverage)x_seed$(seed)_contigs.fasta")
    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)  # contigs are String type
        end
    end

    # -- QUAST aligned-truth scoring --
    nga50 = missing
    num_misassemblies = missing
    dup_ratio = missing
    largest_alignment = missing
    genome_fraction = missing
    if RUN_QUAST && n_contigs > 0 && filesize(contigs_path) > 0
        quast_dir = joinpath(cell_dir, "quast")
        Mycelia.run_quast(
            contigs_path;
            outdir = quast_dir,
            reference = reference_path,
            min_contig = QUAST_MIN_CONTIG
        )
        report_tsv = joinpath(quast_dir, "report.tsv")
        nga50 = parse_quast_metric(report_tsv, "NGA50")
        num_misassemblies = parse_quast_metric(report_tsv, "# misassemblies")
        dup_ratio = parse_quast_metric(report_tsv, "Duplication ratio")
        largest_alignment = parse_quast_metric(report_tsv, "Largest alignment")
        genome_fraction = parse_quast_metric(report_tsv, "Genome fraction (%)")
    end

    return (
        contigs = n_contigs,
        NGA50 = nga50,
        num_misassemblies = num_misassemblies,
        dup_ratio = dup_ratio,
        largest_alignment = largest_alignment,
        genome_fraction_quast = genome_fraction,
        runtime_s = round(runtime; digits = 3),
        peak_rss_kb = peak_rss_kb
    )
end

println("\n--- Phase 2: coverage x seed sweep ---")
for seed in SEEDS
    for coverage in COVERAGES
        println("  cell: coverage=$(coverage)x seed=$(seed) ...")
        cell_dir = joinpath(benchmark_dir, "cov$(coverage)x_seed$(seed)")
        try
            cell = run_cell(coverage, seed, cell_dir)
            if cell === nothing
                continue
            end
            push!(results,
                (
                    organism = ORGANISM,
                    technology = TECHNOLOGY,
                    coverage = coverage,
                    seed = seed,
                    contigs = cell.contigs,
                    NGA50 = cell.NGA50,
                    num_misassemblies = cell.num_misassemblies,
                    dup_ratio = cell.dup_ratio,
                    largest_alignment = cell.largest_alignment,
                    genome_fraction_quast = cell.genome_fraction_quast,
                    runtime_s = cell.runtime_s,
                    peak_rss_kb = cell.peak_rss_kb
                ))
            println("    -> $(cell.contigs) contigs, NGA50=$(cell.NGA50), " *
                    "misasm=$(cell.num_misassemblies), GF=$(cell.genome_fraction_quast)%, " *
                    "$(cell.runtime_s)s")
        catch e
            @warn "Cell failed (coverage=$coverage seed=$seed) — skipping" exception = (e, catch_backtrace())
            continue
        end
    end
end

println("\nCompleted $(DataFrames.nrow(results))/$(length(COVERAGES) * length(SEEDS)) cells")

# === Phase 3: write results CSV ===

# NOTE: format with no 'Z' directive — 'Z' is a TimeZones timezone token that
# errors on a plain Dates.DateTime ("type DateTime has no field zone"). The
# timestamp is UTC (Dates.now(Dates.UTC)); matches the repo's yyyymmdd_HHMMSS convention.
timestamp = Dates.format(Dates.now(Dates.UTC), "yyyymmdd_HHMMSS")
csv_path = joinpath(results_dir, "track_a_lambda_pilot_$(timestamp).csv")
CSV.write(csv_path, results)
println("\nResults written to: $csv_path")
println(results)

# === Phase 4: NGA50 coefficient of variation (the pre-reg deliverable) ===

"""
    nga50_cv(values) -> Union{Missing, Float64}

Coefficient of variation (std/mean) over the non-missing NGA50 values. Returns
`missing` if fewer than two numeric values are available.
"""
function nga50_cv(values)::Union{Missing, Float64}
    numeric = Float64[v for v in values if !ismissing(v)]
    length(numeric) < 2 && return missing
    m = Statistics.mean(numeric)
    m == 0 && return missing
    return Statistics.std(numeric) / m
end

println("\n--- Phase 4: NGA50 coefficient of variation ---")
overall_cv = nga50_cv(results.NGA50)
n_numeric = count(!ismissing, results.NGA50)
println("Overall NGA50 CV across $(n_numeric) scored cells: " *
        (ismissing(overall_cv) ? "n/a (insufficient numeric NGA50 values)" :
         string(round(overall_cv; digits = 4))))

println("\nPer-coverage NGA50 CV:")
for coverage in COVERAGES
    sub = results[results.coverage .== coverage, :]
    cv = nga50_cv(sub.NGA50)
    vals = [ismissing(v) ? "NA" : string(v) for v in sub.NGA50]
    println("  $(coverage)x: CV=" *
            (ismissing(cv) ? "n/a" : string(round(cv; digits = 4))) *
            "   NGA50=[$(join(vals, ", "))]")
end

println("\n=== Track-A Lambda Pilot complete (UTC): $(Dates.now(Dates.UTC)) ===")
