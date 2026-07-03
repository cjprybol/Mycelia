# Track-A Headroom Pilot Benchmark — greedy-only "distance-from-perfect" sweep
#
# Purpose
# -------
# The base Track-A Lambda pilot (`benchmarking/track_a_lambda_pilot.jl`) showed
# a CEILING: the current greedy assembler already recovers ~99.8% of Lambda
# (NGA50 ~48.4k / 48.5k, 0 misassemblies) at >=30x coverage. At the ceiling
# there is NO room for a future Viterbi-DP strategy to help — greedy is already
# effectively perfect. The only fragmented base-pilot cell was 10x
# (NGA50 ~1467). So the DISCRIMINATING regime — where a DP strategy could
# plausibly beat greedy — is (i) LOW coverage and (ii) HARDER / LARGER genomes.
#
# This extended pilot MEASURES that headroom with the greedy baseline alone.
# It does NOT need a DP strategy: the DP/greedy strategy pair does not exist
# yet. The point is to quantify the OPPORTUNITY = how far the current greedy
# assembly is from perfect, per (organism, coverage) cell, so the Track-A
# program can locate the regime worth building a DP strategy for.
#
# What "headroom" means (distance-from-perfect; higher = more opportunity)
# ------------------------------------------------------------------------
# A perfect assembly of a single-replicon genome is ONE contig spanning the
# whole reference with 100% genome fraction and 0 misassemblies. We report
# three complementary distance-from-perfect measures per cell:
#
#   headroom_nga50 = 1 - NGA50 / ref_size
#       0.0  => NGA50 equals the genome size (one perfect contig, no headroom).
#       ->1  => heavily fragmented aligned blocks (large headroom for DP).
#       ref_size is the reference length in bp (summed from the reference FASTA).
#
#   headroom_gf = 1 - genome_fraction_quast / 100
#       0.0  => QUAST genome fraction is 100% (whole genome recovered).
#       ->1  => most of the genome is unrecovered (large headroom).
#
#   num_misassemblies (absolute count, reported as-is)
#       0    => no structural errors. >0 => structural headroom: cells where
#       greedy makes join errors a DP strategy might avoid.
#
# A cell is "materially below perfect" (a DP-opportunity cell) when
# headroom_nga50 exceeds a small threshold (default 0.05, i.e. NGA50 < 95% of
# ref_size) OR num_misassemblies > 0. Those cells are highlighted in the summary.
#
# Design
# ------
# Two organism blocks share ONE download -> simulate -> greedy-assemble ->
# QUAST -> parse pipeline (reused verbatim from the base Lambda pilot):
#
#   1. LOW-COVERAGE LAMBDA CURVE (cheap): Lambda (NC_001416, ~48.5 kb),
#      Illumina, coverage {5,10,15,20,30}x x seed {42,123,456} = 15 cells.
#      Maps the fragmentation -> saturation transition precisely.
#
#   2. BACTERIAL FEASIBILITY PROBE (BOUNDED — can be compute-heavy, guarded):
#      E. coli K-12 MG1655 (U00096.3, ~4.64 Mb), Illumina, coverage {10,30}x,
#      seed {42} = 2 cells. This is a bounded probe, not a full grid. The
#      greedy assembler emitted ~17k contigs for the 48 kb Lambda genome; on a
#      4.6 Mb genome that could produce very large contig counts / high memory.
#      Each bacterial cell is wrapped in try/catch (@warn + continue) and prints
#      elapsed seconds + Sys.maxrss() so a blow-up is visible rather than fatal.
#      Bacterial cells are gated behind RUN_BACTERIAL (see constants) so the run
#      can be scoped down or disabled. 100x bacterial is deliberately NOT run.
#
# QUAST NGA50 here is the HONEST aligned N50 from QUAST's aligned-truth
# report.tsv — NOT the size-ratio "genome fraction" proxy used elsewhere.
#
# External tools
# --------------
# The whole pilot requires Conda-backed external tools (art_illumina for
# simulation, QUAST for scoring); these are set up on Lovelace, not the laptop.
# QUAST is gated behind MYCELIA_RUN_QUAST (default "true" — QUAST is the point
# of the pilot). Set MYCELIA_RUN_QUAST=false to run assembly-only (no
# aligned-truth metrics; headroom_nga50 / headroom_gf will be missing).
# LD_LIBRARY_PATH / thread caps are the caller's job (the run command sets them).
#
# Env knobs (all optional):
#   MYCELIA_RUN_QUAST         default "true"  — run QUAST scoring
#   MYCELIA_PILOT_K           default "31"    — k-mer size for greedy assembly
#   MYCELIA_QUAST_MIN_CONTIG  default "100"   — QUAST --min-contig (fixed across
#                                               cells so it does not add noise)
#   MYCELIA_RUN_BACTERIAL     default "true"  — include the E. coli probe;
#                                               set "false" for a Lambda-only run
#
# Usage (on Lovelace, see PR for the exact launch command):
#   julia --project=. benchmarking/headroom_pilot.jl

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

# Low-coverage Lambda curve (cheap; full seed replication).
const LAMBDA_ORGANISM::String = "Lambda"
const LAMBDA_ACCESSION::String = "NC_001416"
const LAMBDA_COVERAGES::Vector{Int} = [5, 10, 15, 20, 30]
const LAMBDA_SEEDS::Vector{Int} = [42, 123, 456]

# Bacterial feasibility probe (BOUNDED — heavy; scope with these constants).
# Set RUN_BACTERIAL=false (or MYCELIA_RUN_BACTERIAL=false) for a Lambda-only run.
# U00096.3 is the canonical E. coli K-12 MG1655 GenBank accession (~4.64 Mb);
# NC_000913.3 is the equivalent RefSeq accession if U00096.3 ever fails to fetch.
const RUN_BACTERIAL::Bool =
    get(ENV, "MYCELIA_RUN_BACTERIAL", "true") in ["1", "true", "yes"]
const BACTERIAL_ORGANISM::String = "Ecoli_K12_MG1655"
const BACTERIAL_ACCESSION::String = "U00096.3"
const BACTERIAL_COVERAGES::Vector{Int} = [10, 30]
const BACTERIAL_SEEDS::Vector{Int} = [42]

const TECHNOLOGY::String = "illumina"

const RUN_QUAST::Bool = get(ENV, "MYCELIA_RUN_QUAST", "true") in ["1", "true", "yes"]
const K::Int = parse(Int, get(ENV, "MYCELIA_PILOT_K", "31"))
const QUAST_MIN_CONTIG::Int = parse(Int, get(ENV, "MYCELIA_QUAST_MIN_CONTIG", "100"))

# A cell is flagged as a DP-opportunity cell when its NGA50 falls below this
# fraction of the reference size (headroom_nga50 > threshold) or it has any
# misassemblies. 0.05 => NGA50 must reach >=95% of ref_size to count as "no
# headroom" (matches the base pilot's ~99.8% Lambda ceiling).
const HEADROOM_FLAG_THRESHOLD::Float64 = 0.05

println("=== Track-A Headroom Pilot ===")
println("Start time (UTC): $(Dates.now(Dates.UTC))")
println("Lambda curve: cov=$LAMBDA_COVERAGES x seeds=$LAMBDA_SEEDS " *
        "= $(length(LAMBDA_COVERAGES) * length(LAMBDA_SEEDS)) cells")
if RUN_BACTERIAL
    println("Bacterial probe ($BACTERIAL_ORGANISM, $BACTERIAL_ACCESSION): " *
            "cov=$BACTERIAL_COVERAGES x seeds=$BACTERIAL_SEEDS " *
            "= $(length(BACTERIAL_COVERAGES) * length(BACTERIAL_SEEDS)) cells (BOUNDED)")
else
    println("Bacterial probe: DISABLED (RUN_BACTERIAL=false)")
end
println("Greedy assembly k: $K")
println("QUAST scoring: $(RUN_QUAST ? "ENABLED (min_contig=$QUAST_MIN_CONTIG)" : "DISABLED (MYCELIA_RUN_QUAST=false)")")

# Working + results directories
benchmark_dir = mktempdir(prefix = "track_a_headroom_pilot_")
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
println("Working directory: $benchmark_dir")

# === QUAST report.tsv parsing ===
# QUAST writes a tab-separated report.tsv whose first column is the metric name
# and second column is the value for the (single) assembly we pass. Reused
# verbatim from the base Lambda pilot.

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
    reference_length_bp(reference_path) -> Int

Total reference length in bp, summed over all records in the reference FASTA.
Used as the denominator for headroom_nga50 (distance-from-perfect). Independent
of QUAST so it is available even when MYCELIA_RUN_QUAST=false.
"""
function reference_length_bp(reference_path::String)::Int
    total = 0
    reader = Mycelia.open_fastx(reference_path)
    for record in reader
        total += length(FASTX.sequence(String, record))
    end
    close(reader)
    return total
end

# === Results table ===

results = DataFrames.DataFrame(
    organism = String[],
    technology = String[],
    ref_size = Int[],
    coverage = Int[],
    seed = Int[],
    contigs = Union{Missing, Int}[],
    NGA50 = Union{Missing, Float64}[],
    num_misassemblies = Union{Missing, Float64}[],
    dup_ratio = Union{Missing, Float64}[],
    largest_alignment = Union{Missing, Float64}[],
    genome_fraction_quast = Union{Missing, Float64}[],
    runtime_s = Union{Missing, Float64}[],
    peak_rss_kb = Union{Missing, Int}[],
    headroom_nga50 = Union{Missing, Float64}[],
    headroom_gf = Union{Missing, Float64}[]
)

"""
    load_reads(sim) -> Vector{FASTX.FASTQ.Record}

Load the forward and (paired) reverse simulated reads into a single flat record
vector, matching how the greedy assembler consumes reads (no pairing metadata).
Reused verbatim from the base Lambda pilot.
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
    headroom_from(nga50, genome_fraction, ref_size) -> (headroom_nga50, headroom_gf)

Distance-from-perfect measures for one cell. Each is `missing` when its input
metric is missing. See the file header for full definitions.
"""
function headroom_from(nga50::Union{Missing, Float64},
        genome_fraction::Union{Missing, Float64},
        ref_size::Int)
    hr_nga50 = (ismissing(nga50) || ref_size == 0) ? missing : 1.0 - nga50 / ref_size
    hr_gf = ismissing(genome_fraction) ? missing : 1.0 - genome_fraction / 100.0
    return (hr_nga50, hr_gf)
end

"""
    run_cell(organism, reference_path, ref_size, coverage, seed, cell_dir) -> NamedTuple

Simulate reads, run the current greedy assembly, and (optionally) score with
QUAST for one (organism, coverage, seed) cell. Returns the per-cell metric row,
including derived headroom. Structure reused from the base Lambda pilot.
"""
function run_cell(organism::String, reference_path::String, ref_size::Int,
        coverage::Int, seed::Int, cell_dir::String)
    mkpath(cell_dir)

    # -- Simulate Illumina reads directly at this coverage + seed --
    outbase = joinpath(cell_dir, "$(organism)_cov$(coverage)x_seed$(seed)")
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
        @warn "No reads simulated for organism=$organism coverage=$coverage seed=$seed"
        return nothing
    end

    # -- Current greedy assembly (quality-aware qualmer path auto-detected for FASTQ) --
    t0 = time()
    result = Mycelia.Rhizomorph.assemble_genome(records; k = K, verbose = false)
    runtime = time() - t0
    peak_rss_kb = Int(Sys.maxrss() ÷ 1024)  # process-cumulative peak RSS (not per-cell isolated)

    n_contigs = length(result.contigs)
    contigs_path = joinpath(cell_dir, "$(organism)_cov$(coverage)x_seed$(seed)_contigs.fasta")
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

    hr_nga50, hr_gf = headroom_from(nga50, genome_fraction, ref_size)

    return (
        contigs = n_contigs,
        NGA50 = nga50,
        num_misassemblies = num_misassemblies,
        dup_ratio = dup_ratio,
        largest_alignment = largest_alignment,
        genome_fraction_quast = genome_fraction,
        runtime_s = round(runtime; digits = 3),
        peak_rss_kb = peak_rss_kb,
        headroom_nga50 = hr_nga50,
        headroom_gf = hr_gf
    )
end

"""
    run_block(organism, accession, coverages, seeds; bacterial=false)

Download the reference for `organism` once, then sweep `coverages` x `seeds`,
pushing one row per completed cell into the shared `results` DataFrame. When
`bacterial=true`, each cell prints elapsed seconds + peak RSS and is wrapped so
a heavy/failing bacterial cell @warns and continues rather than aborting the run.
"""
function run_block(organism::String, accession::String,
        coverages::Vector{Int}, seeds::Vector{Int}; bacterial::Bool = false)
    println("\n--- Block: $organism ($accession) ---")

    # Download reference once for this organism.
    ref_dir = joinpath(benchmark_dir, "reference_$(organism)")
    mkpath(ref_dir)
    reference_path = Mycelia.download_genome_by_accession(
        accession = accession,
        outdir = ref_dir,
        compressed = false
    )
    if !isfile(reference_path) || filesize(reference_path) == 0
        @warn "Failed to download reference ($accession) for $organism — skipping block"
        return
    end
    ref_size = reference_length_bp(reference_path)
    println("Reference: $reference_path ($(filesize(reference_path)) bytes, ref_size=$ref_size bp)")

    for seed in seeds
        for coverage in coverages
            println("  cell: $organism coverage=$(coverage)x seed=$(seed) ...")
            cell_dir = joinpath(benchmark_dir, "$(organism)_cov$(coverage)x_seed$(seed)")
            wall0 = time()
            try
                cell = run_cell(organism, reference_path, ref_size, coverage, seed, cell_dir)
                if cell === nothing
                    continue
                end
                push!(results,
                    (
                        organism = organism,
                        technology = TECHNOLOGY,
                        ref_size = ref_size,
                        coverage = coverage,
                        seed = seed,
                        contigs = cell.contigs,
                        NGA50 = cell.NGA50,
                        num_misassemblies = cell.num_misassemblies,
                        dup_ratio = cell.dup_ratio,
                        largest_alignment = cell.largest_alignment,
                        genome_fraction_quast = cell.genome_fraction_quast,
                        runtime_s = cell.runtime_s,
                        peak_rss_kb = cell.peak_rss_kb,
                        headroom_nga50 = cell.headroom_nga50,
                        headroom_gf = cell.headroom_gf
                    ))
                hr = ismissing(cell.headroom_nga50) ? "NA" : string(round(cell.headroom_nga50; digits = 4))
                println("    -> $(cell.contigs) contigs, NGA50=$(cell.NGA50), " *
                        "misasm=$(cell.num_misassemblies), GF=$(cell.genome_fraction_quast)%, " *
                        "headroom_nga50=$(hr), $(cell.runtime_s)s")
            catch e
                @warn "Cell failed (organism=$organism coverage=$coverage seed=$seed) — skipping" exception = (e, catch_backtrace())
                continue
            finally
                if bacterial
                    # Per-cell compute-risk telemetry: watch for blow-ups on the
                    # large bacterial genome (huge contig counts / high memory).
                    elapsed = round(time() - wall0; digits = 1)
                    println("    [bacterial cell telemetry] elapsed=$(elapsed)s " *
                            "peak_rss=$(Int(Sys.maxrss() ÷ 1024)) KB")
                end
            end
        end
    end
end

# === Run blocks ===

run_block(LAMBDA_ORGANISM, LAMBDA_ACCESSION, LAMBDA_COVERAGES, LAMBDA_SEEDS)

if RUN_BACTERIAL
    run_block(BACTERIAL_ORGANISM, BACTERIAL_ACCESSION,
        BACTERIAL_COVERAGES, BACTERIAL_SEEDS; bacterial = true)
end

expected_cells = length(LAMBDA_COVERAGES) * length(LAMBDA_SEEDS) +
                 (RUN_BACTERIAL ? length(BACTERIAL_COVERAGES) * length(BACTERIAL_SEEDS) : 0)
println("\nCompleted $(DataFrames.nrow(results))/$(expected_cells) cells")

# === Write results CSV ===

# NOTE: format with no 'Z' directive — 'Z' is a TimeZones timezone token that
# errors on a plain Dates.DateTime ("type DateTime has no field zone"). The
# timestamp is UTC (Dates.now(Dates.UTC)); matches the repo's yyyymmdd_HHMMSS convention.
timestamp = Dates.format(Dates.now(Dates.UTC), "yyyymmdd_HHMMSS")
csv_path = joinpath(results_dir, "headroom_pilot_$(timestamp).csv")
CSV.write(csv_path, results)
println("\nResults written to: $csv_path")
println(results)

# === Summary: headroom per organism/coverage + DP-opportunity cells ===

"""
    mean_or_missing(values) -> Union{Missing, Float64}

Mean over non-missing values, or `missing` if none are numeric.
"""
function mean_or_missing(values)::Union{Missing, Float64}
    numeric = Float64[v for v in values if !ismissing(v)]
    isempty(numeric) && return missing
    return Statistics.mean(numeric)
end

println("\n--- Headroom summary (distance-from-perfect; higher = more DP opportunity) ---")
println("Per organism x coverage (mean over seeds):")
for organism in unique(results.organism)
    org_rows = results[results.organism .== organism, :]
    for coverage in sort(unique(org_rows.coverage))
        sub = org_rows[org_rows.coverage .== coverage, :]
        nga50_mean = mean_or_missing(sub.NGA50)
        gf_mean = mean_or_missing(sub.genome_fraction_quast)
        hr_mean = mean_or_missing(sub.headroom_nga50)
        misasm_mean = mean_or_missing(sub.num_misassemblies)
        fmt(x) = ismissing(x) ? "NA" : string(round(x; digits = 4))
        println("  $(organism) $(coverage)x: NGA50=$(fmt(nga50_mean)) " *
                "GF=$(fmt(gf_mean))% misasm=$(fmt(misasm_mean)) " *
                "headroom_nga50=$(fmt(hr_mean))")
    end
end

println("\nDP-opportunity cells (headroom_nga50 > $(HEADROOM_FLAG_THRESHOLD) OR misassemblies > 0):")
n_flagged = 0
for row in eachrow(results)
    flagged = (!ismissing(row.headroom_nga50) && row.headroom_nga50 > HEADROOM_FLAG_THRESHOLD) ||
              (!ismissing(row.num_misassemblies) && row.num_misassemblies > 0)
    if flagged
        n_flagged += 1
        hr = ismissing(row.headroom_nga50) ? "NA" : string(round(row.headroom_nga50; digits = 4))
        println("  $(row.organism) cov=$(row.coverage)x seed=$(row.seed): " *
                "NGA50=$(row.NGA50)/$(row.ref_size) headroom_nga50=$(hr) " *
                "misasm=$(row.num_misassemblies) contigs=$(row.contigs)")
    end
end
if n_flagged == 0
    println("  (none — every scored cell is at/near the greedy ceiling; no DP headroom found)")
else
    println("  => $(n_flagged) cell(s) where greedy is materially below perfect — the regime to target for a DP strategy.")
end

println("\n=== Track-A Headroom Pilot complete (UTC): $(Dates.now(Dates.UTC)) ===")
