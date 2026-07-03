# E. coli Lightweight Crash-Retest — does :lightweight COMPLETE where :full crashed?
#
# Purpose
# -------
# The greedy Rhizomorph assembler CRASHED on a single 10x E. coli (K-12 MG1655,
# ~4.64 Mb) cell under the default memory_profile=:full: >80 min wall, ~1.27
# BILLION allocations, and an out-of-memory crash. Root cause: :full keeps
# per-vertex evidence that is O(coverage x genome), which is intractable at
# bacterial scale. The lightweight edge-support fix (branch
# cp/lightweight-genomefrac) makes per-vertex memory O(distinct-kmers) and
# restored accuracy on Lambda (genome fraction 27.8% -> 99.9%) at ~2.3x less
# memory. This script answers the single feasibility question that fix opens up:
#
#   Does E. coli 10x now COMPLETE with memory_profile=:lightweight instead of
#   crashing?
#
# It assembles the SINGLE 10x / seed-42 E. coli cell ONCE with
# memory_profile=:lightweight, wraps the assemble in try/catch so an OOM/crash is
# REPORTED (not silently fatal), records whether it completed plus elapsed time,
# whole-process peak Sys.maxrss(), and contig count, and — if it completes —
# scores the contigs against the reference with QUAST for NGA50 / misassemblies /
# duplication ratio / genome fraction. A single PASS/FAIL verdict line summarizes.
#
# Why the NON-quality k-mer path (use_quality_scores = false)
# -----------------------------------------------------------
# memory_profile is threaded ONLY into the k-mer assembly path
# (`_assemble_kmer_graph`, the use_quality_scores = false arm). The quality-aware
# qualmer path (auto-selected for FASTQ) does NOT honor memory_profile. To make
# :lightweight actually take effect, the AssemblyConfig here sets
# use_quality_scores = false, exactly matching benchmarking/mode_comparison.jl.
# DoubleStrand is the default graph mode; dedup_revcomp / compact_unitigs are left
# at their defaults (false), so :lightweight is the only efficiency lever engaged.
#
# External tools
# --------------
# Requires Conda-backed external tools (art_illumina for simulation, QUAST for
# scoring) — set up on Lovelace, NOT the laptop. This script CANNOT be run
# locally; it is static-verified (parse + field-name checks) only. QUAST is gated
# behind MYCELIA_RUN_QUAST (default "true"). LD_LIBRARY_PATH / thread caps /
# JULIA_DEPOT_PATH are the caller's job (the run command sets them).
#
# Env knobs (all optional):
#   MYCELIA_RUN_QUAST        default "true"  — run QUAST scoring
#   MYCELIA_MODE_K           default "31"    — k-mer size for assembly
#   MYCELIA_QUAST_MIN_CONTIG default "100"   — QUAST --min-contig
#   MYCELIA_ECOLI_COVERAGE   default "10"    — single coverage tier (the crash cell)
#   MYCELIA_ECOLI_SEED       default "42"    — single ART seed
#
# Usage (on Lovelace, see PR for the exact launch command):
#   julia --project=. benchmarking/ecoli_lightweight_crashtest.jl

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import DataFrames
import CSV
import Dates

# === Configuration ===

# E. coli K-12 MG1655. U00096.3 is the canonical GenBank accession (~4.64 Mb);
# NC_000913.3 is the equivalent RefSeq accession if U00096.3 ever fails to fetch.
const ORGANISM::String = "Ecoli_K12_MG1655"
const ACCESSION::String = "U00096.3"
const TECHNOLOGY::String = "illumina"
const MEMORY_PROFILE::Symbol = :lightweight

const COVERAGE::Int = parse(Int, get(ENV, "MYCELIA_ECOLI_COVERAGE", "10"))
const SEED::Int = parse(Int, get(ENV, "MYCELIA_ECOLI_SEED", "42"))

const RUN_QUAST::Bool = get(ENV, "MYCELIA_RUN_QUAST", "true") in ["1", "true", "yes"]
const K::Int = parse(Int, get(ENV, "MYCELIA_MODE_K", "31"))
const QUAST_MIN_CONTIG::Int = parse(Int, get(ENV, "MYCELIA_QUAST_MIN_CONTIG", "100"))

println("=== E. coli Lightweight Crash-Retest ===")
println("Start time (UTC): $(Dates.now(Dates.UTC))")
println("Organism: $ORGANISM ($ACCESSION), technology: $TECHNOLOGY")
println("Single cell: coverage=$(COVERAGE)x seed=$(SEED), k=$K")
println("memory_profile=:$(MEMORY_PROFILE), use_quality_scores=false (k-mer path — required for memory_profile to take effect)")
println("QUAST scoring: $(RUN_QUAST ? "ENABLED (min_contig=$QUAST_MIN_CONTIG)" : "DISABLED (MYCELIA_RUN_QUAST=false)")")

# Working + results directories
benchmark_dir = mktempdir(prefix = "ecoli_lightweight_crashtest_")
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
println("Working directory: $benchmark_dir")

# === QUAST report.tsv parsing ===
# QUAST writes a tab-separated report.tsv whose first column is the metric name
# and second column is the value for the (single) assembly we pass. Reused
# verbatim from the Track-A pilots.

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
    load_reads(sim) -> Vector{FASTX.FASTQ.Record}

Load the forward and (paired) reverse simulated reads into a single flat record
vector, matching how the k-mer assembler consumes reads (no pairing metadata).
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

# === Phase 1: download E. coli reference ===

println("\n--- Phase 1: download E. coli reference ---")
ref_dir = joinpath(benchmark_dir, "reference")
mkpath(ref_dir)
reference_path = Mycelia.download_genome_by_accession(
    accession = ACCESSION,
    outdir = ref_dir,
    compressed = false
)
if !isfile(reference_path) || filesize(reference_path) == 0
    error("Failed to download E. coli reference ($ACCESSION) to $ref_dir")
end
println("Reference: $reference_path ($(filesize(reference_path)) bytes)")

# === Phase 2: simulate reads for the single crash cell ===

println("\n--- Phase 2: simulate reads (coverage=$(COVERAGE)x seed=$(SEED)) ---")
reads_dir = joinpath(benchmark_dir, "reads")
mkpath(reads_dir)
outbase = joinpath(reads_dir, "$(ORGANISM)_cov$(COVERAGE)x_seed$(SEED)")
sim = Mycelia.simulate_illumina_reads(
    fasta = reference_path,
    coverage = COVERAGE,
    rndSeed = SEED,
    read_length = 150,
    paired = true,
    quiet = true,
    outbase = outbase
)
records = load_reads(sim)
if isempty(records)
    error("No reads simulated for coverage=$(COVERAGE)x seed=$(SEED)")
end
println("Loaded $(length(records)) reads")

# === Phase 3: assemble ONCE with memory_profile=:lightweight (crash-guarded) ===

println("\n--- Phase 3: assemble (memory_profile=:$(MEMORY_PROFILE), crash-guarded) ---")

# Build the AssemblyConfig for the k-mer path. use_quality_scores=false routes
# through `_assemble_kmer_graph`, the only path that honors memory_profile.
# DoubleStrand is the default graph mode; dedup_revcomp / compact_unitigs default
# to false, so :lightweight is the sole efficiency lever under test.
seqtype = Mycelia.Rhizomorph._detect_sequence_type(records)
config = Mycelia.Rhizomorph.AssemblyConfig(;
    k = K,
    sequence_type = seqtype,
    graph_mode = Mycelia.Rhizomorph.DoubleStrand,
    use_quality_scores = false,
    verbose = false,
    memory_profile = MEMORY_PROFILE
)

completed::Bool = false
elapsed_s::Union{Missing, Float64} = missing
n_contigs::Union{Missing, Int} = missing
contigs_path::String = joinpath(benchmark_dir, "$(ORGANISM)_cov$(COVERAGE)x_seed$(SEED)_contigs.fasta")

GC.gc(true)
t0 = time()
try
    result = Mycelia.Rhizomorph.assemble_genome(records, config)
    global elapsed_s = round(time() - t0; digits = 3)
    global n_contigs = length(result.contigs)
    global completed = true
    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)  # contigs are String type
        end
    end
    println("  Assembly COMPLETED: $(n_contigs) contigs in $(elapsed_s)s")
catch e
    global elapsed_s = round(time() - t0; digits = 3)
    global completed = false
    @warn "Assembly FAILED/crashed for E. coli $(COVERAGE)x seed=$(SEED) after $(elapsed_s)s " *
          "(memory_profile=:$(MEMORY_PROFILE)) — recording as a crash, not aborting" exception = (e, catch_backtrace())
end

# Whole-process peak RSS (monotonic; captured after the assemble attempt).
peak_rss_kb::Int = Int(Sys.maxrss() ÷ 1024)
peak_rss_gb::Float64 = round(peak_rss_kb / (1024 * 1024); digits = 3)
println("  Peak Sys.maxrss(): $(peak_rss_gb) GB ($(peak_rss_kb) KB)")

# === Phase 4: QUAST aligned-truth scoring (only if assembly completed) ===

nga50::Union{Missing, Float64} = missing
num_misassemblies::Union{Missing, Float64} = missing
dup_ratio::Union{Missing, Float64} = missing
genome_fraction::Union{Missing, Float64} = missing

if completed && RUN_QUAST && !ismissing(n_contigs) && n_contigs > 0 && filesize(contigs_path) > 0
    println("\n--- Phase 4: QUAST aligned-truth scoring ---")
    quast_dir = joinpath(benchmark_dir, "quast")
    try
        Mycelia.run_quast(
            contigs_path;
            outdir = quast_dir,
            reference = reference_path,
            min_contig = QUAST_MIN_CONTIG
        )
        report_tsv = joinpath(quast_dir, "report.tsv")
        global nga50 = parse_quast_metric(report_tsv, "NGA50")
        global num_misassemblies = parse_quast_metric(report_tsv, "# misassemblies")
        global dup_ratio = parse_quast_metric(report_tsv, "Duplication ratio")
        global genome_fraction = parse_quast_metric(report_tsv, "Genome fraction (%)")
        println("  NGA50=$(nga50), misasm=$(num_misassemblies), dup_ratio=$(dup_ratio), GF=$(genome_fraction)%")
    catch e
        @warn "QUAST failed — leaving accuracy metrics missing" exception = (e,)
    end
elseif completed && !RUN_QUAST
    println("\n--- Phase 4: QUAST scoring DISABLED (MYCELIA_RUN_QUAST=false) ---")
elseif !completed
    println("\n--- Phase 4: skipped (assembly did not complete) ---")
end

# === Phase 5: write results CSV ===

# NOTE: format with no 'Z' directive — 'Z' is a TimeZones timezone token that
# errors on a plain Dates.DateTime. The timestamp is UTC (Dates.now(Dates.UTC));
# matches the repo's yyyymmdd_HHMMSS convention.
timestamp = Dates.format(Dates.now(Dates.UTC), "yyyymmdd_HHMMSS")
csv_path = joinpath(results_dir, "ecoli_lightweight_crashtest_$(timestamp).csv")

results = DataFrames.DataFrame(
    organism = String[ORGANISM],
    memory_profile = String[string(MEMORY_PROFILE)],
    coverage = Int[COVERAGE],
    seed = Int[SEED],
    completed = Bool[completed],
    elapsed_s = Union{Missing, Float64}[elapsed_s],
    peak_rss_kb = Int[peak_rss_kb],
    n_contigs = Union{Missing, Int}[n_contigs],
    NGA50 = Union{Missing, Float64}[nga50],
    num_misassemblies = Union{Missing, Float64}[num_misassemblies],
    dup_ratio = Union{Missing, Float64}[dup_ratio],
    genome_fraction_quast = Union{Missing, Float64}[genome_fraction]
)
CSV.write(csv_path, results)
println("\nResults written to: $csv_path")
println(results)

# === Phase 6: PASS/FAIL verdict ===

println("\n--- Verdict ---")
if completed
    gf_str = ismissing(genome_fraction) ? "NA" : string(round(genome_fraction; digits = 3))
    println("PASS: E. coli $(COVERAGE)x :lightweight COMPLETED in $(elapsed_s)s at $(peak_rss_gb)GB, " *
            "GF=$(gf_str)% ($(ismissing(n_contigs) ? "NA" : n_contigs) contigs) — " *
            "no crash where memory_profile=:full previously OOM'd.")
else
    println("FAIL: E. coli $(COVERAGE)x :lightweight FAILED/crashed after $(elapsed_s)s at $(peak_rss_gb)GB peak — " *
            "did NOT complete (see the warning above for the exception).")
end

println("\n=== E. coli Lightweight Crash-Retest complete (UTC): $(Dates.now(Dates.UTC)) ===")
