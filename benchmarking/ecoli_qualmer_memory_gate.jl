# E. coli Qualmer-Path Memory Gate — does the CORRECTOR complete under a
# SPAdes-class memory ceiling with the composed prefilter + aggregate storage?
#
# Purpose
# -------
# FASTQ input routes Rhizomorph assembly/correction through the QUALMER graph,
# whose per-occurrence quality evidence is O(reads x read_len) ~ O(genome x
# coverage). E. coli @30-50x hit 74-200 GB+ and OOM'd on the lr6 95 GB node. The
# PR #425 coverage prefilter (qualmer_prefilter_min_count) removes singleton
# error-k-mer vertices (~2x) but leaves the true-k-mer per-occurrence residual.
# td-n8ax composes the prefilter WITH aggregate quality storage
# (qualmer_memory_profile=:ultralight_quality — a coverage counter + exact
# per-position mean Phred per DISTINCT k-mer, O(distinct), no per-observation
# tracking) so the residual collapses.
#
# This is the qualmer-path analogue of ecoli_lightweight_crashtest.jl (which
# gates only the NON-quality k-mer path's :lightweight profile). It assembles the
# single E. coli cell ONCE through corrector=:iterative with the COMPOSED config,
# crash-guarded, records whole-process peak Sys.maxrss(), and applies a hard
# memory GATE.
#
# GATE (td-n8ax acceptance)
# -------------------------
#   PASS  iff  the run COMPLETES  AND  peak RSS <= MEMORY_CEILING_GB (default 32).
#   Also reports whether it cleared the ~16 GB SPAdes-class aspiration.
# Was 74-200 GB+ (non-terminating). QUAST scoring (if it completes) confirms the
# correction quality is preserved (genome fraction / NGA50).
#
# External tools + host
# ---------------------
# Requires Conda-backed art_illumina (simulation) + QUAST (scoring), set up on
# Lawrencium, NOT the laptop. This script CANNOT run locally — it is static-
# verified (parse + field-name checks) here and DISPATCHED to lr6 for the real
# gate. LD_LIBRARY_PATH / thread caps / JULIA_DEPOT_PATH are the caller's job.
#
# Env knobs (all optional):
#   MYCELIA_ECOLI_COVERAGE   default "30"   — coverage tier (the gate cell)
#   MYCELIA_ECOLI_SEED       default "42"   — ART seed
#   MYCELIA_MODE_K           default "31"   — base k (the iterative corrector ladders from here)
#   MYCELIA_QUALMER_PROFILE  default "ultralight_quality" — aggregate storage profile
#   MYCELIA_PREFILTER_MIN    default "2"    — coverage prefilter floor (composed lever)
#   MYCELIA_MEMORY_CEILING_GB default "32"  — hard gate ceiling
#   MYCELIA_RUN_QUAST        default "true" — QUAST accuracy scoring
#   MYCELIA_QUAST_MIN_CONTIG default "500"  — QUAST --min-contig
#   MYCELIA_GATE_HEARTBEAT_S default "300"  — seconds between flushed memory lines (0 = off)
#   MYCELIA_GATE_VERBOSE     default "false" — AssemblyConfig verbose (assembler stage logging)
#
# Usage (on Lawrencium lr6, see PR for the exact SLURM launch):
#   julia --project=. benchmarking/ecoli_qualmer_memory_gate.jl

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

const ORGANISM::String = "Ecoli_K12_MG1655"
const ACCESSION::String = "U00096.3"  # RefSeq equivalent: NC_000913.3

const COVERAGE::Int = parse(Int, get(ENV, "MYCELIA_ECOLI_COVERAGE", "30"))
const SEED::Int = parse(Int, get(ENV, "MYCELIA_ECOLI_SEED", "42"))
const K::Int = parse(Int, get(ENV, "MYCELIA_MODE_K", "31"))
const QUALMER_PROFILE::Symbol = Symbol(get(ENV, "MYCELIA_QUALMER_PROFILE", "ultralight_quality"))
const PREFILTER_MIN::Int = parse(Int, get(ENV, "MYCELIA_PREFILTER_MIN", "2"))
const MEMORY_CEILING_GB::Float64 = parse(Float64, get(ENV, "MYCELIA_MEMORY_CEILING_GB", "32"))
const ASPIRATION_GB::Float64 = 16.0

const RUN_QUAST::Bool = get(ENV, "MYCELIA_RUN_QUAST", "true") in ["1", "true", "yes"]
const QUAST_MIN_CONTIG::Int = parse(Int, get(ENV, "MYCELIA_QUAST_MIN_CONTIG", "500"))
const HEARTBEAT_S::Int = parse(Int, get(ENV, "MYCELIA_GATE_HEARTBEAT_S", "300"))
const ASSEMBLY_VERBOSE::Bool = get(ENV, "MYCELIA_GATE_VERBOSE", "false") in [
    "1", "true", "yes"]

# === Progress that survives a Slurm kill ===
# Julia block-buffers stdout when it is redirected to a file, so a job killed for
# time or memory loses everything it printed since the last buffer flush: on
# Lawrencium (jobs 28071705-07) the .out files stopped at the launcher's RUN line
# and the phase reached was unknowable. Every progress line is therefore flushed,
# and a background task prints peak and live memory every HEARTBEAT_S seconds
# (0 disables it). The task waits with Libc.systemsleep, an OS sleep, NOT sleep():
# sleep() is a libuv timer, and while the main thread is inside the single
# blocking assemble_genome call the event loop is not serviced, so a sleep()-based
# heartbeat never fires (checked: zero lines in 6 s against a busy main thread,
# versus one line per second with systemsleep, both surviving SIGKILL). Because
# systemsleep blocks its thread, the heartbeat only starts with
# Threads.nthreads() > 1; the gate sbatch sets JULIA_NUM_THREADS from
# --cpus-per-task=16.

function gate_println(args...)::Nothing
    println(args...)
    flush(stdout)
    return nothing
end

function memory_line(tag::AbstractString)::String
    maxrss_gb = round(Sys.maxrss() / 1024^3; digits = 3)
    live_gb = round(Base.gc_live_bytes() / 1024^3; digits = 3)
    return "[$(tag) $(Dates.now(Dates.UTC))] maxrss=$(maxrss_gb) GB gc_live=$(live_gb) GB"
end

const HEARTBEAT_TASK = if HEARTBEAT_S > 0 && Threads.nthreads() > 1
    Threads.@spawn begin
        while true
            Libc.systemsleep(HEARTBEAT_S)
            gate_println(memory_line("heartbeat"))
            flush(stderr)
        end
    end
else
    HEARTBEAT_S > 0 &&
        @warn "Memory heartbeat disabled: needs Threads.nthreads() > 1 " *
              "(got $(Threads.nthreads())); phase lines are still flushed."
    nothing
end

gate_println("=== E. coli Qualmer-Path Memory Gate (td-n8ax) ===")
gate_println("Start time (UTC): $(Dates.now(Dates.UTC))")
gate_println("Organism: $ORGANISM ($ACCESSION)")
gate_println("Cell: coverage=$(COVERAGE)x seed=$(SEED), base k=$K")
gate_println("COMPOSED config: corrector=:iterative, qualmer_memory_profile=:$(QUALMER_PROFILE), " *
             "qualmer_prefilter_min_count=$(PREFILTER_MIN), graph_mode=DoubleStrand")
gate_println("GATE: completes AND peak RSS <= $(MEMORY_CEILING_GB) GB (aspiration <= $(ASPIRATION_GB) GB); was 74-200 GB+")

# Under the sbatch wrapper (which sets QUALMER_GATE_HEAD) this directory survives
# process exit so the wrapper can save the contigs FASTA and QUAST report; the
# wrapper owns its lifetime. Run by hand, Julia removes it at exit as usual.
benchmark_dir = mktempdir(prefix = "ecoli_qualmer_memory_gate_";
    cleanup = !haskey(ENV, "QUALMER_GATE_HEAD"))
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
gate_println("Working directory: $benchmark_dir")

# === QUAST report.tsv parsing (reused verbatim from the crashtest) ===

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

gate_println("\n--- Phase 1: download E. coli reference ---")
ref_dir = joinpath(benchmark_dir, "reference")
mkpath(ref_dir)
reference_path = Mycelia.download_genome_by_accession(
    accession = ACCESSION, outdir = ref_dir, compressed = false)
if !isfile(reference_path) || filesize(reference_path) == 0
    error("Failed to download E. coli reference ($ACCESSION) to $ref_dir")
end
gate_println("Reference: $reference_path ($(filesize(reference_path)) bytes)")

# === Phase 2: simulate reads ===

gate_println("\n--- Phase 2: simulate reads (coverage=$(COVERAGE)x seed=$(SEED)) ---")
reads_dir = joinpath(benchmark_dir, "reads")
mkpath(reads_dir)
outbase = joinpath(reads_dir, "$(ORGANISM)_cov$(COVERAGE)x_seed$(SEED)")
sim = Mycelia.simulate_illumina_reads(
    fasta = reference_path, coverage = COVERAGE, rndSeed = SEED,
    read_length = 150, paired = true, quiet = true, outbase = outbase)
records = load_reads(sim)
if isempty(records)
    error("No reads simulated for coverage=$(COVERAGE)x seed=$(SEED)")
end
gate_println("Loaded $(length(records)) reads")

# === Phase 3: assemble ONCE with the COMPOSED corrector config (crash-guarded) ===

gate_println("\n--- Phase 3: assemble (corrector=:iterative + prefilter + aggregate, crash-guarded) ---")

# The corrector path: use_quality_scores=true (FASTQ default) routes through the
# qualmer graph; qualmer_memory_profile + qualmer_prefilter_min_count are the two
# composed memory levers under test. DoubleStrand is the corrector's default mode.
config = Mycelia.Rhizomorph.AssemblyConfig(;
    k = K,
    graph_mode = Mycelia.Rhizomorph.DoubleStrand,
    use_quality_scores = true,
    corrector = :iterative,
    qualmer_memory_profile = QUALMER_PROFILE,
    qualmer_prefilter_min_count = PREFILTER_MIN,
    verbose = ASSEMBLY_VERBOSE)

completed::Bool = false
# failure_kind distinguishes a genuine OOM (the memory result this gate exists to
# measure) from a plain code error (MethodError/etc). Both would otherwise print the
# identical "did NOT complete" verdict and mislead an operator into chasing memory
# when the fault is code — so we branch on OutOfMemoryError and use a distinct exit.
failure_kind::Symbol = :none
elapsed_s::Union{Missing, Float64} = missing
n_contigs::Union{Missing, Int} = missing
contigs_path::String = joinpath(benchmark_dir, "$(ORGANISM)_cov$(COVERAGE)x_seed$(SEED)_contigs.fasta")

GC.gc(true)
gate_println(memory_line("before-assemble"))
t0 = time()
try
    result = Mycelia.Rhizomorph.assemble_genome(records, config)
    global elapsed_s = round(time() - t0; digits = 3)
    global n_contigs = length(result.contigs)
    global completed = true
    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)
        end
    end
    gate_println("  Assembly COMPLETED: $(n_contigs) contigs in $(elapsed_s)s")
catch e
    global elapsed_s = round(time() - t0; digits = 3)
    global completed = false
    global failure_kind = e isa OutOfMemoryError ? :oom : :error
    @warn "Assembly FAILED ($(failure_kind)) for E. coli $(COVERAGE)x seed=$(SEED) after " *
          "$(elapsed_s)s — recording, not aborting" exception = (e, catch_backtrace())
end

# Whole-process peak RSS (monotonic; captured after the assemble attempt).
peak_rss_kb::Int = Int(Sys.maxrss() ÷ 1024)
peak_rss_gb::Float64 = round(peak_rss_kb / (1024 * 1024); digits = 3)
gate_println("  Peak Sys.maxrss(): $(peak_rss_gb) GB ($(peak_rss_kb) KB)")

# === Phase 4: QUAST aligned-truth scoring (only if it completed) ===

nga50::Union{Missing, Float64} = missing
num_misassemblies::Union{Missing, Float64} = missing
dup_ratio::Union{Missing, Float64} = missing
genome_fraction::Union{Missing, Float64} = missing

if completed && RUN_QUAST && !ismissing(n_contigs) && n_contigs > 0 &&
   filesize(contigs_path) > 0
    gate_println("\n--- Phase 4: QUAST aligned-truth scoring ---")
    quast_dir = joinpath(benchmark_dir, "quast")
    try
        Mycelia.run_quast(contigs_path; outdir = quast_dir,
            reference = reference_path, min_contig = QUAST_MIN_CONTIG)
        report_tsv = joinpath(quast_dir, "report.tsv")
        global nga50 = parse_quast_metric(report_tsv, "NGA50")
        global num_misassemblies = parse_quast_metric(report_tsv, "# misassemblies")
        global dup_ratio = parse_quast_metric(report_tsv, "Duplication ratio")
        global genome_fraction = parse_quast_metric(report_tsv, "Genome fraction (%)")
        gate_println("  NGA50=$(nga50), misasm=$(num_misassemblies), dup_ratio=$(dup_ratio), GF=$(genome_fraction)%")
    catch e
        @warn "QUAST failed — leaving accuracy metrics missing" exception = (
            e, catch_backtrace())
    end
elseif !completed
    gate_println("\n--- Phase 4: skipped (assembly did not complete) ---")
end

# === Phase 5: write results CSV ===

timestamp = Dates.format(Dates.now(Dates.UTC), "yyyymmdd_HHMMSS")
csv_path = joinpath(results_dir, "ecoli_qualmer_memory_gate_$(timestamp).csv")
results = DataFrames.DataFrame(
    organism = String[ORGANISM],
    qualmer_memory_profile = String[string(QUALMER_PROFILE)],
    prefilter_min_count = Int[PREFILTER_MIN],
    coverage = Int[COVERAGE],
    seed = Int[SEED],
    completed = Bool[completed],
    elapsed_s = Union{Missing, Float64}[elapsed_s],
    peak_rss_kb = Int[peak_rss_kb],
    peak_rss_gb = Float64[peak_rss_gb],
    memory_ceiling_gb = Float64[MEMORY_CEILING_GB],
    n_contigs = Union{Missing, Int}[n_contigs],
    NGA50 = Union{Missing, Float64}[nga50],
    num_misassemblies = Union{Missing, Float64}[num_misassemblies],
    dup_ratio = Union{Missing, Float64}[dup_ratio],
    genome_fraction_quast = Union{Missing, Float64}[genome_fraction],
    # Provenance set by the sbatch wrapper (empty when run by hand): the checkout
    # HEAD, the sha256 of the env Manifest.toml the run actually used (includes the
    # Mycelia develop-path entry), and the sha256 of the supplied source Manifest
    # (comparable with a checkout's file; empty in resolve mode).
    git_head = String[get(ENV, "QUALMER_GATE_HEAD", "")],
    manifest_sha256 = String[get(ENV, "QUALMER_GATE_MANIFEST_SHA", "")],
    manifest_src_sha256 = String[get(ENV, "QUALMER_GATE_MANIFEST_SRC_SHA", "")])
CSV.write(csv_path, results)
gate_println("\nResults written to: $csv_path")
gate_println(results)

# === Phase 6: PASS/FAIL memory GATE verdict ===

gate_println("\n--- GATE Verdict ---")
under_ceiling = peak_rss_gb <= MEMORY_CEILING_GB
under_aspiration = peak_rss_gb <= ASPIRATION_GB
gate_pass = completed && under_ceiling

if gate_pass
    gf_str = ismissing(genome_fraction) ? "NA" : string(round(genome_fraction; digits = 3))
    aspiration_str = under_aspiration ?
                     "and CLEARED the $(ASPIRATION_GB) GB SPAdes-class aspiration" :
                     "(above the $(ASPIRATION_GB) GB aspiration, but under the ceiling)"
    gate_println("PASS: E. coli $(COVERAGE)x COMPLETED in $(elapsed_s)s at $(peak_rss_gb) GB peak " *
                 "<= $(MEMORY_CEILING_GB) GB ceiling $(aspiration_str). GF=$(gf_str)% " *
                 "($(ismissing(n_contigs) ? "NA" : n_contigs) contigs). Was 74-200 GB+.")
elseif completed && !under_ceiling
    gate_println("FAIL: E. coli $(COVERAGE)x COMPLETED but peak $(peak_rss_gb) GB EXCEEDS the " *
                 "$(MEMORY_CEILING_GB) GB ceiling — memory not yet bounded to target.")
elseif failure_kind == :error
    gate_println("ERROR: E. coli $(COVERAGE)x assembly threw a NON-OOM exception after $(elapsed_s)s " *
                 "at $(peak_rss_gb) GB — this is a CODE failure, NOT a memory-gate result. " *
                 "See the backtrace above. (Do not read this as 'memory unbounded'.)")
else
    gate_println("FAIL: E. coli $(COVERAGE)x did NOT complete — OUT OF MEMORY after $(elapsed_s)s at " *
                 "$(peak_rss_gb) GB peak. Memory not bounded to target.")
end

gate_println("\n=== E. coli Qualmer-Path Memory Gate complete (UTC): $(Dates.now(Dates.UTC)) ===")
# Exit codes: 0 = gate PASS; 1 = memory FAIL (OOM or over-ceiling, a real gate result);
# 2 = CODE error (non-OOM exception — NOT a memory result, needs a code fix).
exit(gate_pass ? 0 : (failure_kind == :error ? 2 : 1))
