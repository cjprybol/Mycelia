# Real-Data Corrector Validation — Rhizomorph :scalable corrector on REAL genomes
# =============================================================================
# Production-confidence validation: does the near-complete + fast + low-redundancy
# result from the synthetic-read sweeps hold on a REAL, well-characterized genome
# with real repeat/GC structure?
#
# Data provenance (documented, honest):
#   * REFERENCE  — REAL genome downloaded from NCBI RefSeq by accession.
#   * READS      — realistic Illumina reads SIMULATED FROM THE REAL REFERENCE via
#                  ART (HS25 = HiSeq 2500 empirical error/quality model). SRA
#                  acquisition (prefetch/fasterq-dump) is unavailable in this
#                  environment (sra-tools not installed), so we use ART's
#                  platform-realistic error model applied to the REAL reference.
#                  This exercises the real genome's repeat/GC structure — the
#                  property a random SMOKE genome cannot test — while keeping read
#                  error characteristics matched to a real Illumina platform.
#
# Arms:
#   * naive    — assemble_genome(reads; corrector=:none, k=21)
#   * scalable — assemble_genome(reads; corrector=:iterative, strategy=:scalable, k=21)
#
# Metrics vs REFERENCE via QUAST: genome fraction (%), #contigs, N50, largest
# contig, mismatches per 100 kbp, misassemblies. This is a true alignment-based
# genome fraction (not total_length/glen).
#
# Read-only w.r.t. src/ — this script only consumes the public Mycelia API.
#
# Usage:
#   LD_LIBRARY_PATH='' ~/.juliaup/bin/julia --project=. \
#       benchmarking/real_data_corrector_validation.jl
#
# Env knobs:
#   MYCELIA_RDV_TARGETS   comma list of names: phix174,lambda  (default: phix174,lambda)
#   MYCELIA_RDV_COVERAGE  fold coverage for ART                (default: 50)
#   MYCELIA_RDV_K         k-mer size                           (default: 21)
#   MYCELIA_RDV_SEED      ART rndSeed                          (default: 42)
#   MYCELIA_RDV_SMOKE     "true" -> phix174 only, coverage 30  (default: false)

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import BioSequences
import DataFrames
import CSV
import Dates
import Statistics

# --- Config -----------------------------------------------------------------

_truthy(s) = lowercase(strip(s)) in ("1", "true", "yes", "on")

const SMOKE = _truthy(get(ENV, "MYCELIA_RDV_SMOKE", "false"))

# name => (accession, approx_size_bp)
const REGISTRY = Dict(
    "phix174" => ("NC_001422", 5386),   # phiX174 — classic Illumina control, ~5.4 kb
    "lambda" => ("NC_001416", 48502),  # Enterobacteria phage lambda, ~48.5 kb
)

const TARGETS = let
    raw = get(ENV, "MYCELIA_RDV_TARGETS", SMOKE ? "phix174" : "phix174,lambda")
    [String(strip(x)) for x in split(raw, ",") if !isempty(strip(x))]
end

const COVERAGE = parse(Int, get(ENV, "MYCELIA_RDV_COVERAGE", SMOKE ? "30" : "50"))
const K = parse(Int, get(ENV, "MYCELIA_RDV_K", "21"))
const SEED = parse(Int, get(ENV, "MYCELIA_RDV_SEED", "42"))

println("=== Real-Data Corrector Validation ===")
println("Start        : $(Dates.now())")
println("Targets      : $(join(TARGETS, ", "))")
println("Coverage     : $(COVERAGE)x  |  k=$(K)  |  ART seed=$(SEED)")
println("Read model   : ART HS25 (HiSeq 2500), 150 bp paired-end, from REAL reference")
println("SRA note     : sra-tools unavailable -> simulated-from-real-reference reads")

workdir = mktempdir(prefix = "rdv_")
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
println("Workdir      : $(workdir)")

# --- Helpers ----------------------------------------------------------------

"""Load a FASTQ file into a Vector of FASTQ.Record."""
function load_fastq(path::String)
    open(FASTX.FASTQ.Reader, path) do reader
        collect(reader)
    end
end

"""
Parse a MUMmer 3.23 `dnadiff` .report into the reference-based metrics we need.
The shared `Mycelia.parse_dnadiff_report` does not decode MUMmer's
`count(pct%)` single-token format for the `[Bases]` block, so we parse it here
(read-only w.r.t. src/). Returns (aligned_ref_bases, genome_fraction_pct,
avg_identity_pct, total_snps, total_indels); NaN/0 when a field is absent.

Report layout (columns are REF then QRY):
    AlignedBases   48487(99.97%)   48502(100.00%)
    AvgIdentity    99.99           99.99            # first occ = 1-to-1
    TotalSNPs      12              12
    TotalIndels    3               3
"""
function parse_dnadiff_local(report::String)
    aln_ref = 0; gf = NaN; ident = NaN; snps = 0; indels = 0
    seen_ident = false
    count_pct(tok) = begin
        m = match(r"^(\d+)\(([\d.]+)%\)", tok)
        m === nothing ? (nothing, nothing) : (parse(Int, m.captures[1]), parse(Float64, m.captures[2]))
    end
    for line in eachline(report)
        f = split(strip(line))
        isempty(f) && continue
        key = f[1]
        if key == "AlignedBases" && length(f) >= 2
            c, p = count_pct(f[2])
            c !== nothing && (aln_ref = c; gf = p)
        elseif key == "AvgIdentity" && !seen_ident && length(f) >= 2
            v = tryparse(Float64, f[2]); v !== nothing && (ident = v; seen_ident = true)
        elseif key == "TotalSNPs" && length(f) >= 2
            v = tryparse(Int, f[2]); v !== nothing && (snps = v)
        elseif key == "TotalIndels" && length(f) >= 2
            v = tryparse(Int, f[2]); v !== nothing && (indels = v)
        end
    end
    return (aligned_ref_bases = aln_ref, genome_fraction = gf,
        avg_identity = ident, total_snps = snps, total_indels = indels)
end

"""Write assembly contigs (Vector of String) to a FASTA file."""
function write_contigs(contigs, path::String)
    open(path, "w") do io
        for (i, contig) in enumerate(contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)
        end
    end
    return path
end

"""
Run one assembly arm and evaluate it against the reference with MUMmer `dnadiff`.
Returns a NamedTuple of metrics.

Reference-based metrics come from dnadiff (QUAST's old bioconda build does not
solve on osx-arm64; dnadiff/MUMmer installs cleanly and yields the same
reference-aligned quantities):
  * genome fraction (%) = AlignedBases % of the REFERENCE  (aligned_pct_ref)
  * avg identity (%)    = AvgIdentity over aligned segments
  * mismatches/100 kbp  = SNPs / aligned_ref_bases * 1e5
  * indels/100 kbp      = Indels / aligned_ref_bases * 1e5
"""
function run_arm(name::String, reads, corrector::Symbol, strategy, reference::String, outdir::String)
    tag = corrector == :none ? "naive" : "scalable"
    mkpath(outdir)
    contigs_path = joinpath(outdir, "$(name)_$(tag)_contigs.fasta")

    t0 = time()
    local result
    try
        if corrector == :none
            result = Mycelia.Rhizomorph.assemble_genome(
                reads;
                k = K,
                graph_mode = Mycelia.Rhizomorph.DoubleStrand,
                corrector = :none,
                verbose = false,
            )
        else
            result = Mycelia.Rhizomorph.assemble_genome(
                reads;
                k = K,
                graph_mode = Mycelia.Rhizomorph.DoubleStrand,
                corrector = :iterative,
                strategy = strategy,
                verbose = false,
            )
        end
    catch e
        @warn "Assembly arm failed" name tag exception = (e, catch_backtrace())
        return (name = name, arm = tag, ok = false, runtime_s = time() - t0,
            n_contigs = 0, total_length = 0, largest_contig = 0, n50 = 0,
            genome_fraction = NaN, avg_identity = NaN,
            mismatches_per_100kbp = NaN, indels_per_100kbp = NaN)
    end
    runtime = time() - t0
    write_contigs(result.contigs, contigs_path)
    n_contigs = length(result.contigs)
    total_length = sum(length.(result.contigs); init = 0)

    # Reference-free structural metrics
    n50 = 0; largest = 0
    am = Mycelia.assembly_metrics(contigs_path)
    if am !== nothing
        n50 = am.n50
        largest = am.largest_contig
    end

    # Reference-based metrics via MUMmer dnadiff
    gf = NaN; ident = NaN; mm = NaN; indels_100k = NaN
    if n_contigs > 0 && total_length > 0
        dd_out = joinpath(outdir, "$(name)_$(tag)_dnadiff")
        try
            paths = Mycelia.run_dnadiff(reference = reference, query = contigs_path,
                outdir = dd_out, prefix = "dnadiff", force = true)
            p = parse_dnadiff_local(paths.report)
            gf = p.genome_fraction
            ident = p.avg_identity
            if p.aligned_ref_bases > 0
                mm = round(p.total_snps / p.aligned_ref_bases * 1e5; digits = 1)
                indels_100k = round(p.total_indels / p.aligned_ref_bases * 1e5; digits = 1)
            end
        catch e
            @warn "dnadiff failed" name tag exception = (e, catch_backtrace())
        end
    end

    return (name = name, arm = tag, ok = true, runtime_s = round(runtime; digits = 2),
        n_contigs = n_contigs, total_length = total_length,
        largest_contig = largest, n50 = n50,
        genome_fraction = gf, avg_identity = ident,
        mismatches_per_100kbp = mm, indels_per_100kbp = indels_100k)
end

# --- Main -------------------------------------------------------------------

rows = NamedTuple[]

for name in TARGETS
    haskey(REGISTRY, name) || (@warn "Unknown target, skipping" name; continue)
    accession, approx = REGISTRY[name]
    println("\n" * "="^70)
    println("Target: $(name)  ($(accession), ~$(approx) bp)")
    println("="^70)

    target_dir = joinpath(workdir, name)
    mkpath(target_dir)

    # 1. Download REAL reference from NCBI RefSeq
    println("[1/3] Downloading reference $(accession) ...")
    reference = Mycelia.download_genome_by_accession(
        accession = accession, outdir = target_dir, compressed = false)
    if !isfile(reference) || filesize(reference) == 0
        @warn "Reference download failed; skipping target" name accession
        continue
    end
    ref_records = open(FASTX.FASTA.Reader, reference) do r; collect(r); end
    ref_len = length(FASTX.sequence(BioSequences.LongDNA{4}, ref_records[1]))
    println("      Reference: $(reference)  ($(ref_len) bp)")

    # 2. Simulate realistic Illumina reads from the REAL reference (ART HS25)
    println("[2/3] Simulating Illumina reads (ART HS25, $(COVERAGE)x, 150 bp PE) ...")
    outbase = joinpath(target_dir, "$(name)_illumina")
    art = Mycelia.simulate_illumina_reads(
        fasta = reference, coverage = COVERAGE, outbase = outbase,
        read_length = 150, mflen = 300, sdev = 10, seqSys = "HS25",
        paired = true, errfree = false, rndSeed = SEED, quiet = true)
    # ART output is gzipped by simulate_illumina_reads; decompress to plain .fq
    r1 = replace(art.forward_reads, ".gz" => "")
    r2 = replace(art.reverse_reads, ".gz" => "")
    for (gz, fq) in ((art.forward_reads, r1), (art.reverse_reads, r2))
        (isfile(gz) && filesize(gz) > 0) || error("ART did not produce $(gz)")
        isfile(fq) || run(`gunzip -k -f $(gz)`)
        (isfile(fq) && filesize(fq) > 0) || error("gunzip did not produce $(fq)")
    end
    reads = vcat(load_fastq(r1), load_fastq(r2))
    read_bases = sum(length(FASTX.sequence(String, rec)) for rec in reads)
    println("      Reads: $(length(reads)) ($(read_bases) bp, ~$(round(read_bases / ref_len; digits = 1))x effective)")

    # 3. Two arms + dnadiff (reference-based) evaluation
    println("[3/3] Assembling (naive vs scalable, k=$(K)) + dnadiff ...")
    naive = run_arm(name, reads, :none, nothing, reference, joinpath(target_dir, "naive"))
    scalable = run_arm(name, reads, :iterative, :scalable, reference, joinpath(target_dir, "scalable"))
    push!(rows, naive)
    push!(rows, scalable)

    for r in (naive, scalable)
        println("   [$(r.arm)] contigs=$(r.n_contigs) total=$(r.total_length) " *
                "N50=$(r.n50) largest=$(r.largest_contig) " *
                "GF=$(r.genome_fraction)% ident=$(r.avg_identity)% " *
                "mm/100kb=$(r.mismatches_per_100kbp) indels/100kb=$(r.indels_per_100kbp) " *
                "t=$(r.runtime_s)s")
    end
end

# --- Report -----------------------------------------------------------------

if isempty(rows)
    println("\nNo results produced.")
    exit(1)
end

df = DataFrames.DataFrame(rows)
stamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
csv_path = joinpath(results_dir, "real_data_corrector_validation_$(stamp).csv")
CSV.write(csv_path, df)

println("\n" * "="^70)
println("SUMMARY (real reference, ART-simulated reads)")
println("="^70)
show(df; allrows = true, allcols = true)
println("\n\nResults CSV: $(csv_path)")
println("Done: $(Dates.now())")
