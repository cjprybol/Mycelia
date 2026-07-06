# Rhizomorph Correction-Core Validation Sweep
# =============================================
#
# Assembly-vs-assembly, error-rate x read-regime sweep that validates the
# corrected Rhizomorph correction core by comparing two assembly ARMS on the
# SAME simulated reads:
#
#   * naive     — assemble_genome(reads; corrector=:none)      (baseline)
#   * iterative — assemble_genome(reads; corrector=:iterative) (correction core)
#
# Both arms run on the DoubleStrand graph mode (NOT :canonical). Rationale
# (see src/rhizomorph/assembly.jl::_assemble_with_iterative_corrector): the
# naive contig path is INVALID under :canonical, so a canonical baseline would
# not be apples-to-apples. DoubleStrand is the auto-detected DNA mode for both
# the naive path and the iterative corrector's re-assembly, so pinning both arms
# to DoubleStrand keeps the comparison honest.
#
# Axes swept:
#   * error rate  — MYCELIA_RGV_ERR      (default 0.01,0.05,0.10)
#   * read regime — MYCELIA_RGV_READLEN  (default 150,5000)
#       readlen <= 500  -> "short-low-error"  regime, Illumina error model
#       readlen  > 500  -> "long-high-error"  regime, Nanopore error model
#
# SCALE-ASSERTION GUARD (rhizomorph_scale_guard.jl): a VERDICT is only printed
# when effective coverage x effective genome length exceeds a floor. Below the
# floor the harness emits a SMOKE-ONLY notice, so a toy run can never be quoted
# as validation. This is what makes the harness safe to smoke-test locally at
# tiny scale while the real validation runs on Lovelace.
#
# Usage:
#   # Local smoke test (synthetic 2 kb genome, no network, SMOKE-ONLY):
#   MYCELIA_RGV_SMOKE=true LD_LIBRARY_PATH='' \
#     julia --project=. benchmarking/rhizomorph_correction_validation_sweep.jl
#
#   # Full validation run (Lambda phage, real coverage, VERDICT) on Lovelace:
#   MYCELIA_RGV_COVERAGE=30 MYCELIA_RUN_EXTERNAL=true LD_LIBRARY_PATH='' \
#     julia --project=. benchmarking/rhizomorph_correction_validation_sweep.jl
#
# Environment variables:
#   MYCELIA_RGV_SMOKE          truthy -> synthetic genome, no download (default false)
#   MYCELIA_RGV_ACCESSION      reference accession (default NC_001416, Lambda phage)
#   MYCELIA_RGV_SMOKE_GENOME_LEN synthetic genome length in smoke mode (default 2000)
#   MYCELIA_RGV_ERR            comma-separated error rates (default 0.01,0.05,0.10)
#   MYCELIA_RGV_READLEN        comma-separated read lengths (default 150,5000)
#   MYCELIA_RGV_COVERAGE       target fold coverage per cell (default 30; smoke default 10)
#   MYCELIA_RGV_K              assembly k-mer size (default 21)
#   MYCELIA_RGV_SEED           RNG seed for reproducibility (default 42)
#   MYCELIA_RGV_SCALE_FLOOR    override the scale-guard floor in bases
#   MYCELIA_RUN_EXTERNAL       truthy -> run QUAST (degrades gracefully if absent)

import Pkg
if isinteractive()
    Pkg.activate(joinpath(@__DIR__, ".."))
end

import Mycelia
import FASTX
import BioSequences
import DataFrames
import CSV
import Dates
import Statistics
import Random

# Pure, dependency-free scale guard (shared with the unit test).
include(joinpath(@__DIR__, "rhizomorph_scale_guard.jl"))

# === Configuration parsing =================================================

_truthy(s) = lowercase(strip(s)) in ("1", "true", "yes", "on")

function _parse_float_list(s, default::Vector{Float64})
    isempty(strip(s)) && return default
    return Float64[parse(Float64, strip(x)) for x in split(s, ",") if !isempty(strip(x))]
end

function _parse_int_list(s, default::Vector{Int})
    isempty(strip(s)) && return default
    return Int[parse(Int, strip(x)) for x in split(s, ",") if !isempty(strip(x))]
end

"""
Map a read length to a (label, technology) read-regime pair. Short reads model
low-error Illumina chemistry; long reads model high-error Nanopore chemistry.
The per-base error is still overridden by the swept error rate — the technology
selects the error-TYPE mix (substitution- vs indel-dominated) and quality model.
"""
function regime_for_readlen(readlen::Int)
    return readlen <= 500 ? ("short-low-error", :illumina) : ("long-high-error", :nanopore)
end

# === Reference acquisition =================================================

"""
Acquire the reference sequence for the sweep. In smoke mode a synthetic random
genome is generated in-process (no network, fully self-contained). Otherwise the
reference is downloaded by accession, reusing real_genome_benchmark.jl's
`download_genome_by_accession` path.

Returns `(refseq::BioSequences.LongDNA{4}, ref_fasta_path::String, label::String)`.
"""
function acquire_reference(; smoke::Bool, accession::String, smoke_len::Int, seed::Int, workdir::String)
    if smoke
        rec = Mycelia.random_fasta_record(moltype = :DNA, seed = seed, L = smoke_len)
        ref_path = joinpath(workdir, "synthetic_reference.fasta")
        open(FASTX.FASTA.Writer, ref_path) do w
            write(w, rec)
        end
        refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
        return refseq, ref_path, "synthetic_$(smoke_len)bp"
    end
    genome_dir = joinpath(workdir, accession)
    mkpath(genome_dir)
    ref_path = Mycelia.download_genome_by_accession(
        accession = accession, outdir = genome_dir, compressed = false)
    if !isfile(ref_path) || filesize(ref_path) == 0
        error("Failed to download reference $accession")
    end
    # Concatenate all records into one reference sequence (Lambda is single-record;
    # multi-record references are joined so coverage/length math is well defined).
    records = open(FASTX.FASTA.Reader, ref_path) do reader
        collect(reader)
    end
    isempty(records) && error("No records in downloaded reference $ref_path")
    refseq = BioSequences.LongDNA{4}()
    for rec in records
        append!(refseq, FASTX.sequence(BioSequences.LongDNA{4}, rec))
    end
    return refseq, ref_path, accession
end

# === Read simulation =======================================================

"""
Simulate fixed-length reads for one (error_rate, regime) cell by sampling random
fragments across both strands of the reference and passing each through
`Mycelia.observe` with the swept error rate and regime technology. Returns
`(records::Vector{FASTX.FASTQ.Record}, sampled_bases::Int)` where `sampled_bases`
is the pre-error total used to compute effective coverage.
"""
function simulate_regime_reads(refseq::BioSequences.LongDNA{4}, readlen::Int,
        coverage::Real, error_rate::Float64, tech::Symbol, rng::Random.AbstractRNG)
    glen = length(refseq)
    effective_readlen = min(readlen, glen)
    n_reads = max(1, ceil(Int, coverage * glen / effective_readlen))
    records = Vector{FASTX.FASTQ.Record}()
    sampled_bases = 0
    for i in 1:n_reads
        start = glen == effective_readlen ? 1 : rand(rng, 1:(glen - effective_readlen + 1))
        frag = refseq[start:(start + effective_readlen - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        sampled_bases += effective_readlen
        obs_seq, quals = Mycelia.observe(frag; error_rate = error_rate, tech = tech)
        isempty(obs_seq) && continue
        qstr = String([Char(q + 33) for q in quals])
        push!(records, FASTX.FASTQ.Record("read_$(i)", string(obs_seq), qstr))
    end
    return records, sampled_bases
end

# === One assembly arm ======================================================

"""
Run one assembly ARM (`corrector` = :none or :iterative) on `reads`, both pinned
to DoubleStrand, write contigs to FASTA, and return a metrics NamedTuple. On
failure the arm is recorded with `ok=false` rather than aborting the whole sweep.
"""
function run_arm(reads, corrector::Symbol, k::Int, glen::Int, outdir::String, tag::String)
    contigs_path = joinpath(outdir, "$(tag)_$(corrector)_contigs.fasta")
    t0 = time()
    local result
    try
        result = Mycelia.Rhizomorph.assemble_genome(
            reads;
            k = k,
            graph_mode = Mycelia.Rhizomorph.DoubleStrand,
            corrector = corrector,
            verbose = false
        )
    catch e
        @warn "Assembly arm failed" corrector tag exception = (e, catch_backtrace())
        return (ok = false, corrector = corrector, n_contigs = 0, total_length = 0,
            largest_contig = 0, n50 = 0, genome_fraction = 0.0, runtime_s = time() - t0,
            contigs_path = "")
    end
    runtime = time() - t0

    open(contigs_path, "w") do io
        for (i, contig) in enumerate(result.contigs)
            println(io, ">contig_$(i) length=$(length(contig))")
            println(io, contig)
        end
    end

    n_contigs = length(result.contigs)
    total_length = sum(length.(result.contigs); init = 0)
    largest_contig = 0
    n50 = 0
    if isfile(contigs_path) && filesize(contigs_path) > 0
        try
            m = Mycelia.assembly_metrics(contigs_path)
            if m !== nothing
                largest_contig = m.largest_contig
                n50 = m.n50
            end
        catch e
            @warn "assembly_metrics failed" tag corrector exception = e
        end
    end
    genome_fraction = glen == 0 ? 0.0 : round(total_length / glen * 100; digits = 1)

    return (ok = true, corrector = corrector, n_contigs = n_contigs,
        total_length = total_length, largest_contig = largest_contig, n50 = n50,
        genome_fraction = genome_fraction, runtime_s = round(runtime; digits = 3),
        contigs_path = contigs_path)
end

# === Main sweep ============================================================

function run_sweep()
    smoke = _truthy(get(ENV, "MYCELIA_RGV_SMOKE", "false"))
    accession = get(ENV, "MYCELIA_RGV_ACCESSION", "NC_001416")  # Lambda phage
    smoke_len = parse(Int, get(ENV, "MYCELIA_RGV_SMOKE_GENOME_LEN", "2000"))
    errs = _parse_float_list(get(ENV, "MYCELIA_RGV_ERR", ""), [0.01, 0.05, 0.10])
    readlens = _parse_int_list(get(ENV, "MYCELIA_RGV_READLEN", ""), [150, 5000])
    coverage = parse(Float64, get(ENV, "MYCELIA_RGV_COVERAGE", smoke ? "10" : "30"))
    k = parse(Int, get(ENV, "MYCELIA_RGV_K", "21"))
    seed = parse(Int, get(ENV, "MYCELIA_RGV_SEED", "42"))
    scale_floor = parse(Float64, get(ENV, "MYCELIA_RGV_SCALE_FLOOR", string(SCALE_FLOOR_BASES)))
    run_external = _truthy(get(ENV, "MYCELIA_RUN_EXTERNAL", "false"))

    println("=== Rhizomorph Correction-Core Validation Sweep ===")
    println("Start time     : $(Dates.now())")
    println("Mode           : $(smoke ? "SMOKE (synthetic genome, no network)" : "FULL (download $accession)")")
    println("Error rates    : $errs")
    println("Read lengths   : $readlens")
    println("Coverage       : $(coverage)x")
    println("k              : $k")
    println("Arms           : naive (corrector=:none) vs iterative (corrector=:iterative), both DoubleStrand")
    println("Scale floor    : $(scale_floor) bases")
    println("QUAST          : $(run_external ? "enabled" : "disabled (set MYCELIA_RUN_EXTERNAL=true)")")

    workdir = mktempdir(prefix = "rgv_sweep_")
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)

    println("\n--- Acquiring reference ---")
    refseq, ref_path, ref_label = acquire_reference(
        smoke = smoke, accession = accession, smoke_len = smoke_len, seed = seed, workdir = workdir)
    glen = length(refseq)
    println("Reference: $ref_label ($glen bp)")

    # k must not exceed the shortest read; clamp and warn rather than fail.
    min_readlen = minimum(min.(readlens, glen))
    if k > min_readlen
        @warn "k=$k exceeds shortest effective read length $min_readlen; clamping k to $min_readlen"
        k = min_readlen
    end

    rng = Random.MersenneTwister(seed)

    rows = DataFrames.DataFrame(
        reference = String[], genome_len = Int[], error_rate = Float64[],
        regime = String[], readlen = Int[], tech = String[], target_coverage = Float64[],
        effective_coverage = Float64[], k = Int[], arm = String[], ok = Bool[],
        n_contigs = Int[], total_length = Int[], largest_contig = Int[], n50 = Int[],
        genome_fraction = Float64[], runtime_s = Float64[]
    )

    # Track the minimum effective coverage across cells for the scale guard: the
    # weakest cell governs whether the whole sweep earns a VERDICT.
    min_effective_coverage = Inf

    println("\n--- Sweeping (error_rate x read-regime) x {naive, iterative} ---")
    for err in errs
        for readlen in readlens
            regime, tech = regime_for_readlen(readlen)
            cell_dir = joinpath(workdir, "err$(err)_len$(readlen)")
            mkpath(cell_dir)
            println("\n[cell] err=$err  regime=$regime  readlen=$readlen  tech=$tech")

            reads, sampled_bases = simulate_regime_reads(refseq, readlen, coverage, err, tech, rng)
            eff_cov = glen == 0 ? 0.0 : round(sampled_bases / glen; digits = 2)
            min_effective_coverage = min(min_effective_coverage, eff_cov)
            println("  simulated $(length(reads)) reads, effective coverage $(eff_cov)x")

            for corrector in (:none, :iterative)
                tag = "$(ref_label)_err$(err)_len$(readlen)"
                arm_name = corrector == :none ? "naive" : "iterative"
                res = run_arm(reads, corrector, k, glen, cell_dir, tag)
                push!(rows, (
                    reference = ref_label, genome_len = glen, error_rate = err,
                    regime = regime, readlen = readlen, tech = String(tech),
                    target_coverage = coverage, effective_coverage = eff_cov, k = k,
                    arm = arm_name, ok = res.ok, n_contigs = res.n_contigs,
                    total_length = res.total_length, largest_contig = res.largest_contig,
                    n50 = res.n50, genome_fraction = res.genome_fraction,
                    runtime_s = res.runtime_s))
                println("    $(rpad(arm_name, 9)) -> ok=$(res.ok) contigs=$(res.n_contigs) " *
                        "total=$(res.total_length)bp largest=$(res.largest_contig) " *
                        "n50=$(res.n50) frac=$(res.genome_fraction)% $(res.runtime_s)s")
            end

            # Optional QUAST validation for this cell (both arms vs reference).
            if run_external
                cell_contigs = String[joinpath(cell_dir, f) for f in readdir(cell_dir)
                                      if endswith(f, "_contigs.fasta") && filesize(joinpath(cell_dir, f)) > 0]
                if !isempty(cell_contigs)
                    try
                        Mycelia.run_quast(cell_contigs;
                            outdir = joinpath(cell_dir, "quast"),
                            reference = ref_path,
                            min_contig = max(50, glen ÷ 10))
                        println("  QUAST: $(joinpath(cell_dir, "quast"))")
                    catch e
                        @warn "QUAST unavailable / failed for this cell — continuing without it" exception = e
                    end
                end
            end
        end
    end

    # === Tidy per-(err, regime) table ===
    println("\n=== Per-(error_rate, regime) summary: naive vs iterative ===")
    _fmt = (v) -> rpad(string(v), 10)
    println(join(_fmt.(["err", "regime", "readlen", "arm", "contigs", "total_bp",
        "largest", "n50", "frac_%", "runtime_s"]), ""))
    for r in eachrow(sort(rows, [:error_rate, :readlen, :arm]))
        println(join(_fmt.([r.error_rate, r.regime, r.readlen, r.arm, r.n_contigs,
            r.total_length, r.largest_contig, r.n50, r.genome_fraction, r.runtime_s]), ""))
    end

    # === Scale guard ===
    eff_cov_for_guard = isfinite(min_effective_coverage) ? min_effective_coverage : 0.0
    metric = scale_metric_bases(eff_cov_for_guard, glen)
    verdict_allowed = scale_verdict_allowed(eff_cov_for_guard, glen; floor = scale_floor)
    mode = verdict_allowed ? "VERDICT" : "SMOKE-ONLY"

    # === Machine-readable CSV ===
    rows[!, :mode] = fill(mode, DataFrames.nrow(rows))
    rows[!, :scale_metric_bases] = fill(metric, DataFrames.nrow(rows))
    rows[!, :scale_floor_bases] = fill(scale_floor, DataFrames.nrow(rows))
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    csv_path = joinpath(results_dir, "rhizomorph_correction_validation_sweep_$(timestamp).csv")
    CSV.write(csv_path, rows)
    println("\nCSV written: $csv_path")

    if verdict_allowed
        println("\n=== VERDICT (scale metric $(round(metric; digits=0)) bases >= floor $(scale_floor)) ===")
        println("Sweep meets the scale floor; the naive-vs-iterative comparison above is a")
        println("validation-grade result. Compare `iterative` against `naive` per cell:")
        println("  - genome_fraction / n50 higher for iterative => correction core helps at that (err, regime)")
        println("  - lower or equal => correction core neutral or harmful at that (err, regime)")
    else
        println("\n=== SMOKE-ONLY (scale metric $(round(metric; digits=0)) bases < floor $(scale_floor)) ===")
        println("This run is BELOW the scale floor and MUST NOT be quoted as validation.")
        println("It confirms the harness parses and runs end-to-end only. Raise coverage,")
        println("genome size, or read count (or run the full Lambda-phage config on Lovelace)")
        println("so effective coverage x genome length exceeds $(scale_floor) bases for a VERDICT.")
    end

    println("\n=== Sweep complete: $(Dates.now()) ===")
    return (csv_path = csv_path, mode = mode, rows = rows)
end

# Only run the sweep when invoked as a script; `include`-ing this file (e.g. from
# the unit test) just defines the functions.
if abspath(PROGRAM_FILE) == @__FILE__
    run_sweep()
end
