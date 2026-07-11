# Rhizomorph Correction CALIBRATION Benchmark (Tier-1)
# ====================================================
#
# "Reference-free correctness prediction": can a signal available WITHOUT a
# reference predict whether a corrected read is actually right? This is the
# infrastructure for td-pw7p Priority #2 (the flagship differentiator that
# justifies the Nature-Methods ambition).
#
# TIER-1 (this harness): calibrate the reference-free MINIMUM K-MER COVERAGE of a
# corrected read — computed from a k-mer count table over the corrected read set
# ALONE (no reference, no decoder internals). A residual error spawns
# low-coverage k-mers, so min-coverage is the canonical "solid vs weak" signal
# Stage 0 already trusts. We measure how well it DISCRIMINATES correct from
# incorrect reads (AUROC) and how CALIBRATED a naive monotone map of it is
# (reliability / ECE / Brier).
#
# TIER-2 (follow-on, needs decoder instrumentation): the per-BASE best-vs-2nd-best
# path log-probability GAP. Both Viterbi loops (viterbi-next.jl, path-finding.jl)
# already visit every competing successor but keep only the argmax; capturing the
# runner-up is a localized but hot-loop change with a perf-neutral requirement, so
# it is deliberately deferred. This Tier-1 harness establishes the calibration
# toolkit (calibration_metrics.jl) that Tier-2 will reuse verbatim on the gap
# signal; only the confidence-extraction changes.
#
# IMPORTANT (why not the corrected read's Phred): adjust_quality_scores
# (src/iterative-assembly.jl) is a NO-OP on same-length reads — it returns the
# input quality unchanged — so the corrected FASTQ's Phred carries zero
# correction-confidence and is NOT a usable signal. Min-k-mer-coverage is.
#
# Signal -> probability: AUROC is rank-based and uses the raw min-coverage
# (discrimination, no mapping needed) — it is the HONEST headline. ECE / Brier /
# reliability need a probability, so we apply a naive monotone map
# p = min_cov / (min_cov + C0). A PROPERLY FITTED map (Platt / isotonic on a
# held-out split) is a follow-on; the naive map is illustrative and labelled as
# such.
#
# Usage:
#   MYCELIA_RCC_SMOKE=true LD_LIBRARY_PATH='' \
#     julia --project=. benchmarking/rhizomorph_correction_calibration_benchmark.jl
#
# Env (RCC prefix; error/coverage/k mirror the accuracy harness):
#   MYCELIA_RCC_SMOKE, MYCELIA_RCC_SMOKE_GENOME_LEN, MYCELIA_RCC_ACCESSION,
#   MYCELIA_RCC_ERR, MYCELIA_RCC_READLEN, MYCELIA_RCC_COVERAGE, MYCELIA_RCC_K,
#   MYCELIA_RCC_ASSIGNED_Q, MYCELIA_RCC_SEED,
#   MYCELIA_RCC_KCAL   (k for the coverage signal, default 15),
#   MYCELIA_RCC_C0     (naive map denominator constant, default 2.0),
#   MYCELIA_RCC_NBINS  (reliability bins, default 10)

import Pkg
if isinteractive()
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# Reuse the accuracy harness (simulate_substitution_reads, correct_reads_scalable,
# acquire_reference) + the pure calibration metrics. Including the accuracy
# benchmark also pulls the validation-sweep helpers + scale guard it includes.
include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_benchmark.jl"))
include(joinpath(@__DIR__, "calibration_metrics.jl"))

# === Reference-free min-k-mer-coverage signal ===============================

"""Canonical (strand-merged) k-mer string: min of the k-mer and its rev-comp."""
function _canonical_kmer(s::AbstractString)
    comp = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'N' => 'N')
    rc = String([get(comp, c, c) for c in reverse(collect(s))])
    return s <= rc ? String(s) : rc
end

"""
Count canonical k_cal-mers across ALL corrected reads (a reference-free coverage
table). Returns a Dict{String,Int}.
"""
function corrected_kmer_coverage(corrected_by_id::Dict{String, String}, k_cal::Int)
    counts = Dict{String, Int}()
    for (_id, seq) in corrected_by_id
        length(seq) < k_cal && continue
        cs = collect(seq)
        for i in 1:(length(cs) - k_cal + 1)
            kmer = _canonical_kmer(String(cs[i:(i + k_cal - 1)]))
            counts[kmer] = get(counts, kmer, 0) + 1
        end
    end
    return counts
end

"""Minimum canonical-k_cal-mer coverage along `seq` (Inf-proxy -1 if too short)."""
function min_kmer_coverage(seq::AbstractString, counts::Dict{String, Int}, k_cal::Int)
    length(seq) < k_cal && return -1
    cs = collect(seq)
    m = typemax(Int)
    for i in 1:(length(cs) - k_cal + 1)
        kmer = _canonical_kmer(String(cs[i:(i + k_cal - 1)]))
        m = min(m, get(counts, kmer, 0))
    end
    return m
end

# === Main =====================================================================

function run_calibration_benchmark()
    smoke = _truthy(get(ENV, "MYCELIA_RCC_SMOKE", "false"))
    accession = get(ENV, "MYCELIA_RCC_ACCESSION", "NC_001416")
    smoke_len = parse(Int, get(ENV, "MYCELIA_RCC_SMOKE_GENOME_LEN", "2000"))
    errs = _parse_float_list(get(ENV, "MYCELIA_RCC_ERR", ""), [0.02, 0.05, 0.10])
    readlens = _parse_int_list(get(ENV, "MYCELIA_RCC_READLEN", ""), [150])
    coverage = parse(Float64, get(ENV, "MYCELIA_RCC_COVERAGE", smoke ? "15" : "30"))
    k = parse(Int, get(ENV, "MYCELIA_RCC_K", "21"))
    assigned_q = parse(Int, get(ENV, "MYCELIA_RCC_ASSIGNED_Q", "20"))
    seed = parse(Int, get(ENV, "MYCELIA_RCC_SEED", "42"))
    k_cal = parse(Int, get(ENV, "MYCELIA_RCC_KCAL", "15"))
    c0 = parse(Float64, get(ENV, "MYCELIA_RCC_C0", "2.0"))
    nbins = parse(Int, get(ENV, "MYCELIA_RCC_NBINS", "10"))

    println("=== Rhizomorph Correction CALIBRATION Benchmark (Tier-1: min-k-mer-coverage) ===")
    println("Start time : $(Dates.now())")
    println("Mode       : $(smoke ? "SMOKE (synthetic, no network)" : "FULL (download $accession)")")
    println("Error rates: $errs   coverage: $(coverage)x   k: $k   k_cal: $k_cal   C0: $c0")

    workdir = mktempdir(prefix = "rcc_bench_")
    refseq, _ref_path,
    ref_label = acquire_reference(
        smoke = smoke, accession = accession, smoke_len = smoke_len, seed = seed, workdir = workdir)
    glen = length(refseq)
    println("Reference  : $ref_label ($glen bp)")
    if k > min(minimum(readlens), glen)
        k = min(minimum(readlens), glen)
    end
    rng = Random.MersenneTwister(seed)

    # Pool per-read (confidence, correct) across all error-rate cells so both
    # classes (correct / incorrect reads) are represented.
    min_covs = Int[]
    corrects = Bool[]
    n_dropped = 0
    n_short = 0

    for err in errs
        for readlen in readlens
            println("\n[cell] err=$err readlen=$readlen")
            records, truth_by_id,
            _obs,
            injected_total,
            _sb = simulate_substitution_reads(
                refseq, readlen, coverage, err, rng; assigned_q = assigned_q)
            println("  $(length(records)) reads, injected $injected_total substitutions")
            local corrected_by_id
            try
                corrected_by_id, _r = correct_reads_scalable(records, k)
            catch e
                @warn "corrector failed; skipping cell" err readlen exception = (
                    e, catch_backtrace())
                continue
            end
            counts = corrected_kmer_coverage(corrected_by_id, k_cal)
            for (rid, T) in truth_by_id
                if !haskey(corrected_by_id, rid)
                    n_dropped += 1
                    continue
                end
                C = corrected_by_id[rid]
                mc = min_kmer_coverage(C, counts, k_cal)
                if mc < 0
                    n_short += 1
                    continue
                end
                push!(min_covs, mc)
                push!(corrects, C == T)   # reference-free signal vs TRUE correctness
            end
        end
    end

    n = length(min_covs)
    println("\n--- Pooled calibration set ---")
    println("scored reads: $n   dropped: $n_dropped   too-short(<k_cal): $n_short")
    if n == 0
        println("No scored reads — cannot calibrate.");
        return nothing
    end
    n_correct = count(corrects)
    println("fully-correct reads: $n_correct  ($(round(100n_correct/n; digits=1))%)   incorrect: $(n - n_correct)")

    # Rank-based discrimination on the RAW signal (no probability mapping needed).
    au = auroc(min_covs, corrects)
    println("\nAUROC (min-k-mer-coverage discriminates correct vs incorrect reads): $(round(au; digits=4))")
    if isnan(au)
        println("  (AUROC undefined — only one correctness class present; vary error rate/coverage)")
    end

    # Naive monotone map to a probability for reliability/ECE/Brier (a fitted
    # Platt/isotonic map is the follow-on; this is illustrative).
    probs = Float64[mc / (mc + c0) for mc in min_covs]
    ece = expected_calibration_error(probs, corrects; nbins = nbins)
    br = brier_score(probs, corrects)
    println("\nNaive map p = min_cov/(min_cov + $c0):  ECE=$(round(ece; digits=4))  Brier=$(round(br; digits=4))")
    println("Reliability (naive map):")
    println("  ", rpad("bin", 12), rpad("mean_conf", 12), rpad("accuracy", 12), "count")
    for b in reliability_bins(probs, corrects; nbins = nbins)
        println("  ", rpad("[$(round(b.lo;digits=2)),$(round(b.hi;digits=2)))", 12),
            rpad(string(round(b.mean_conf; digits = 4)), 12),
            rpad(string(round(b.accuracy; digits = 4)), 12), b.count)
    end

    println("\nInterpretation:")
    println("  AUROC >> 0.5 => min-k-mer-coverage is a usable reference-free correctness signal.")
    println("  AUROC ~ 0.5  => this signal does not discriminate here (try Tier-2 log-prob gap).")
    println("  ECE near 0   => the naive map is already well-calibrated; large => needs a fitted map.")
    println("NOTE: Tier-2 (per-base best-vs-2nd-best path log-prob gap) reuses calibration_metrics.jl")
    println("      verbatim; only the confidence-extraction changes. That is the flagship signal.")
    println("\n=== Calibration benchmark complete: $(Dates.now()) ===")
    return (; n, auroc = au, ece, brier = br)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_calibration_benchmark()
end
