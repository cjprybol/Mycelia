# ==============================================================================
# Stage 0 k-mer classification — score-threshold calibration
#
# The comparison matrix (stage0_classification_matrix.jl) surfaces a key finding:
# threshold-independent AUC is nearly identical across arms (~0.95), yet their
# DEFAULT decision thresholds are badly miscalibrated — QualityThreshold's fixed
# Q20 cutoff yields specificity 0.0 (every k-mer called solid), and
# LogisticFusion's 0.5 posterior cutoff sits far off the class-imbalanced optimum
# (balanced accuracy ≈ 0.76). Since AUC says the SCORES separate the classes well,
# the loss is purely in where each arm draws its line.
#
# This script quantifies that gap. For one base config it builds the same
# simulated truth + per-k-mer scores as the matrix, then for EACH arm sweeps its
# continuous score over the observed range, finds the threshold that MAXIMIZES
# balanced accuracy, and reports the recovered balanced accuracy against the
# arm's default. The delta is exactly what calibration buys — with the AUC held
# fixed, it isolates "better cutoff" from "better score".
#
# Config (env-overridable, matched to the matrix harness):
#   MYCELIA_S0_REFLEN   reference length            (default 10000)
#   MYCELIA_S0_COVERAGE  mean coverage               (default 30)
#   MYCELIA_S0_READLEN   read length                 (default 100)
#   MYCELIA_S0_ERRORRATE substitution error rate     (default 0.01)
#   MYCELIA_S0_K         k-mer length                (default 21)
#   MYCELIA_S0_SEED      RNG seed                    (default 42)
# ==============================================================================

import Mycelia
import FASTX
import BioSequences
import Kmers
import MetaGraphsNext
import Random
import Statistics
import JSON

const R = Mycelia.Rhizomorph

getenv(k, default) = haskey(ENV, k) ? parse(typeof(default), ENV[k]) : default

const REFLEN     = getenv("MYCELIA_S0_REFLEN", 10000)
const COVERAGE   = getenv("MYCELIA_S0_COVERAGE", 30)
const READLEN    = getenv("MYCELIA_S0_READLEN", 100)
const ERROR_RATE = getenv("MYCELIA_S0_ERRORRATE", 0.01)
const K          = getenv("MYCELIA_S0_K", 21)
const SEED       = getenv("MYCELIA_S0_SEED", 42)
const DATASET    = "s0"

const BASES = ['A', 'C', 'G', 'T']

# ------------------------------------------------------------------------------
# Simulate + truth — identical model to stage0_classification_matrix.jl so the
# calibration is measured on the SAME data the matrix scores against.
# ------------------------------------------------------------------------------

"""Phred quality model: correct bases ~Q35, error bases ~Q22 (overlapping), so
quality carries partial signal. Matches the matrix harness."""
function sample_quality(rng, is_error::Bool)::Char
    q = is_error ? round(Int, 22 + 8 * randn(rng)) : round(Int, 35 + 3 * randn(rng))
    q = clamp(q, 2, 40)
    return Char(q + 33)
end

function simulate(rng)
    ref = join(rand(rng, BASES, REFLEN))
    n_reads = max(1, round(Int, COVERAGE * REFLEN / READLEN))
    recs = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        start = rand(rng, 1:(REFLEN - READLEN + 1))
        seq = collect(ref[start:(start + READLEN - 1)])
        quals = Vector{Char}(undef, READLEN)
        for j in 1:READLEN
            is_err = rand(rng) < ERROR_RATE
            if is_err
                alt = rand(rng, filter(!=(seq[j]), BASES))
                seq[j] = alt
            end
            quals[j] = sample_quality(rng, is_err)
        end
        push!(recs, FASTX.FASTQ.Record("r$i", join(seq), join(quals)))
    end
    return ref, recs
end

"""Truth set: canonical k-mers of the error-free reference = genomic k-mers."""
function reference_kmer_set(ref::AbstractString)
    truth = Set{Kmers.DNAKmer{K}}()
    seq = BioSequences.LongDNA{4}(ref)
    for (kmer, _) in Kmers.UnambiguousDNAMers{K}(seq)
        push!(truth, BioSequences.canonical(kmer))
    end
    return truth
end

# ------------------------------------------------------------------------------
# Threshold sweep — the calibration core.
# ------------------------------------------------------------------------------

"""
    balanced_accuracy_at(scores, ys, threshold) -> Float64

Balanced accuracy `(sensitivity + specificity) / 2` when a k-mer is called solid
iff `score >= threshold`. `ys[i]` is `true` for a genomic k-mer. NaN if either
class is empty.
"""
function balanced_accuracy_at(scores::AbstractVector{Float64},
                              ys::AbstractVector{Bool},
                              threshold::Float64)::Float64
    tp = fp = tn = fn = 0
    for (s, y) in zip(scores, ys)
        called = s >= threshold
        if y
            called ? (tp += 1) : (fn += 1)
        else
            called ? (fp += 1) : (tn += 1)
        end
    end
    sens = tp + fn > 0 ? tp / (tp + fn) : NaN
    spec = tn + fp > 0 ? tn / (tn + fp) : NaN
    return (sens + spec) / 2
end

"""
    calibrate(scores, ys) -> (best_threshold, best_bal_acc)

Sweep candidate thresholds spanning the observed score range and return the one
that maximizes balanced accuracy. Candidates are the midpoints between adjacent
unique scores (plus an above-max cutoff that rejects everything), so every
distinct partition of the k-mers into solid/weak is evaluated.
"""
function calibrate(scores::AbstractVector{Float64}, ys::AbstractVector{Bool})
    uniq = sort!(unique(scores))
    # midpoints between adjacent unique scores + a below-min and above-max bound,
    # covering "accept all" through "reject all".
    candidates = Float64[]
    push!(candidates, nextfloat(uniq[end]))            # reject all
    for i in 1:(length(uniq) - 1)
        push!(candidates, (uniq[i] + uniq[i + 1]) / 2)  # every class-partition boundary
    end
    push!(candidates, uniq[1])                          # accept all (score >= min)
    best_t = candidates[1]
    best_ba = -Inf
    for t in candidates
        ba = balanced_accuracy_at(scores, ys, t)
        if !isnan(ba) && ba > best_ba
            best_ba = ba
            best_t = t
        end
    end
    return best_t, best_ba
end

"""Default-threshold balanced accuracy straight from the arm's own verdicts."""
function default_balanced_accuracy(classification, truth, labels)::Float64
    tp = fp = tn = fn = 0
    for label in labels
        y = label in truth
        called = classification.verdicts[label]
        if y
            called ? (tp += 1) : (fn += 1)
        else
            called ? (fp += 1) : (tn += 1)
        end
    end
    sens = tp + fn > 0 ? tp / (tp + fn) : NaN
    spec = tn + fp > 0 ? tn / (tn + fp) : NaN
    return (sens + spec) / 2
end

function main()
    rng = Random.MersenneTwister(SEED)
    println("=== Stage 0 threshold calibration ===")
    println("ref=$REFLEN cov=$COVERAGE readlen=$READLEN err=$ERROR_RATE k=$K seed=$SEED")
    ref, recs = simulate(rng)
    truth = reference_kmer_set(ref)
    # :full retains per-observation quality for the depth-independent per-base mean
    # the fused arms need (review F2/C1); :lightweight_quality would silently zero
    # the quality signal and corrupt QualityThreshold / the fused arms.
    graph = R.build_qualmer_graph(recs, K; dataset_id = DATASET, mode = :canonical,
                                  memory_profile = :full)
    labels = collect(MetaGraphsNext.labels(graph))
    ys = Bool[label in truth for label in labels]
    n_genomic = count(ys)
    println("graph k-mers: $(length(labels))  (genomic=$n_genomic  error=$(length(labels) - n_genomic))\n")

    # fit the learned arm on the graph's own (coverage, phred, truth), matching
    # the matrix harness so the LogisticFusion default is the same one it reports.
    covs = Int[]; phreds = Float64[]
    for label in labels
        cov, jq = R._kmer_evidence(graph, label, DATASET)
        push!(covs, cov); push!(phreds, R._mean_phred(jq))
    end
    logistic = R.fit_logistic_fusion(covs, phreds, ys)

    arms = [
        R.FixedCoverageThreshold(2),
        R.QualityThreshold(20.0),
        R.MixtureModelClassifier(),
        R.BloomFilterClassifier(),
        R.BayesianMixtureClassifier(genomic_mean_coverage = Float64(COVERAGE)),
        R.EffectiveCoverageClassifier(genomic_mean_coverage = Float64(COVERAGE)),
        logistic,
    ]

    println(rpad("strategy", 22), lpad("default_ba", 12), lpad("calib_ba", 10),
            lpad("gain", 8), lpad("best_thr", 12))
    println("-"^64)
    results = Dict{String, Any}()
    for clf in arms
        c = R.classify_kmers(clf, graph; dataset_id = DATASET)
        scores = Float64[c.scores[label] for label in labels]
        default_ba = default_balanced_accuracy(c, truth, labels)
        best_t, best_ba = calibrate(scores, ys)
        gain = best_ba - default_ba
        results[c.strategy] = Dict(
            "default_bal_acc" => default_ba,
            "calibrated_bal_acc" => best_ba,
            "gain" => gain,
            "best_threshold" => best_t,
            "score_min" => minimum(scores),
            "score_max" => maximum(scores),
        )
        r3(x) = string(round(x; digits = 3))
        println(rpad(c.strategy, 22), lpad(r3(default_ba), 12), lpad(r3(best_ba), 10),
                lpad(r3(gain), 8), lpad(string(round(best_t; digits = 4)), 12))
    end

    outdir = joinpath(@__DIR__, "results"); mkpath(outdir)
    config = Dict("reflen" => REFLEN, "coverage" => COVERAGE, "readlen" => READLEN,
                  "error_rate" => ERROR_RATE, "k" => K, "seed" => SEED)
    open(joinpath(outdir, "stage0_calibration_seed$(SEED).json"), "w") do io
        JSON.print(io, Dict("config" => config, "calibration" => results), 2)
    end
    println("\npersisted calibration: ",
            joinpath(outdir, "stage0_calibration_seed$(SEED).json"))
    println("=== DONE ===")
end

main()
