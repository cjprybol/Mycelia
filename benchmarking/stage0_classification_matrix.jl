# ==============================================================================
# Stage 0 k-mer classification — ground-truth comparison matrix
#
# Simulates reads from a KNOWN reference at a controlled coverage + error rate, so
# the truth (which k-mers are genomic vs sequencing-error) is exact. Runs every
# KmerClassifier arm and scores it against truth with threshold-independent AUC
# plus sensitivity/specificity/precision. Emits the comparison matrix (the paper's
# results table) AND persists the per-k-mer (coverage, mean_phred, arm-scores,
# truth) table — the input→workflow router's training data (td-26tt).
#
# Toy-scale-first (feedback_prove-completion-at-toy-scale-first): defaults are a
# small reference so this COMPLETES; env vars scale it up.
#   MYCELIA_S0_REFLEN   reference length            (default 2000)
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

const REFLEN     = getenv("MYCELIA_S0_REFLEN", 2000)
const COVERAGE   = getenv("MYCELIA_S0_COVERAGE", 30)
const READLEN    = getenv("MYCELIA_S0_READLEN", 100)
const ERROR_RATE = getenv("MYCELIA_S0_ERRORRATE", 0.01)
const K          = getenv("MYCELIA_S0_K", 21)
const SEED       = getenv("MYCELIA_S0_SEED", 42)
const DATASET    = "s0"

const BASES = ['A', 'C', 'G', 'T']

"""Phred quality model: correct bases ~Q35, error bases ~Q22 (lower but
overlapping), so quality carries PARTIAL signal — not a perfect oracle. This is
what lets the matrix show where fusing quality with coverage actually helps."""
function sample_quality(rng, is_error::Bool)::Char
    q = is_error ? round(Int, 22 + 8 * randn(rng)) : round(Int, 35 + 3 * randn(rng))
    q = clamp(q, 2, 40)
    return Char(q + 33)  # Phred+33
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

"""Truth set: canonical k-mers of the ERROR-FREE reference = genomic k-mers."""
function reference_kmer_set(ref::AbstractString)
    truth = Set{Kmers.DNAKmer{K}}()
    seq = BioSequences.LongDNA{4}(ref)
    for (kmer, _) in Kmers.UnambiguousDNAMers{K}(seq)
        push!(truth, BioSequences.canonical(kmer))
    end
    return truth
end

"""Mann–Whitney AUC = P(score(genomic) > score(error)); 0.5 = chance, 1.0 = perfect."""
function auc(scores_pos::Vector{Float64}, scores_neg::Vector{Float64})::Float64
    (isempty(scores_pos) || isempty(scores_neg)) && return NaN
    wins = 0.0
    for sp in scores_pos, sn in scores_neg
        wins += sp > sn ? 1.0 : (sp == sn ? 0.5 : 0.0)
    end
    return wins / (length(scores_pos) * length(scores_neg))
end

function evaluate(classification, truth, labels)
    tp = fp = tn = fn = 0
    pos_scores = Float64[]; neg_scores = Float64[]
    for label in labels
        is_genomic = label in truth
        called_solid = classification.verdicts[label]
        score = classification.scores[label]
        is_genomic ? push!(pos_scores, score) : push!(neg_scores, score)
        if is_genomic
            called_solid ? (tp += 1) : (fn += 1)
        else
            called_solid ? (fp += 1) : (tn += 1)
        end
    end
    sens = tp + fn > 0 ? tp / (tp + fn) : NaN          # recall of genomic
    spec = tn + fp > 0 ? tn / (tn + fp) : NaN          # rejection of error
    prec = tp + fp > 0 ? tp / (tp + fp) : NaN
    return (; sens, spec, prec, bal_acc = (sens + spec) / 2,
            auc = auc(pos_scores, neg_scores), tp, fp, tn, fn)
end

function main()
    rng = Random.MersenneTwister(SEED)
    println("=== Stage 0 classification matrix ===")
    println("ref=$REFLEN cov=$COVERAGE readlen=$READLEN err=$ERROR_RATE k=$K seed=$SEED")
    ref, recs = simulate(rng)
    truth = reference_kmer_set(ref)
    graph = R.build_qualmer_graph(recs, K; dataset_id = DATASET, mode = :canonical,
                                  memory_profile = :lightweight_quality)
    labels = collect(MetaGraphsNext.labels(graph))
    n_genomic = count(l -> l in truth, labels)
    println("graph k-mers: $(length(labels))  (genomic=$n_genomic  error=$(length(labels) - n_genomic))\n")

    # fit the learned arm on the graph's own (coverage, phred, truth) — same data
    # the matrix scores against; a production router would use held-out data.
    covs = Int[]; phreds = Float64[]; ys = Bool[]
    for label in labels
        cov, jq = R._kmer_evidence(graph, label, DATASET)
        push!(covs, cov); push!(phreds, R._mean_phred(jq)); push!(ys, label in truth)
    end
    logistic = R.fit_logistic_fusion(covs, phreds, ys)

    arms = [
        R.FixedCoverageThreshold(2),
        R.QualityThreshold(20.0),
        R.BayesianMixtureClassifier(genomic_mean_coverage = Float64(COVERAGE)),
        R.EffectiveCoverageClassifier(genomic_mean_coverage = Float64(COVERAGE)),
        logistic,
    ]

    println(rpad("strategy", 22), lpad("sens", 7), lpad("spec", 8), lpad("prec", 8),
            lpad("bal_acc", 9), lpad("AUC", 8))
    println("-"^62)
    results = Dict{String, Any}()
    classifications = Dict{String, Any}()
    for clf in arms
        c = R.classify_kmers(clf, graph; dataset_id = DATASET)
        m = evaluate(c, truth, labels)
        classifications[c.strategy] = c
        results[c.strategy] = Dict("sens" => m.sens, "spec" => m.spec, "prec" => m.prec,
                                   "bal_acc" => m.bal_acc, "auc" => m.auc,
                                   "tp" => m.tp, "fp" => m.fp, "tn" => m.tn, "fn" => m.fn,
                                   "params" => c.params)
        f(x) = isnan(x) ? "  n/a" : string(round(x; digits = 3))
        println(rpad(c.strategy, 22), lpad(f(m.sens), 7), lpad(f(m.spec), 8),
                lpad(f(m.prec), 8), lpad(f(m.bal_acc), 9), lpad(f(m.auc), 8))
    end

    # persist the router's training table: per-kmer coverage, quality, arm scores, truth
    outdir = joinpath(@__DIR__, "results"); mkpath(outdir)
    per_kmer = Vector{Dict{String, Any}}()
    for (i, label) in enumerate(labels)
        row = Dict{String, Any}("coverage" => covs[i], "mean_phred" => phreds[i],
                                "genomic" => ys[i])
        for (name, c) in classifications
            row["score_$name"] = c.scores[label]
        end
        push!(per_kmer, row)
    end
    open(joinpath(outdir, "stage0_matrix_seed$(SEED).json"), "w") do io
        JSON.print(io, Dict("config" => Dict("reflen" => REFLEN, "coverage" => COVERAGE,
                            "readlen" => READLEN, "error_rate" => ERROR_RATE, "k" => K, "seed" => SEED),
                            "matrix" => results, "per_kmer" => per_kmer), 2)
    end
    println("\npersisted: ", joinpath(outdir, "stage0_matrix_seed$(SEED).json"),
            " ($(length(per_kmer)) k-mer rows — router training data)")
    println("=== DONE ===")
end

main()
