# ==============================================================================
# K-mer classification (Stage 0) — à la carte comparison ecosystem
#
# Distinguishes real-genomic k-mers from sequencing-error k-mers so that clean
# reads/positions can be skipped before expensive per-read correction. Unlike the
# speed-first heuristics in most assemblers (which discard low-coverage k-mers and
# throw out real biological variation alongside error), this is an EXTENSIBLE
# family of strategies benchmarked head-to-head, with a probabilistic,
# quality-aware default that separates error from real low-frequency variation
# with statistical grounding.
#
# Each strategy is a `KmerClassifier` subtype implementing `classify_kmers`. All
# return a `KmerClassification`: a per-k-mer boolean verdict PLUS a continuous
# score (higher = more likely real), so strategies can be compared by ROC/PR
# against a known-truth error set — the comparison matrix that IS the manuscript
# evidence and the input→workflow router's training data.
#
# Coverage is read from `count_total_observations(vertex_data)`; quality from
# `get_vertex_joint_quality(vertex_data, dataset_id)` (per-position combined Phred,
# where combining = adding Phred in log space = accumulating observation
# log-likelihoods). The principled default fuses BOTH.
# ==============================================================================

"""
    KmerClassifier

Abstract supertype for k-mer classification strategies (Stage 0). Mirrors the
`PathFindingStrategy` dispatch pattern: each concrete strategy implements
[`classify_kmers`](@ref) and is one arm of the comparison matrix.
"""
abstract type KmerClassifier end

"""
    KmerClassification{K}

Result of classifying every k-mer (vertex label type `K`) in a graph.

- `verdicts`  : `label => true` (solid / real-genomic) or `false` (weak / error)
- `scores`    : `label => Float64` continuous confidence (higher ⇒ more likely
                real); used for threshold-independent ROC/PR comparison
- `strategy`  : human-readable strategy name (for the comparison matrix / paper)
- `params`    : the strategy's resolved parameters (e.g. the chosen threshold),
                so an auto-selected cutoff is recorded, not hidden
"""
struct KmerClassification{K}
    verdicts::Dict{K, Bool}
    scores::Dict{K, Float64}
    strategy::String
    params::Dict{Symbol, Any}
end

"""
    classify_kmers(classifier::KmerClassifier, graph; dataset_id) -> KmerClassification

Classify every vertex (k-mer) of a qualmer `graph` as solid (real-genomic) or weak
(sequencing error). `dataset_id` selects which observation set's quality/coverage
evidence to use. Each `KmerClassifier` subtype provides its own method.
"""
function classify_kmers end

# ------------------------------------------------------------------------------
# Shared evidence extraction — every strategy reads (coverage, joint-quality)
# through this one accessor so the arms differ only in their DECISION RULE, not
# in how they touch the graph.
# ------------------------------------------------------------------------------

"""
    _kmer_evidence(graph, label, dataset_id) -> (coverage::Int, joint_quality::Vector{UInt8})

Pull the two evidence signals for one k-mer: `coverage` = number of observations
supporting it; `joint_quality` = per-position combined Phred across those
observations (empty if none for this dataset).
"""
function _kmer_evidence(graph, label, dataset_id::AbstractString)
    vertex_data = graph[label]
    coverage = count_total_observations(vertex_data)
    joint_quality = get_vertex_joint_quality(vertex_data, dataset_id)
    return coverage, joint_quality
end

"""
    _mean_phred(joint_quality) -> Float64

Mean combined-Phred across the k-mer's positions (0.0 if no quality evidence).
A convenience summary of the per-position joint-quality vector.
"""
function _mean_phred(joint_quality::AbstractVector{UInt8})::Float64
    isempty(joint_quality) && return 0.0
    return Statistics.mean(Float64.(joint_quality))
end

# ------------------------------------------------------------------------------
# Arm 1 (comparison baseline) — FixedCoverageThreshold
# The classic heuristic: a k-mer is solid iff its multiplicity ≥ a fixed cutoff.
# Coverage ONLY; ignores quality. This is what most assemblers do — included so
# the matrix can quantify exactly what fusing quality buys over it.
# ------------------------------------------------------------------------------

"""
    FixedCoverageThreshold(min_count::Int = 2)

Solid iff observation count ≥ `min_count`. Score = raw coverage. Coverage-only
heuristic baseline (the "bag of tricks" reference point).
"""
struct FixedCoverageThreshold <: KmerClassifier
    min_count::Int
end
FixedCoverageThreshold() = FixedCoverageThreshold(2)

function classify_kmers(clf::FixedCoverageThreshold, graph; dataset_id::AbstractString)
    labels = collect(MetaGraphsNext.labels(graph))
    K = eltype(labels)
    verdicts = Dict{K, Bool}()
    scores = Dict{K, Float64}()
    for label in labels
        coverage, _ = _kmer_evidence(graph, label, dataset_id)
        verdicts[label] = coverage >= clf.min_count
        scores[label] = Float64(coverage)
    end
    return KmerClassification{K}(
        verdicts, scores, "FixedCoverageThreshold",
        Dict{Symbol, Any}(:min_count => clf.min_count),
    )
end

# ------------------------------------------------------------------------------
# Arm 2 (comparison baseline) — QualityThreshold
# The mirror image: a k-mer is solid iff its mean combined-Phred ≥ a cutoff.
# Quality ONLY; ignores coverage. Isolating the quality signal lets the matrix
# separate "quality helps" from "coverage helps" before the fused model claims
# credit for either.
# ------------------------------------------------------------------------------

"""
    QualityThreshold(min_mean_phred::Float64 = 20.0)

Solid iff mean combined-Phred ≥ `min_mean_phred`. Score = mean combined-Phred.
Quality-only baseline (contrast against coverage-only).
"""
struct QualityThreshold <: KmerClassifier
    min_mean_phred::Float64
end
QualityThreshold() = QualityThreshold(20.0)

function classify_kmers(clf::QualityThreshold, graph; dataset_id::AbstractString)
    labels = collect(MetaGraphsNext.labels(graph))
    K = eltype(labels)
    verdicts = Dict{K, Bool}()
    scores = Dict{K, Float64}()
    for label in labels
        _, joint_quality = _kmer_evidence(graph, label, dataset_id)
        mean_phred = _mean_phred(joint_quality)
        verdicts[label] = mean_phred >= clf.min_mean_phred
        scores[label] = mean_phred
    end
    return KmerClassification{K}(
        verdicts, scores, "QualityThreshold",
        Dict{Symbol, Any}(:min_mean_phred => clf.min_mean_phred),
    )
end

# ------------------------------------------------------------------------------
# Shared machinery for the FUSED arms (coverage + quality together).
# The user's directive: implement all three fusion rules and let the comparison
# matrix rank them honestly against ground truth — do not pick a priori.
# A `FusedClassifier` differs from its siblings ONLY in `_posterior`; the graph
# walk, verdict, and result packaging are shared.
# ------------------------------------------------------------------------------

"""
    FusedClassifier <: KmerClassifier

Supertype for strategies that fuse coverage AND quality into a per-k-mer
posterior. Concrete subtypes implement `_posterior(clf, coverage, mean_phred)`,
`_strategy_name(clf)`, and `_params(clf)`; all share one `classify_kmers` method.
"""
abstract type FusedClassifier <: KmerClassifier end

_threshold(clf::FusedClassifier)::Float64 = clf.decision_threshold

"""
    _quality_correct_prob(mean_phred) -> Float64

P(the k-mer's bases are correct) from mean combined-Phred: `1 - 10^(-Q/10)`.
Strong accumulated Phred ⇒ ≈1. The quality signal every fused arm consumes.
"""
function _quality_correct_prob(mean_phred::Float64)::Float64
    return 1.0 - 10.0^(-mean_phred / 10.0)
end

function classify_kmers(clf::FusedClassifier, graph; dataset_id::AbstractString)
    labels = collect(MetaGraphsNext.labels(graph))
    K = eltype(labels)
    verdicts = Dict{K, Bool}()
    scores = Dict{K, Float64}()
    thr = _threshold(clf)
    for label in labels
        coverage, joint_quality = _kmer_evidence(graph, label, dataset_id)
        posterior = _posterior(clf, coverage, _mean_phred(joint_quality))
        verdicts[label] = posterior >= thr
        scores[label] = posterior
    end
    return KmerClassification{K}(verdicts, scores, _strategy_name(clf), _params(clf))
end

# ------------------------------------------------------------------------------
# Fusion A — BayesianMixtureClassifier (quality as INDEPENDENT evidence)
# P(real) ∝ prior · Poisson(cov; λ_genomic) · q ;  q = P(bases correct).
# A clean single-coverage k-mer can still be called real (true rare variant).
# ------------------------------------------------------------------------------

"""
    BayesianMixtureClassifier(; genomic_mean_coverage, error_mean_coverage = 1.0,
                                prior_real = 0.5, decision_threshold = 0.5)

Fuse coverage and quality as independent evidence via a Poisson coverage mixture
weighted by the quality-correct probability. Quality votes separately from
coverage, so strong quality can rescue a low-coverage real k-mer.
"""
struct BayesianMixtureClassifier <: FusedClassifier
    genomic_mean_coverage::Float64
    error_mean_coverage::Float64
    prior_real::Float64
    decision_threshold::Float64
end
function BayesianMixtureClassifier(; genomic_mean_coverage::Float64,
                                     error_mean_coverage::Float64 = 1.0,
                                     prior_real::Float64 = 0.5,
                                     decision_threshold::Float64 = 0.5)
    return BayesianMixtureClassifier(genomic_mean_coverage, error_mean_coverage,
                                     prior_real, decision_threshold)
end

function _posterior(clf::BayesianMixtureClassifier, coverage::Int, mean_phred::Float64)::Float64
    q = _quality_correct_prob(mean_phred)
    l_real = Distributions.pdf(Distributions.Poisson(clf.genomic_mean_coverage), coverage) * q
    l_error = Distributions.pdf(Distributions.Poisson(clf.error_mean_coverage), coverage) * (1.0 - q)
    num = clf.prior_real * l_real
    den = num + (1.0 - clf.prior_real) * l_error
    return den > 0.0 ? num / den : clf.prior_real
end
_strategy_name(::BayesianMixtureClassifier) = "BayesianMixture"
_params(clf::BayesianMixtureClassifier) = Dict{Symbol, Any}(
    :genomic_mean_coverage => clf.genomic_mean_coverage,
    :error_mean_coverage => clf.error_mean_coverage,
    :prior_real => clf.prior_real,
    :decision_threshold => clf.decision_threshold,
)

# ------------------------------------------------------------------------------
# Fusion B — EffectiveCoverageClassifier (quality DISCOUNTS coverage)
# eff_cov = cov · q, then a Poisson likelihood ratio on that single axis.
# Fewer independence assumptions; quality shrinks coverage rather than voting
# separately (closest to Quake's q-mer idea).
# ------------------------------------------------------------------------------

"""
    EffectiveCoverageClassifier(; genomic_mean_coverage, error_mean_coverage = 1.0,
                                  prior_real = 0.5, decision_threshold = 0.5)

Quality-adjust coverage (`eff_cov = coverage · P(correct)`), then a Poisson
likelihood ratio between genomic and error components on the effective coverage.
"""
struct EffectiveCoverageClassifier <: FusedClassifier
    genomic_mean_coverage::Float64
    error_mean_coverage::Float64
    prior_real::Float64
    decision_threshold::Float64
end
function EffectiveCoverageClassifier(; genomic_mean_coverage::Float64,
                                       error_mean_coverage::Float64 = 1.0,
                                       prior_real::Float64 = 0.5,
                                       decision_threshold::Float64 = 0.5)
    return EffectiveCoverageClassifier(genomic_mean_coverage, error_mean_coverage,
                                       prior_real, decision_threshold)
end

function _posterior(clf::EffectiveCoverageClassifier, coverage::Int, mean_phred::Float64)::Float64
    # Quality discounts coverage, then a Poisson likelihood ratio on the single
    # effective-coverage axis (rounded to an integer count for the discrete model).
    eff_cov = round(Int, coverage * _quality_correct_prob(mean_phred))
    l_real = clf.prior_real * Distributions.pdf(Distributions.Poisson(clf.genomic_mean_coverage), eff_cov)
    l_error = (1.0 - clf.prior_real) * Distributions.pdf(Distributions.Poisson(clf.error_mean_coverage), eff_cov)
    den = l_real + l_error
    return den > 0.0 ? l_real / den : clf.prior_real
end
_strategy_name(::EffectiveCoverageClassifier) = "EffectiveCoverage"
_params(clf::EffectiveCoverageClassifier) = Dict{Symbol, Any}(
    :genomic_mean_coverage => clf.genomic_mean_coverage,
    :error_mean_coverage => clf.error_mean_coverage,
    :prior_real => clf.prior_real,
    :decision_threshold => clf.decision_threshold,
)

# ------------------------------------------------------------------------------
# Fusion C — LogisticFusionClassifier (LEARNED boundary)
# post = σ(w0 + w1·cov + w2·mean_phred), weights fit from a labeled k-mer set
# (ground truth from simulated reads). Most flexible; least interpretable.
# ------------------------------------------------------------------------------

"""
    LogisticFusionClassifier(w0, w_cov, w_phred; decision_threshold = 0.5)

Learned logistic fusion of coverage and mean-Phred. Fit weights with
[`fit_logistic_fusion`](@ref) from a labeled (coverage, mean_phred) → {real,error}
training set, then classify.
"""
struct LogisticFusionClassifier <: FusedClassifier
    w0::Float64
    w_cov::Float64
    w_phred::Float64
    decision_threshold::Float64
end
LogisticFusionClassifier(w0, w_cov, w_phred; decision_threshold::Float64 = 0.5) =
    LogisticFusionClassifier(w0, w_cov, w_phred, decision_threshold)

function _posterior(clf::LogisticFusionClassifier, coverage::Int, mean_phred::Float64)::Float64
    z = clf.w0 + clf.w_cov * coverage + clf.w_phred * mean_phred
    return 1.0 / (1.0 + exp(-z))
end
_strategy_name(::LogisticFusionClassifier) = "LogisticFusion"
_params(clf::LogisticFusionClassifier) = Dict{Symbol, Any}(
    :w0 => clf.w0, :w_cov => clf.w_cov, :w_phred => clf.w_phred,
    :decision_threshold => clf.decision_threshold,
)

"""
    fit_logistic_fusion(coverages, mean_phreds, labels; iters = 500, lr = 0.05,
                        decision_threshold = 0.5) -> LogisticFusionClassifier

Fit logistic-fusion weights by gradient descent on standardized features.
`labels[i]` is `true` for a real-genomic k-mer, `false` for error. Features are
standardized internally and the learned weights are folded back to raw
(coverage, mean_phred) space so the returned classifier applies directly.
"""
function fit_logistic_fusion(coverages::AbstractVector{<:Real},
                             mean_phreds::AbstractVector{<:Real},
                             labels::AbstractVector{Bool};
                             iters::Int = 500, lr::Float64 = 0.05,
                             decision_threshold::Float64 = 0.5)
    n = length(labels)
    @assert n == length(coverages) == length(mean_phreds) "feature/label length mismatch"
    x1 = Float64.(coverages); x2 = Float64.(mean_phreds); y = Float64.(labels)
    m1, s1 = Statistics.mean(x1), max(Statistics.std(x1), 1e-9)
    m2, s2 = Statistics.mean(x2), max(Statistics.std(x2), 1e-9)
    z1 = (x1 .- m1) ./ s1; z2 = (x2 .- m2) ./ s2   # standardized features
    b0, b1, b2 = 0.0, 0.0, 0.0
    for _ in 1:iters
        p = 1.0 ./ (1.0 .+ exp.(-(b0 .+ b1 .* z1 .+ b2 .* z2)))
        g = p .- y
        b0 -= lr * Statistics.mean(g)
        b1 -= lr * Statistics.mean(g .* z1)
        b2 -= lr * Statistics.mean(g .* z2)
    end
    # fold standardization back into raw-space weights
    w_cov = b1 / s1
    w_phred = b2 / s2
    w0 = b0 - b1 * m1 / s1 - b2 * m2 / s2
    return LogisticFusionClassifier(w0, w_cov, w_phred; decision_threshold = decision_threshold)
end
