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
# `get_vertex_mean_quality(vertex_data, dataset_id)` (per-position MEAN Phred —
# depth-INDEPENDENT, so quality is an evidence axis distinct from coverage rather
# than a proxy for it). The fused arms combine BOTH independent signals.
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
# Shared evidence extraction — every strategy reads (coverage, mean-quality)
# through this one accessor so the arms differ only in their DECISION RULE, not
# in how they touch the graph.
# ------------------------------------------------------------------------------

"""
    _kmer_evidence(graph, label, dataset_id) -> (coverage::Int, mean_quality)

Pull the two evidence signals for one k-mer: `coverage` = number of observations
supporting it; `mean_quality` = per-position MEAN Phred across observations
(`Vector{Float64}`, raw Phred 0–60), or `nothing` when the vertex has no evidence
for this dataset.

The quality signal is the per-position MEAN, which is depth-INDEPENDENT. Using the
accumulated joint quality (`get_vertex_joint_quality`, a sum of Phred over
observations that grows with coverage and saturates at 255) would make quality a
coverage proxy and defeat the independence the fused arms assume (review C1).
"""
function _kmer_evidence(graph, label, dataset_id::AbstractString)
    vertex_data = graph[label]
    coverage = count_total_observations(vertex_data)
    mean_quality = get_vertex_mean_quality(vertex_data, dataset_id)
    return coverage, mean_quality
end

"""
    _mean_phred(quality) -> Float64

Mean per-base Phred across the k-mer's positions (0.0 for `nothing` or empty).
Accepts the `nothing`/`Vector{Float64}` returned by `get_vertex_mean_quality`
(review C2: the accessor returns `nothing`, not an empty vector, on no evidence).
"""
function _mean_phred(quality)::Float64
    (quality === nothing || isempty(quality)) && return 0.0
    return Statistics.mean(Float64.(quality))
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
# The mirror image: a k-mer is solid iff its mean per-base Phred ≥ a cutoff.
# Quality ONLY; ignores coverage. Isolating the quality signal lets the matrix
# separate "quality helps" from "coverage helps" before the fused model claims
# credit for either.
# ------------------------------------------------------------------------------

"""
    QualityThreshold(min_mean_phred::Float64 = 20.0)

Solid iff mean per-base Phred ≥ `min_mean_phred`. Score = mean per-base Phred.
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
        _, mean_quality = _kmer_evidence(graph, label, dataset_id)
        mean_phred = _mean_phred(mean_quality)
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

P(the k-mer's bases are correct) from mean per-base Phred: `1 - 10^(-Q/10)`.
Strong per-base Phred ⇒ ≈1. The quality signal every fused arm consumes.
"""
function _quality_correct_prob(mean_phred::Float64)::Float64
    # Clamp strictly below 1.0 so the (1 - q) error-likelihood term in the fused
    # arms never underflows to exactly 0 and collapses the posterior (review C1).
    return min(1.0 - 10.0^(-mean_phred / 10.0), 1.0 - 1e-6)
end

function classify_kmers(clf::FusedClassifier, graph; dataset_id::AbstractString)
    labels = collect(MetaGraphsNext.labels(graph))
    K = eltype(labels)
    verdicts = Dict{K, Bool}()
    scores = Dict{K, Float64}()
    thr = _threshold(clf)
    for label in labels
        coverage, mean_quality = _kmer_evidence(graph, label, dataset_id)
        posterior = _posterior(clf, coverage, _mean_phred(mean_quality))
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

# ------------------------------------------------------------------------------
# Arm — MixtureModelClassifier (classic auto-threshold; BLESS/Quake/KmerGenie)
# Build the k-mer coverage histogram over ALL vertices, model it as a mixture of a
# low-count error component and a higher-count genomic component, and auto-select
# the solid/weak boundary at the histogram VALLEY (the first local minimum between
# the error peak and the genomic peak). Coverage-only, quality-agnostic. Because
# it needs the whole histogram, it computes the valley once inside
# `classify_kmers` rather than fitting the fused `_posterior(cov, phred)` contract.
# ------------------------------------------------------------------------------

"""
    MixtureModelClassifier(; min_valley = 2, steepness = 1.0)

Classic auto-threshold coverage classifier (BLESS/Quake/KmerGenie lineage). Builds
the k-mer coverage histogram, treats it as an error-plus-genomic mixture, and picks
the solid/weak cutoff at the histogram VALLEY — the first local minimum separating
the error peak (near coverage 1) from the genomic peak. A k-mer is solid iff its
coverage exceeds the valley. Coverage-only; quality is ignored.

- `min_valley` : floor for the auto-selected valley when no clear valley exists
                 (degenerate/monotonic histograms fall back to this cutoff)
- `steepness`  : logistic steepness for the confidence score around the valley

Score is a monotone confidence `σ(steepness · (coverage − valley))`: increasing in
coverage and crossing 0.5 at the chosen valley.
"""
struct MixtureModelClassifier <: KmerClassifier
    min_valley::Int
    steepness::Float64
end
MixtureModelClassifier(; min_valley::Int = 2, steepness::Float64 = 1.0) =
    MixtureModelClassifier(min_valley, steepness)

"""
    _coverage_histogram(coverages) -> Vector{Int}

Histogram of k-mer coverages as counts indexed by coverage: `hist[c]` = number of
k-mers observed exactly `c` times (index 1 ⇒ coverage 1). Length is the maximum
observed coverage; empty input yields an empty histogram.
"""
function _coverage_histogram(coverages::AbstractVector{<:Integer})::Vector{Int}
    isempty(coverages) && return Int[]
    maxcov = maximum(coverages)
    hist = zeros(Int, maxcov)
    for c in coverages
        c >= 1 && (hist[c] += 1)
    end
    return hist
end

"""
    _histogram_valley(hist, min_valley) -> Int

Auto-select the solid/weak coverage boundary as the first local minimum (valley)
of the coverage histogram that lies above the initial error peak. Walks up from the
low-count error mode until the histogram stops decreasing (the valley), returning
that coverage. Falls back to `min_valley` for monotonic or degenerate histograms.
"""
function _histogram_valley(hist::AbstractVector{Int}, min_valley::Int)::Int
    n = length(hist)
    n < 2 && return min_valley
    # Descend from coverage 1 while the histogram keeps falling (the error mode),
    # then the first point where it turns back up is the valley.
    valley = 1
    while valley < n && hist[valley + 1] < hist[valley]
        valley += 1
    end
    # A genuine valley requires a subsequent genomic peak (histogram rises again).
    # If we ran to the end without rising, there is no clear bimodal split.
    if valley >= n || valley == 1
        return min_valley
    end
    return max(valley, min_valley)
end

function classify_kmers(clf::MixtureModelClassifier, graph; dataset_id::AbstractString)
    labels = collect(MetaGraphsNext.labels(graph))
    K = eltype(labels)
    coverages = Dict{K, Int}()
    for label in labels
        coverage, _ = _kmer_evidence(graph, label, dataset_id)
        coverages[label] = coverage
    end
    hist = _coverage_histogram(collect(values(coverages)))
    valley = _histogram_valley(hist, clf.min_valley)
    verdicts = Dict{K, Bool}()
    scores = Dict{K, Float64}()
    for label in labels
        coverage = coverages[label]
        verdicts[label] = coverage > valley
        # Monotone confidence: logistic centred at the valley (0.5 at the cutoff).
        scores[label] = 1.0 / (1.0 + exp(-clf.steepness * (coverage - valley)))
    end
    return KmerClassification{K}(
        verdicts, scores, "MixtureModel",
        Dict{Symbol, Any}(:valley => valley, :min_valley => clf.min_valley,
                          :steepness => clf.steepness),
    )
end

# ------------------------------------------------------------------------------
# Arm — BloomFilterClassifier (memory-efficient tier)
# Insert every k-mer whose coverage ≥ a cutoff into a minimal Bloom filter (a
# BitVector with a few hash functions over the k-mer's hash), then classify solid =
# membership. Coverage is still the decision signal; the Bloom filter is the
# compact set representation the classic memory-lean pipelines use in place of an
# exact solid-k-mer set. Score = raw coverage. No external Bloom dependency exists,
# so a small self-contained one lives here.
# ------------------------------------------------------------------------------

"""
    _BloomFilter(nbits, nhashes)

Minimal Bloom filter: a `BitVector` of `nbits` slots probed by `nhashes` hash
functions derived from an item's `hash`. Approximate-membership only — no false
negatives, tunable false-positive rate via `nbits`/`nhashes`.
"""
struct _BloomFilter
    bits::BitVector
    nhashes::Int
end
_BloomFilter(nbits::Int, nhashes::Int) = _BloomFilter(falses(nbits), nhashes)

"""
    _bloom_indices(bf, item) -> Vector{Int}

The `nhashes` bit positions an `item` maps to, using a seeded `hash` per function
folded into `[1, length(bf.bits)]`.
"""
function _bloom_indices(bf::_BloomFilter, item)::Vector{Int}
    nbits = length(bf.bits)
    return [Int(mod(hash(item, UInt(i)), nbits)) + 1 for i in 1:bf.nhashes]
end

"""
    _bloom_insert!(bf, item)

Set every bit position `item` maps to. Idempotent.
"""
function _bloom_insert!(bf::_BloomFilter, item)
    for idx in _bloom_indices(bf, item)
        bf.bits[idx] = true
    end
    return bf
end

"""
    _bloom_contains(bf, item) -> Bool

`true` if every bit position `item` maps to is set (possible member; may be a false
positive). `false` is definitive (never a member).
"""
function _bloom_contains(bf::_BloomFilter, item)::Bool
    return all(bf.bits[idx] for idx in _bloom_indices(bf, item))
end

"""
    BloomFilterClassifier(; min_count = 2, nbits = 1 << 16, nhashes = 4)

Memory-efficient solid-k-mer tier. Inserts every k-mer with coverage ≥ `min_count`
into a Bloom filter (`nbits` slots, `nhashes` hash functions), then classifies a
k-mer solid iff the filter reports membership. Coverage-only decision signal;
score = raw coverage. The Bloom filter trades a small false-positive rate for a
compact set representation in place of an exact solid-k-mer table.
"""
struct BloomFilterClassifier <: KmerClassifier
    min_count::Int
    nbits::Int
    nhashes::Int
end
BloomFilterClassifier(; min_count::Int = 2, nbits::Int = 1 << 16, nhashes::Int = 4) =
    BloomFilterClassifier(min_count, nbits, nhashes)

function classify_kmers(clf::BloomFilterClassifier, graph; dataset_id::AbstractString)
    labels = collect(MetaGraphsNext.labels(graph))
    K = eltype(labels)
    coverages = Dict{K, Int}()
    bf = _BloomFilter(clf.nbits, clf.nhashes)
    for label in labels
        coverage, _ = _kmer_evidence(graph, label, dataset_id)
        coverages[label] = coverage
        coverage >= clf.min_count && _bloom_insert!(bf, label)
    end
    verdicts = Dict{K, Bool}()
    scores = Dict{K, Float64}()
    for label in labels
        verdicts[label] = _bloom_contains(bf, label)
        scores[label] = Float64(coverages[label])
    end
    return KmerClassification{K}(
        verdicts, scores, "BloomFilter",
        Dict{Symbol, Any}(:min_count => clf.min_count, :nbits => clf.nbits,
                          :nhashes => clf.nhashes),
    )
end
