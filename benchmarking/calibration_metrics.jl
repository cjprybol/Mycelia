# Pure calibration metrics — NO Mycelia dependency.
# ================================================
#
# Reliability curve, Expected Calibration Error (ECE), AUROC, and Brier score for
# a set of confidence scores vs boolean ground-truth labels. These did not exist
# in the repo (only a Mann-Whitney AUC buried in benchmarking/stage0_*); this is
# the reusable calibration toolkit the Rhizomorph "reference-free correctness
# prediction" work (td-pw7p Priority #2) is built on. Depends only on Julia Base.
#
# Terminology: `confidences` are predicted P(label=true) in [0,1]; `labels` are
# the observed booleans. A well-calibrated predictor has, among all items it
# assigned confidence ~p, a fraction ~p that are actually true.
#
#   reliability_bins(conf, labels; nbins)  -> per-bin (mean_conf, accuracy, count)
#   expected_calibration_error(...)        -> ECE = sum_b (n_b/N) |acc_b - conf_b|
#   auroc(scores, labels)                  -> P(score(true) > score(false)) (rank-based; scores need NOT be in [0,1])
#   brier_score(conf, labels)              -> mean((conf - label)^2)
#
# AUROC measures DISCRIMINATION (does the signal separate true from false) and is
# rank-based, so it accepts any real-valued score. ECE / reliability / Brier
# measure CALIBRATION and require `confidences` already mapped to [0,1]
# probabilities (fit a Platt/isotonic map first for a raw signal).

"""
    reliability_bins(confidences, labels; nbins=10) -> Vector{NamedTuple}

Equal-width binning of `confidences` in [0,1]. Returns one entry per NON-EMPTY
bin: `(lo, hi, mean_conf, accuracy, count)` where `accuracy` is the fraction of
`labels` that are true among items whose confidence fell in `[lo, hi)` (the last
bin is closed, `[lo, hi]`). A calibrated predictor has `accuracy ~ mean_conf` in
every bin.
"""
function reliability_bins(confidences::AbstractVector{<:Real},
        labels::AbstractVector{Bool}; nbins::Int = 10)
    length(confidences) == length(labels) ||
        throw(ArgumentError("confidences and labels must have equal length"))
    nbins >= 1 || throw(ArgumentError("nbins must be >= 1"))
    edges = range(0.0, 1.0; length = nbins + 1)
    out = NamedTuple[]
    for b in 1:nbins
        lo = edges[b]
        hi = edges[b + 1]
        # last bin is closed on the right so confidence==1.0 lands somewhere
        idx = b < nbins ? findall(c -> lo <= c < hi, confidences) :
              findall(c -> lo <= c <= hi, confidences)
        isempty(idx) && continue
        cnt = length(idx)
        mean_conf = sum(confidences[idx]) / cnt
        accuracy = count(labels[idx]) / cnt
        push!(out, (
            lo = lo, hi = hi, mean_conf = mean_conf, accuracy = accuracy, count = cnt))
    end
    return out
end

"""
    expected_calibration_error(confidences, labels; nbins=10) -> Float64

ECE = sum over bins of (bin_count / N) * |bin_accuracy - bin_mean_confidence|.
0 is perfectly calibrated; 1 is worst. Requires `confidences` in [0,1].
"""
function expected_calibration_error(confidences::AbstractVector{<:Real},
        labels::AbstractVector{Bool}; nbins::Int = 10)
    n = length(confidences)
    n == 0 && return NaN
    bins = reliability_bins(confidences, labels; nbins = nbins)
    return sum(bin.count / n * abs(bin.accuracy - bin.mean_conf) for bin in bins; init = 0.0)
end

"""
    auroc(scores, labels) -> Float64

Area under the ROC curve via the Mann-Whitney U statistic:
P(score(a true item) > score(a false item)), with ties counted as 0.5. Returns
NaN if either class is empty (AUROC undefined). Rank-based, so `scores` need not
be probabilities — any real-valued discrimination signal works. 0.5 = no
discrimination, 1.0 = perfect, <0.5 = anti-correlated.
"""
function auroc(scores::AbstractVector{<:Real}, labels::AbstractVector{Bool})
    length(scores) == length(labels) ||
        throw(ArgumentError("scores and labels must have equal length"))
    pos = scores[labels]
    neg = scores[.!labels]
    (isempty(pos) || isempty(neg)) && return NaN
    wins = 0.0
    for p in pos, q in neg

        wins += p > q ? 1.0 : (p == q ? 0.5 : 0.0)
    end
    return wins / (length(pos) * length(neg))
end

"""
    brier_score(confidences, labels) -> Float64

Mean squared error of the predicted probabilities: mean((conf - label)^2), where
label is 1.0 for true / 0.0 for false. 0 is best; combines calibration and
refinement. Requires `confidences` in [0,1].
"""
function brier_score(confidences::AbstractVector{<:Real}, labels::AbstractVector{Bool})
    n = length(confidences)
    n == 0 && return NaN
    length(confidences) == length(labels) ||
        throw(ArgumentError("confidences and labels must have equal length"))
    return sum((confidences[i] - (labels[i] ? 1.0 : 0.0))^2 for i in 1:n) / n
end

"""
    fit_isotonic_map(scores, labels) -> NamedTuple

Fit a non-decreasing, piecewise-constant probability map with the pooled
adjacent violators algorithm. The returned `thresholds` are sorted score upper
bounds and `probabilities` are the corresponding empirical probabilities.
"""
function fit_isotonic_map(scores::AbstractVector{<:Real},
        labels::AbstractVector{Bool})::NamedTuple
    length(scores) == length(labels) ||
        throw(ArgumentError("scores and labels must have equal length"))
    isempty(scores) && throw(ArgumentError("scores and labels must be non-empty"))
    all(isfinite, scores) || throw(ArgumentError("scores must be finite"))

    order = sortperm(scores)
    sorted_scores = Float64.(scores[order])
    sorted_labels = labels[order]
    # Pool TIED scores into one point BEFORE PAVA. A calibration map is a function
    # of the score, so all observations with an identical score must share one
    # probability; without this, equal scores whose labels happen to be
    # non-decreasing land in separate blocks with duplicate thresholds, and
    # `searchsortedfirst` then maps a repeated score to the first block's
    # probability instead of the shared empirical one (repeated Viterbi gaps
    # tie routinely, so this is not a corner case).
    point_scores = Float64[]
    point_sums = Float64[]
    point_counts = Int[]
    for i in eachindex(sorted_scores)
        if !isempty(point_scores) && sorted_scores[i] == point_scores[end]
            point_sums[end] += sorted_labels[i] ? 1.0 : 0.0
            point_counts[end] += 1
        else
            push!(point_scores, sorted_scores[i])
            push!(point_sums, sorted_labels[i] ? 1.0 : 0.0)
            push!(point_counts, 1)
        end
    end
    # PAVA over the unique-score points → non-decreasing, unique-threshold blocks.
    block_scores = Float64[]
    block_sums = Float64[]
    block_counts = Int[]
    for j in eachindex(point_scores)
        push!(block_scores, point_scores[j])
        push!(block_sums, point_sums[j])
        push!(block_counts, point_counts[j])
        while length(block_sums) >= 2 &&
            block_sums[end - 1] / block_counts[end - 1] >
            block_sums[end] / block_counts[end]
            block_scores[end - 1] = block_scores[end]   # block threshold = its max score
            block_sums[end - 1] += block_sums[end]
            block_counts[end - 1] += block_counts[end]
            pop!(block_scores)
            pop!(block_sums)
            pop!(block_counts)
        end
    end
    return (thresholds = block_scores, probabilities = block_sums ./ block_counts)
end

"""Map raw `scores` through a fitted isotonic probability map."""
function predict_isotonic(model::NamedTuple,
        scores::AbstractVector{<:Real})::Vector{Float64}
    isempty(model.thresholds) && throw(ArgumentError("isotonic model is empty"))
    return [model.probabilities[clamp(
                searchsortedfirst(model.thresholds, Float64(score)),
                1,
                length(model.thresholds)
            )] for score in scores]
end

"""
    fit_logistic_map(scores, labels; iters=500, lr=0.1) -> NamedTuple

Fit a Platt-style single-feature logistic probability map
`P(true | score) = 1 / (1 + exp(-(a + b*score)))` by gradient descent on the
standardized score, folded back to raw-score space. Returns `(a, b)` — the
raw-space intercept and slope. Mirrors the standardize→fit→unfold structure of
the repo's `fit_logistic_fusion`, but with a single feature and no Mycelia
dependency. Two parameters, strictly monotone, so it inverts cleanly to a single
raw gap cutoff (see `logistic_gap_for_probability`).
"""
function fit_logistic_map(scores::AbstractVector{<:Real},
        labels::AbstractVector{Bool}; iters::Int = 500, lr::Float64 = 0.1)::NamedTuple
    length(scores) == length(labels) ||
        throw(ArgumentError("scores and labels must have equal length"))
    isempty(scores) && throw(ArgumentError("scores and labels must be non-empty"))
    all(isfinite, scores) || throw(ArgumentError("scores must be finite"))
    x = Float64.(scores)
    y = [l ? 1.0 : 0.0 for l in labels]
    n = length(x)
    μ = sum(x) / n
    σ = sqrt(max(sum((xi - μ)^2 for xi in x) / n, eps()))
    z = (x .- μ) ./ σ
    α = 0.0
    β = 0.0
    for _ in 1:iters
        gα = 0.0
        gβ = 0.0
        for i in 1:n
            p = 1.0 / (1.0 + exp(-(α + β * z[i])))
            err = p - y[i]
            gα += err
            gβ += err * z[i]
        end
        α -= lr * gα / n
        β -= lr * gβ / n
    end
    # Fold the standardized (α, β) back to raw-score space:
    # α + β*(x-μ)/σ = (α - β*μ/σ) + (β/σ)*x.
    return (a = α - β * μ / σ, b = β / σ)
end

"""
    fit_logistic_map(features, labels; iters=500, lr=0.1) -> NamedTuple

Fit a multi-feature logistic probability map with samples in rows and features
in columns. Each feature is standardized independently during optimization and
the fitted coefficients are folded back to raw-feature space. Returns
`(a, b::Vector{Float64})`, where `a` is the intercept and `b[j]` multiplies
`features[:, j]`. Constant columns receive a zero raw-space coefficient.
"""
function fit_logistic_map(features::AbstractMatrix{<:Real},
        labels::AbstractVector{Bool}; iters::Int = 500, lr::Float64 = 0.1)::NamedTuple
    n_samples, n_features = size(features)
    n_samples == length(labels) ||
        throw(ArgumentError("feature rows and labels must have equal length"))
    n_samples >= 1 || throw(ArgumentError("features and labels must be non-empty"))
    n_features >= 1 || throw(ArgumentError("features must have at least one column"))
    all(isfinite, features) || throw(ArgumentError("features must be finite"))
    iters >= 1 || throw(ArgumentError("iters must be >= 1"))
    isfinite(lr) && lr > 0.0 || throw(ArgumentError("lr must be finite and > 0"))

    x = Float64.(features)
    y = [label ? 1.0 : 0.0 for label in labels]
    means = zeros(Float64, n_features)
    scales = ones(Float64, n_features)
    standardized = zeros(Float64, n_samples, n_features)
    for j in 1:n_features
        means[j] = sum(x[:, j]) / n_samples
        variance = sum((x[i, j] - means[j])^2 for i in 1:n_samples) / n_samples
        if variance > eps(Float64)
            scales[j] = sqrt(variance)
            for i in 1:n_samples
                standardized[i, j] = (x[i, j] - means[j]) / scales[j]
            end
        end
    end

    intercept = 0.0
    coefficients = zeros(Float64, n_features)
    coefficient_gradient = zeros(Float64, n_features)
    for _ in 1:iters
        intercept_gradient = 0.0
        fill!(coefficient_gradient, 0.0)
        for i in 1:n_samples
            linear_predictor = intercept
            for j in 1:n_features
                linear_predictor += coefficients[j] * standardized[i, j]
            end
            probability = 1.0 / (1.0 + exp(-linear_predictor))
            residual = probability - y[i]
            intercept_gradient += residual
            for j in 1:n_features
                coefficient_gradient[j] += residual * standardized[i, j]
            end
        end
        intercept -= lr * intercept_gradient / n_samples
        for j in 1:n_features
            coefficients[j] -= lr * coefficient_gradient[j] / n_samples
        end
    end

    raw_coefficients = coefficients ./ scales
    raw_intercept = intercept -
                    sum(coefficients[j] * means[j] / scales[j] for j in 1:n_features)
    return (a = raw_intercept, b = raw_coefficients)
end

"""Map raw `scores` through a fitted logistic probability map."""
function predict_logistic(model::NamedTuple,
        scores::AbstractVector{<:Real})::Vector{Float64}
    return [1.0 / (1.0 + exp(-(model.a + model.b * Float64(score)))) for score in scores]
end

"""Map raw row-oriented `features` through a fitted multi-feature logistic map."""
function predict_logistic(model::NamedTuple,
        features::AbstractMatrix{<:Real})::Vector{Float64}
    hasproperty(model, :a) && hasproperty(model, :b) ||
        throw(ArgumentError("logistic model must contain a and b"))
    model.b isa AbstractVector ||
        throw(ArgumentError("multi-feature logistic model b must be a vector"))
    size(features, 2) == length(model.b) ||
        throw(ArgumentError("feature columns must match model coefficients"))
    all(isfinite, features) || throw(ArgumentError("features must be finite"))
    isfinite(model.a) && all(isfinite, model.b) ||
        throw(ArgumentError("logistic model coefficients must be finite"))

    probabilities = Vector{Float64}(undef, size(features, 1))
    for i in axes(features, 1)
        linear_predictor = Float64(model.a)
        for j in axes(features, 2)
            linear_predictor += Float64(model.b[j]) * Float64(features[i, j])
        end
        probabilities[i] = 1.0 / (1.0 + exp(-linear_predictor))
    end
    return probabilities
end

"""
    logistic_gap_for_probability(model, p) -> Float64

Invert the logistic map: the raw score at which `predict_logistic` equals `p`,
usable as a calibrated gap cutoff. `logit(p) = a + b*score` ⇒
`score = (logit(p) - a) / b`. Requires a non-degenerate slope.
"""
function logistic_gap_for_probability(model::NamedTuple, p::Real)::Float64
    0.0 < p < 1.0 || throw(ArgumentError("p must be in the open interval (0, 1)"))
    model.b == 0.0 &&
        throw(ArgumentError("degenerate logistic map (zero slope) cannot be inverted"))
    return (log(p / (1 - p)) - model.a) / model.b
end

"""
    fit_binned_map(scores, labels; nbins=10) -> NamedTuple

Fit a coarse probability map by partitioning the score range into `nbins`
equal-width bins and taking each bin's empirical fraction of true labels. Returns
`edges` (nbins+1 boundaries) and `probabilities` (one per bin; an empty bin
inherits the global mean so the map is always defined). The coarsest of the three
calibration models — zero model assumptions, but edge-unstable on sparse bins.
"""
function fit_binned_map(scores::AbstractVector{<:Real},
        labels::AbstractVector{Bool}; nbins::Int = 10)::NamedTuple
    length(scores) == length(labels) ||
        throw(ArgumentError("scores and labels must have equal length"))
    isempty(scores) && throw(ArgumentError("scores and labels must be non-empty"))
    all(isfinite, scores) || throw(ArgumentError("scores must be finite"))
    nbins >= 1 || throw(ArgumentError("nbins must be >= 1"))
    x = Float64.(scores)
    lo, hi = extrema(x)
    global_mean = sum(l ? 1.0 : 0.0 for l in labels) / length(labels)
    # Degenerate range (all scores equal): a single bin at the global mean.
    hi <= lo && return (edges = [lo, hi], probabilities = [global_mean])
    width = (hi - lo) / nbins
    sums = zeros(Float64, nbins)
    counts = zeros(Int, nbins)
    for i in eachindex(x)
        b = clamp(Int(floor((x[i] - lo) / width)) + 1, 1, nbins)
        sums[b] += labels[i] ? 1.0 : 0.0
        counts[b] += 1
    end
    probs = [counts[b] == 0 ? global_mean : sums[b] / counts[b] for b in 1:nbins]
    return (edges = [lo + b * width for b in 0:nbins], probabilities = probs)
end

"""Map raw `scores` through a fitted equal-width binned probability map."""
function predict_binned(model::NamedTuple,
        scores::AbstractVector{<:Real})::Vector{Float64}
    nbins = length(model.probabilities)
    lo = model.edges[1]
    hi = model.edges[end]
    width = (hi - lo) / nbins
    return [model.probabilities[width <= 0 ? 1 :
                                clamp(Int(floor((Float64(score) - lo) / width)) + 1, 1, nbins)]
            for score in scores]
end
