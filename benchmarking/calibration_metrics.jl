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
    block_starts = Int[]
    block_ends = Int[]
    block_sums = Float64[]
    block_counts = Int[]
    for i in eachindex(sorted_scores)
        push!(block_starts, i)
        push!(block_ends, i)
        push!(block_sums, sorted_labels[i] ? 1.0 : 0.0)
        push!(block_counts, 1)
        while length(block_sums) >= 2 &&
              block_sums[end - 1] / block_counts[end - 1] >
              block_sums[end] / block_counts[end]
            block_ends[end - 1] = block_ends[end]
            block_sums[end - 1] += block_sums[end]
            block_counts[end - 1] += block_counts[end]
            pop!(block_starts)
            pop!(block_ends)
            pop!(block_sums)
            pop!(block_counts)
        end
    end
    return (
        thresholds = [sorted_scores[i] for i in block_ends],
        probabilities = block_sums ./ block_counts,
    )
end

"""Map raw `scores` through a fitted isotonic probability map."""
function predict_isotonic(model::NamedTuple,
        scores::AbstractVector{<:Real})::Vector{Float64}
    isempty(model.thresholds) && throw(ArgumentError("isotonic model is empty"))
    return [model.probabilities[clamp(
                searchsortedfirst(model.thresholds, Float64(score)),
                1,
                length(model.thresholds),
            )] for score in scores]
end
