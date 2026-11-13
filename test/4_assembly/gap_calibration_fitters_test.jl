# Gap-calibration fitters (td-4osf Stage 2): the three probability maps compared
# head-to-head — isotonic (PAVA), logistic (Platt), and equal-width binning.
#
# These are pure (Base-only) helpers in benchmarking/calibration_metrics.jl; this
# CI test pins their contract (monotone recovery on separable data, probability
# range, logistic inverse round-trip, degenerate/short-input guards) inside the
# Pkg.test() sweep, since benchmarking/ is not otherwise run in CI.
#
# Run:
#   julia --project=. test/4_assembly/gap_calibration_fitters_test.jl

import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "calibration_metrics.jl"))

Test.@testset "gap calibration fitters (isotonic / logistic / binned)" begin
    # Separable synthetic signal: a large gap ⇒ the corrected base is right.
    gaps = collect(range(0.0, 5.0; length = 240))
    labels = [g > 2.5 for g in gaps]

    # Discrimination is a property of the raw ranking — identical across models.
    Test.@test auroc(gaps, labels) == 1.0

    for (fit, predict) in ((fit_isotonic_map, predict_isotonic),
        (fit_logistic_map, predict_logistic),
        (fit_binned_map, predict_binned))
        model = fit(gaps, labels)
        p = predict(model, gaps)
        Test.@test length(p) == length(gaps)
        Test.@test all(x -> 0.0 <= x <= 1.0, p)          # valid probabilities
        Test.@test issorted(p)                           # monotone non-decreasing in gap
        Test.@test p[1] < 0.5 < p[end]                   # separates the two classes
        # Every model should be reasonably calibrated on cleanly separable data.
        Test.@test brier_score(p, labels) < 0.2
    end

    # --- Logistic inverse round-trips predict_logistic ----------------------
    lm = fit_logistic_map(gaps, labels)
    for target in (0.25, 0.5, 0.75, 0.9)
        g = logistic_gap_for_probability(lm, target)
        Test.@test isapprox(only(predict_logistic(lm, [g])), target; atol = 1e-6)
    end
    # The decision boundary (P=0.5) lands near the true 2.5 threshold.
    Test.@test isapprox(logistic_gap_for_probability(lm, 0.5), 2.5; atol = 0.3)

    # --- Tied scores share one isotonic probability (PAVA tie-pooling) -------
    # Repeated Viterbi gaps tie; all observations at an identical score must map
    # to their shared empirical probability, not to a duplicate-threshold artifact.
    tie_model = fit_isotonic_map([1.0, 1.0, 1.0, 1.0], [false, true, false, true])
    Test.@test all(isapprox(0.5), predict_isotonic(tie_model, [1.0]))
    # Mixed tied + separable: the two 2.0-scored obs (one T one F) pool to 0.5,
    # strictly above the all-false 1.0 group and below the all-true 3.0 group.
    mix = fit_isotonic_map([1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
        [false, false, false, true, true, true])
    Test.@test predict_isotonic(mix, [1.0])[1] < predict_isotonic(mix, [2.0])[1] <
               predict_isotonic(mix, [3.0])[1]
    Test.@test isapprox(predict_isotonic(mix, [2.0])[1], 0.5)

    # --- Guards (assert both the type AND an identifying message fragment) ----
    throws_msg(f, frag) = begin
        threw = false
        try
            f()
        catch e
            threw = e isa ArgumentError && occursin(frag, e.msg)
        end
        threw
    end
    Test.@test throws_msg(() -> fit_isotonic_map([0.1, Inf], [false, true]), "finite")
    Test.@test throws_msg(() -> fit_logistic_map(Float64[], Bool[]), "non-empty")
    Test.@test throws_msg(() -> fit_binned_map([0.1, 0.2], [true]), "equal length")
    Test.@test throws_msg(() -> logistic_gap_for_probability(lm, 1.0), "open interval")

    # --- Degenerate range: all-equal scores ⇒ single-bin global mean --------
    flat = fill(1.0, 10)
    flat_labels = [true, false, true, false, true, false, true, false, true, false]
    bm = fit_binned_map(flat, flat_labels)
    Test.@test all(isapprox(0.5), predict_binned(bm, flat))

    # --- L2 regularization: slope is well-defined on separable data (td-ufsu) -
    # Without L2, gradient descent on a linearly separable problem lets |b|
    # grow with the iteration cap (measured: 1.16 at 200 iters, 2.11 at 2000).
    # Regularization must keep the 200-iter and 2000-iter slopes close.
    sep_scores = vcat(fill(0.0, 50), fill(5.0, 50))
    sep_labels = vcat(fill(false, 50), fill(true, 50))
    b_short = fit_logistic_map(sep_scores, sep_labels; iters = 200).b
    b_long = fit_logistic_map(sep_scores, sep_labels; iters = 2000).b
    Test.@test isfinite(b_short) && isfinite(b_long)
    Test.@test isapprox(b_short, b_long; rtol = 0.1)

    # 90/10 labels. Regularizing the intercept would pull the fitted mean
    # toward 0.5. With the penalty on the slope only, it stays near the
    # label mean.
    n_majority = 180
    n_minority = 20
    imb_scores = vcat(fill(0.0, n_majority), fill(1.0, n_minority))
    imb_labels = vcat(fill(false, n_majority), fill(true, n_minority))
    imb_model = fit_logistic_map(imb_scores, imb_labels; iters = 4000)
    imb_mean = sum(predict_logistic(imb_model, imb_scores)) / length(imb_scores)
    label_mean = n_minority / (n_majority + n_minority)
    Test.@test isapprox(imb_mean, label_mean; atol = 0.02)
end
