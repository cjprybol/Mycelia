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

    # --- Guards -------------------------------------------------------------
    Test.@test_throws ArgumentError fit_isotonic_map([0.1, Inf], [false, true])
    Test.@test_throws ArgumentError fit_logistic_map(Float64[], Bool[])
    Test.@test_throws ArgumentError fit_binned_map([0.1, 0.2], [true])   # length mismatch
    Test.@test_throws ArgumentError logistic_gap_for_probability(lm, 1.0)  # p not in (0,1)

    # --- Degenerate range: all-equal scores ⇒ single-bin global mean --------
    flat = fill(1.0, 10)
    flat_labels = [true, false, true, false, true, false, true, false, true, false]
    bm = fit_binned_map(flat, flat_labels)
    Test.@test all(isapprox(0.5), predict_binned(bm, flat))
end
