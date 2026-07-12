# Pure unit test for calibration_metrics.jl — NO Mycelia dependency, runs in ms.
#
# Run:
#   julia benchmarking/calibration_metrics_test.jl

import Test

include(joinpath(@__DIR__, "calibration_metrics.jl"))

Test.@testset "calibration metrics" begin
    # --- reliability_bins ----------------------------------------------------
    conf = [0.1, 0.2, 0.8, 0.9]
    labels = [false, false, true, true]
    bins = reliability_bins(conf, labels; nbins = 2)
    Test.@test length(bins) == 2
    Test.@test bins[1].count == 2
    Test.@test bins[1].mean_conf ≈ 0.15
    Test.@test bins[1].accuracy == 0.0
    Test.@test bins[2].count == 2
    Test.@test bins[2].mean_conf ≈ 0.85
    Test.@test bins[2].accuracy == 1.0
    # empty bins are dropped (nothing in the middle range here with nbins=2)
    Test.@test all(b -> b.count > 0, bins)

    # --- expected_calibration_error ------------------------------------------
    # bin1: |0 - 0.15|=0.15 weight 0.5; bin2: |1 - 0.85|=0.15 weight 0.5 => 0.15
    Test.@test expected_calibration_error(conf, labels; nbins = 2) ≈ 0.15
    # perfectly-calibrated extremes => ECE 0
    Test.@test expected_calibration_error([0.0, 1.0], [false, true]; nbins = 2) ≈ 0.0
    # a confidently-wrong predictor => ECE 1.0 (conf 1.0 on all-false, 0.0 on all-true)
    Test.@test expected_calibration_error([1.0, 0.0], [false, true]; nbins = 2) ≈ 1.0

    # --- auroc (rank-based discrimination) -----------------------------------
    Test.@test auroc([0.1, 0.2, 0.8, 0.9], [false, false, true, true]) == 1.0   # perfect
    Test.@test auroc([0.9, 0.1], [false, true]) == 0.0                          # reversed
    Test.@test auroc([0.5, 0.5], [true, false]) == 0.5                          # tie => 0.5
    Test.@test isnan(auroc([0.3, 0.7], [true, true]))                           # one class => undefined
    # accepts non-probability scores (rank-based): coverage-like integers
    Test.@test auroc([2, 5, 9, 1], [false, true, true, false]) == 1.0

    # --- brier_score ----------------------------------------------------------
    # mean((0.1)^2,(0.2)^2,(0.2)^2,(0.1)^2) = (0.01+0.04+0.04+0.01)/4 = 0.025
    Test.@test brier_score(conf, labels) ≈ 0.025
    Test.@test brier_score([0.0, 1.0], [false, true]) == 0.0                    # perfect
    Test.@test brier_score([1.0, 0.0], [false, true]) == 1.0                    # worst

    # --- degenerate / empty inputs -------------------------------------------
    Test.@test isnan(expected_calibration_error(Float64[], Bool[]))
    Test.@test isnan(brier_score(Float64[], Bool[]))
    Test.@test isnan(auroc(Float64[], Bool[]))
    Test.@test isempty(reliability_bins(Float64[], Bool[]))

    # --- length-mismatch guards ----------------------------------------------
    Test.@test_throws ArgumentError auroc([0.1, 0.2], [true])
    Test.@test_throws ArgumentError reliability_bins([0.1], [true, false])
end
