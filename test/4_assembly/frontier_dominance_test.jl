# Unit test for assess_frontier_dominance (td-4osf Stage 3): the programmatic
# verdict that replaces eyeballing the PR-curve CSV. Exercised here on synthetic
# DataFrames — the real sweep is not run under Pkg.test(). Guards the load-bearing
# STRICT-improvement rule: `calibrated-gap-0.0` ties the baseline (threshold 0
# reverts nothing) and must NOT be reported as dominating.
#
# Run:
#   julia --project=. test/4_assembly/frontier_dominance_test.jl

import Test
import DataFrames

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rhizomorph_correction_pr_curve.jl"))

# Minimal row set the verdict reads: error_rate, point, recall, precision,
# over_correction_rate. Baseline is always "noskip+nogate".
function _frontier(rows)
    DataFrames.DataFrame(
        error_rate = [r[1] for r in rows],
        point = [r[2] for r in rows],
        recall = [r[3] for r in rows],
        precision = [r[4] for r in rows],
        over_correction_rate = [r[5] for r in rows])
end

Test.@testset "assess_frontier_dominance" begin
    # --- A genuine dominating point: recall retained, over strictly down, precision up
    df = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-gap-0.0", 0.90, 0.70, 0.030),   # ties baseline (must NOT win)
        (0.10, "calibrated-gap-1.0", 0.89, 0.80, 0.010)   # strict improvement
    ])
    v = assess_frontier_dominance(df; err = 0.10)
    Test.@test v.dominates
    Test.@test v.best.point == "calibrated-gap-1.0"
    Test.@test all(c -> c.point != "calibrated-gap-0.0", v.candidates)   # tie excluded

    # --- gap-0.0 alone (ties baseline) does NOT dominate (strict over-correction)
    df_tie = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-gap-0.0", 0.90, 0.70, 0.030)
    ])
    Test.@test !assess_frontier_dominance(df_tie; err = 0.10).dominates

    # --- recall collapse ⇒ not dominating even though over drops
    df_crash = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-gap-4.0", 0.40, 0.85, 0.005)
    ])
    Test.@test !assess_frontier_dominance(df_crash; err = 0.10).dominates

    # --- NaN candidate metrics are excluded (never a false positive)
    df_nan = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-gap-1.0", NaN, NaN, NaN)
    ])
    Test.@test !assess_frontier_dominance(df_nan; err = 0.10).dominates

    # --- NaN BASELINE ⇒ non-dominating with an explicit reason (not a false pass)
    df_nanbase = _frontier([
        (0.10, "noskip+nogate", NaN, NaN, NaN),
        (0.10, "calibrated-gap-1.0", 0.89, 0.80, 0.010)
    ])
    vnb = assess_frontier_dominance(df_nanbase; err = 0.10)
    Test.@test !vnb.dominates
    Test.@test occursin("non-finite", vnb.reason)

    # --- missing baseline ⇒ non-dominating with a reason
    df_nobase = _frontier([(0.10, "calibrated-gap-1.0", 0.89, 0.80, 0.010)])
    vmb = assess_frontier_dominance(df_nobase; err = 0.10)
    Test.@test !vmb.dominates
    Test.@test occursin("no noskip+nogate baseline", vmb.reason)
end
