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
function _frontier(rows::AbstractVector{<:Tuple})::DataFrames.DataFrame
    points = [r[2] for r in rows]
    gate_requested = [startswith(point, "calibrated-") for point in points]
    DataFrames.DataFrame(
        error_rate = [r[1] for r in rows],
        point = points,
        ok = [length(r) >= 6 ? r[6] : true for r in rows],
        n_reads = [length(r) >= 7 ? r[7] : 100 for r in rows],
        reads_scored = [length(r) >= 8 ? r[8] : 100 for r in rows],
        injected = [length(r) >= 9 ? r[9] : 100 for r in rows],
        recall = [r[3] for r in rows],
        precision = [r[4] for r in rows],
        over_correction_rate = [r[5] for r in rows],
        gate_requested = gate_requested,
        gate_effective = gate_requested,
        gate_executed = gate_requested,
        feature_contract_skips = zeros(Int, length(rows)),
        candidate_evaluations = Int.(gate_requested),
        feature_events = Int.(gate_requested),
        gate_reverts = zeros(Int, length(rows)),
        structural_errors = zeros(Int, length(rows)),
        unkmerizable_errors = zeros(Int, length(rows)),
        model_digest = fill("", length(rows)),
        artifact_digest = fill("", length(rows)))
end

function _frontier_throws_message(f::Function, fragment::AbstractString)::Bool
    try
        f()
    catch error
        return error isa ArgumentError && occursin(fragment, error.msg)
    end
    return false
end

Test.@testset "calibrated probability operating points" begin
    defaults = withenv("MYCELIA_RPC_PROB_THRESHOLDS" => nothing) do
        calibrated_probability_points()
    end
    Test.@test [point.calibrated_probability_threshold for point in defaults] ==
               [0.0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 0.9]
    Test.@test all(point -> startswith(point.name, "calibrated-prob-"), defaults)

    overridden = withenv("MYCELIA_RPC_PROB_THRESHOLDS" => "0.4,0.85") do
        calibrated_probability_points()
    end
    Test.@test [point.calibrated_probability_threshold for point in overridden] ==
               [0.4, 0.85]
    Test.@test _frontier_throws_message(
        () -> withenv("MYCELIA_RPC_PROB_THRESHOLDS" => "-0.1,0.5") do
            calibrated_probability_points()
        end,
        "must be in [0, 1]")
    Test.@test _frontier_throws_message(
        () -> withenv("MYCELIA_RPC_ASSIGNED_Q" => "94") do
            run_pr_curve()
        end,
        "supported Phred range 0:93")

    model = Mycelia.CorrectionConfidenceModel(
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    mktempdir() do dir
        model_path = persist_runtime_probability_model(model, dir)
        Test.@test isfile(model_path)
        lines = readlines(model_path)
        Test.@test first(lines) == "term,coefficient"
        Test.@test "collapsed_frontier,6.0" in lines
        digest = artifact_digest(dir)
        Test.@test length(digest) == 64
        Test.@test digest == artifact_digest(dir)
        write(joinpath(dir, "extra.txt"), "digest changes")
        Test.@test digest != artifact_digest(dir)
    end

    selected = select_operating_points(PR_OPERATING_POINTS, "noskip,scalable")
    Test.@test [point.name for point in selected] == ["scalable", "noskip"]
    Test.@test _frontier_throws_message(
        () -> select_operating_points(PR_OPERATING_POINTS, "not-a-point"),
        "unknown MYCELIA_RPC_POINTS")
    Test.@test _frontier_throws_message(
        () -> select_operating_points(PR_OPERATING_POINTS, "noskip,"),
        "empty point name")
    Test.@test _frontier_throws_message(
        () -> select_operating_points(NamedTuple[], ""),
        "no precision-recall operating points")
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

    # --- NaN candidate metrics make the requested sweep explicitly incomplete
    df_nan = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-gap-1.0", NaN, NaN, NaN)
    ])
    nan_candidate = assess_frontier_dominance(df_nan; err = 0.10)
    Test.@test !nan_candidate.dominates
    Test.@test occursin("incomplete candidate sweep", nan_candidate.reason)
    Test.@test occursin("non-finite metrics", nan_candidate.reason)

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

    # --- candidate-prefix parameter selects the probability family independently
    df_probability = _frontier([
        (0.05, "noskip+nogate", 0.96, 0.82, 0.020),
        (0.05, "calibrated-gap-1.0", 0.95, 0.90, 0.005),
        (0.05, "calibrated-prob-0.0", 0.96, 0.82, 0.020),
        (0.05, "calibrated-prob-0.8", 0.95, 0.91, 0.004),
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.0", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009)
    ])
    default_gap = assess_frontier_dominance(df_probability; err = 0.05)
    probability_05 = assess_frontier_dominance(
        df_probability; err = 0.05, candidate_prefix = "calibrated-prob-")
    probability_10 = assess_frontier_dominance(
        df_probability; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test default_gap.dominates
    Test.@test default_gap.best.point == "calibrated-gap-1.0"
    Test.@test probability_05.dominates
    Test.@test probability_05.best.point == "calibrated-prob-0.8"
    Test.@test probability_10.dominates
    Test.@test probability_10.best.point == "calibrated-prob-0.8"
    Test.@test all(candidate -> candidate.point != "calibrated-prob-0.0",
        probability_10.candidates)

    # --- probability-family NaNs and ties remain non-dominating
    df_probability_bad = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.0", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.9", NaN, NaN, NaN)
    ])
    Test.@test !assess_frontier_dominance(
        df_probability_bad; err = 0.10,
        candidate_prefix = "calibrated-prob-").dominates

    # --- a superficially winning candidate is ineligible if its run failed
    df_failed_candidate = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009, false, 100, 100)
    ])
    Test.@test !assess_frontier_dominance(
        df_failed_candidate; err = 0.10,
        candidate_prefix = "calibrated-prob-").dominates

    # --- a superficially winning candidate is ineligible if not all reads scored
    df_partial_candidate = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009, true, 100, 99)
    ])
    Test.@test !assess_frontier_dominance(
        df_partial_candidate; err = 0.10,
        candidate_prefix = "calibrated-prob-").dominates

    # --- the baseline itself must be a complete, successful run
    df_failed_baseline = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030, false, 100, 100),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009)
    ])
    failed_baseline = assess_frontier_dominance(
        df_failed_baseline; err = 0.10,
        candidate_prefix = "calibrated-prob-")
    Test.@test !failed_baseline.dominates
    Test.@test occursin("failed or partial", failed_baseline.reason)

    df_partial_baseline = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030, true, 100, 99),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009)
    ])
    partial_baseline = assess_frontier_dominance(
        df_partial_baseline; err = 0.10,
        candidate_prefix = "calibrated-prob-")
    Test.@test !partial_baseline.dominates
    Test.@test occursin("failed or partial", partial_baseline.reason)

    # --- no requested family is not a valid negative result
    df_empty_family = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "scalable", 0.40, 1.00, 0.000)
    ])
    empty_family = assess_frontier_dominance(
        df_empty_family; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !empty_family.dominates
    Test.@test occursin("no calibrated-prob-", empty_family.reason)

    # --- all-invalid and mixed valid+invalid sweeps both fail closed
    df_all_invalid = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.5", NaN, NaN, NaN),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009, false, 100, 100)
    ])
    all_invalid = assess_frontier_dominance(
        df_all_invalid; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !all_invalid.dominates
    Test.@test occursin("incomplete candidate sweep", all_invalid.reason)
    Test.@test occursin("calibrated-prob-0.5", all_invalid.reason)
    Test.@test occursin("calibrated-prob-0.8", all_invalid.reason)

    df_mixed = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009),
        (0.10, "calibrated-prob-0.9", NaN, NaN, NaN)
    ])
    mixed = assess_frontier_dominance(
        df_mixed; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !mixed.dominates
    Test.@test isnothing(mixed.best)
    Test.@test occursin("calibrated-prob-0.9", mixed.reason)

    # --- zero-read rows cannot establish either a baseline or candidate verdict
    df_zero_baseline = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030, true, 0, 0),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009)
    ])
    zero_baseline = assess_frontier_dominance(
        df_zero_baseline; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !zero_baseline.dominates
    Test.@test occursin("failed or partial", zero_baseline.reason)

    df_zero_candidate = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009, true, 0, 0)
    ])
    zero_candidate = assess_frontier_dominance(
        df_zero_candidate; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !zero_candidate.dominates
    Test.@test occursin("incomplete candidate sweep", zero_candidate.reason)

    # --- internally complete candidates must still match the baseline cohort
    df_smaller_cohort = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030, true, 200, 200, 1_000),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009, true, 20, 20, 100)
    ])
    smaller_cohort = assess_frontier_dominance(
        df_smaller_cohort; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !smaller_cohort.dominates
    Test.@test occursin("cohort-size mismatch", smaller_cohort.reason)

    df_injected_mismatch = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030, true, 100, 100, 1_000),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009, true, 100, 100, 999)
    ])
    injected_mismatch = assess_frontier_dominance(
        df_injected_mismatch; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !injected_mismatch.dominates
    Test.@test occursin("injected-error mismatch", injected_mismatch.reason)

    # --- configured-but-inactive gates and contract skips fail closed
    df_inactive = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009)
    ])
    df_inactive.candidate_evaluations[2] = 0
    df_inactive.gate_executed[2] = false
    inactive = assess_frontier_dominance(
        df_inactive; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !inactive.dominates
    Test.@test occursin("zero candidates", inactive.reason)

    df_skips = _frontier([
        (0.10, "noskip+nogate", 0.90, 0.70, 0.030),
        (0.10, "calibrated-prob-0.8", 0.89, 0.81, 0.009)
    ])
    df_skips.feature_contract_skips[2] = 1
    strict_skips = assess_frontier_dominance(
        df_skips; err = 0.10, candidate_prefix = "calibrated-prob-")
    Test.@test !strict_skips.dominates
    Test.@test occursin("exceed bound 0", strict_skips.reason)
    Test.@test assess_frontier_dominance(
        df_skips; err = 0.10, candidate_prefix = "calibrated-prob-",
        max_contract_skips = 1).dominates
end
