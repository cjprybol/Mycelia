# VERDICT gate for the correction-accuracy benchmark (td-womt).
#
# Pins the pure decision in benchmarking/rhizomorph_correction_accuracy_verdict.jl
# inside the default Pkg.test() sweep. benchmarking/ is not otherwise included.
# Fixtures only: no corrector, no download, no Lambda / Lovelace / cluster run.
#
#   (a) an unhealthy null arm (crash, mis-join, or mass drop) cannot emit VERDICT
#   (b) a healthy null arm whose scale and scored checks pass still emits VERDICT
#   (c) the existing local smoke scale (synthetic ~2 kb at 10x) stays SMOKE-ONLY
#   The harness source is pinned too: run_accuracy_benchmark must call
#   correction_accuracy_verdict and record null_scored_fraction. Restoring
#   `verdict_allowed = scale_ok && scored_ok` fails that contract.
#
# Run:
#   julia --project=. test/4_assembly/rhizomorph_correction_accuracy_verdict_test.jl

module CorrectionAccuracyVerdictTest

import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rhizomorph_scale_guard.jl"))
include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_correction_accuracy_verdict.jl"))

# Lambda phage at 30x is the harness's documented VERDICT-scale configuration
# (48,502 bp x 30x ~ 1.46 Mbase, above SCALE_FLOOR_BASES).
const LAMBDA_COVERAGE = 30.0
const LAMBDA_LENGTH = 48_502
# Local smoke, from the benchmark header and the scale-guard test: synthetic
# ~2 kb genome at 10x = 20,000 bases, below the 1 Mbase floor.
const SMOKE_COVERAGE = 10.0
const SMOKE_LENGTH = 2_000

function _cell(; injected_errors::Int, n_reads::Int, null_ok::Bool,
        null_reads_scored::Int, null_corrected_unjoined::Int)
    return NullArmCell(;
        injected_errors = injected_errors,
        n_reads = n_reads,
        null_ok = null_ok,
        null_reads_scored = null_reads_scored,
        null_corrected_unjoined = null_corrected_unjoined)
end

function _healthy_injected(; n_reads::Int = 100, injected_errors::Int = 40)
    return _cell(;
        injected_errors = injected_errors,
        n_reads = n_reads,
        null_ok = true,
        null_reads_scored = n_reads,
        null_corrected_unjoined = 0)
end

function _decide(scale_ok::Bool, agg_scored_fraction::Real, cells::Vector{NullArmCell};
        metric::Real, scale_floor::Real)
    return correction_accuracy_verdict(;
        scale_ok = scale_ok,
        agg_scored_fraction = agg_scored_fraction,
        cells = cells,
        metric = metric,
        scale_floor = scale_floor)
end

function _lambda()
    return (
        scale_ok = scale_verdict_allowed(LAMBDA_COVERAGE, LAMBDA_LENGTH),
        metric = scale_metric_bases(LAMBDA_COVERAGE, LAMBDA_LENGTH),
        scale_floor = Float64(SCALE_FLOOR_BASES))
end

function _smoke()
    return (
        scale_ok = scale_verdict_allowed(SMOKE_COVERAGE, SMOKE_LENGTH),
        metric = scale_metric_bases(SMOKE_COVERAGE, SMOKE_LENGTH),
        scale_floor = Float64(SCALE_FLOOR_BASES))
end

Test.@testset "correction-accuracy VERDICT requires null-arm health" begin
    lambda = _lambda()
    Test.@test lambda.scale_ok
    healthy = [_healthy_injected()]

    Test.@testset "unhealthy null arm cannot emit VERDICT" begin
        # Crash: corrector threw, but every read that did exist would have
        # scored. The scored fraction alone would pass; null_ok must still block.
        crashed = [_cell(;
            injected_errors = 40, n_reads = 100, null_ok = false,
            null_reads_scored = 100, null_corrected_unjoined = 0)]
        crash = _decide(lambda.scale_ok, 1.0, crashed;
            metric = lambda.metric, scale_floor = lambda.scale_floor)
        Test.@test crash.scale_ok
        Test.@test crash.scored_ok
        Test.@test !crash.null_ok
        Test.@test !crash.verdict_allowed
        Test.@test crash.mode == "SMOKE-ONLY"
        Test.@test !occursin("=== VERDICT", crash.line)
        Test.@test occursin("=== SMOKE-ONLY", crash.line)
        Test.@test occursin("null_scored_fraction=1.0", crash.line)

        # Mis-join: null corrector ran, emitted reads, none joined to truth.
        # A sibling cell keeps the pooled scored fraction above the floor, so
        # only the mis-join (not the fraction) can downgrade the run.
        misjoined = [
            _healthy_injected(; n_reads = 100),
            _cell(;
                injected_errors = 10, n_reads = 10, null_ok = true,
                null_reads_scored = 0, null_corrected_unjoined = 10),
        ]
        misjoin = _decide(lambda.scale_ok, 1.0, misjoined;
            metric = lambda.metric, scale_floor = lambda.scale_floor)
        Test.@test misjoin.null_scored_fraction > 0.5
        Test.@test !misjoin.null_ok
        Test.@test !misjoin.verdict_allowed
        Test.@test misjoin.mode == "SMOKE-ONLY"
        Test.@test !occursin("=== VERDICT", misjoin.line)
        Test.@test occursin("null_scored_fraction=", misjoin.line)

        # Mass drop: corrector ran and ids joined, but almost every row was
        # excluded. null_ok stays true; the scored fraction is the gate.
        dropped = [_cell(;
            injected_errors = 40, n_reads = 100, null_ok = true,
            null_reads_scored = 10, null_corrected_unjoined = 0)]
        drop = _decide(lambda.scale_ok, 1.0, dropped;
            metric = lambda.metric, scale_floor = lambda.scale_floor)
        Test.@test drop.null_ok
        Test.@test !drop.null_scored_ok
        Test.@test drop.null_scored_fraction == 0.1
        Test.@test !drop.verdict_allowed
        Test.@test drop.mode == "SMOKE-ONLY"
        Test.@test !occursin("=== VERDICT", drop.line)
        Test.@test occursin("=== SMOKE-ONLY", drop.line)
        Test.@test occursin("null_scored_fraction=0.1", drop.line)
    end

    Test.@testset "healthy null arm with passing scale and scored checks emits VERDICT" begin
        # Same scale and real-arm scored fraction as the unhealthy fixtures.
        # The previous gate was `scale_ok && scored_ok`; a healthy null arm
        # must leave that decision unchanged.
        decision = _decide(lambda.scale_ok, 1.0, healthy;
            metric = lambda.metric, scale_floor = lambda.scale_floor)
        Test.@test decision.scale_ok && decision.scored_ok && decision.null_arm_ok
        Test.@test decision.verdict_allowed == (decision.scale_ok && decision.scored_ok)
        Test.@test decision.verdict_allowed
        Test.@test decision.mode == "VERDICT"
        Test.@test occursin("=== VERDICT", decision.line)
        Test.@test !occursin("SMOKE-ONLY", decision.line)
        Test.@test occursin("null_scored_fraction=1.0", decision.line)
        Test.@test decision.null_scored_fraction == 1.0

        # Control A (no injected errors) is outside the null gate, even when
        # that cell's null arm crashes or mis-joins.
        control_a_only_failure = [
            _cell(;
                injected_errors = 0, n_reads = 100, null_ok = false,
                null_reads_scored = 0, null_corrected_unjoined = 100),
            _healthy_injected(),
        ]
        spared = _decide(lambda.scale_ok, 1.0, control_a_only_failure;
            metric = lambda.metric, scale_floor = lambda.scale_floor)
        Test.@test spared.n_injected == 1
        Test.@test spared.null_arm_ok
        Test.@test spared.mode == "VERDICT"
    end

    Test.@testset "existing smoke path stays SMOKE-ONLY" begin
        smoke = _smoke()
        Test.@test smoke.metric == 20_000.0
        Test.@test !smoke.scale_ok
        # Healthy null and a fully scored real arm: the label comes from the
        # scale floor, which is what a local MYCELIA_RCA_SMOKE run hits.
        decision = _decide(smoke.scale_ok, 1.0, healthy;
            metric = smoke.metric, scale_floor = smoke.scale_floor)
        Test.@test decision.scored_ok
        Test.@test decision.null_arm_ok
        Test.@test !decision.verdict_allowed
        Test.@test decision.mode == "SMOKE-ONLY"
        Test.@test occursin("=== SMOKE-ONLY", decision.line)
        Test.@test !occursin("=== VERDICT", decision.line)
        Test.@test occursin("null_scored_fraction=1.0", decision.line)
        Test.@test occursin("scale metric", decision.line)
    end

    Test.@testset "boolean decision matches the fixture gate" begin
        Test.@test correction_accuracy_verdict_allowed(true, true, true, 1.0)
        Test.@test !correction_accuracy_verdict_allowed(true, true, false, 1.0)
        Test.@test !correction_accuracy_verdict_allowed(true, true, true, 0.1)
        Test.@test !correction_accuracy_verdict_allowed(false, true, true, 1.0)
        Test.@test !correction_accuracy_verdict_allowed(true, false, true, 1.0)
        # No injected-error cells -> null scored fraction is NaN and does not block.
        Test.@test correction_accuracy_verdict_allowed(true, true, true, NaN)
        Test.@test null_scored_fraction_ok(0.5, 0.5)
        Test.@test !null_scored_fraction_ok(0.5 - eps(0.5), 0.5)
    end

    # The fixtures above include only the pure helper. This pins the harness
    # call site, so reverting run_accuracy_benchmark to `scale_ok && scored_ok`
    # fails in the default suite without loading Mycelia.
    Test.@testset "harness source uses the verdict helper" begin
        src = read(joinpath(@__DIR__, "..", "..", "benchmarking",
            "rhizomorph_correction_accuracy_benchmark.jl"), String)
        call_at = findfirst("decision = correction_accuracy_verdict(", src)
        use_at = findfirst("verdict_allowed = decision.verdict_allowed", src)
        Test.@test call_at !== nothing
        Test.@test use_at !== nothing
        Test.@test first(call_at) < first(use_at)
        Test.@test !occursin("verdict_allowed = scale_ok && scored_ok", src)
        Test.@test occursin("null_scored_fraction = Float64[]", src)
        Test.@test occursin("null_scored_fraction = mn.scored_fraction", src)
        Test.@test occursin("rows[!, :agg_null_scored_fraction]", src)
        Test.@test occursin("round(decision.null_scored_fraction; digits = 4)", src)
        Test.@test occursin("null_scored_fraction = decision.null_scored_fraction", src)
        summary_prints = count(l -> occursin("println(\"\\n\$(decision.line)\")", l),
            eachline(IOBuffer(src)))
        Test.@test summary_prints == 2
    end
end

end # module CorrectionAccuracyVerdictTest
