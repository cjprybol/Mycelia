# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/qualmer_profile_quality_verdict_test.jl")'
# ```
#
# Regression cover for the qualmer profile quality GATE decision (td-5olt).
#
# The gate in `benchmarking/qualmer_profile_quality_compare.jl` is the decision-2
# authoritative check for "does aggregate qualmer storage preserve correction
# quality at real genome scale". Before td-5olt its error branch and its
# "DEGRADED" branch only printed, so the script exited 0 no matter what — an
# unfalsifiable gate. These tests pin the three ways it must now go red, and the
# one way it must stay green (a gate that always fails is equally useless).
#
# The decision is deliberately dependency-free, so this runs in milliseconds and
# needs no assembly, no ART, and no HPC.

import Test

include(joinpath(dirname(dirname(@__DIR__)), "benchmarking",
    "qualmer_profile_quality_verdict.jl"))

Test.@testset "qualmer profile quality gate decision" begin
    requested = [:full, :ultralight_quality]

    Test.@testset "healthy run PASSES" begin
        results = [(:full, 1, 0.9870, 5316, 12.0),
            (:ultralight_quality, 1, 0.9870, 5316, 9.0)]
        verdicts, failures = qualmer_qc_gate_verdict(results, requested)
        Test.@test isempty(failures)
        Test.@test length(verdicts) == 1
        Test.@test occursin("PARITY", only(verdicts))
        Test.@test qualmer_qc_gate_check(results, requested; io = IOBuffer())
    end

    Test.@testset "within-tolerance shortfall still PASSES" begin
        # 0.005 below :full is inside the 0.01 tolerance — must not go red, or the
        # gate becomes noise and gets ignored.
        results = [(:full, 1, 0.9870, 5316, 12.0),
            (:ultralight_quality, 1, 0.9820, 5300, 9.0)]
        _, failures = qualmer_qc_gate_verdict(results, requested)
        Test.@test isempty(failures)
    end

    Test.@testset "DEGRADED profile FAILS" begin
        results = [(:full, 1, 0.9870, 5316, 12.0),
            (:ultralight_quality, 7, 0.6100, 1200, 9.0)]
        verdicts, failures = qualmer_qc_gate_verdict(results, requested)
        Test.@test occursin("DEGRADED", only(verdicts))
        Test.@test length(failures) == 1
        Test.@test occursin("DEGRADED", only(failures))
        Test.@test !qualmer_qc_gate_check(results, requested; io = IOBuffer())
    end

    Test.@testset "errored profile (absent from results) FAILS" begin
        # The benchmark's catch branch `continue`s, so an errored profile simply
        # never lands in `results`. That absence is the failure signal.
        results = [(:full, 1, 0.9870, 5316, 12.0)]
        _, failures = qualmer_qc_gate_verdict(results, requested)
        Test.@test length(failures) == 1
        Test.@test occursin("NO result", only(failures))
        Test.@test !qualmer_qc_gate_check(results, requested; io = IOBuffer())
    end

    Test.@testset "missing :full baseline FAILS" begin
        # Without the baseline there is nothing to claim parity against, so
        # reporting PASS would be a claim the run cannot support.
        results = [(:ultralight_quality, 1, 0.9870, 5316, 9.0)]
        _, failures = qualmer_qc_gate_verdict(results, requested)
        Test.@test any(f -> occursin(":full is absent", f), failures)
        Test.@test !qualmer_qc_gate_check(results, requested; io = IOBuffer())
    end

    Test.@testset "total assembly failure FAILS" begin
        results = Tuple{Symbol, Int, Float64, Int, Float64}[]
        _, failures = qualmer_qc_gate_verdict(results, requested)
        Test.@test !isempty(failures)
        Test.@test !qualmer_qc_gate_check(results, requested; io = IOBuffer())
    end
end

println("✓ qualmer profile quality gate decision tests completed")
