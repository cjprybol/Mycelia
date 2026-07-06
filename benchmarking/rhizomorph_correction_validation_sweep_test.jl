# Unit test for the Rhizomorph correction-validation SCALE GUARD.
#
# Verifies the core safety property of the sweep: a toy-scale configuration
# triggers SMOKE-ONLY (verdict NOT allowed), while a realistic large-scale
# configuration would pass (verdict allowed). This test is intentionally
# dependency-free — it includes only the pure guard helpers, so it runs in
# milliseconds without loading Mycelia, hitting the network, or assembling.
#
# Run:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/rhizomorph_correction_validation_sweep_test.jl

import Test

include(joinpath(@__DIR__, "rhizomorph_scale_guard.jl"))

Test.@testset "rhizomorph correction validation scale guard" begin
    # --- Toy-scale config triggers SMOKE-ONLY --------------------------------
    # Synthetic 2 kb genome at 10x effective coverage = 20,000 sequenced bases,
    # far below the 1 Mbase floor. This is exactly the local smoke configuration.
    toy_coverage = 10.0
    toy_genome_len = 2_000
    Test.@test scale_metric_bases(toy_coverage, toy_genome_len) == 20_000.0
    Test.@test scale_verdict_allowed(toy_coverage, toy_genome_len) == false

    # A real genome at trivially low coverage is also SMOKE-ONLY: Lambda (48,502 bp)
    # at 1x = 48,502 bases < floor. Guards against "big genome" being mistaken for
    # "enough data".
    Test.@test scale_verdict_allowed(1.0, 48_502) == false

    # --- Large-scale config passes -------------------------------------------
    # Lambda phage (48,502 bp) at 30x = 1,455,060 sequenced bases, above the floor.
    lambda_len = 48_502
    Test.@test scale_metric_bases(30.0, lambda_len) > SCALE_FLOOR_BASES
    Test.@test scale_verdict_allowed(30.0, lambda_len) == true

    # --- Boundary behaviour ---------------------------------------------------
    # Exactly at the floor is allowed (>=), just below is not.
    Test.@test scale_verdict_allowed(1.0, SCALE_FLOOR_BASES) == true
    Test.@test scale_verdict_allowed(1.0, SCALE_FLOOR_BASES - 1) == false

    # --- Custom floor override respected -------------------------------------
    # The harness exposes MYCELIA_RGV_SCALE_FLOOR; the guard must honor an
    # explicit floor argument in both directions.
    Test.@test scale_verdict_allowed(10.0, 2_000; floor = 10_000) == true
    Test.@test scale_verdict_allowed(10.0, 2_000; floor = 100_000) == false
end
