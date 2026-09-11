# CI entry point for benchmarking/run_ont_k_sweep_shards_test.jl.
#
# The concurrency-cap suite lives next to the shell driver it tests, but
# runtests.jl's include_all_tests only walks test/1_* .. test/8_*, so nothing
# under benchmarking/ is ever executed by Pkg.test(). This file pins it into the
# CI sweep -- the same pattern as track_a_baseline_benchmark_test.jl for
# benchmarking/track_a_baseline_benchmark.jl and gap_calibration_fitters_test.jl
# for benchmarking/calibration_metrics.jl.
#
# Without this the cap has no CI guard at all. That matters more than usual
# here: the suite is the only thing standing between a regression and a repeat
# of the 2026-08-04 thrashing incident, and a cap is exactly the kind of code
# whose failure is silent -- a broken cap and a working one produce identical
# output on a host with room to spare.
#
# Wrapped in a module for the same reason track_a_baseline_benchmark_test.jl is:
# include_all_tests pulls every test file into one shared Main, and this suite
# defines top-level names (DRIVER, STUB_JULIA, stage_sandbox, run_driver) that
# would otherwise collide with a sibling.
#
# Run:
#   julia --project=. test/4_assembly/ont_k_sweep_shards_driver_test.jl

module ONTKSweepShardsDriverTests

include(joinpath(@__DIR__, "..", "..", "benchmarking", "run_ont_k_sweep_shards_test.jl"))

end
