# Default-CI pin for the dependency-free QUAST report parser (bead td-i2bg,
# the #406 slice). benchmarking/ is not otherwise included by Pkg.test(), so
# this file runs the sweep's standalone test script: scale-guard checks, the
# shared parser (including the deferred hardening fixtures), the
# mode_comparison.jl source contract, and the DataFrame schema guard.
#
# The script imports DataFrames and does not import Mycelia, start QUAST, or
# touch the network.
#
# Run:
#   julia --project=. test/4_assembly/quast_report_parsing_test.jl

include(joinpath(
    @__DIR__, "..", "..", "benchmarking",
    "rhizomorph_correction_validation_sweep_test.jl"))
