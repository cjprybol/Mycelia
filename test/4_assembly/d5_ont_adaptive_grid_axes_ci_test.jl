# Runs the D5 axis-registry pin (benchmarking/d5_ont_adaptive_grid_axes_test.jl)
# in CI. test/runtests.jl walks only test/<stage>/ directories, so the registry
# test under benchmarking/ never ran on its own. The registry and its test are
# included unchanged; the module keeps their names out of the shared test scope.

module D5OntAdaptiveGridAxesCITest
include(joinpath(@__DIR__, "..", "..", "benchmarking", "d5_ont_adaptive_grid_axes_test.jl"))
end
