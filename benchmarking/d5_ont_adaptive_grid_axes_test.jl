# Pins the D5 axis registry to the 2026-09-11 registration (td-h3yc).
# A mutant that drops 15/21, shrinks the comparator, or changes the dropped
# set must fail HERE — before any HPC job is named 1248 cells.

import Test
include(joinpath(@__DIR__, "d5_ont_adaptive_grid_axes.jl"))
import .D5OntAdaptiveGridAxes

Test.@testset "D5 axis registry matches the 2026-09-11 lock arithmetic" begin
    Test.@test length(D5OntAdaptiveGridAxes.D5_K_POLICIES) == 8
    Test.@test D5OntAdaptiveGridAxes.D5_K_POLICIES[1] === :adaptive
    Test.@test D5OntAdaptiveGridAxes.D5_COMPARATOR_RUNGS == (13, 15, 17, 19, 21, 23, 31)
    Test.@test 15 in D5OntAdaptiveGridAxes.D5_CANDIDATE_LADDER &&
               21 in D5OntAdaptiveGridAxes.D5_CANDIDATE_LADDER
    Test.@test length(D5OntAdaptiveGridAxes.D5_CANDIDATE_LADDER) == 23
    Test.@test length(D5OntAdaptiveGridAxes.D5_DROPPED_RUNGS) == 16
    Test.@test D5OntAdaptiveGridAxes.D5_DROPPED_RUNGS ==
               (29, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)
    Test.@test isempty(intersect(D5OntAdaptiveGridAxes.D5_DROPPED_RUNGS,
        D5OntAdaptiveGridAxes.D5_COMPARATOR_RUNGS))
    Test.@test D5OntAdaptiveGridAxes.n_ont_grid_cells() == 1152
    Test.@test D5OntAdaptiveGridAxes.n_validation_cells() == 96
    Test.@test D5OntAdaptiveGridAxes.n_d5_cells() == 1248
    # 3×4×2×8×2×3 — if you "fix" a count without changing an axis, this product fails.
    Test.@test D5OntAdaptiveGridAxes.n_ont_grid_cells() == 3 * 4 * 2 * 8 * 2 * 3
    v = D5OntAdaptiveGridAxes.validation_cell(29, :Lambda, 42)
    Test.@test v.identity == 99 && v.coverage == 100 && v.corrector === :iterative
    Test.@test v.stratum === :validation
    g = D5OntAdaptiveGridAxes.ont_grid_cell(95, 30, :none, :adaptive, :T4, 123)
    Test.@test g.stratum === :grid && g.k_policy === :adaptive
end
