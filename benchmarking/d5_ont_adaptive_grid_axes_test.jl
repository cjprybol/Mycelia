# Pins the D5 axis registry to the 2026-09-11 registration (td-h3yc).
# A mutant that drops 15/21, shrinks the comparator, or changes the dropped
# set must fail HERE — before any HPC job is named 1248 cells.

using Test
include(joinpath(@__DIR__, "d5_ont_adaptive_grid_axes.jl"))
using .D5OntAdaptiveGridAxes

@testset "D5 axis registry matches the 2026-09-11 lock arithmetic" begin
    @test length(D5_K_POLICIES) == 8
    @test D5_K_POLICIES[1] === :adaptive
    @test D5_COMPARATOR_RUNGS == (13, 15, 17, 19, 21, 23, 31)
    @test 15 in D5_CANDIDATE_LADDER && 21 in D5_CANDIDATE_LADDER
    @test length(D5_CANDIDATE_LADDER) == 23
    @test length(D5_DROPPED_RUNGS) == 16
    @test D5_DROPPED_RUNGS == (29, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)
    @test isempty(intersect(D5_DROPPED_RUNGS, D5_COMPARATOR_RUNGS))
    @test n_ont_grid_cells() == 1152
    @test n_validation_cells() == 96
    @test n_d5_cells() == 1248
    # 3×4×2×8×2×3 — if you "fix" a count without changing an axis, this product fails.
    @test n_ont_grid_cells() == 3 * 4 * 2 * 8 * 2 * 3
    v = validation_cell(29, :Lambda, 42)
    @test v.identity == 99 && v.coverage == 100 && v.corrector === :iterative
    @test v.stratum === :validation
    g = ont_grid_cell(95, 30, :none, :adaptive, :T4, 123)
    @test g.stratum === :grid && g.k_policy === :adaptive
end
