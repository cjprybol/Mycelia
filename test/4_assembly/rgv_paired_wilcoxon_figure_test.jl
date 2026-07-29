# Unit test for the paired-Wilcoxon FIGURE caption (beads td-59o7, td-9p91).
#
# `paired_figure_caption` was made pure — no CairoMakie call, no IO — precisely so
# the integrity caveats it stamps could be asserted without rendering anything, and
# it shipped with no test. This file supplies one.
#
# Why the caveats matter enough to test: under this repo's
# `1 notebook = 1 figure = 1 slide` model the SVG/PNG is the artifact that reaches
# a deck or a manuscript, detached from `report.md`. A run whose guard was
# overridden, or whose metric definition was asserted by an operator rather than
# observed, previously rendered a caption bit-for-bit identical to a measured run,
# so the artifact that travels furthest carried the least disclosure.
#
# The script is loaded into its own module: it defines `main` / `_fig_arg` at top
# level, as do the other RGV scripts, and the full suite includes them all into the
# same `Main`.

import Test

module _RGVFigure
include(joinpath(@__DIR__, "..", "..", "benchmarking", "rgv_paired_wilcoxon_figure.jl"))
end

# The minimum a caption needs. `metric_definition` values are vectors because the
# analysis emits one entry per distinct value on the axis.
function _figt_payload(; extra...)
    payload = Dict{String, Any}(
        "metric_definition" => Dict("metric_source" => ["quast"],
            "quast_min_contig" => [500]),
        "seeds" => [42, 123, 456]
    )
    for (k, v) in extra
        payload[string(k)] = v
    end
    return payload
end

Test.@testset "RGV paired-Wilcoxon figure caption" begin
    Test.@testset "a clean run carries no warning glyph" begin
        caption = _RGVFigure.paired_figure_caption(
            _figt_payload(metric_definition_operator_asserted = false,
            metric_definition_override_bound = false))
        Test.@test !occursin("⚠", caption)
        Test.@test occursin("paired Wilcoxon", caption)
        # The definition itself still belongs on the figure: a paired plot of QUAST
        # NGA50 and one of an internal size-ratio proxy look identical.
        Test.@test occursin("metric_source=quast", caption)
        Test.@test occursin("42, 123, 456", caption)
    end

    Test.@testset "an operator-asserted definition is stamped on the figure" begin
        caption = _RGVFigure.paired_figure_caption(
            _figt_payload(metric_definition_operator_asserted = true,
            metric_definition_override_bound = false))
        Test.@test occursin("⚠", caption)
        Test.@test occursin("OPERATOR-ASSERTED", caption)
        Test.@test !occursin("GUARD OVERRIDDEN", caption)
    end

    Test.@testset "a bound guard override is stamped on the figure" begin
        caption = _RGVFigure.paired_figure_caption(
            _figt_payload(metric_definition_operator_asserted = false,
            metric_definition_override_bound = true))
        Test.@test occursin("GUARD OVERRIDDEN", caption)
        Test.@test occursin("NOT validation-grade", caption)
        Test.@test !occursin("OPERATOR-ASSERTED", caption)
    end

    Test.@testset "the two caveats are independent, not exclusive" begin
        # They answer different questions — whether the guard was overridden, and
        # whether the values it compared were observed — so a run in both states
        # must show both. Chaining them would let the weaker one shadow the stronger.
        caption = _RGVFigure.paired_figure_caption(
            _figt_payload(metric_definition_operator_asserted = true,
            metric_definition_override_bound = true))
        Test.@test occursin("GUARD OVERRIDDEN", caption)
        Test.@test occursin("OPERATOR-ASSERTED", caption)
        Test.@test count("⚠", caption) == 2
    end

    Test.@testset "an older results.json renders exactly the clean caption" begin
        # Both flags are read with a `false` default so a payload written before the
        # fields existed does not gain a caveat it has no basis for — and does not
        # raise a KeyError either.
        clean = _RGVFigure.paired_figure_caption(
            _figt_payload(metric_definition_operator_asserted = false,
            metric_definition_override_bound = false))
        legacy = _RGVFigure.paired_figure_caption(_figt_payload())
        Test.@test legacy == clean
        Test.@test !occursin("⚠", legacy)
    end
end
