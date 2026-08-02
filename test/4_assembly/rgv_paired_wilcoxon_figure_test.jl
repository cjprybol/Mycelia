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

    Test.@testset "undeterminable provenance is stamped on the figure" begin
        # The third state `report.md` can report. Carrying two of three would
        # recreate the partial-propagation defect this series has been closing: a
        # fact computed in one file, routed to two of its three consumers, dropped
        # by the third. This is the state EVERY pre-provenance CSV lands in.
        caption = _RGVFigure.paired_figure_caption(
            _figt_payload(metric_definition_provenance_undeterminable = true))
        Test.@test occursin("PROVENANCE UNDETERMINABLE", caption)
        Test.@test occursin("measured-vs-asserted is unknown", caption)
        Test.@test !occursin("OPERATOR-ASSERTED", caption)
        Test.@test !occursin("GUARD OVERRIDDEN", caption)
    end

    Test.@testset "the caveats are independent, not exclusive" begin
        # They answer different questions — whether the guard was overridden,
        # whether the values it compared were observed, and whether that is knowable
        # at all — so a run in several states must show each. Chaining them would let
        # a weaker one shadow a stronger one, which is exactly how the report's own
        # branch chain let a non-binding flag suppress the asserted qualifier.
        caption = _RGVFigure.paired_figure_caption(
            _figt_payload(metric_definition_operator_asserted = true,
            metric_definition_override_bound = true))
        Test.@test occursin("GUARD OVERRIDDEN", caption)
        Test.@test occursin("OPERATOR-ASSERTED", caption)
        Test.@test count("⚠", caption) == 2

        all_three = _RGVFigure.paired_figure_caption(
            _figt_payload(metric_definition_operator_asserted = true,
            metric_definition_override_bound = true,
            metric_definition_provenance_undeterminable = true))
        Test.@test count("⚠", all_three) == 3
    end

    Test.@testset "the definition axes are named in a deterministic order" begin
        # `metric_definition` is a Dict, so iteration order is unspecified: an
        # unsorted comprehension could name the axes differently across two renders
        # of the SAME results.json, making the committed SVG/PNG non-reproducible.
        # Asserting the sorted order pins the contract; building the same payload
        # from two different insertion orders pins that it does not depend on how
        # the Dict happened to be constructed.
        a = Dict{String, Any}(
            "metric_definition" => Dict("metric_source" => ["quast"],
                "quast_min_contig" => [500]),
            "seeds" => [42])
        b = Dict{String, Any}(
            "metric_definition" => Dict("quast_min_contig" => [500],
                "metric_source" => ["quast"]),
            "seeds" => [42])
        Test.@test _RGVFigure.paired_figure_caption(a) ==
                   _RGVFigure.paired_figure_caption(b)
        # Sorted: `metric_source` precedes `quast_min_contig`.
        cap = _RGVFigure.paired_figure_caption(a)
        Test.@test findfirst("metric_source=", cap).start <
                   findfirst("quast_min_contig=", cap).start
    end

    Test.@testset "the caption does not call an exploratory run pre-registered" begin
        # `report.md` opens by stating the run is exploratory: the pre-registration's
        # H1 is Viterbi DP vs greedy, while this sweep compares corrector arms. A
        # figure captioned "pre-registered" contradicts its own report, in the
        # artifact least likely to be read alongside the correction.
        caption = _RGVFigure.paired_figure_caption(_figt_payload())
        Test.@test occursin("EXPLORATORY", caption)
        Test.@test !occursin("pre-registered", caption)
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
