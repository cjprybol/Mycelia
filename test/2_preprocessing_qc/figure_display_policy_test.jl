import Test
import Mycelia
import CairoMakie
import StatsPlots

# Locks the figure-display policy that keeps the suite headless: library plotting
# functions must not open windows unless the user opts in, and must never open
# them under CI. All ENV mutations are scoped with `withenv` (`=> nothing` unsets)
# so this test cannot leak display-on into the rest of the suite.
Test.@testset "figure display policy" begin
    Test.@testset "should_display_plots" begin
        # Off by default (no env set).
        withenv("MYCELIA_SHOW_PLOTS" => nothing, "CI" => nothing,
            "GITHUB_ACTIONS" => nothing) do
            Test.@test Mycelia.should_display_plots() == false
        end
        # Opt-in on.
        withenv("MYCELIA_SHOW_PLOTS" => "true", "CI" => nothing,
            "GITHUB_ACTIONS" => nothing) do
            Test.@test Mycelia.should_display_plots() == true
        end
        # CI forces off even with opt-in (both CI and GITHUB_ACTIONS signals).
        withenv("MYCELIA_SHOW_PLOTS" => "true", "CI" => "true") do
            Test.@test Mycelia.should_display_plots() == false
        end
        withenv("MYCELIA_SHOW_PLOTS" => "true", "GITHUB_ACTIONS" => "true",
            "CI" => nothing) do
            Test.@test Mycelia.should_display_plots() == false
        end
        # Any non-"true" value is off.
        withenv("MYCELIA_SHOW_PLOTS" => "false", "CI" => nothing,
            "GITHUB_ACTIONS" => nothing) do
            Test.@test Mycelia.should_display_plots() == false
        end
        # Opt-in is case-insensitive (locks the lowercase() normalization).
        withenv("MYCELIA_SHOW_PLOTS" => "TRUE", "CI" => nothing,
            "GITHUB_ACTIONS" => nothing) do
            Test.@test Mycelia.should_display_plots() == true
        end
    end

    Test.@testset "present_figure" begin
        fig = CairoMakie.Figure()
        ax = CairoMakie.Axis(fig[1, 1])
        CairoMakie.lines!(ax, [0.0, 1.0], [0.0, 1.0])

        # No-op (no display, no save) + returns the figure when nothing configured.
        withenv("MYCELIA_PLOT_ARTIFACTS" => nothing, "MYCELIA_SHOW_PLOTS" => nothing) do
            Test.@test Mycelia.present_figure(fig; name = "x") === fig
        end

        # Saves png + svg when an artifacts dir AND a name are given.
        dir = mktempdir()
        withenv("MYCELIA_PLOT_ARTIFACTS" => dir, "MYCELIA_SHOW_PLOTS" => nothing) do
            Mycelia.present_figure(fig; name = "unit_smoke")
        end
        Test.@test isfile(joinpath(dir, "unit_smoke.png"))
        Test.@test isfile(joinpath(dir, "unit_smoke.svg"))

        # Also exercise the StatsPlots branch of save_plot — most routed sites
        # (assess_dnamer_saturation, analyze_kmer_spectra, sankey, cluster
        # assessment) emit StatsPlots plots, and save_plot handles them on a
        # separate path. Assert non-empty to defeat a doubly-swallowed save
        # failure (save_plot @warns internally; present_figure @warns on top).
        sp = StatsPlots.plot([1, 2, 3], [3, 1, 2])
        dir_sp = mktempdir()
        withenv("MYCELIA_PLOT_ARTIFACTS" => dir_sp, "MYCELIA_SHOW_PLOTS" => nothing) do
            Mycelia.present_figure(sp; name = "statsplots_smoke")
        end
        for ext in ("png", "svg")
            path = joinpath(dir_sp, "statsplots_smoke.$(ext)")
            Test.@test isfile(path)
            Test.@test filesize(path) > 0
        end

        # No artifact written when the dir is set but no name is given.
        dir2 = mktempdir()
        withenv("MYCELIA_PLOT_ARTIFACTS" => dir2, "MYCELIA_SHOW_PLOTS" => nothing) do
            Mycelia.present_figure(fig)
        end
        Test.@test isempty(readdir(dir2))
    end
end
