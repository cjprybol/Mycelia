import Test
import DataFrames
import Mycelia

function make_batch_qc_dataframe(; include_n50 = true, stage = "after_filtering")
    df = DataFrames.DataFrame(
        stage = [stage, stage, "before_filtering"],
        yield_gb = [2.4, 1.8, 2.9],
        total_reads = [2400, 1800, 2900],
        mean_length = [1200.0, 900.0, 1500.0],
        q20_percent = [95.0, 93.0, 96.0],
        q30_percent = [91.0, 88.0, 92.0]
    )

    if include_n50
        df[!, :n50] = [1500.0, 1100.0, 1700.0]
    end

    return df
end

function make_fastplong_summary(; n50_values = [2000.0, 1800.0])
    return DataFrames.DataFrame(
        total_reads = [1200, 900],
        yield_gb = [2.2, 1.7],
        mean_length = [1500.0, 1700.0],
        n50 = n50_values,
        q20_percent = [94.0, 96.0],
        q30_percent = [88.0, 91.0]
    )
end

Test.@testset "Plotting Utilities" begin
    Test.@testset "Optimal subsequence length" begin
        Test.@test Mycelia.optimal_subsequence_length(error_rate = 0.0) == typemax(Int)
        Test.@test Mycelia.optimal_subsequence_length(error_rate = 1.0) == 1
        Test.@test Mycelia.optimal_subsequence_length(error_rate = 0.01, threshold = 0.99) ==
                   1
    end

    Test.@testset "Jitter" begin
        values = Mycelia.jitter(5.0, 20)
        Test.@test length(values) == 20
        Test.@test maximum(abs.(values .- 5.0)) <= (1.0 / 3.0)
    end

    Test.@testset "Color helpers" begin
        colors = Mycelia.n_maximally_distinguishable_colors(4)
        Test.@test length(colors) == 4

        red = Mycelia.Colors.RGB(1, 0, 0)
        Test.@test Mycelia.merge_colors(red, red) == red
    end

    Test.@testset "Unit conversion helpers" begin
        Test.@test Mycelia.pixels_to_points(12) == 9
        Test.@test Mycelia.points_to_pixels(9) == 12
    end

    Test.@testset "Batch QC plotting" begin
        no_filtered_df = make_batch_qc_dataframe(stage = "before_filtering")
        fig = Test.@test_logs (:warn, r"No 'after_filtering' data found.") Mycelia.plot_batch_qc_distributions(no_filtered_df)
        Test.@test fig isa Mycelia.CairoMakie.Figure

        fig = Mycelia.plot_batch_qc_distributions(make_batch_qc_dataframe())
        Test.@test fig isa Mycelia.CairoMakie.Figure

        fig = Mycelia.plot_batch_qc_distributions(make_batch_qc_dataframe(include_n50 = false))
        Test.@test fig isa Mycelia.CairoMakie.Figure
    end

    Test.@testset "Single-sample QC plotting" begin
        Test.@test Mycelia.visualize_fastplong_single(nothing) isa Mycelia.CairoMakie.Figure

        data_with_drops = (
            summary = make_fastplong_summary(),
            sample_id = "sample-a",
            filtering_stats = Dict(
                "passed_filter_reads" => 900,
                "too_short_reads" => 120,
                "low_quality_reads" => 80
            )
        )
        Test.@test Mycelia.visualize_fastplong_single(data_with_drops) isa Mycelia.CairoMakie.Figure

        data_with_no_drops = (
            summary = make_fastplong_summary(n50_values = [0.0, 0.0]),
            sample_id = "sample-b",
            filtering_stats = Dict("passed_filter_reads" => 900)
        )
        Test.@test Mycelia.visualize_fastplong_single(
            data_with_no_drops;
            title = "Custom title"
        ) isa Mycelia.CairoMakie.Figure

        data_without_filter_stats = (
            summary = make_fastplong_summary(),
            sample_id = "sample-c",
            filtering_stats = Dict{String, Int}()
        )
        Test.@test Mycelia.visualize_fastplong_single(
            data_without_filter_stats
        ) isa Mycelia.CairoMakie.Figure
    end

    Test.@testset "Taxa abundance plotting" begin
        missing_sample_id_df = DataFrames.DataFrame(genus = ["Alpha"])
        Test.@test_throws ErrorException Mycelia.plot_taxa_abundances(
            missing_sample_id_df,
            "genus"
        )

        missing_taxa_level_df = DataFrames.DataFrame(sample_id = ["sample-a"])
        Test.@test_throws ErrorException Mycelia.plot_taxa_abundances(
            missing_taxa_level_df,
            "genus"
        )

        taxa_df = DataFrames.DataFrame(
            sample_id = [
                "sample-b",
                "sample-b",
                "sample-a",
                "sample-a",
                "sample-a",
                "sample-a"
            ],
            genus = [
                "Alpha",
                "Beta",
                "Alpha",
                missing,
                "missing",
                "Gamma"
            ]
        )

        fig, ax, taxa_colors = Mycelia.plot_taxa_abundances(
            taxa_df,
            "genus";
            top_n = 1,
            filter_taxa = ["Gamma"],
            sort_samples = true,
            legend_nbanks = 2
        )
        Test.@test fig isa Mycelia.CairoMakie.Figure
        Test.@test ax isa Mycelia.CairoMakie.Axis
        Test.@test Set(keys(taxa_colors)) == Set(["Alpha", "Other", "Missing"])

        mktempdir() do dir
            save_path = joinpath(dir, "taxa_abundance.png")
            fig, ax, taxa_colors = Mycelia.generate_taxa_abundances_plot(
                taxa_df;
                taxa_level = "genus",
                top_n = 31,
                save_path = save_path
            )

            Test.@test fig isa Mycelia.CairoMakie.Figure
            Test.@test ax isa Mycelia.CairoMakie.Axis
            Test.@test haskey(taxa_colors, "Missing")
            Test.@test isfile(save_path)
        end
    end
end
