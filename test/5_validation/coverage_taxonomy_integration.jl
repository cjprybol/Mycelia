# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/5_validation/coverage_taxonomy_integration.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/5_validation/coverage_taxonomy_integration.jl", "test/5_validation", execute=false)'
# ```

import Test
import Mycelia
import DataFrames

Test.@testset "Coverage-Weighted Taxonomy Integration" begin
    ## Setup: Create mock data for testing
    test_dir = mktempdir()

    Test.@testset "merge_coverage_with_taxonomy" begin
        ## Create mock mosdepth summary data
        coverage_df = DataFrames.DataFrame(
            chrom = ["contig_001", "contig_002", "contig_003", "contig_004", "total"],
            length = [5000, 3000, 4000, 2000, 14000],
            bases = [5000, 3000, 4000, 2000, 14000],
            mean = [15.2, 8.5, 22.1, 3.0, 12.5],
            min = [0, 0, 0, 0, 0],
            max = [45, 30, 55, 10, 55]
        )

        ## Create mock taxonomy data
        taxonomy_df = DataFrames.DataFrame(
            contig_id = ["contig_001", "contig_002", "contig_003"],  ## contig_004 has no taxonomy
            domain = ["Bacteria", "Viruses", "Bacteria"],
            kingdom = [missing, missing, missing],
            phylum = ["Proteobacteria", missing, "Firmicutes"],
            class = ["Gammaproteobacteria", missing, "Bacilli"],
            order = ["Enterobacterales", missing, "Lactobacillales"],
            family = ["Enterobacteriaceae", "Siphoviridae", "Streptococcaceae"],
            genus = ["Escherichia", "Lambdavirus", "Streptococcus"],
            species = ["Escherichia coli", "Lambda phage", "Streptococcus pyogenes"]
        )

        Test.@testset "Basic merge with default parameters" begin
            merged = Mycelia.merge_coverage_with_taxonomy(
                coverage_df,
                taxonomy_df,
                contig_col = :contig_id
            )

            Test.@test merged isa DataFrames.DataFrame
            ## Should exclude "total" row and only have contigs
            Test.@test DataFrames.nrow(merged) == 4  ## 4 contigs, 3 with taxonomy
            Test.@test "chrom" in names(merged)
            Test.@test "mean" in names(merged)
        end

        Test.@testset "Merge with minimum coverage filter" begin
            merged = Mycelia.merge_coverage_with_taxonomy(
                coverage_df,
                taxonomy_df,
                contig_col = :contig_id,
                min_coverage = 5.0
            )

            ## contig_004 has mean=3.0, should be filtered out
            Test.@test DataFrames.nrow(merged) == 3
            Test.@test !("contig_004" in merged.chrom)
        end

        Test.@testset "Merge with minimum length filter" begin
            merged = Mycelia.merge_coverage_with_taxonomy(
                coverage_df,
                taxonomy_df,
                contig_col = :contig_id,
                min_length = 3000
            )

            ## contig_004 has length=2000, should be filtered out
            Test.@test DataFrames.nrow(merged) == 3
            Test.@test !("contig_004" in merged.chrom)
        end

        Test.@testset "Merge with both filters" begin
            merged = Mycelia.merge_coverage_with_taxonomy(
                coverage_df,
                taxonomy_df,
                contig_col = :contig_id,
                min_coverage = 10.0,
                min_length = 3000
            )

            ## Only contig_001 (mean=15.2, length=5000) and contig_003 (mean=22.1, length=4000) pass
            Test.@test DataFrames.nrow(merged) == 2
            Test.@test all(merged.mean .>= 10.0)
            Test.@test all(merged.length .>= 3000)
        end
    end

    Test.@testset "compute_coverage_weighted_abundance" begin
        ## Create mock merged data
        merged_df = DataFrames.DataFrame(
            chrom = ["contig_001", "contig_002", "contig_003", "contig_004"],
            mean = [15.2, 8.5, 22.1, 5.0],
            domain = ["Bacteria", "Viruses", "Bacteria", missing],
            phylum = ["Proteobacteria", missing, "Firmicutes", missing],
            genus = ["Escherichia", "Lambdavirus", "Streptococcus", missing]
        )

        Test.@testset "Abundance at genus level" begin
            abundance = Mycelia.compute_coverage_weighted_abundance(
                merged_df,
                "Sample_001",
                rank = :genus
            )

            Test.@test abundance isa DataFrames.DataFrame
            Test.@test "sample" in names(abundance)
            Test.@test "taxon" in names(abundance)
            Test.@test "relative_abundance" in names(abundance)
            Test.@test "total_coverage" in names(abundance)
            Test.@test "n_contigs" in names(abundance)

            ## All samples should be "Sample_001"
            Test.@test all(abundance.sample .== "Sample_001")

            ## Relative abundances should sum to 1
            Test.@test isapprox(sum(abundance.relative_abundance), 1.0, atol=1e-10)
        end

        Test.@testset "Abundance at domain level" begin
            abundance = Mycelia.compute_coverage_weighted_abundance(
                merged_df,
                "Sample_001",
                rank = :domain
            )

            ## Should have Bacteria, Viruses, and Unclassified
            Test.@test DataFrames.nrow(abundance) == 3
            Test.@test "Bacteria" in abundance.taxon
            Test.@test "Viruses" in abundance.taxon
            Test.@test "Unclassified" in abundance.taxon

            ## Check coverage sums
            bacteria_row = abundance[abundance.taxon .== "Bacteria", :]
            Test.@test bacteria_row.total_coverage[1] == 15.2 + 22.1  ## contig_001 + contig_003
            Test.@test bacteria_row.n_contigs[1] == 2
        end

        Test.@testset "Abundance excluding unclassified" begin
            abundance = Mycelia.compute_coverage_weighted_abundance(
                merged_df,
                "Sample_001",
                rank = :domain,
                include_unclassified = false
            )

            ## Should only have Bacteria and Viruses
            Test.@test DataFrames.nrow(abundance) == 2
            Test.@test "Bacteria" in abundance.taxon
            Test.@test "Viruses" in abundance.taxon
            Test.@test !("Unclassified" in abundance.taxon)

            ## Relative abundances should still sum to 1 (among classified only)
            Test.@test isapprox(sum(abundance.relative_abundance), 1.0, atol=1e-10)
        end

        Test.@testset "Custom unclassified label" begin
            abundance = Mycelia.compute_coverage_weighted_abundance(
                merged_df,
                "Sample_001",
                rank = :genus,
                include_unclassified = true,
                unclassified_label = "Unknown"
            )

            Test.@test "Unknown" in abundance.taxon
            Test.@test !("Unclassified" in abundance.taxon)
        end

        Test.@testset "Abundance is sorted by relative abundance" begin
            abundance = Mycelia.compute_coverage_weighted_abundance(
                merged_df,
                "Sample_001",
                rank = :genus
            )

            ## Should be sorted in descending order
            for i in 1:(DataFrames.nrow(abundance)-1)
                Test.@test abundance.relative_abundance[i] >= abundance.relative_abundance[i+1]
            end
        end
    end

    Test.@testset "Contig column detection" begin
        ## Test _detect_contig_column helper

        Test.@testset "Detect :contig_id" begin
            df = DataFrames.DataFrame(contig_id = ["a", "b"], value = [1, 2])
            col = Mycelia._detect_contig_column(df)
            Test.@test col == :contig_id
        end

        Test.@testset "Detect :query_id" begin
            df = DataFrames.DataFrame(query_id = ["a", "b"], value = [1, 2])
            col = Mycelia._detect_contig_column(df)
            Test.@test col == :query_id
        end

        Test.@testset "Detect Symbol(\"query id\")" begin
            df = DataFrames.DataFrame()
            df[!, "query id"] = ["a", "b"]
            df[!, "value"] = [1, 2]
            col = Mycelia._detect_contig_column(df)
            Test.@test col == Symbol("query id")
        end

        Test.@testset "Fallback to first column" begin
            df = DataFrames.DataFrame(unknown_col = ["a", "b"], value = [1, 2])
            col = Mycelia._detect_contig_column(df)
            Test.@test col == :unknown_col
        end
    end

    Test.@testset "MicrobiomePlotConfig and AxisOrdering structs" begin
        ## Test that configuration structs are correctly defined

        Test.@testset "MicrobiomePlotConfig defaults" begin
            config = Mycelia.MicrobiomePlotConfig()

            Test.@test config.max_width == 1600
            Test.@test config.min_height == 800
            Test.@test config.pixels_per_sample == 12
            Test.@test config.orientation == :auto
            Test.@test config.min_label_fontsize == 6.0
            Test.@test config.max_label_fontsize == 12.0
            Test.@test config.top_n_taxa == 25
            Test.@test config.large_dataset_view == :auto
        end

        Test.@testset "MicrobiomePlotConfig custom values" begin
            config = Mycelia.MicrobiomePlotConfig(
                max_width = 2000,
                top_n_taxa = 15,
                orientation = :portrait
            )

            Test.@test config.max_width == 2000
            Test.@test config.top_n_taxa == 15
            Test.@test config.orientation == :portrait
        end

        Test.@testset "AxisOrdering defaults" begin
            ordering = Mycelia.AxisOrdering()

            Test.@test ordering.method == :hclust
            Test.@test ordering.distance_metric == :braycurtis
            Test.@test ordering.linkage == :average
            Test.@test ordering.preordered_labels === nothing
        end

        Test.@testset "AxisOrdering alphabetical" begin
            ordering = Mycelia.AxisOrdering(method = :alphabetical)
            Test.@test ordering.method == :alphabetical
        end

        Test.@testset "AxisOrdering preordered" begin
            labels = ["Sample_C", "Sample_A", "Sample_B"]
            ordering = Mycelia.AxisOrdering(
                method = :preordered,
                preordered_labels = labels
            )
            Test.@test ordering.method == :preordered
            Test.@test ordering.preordered_labels == labels
        end
    end

    Test.@testset "Adaptive sizing functions" begin
        Test.@testset "adaptive_label_fontsize" begin
            ## Few samples should get max fontsize
            fontsize_small = Mycelia.adaptive_label_fontsize(30)
            Test.@test fontsize_small == 12.0

            ## Many samples should get min fontsize
            fontsize_large = Mycelia.adaptive_label_fontsize(600)
            Test.@test fontsize_large == 6.0

            ## Medium samples should be interpolated
            fontsize_medium = Mycelia.adaptive_label_fontsize(275)
            Test.@test 6.0 < fontsize_medium < 12.0
        end

        Test.@testset "calculate_figure_size" begin
            config = Mycelia.MicrobiomePlotConfig()

            ## Small sample count (landscape)
            width, height = Mycelia.calculate_figure_size(50, 25, config=config)
            Test.@test width >= 800
            Test.@test height >= config.min_height

            ## Large sample count (portrait)
            width_large, height_large = Mycelia.calculate_figure_size(300, 25, config=config)
            Test.@test width_large <= config.max_width
            Test.@test height_large > height  ## Should be taller for more samples
        end

        Test.@testset "determine_view_types" begin
            config = Mycelia.MicrobiomePlotConfig()

            ## Few samples - just barplot
            views_small = Mycelia.determine_view_types(50, config=config)
            Test.@test :barplot in views_small
            Test.@test length(views_small) == 1

            ## Medium samples - barplot and heatmap
            views_medium = Mycelia.determine_view_types(200, config=config)
            Test.@test :barplot in views_medium
            Test.@test :heatmap in views_medium

            ## Large samples - heatmap and paginated
            views_large = Mycelia.determine_view_types(500, config=config)
            Test.@test :heatmap in views_large
            Test.@test :paginated in views_large
        end

        Test.@testset "calculate_tick_step" begin
            ## Few samples - show all
            step_small = Mycelia.calculate_tick_step(30, config=Mycelia.MicrobiomePlotConfig())
            Test.@test step_small == 1

            ## Many samples - skip some
            step_large = Mycelia.calculate_tick_step(300, config=Mycelia.MicrobiomePlotConfig())
            Test.@test step_large > 1
        end
    end

    Test.@testset "Integration: Full abundance workflow with mock data" begin
        ## Create a realistic mock dataset
        n_samples = 5
        n_contigs_per_sample = 20

        all_abundances = DataFrames.DataFrame()

        for s in 1:n_samples
            sample_id = "Sample_$(lpad(s, 3, '0'))"

            ## Create mock coverage data
            coverage_df = DataFrames.DataFrame(
                chrom = ["contig_$(lpad(i, 3, '0'))" for i in 1:n_contigs_per_sample],
                length = rand(1000:10000, n_contigs_per_sample),
                bases = rand(1000:10000, n_contigs_per_sample),
                mean = rand(1.0:50.0, n_contigs_per_sample),
                min = zeros(Int, n_contigs_per_sample),
                max = rand(10:100, n_contigs_per_sample)
            )

            ## Create mock taxonomy data (some contigs unclassified)
            taxa = ["Bacteria", "Viruses", "Archaea", missing]
            genera = ["Escherichia", "Bacillus", "Streptococcus", "Lambdavirus", "Methanobacterium", missing]

            taxonomy_df = DataFrames.DataFrame(
                contig_id = coverage_df.chrom,
                domain = rand(taxa, n_contigs_per_sample),
                genus = rand(genera, n_contigs_per_sample)
            )

            ## Merge and compute abundance
            merged = Mycelia.merge_coverage_with_taxonomy(
                coverage_df,
                taxonomy_df,
                contig_col = :contig_id,
                min_coverage = 1.0
            )

            abundance = Mycelia.compute_coverage_weighted_abundance(
                merged,
                sample_id,
                rank = :domain
            )

            all_abundances = DataFrames.vcat(all_abundances, abundance, cols = :union)
        end

        Test.@test DataFrames.nrow(all_abundances) > 0
        Test.@test length(unique(all_abundances.sample)) == n_samples

        ## Check that each sample's abundances sum to 1
        for sample_df in DataFrames.groupby(all_abundances, :sample)
            total_abundance = sum(sample_df.relative_abundance)
            Test.@test isapprox(total_abundance, 1.0, atol=1e-10)
        end

        ## Test that the data format is compatible with plot_microbiome_abundance
        ## (don't actually generate plot in tests to avoid graphics dependencies)
        Test.@test "sample" in names(all_abundances)
        Test.@test "taxon" in names(all_abundances)
        Test.@test "relative_abundance" in names(all_abundances)
    end

    ## Cleanup
    rm(test_dir, recursive=true, force=true)
end
