# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/environmental_metagenome_analysis_test.jl")'
# ```

import Test
import Mycelia
import DataFrames

Test.@testset "Environmental Metagenome Analysis" begin
    catalog = Mycelia.environmental_metagenome_dataset_catalog()

    Test.@test DataFrames.nrow(catalog) == 1
    Test.@test catalog.dataset_id[1] == "tara_oceans_surface_prokaryotes"
    Test.@test catalog.public_accession[1] == "PRJEB1787"

    example_data = Mycelia.environmental_metagenome_example_data()

    Test.@test DataFrames.nrow(example_data.dataset_info) == 1
    Test.@test DataFrames.nrow(example_data.sample_metadata) == 6
    Test.@test DataFrames.nrow(example_data.abundance_table) == 48

    prepared = Mycelia.prepare_environmental_abundance_matrix(example_data.abundance_table)
    Test.@test size(prepared.matrix) == (8, 6)
    Test.@test prepared.sample_names == example_data.sample_metadata.sample
    Test.@test prepared.taxon_names[1] == "Pelagibacter"

    analysis = Mycelia.analyze_environmental_metagenome(
        example_data.abundance_table,
        example_data.sample_metadata
    )

    Test.@test size(analysis.abundance_matrix) == (8, 6)
    Test.@test DataFrames.nrow(analysis.alpha_diversity) == 6
    Test.@test "region" in names(analysis.alpha_diversity)
    Test.@test size(analysis.beta_diversity.distance_matrix) == (6, 6)
    Test.@test "PC1" in names(analysis.beta_diversity.pcoa_df)
    Test.@test "PC2" in names(analysis.beta_diversity.pcoa_df)
    Test.@test "region" in names(analysis.beta_diversity.pcoa_df)
    Test.@test DataFrames.nrow(analysis.top_taxa) == 8

    relative = analysis.relative_abundance
    grouped_relative = DataFrames.combine(
        DataFrames.groupby(relative, :sample),
        :relative_abundance => sum => :relative_sum
    )
    Test.@test all(isapprox.(grouped_relative.relative_sum, 1.0; atol = 1e-10))

    output_dir = mktempdir()
    figures = Mycelia.generate_environmental_metagenome_figures(
        analysis;
        title_prefix = "Environmental Metagenome Test",
        output_dir = output_dir
    )

    Test.@test haskey(figures.abundance_views, :barplot)
    Test.@test haskey(figures.saved_paths, :abundance_barplot)
    Test.@test haskey(figures.saved_paths, :alpha_diversity)
    Test.@test haskey(figures.saved_paths, :beta_diversity_pcoa)

    for paths in values(figures.saved_paths)
        for filepath in values(paths)
            Test.@test isfile(filepath)
        end
    end

    case_study = Mycelia.run_environmental_metagenome_case_study(output_dir = mktempdir())
    Test.@test case_study.dataset_info.public_accession[1] == "PRJEB1787"
    Test.@test haskey(case_study.figures.saved_paths, :beta_diversity_pcoa)
end

println("Environmental metagenome analysis tests passed!")
