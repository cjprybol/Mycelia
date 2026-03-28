# # Tutorial 22: Environmental Metagenome Analysis
#
# This tutorial shows how to turn a public environmental metagenome study into a
# compact, reproducible Mycelia case study. The selected public dataset is the
# Tara Oceans marine metagenome collection (`PRJEB1787`), while the executable
# example uses a bundled lightweight abundance table inspired by the study so it
# can run without network access.
#
# ## Learning Objectives
# - Inspect the curated environmental metagenome dataset catalog
# - Load the bundled Tara Oceans-inspired environmental abundance table
# - Run alpha and beta diversity analysis with Mycelia
# - Save figure outputs for a narrative-ready environmental case study
#
# ## Setup
#
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/22_environmental_metagenome_analysis.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

# Keep GR/StatsPlots in a file-saving mode when the tutorial is run headlessly.
ENV["GKS_WSTYPE"] = get(ENV, "GKS_WSTYPE", "100")

import Mycelia
import DataFrames

println("=== Environmental Metagenome Dataset Catalog ===")
catalog = Mycelia.environmental_metagenome_dataset_catalog()
display(catalog)

selected_dataset = "tara_oceans_surface_prokaryotes"
println("\nSelected dataset id: $selected_dataset")

example_data = Mycelia.environmental_metagenome_example_data(dataset_id = selected_dataset)

println("\n=== Sample Metadata ===")
display(example_data.sample_metadata)

println("\n=== Abundance Table Preview ===")
display(first(example_data.abundance_table, 12))

println("\n=== Running Environmental Metagenome Analysis ===")
analysis = Mycelia.analyze_environmental_metagenome(
    example_data.abundance_table,
    example_data.sample_metadata
)

println("Top taxa across samples:")
display(analysis.top_taxa)

println("\nAlpha diversity:")
display(analysis.alpha_diversity)

println("\nPCoA coordinates:")
display(analysis.beta_diversity.pcoa_df)

output_dir = get(ENV, "MYCELIA_ENV_METAGENOME_OUTPUT_DIR", mktempdir())
println("\nSaving figures to: $output_dir")

figures = Mycelia.generate_environmental_metagenome_figures(
    analysis;
    title_prefix = example_data.dataset_info.title[1],
    output_dir = output_dir,
    top_n_taxa = 8
)

println("\nSaved figure files:")
for (figure_name, paths) in sort(collect(figures.saved_paths), by = first)
    println("  $(figure_name)")
    for (_, filepath) in sort(collect(paths), by = first)
        println("    $filepath")
    end
end

println("\nTutorial complete.")
