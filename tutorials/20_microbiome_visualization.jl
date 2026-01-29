# # Tutorial 20: Microbiome Visualization
#
# This tutorial demonstrates Mycelia's unified microbiome visualization system,
# which produces publication-quality figures with adaptive sizing for datasets
# of any size (from a few samples to 600+).
#
# ## Learning Objectives
# - Create relative abundance visualizations using `plot_microbiome_abundance`
# - Configure visualization options with `MicrobiomePlotConfig`
# - Control sample ordering with `AxisOrdering`
# - Save plots in multiple formats
#
# ## Prerequisites
# - Basic familiarity with DataFrames
# - Microbiome abundance data in long format

# ## Setup

import Mycelia
import DataFrames
import Random

# Set seed for reproducible example data
Random.seed!(42)

# ## Generate Example Data
#
# We'll create synthetic microbiome abundance data to demonstrate the visualization system.

function generate_example_data(n_samples::Int, n_taxa::Int; seed::Int=42)
    Random.seed!(seed)

    # Generate sample names
    samples = ["Sample_$(lpad(i, 3, '0'))" for i in 1:n_samples]

    # Generate taxa names (common gut bacteria genera)
    taxa_names = [
        "Bacteroides", "Prevotella", "Faecalibacterium", "Ruminococcus",
        "Blautia", "Coprococcus", "Dorea", "Roseburia", "Lachnospira",
        "Eubacterium", "Clostridium", "Akkermansia", "Bifidobacterium",
        "Lactobacillus", "Streptococcus", "Enterococcus", "Escherichia",
        "Veillonella", "Dialister", "Megasphaera", "Alistipes", "Parabacteroides",
        "Oscillospira", "Sutterella", "Bilophila"
    ]

    # Use first n_taxa
    taxa = taxa_names[1:min(n_taxa, length(taxa_names))]

    # Generate abundance matrix with some structure
    rows = DataFrames.DataFrame[]

    for sample in samples
        # Random abundances that sum to 1
        abundances = rand(length(taxa))
        abundances ./= sum(abundances)

        for (taxon, abund) in zip(taxa, abundances)
            push!(rows, DataFrames.DataFrame(
                sample = sample,
                taxon = taxon,
                relative_abundance = abund
            ))
        end
    end

    return DataFrames.vcat(rows...)
end

# Generate a small dataset (30 samples, 15 taxa)
small_df = generate_example_data(30, 15)
println("Small dataset: $(DataFrames.nrow(small_df)) rows")
display(first(small_df, 10))

# ## Basic Visualization
#
# The simplest usage requires just a DataFrame with sample, taxon, and abundance columns.

# Basic plot with defaults
results_basic = Mycelia.plot_microbiome_abundance(
    small_df,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus"
)

# The results dictionary contains the generated visualizations
println("\nGenerated views: $(keys(results_basic))")

# Display the barplot
display(results_basic[:barplot])

# ## Customizing Configuration
#
# Use `MicrobiomePlotConfig` to control visualization parameters.

# Create custom configuration
custom_config = Mycelia.MicrobiomePlotConfig(
    top_n_taxa = 10,                    # Show top 10 taxa
    orientation = :landscape,           # Force landscape orientation
    legend_fontsize = 10.0,             # Larger legend text
    show_sample_dendrogram = true,      # Add dendrogram
    output_formats = [:png, :svg]       # Output formats
)

results_custom = Mycelia.plot_microbiome_abundance(
    small_df,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus",
    config = custom_config,
    title = "Custom Configuration Example"
)

display(results_custom[:barplot])

# ## Controlling Sample Order
#
# Use `AxisOrdering` to specify how samples should be ordered.

# Alphabetical ordering
alpha_config = Mycelia.MicrobiomePlotConfig(
    sample_ordering = Mycelia.AxisOrdering(method = :alphabetical),
    top_n_taxa = 10
)

results_alpha = Mycelia.plot_microbiome_abundance(
    small_df,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus",
    config = alpha_config,
    title = "Alphabetical Sample Order"
)

display(results_alpha[:barplot])

# Hierarchical clustering with different distance metrics
hclust_config = Mycelia.MicrobiomePlotConfig(
    sample_ordering = Mycelia.AxisOrdering(
        method = :hclust,
        distance_metric = :braycurtis,  # Bray-Curtis (standard for microbiome)
        linkage = :average
    ),
    show_sample_dendrogram = true,
    top_n_taxa = 10
)

results_hclust = Mycelia.plot_microbiome_abundance(
    small_df,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus",
    config = hclust_config,
    title = "Hierarchical Clustering (Bray-Curtis)"
)

display(results_hclust[:barplot])

# ## Medium Dataset (150 samples)
#
# For larger datasets, the system automatically adjusts sizing and may generate
# multiple view types.

medium_df = generate_example_data(150, 20)
println("\nMedium dataset: $(DataFrames.nrow(medium_df)) rows, $(length(unique(medium_df.sample))) samples")

results_medium = Mycelia.plot_microbiome_abundance(
    medium_df,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus",
    title = "Medium Cohort (150 samples)"
)

println("Generated views: $(keys(results_medium))")

# Display barplot (may be in portrait orientation)
if haskey(results_medium, :barplot)
    display(results_medium[:barplot])
end

# ## Large Dataset (400 samples)
#
# For very large datasets, heatmaps and paginated views are generated automatically.

large_df = generate_example_data(400, 25)
println("\nLarge dataset: $(DataFrames.nrow(large_df)) rows, $(length(unique(large_df.sample))) samples")

large_config = Mycelia.MicrobiomePlotConfig(
    large_dataset_view = :auto,  # Let system decide
    samples_per_page = 100       # For paginated view
)

results_large = Mycelia.plot_microbiome_abundance(
    large_df,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus",
    config = large_config,
    title = "Large Cohort Study"
)

println("Generated views: $(keys(results_large))")

# Display heatmap if generated
if haskey(results_large, :heatmap)
    println("\nHeatmap view:")
    display(results_large[:heatmap])
end

# Display first page of paginated view if generated
if haskey(results_large, :paginated)
    println("\nPaginated view ($(length(results_large[:paginated])) pages):")
    display(results_large[:paginated][1])
end

# ## Saving Plots
#
# Use `save_plot` to export figures in multiple formats.

# Create output directory
output_dir = mktempdir()
println("\nSaving to: $output_dir")

# Save a figure
if haskey(results_basic, :barplot)
    paths = Mycelia.save_plot(
        results_basic[:barplot],
        joinpath(output_dir, "basic_barplot"),
        formats = [:png, :svg],
        dpi = 300
    )
    println("Saved files: $(values(paths))")
end

# ## Automatic Output Directory
#
# Pass `output_dir` to automatically save all generated views.

auto_output = mktempdir()
results_auto = Mycelia.plot_microbiome_abundance(
    small_df,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus",
    output_dir = auto_output  # Automatic saving
)

println("\nAuto-saved files:")
for f in readdir(auto_output)
    println("  $f")
end

# ## Summary
#
# Key takeaways:
# 1. `plot_microbiome_abundance` is the main entry point
# 2. `MicrobiomePlotConfig` controls all visualization parameters
# 3. `AxisOrdering` specifies sample/taxa ordering (hclust, alphabetical, preordered)
# 4. The system automatically selects appropriate views based on sample count
# 5. `save_plot` exports figures in multiple formats
#
# For more details, see the [Microbiome Visualization](../docs/src/microbiome-visualization.md) documentation.
