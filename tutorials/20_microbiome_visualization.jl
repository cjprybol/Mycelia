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

# ## Coverage-Weighted Abundance from Sequencing Data
#
# A common workflow is computing abundance from sequencing coverage combined with
# taxonomic classification. Mycelia provides functions to integrate mosdepth coverage
# data with BLAST or MMseqs2 taxonomy assignments.

# ### Simulated Coverage + Taxonomy Data
#
# Let's create mock data representing a typical metagenomics workflow output.

function generate_coverage_taxonomy_data(n_contigs::Int; seed::Int=42)
    Random.seed!(seed)

    ## Mock mosdepth summary data
    coverage_df = DataFrames.DataFrame(
        chrom = ["contig_$(lpad(i, 4, '0'))" for i in 1:n_contigs],
        length = rand(500:10000, n_contigs),
        bases = rand(500:10000, n_contigs),
        mean = rand(1.0:50.0, n_contigs),
        min = zeros(Int, n_contigs),
        max = rand(10:100, n_contigs)
    )

    ## Mock taxonomy data (some contigs unclassified)
    domains = ["Bacteria", "Viruses", "Archaea", missing]
    phyla = ["Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", missing]
    genera = ["Escherichia", "Bacillus", "Bacteroides", "Streptococcus",
              "Lactobacillus", "Clostridium", "Prevotella", missing]

    taxonomy_df = DataFrames.DataFrame(
        contig_id = coverage_df.chrom,
        domain = rand(domains, n_contigs),
        phylum = rand(phyla, n_contigs),
        genus = rand(genera, n_contigs)
    )

    return coverage_df, taxonomy_df
end

## Generate mock data for 100 contigs
coverage_df, taxonomy_df = generate_coverage_taxonomy_data(100)
println("Coverage data: $(DataFrames.nrow(coverage_df)) contigs")
println("Taxonomy data: $(DataFrames.nrow(taxonomy_df)) assignments")

# ### Merge Coverage with Taxonomy
#
# The `merge_coverage_with_taxonomy` function joins coverage and taxonomy data,
# optionally filtering by minimum coverage or contig length.

merged = Mycelia.merge_coverage_with_taxonomy(
    coverage_df,
    taxonomy_df,
    contig_col = :contig_id,
    min_coverage = 1.0,    ## Require at least 1x mean coverage
    min_length = 500       ## Require at least 500bp
)

println("\nMerged data: $(DataFrames.nrow(merged)) contigs after filtering")
display(first(merged, 5))

# ### Compute Coverage-Weighted Abundance
#
# Calculate relative abundance by summing coverage across taxonomic groups
# and normalizing to proportions.

abundance = Mycelia.compute_coverage_weighted_abundance(
    merged,
    "Sample_001",
    rank = :genus,
    include_unclassified = true
)

println("\nGenus-level abundance:")
display(abundance)

## Verify abundances sum to 1
println("\nTotal relative abundance: $(sum(abundance.relative_abundance))")

# ### Visualize Coverage-Based Abundance

results_coverage = Mycelia.plot_microbiome_abundance(
    abundance,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus",
    title = "Coverage-Weighted Genus Abundance"
)

display(results_coverage[:barplot])

# ### Multi-Sample Coverage Workflow
#
# Process multiple samples in batch for cohort-level visualization.

function generate_multi_sample_coverage_data(n_samples::Int, n_contigs::Int; seed::Int=42)
    Random.seed!(seed)

    all_abundances = DataFrames.DataFrame()

    for s in 1:n_samples
        sample_id = "Sample_$(lpad(s, 3, '0'))"

        ## Each sample has its own coverage/taxonomy data
        cov_df, tax_df = generate_coverage_taxonomy_data(n_contigs, seed=seed+s)

        ## Merge and compute abundance
        merged = Mycelia.merge_coverage_with_taxonomy(
            cov_df, tax_df,
            contig_col = :contig_id,
            min_coverage = 1.0
        )

        abundance = Mycelia.compute_coverage_weighted_abundance(
            merged, sample_id,
            rank = :genus
        )

        all_abundances = DataFrames.vcat(all_abundances, abundance, cols=:union)
    end

    return all_abundances
end

## Generate data for 20 samples
multi_sample_abundance = generate_multi_sample_coverage_data(20, 50)
println("\nMulti-sample data: $(DataFrames.nrow(multi_sample_abundance)) rows")
println("Samples: $(length(unique(multi_sample_abundance.sample)))")
println("Taxa: $(length(unique(multi_sample_abundance.taxon)))")

# ### Visualize Multi-Sample Cohort

cohort_config = Mycelia.MicrobiomePlotConfig(
    top_n_taxa = 10,
    sample_ordering = Mycelia.AxisOrdering(
        method = :hclust,
        distance_metric = :braycurtis
    ),
    show_sample_dendrogram = true
)

results_cohort = Mycelia.plot_microbiome_abundance(
    multi_sample_abundance,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus",
    config = cohort_config,
    title = "Cohort Coverage-Weighted Abundance (N=20)"
)

display(results_cohort[:barplot])

# ### Abundance at Different Taxonomic Ranks
#
# Compare community composition at different taxonomic levels.

for rank in [:domain, :phylum, :genus]
    ## Recompute abundance at each rank
    cov_df, tax_df = generate_coverage_taxonomy_data(100, seed=123)
    merged = Mycelia.merge_coverage_with_taxonomy(cov_df, tax_df, contig_col=:contig_id)

    abundance = Mycelia.compute_coverage_weighted_abundance(
        merged, "Sample_001",
        rank = rank
    )

    println("\n$(titlecase(string(rank)))-level abundance (top 5):")
    display(first(abundance, 5))
end

# ## Summary
#
# Key takeaways:
# 1. `plot_microbiome_abundance` is the main entry point
# 2. `MicrobiomePlotConfig` controls all visualization parameters
# 3. `AxisOrdering` specifies sample/taxa ordering (hclust, alphabetical, preordered)
# 4. The system automatically selects appropriate views based on sample count
# 5. `save_plot` exports figures in multiple formats
# 6. `merge_coverage_with_taxonomy` combines mosdepth coverage with taxonomy data
# 7. `compute_coverage_weighted_abundance` calculates relative abundance from coverage
# 8. Coverage-based abundance works at any taxonomic rank (domain to species)
#
# For more details, see the [Microbiome Visualization](../docs/src/microbiome-visualization.md) documentation.
