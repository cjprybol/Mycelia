# Microbiome Visualization

Mycelia provides a unified visualization system for microbiome abundance data, designed to handle datasets ranging from a few samples to 600+ samples with adaptive sizing and automatic view selection.

## Overview

The visualization system addresses common challenges in microbiome data presentation:

- **Large sample counts**: Adaptive sizing prevents label overlap with 300-600+ samples
- **Consistent output**: Publication-quality figures in PNG, SVG, and PDF formats
- **Automatic clustering**: Samples ordered by Bray-Curtis hierarchical clustering
- **Multiple view types**: Barplots, heatmaps, and paginated views auto-selected by sample count

## Quick Start

```julia
import Mycelia
import DataFrames

# Your abundance data in long format
abundance_df = DataFrames.DataFrame(
    sample = ["S1", "S1", "S2", "S2", "S3", "S3"],
    taxon = ["Bacteroides", "Prevotella", "Bacteroides", "Prevotella", "Bacteroides", "Prevotella"],
    relative_abundance = [0.4, 0.6, 0.7, 0.3, 0.5, 0.5]
)

# Generate visualization
results = Mycelia.plot_microbiome_abundance(
    abundance_df,
    sample_col = :sample,
    taxon_col = :taxon,
    abundance_col = :relative_abundance,
    rank = "genus"
)

# Access the barplot
display(results[:barplot])
```

## Configuration

Use `MicrobiomePlotConfig` to customize visualization behavior:

```julia
config = Mycelia.MicrobiomePlotConfig(
    # Figure dimensions
    max_width = 1600,
    min_height = 800,
    pixels_per_sample = 12,
    orientation = :auto,  # :auto, :landscape, :portrait

    # Label sizing
    min_label_fontsize = 6.0,
    max_label_fontsize = 12.0,
    label_rotation = 90.0,

    # Taxa display
    top_n_taxa = 25,
    legend_fontsize = 8.0,

    # Clustering
    sample_ordering = Mycelia.AxisOrdering(
        method = :hclust,
        distance_metric = :braycurtis,
        linkage = :average
    ),

    # Dendrograms
    show_sample_dendrogram = true,
    color_branches_by_cluster = true,
    n_clusters = 5,

    # Large dataset handling
    large_dataset_view = :auto,
    samples_per_page = 150,
    heatmap_threshold = 300,

    # Output
    output_formats = [:png, :svg],
    dpi = 300
)

results = Mycelia.plot_microbiome_abundance(df, config=config, output_dir="figures/")
```

## Ordering Options

Control how samples and taxa are ordered using `AxisOrdering`:

### Hierarchical Clustering (Default)

```julia
ordering = Mycelia.AxisOrdering(
    method = :hclust,
    distance_metric = :braycurtis,  # :braycurtis, :euclidean, :cosine
    linkage = :average              # :single, :complete, :average, :ward
)
```

### Pre-specified Order

```julia
ordering = Mycelia.AxisOrdering(
    method = :preordered,
    preordered_labels = ["Sample_A", "Sample_B", "Sample_C"]
)
```

### Alphabetical

```julia
ordering = Mycelia.AxisOrdering(method = :alphabetical)
```

### Sort by Value

```julia
ordering = Mycelia.AxisOrdering(
    method = :sort,
    sort_by = :mean_abundance  # Or provide a custom function
)
```

## Automatic View Selection

The visualization system automatically selects appropriate views based on sample count:

| Sample Count | Views Generated |
|--------------|-----------------|
| ≤150 | Barplot only |
| 150-300 | Barplot + Heatmap |
| >300 | Heatmap + Paginated barplots |

Override with `large_dataset_view`:

```julia
config = Mycelia.MicrobiomePlotConfig(
    large_dataset_view = :all  # Generate all view types
)
```

Options: `:auto`, `:barplot`, `:heatmap`, `:paginated`, `:all`

## Sizing Strategy

The system uses adaptive sizing to prevent label overlap:

| Samples | Orientation | Label Font | Strategy |
|---------|-------------|------------|----------|
| 1-50 | Landscape | 12pt | Standard |
| 50-150 | Square | 10pt | Compact |
| 150-300 | Portrait | 8pt | Tall |
| 300-600 | Portrait | 6-7pt | Very tall |
| 600+ | Paginated | 8pt | Multi-page |

## Saving Plots

Use `save_plot` to export figures:

```julia
# Save CairoMakie figure
Mycelia.save_plot(results[:barplot], "output/genus_barplot")

# Save StatsPlots figure
Mycelia.save_plot(statsplots_figure, "output/comparison")

# Custom formats and DPI
Mycelia.save_plot(fig, "output/figure", formats=[:png, :svg, :pdf], dpi=600)
```

## Complete Example

```julia
import Mycelia
import DataFrames
import CSV

# Load your abundance data
df = CSV.read("relative_abundance.csv", DataFrames.DataFrame)

# Configure for a large cohort study
config = Mycelia.MicrobiomePlotConfig(
    top_n_taxa = 30,
    orientation = :portrait,
    show_sample_dendrogram = true,
    sample_ordering = Mycelia.AxisOrdering(
        method = :hclust,
        distance_metric = :braycurtis
    ),
    output_formats = [:png, :svg]
)

# Generate all visualizations
results = Mycelia.plot_microbiome_abundance(
    df,
    sample_col = :sample_id,
    taxon_col = :genus,
    abundance_col = :relative_abundance,
    rank = "genus",
    config = config,
    title = "Genus-Level Composition (N=500 samples)",
    output_dir = "figures/genus/"
)

# Access individual views
barplot_fig = results[:barplot]      # CairoMakie.Figure
heatmap_fig = results[:heatmap]      # CairoMakie.Figure (if generated)
pages = results[:paginated]          # Vector{CairoMakie.Figure} (if generated)
prepared_data = results[:data]       # NamedTuple with processed matrices
```

## API Reference

### Main Functions

- [`plot_microbiome_abundance`](@ref) - Main entry point for visualization
- [`save_plot`](@ref) - Save figures in multiple formats

### Configuration Types

- [`MicrobiomePlotConfig`](@ref) - Comprehensive configuration struct
- [`AxisOrdering`](@ref) - Ordering specification for axes

### Utility Functions

- [`adaptive_label_fontsize`](@ref) - Calculate font size based on sample count
- [`calculate_figure_size`](@ref) - Calculate figure dimensions
- [`compute_axis_ordering`](@ref) - Compute sample/taxa ordering
- [`determine_view_types`](@ref) - Select view types based on sample count
- [`calculate_tick_step`](@ref) - Calculate label spacing

## See Also

- [Metagenomic Workflow](metagenomic-workflow.md) - Full workflow including classification
- [API Reference](api/all-functions.md) - Complete function documentation
