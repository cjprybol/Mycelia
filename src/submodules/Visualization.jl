"""
Visualization Submodule for Mycelia.jl

Handles plotting and data visualization for biological data.
"""
module Visualization

import AlgebraOfGraphics
import CairoMakie
import Colors
import ColorSchemes
import GraphMakie
import GraphPlot
import Karnak
import Luxor
import Makie
import Plots
import SankeyPlots
import StatsPlots

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# Sequence visualization
export plot_sequence_quality, plot_length_distribution
export plot_gc_content, visualize_alignment

# Assembly visualization
export plot_assembly_stats, contig_length_plot
export assembly_graph_plot, coverage_plot

# Phylogenetic visualization
export plot_phylogeny, plot_tree, tree_heatmap
export dendrogram_plot

# General plotting utilities
export save_plot, multi_panel_plot
export set_plot_theme, customize_colors

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _prepare_plot_data, _apply_plot_styling
export _validate_plot_parameters

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from plotting-and-visualization.jl to here

end # module Visualization