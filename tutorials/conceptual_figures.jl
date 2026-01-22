# # Conceptual Figures for Mycelia Publication
#
# This script generates conceptual figures for the Mycelia manuscript.
# These figures illustrate core concepts without requiring data generation.
#
# ## Figures Overview
#
# 1. **Figure 1: Graph Type Hierarchy** - The 6-graph progression from n-gram to FASTQ
# 2. **Figure 2: Quality-Aware Assembly** - How per-base quality scores are preserved
# 3. **Figure 3: Workflow Comparison** - Mycelia vs. traditional assembly pipelines
#
# ## Usage
#
# From the Mycelia base directory:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/conceptual_figures.jl", "tutorials/notebooks", execute=false)'
# ```

# ## Setup

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import Plots
import Colors

## Output directory for figures
const FIGURE_DIR = joinpath(@__DIR__, "..", "results", "figures")
mkpath(FIGURE_DIR)

# ## Figure 1: Graph Type Hierarchy
#
# Illustrates the 6-graph type progression in Mycelia's assembly framework:
#
# 1. **N-gram Graph**: Fixed-length character substrings
# 2. **K-mer Graph**: Canonical k-mers (De Bruijn graph)
# 3. **Qualmer Graph**: K-mers with quality score annotations
# 4. **String Graph**: Variable-length overlap graph
# 5. **FASTA Graph**: Sequence records with headers
# 6. **FASTQ Graph**: Sequences with per-base quality scores
#
# The progression represents increasing information content and biological meaning.

function figure1_graph_hierarchy()
    ## TODO: Create diagram showing 6-graph hierarchy
    ## - Boxes for each graph type
    ## - Arrows showing progression/transformation
    ## - Key properties listed for each type
    ## - Color coding: blue (fixed-length) → green (variable-length) → orange (with quality)

    println("Figure 1: Graph Type Hierarchy")
    println("==============================")
    println("N-gram → K-mer → Qualmer → String → FASTA → FASTQ")
    println("")
    println("Fixed-length graphs (information compression):")
    println("  - N-gram: character substrings, language-agnostic")
    println("  - K-mer: canonical DNA/RNA k-mers, De Bruijn structure")
    println("  - Qualmer: k-mers + quality annotations")
    println("")
    println("Variable-length graphs (information preservation):")
    println("  - String: overlap graph, variable-length sequences")
    println("  - FASTA: sequences with identifiers")
    println("  - FASTQ: full sequencing records with quality scores")

    ## Placeholder for actual figure generation
    ## p = Plots.plot(...)
    ## Plots.savefig(p, joinpath(FIGURE_DIR, "figure1_graph_hierarchy.pdf"))
end

figure1_graph_hierarchy()

# ## Figure 2: Quality-Aware Assembly
#
# Illustrates how Mycelia preserves per-base quality information through assembly:
#
# 1. **Input reads** with quality scores (FASTQ)
# 2. **Qualmer graph** nodes carry quality distributions
# 3. **Path selection** uses quality weights
# 4. **Output contigs** retain quality annotations
#
# Traditional assemblers discard quality after initial filtering.
# Mycelia propagates quality through the entire assembly process.

function figure2_quality_aware_assembly()
    ## TODO: Create diagram showing quality preservation
    ## - Left: Input reads with quality heatmap
    ## - Center: Qualmer graph with node quality annotations
    ## - Right: Output contigs with quality tracks
    ## - Arrows showing information flow
    ## - Comparison panel: traditional (quality discarded) vs Mycelia (quality preserved)

    println("Figure 2: Quality-Aware Assembly")
    println("================================")
    println("")
    println("Traditional Pipeline:")
    println("  FASTQ → (discard quality) → K-mer graph → Contigs (no quality)")
    println("")
    println("Mycelia Pipeline:")
    println("  FASTQ → Qualmer graph → Quality-weighted paths → Contigs + quality tracks")
    println("")
    println("Key insight: Quality scores inform which paths through the graph")
    println("represent true biological sequence vs. sequencing errors.")

    ## Placeholder for actual figure generation
end

figure2_quality_aware_assembly()

# ## Figure 3: Workflow Comparison
#
# Compares Mycelia's integrated approach to traditional pipelines:
#
# **Traditional Pipeline:**
# - Separate tools for each step (QC, assembly, annotation, analysis)
# - Data format conversions between steps
# - Quality information lost early
# - Parameter tuning per tool
#
# **Mycelia Approach:**
# - Unified graph framework
# - Native Julia implementation + external tool integration
# - Quality preservation throughout
# - Consistent API and data structures

function figure3_workflow_comparison()
    ## TODO: Create side-by-side workflow diagrams
    ## - Left: Traditional (multiple disconnected tools)
    ## - Right: Mycelia (integrated framework)
    ## - Highlight: data flow, format conversions, quality preservation
    ## - Color coding for tool categories

    println("Figure 3: Workflow Comparison")
    println("=============================")
    println("")
    println("Traditional Pipeline:")
    println("  FastQC → Trimmomatic → SPAdes → Prokka → BLAST → ...")
    println("  (6+ tools, multiple formats, quality discarded)")
    println("")
    println("Mycelia Pipeline:")
    println("  Mycelia.qc() → Mycelia.assemble() → Mycelia.annotate()")
    println("  (unified framework, quality preserved, consistent API)")
    println("")
    println("Mycelia integrates external tools when beneficial while")
    println("maintaining a consistent interface and data model.")

    ## Placeholder for actual figure generation
end

figure3_workflow_comparison()

# ## Export All Figures
#
# Generate publication-quality versions of all figures.

function export_all_figures()
    println("\n" * "="^60)
    println("Generating publication-quality figures...")
    println("="^60)

    ## TODO: Implement actual figure generation with:
    ## - Vector graphics (PDF/SVG) for publication
    ## - Consistent styling (fonts, colors, line weights)
    ## - Appropriate dimensions for journal requirements
    ## - Color-blind friendly palette

    figure1_graph_hierarchy()
    figure2_quality_aware_assembly()
    figure3_workflow_comparison()

    println("\nFigures saved to: $FIGURE_DIR")
    println("Files: figure1_graph_hierarchy.pdf, figure2_quality_aware_assembly.pdf, figure3_workflow_comparison.pdf")
end

## Uncomment to generate all figures:
## export_all_figures()

# ## Summary
#
# This script provides the skeleton for generating conceptual figures.
# The actual visualization code should be added to each function.
#
# Recommended visualization approach:
# - Use Plots.jl or Makie.jl for programmatic figure generation
# - Consider TikZ/LaTeX for precise control over diagram elements
# - Luxor.jl is available in Mycelia for custom graphics
#
# Next steps:
# 1. Decide on specific visual design for each figure
# 2. Implement drawing code in each function
# 3. Add styling consistent with target journal
# 4. Generate final outputs in required formats
