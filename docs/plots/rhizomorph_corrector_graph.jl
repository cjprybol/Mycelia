#!/usr/bin/env julia
# Assembly-redundancy ("bloat") figure for the Rhizomorph corrector.
#
# Complements the contigs/N50 bar chart with a metric it does not show:
# fold-redundancy = total assembled length / reference genome length. A perfect
# assembly reconstructs the genome once (1x); the naive short-read assembly emits
# each strand plus error debris, inflating total length to ~9x the genome, while
# the corrected (:scalable) arm returns to ~1x (a single near-full-length contig).
#
# Data source is the SAME committed CSV as the bar chart and the benchmarks table
# (real_data_corrector_validation_20260709_192624.csv), so all three are
# consistent by construction. Genome lengths are the RefSeq reference sizes.
#
# Usage (from the Mycelia base directory):
#   julia --project=. docs/plots/rhizomorph_corrector_graph.jl
#
# Outputs: docs/src/assets/rhizomorph/corrector_graph.{png,svg}

import CSV
import DataFrames
import CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CSV_PATH = joinpath(REPO_ROOT, "benchmarking", "results",
    "real_data_corrector_validation_20260709_192624.csv")
const OUT_DIR = joinpath(REPO_ROOT, "docs", "src", "assets", "rhizomorph")

const COLOR_NAIVE = CairoMakie.RGBf(0.835, 0.369, 0.000)     # #D55E00
const COLOR_SCALABLE = CairoMakie.RGBf(0.000, 0.447, 0.698)  # #0072B2

# RefSeq reference lengths (bp) for the fold-redundancy denominator.
const GENOME_LEN = Dict("phix174" => 5386, "lambda" => 48502)

isfile(CSV_PATH) || error("Validation CSV not found: $(CSV_PATH)")
df = CSV.read(CSV_PATH, DataFrames.DataFrame)

genomes = ["phix174", "lambda"]
genome_labels = ["phiX174\n(5.4 kb)", "lambda\n(48.5 kb)"]
arms = ["naive", "scalable"]
arm_labels = ["naive (:none)", "corrected (:scalable)"]
arm_colors = [COLOR_NAIVE, COLOR_SCALABLE]

"""Fold-redundancy (total_length / genome_length) for one genome/arm row."""
function fold_redundancy(df, genome, arm)
    row = df[(df.name .== genome) .& (df.arm .== arm), :]
    DataFrames.nrow(row) == 1 || error("Expected one row for $(genome)/$(arm)")
    return Float64(row[1, :total_length]) / GENOME_LEN[genome]
end

fold = [fold_redundancy(df, g, a) for g in genomes, a in arms]  # genomes × arms

# --- Grouped bars ---------------------------------------------------------
CairoMakie.set_theme!(CairoMakie.theme_minimal())
fig = CairoMakie.Figure(size = (720, 470), figure_padding = 18)

CairoMakie.Label(fig[0, 1],
    "Assembly redundancy: total length ÷ genome length";
    fontsize = 16, font = :bold, padding = (0, 0, 6, 0))

ax = CairoMakie.Axis(fig[1, 1];
    ylabel = "fold-redundancy (× genome; 1× = ideal)",
    titlesize = 14)

xs = Int[];
ys = Float64[];
dodge = Int[];
colors = eltype(arm_colors)[]
for gi in 1:length(genomes), ai in 1:length(arms)

    push!(xs, gi);
    push!(ys, fold[gi, ai]);
    push!(dodge, ai)
    push!(colors, arm_colors[ai])
end
CairoMakie.barplot!(ax, xs, ys;
    dodge = dodge, color = colors, strokecolor = :white, strokewidth = 1.5,
    bar_labels = [string(round(y, digits = 1)) * "×" for y in ys],
    label_size = 12, label_offset = 4)

# 1× reference line: a perfect single-copy reconstruction.
CairoMakie.hlines!(ax, [1.0]; color = :gray45, linestyle = :dash, linewidth = 1.3)
CairoMakie.text!(ax, 0.62, 1.0; text = "1× (ideal)", align = (:left, :bottom),
    fontsize = 10, color = :gray40)

ax.xticks = (1:length(genomes), genome_labels)
ax.xgridvisible = false
CairoMakie.ylims!(ax, 0, maximum(fold) * 1.18)

legend_elems = [CairoMakie.PolyElement(color = c) for c in arm_colors]
CairoMakie.Legend(fig[2, 1], legend_elems, arm_labels;
    orientation = :horizontal, framevisible = false)

CairoMakie.Label(fig[3, 1],
    "Naive over-produces ≈9× the genome (both strands + error debris); " *
    "the corrector returns to ≈1×. Source: committed real-data validation CSV.";
    fontsize = 10, color = :gray40, padding = (0, 0, 0, 0))

CairoMakie.rowsize!(fig.layout, 1, CairoMakie.Relative(0.8))

mkpath(OUT_DIR)
png_path = joinpath(OUT_DIR, "corrector_graph.png")
svg_path = joinpath(OUT_DIR, "corrector_graph.svg")
CairoMakie.save(png_path, fig; px_per_unit = 2)
CairoMakie.save(svg_path, fig)
println("Wrote:\n  $(png_path)\n  $(svg_path)")
