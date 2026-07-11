#!/usr/bin/env julia
# Before/after metrics figure for the Rhizomorph corrector.
#
# Reads the committed real-genome validation CSV (Session A, PR #390) and renders
# a grouped bar chart contrasting the naive (corrector=:none) and scalable
# (corrector=:iterative, strategy=:scalable) arms on phiX174 + lambda.
#
# Two SEPARATE log-scale panels (contigs, N50) — never a dual axis — because the
# two metrics live on different scales. Colors are the Okabe-Ito colorblind-safe
# pair (validated: adjacent-pair CVD dE ~92). Fewer contigs and larger N50 are
# both better; the corrector collapses thousands of fragments into one contig.
#
# Usage (from the Mycelia base directory):
#   julia --project=. docs/plots/rhizomorph_corrector_before_after.jl
#
# Outputs: docs/src/assets/rhizomorph/corrector_before_after.{png,svg}

import CSV
import DataFrames
import CairoMakie

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CSV_PATH = joinpath(REPO_ROOT, "benchmarking", "results",
    "real_data_corrector_validation_20260709_192624.csv")
const OUT_DIR = joinpath(REPO_ROOT, "docs", "src", "assets", "rhizomorph")

# Okabe-Ito, validated colorblind-safe pair. Order is fixed: naive first.
const COLOR_NAIVE = CairoMakie.RGBf(0.835, 0.369, 0.000)     # #D55E00 vermillion
const COLOR_SCALABLE = CairoMakie.RGBf(0.000, 0.447, 0.698)  # #0072B2 blue

isfile(CSV_PATH) || error("Validation CSV not found: $(CSV_PATH)")
df = CSV.read(CSV_PATH, DataFrames.DataFrame)

# Fixed genome + arm order so color follows the entity, never its rank.
genomes = ["phix174", "lambda"]
genome_labels = ["phiX174\n(5.4 kb)", "lambda\n(48.5 kb)"]
arms = ["naive", "scalable"]
arm_labels = ["naive (:none)", "scalable (:iterative)"]
arm_colors = [COLOR_NAIVE, COLOR_SCALABLE]

"""Pull one metric column into a genome x arm matrix in the fixed order."""
function metric_matrix(df, col::Symbol)
    M = Matrix{Float64}(undef, length(genomes), length(arms))
    for (gi, g) in enumerate(genomes), (ai, a) in enumerate(arms)

        row = df[(df.name .== g) .& (df.arm .== a), :]
        DataFrames.nrow(row) == 1 || error("Expected one row for $(g)/$(a)")
        M[gi, ai] = Float64(row[1, col])
    end
    return M
end

contigs = metric_matrix(df, :n_contigs)
n50 = metric_matrix(df, :n50)

"""Grouped, log-scaled bars for one metric with direct value labels."""
function panel!(ax, M; fillto = 0.8)
    xs = Int[];
    ys = Float64[];
    dodge = Int[];
    colors = eltype(arm_colors)[]
    for gi in 1:length(genomes), ai in 1:length(arms)

        push!(xs, gi);
        push!(ys, M[gi, ai]);
        push!(dodge, ai)
        push!(colors, arm_colors[ai])
    end
    CairoMakie.barplot!(ax, xs, ys;
        dodge = dodge, color = colors, fillto = fillto,
        strokecolor = :white, strokewidth = 1.5,
        bar_labels = [string(round(Int, y)) for y in ys],
        label_size = 12, label_offset = 4)
    ax.xticks = (1:length(genomes), genome_labels)
    ax.yscale = log10
    ax.ygridvisible = true
    ax.xgridvisible = false
    ## Headroom so the tallest bars' value labels are not clipped at the top.
    CairoMakie.ylims!(ax, 0.7, maximum(M) * 6)
    return ax
end

CairoMakie.set_theme!(CairoMakie.theme_minimal())
fig = CairoMakie.Figure(size = (860, 460), figure_padding = 18)

CairoMakie.Label(fig[0, 1:2],
    "Rhizomorph corrector: naive vs :iterative/:scalable (real-genome validation)";
    fontsize = 17, font = :bold, padding = (0, 0, 6, 0))

ax1 = CairoMakie.Axis(fig[1, 1];
    title = "Contigs  (fewer is better)", ylabel = "contigs (log scale)",
    titlesize = 14)
panel!(ax1, contigs)

ax2 = CairoMakie.Axis(fig[1, 2];
    title = "N50  (larger is better)", ylabel = "N50, bp (log scale)",
    titlesize = 14)
panel!(ax2, n50)

# Single shared legend (identity is never color-alone).
legend_elems = [CairoMakie.PolyElement(color = c) for c in arm_colors]
CairoMakie.Legend(fig[2, 1:2], legend_elems, arm_labels;
    orientation = :horizontal, framevisible = false, padding = (0, 0, 0, 0))

CairoMakie.Label(fig[3, 1:2],
    "Source: benchmarking/results/real_data_corrector_validation_20260709_192624.csv " *
    "(ART HS25 150 bp PE, 50x, k=21). Corrected arm resolves each genome to a single contig.";
    fontsize = 10, color = :gray40, padding = (0, 0, 0, 0))

CairoMakie.rowsize!(fig.layout, 1, CairoMakie.Relative(0.78))

mkpath(OUT_DIR)
png_path = joinpath(OUT_DIR, "corrector_before_after.png")
svg_path = joinpath(OUT_DIR, "corrector_before_after.svg")
CairoMakie.save(png_path, fig; px_per_unit = 2)
CairoMakie.save(svg_path, fig)
println("Wrote:\n  $(png_path)\n  $(svg_path)")
