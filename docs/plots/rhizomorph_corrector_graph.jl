#!/usr/bin/env julia
# Assembly-contiguity ("graph structure") figure for the Rhizomorph corrector.
#
# Runs a small live phiX174 assembly through both arms and renders the resulting
# contig structure as two horizontal tracks: each track is the assembled length
# partitioned into its contigs (white separators mark contig boundaries). The
# naive (:none) track is densely striped — many short, redundant fragments — while
# the scalable (:iterative) track is a single solid block spanning the genome.
# This is the "fragmented graph -> collapsed graph" story in a legible form (a raw
# k-mer node-link diagram at k=21 would be thousands of nodes and unreadable).
#
# Live-run so the fragment structure is real; the headline metrics are also
# reported by the committed CSV / benchmarks page. Uses 30x to run quickly.
#
# Usage (from the Mycelia base directory):
#   julia --project=. docs/plots/rhizomorph_corrector_graph.jl
#
# Outputs: docs/src/assets/rhizomorph/corrector_graph.{png,svg}

import Mycelia
import FASTX
import Random
import CairoMakie

Random.seed!(42)

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT_DIR = joinpath(REPO_ROOT, "docs", "src", "assets", "rhizomorph")

const COLOR_NAIVE = CairoMakie.RGBf(0.835, 0.369, 0.000)     # #D55E00
const COLOR_SCALABLE = CairoMakie.RGBf(0.000, 0.447, 0.698)  # #0072B2

println("Downloading phiX174 + simulating reads (30x Illumina)...")
reference_gz = Mycelia.download_genome_by_accession(accession = "NC_001422.1")
# Decompress to a clean .gz-free path so ART output filenames stay clean.
work_dir = mktempdir()
reference_file = joinpath(work_dir, "phix174.fna")
run(pipeline(`gzip -dc $(reference_gz)`, stdout = reference_file))
ref_len = length(FASTX.sequence(String, first(Mycelia.open_fastx(reference_file))))

art = Mycelia.simulate_illumina_reads(fasta = reference_file, coverage = 30,
    read_length = 150, seqSys = "HS25", paired = true, rndSeed = 42)
read_paths = String[]
for gz in (art.forward_reads, art.reverse_reads)
    fq = replace(gz, r"\.gz$" => "")  # strip only the trailing .gz (reference is .fna.gz)
    isfile(fq) || run(`gunzip -k -f $(gz)`)
    push!(read_paths, fq)
end
reads = reduce(vcat, [collect(Mycelia.open_fastx(p)) for p in read_paths])

println("Assembling both arms...")
naive = Mycelia.Rhizomorph.assemble_genome(reads; k = 21, corrector = :none)
scalable = Mycelia.Rhizomorph.assemble_genome(reads; k = 21,
    corrector = :iterative, strategy = :scalable)

"""Contig lengths sorted descending."""
sorted_lengths(assembly) = sort(length.(assembly.contigs); rev = true)

naive_lengths = sorted_lengths(naive)
scalable_lengths = sorted_lengths(scalable)
println("naive: $(length(naive_lengths)) contigs; scalable: $(length(scalable_lengths)) contigs")

# --- Render two contiguity tracks -----------------------------------------
CairoMakie.set_theme!(CairoMakie.theme_minimal())
fig = CairoMakie.Figure(size = (900, 340), figure_padding = 18)

CairoMakie.Label(fig[0, 1],
    "Assembly contiguity: naive vs :iterative/:scalable (phiX174, live 30x run)";
    fontsize = 16, font = :bold, padding = (0, 0, 8, 0))

ax = CairoMakie.Axis(fig[1, 1];
    xlabel = "assembled length (bp) — each block is one contig, white lines are contig boundaries",
    yticks = ([1, 2], ["scalable\n(:iterative)", "naive\n(:none)"]),
    xlabelsize = 12)
CairoMakie.hidespines!(ax, :t, :r, :l)
ax.ygridvisible = false

"""Draw one horizontal track at row y: contigs laid end-to-end with white gaps."""
function draw_track!(ax, y, lengths, color)
    x = 0.0
    gap = maximum(sum(lengths), init = 1) * 0.0  # solid; separators drawn as strokes
    for L in lengths
        CairoMakie.poly!(ax,
            CairoMakie.Rect(x, y - 0.32, L, 0.64);
            color = color, strokecolor = :white, strokewidth = 1.2)
        x += L
    end
    return x  # total assembled length
end

# Reference genome extent as a faint backdrop on both rows.
for y in (1, 2)
    CairoMakie.poly!(ax, CairoMakie.Rect(0.0, y - 0.34, Float64(ref_len), 0.68);
        color = (:gray70, 0.18), strokecolor = :gray60, strokewidth = 0.8)
end

total_scalable = draw_track!(ax, 1, scalable_lengths, COLOR_SCALABLE)
total_naive = draw_track!(ax, 2, naive_lengths, COLOR_NAIVE)

CairoMakie.text!(ax, total_naive, 2; text = "  $(length(naive_lengths)) contigs",
    align = (:left, :center), fontsize = 12, color = :gray30)
CairoMakie.text!(ax, total_scalable, 1; text = "  $(length(scalable_lengths)) contig",
    align = (:left, :center), fontsize = 12, color = :gray30)

CairoMakie.xlims!(ax, 0, max(total_naive, total_scalable, ref_len) * 1.16)
CairoMakie.ylims!(ax, 0.4, 2.6)

CairoMakie.Label(fig[2, 1],
    "Faint gray band = reference genome length ($(ref_len) bp). The naive assembly " *
    "over-produces short, redundant fragments; the corrector collapses them into one " *
    "near-full-length contig.";
    fontsize = 10, color = :gray40, padding = (0, 0, 0, 0))

mkpath(OUT_DIR)
png_path = joinpath(OUT_DIR, "corrector_graph.png")
svg_path = joinpath(OUT_DIR, "corrector_graph.svg")
CairoMakie.save(png_path, fig; px_per_unit = 2)
CairoMakie.save(svg_path, fig)
println("Wrote:\n  $(png_path)\n  $(svg_path)")
