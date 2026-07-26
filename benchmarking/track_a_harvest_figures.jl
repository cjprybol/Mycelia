#!/usr/bin/env julia
# Track A baseline harvest — figures for the 2026-07-25 Lovelace/LRC runs.
#
# Produces four figures (SVG + PNG each) from already-computed benchmark outputs.
# This script does NOT run assemblies; it only reads result tables, so it is cheap
# to re-run and safe on a laptop.
#
#   Figure 1  NGA50 coefficient of variation vs the pre-registration's assumed 0.15.
#             -> answers td-bblmi: is CV ~= 0.15 supportable?
#   Figure 2  Per-technology assembly envelope (genome fraction + NGA50/reference).
#             -> explains WHY 8 configuration groups have no computable CV.
#   Figure 3  Cross-host reproducibility (Lovelace vs LRC, identical seeds).
#             -> integrity caveat: which technologies replicate bit-for-bit.
#   Figure 4  Rhizomorph iterative correction vs naive, by error rate x read regime.
#             -> answers td-4e19d.1: where does correction pay?
#
# Usage:
#   julia benchmarking/track_a_harvest_figures.jl \
#       --track-a <track_a_results.tsv> \
#       [--track-a-replicate <second-host track_a_results.tsv>] \
#       [--rgv <rhizomorph_correction_validation_sweep_*.csv>] \
#       [--outdir <plots dir>]
#
# Run with the default environment (CairoMakie/CSV/DataFrames); no Mycelia import
# needed, which keeps it runnable where the Mycelia depot is unavailable.

import CairoMakie
import CSV
import DataFrames
import Statistics
import Printf

# The value the power analysis assumed. Mirrors CV_THRESHOLD in
# benchmarking/track_a_baseline_benchmark.jl — keep the two in sync.
const CV_THRESHOLD = 0.15

# Reference genome lengths (bp), matching ORGANISMS in the Track A driver.
const REFERENCE_LENGTHS = Dict(
    "Lambda" => 48_502,
    "T4" => 168_903,
    "phi29" => 19_282,
    "SARS-CoV-2" => 29_903
)

const TECHNOLOGY_ORDER = ["illumina", "pacbio", "ont"]
const TECHNOLOGY_COLORS = Dict(
    "illumina" => :dodgerblue3,
    "pacbio" => :darkorange3,
    "ont" => :firebrick3
)

# === Argument parsing ===

function arg_value(flag::AbstractString, default = nothing)
    index = findfirst(==(flag), ARGS)
    return (index !== nothing && index < length(ARGS)) ? ARGS[index + 1] : default
end

# === Loading ===

function load_track_a(path::AbstractString)::DataFrames.DataFrame
    return CSV.read(path, DataFrames.DataFrame; delim = '\t')
end

# === Coefficient of variation ===

"""
Per (organism, technology, coverage, decoder_arm) group, compute the NGA50 CV
across seeds. Groups whose mean NGA50 is zero are retained with `computable =
false` so a downstream figure can show them explicitly rather than dropping them
silently — an all-zero group means every seed produced an assembly with no
aligned contig, which is a result, not missing data.
"""
function coefficient_of_variation_table(df::DataFrames.DataFrame)::DataFrames.DataFrame
    grouped = DataFrames.groupby(df, [:organism, :technology, :coverage, :decoder_arm])
    rows = NamedTuple[]
    for group in grouped
        values = Float64.(group.NGA50)
        mean_value = Statistics.mean(values)
        n = length(values)
        computable = n >= 2 && mean_value > 0
        sd_value = n >= 2 ? Statistics.std(values) : NaN
        push!(rows,
            (organism = String(group.organism[1]),
                technology = String(group.technology[1]),
                coverage = Int(group.coverage[1]),
                decoder_arm = String(group.decoder_arm[1]),
                n_seeds = n,
                mean_nga50 = mean_value,
                sd_nga50 = sd_value,
                cv = computable ? sd_value / mean_value : NaN,
                computable = computable))
    end
    return DataFrames.DataFrame(rows)
end

# === Figure 1: CV vs the assumed threshold ===

function figure_cv_vs_threshold(cv_table::DataFrames.DataFrame, outdir::AbstractString)
    # One arm only: the arms are analytically identical on every quality endpoint
    # (contig paths are chosen topologically; quality is applied afterwards), so
    # plotting both would double every point on top of itself.
    single_arm = cv_table[cv_table.decoder_arm .== "kmer", :]
    computable = single_arm[single_arm.computable, :]
    undefined_groups = single_arm[.!single_arm.computable, :]

    figure = CairoMakie.Figure(size = (1450, 620), fontsize = 14)

    axis_a = CairoMakie.Axis(figure[1, 1],
        title = "A. NGA50 CV per configuration, by technology and coverage",
        xlabel = "coverage (x)", ylabel = "NGA50 coefficient of variation (n = 3 seeds)",
        xticks = ([10, 30, 50, 100], ["10", "30", "50", "100"]))

    for (offset, technology) in zip((-4.0, 0.0, 4.0), TECHNOLOGY_ORDER)
        subset = computable[computable.technology .== technology, :]
        isempty(subset) && continue
        CairoMakie.scatter!(axis_a,
            Float64.(subset.coverage) .+ offset, Float64.(subset.cv);
            color = TECHNOLOGY_COLORS[technology], markersize = 13,
            strokewidth = 0.5, strokecolor = :white, label = technology)
    end
    CairoMakie.hlines!(axis_a, [CV_THRESHOLD];
        color = :black, linestyle = :dash, linewidth = 2.5)
    CairoMakie.text!(axis_a, 104, CV_THRESHOLD + 0.008;
        text = "assumed CV = $(CV_THRESHOLD)", align = (:right, :bottom), fontsize = 12)
    CairoMakie.axislegend(axis_a; position = :rt, framevisible = true)

    exceed = sum(computable.cv .> CV_THRESHOLD)
    total = DataFrames.nrow(computable)
    axis_b = CairoMakie.Axis(figure[1, 2],
        title = "B. Distribution of CV by coverage",
        xlabel = "coverage (x)", ylabel = "NGA50 coefficient of variation",
        xticks = ([10, 30, 50, 100], ["10", "30", "50", "100"]))
    for coverage in (10, 30, 50, 100)
        subset = computable[computable.coverage .== coverage, :]
        isempty(subset) && continue
        jitter = range(-3.0, 3.0; length = DataFrames.nrow(subset))
        CairoMakie.scatter!(axis_b,
            fill(Float64(coverage), DataFrames.nrow(subset)) .+ collect(jitter),
            Float64.(subset.cv);
            color = [TECHNOLOGY_COLORS[t] for t in subset.technology], markersize = 12)
        median_cv = Statistics.median(subset.cv)
        CairoMakie.lines!(axis_b, [coverage - 4.5, coverage + 4.5], [median_cv, median_cv];
            color = :black, linewidth = 3)
    end
    CairoMakie.hlines!(axis_b, [CV_THRESHOLD];
        color = :black, linestyle = :dash, linewidth = 2.5)

    CairoMakie.Label(figure[0, 1:2],
        "Track A baseline: the assumed NGA50 CV of $(CV_THRESHOLD) is exceeded in " *
        "$(exceed) of $(total) computable configurations (max $(Printf.@sprintf("%.3f", maximum(computable.cv))))." *
        "\nBlack bars in B are per-coverage medians. " *
        "$(DataFrames.nrow(undefined_groups)) further configurations " *
        "(ONT at 10x and 30x) had NGA50 = 0 for every seed, so no CV is defined — see Figure 2.";
        fontsize = 15, font = :bold, justification = :left, lineheight = 1.15,
        tellwidth = false)

    save_figure(figure, outdir, "track_a_nga50_cv_vs_threshold")
    return (exceed = exceed, total = total, undefined = DataFrames.nrow(undefined_groups))
end

# === Figure 2: per-technology assembly envelope ===

function figure_technology_envelope(df::DataFrames.DataFrame, outdir::AbstractString)
    single_arm = df[df.decoder_arm .== "kmer", :]
    figure = CairoMakie.Figure(size = (1450, 620), fontsize = 14)

    axis_a = CairoMakie.Axis(figure[1, 1],
        title = "A. Reference coverage recovered (QUAST genome fraction)",
        xlabel = "sequencing coverage (x)", ylabel = "genome fraction (%)",
        xticks = ([10, 30, 50, 100], ["10", "30", "50", "100"]))
    axis_b = CairoMakie.Axis(figure[1, 2],
        title = "B. Contiguity: NGA50 as a fraction of the reference length",
        xlabel = "sequencing coverage (x)", ylabel = "NGA50 / reference length",
        xticks = ([10, 30, 50, 100], ["10", "30", "50", "100"]))

    for technology in TECHNOLOGY_ORDER
        subset = single_arm[single_arm.technology .== technology, :]
        isempty(subset) && continue
        coverages = sort(unique(Int.(subset.coverage)))
        fractions = Float64[]
        ratios = Float64[]
        for coverage in coverages
            cells = subset[subset.coverage .== coverage, :]
            push!(fractions, Statistics.mean(Float64.(cells.genome_fraction)))
            push!(ratios,
                Statistics.mean([Float64(row.NGA50) / REFERENCE_LENGTHS[row.organism]
                                 for row in eachrow(cells)]))
        end
        CairoMakie.scatterlines!(axis_a, coverages, fractions;
            color = TECHNOLOGY_COLORS[technology], linewidth = 3, markersize = 12,
            label = technology)
        CairoMakie.scatterlines!(axis_b, coverages, ratios;
            color = TECHNOLOGY_COLORS[technology], linewidth = 3, markersize = 12,
            label = technology)
    end
    CairoMakie.ylims!(axis_a, -3, 105)
    CairoMakie.ylims!(axis_b, -0.03, 1.05)
    CairoMakie.axislegend(axis_a; position = :rb)
    CairoMakie.axislegend(axis_b; position = :lt)

    CairoMakie.Label(figure[0, 1:2],
        "At k = 31, ONT reads recover the reference but never assemble it: genome " *
        "fraction climbs 0.7% -> 98.6% from 10x to 100x,\nyet NGA50 reaches only " *
        "14.9% of the reference length, versus ~83% for Illumina and PacBio. " *
        "PacBio is already at\n99.7% genome fraction and 59% contiguity by 10x. " *
        "Means over 4 organisms x 3 seeds; k-mer arm.";
        fontsize = 15, font = :bold, justification = :left, lineheight = 1.15,
        tellwidth = false)

    save_figure(figure, outdir, "track_a_technology_envelope")
    return nothing
end

# === Figure 3: cross-host reproducibility ===

function figure_cross_host(primary::DataFrames.DataFrame,
        replicate::DataFrames.DataFrame, outdir::AbstractString)
    cell_key(row) = (row.organism, row.technology, row.coverage, row.seed, row.decoder_arm)
    replicate_index = Dict(cell_key(row) => row for row in eachrow(replicate))

    metrics = [:n_reads, :n_contigs, :NGA50, :genome_fraction,
        :duplication_ratio, :largest_contig, :misassemblies]
    shared = 0
    differing = Dict(t => 0 for t in TECHNOLOGY_ORDER)
    totals = Dict(t => 0 for t in TECHNOLOGY_ORDER)
    reads_differ = Dict(t => 0 for t in TECHNOLOGY_ORDER)
    for row in eachrow(primary)
        key = cell_key(row)
        haskey(replicate_index, key) || continue
        other = replicate_index[key]
        shared += 1
        technology = String(row.technology)
        totals[technology] += 1
        if any(row[m] != other[m] for m in metrics)
            differing[technology] += 1
        end
        if row.n_reads != other.n_reads
            reads_differ[technology] += 1
        end
    end

    figure = CairoMakie.Figure(size = (1250, 560), fontsize = 14)
    axis = CairoMakie.Axis(figure[1, 1],
        title = "Cross-host disagreement on identical (organism, coverage, seed, arm) cells",
        xlabel = "sequencing technology", ylabel = "share of shared cells that disagree (%)",
        xticks = (1:length(TECHNOLOGY_ORDER), TECHNOLOGY_ORDER))
    present = [t for t in TECHNOLOGY_ORDER if totals[t] > 0]
    percentages = [100 * differing[t] / totals[t] for t in present]
    CairoMakie.barplot!(axis, 1:length(present), percentages;
        color = [TECHNOLOGY_COLORS[t] for t in present], width = 0.55)
    for (index, technology) in enumerate(present)
        CairoMakie.text!(axis, index, percentages[index] + 1.5;
            text = "$(differing[technology])/$(totals[technology]) cells\n" *
                   "$(reads_differ[technology]) differ in read count",
            align = (:center, :bottom), fontsize = 12)
    end
    CairoMakie.ylims!(axis, 0, 108)

    CairoMakie.Label(figure[0, 1],
        "Same seed, two hosts: Illumina and ONT replicate exactly; PacBio does not." *
        "\nBadread's PacBio arm alone passes --identity, and its bioconda env is " *
        "unpinned, so PacBio\n'seed variance' partly reflects simulator " *
        "nondeterminism rather than assembler behaviour.\n($(shared) shared cells; " *
        "Lovelace 218bfa378 vs LRC 250f7f26.)";
        fontsize = 15, font = :bold, justification = :left, lineheight = 1.15,
        tellwidth = false)

    save_figure(figure, outdir, "track_a_cross_host_reproducibility")
    return (shared = shared, differing = differing, totals = totals)
end

# === Figure 4: iterative correction vs naive ===

function figure_correction_sweep(path::AbstractString, outdir::AbstractString)
    df = CSV.read(path, DataFrames.DataFrame)
    error_rates = sort(unique(Float64.(df.error_rate)))
    regimes = sort(unique(String.(df.regime)))
    arm_colors = Dict("naive" => :grey45, "iterative" => :seagreen4)

    figure = CairoMakie.Figure(size = (1500, 640), fontsize = 14)
    axis_a = CairoMakie.Axis(figure[1, 1],
        title = "A. Largest contig (all cells; internal metric)",
        xlabel = "simulated per-base error rate", ylabel = "largest contig (bp, log scale)",
        yscale = log10, xticks = (1:length(error_rates), string.(error_rates)))
    axis_b = CairoMakie.Axis(figure[1, 2],
        title = "B. Wall-clock cost of correction",
        xlabel = "simulated per-base error rate", ylabel = "runtime (s, log scale)",
        yscale = log10, xticks = (1:length(error_rates), string.(error_rates)))

    # The sweep's own regime names ("short-low-error", "long-high-error") describe the
    # ORIGINAL sweep design and contradict this figure's x-axis, which walks the error
    # rate independently — a series labelled "low-error" plotted at err = 0.10 reads as
    # a mistake. Relabel by what actually varies between the regimes: read length and
    # platform. Unrecognised regimes fall through to their raw name.
    regime_label(regime) =
        if occursin("short", regime)
            "short reads (150 bp, Illumina)"
        elseif occursin("long", regime)
            "long reads (5 kb, ONT)"
        else
            regime
        end

    marker_for(regime) = regime == regimes[1] ? :circle : :utriangle
    for regime in regimes, arm in ("naive", "iterative")

        xs = Float64[]
        largest = Float64[]
        runtimes = Float64[]
        for (index, rate) in enumerate(error_rates)
            cells = df[(df.regime .== regime) .& (df.arm .== arm) .& (Float64.(df.error_rate) .== rate), :]
            isempty(cells) && continue
            push!(xs, index)
            push!(largest, Statistics.mean(Float64.(cells.largest_contig)))
            push!(runtimes, Statistics.mean(Float64.(cells.runtime_s)))
        end
        CairoMakie.scatterlines!(axis_a, xs, largest;
            color = arm_colors[arm], marker = marker_for(regime), markersize = 14,
            linewidth = 3, linestyle = regime == regimes[1] ? :solid : :dash,
            label = "$(arm) / $(regime_label(regime))")
        CairoMakie.scatterlines!(axis_b, xs, runtimes;
            color = arm_colors[arm], marker = marker_for(regime), markersize = 14,
            linewidth = 3, linestyle = regime == regimes[1] ? :solid : :dash,
            label = "$(arm) / $(regime_label(regime))")
    end

    reference_length = Float64(first(df.genome_len))
    CairoMakie.hlines!(axis_a, [reference_length];
        color = :black, linestyle = :dot, linewidth = 2)
    CairoMakie.text!(axis_a, length(error_rates) + 0.02, reference_length;
        text = "reference $(Int(reference_length)) bp",
        align = (:right, :bottom), fontsize = 12)
    CairoMakie.axislegend(axis_a; position = :lb, labelsize = 11)
    CairoMakie.axislegend(axis_b; position = :rb, labelsize = 11)

    quast_ok = sum(.!ismissing.(df.quast_nga50))
    CairoMakie.Label(figure[0, 1:2],
        "Iterative correction pays at 1% error and reverses by 10%: at err = 0.01 it lifts " *
        "NGA50\n12,720 -> 31,181 (short reads) and 28,529 -> 44,523 (long reads), but at " *
        "err = 0.10 it yields MORE and\nSHORTER contigs than naive. Reference-based QUAST " *
        "metrics exist for only $(quast_ok)/$(DataFrames.nrow(df)) cells: above 1% error\nno " *
        "contig clears the min-contig filter, so both panels use internal metrics.";
        fontsize = 15, font = :bold, justification = :left, lineheight = 1.15,
        tellwidth = false)

    save_figure(figure, outdir, "rhizomorph_correction_by_error_regime")
    return nothing
end

# === Output ===

function save_figure(figure, outdir::AbstractString, stem::AbstractString)
    mkpath(outdir)
    png_path = joinpath(outdir, "$(stem).png")
    svg_path = joinpath(outdir, "$(stem).svg")
    CairoMakie.save(png_path, figure)
    CairoMakie.save(svg_path, figure)
    println("  wrote $(png_path)")
    println("  wrote $(svg_path)")
    return (png = png_path, svg = svg_path)
end

# === Entry point ===

function main()
    track_a_path = arg_value("--track-a")
    if track_a_path === nothing
        error("--track-a <track_a_results.tsv> is required")
    end
    outdir = arg_value("--outdir", joinpath(dirname(track_a_path), "plots"))
    replicate_path = arg_value("--track-a-replicate")
    rgv_path = arg_value("--rgv")

    println("=== Track A harvest figures ===")
    track_a = load_track_a(track_a_path)
    println("Track A cells: $(DataFrames.nrow(track_a)) from $(track_a_path)")

    cv_table = coefficient_of_variation_table(track_a)
    summary = figure_cv_vs_threshold(cv_table, outdir)
    println("CV verdict: $(summary.exceed)/$(summary.total) computable groups exceed " *
            "$(CV_THRESHOLD); $(summary.undefined) groups undefined (NGA50 = 0).")

    figure_technology_envelope(track_a, outdir)

    if replicate_path !== nothing
        replicate = load_track_a(replicate_path)
        result = figure_cross_host(track_a, replicate, outdir)
        println("Cross-host: $(result.shared) shared cells; " *
                "differing by technology = $(result.differing)")
    else
        println("(no --track-a-replicate given; skipping cross-host figure)")
    end

    if rgv_path !== nothing
        figure_correction_sweep(rgv_path, outdir)
    else
        println("(no --rgv given; skipping correction-sweep figure)")
    end

    println("Done. Figures in $(outdir)")
    return nothing
end

main()
