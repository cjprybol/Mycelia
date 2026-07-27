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

function require_columns(df::DataFrames.DataFrame, path::AbstractString,
        required::Vector{Symbol})
    absent = [c for c in required if !hasproperty(df, c)]
    isempty(absent) && return nothing
    error("$(path) is missing required column(s): $(join(absent, ", ")). " *
          "Present: $(join(DataFrames.names(df), ", ")). " *
          "This table is likely from an older schema than the figures expect.")
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
    # Group by :k too. A result tree can now hold several k (--k), and pooling them
    # turns a between-SEED CV into a between-K CV — measured at up to ~39x inflation,
    # enough to flip the pre-registration verdict under a y-axis still labelled
    # "n = 3 seeds". Also drop non-ok cells: one crashed seed moved a group's CV from
    # 0.005 to 0.87, and an error row is not a measurement.
    if hasproperty(df, :status)
        df = df[df.status .== "ok", :]
    end
    group_keys = [:organism, :technology, :coverage, :decoder_arm]
    hasproperty(df, :k) && push!(group_keys, :k)
    grouped = DataFrames.groupby(df, group_keys)
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
                k = hasproperty(group, :k) ? Int(group.k[1]) : -1,
                n_seeds = n,
                mean_nga50 = mean_value,
                sd_nga50 = sd_value,
                cv = computable ? sd_value / mean_value : NaN,
                computable = computable))
    end
    return DataFrames.DataFrame(rows)
end

# Report the k actually present in the table. The caption used to open "At k = 31"
# as a literal, which the --k flag turned into a claim the script cannot honour: point
# it at a k-sweep and it would print "At k = 31" over k=11 data.
function describe_k(df::DataFrames.DataFrame)::String
    hasproperty(df, :k) || return ""
    ks = sort(unique(Int.(df.k)))
    length(ks) == 1 && return "k = $(only(ks))."
    return "MIXED k = $(join(ks, ", ")) — series pool incomparable k; facet before reading."
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
    hasproperty(single_arm, :status) && (single_arm = single_arm[single_arm.status .== "ok", :])
    n_organisms = length(unique(single_arm.organism))
    n_seeds = length(unique(single_arm.seed))
    figure = CairoMakie.Figure(size = (1450, 620), fontsize = 14)
    series = Dict{String, NamedTuple}()

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
        series[technology] = (coverages = coverages, fractions = fractions, ratios = ratios)
    end
    CairoMakie.ylims!(axis_a, -3, 105)
    CairoMakie.ylims!(axis_b, -0.03, 1.05)
    CairoMakie.axislegend(axis_a; position = :rb)
    CairoMakie.axislegend(axis_b; position = :lt)

    # Every number below is computed from `series`, not written as a literal. The
    # previous caption hard-coded all eleven of them, so pointing the script at a new
    # table (an ONT k-sweep, say) would redraw the marks while the bold text above them
    # kept asserting the old run's values — and a reader trusts the caption.
    pct(x) = Printf.@sprintf("%.1f%%", x)
    k_note = describe_k(single_arm)
    ont = get(series, "ont", nothing)
    caption = if ont === nothing
        "Per-technology assembly envelope. Means over $(n_organisms) organisms x " *
        "$(n_seeds) seeds; k-mer arm. $(k_note)"
    else
        others = [t for t in ("illumina", "pacbio") if haskey(series, t)]
        plateau = isempty(others) ? "n/a" :
                  pct(100 * Statistics.mean([last(series[t].ratios) for t in others]))
        pac = get(series, "pacbio", nothing)
        pac_note = pac === nothing ? "" :
            " PacBio is already at $(pct(first(pac.fractions))) genome fraction and " *
            "$(pct(100 * first(pac.ratios))) contiguity by $(first(pac.coverages))x."
        "ONT reads recover the reference but never assemble it: genome fraction " *
        "climbs $(pct(first(ont.fractions))) -> $(pct(last(ont.fractions))) from " *
        "$(first(ont.coverages))x to $(last(ont.coverages))x,\nyet NGA50 reaches only " *
        "$(pct(100 * last(ont.ratios))) of the reference length, versus ~$(plateau) " *
        "for Illumina and PacBio.$(pac_note)\nMeans over $(n_organisms) organisms x " *
        "$(n_seeds) seeds; k-mer arm. $(k_note)"
    end
    CairoMakie.Label(figure[0, 1:2], caption;
        fontsize = 15, font = :bold, justification = :left, lineheight = 1.15,
        tellwidth = false)

    save_figure(figure, outdir, "track_a_technology_envelope")
    return nothing
end

# === Figure 3: cross-host reproducibility ===

function figure_cross_host(primary::DataFrames.DataFrame,
        replicate::DataFrames.DataFrame, outdir::AbstractString)
    # k is part of the identity: without it a mixed-k replicate table collapses to a
    # last-write-wins Dict and a k=31 primary is silently compared against a k=19
    # replicate, reporting the difference as cross-host non-reproducibility — exactly
    # the conclusion this figure exists to draw.
    cell_key(row) = (row.organism, row.technology, row.coverage, row.seed,
        row.decoder_arm, hasproperty(row, :k) ? row.k : -1)
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

    figure = CairoMakie.Figure(size = (1250, 620), fontsize = 14)
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

    # Derive the verdict from the computed counts. It used to be the literal string
    # "Illumina and ONT replicate exactly; PacBio does not", which would print
    # unchanged over any pair of inputs — including a pair that disagreed everywhere.
    clean = [t for t in present if differing[t] == 0]
    dirty = [t for t in present if differing[t] > 0]
    agrees(xs) = length(xs) == 1 ? "$(only(xs)) replicates exactly" :
                 "$(join(xs, " and ")) replicate exactly"
    verdict = isempty(dirty) ? "all technologies replicate exactly" :
              isempty(clean) ? "no technology replicates exactly" :
              "$(agrees(clean)); $(join(dirty, " and ")) $(length(dirty) == 1 ? "does" : "do") not"
    organisms = sort(unique(String.(primary.organism)))
    shared_organisms = sort(unique(String.(replicate.organism)))
    scope = length(shared_organisms) < length(organisms) ?
        "\nScope: only $(join(shared_organisms, ", ")) were run on both hosts, so this " *
        "says nothing about $(join(setdiff(organisms, shared_organisms), ", "))." : ""
    CairoMakie.Label(figure[0, 1],
        "Same seed, two runs: $(verdict). ($(shared) shared cells.)" *
        "\nBadread's PacBio arm alone passes --identity and its bioconda env is " *
        "unpinned, which is CONSISTENT WITH\nsimulator nondeterminism — but the " *
        "per-host badread versions were not compared, so the mechanism is inferred," *
        "\nnot established. Note some PacBio cells disagree despite identical read " *
        "counts.$(scope)\nThe two runs may also differ in code revision; this is a " *
        "host+revision comparison unless both are pinned.";
        fontsize = 14, font = :bold, justification = :left, lineheight = 1.15,
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

    headline = Dict{String, Tuple{Float64, Float64}}()
    cell_counts = Dict{Tuple{String, String, Float64}, Int}()
    marker_for(regime) = regime == regimes[1] ? :circle : :utriangle
    for regime in regimes, arm in ("naive", "iterative")

        xs = Float64[]
        largest = Float64[]
        runtimes = Float64[]
        for (index, rate) in enumerate(error_rates)
            cells = df[(df.regime .== regime) .& (df.arm .== arm) .& (Float64.(df.error_rate) .== rate), :]
            isempty(cells) && continue
            push!(xs, index)
            mean_largest = Statistics.mean(Float64.(cells.largest_contig))
            push!(largest, mean_largest)
            push!(runtimes, Statistics.mean(Float64.(cells.runtime_s)))
            cell_counts[(regime, arm, rate)] = DataFrames.nrow(cells)
            if rate == minimum(error_rates)
                prev = get(headline, regime, (0.0, 0.0))
                headline[regime] = arm == "naive" ? (mean_largest, prev[2]) :
                                   (prev[1], mean_largest)
            end
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

    # The only sweep CSV committed to the repo predates the quast_* columns (they
    # landed ~8 h after it was written), so an unguarded `df.quast_nga50` threw
    # ArgumentError here — after both panels were built and before save_figure, so
    # figures 1-3 landed and figure 4 silently did not.
    quast_ok = hasproperty(df, :quast_nga50) ? sum(.!ismissing.(df.quast_nga50)) : 0
    # Derive the headline deltas, and name the metric they come from. The caption used
    # to quote quast_nga50 while BOTH panels plot largest_contig — at err=0.01 short
    # reads panel A shows 13,891 -> 31,195 against a caption saying 12,720 -> 31,181,
    # so a reader taking values off the plot got different numbers with no explanation.
    lowest = minimum(error_rates)
    seeds_per_cell = maximum(values(cell_counts); init = 1)
    replication = seeds_per_cell > 1 ? "$(seeds_per_cell) seeds/cell" :
        "n = 1 SEED PER CELL — no variance estimate; treat the deltas as provisional"
    deltas = String[]
    for regime in regimes
        pair = get(headline, regime, nothing)
        pair === nothing && continue
        push!(deltas, "$(regime_label(regime)) $(Int(round(pair[1]))) -> $(Int(round(pair[2])))")
    end
    CairoMakie.Label(figure[0, 1:2],
        "Iterative correction vs naive at the lowest simulated error rate " *
        "(err = $(lowest)), by largest contig — the\nquantity both panels plot: " *
        "$(join(deltas, "; ")).\nReference-based QUAST metrics exist for only " *
        "$(quast_ok)/$(DataFrames.nrow(df)) cells: above err = $(lowest) no contig " *
        "clears the min-contig\nfilter, so both panels use internal metrics. " *
        "Replication: $(replication).";
        fontsize = 14, font = :bold, justification = :left, lineheight = 1.15,
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
    # Preflight: assert every column the figures and captions depend on exists BEFORE
    # any plotting. Previously a missing column surfaced as an ArgumentError deep in a
    # caption, after three figures had already been written — so the run half-succeeded
    # and the reader had to notice figure 4 was absent.
    isfile(track_a_path) || error("--track-a file not found: $(track_a_path)")
    require_columns(load_track_a(track_a_path), track_a_path,
        [:organism, :technology, :coverage, :seed, :decoder_arm, :NGA50,
            :genome_fraction, :n_contigs])
    if replicate_path !== nothing
        isfile(replicate_path) || error("--track-a-replicate not found: $(replicate_path)")
    end
    if rgv_path !== nothing
        isfile(rgv_path) || error("--rgv file not found: $(rgv_path)")
        require_columns(CSV.read(rgv_path, DataFrames.DataFrame), rgv_path,
            [:regime, :arm, :error_rate, :largest_contig, :runtime_s, :genome_len])
    end
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

# Guard so a test can `include()` this file for its pure helpers
# (coefficient_of_variation_table, describe_k, require_columns) without argument
# parsing or rendering. Matches the convention used across benchmarking/.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
