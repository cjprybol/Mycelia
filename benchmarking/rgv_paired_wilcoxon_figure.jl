# Figure for the pre-registered paired-Wilcoxon analysis.
# =======================================================
#
# Reads `results.json` written by `rgv_paired_wilcoxon.jl` and draws, per metric:
#
#   * LEFT  — a paired slope plot: one line per cell from the control arm's value
#             to the treatment arm's value. A paired test's evidence IS the set of
#             within-cell lines, so showing them (rather than two group means) is
#             what makes the figure honest about what was tested.
#   * RIGHT — the paired differences with the median marked and zero as a reference
#             line, which is the quantity the Wilcoxon statistic actually ranks.
#
# Kept in its own script so the analysis stays free of CairoMakie and fast enough
# to unit-test. Emits BOTH SVG and PNG (house convention).
#
# USAGE
# -----
#   julia --project=. benchmarking/rgv_paired_wilcoxon_figure.jl \
#       --json benchmarking/results/rgv_paired_wilcoxon/results.json
#
#   Flags:
#     --json PATH        (required) results.json from rgv_paired_wilcoxon.jl
#     --output-dir PATH  where to write the figure (default: alongside the JSON)
#     --basename NAME    output stem (default "rgv_paired_wilcoxon")

import CairoMakie
import JSON
import Statistics

"""
    paired_figure_caption(payload) -> String

Caption for the paired figure, including any INTEGRITY CAVEAT the payload carries.

Pure (no CairoMakie, no IO) so the caveats can be unit-tested without rendering.

# Why the caveats belong here and not only in `report.md`

The metric definition already belongs on the figure: a paired plot of QUAST NGA50
and a paired plot of an internal size-ratio proxy look identical, so the caption is
what keeps the reader from conflating them (bead td-9p91).

The same argument applies with more force to the two facts that say the definition
is not trustworthy. Under this repo's `1 notebook = 1 figure = 1 slide` model the
SVG/PNG is what reaches a deck or a manuscript, detached from `report.md`. A run
whose guard was overridden, or whose definition was asserted by an operator rather
than observed, previously rendered a caption bit-for-bit identical to a measured
run — so the artifact that travels furthest carried the least disclosure.

Both flags are read with a `false` default so an older `results.json` (written
before the fields existed) renders exactly as it did before.
"""
function paired_figure_caption(payload)
    definition = join(
        [string(k, "=", join(v, "/")) for (k, v) in payload["metric_definition"]], "  ")
    seeds = join(string.(payload["seeds"]), ", ")
    caveats = String[]
    if get(payload, "metric_definition_override_bound", false) === true
        push!(caveats,
            "GUARD OVERRIDDEN — rows span more than one definition — NOT validation-grade")
    end
    if get(payload, "metric_definition_operator_asserted", false) === true
        push!(caveats, "definition OPERATOR-ASSERTED (backfilled), not observed")
    end
    head = "RGV correction-validation sweep — pre-registered paired Wilcoxon\n" *
           "seeds $seeds · definition: $definition"
    return isempty(caveats) ? head : head * "\n⚠ " * join(caveats, "\n⚠ ")
end

"""
    draw_paired_figure(payload; output_dir, basename) -> (svg_path, png_path)

Render the paired slope plot and paired-difference panel for every metric in
`payload` (the parsed `results.json`). Metrics with no usable pairs are still
given a panel, labelled as such, so an empty result is visible rather than absent.
"""
function draw_paired_figure(payload; output_dir::AbstractString,
        basename::AbstractString = "rgv_paired_wilcoxon")
    metrics = payload["metrics"]
    isempty(metrics) && error("results.json contains no metrics")
    treatment = payload["treatment"]
    control = payload["control"]

    nrows = length(metrics)
    fig = CairoMakie.Figure(size = (980, 380 * nrows))
    for (row, m) in enumerate(metrics)
        name = m["metric"]
        pairs = get(m, "pairs", Any[])
        ctrl = Float64[p["control"] for p in pairs]
        trt = Float64[p["treatment"] for p in pairs]
        diffs = Float64[p["difference"] for p in pairs]

        # Only the verdict LABEL goes on the axis; the full pre-registration wording
        # ("SUPPORTED (FDR-adjusted p < 0.05 and ...)") is long enough to overrun
        # into the neighbouring panel, and it is already in report.md verbatim.
        verdict = first(split(m["verdict"], " ("))
        pf = m["pvalue_fdr"]
        subtitle = "n=$(m["n_pairs"]) pairs · FDR p=" *
                   (pf === nothing || (pf isa Real && isnan(pf)) ? "n/a" :
                    string(round(Float64(pf); digits = 5))) * " · " * verdict

        ax1 = CairoMakie.Axis(fig[row, 1];
            title = "$name — paired by cell\n$subtitle",
            ylabel = name, xticks = ([1, 2], [control, treatment]),
            xticklabelrotation = 0.0)
        CairoMakie.xlims!(ax1, 0.7, 2.3)
        if isempty(pairs)
            CairoMakie.text!(ax1, 0.5, 0.5; text = "no usable pairs",
                align = (:center, :center), space = :relative)
        else
            for i in eachindex(ctrl)
                # Colour by direction of change so the paired structure is readable
                # at a glance without implying a group-level claim.
                improving = m["direction"] == "lower" ? trt[i] < ctrl[i] : trt[i] > ctrl[i]
                CairoMakie.lines!(ax1, [1, 2], [ctrl[i], trt[i]];
                    color = improving ? (:steelblue, 0.75) : (:firebrick, 0.75),
                    linewidth = 1.5)
            end
            CairoMakie.scatter!(ax1, fill(1, length(ctrl)), ctrl;
                color = :black, markersize = 7)
            CairoMakie.scatter!(ax1, fill(2, length(trt)), trt;
                color = :black, markersize = 7)
        end

        ax2 = CairoMakie.Axis(fig[row, 2];
            title = "paired differences ($treatment − $control)",
            ylabel = "Δ $name", xticks = ([1], ["pairs"]))
        CairoMakie.xlims!(ax2, 0.5, 1.5)
        CairoMakie.hlines!(ax2, [0.0]; color = :grey40, linestyle = :dash)
        if !isempty(diffs)
            jitter = range(0.9, 1.1; length = max(length(diffs), 2))
            CairoMakie.scatter!(ax2, collect(jitter)[1:length(diffs)], diffs;
                color = :steelblue, markersize = 9)
            med = Statistics.median(diffs)
            CairoMakie.hlines!(ax2, [med]; color = :black, linewidth = 2)
            CairoMakie.text!(ax2, 0.52, med;
                text = "median $(round(med; digits = 2))",
                align = (:left, :bottom), fontsize = 11)
        else
            CairoMakie.text!(ax2, 0.5, 0.5; text = "no usable pairs",
                align = (:center, :center), space = :relative)
        end
    end

    CairoMakie.Label(fig[0, 1:2], paired_figure_caption(payload);
        fontsize = 15, padding = (0, 0, 8, 0))

    mkpath(output_dir)
    svg_path = joinpath(output_dir, basename * ".svg")
    png_path = joinpath(output_dir, basename * ".png")
    CairoMakie.save(svg_path, fig)
    CairoMakie.save(png_path, fig)
    return svg_path, png_path
end

function _fig_arg(flag)
    i = findfirst(==(flag), ARGS)
    return (i !== nothing && i < length(ARGS)) ? ARGS[i + 1] : nothing
end

function main()
    json_path = _fig_arg("--json")
    if json_path === nothing
        println(stderr,
            "ERROR: --json PATH is required (results.json from rgv_paired_wilcoxon.jl).")
        return 2
    end
    isfile(json_path) || (println(stderr, "ERROR: not found: $json_path"); return 1)
    payload = JSON.parsefile(json_path)
    out = something(_fig_arg("--output-dir"), dirname(abspath(json_path)))
    svg,
    png = draw_paired_figure(payload;
        output_dir = out,
        basename = something(_fig_arg("--basename"), "rgv_paired_wilcoxon"))
    println("Wrote:\n  $svg\n  $png")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
