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

The same argument applies with more force to every fact that says the definition is
not trustworthy. Under this repo's `1 notebook = 1 figure = 1 slide` model the
SVG/PNG is what reaches a deck or a manuscript, detached from `report.md`. A run
whose guard was overridden, whose definition was asserted by an operator rather than
observed, or whose provenance cannot be determined at all, previously rendered a
caption bit-for-bit identical to a measured run — so the artifact that travels
furthest carried the least disclosure.

ALL THREE states `report.md` can report are rendered here. Carrying two of the three
would recreate the partial-propagation defect this series has spent several rounds
closing: a fact computed in one file, routed to two of its three consumers, and
dropped by the third.

The title says EXPLORATORY for the same reason. This analysis applies the
pre-registration's statistical RULE to a comparison the pre-registration does not
describe — its H1 is Viterbi DP vs greedy, while this sweep compares
`corrector=:none` against `corrector=:iterative` — and `report.md` says so in its
first line. A figure captioned "pre-registered" while its own report calls the run
exploratory is the same conflation, in the artifact least likely to be read
alongside the correction.

Every flag is read with a `false` default, so an older `results.json` written before
the fields existed renders with no caveats, exactly as it did before.
"""
function paired_figure_caption(payload)
    # Sorted by axis name: `metric_definition` is a `Dict`, whose iteration order is
    # unspecified, so an unsorted comprehension can name the axes in a different
    # order for two renders of the SAME results.json. That makes the SVG/PNG
    # non-reproducible and every figure diff noisy — which matters here more than
    # usual, because under `1 notebook = 1 figure = 1 slide` these files are
    # committed artifacts that get diffed and re-rendered.
    md = payload["metric_definition"]
    definition = join(
        [string(k, "=", join(md[k], "/")) for k in sort(collect(keys(md)))], "  ")
    seeds = join(string.(payload["seeds"]), ", ")
    caveats = String[]
    # The override caveat lives in this caption, which `draw_paired_figure` stamps
    # on the figure itself. That is the disclosure that travels with the SVG/PNG.
    # A second red band repeating "GUARD OVERRIDDEN — NOT validation-grade" would
    # say the same thing again (td-thq2 item 2); the caption is the record.
    if get(payload, "metric_definition_override_bound", false) === true
        push!(caveats,
            "GUARD OVERRIDDEN — rows span more than one definition — NOT validation-grade")
    end
    if get(payload, "metric_definition_operator_asserted", false) === true
        push!(caveats, "definition OPERATOR-ASSERTED (backfilled), not observed")
    end
    if get(payload, "metric_definition_provenance_undeterminable", false) === true
        push!(caveats,
            "definition PROVENANCE UNDETERMINABLE — measured-vs-asserted is unknown")
    end
    head = "RGV correction-validation sweep — paired Wilcoxon (EXPLORATORY)\n" *
           "seeds $seeds · definition: $definition"
    return isempty(caveats) ? head : head * "\n⚠ " * join(caveats, "\n⚠ ")
end

"""
    paired_line_color(direction, control, treatment)

Colour of one paired slope.

`"none"` metrics get no improvement claim — `decide` refuses one — so the line
is neutral grey. Treating a higher treatment value as "improving" would assert
exactly the directional claim the verdict withholds, including on
`quast_duplication_ratio`.
"""
function paired_line_color(direction, control, treatment)
    if string(direction) == "none"
        return (:grey40, 0.85)
    end
    improving = string(direction) == "lower" ? treatment < control : treatment > control
    return improving ? (:steelblue, 0.75) : (:firebrick, 0.75)
end

"""
    _pair_float(x) -> Float64

A per-pair number. `null` is an error: the analysis emits finite pair values,
and a null here used to become `MethodError` inside `Float64[...]`.
"""
function _pair_float(x)
    x isa Real ||
        error("paired figure expected a finite number, got $(repr(x)); " *
              "null is not a pair value")
    y = Float64(x)
    isfinite(y) || error("paired figure expected a finite number, got $y")
    return y
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
        ctrl = Float64[_pair_float(p["control"]) for p in pairs]
        trt = Float64[_pair_float(p["treatment"]) for p in pairs]
        diffs = Float64[_pair_float(p["difference"]) for p in pairs]

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
                # at a glance. `"none"` stays grey: that is not an improvement claim.
                CairoMakie.lines!(ax1, [1, 2], [ctrl[i], trt[i]];
                    color = paired_line_color(m["direction"], ctrl[i], trt[i]),
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

const _FIG_KNOWN_FLAGS = ("--json", "--output-dir", "--basename")
const _FIG_VALUE_FLAGS = ("--json", "--output-dir", "--basename")

"""
    _cli_problem(args, known, value_flags) -> Union{String, Nothing}

`nothing` when every `--flag` is known and every value-taking flag is followed
by a value. Otherwise a message. Unknown flags used to be ignored, and a value
slot used to accept the next flag (`--output-dir --basename` named a directory
`--basename`).
"""
function _cli_problem(args, known, value_flags)
    known_set = Set(String.(known))
    value_set = Set(String.(value_flags))
    for a in args
        if startswith(a, "--") && !(a in known_set)
            return "unrecognized flag `$a`"
        end
    end
    for (i, a) in enumerate(args)
        a in value_set || continue
        if i == length(args) || startswith(args[i + 1], "--")
            nxt = i == length(args) ? "end of arguments" : "`$(args[i + 1])`"
            return "`$a` requires a value, not $nxt"
        end
    end
    return nothing
end

function _fig_arg(flag, args = ARGS)
    i = findfirst(==(flag), args)
    i === nothing && return nothing
    if i == length(args) || startswith(args[i + 1], "--")
        error("`$flag` requires a value, not " *
              (i == length(args) ? "end of arguments" : "`$(args[i + 1])`"))
    end
    return args[i + 1]
end

function main()
    problem = _cli_problem(ARGS, _FIG_KNOWN_FLAGS, _FIG_VALUE_FLAGS)
    if problem !== nothing
        println(stderr, "ERROR: $problem")
        return 2
    end
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
