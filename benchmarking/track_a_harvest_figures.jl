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

# Colour lookup that tolerates a technology the palette does not know, rather than
# throwing a KeyError partway through rendering.
technology_color(t) = get(TECHNOLOGY_COLORS, String(t), :grey40)

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

# Keep only cells whose metrics are real measurements. Two independent ways a row can
# carry zeros that are NOT measurements, and `status` only detects the first:
#
#   1. status != "ok" — the cell errored, or assembled nothing.
#   2. status == "ok" but QUAST never scored it. On a QUAST exception the harness
#      substitutes empty_metrics() — NGA50, genome_fraction and largest_contig all
#      0.0 — and then derives status from n_contigs ALONE, so a cell QUAST never
#      scored is written as "ok" with a full row of zeros. Filtering on
#      status alone misses it completely, i.e. through the same door the status
#      filter was added to close. Not hypothetical: QUAST was failing on one host
#      (conda ProcessExited(4)), and the merged matrix mixes metric sources by host.
#
# The detector is the implication track_a_merge_hosts.jl's `quast_evidence` already
# uses: n_contigs > 0 AND largest_contig == 0 means QUAST did not score this cell
# (QUAST exits non-zero when nothing clears min_contig, so both routes land in the
# same empty_metrics() fallback). Both columns are already required, so
# the data was loaded and unused.
#
# Filtered HERE rather than by changing `status` at the source, because the 432
# checkpoints already on disk record "ok" — a source-side fix would only protect
# future runs and would leave every existing tree misreported.
#
# Missing-safe throughout: `df.status .== "ok"` yields Union{Missing,Bool} when any
# status is missing and DataFrames throws on indexing with that mask, so a single
# missing field would abort the figure rather than drop a row. coalesce() makes
# missing mean "not a measurement". Shared by the CV table and the envelope figure
# so the two cannot drift apart.
#
# The mask is built by `measurement_mask` and reported by `exclusion_counts`, both of
# which read the SAME predicates, so what the figures drop and what the run reports
# dropping cannot drift apart. Reporting is mandatory rather than nice-to-have: the
# filter is one-directional on this file's headline verdict, and not through the CV
# magnitude — through the DENOMINATOR. Dropping a seed shrinks its group's n, and at
# n == 1 the group leaves `computable` entirely, so a group that would have FAILED the
# threshold is removed from "E of T exceed" instead of counted against it. The sibling
# power analysis (benchmarking/track_a_baseline_benchmark.jl) already prints this
# accounting into power_analysis_summary.md next to the counts; the two consumers of
# the same predicate must not disagree about how much of the matrix went unscored.
function measurement_mask(df::DataFrames.DataFrame)
    n_rows = DataFrames.nrow(df)
    # No status column: nothing is classifiable as a non-measurement, so keep the frame
    # whole rather than emptying it.
    hasproperty(df, :status) ||
        return (keep = trues(n_rows), non_ok = falses(n_rows), unscored = falses(n_rows))
    ok = coalesce.(df.status .== "ok", false)
    unscored = (hasproperty(df, :n_contigs) && hasproperty(df, :largest_contig)) ?
               coalesce.(df.n_contigs .> 0, false) .&
               coalesce.(df.largest_contig .== 0, false) :
               falses(n_rows)
    # Attribute each dropped row to exactly ONE reason, so the two counts sum to the
    # drop total and a reader can add them up: a non-ok row that is ALSO unscored is
    # reported under non-ok, because `status` is the coarser and earlier failure.
    return (keep = ok .& .!unscored, non_ok = .!ok, unscored = ok .& unscored)
end

function ok_cells(df::DataFrames.DataFrame)::DataFrames.DataFrame
    return df[measurement_mask(df).keep, :]
end

"""
    exclusion_counts(df) -> (n_total, n_kept, n_non_ok, n_unscored)

How many rows `ok_cells` drops, split by reason. The two reasons are different
failures and are counted separately: `non_ok` is a cell that errored or assembled
nothing, `unscored` is a cell the harness recorded as "ok" but QUAST never scored.
Conflating them would report a QUAST outage as an assembly failure.

Pure and separately callable so the count can be pinned by a test and printed by
`main` — the drop count previously appeared nowhere at all, so stdout said
"Track A cells: 288" while the figures were computed over however many survived.
"""
function exclusion_counts(df::DataFrames.DataFrame)
    mask = measurement_mask(df)
    return (n_total = DataFrames.nrow(df), n_kept = count(mask.keep),
        n_non_ok = count(mask.non_ok), n_unscored = count(mask.unscored))
end

"""
    exclusion_clause(counts) -> String

The caption sentence disclosing the excluded rows. Emitted on EVERY render, including
the no-exclusions case — a clause that appears only on a shortfall leaves the reader
unable to tell "nothing was excluded" from "exclusions were not accounted for", and
this file already carries a fix for the opposite defect (a clause that always warns).
"""
function exclusion_clause(counts)
    dropped = counts.n_non_ok + counts.n_unscored
    dropped == 0 &&
        return "All $(counts.n_total) cells are measurements; none excluded."
    reasons = String[]
    counts.n_non_ok > 0 && push!(reasons, "$(counts.n_non_ok) status != ok")
    counts.n_unscored > 0 &&
        push!(reasons,
            "$(counts.n_unscored) QUAST-unscored (n_contigs > 0, largest_contig = 0)")
    return "$(dropped) of $(counts.n_total) cells excluded as non-measurements " *
           "($(join(reasons, "; "))); every count above is over the " *
           "$(counts.n_kept) surviving cells."
end

"""
    cv_undefined_reason(n_seeds, mean_nga50) -> String

Why a configuration group has no computable CV — "too_few_seeds", "all_zero", or
"computable". Kept as one function so the CV table's column and figure 1's caption
classify identically.

`n_seeds < 2` is checked FIRST, and the precedence is load-bearing. With a single
surviving seed there is no between-seed spread to estimate, so the CV is undefined
whatever the value is; and that group may have reached n == 1 precisely BECAUSE
`ok_cells` dropped a seed, whose NGA50 is then unknown to this table. Classifying it
as "all_zero" would assert something the surviving data cannot support — the defect
this split exists to remove, where the caption stated "NGA50 = 0 for every seed" over
groups whose surviving seed was a fully assembled genome.
"""
function cv_undefined_reason(n_seeds::Integer, mean_nga50::Real)
    n_seeds < 2 ? "too_few_seeds" : mean_nga50 <= 0 ? "all_zero" : "computable"
end

# "(ont at 10x and 30x)" — name a bucket's configurations from its rows. Empty for an
# empty bucket, so a caption never carries a dangling pair of parentheses.
function configuration_label(rows::DataFrames.DataFrame)::String
    DataFrames.nrow(rows) == 0 && return ""
    techs = sort(unique(String.(rows.technology)))
    covs = sort(unique(Int.(rows.coverage)))
    return "($(join(techs, ", ")) at $(join(string.(covs), "x and "))x)"
end

"""
    split_undefined_groups(undefined_groups) -> (all_zero, too_few_seeds)

Partition the non-computable CV groups by WHY they are non-computable, instead of
asserting one reason for all of them.

The figure-1 caption stated "had NGA50 = 0 for every seed" over every non-computable
group. That is false for any group that reached `n_seeds < 2` because `ok_cells`
dropped a seed. Reproduced against the 2026-07-25 Lovelace table with QUAST made to
fail on a third of its cells: of the 20 groups the old caption covered, 12 had a
SURVIVING seed with NGA50 > 0 — including Lambda/pacbio/100x at NGA50 = 48502, a
fully assembled genome — so the figure asserted the opposite of its own data.
Classified through `cv_undefined_reason`, so this split and the CV table's
`undefined_reason` column cannot disagree.
"""
function split_undefined_groups(undefined_groups::DataFrames.DataFrame)
    reasons = String[cv_undefined_reason(row.n_seeds, row.mean_nga50)
                     for row in eachrow(undefined_groups)]
    return (undefined_groups[reasons .== "all_zero", :],
        undefined_groups[reasons .== "too_few_seeds", :])
end

"""
    undefined_groups_clause(all_zero, too_few_seeds) -> String

The figure-1 caption sentence about configurations with no computable CV. One clause
per reason, each naming only the configurations it actually covers, and each count
derived from its own bucket.

Pure and separately callable for the same reason `quast_coverage_clause` is: the
sentence it replaces was an inline literal asserting a cause, and no assertion could
reach it. "NGA50 = 0 for every SURVIVING seed" is the honest form of that claim — the
CV table only ever saw the seeds `ok_cells` kept, and `exclusion_clause` states how
many it did not see.
"""
function undefined_groups_clause(all_zero::DataFrames.DataFrame,
        too_few_seeds::DataFrames.DataFrame)
    n_zero = DataFrames.nrow(all_zero)
    n_few = DataFrames.nrow(too_few_seeds)
    n_zero + n_few == 0 && return "No configuration is left without a computable CV."
    zero_why = "had NGA50 = 0 for every surviving seed"
    few_why = "have fewer than 2 surviving seeds, so no CV exists whatever their NGA50"
    # With one reason in play the total and the bucket count are the same number, so
    # state it once. Only the mixed case needs a total to add up to.
    n_few == 0 && return "$(n_zero) further configurations " *
           "$(configuration_label(all_zero)) $(zero_why), so no CV is defined — " *
           "see Figure 2."
    n_zero == 0 && return "$(n_few) further configurations " *
           "$(configuration_label(too_few_seeds)) $(few_why) — see Figure 2."
    # One reason per LINE in the mixed case. CairoMakie sizes a Label block by its
    # longest line and centres the block, so a single over-long line pushes the whole
    # caption past the figure edge and clips the other lines with it — verified by
    # rendering the two-reason case, where the unwrapped form ran ~295 characters and
    # truncated all three lines including the headline verdict.
    return "$(n_zero + n_few) further configurations have no defined CV — see Figure 2:" *
           "\n  $(n_zero) $(zero_why) $(configuration_label(all_zero))" *
           "\n  $(n_few) $(few_why) $(configuration_label(too_few_seeds))"
end

"""
    quast_coverage_clause(quast_ok, n_rows, metric_sources) -> String

The figure-4 caption sentence describing how much of the sweep QUAST actually
scored. Pure and separately callable so it can be pinned by a test — the previous
version was an inline literal, and no assertion could reach it.

Two behaviours worth stating, because the literal got both wrong:

  * It is emitted as a SHORTFALL warning only when there is a shortfall. The old
    text printed "exists for only N/N cells: above err = X no contig clears the
    min-contig filter" even when QUAST had scored every row — a caption that always
    warns, which is the defect class this file's tests exist to close.
  * It does not assert WHY. The old clause blamed the min-contig filter; the
    committed sweep records `metric_source = "internal:quast-failed"` for all eight
    unscored rows (QUAST itself failed), and four of those eight have
    `largest_contig > 500`, so they would have cleared the ceiling #439 settled on.
    The sweep carries no `quast_min_contig` column, so no caption can name a
    threshold honestly. Report `metric_source`, which is a fact in the data.
"""
function quast_coverage_clause(quast_ok::Integer, n_rows::Integer,
        metric_sources::AbstractVector{<:AbstractString} = String[])
    if quast_ok >= n_rows
        return "Reference-based QUAST metrics cover all $(quast_ok) cells; both " *
               "panels plot the internal metric so the arms stay comparable."
    end
    others = sort(unique(filter(!=("quast"), metric_sources)))
    why = isempty(others) ? "" : " (metric_source: $(join(others, ", ")))"
    return "Reference-based QUAST metrics exist for only $(quast_ok)/$(n_rows) " *
           "cells$(why), so both panels use internal metrics."
end

# === Coefficient of variation ===

"""
Per (organism, technology, coverage, decoder_arm) group, compute the NGA50 CV
across seeds. Groups whose mean NGA50 is zero are retained with `computable =
false` so a downstream figure can show them explicitly rather than dropping them
silently — an all-zero group means every seed produced an assembly with no
aligned contig, which is a result, not missing data.

`computable` is false for two different reasons, and the group carries which one in
`undefined_reason` rather than leaving a consumer to guess: fewer than 2 SURVIVING
seeds ("too_few_seeds"), or >= 2 surviving seeds that are all zero ("all_zero").
Seeds counted here are those `ok_cells` kept, so a group can reach n == 1 because a
seed was excluded — see `exclusion_counts` for that accounting.
"""
function coefficient_of_variation_table(df::DataFrames.DataFrame)::DataFrames.DataFrame
    # Group by :k too. A result tree can now hold several k (--k), and pooling them
    # turns a between-SEED CV into a between-K CV — measured at up to ~39x inflation,
    # enough to flip the pre-registration verdict under a y-axis still labelled
    # "n = 3 seeds". Also drop non-ok cells: one crashed seed moved a group's CV from
    # 0.005 to 0.87, and an error row is not a measurement.
    df = ok_cells(df)
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
                computable = computable,
                # Carry the REASON, not just the fact, so no consumer has to assert
                # one. "no computable CV" has two disjoint causes with opposite
                # readings — an all-zero group is a result about the assembler, a
                # single-seed group is a result about the harvest (possibly about
                # what ok_cells removed) — and the caption used to state the first
                # over both.
                undefined_reason = cv_undefined_reason(n, mean_value)))
    end
    # An empty `rows` yields a 0x0 DataFrame with no :decoder_arm, so the caller throws a
    # bare ArgumentError from inside plotting instead of saying what is wrong. This is
    # reachable — ok_cells drops every row on a harvest where nothing succeeded — so
    # fail here, naming the cause.
    isempty(rows) && error("no cells with status == \"ok\" to compute a CV over. " *
          "Either the run failed everywhere, or the table's status column is not " *
          "populated as expected.")
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

# `exclusions` is the `exclusion_counts` of the table the CV table was built FROM; the
# CV table itself retains no trace of the rows ok_cells removed, so it cannot be
# recovered here. `main` always supplies it. Defaulting to `nothing` (caption omits the
# clause) rather than to a zero-count tuple is deliberate: a fabricated "none excluded"
# would be the same class of always-reassuring caption this file's other fixes remove.
function figure_cv_vs_threshold(cv_table::DataFrames.DataFrame, outdir::AbstractString;
        exclusions = nothing)
    # One arm only: the arms are analytically identical on every quality endpoint
    # (contig paths are chosen topologically; quality is applied afterwards), so
    # plotting both would double every point on top of itself.
    single_arm = cv_table[cv_table.decoder_arm .== "kmer", :]
    computable = single_arm[single_arm.computable, :]
    undefined_groups = single_arm[.!single_arm.computable, :]

    # Derive the two remaining literals: the seed count on the y-axis, and WHICH
    # configurations are undefined. Both were hard-coded ("n = 3 seeds", "ONT at 10x and
    # 30x") in a file whose whole convention is to compute its captions, so a different
    # matrix would have mislabelled the axis and named the wrong cells.
    # The CV-table guard fires only when EVERY row is dropped, but this figure then
    # narrows to the kmer arm, and that subset can be empty while the table is not — a
    # qualmer-only harvest, say. The reductions below then threw
    # "reducing over an empty collection ... consider supplying `init`", from inside
    # figure construction, which is the uninformative-abort class the guard exists to
    # remove. Note the advice in that message is a trap here: supplying `init` to
    # `minimum` is exactly what produced the caption bug fixed in the previous commit.
    # An explicit emptiness check is the correct instrument.
    isempty(single_arm) && error("no kmer-arm rows in the CV table, so figure 1 has " *
          "nothing to plot. Present arms: " *
          "$(join(sort(unique(String.(cv_table.decoder_arm))), ", ")).")
    seed_label = let n = sort(unique(single_arm.n_seeds))
        length(n) == 1 ? "$(only(n)) seeds" : "$(minimum(n))-$(maximum(n)) seeds"
    end
    all_zero_groups, too_few_groups = split_undefined_groups(undefined_groups)

    figure = CairoMakie.Figure(size = (1450, 620), fontsize = 14)

    axis_a = CairoMakie.Axis(figure[1, 1],
        title = "A. NGA50 CV per configuration, by technology and coverage",
        xlabel = "coverage (x)",
        ylabel = "NGA50 coefficient of variation (n = $(seed_label))",
        xticks = ([10, 30, 50, 100], ["10", "30", "50", "100"]))

    plotted_any = false
    for (offset, technology) in zip((-4.0, 0.0, 4.0), TECHNOLOGY_ORDER)
        subset = computable[computable.technology .== technology, :]
        isempty(subset) && continue
        plotted_any = true
        CairoMakie.scatter!(axis_a,
            Float64.(subset.coverage) .+ offset, Float64.(subset.cv);
            color = technology_color(technology), markersize = 13,
            strokewidth = 0.5, strokecolor = :white, label = technology)
    end
    CairoMakie.hlines!(axis_a, [CV_THRESHOLD];
        color = :black, linestyle = :dash, linewidth = 2.5)
    CairoMakie.text!(axis_a, 104, CV_THRESHOLD + 0.008;
        text = "assumed CV = $(CV_THRESHOLD)", align = (:right, :bottom), fontsize = 12)
    # axislegend errors outright when no plot in the axis carries a label, which happens
    # whenever every technology subset is empty. A missing legend is a cosmetic loss; a
    # thrown figure is not.
    plotted_any && CairoMakie.axislegend(axis_a; position = :rt, framevisible = true)

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
            color = [technology_color(t) for t in subset.technology], markersize = 12)
        median_cv = Statistics.median(subset.cv)
        CairoMakie.lines!(axis_b, [coverage - 4.5, coverage + 4.5], [median_cv, median_cv];
            color = :black, linewidth = 3)
    end
    CairoMakie.hlines!(axis_b, [CV_THRESHOLD];
        color = :black, linestyle = :dash, linewidth = 2.5)

    undefined_sentence = undefined_groups_clause(all_zero_groups, too_few_groups)
    # Never silently omit the exclusion accounting when it is available; when it is not,
    # say so rather than implying nothing was dropped.
    exclusion_sentence = exclusions === nothing ?
                         "Exclusion accounting not supplied to this figure." :
                         exclusion_clause(exclusions)

    CairoMakie.Label(figure[0, 1:2],
        "Track A baseline: the assumed NGA50 CV of $(CV_THRESHOLD) is exceeded in " *
        "$(exceed) of $(total) computable configurations$(total == 0 ? "" : " (max " *
            Printf.@sprintf("%.3f", maximum(computable.cv)) * ")")." *
        "\nBlack bars in B are per-coverage medians. $(undefined_sentence)" *
        "\n$(exclusion_sentence)";
        fontsize = 15, font = :bold, justification = :left, lineheight = 1.15,
        tellwidth = false)

    save_figure(figure, outdir, "track_a_nga50_cv_vs_threshold")
    return (exceed = exceed, total = total,
        undefined = DataFrames.nrow(undefined_groups),
        undefined_all_zero = DataFrames.nrow(all_zero_groups),
        undefined_too_few_seeds = DataFrames.nrow(too_few_groups))
end

# === Figure 2: per-technology assembly envelope ===

function figure_technology_envelope(df::DataFrames.DataFrame, outdir::AbstractString)
    single_arm = df[df.decoder_arm .== "kmer", :]
    single_arm = ok_cells(single_arm)
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

    plotted_any = false
    for technology in TECHNOLOGY_ORDER
        subset = single_arm[single_arm.technology .== technology, :]
        isempty(subset) && continue
        coverages = sort(unique(Int.(subset.coverage)))
        fractions = Float64[]
        ratios = Float64[]
        for coverage in coverages
            cells = subset[subset.coverage .== coverage, :]
            push!(fractions, Statistics.mean(Float64.(cells.genome_fraction)))
            # Skip organisms with no known reference length rather than throwing a
            # KeyError halfway through rendering. A new organism in the table is a
            # reason to extend REFERENCE_LENGTHS, not to lose the whole figure.
            known = [row
                     for row in eachrow(cells) if haskey(REFERENCE_LENGTHS, row.organism)]
            unknown = setdiff(unique(String.(cells.organism)), collect(keys(REFERENCE_LENGTHS)))
            isempty(unknown) ||
                @warn "no reference length known; excluded from contiguity panel" unknown
            push!(ratios,
                isempty(known) ? NaN :
                Statistics.mean([Float64(row.NGA50) / REFERENCE_LENGTHS[row.organism]
                                 for row in known]))
        end
        CairoMakie.scatterlines!(axis_a, coverages, fractions;
            color = technology_color(technology), linewidth = 3, markersize = 12,
            label = technology)
        CairoMakie.scatterlines!(axis_b, coverages, ratios;
            color = technology_color(technology), linewidth = 3, markersize = 12,
            label = technology)
        series[technology] = (coverages = coverages, fractions = fractions, ratios = ratios)
        plotted_any = true
    end
    CairoMakie.ylims!(axis_a, -3, 105)
    CairoMakie.ylims!(axis_b, -0.03, 1.05)
    # Same guard as figure 1. Guarding only figure 1 was worse than guarding neither:
    # figure 1 then survives a technology outside TECHNOLOGY_ORDER while figure 2 dies
    # on it, so the run half-succeeds — exactly the outcome the preflight exists to
    # prevent. Reproduced by renaming every technology to "nanopore".
    if plotted_any
        CairoMakie.axislegend(axis_a; position = :rb)
        CairoMakie.axislegend(axis_b; position = :lt)
    end

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
        color = [technology_color(t) for t in present], width = 0.55)
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
    # `shared == 0` must never read as agreement. An empty comparison makes `dirty`
    # empty, which would otherwise print the strongest possible claim — "all
    # technologies replicate exactly" — from no evidence at all.
    verdict = isempty(dirty) ?
              "all $(length(clean)) compared technologies replicate exactly" :
              isempty(clean) ? "no technology replicates exactly" :
              "$(agrees(clean)); $(join(dirty, " and ")) $(length(dirty) == 1 ? "does" : "do") not"
    # Gate on MEMBERSHIP, not count. Comparing lengths suppressed the caveat exactly
    # where it matters most: two tables with equally many but DISJOINT organisms passed
    # the length test while covering nothing in common. `replicate_organisms` is also
    # the honest name — it is the replicate's list, not an intersection.
    organisms = sort(unique(String.(primary.organism)))
    replicate_organisms = sort(unique(String.(replicate.organism)))
    untested = setdiff(organisms, replicate_organisms)
    covered = intersect(organisms, replicate_organisms)
    # `covered` cannot be empty here: no shared organism forces shared == 0, which routes
    # the caption to the NO SHARED CELLS branch below that never interpolates `scope`.
    # So this branch only ever runs with a non-empty `covered`.
    scope = isempty(untested) ? "" :
            "\nScope: only $(join(covered, ", ")) " *
            "$(length(covered) == 1 ? "appears" : "appear") in both tables, so this says " *
            "nothing about $(join(untested, ", "))."
    # On the zero-shared path the mechanism and scope prose is not merely irrelevant, it
    # reads as commentary on a result that does not exist. Emit a short self-contained
    # caption instead, and never let an empty comparison print an agreement verdict.
    caption = if shared == 0
        "NO SHARED CELLS — the two tables have no cell in common,\nso this figure " *
        "establishes NOTHING about reproducibility.\nCheck that both tables cover the " *
        "same organisms,\ncoverages, seeds, decoder arms and k."
    else
        "Same seed, two runs: $(verdict). ($(shared) shared cells.)" *
        "\nBadread's PacBio arm alone passes --identity and its bioconda env is " *
        "unpinned, which is CONSISTENT WITH\nsimulator nondeterminism — but the " *
        "per-host badread versions were not compared, so the mechanism is inferred," *
        "\nnot established. Note some PacBio cells disagree despite identical read " *
        "counts.$(scope)\nThe two runs may also differ in code revision; this is a " *
        "host+revision comparison unless both are pinned."
    end
    CairoMakie.Label(figure[0, 1], caption;
        fontsize = 14, font = :bold, justification = :left, lineheight = 1.15,
        tellwidth = false)

    save_figure(figure, outdir, "track_a_cross_host_reproducibility")
    return (shared = shared, differing = differing, totals = totals)
end

# === Figure 4: iterative correction vs naive ===

function figure_correction_sweep(df::DataFrames.DataFrame, outdir::AbstractString)
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

    # The element type must admit `nothing`: a regime missing one arm stores a sentinel
    # rather than a fabricated 0.0, and the concrete Tuple{Float64,Float64} rejected it
    # with a convert MethodError on the canonical render.
    headline = Dict{String, Tuple{Union{Float64, Nothing}, Union{Float64, Nothing}}}()
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
                # `nothing`, not (0.0, 0.0): a 0.0 default is indistinguishable from a
                # measurement, and a missing naive arm rendered as "naive produced a
                # 0 bp contig" — the most favourable possible reading of the iterative
                # arm, printed in bold over the figure. Same class as the `init` bug.
                prev = get(headline, regime, (nothing, nothing))
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

    # Only draw the reference line when the sweep has ONE reference. Taking row 1's
    # length and labelling it "reference N bp" across the whole panel silently asserts
    # the wrong length for every other genome in a mixed sweep.
    genome_lengths = unique(Float64.(df.genome_len))
    reference_length = length(genome_lengths) == 1 ? only(genome_lengths) : nothing
    if isempty(genome_lengths)
        @warn "sweep table has no rows; omitting the reference line"
    elseif length(genome_lengths) > 1
        @warn "sweep spans multiple reference lengths; omitting the reference line" genome_lengths
    end
    if reference_length !== nothing
        CairoMakie.hlines!(axis_a, [reference_length];
            color = :black, linestyle = :dot, linewidth = 2)
        CairoMakie.text!(axis_a, length(error_rates) + 0.02, reference_length;
            text = "reference $(Int(reference_length)) bp",
            align = (:right, :bottom), fontsize = 12)
    end
    # Same exposure as figures 1 and 2: an empty sweep table means no labelled series.
    if !isempty(regimes) && !isempty(error_rates)
        CairoMakie.axislegend(axis_a; position = :lb, labelsize = 11)
        CairoMakie.axislegend(axis_b; position = :rb, labelsize = 11)
    end

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
    # MINIMUM, not maximum. This caption exists to warn about thin replication, so
    # taking the best-replicated cell would let one well-replicated cell hide the
    # unreplicated ones it is meant to flag.
    # `init` is the value FOLDED INTO the reduction, not an empty-collection fallback.
    # `minimum(v; init = 0)` is therefore 0 for any non-negative v, which pinned
    # seeds_per_cell to 0 and made this caption print "n = 1 SEED PER CELL" on every
    # sweep including well-replicated ones — inverting the fix it was part of, which
    # moved from maximum to minimum precisely to stop overstating replication.
    # `maximum(...; init = 0)` was safe by luck; the reducer changed and the idiom
    # did not. Handle the empty case explicitly instead.
    counts = collect(values(cell_counts))
    seeds_per_cell = isempty(counts) ? 0 : minimum(counts)
    seeds_max = isempty(counts) ? 0 : maximum(counts)
    replication = isempty(counts) ? "NO CELLS matched — nothing was plotted" :
                  seeds_per_cell > 1 ?
                  (seeds_per_cell == seeds_max ? "$(seeds_per_cell) seeds/cell" :
                   "$(seeds_per_cell)-$(seeds_max) seeds/cell") :
                  "n = 1 SEED PER CELL for at least one cell — no variance estimate there; " *
                  "treat the deltas as provisional"
    deltas = String[]
    for regime in regimes
        pair = get(headline, regime, nothing)
        # Skip a regime missing either arm rather than printing a fabricated 0.
        (pair === nothing || pair[1] === nothing || pair[2] === nothing) && continue
        push!(deltas, "$(regime_label(regime)) $(Int(round(pair[1]))) -> $(Int(round(pair[2])))")
    end
    # Emit the QUAST-shortfall clause ONLY when there is a shortfall, and take the
    # reason from the data rather than asserting one.
    #
    # This was a literal that printed even when quast_ok == nrow(df) — precisely the
    # "caption that always warns" defect this file's tests exist to pin, still live in
    # the delivered figure. Its stated cause was also wrong: the committed sweep
    # records metric_source = "internal:quast-failed" for all 8 unscored rows (QUAST
    # itself failed), not a min-contig shortfall — and 4 of those 8 have
    # largest_contig > 500, so they would have cleared the ceiling #439 settled on.
    # The CSV carries no quast_min_contig column, so the caption cannot name a
    # threshold even in principle; metric_source is the fact it does have.
    sources = hasproperty(df, :metric_source) ?
              String.(collect(skipmissing(df.metric_source))) : String[]
    quast_clause = quast_coverage_clause(quast_ok, DataFrames.nrow(df), sources)
    CairoMakie.Label(figure[0, 1:2],
        "Iterative correction vs naive at the lowest simulated error rate " *
        "(err = $(lowest)), by largest contig — the\nquantity both panels plot: " *
        "$(join(deltas, "; ")).\n$(quast_clause)\n" *
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
    # Load each table ONCE and validate the loaded frame. The first version parsed the
    # primary table twice (once to check, once to use) and never checked the replicate
    # table at all — so a schema problem there still surfaced mid-render in figure 3,
    # which is exactly what the preflight exists to prevent.
    isfile(track_a_path) || error("--track-a file not found: $(track_a_path)")
    track_a = load_track_a(track_a_path)
    # Figure 3 compares on these metric columns, so they are required too.
    cross_host_metrics = [:n_reads, :n_contigs, :NGA50, :genome_fraction,
        :duplication_ratio, :largest_contig, :misassemblies]
    require_columns(track_a, track_a_path,
        vcat(
            [:organism, :technology, :coverage, :seed, :decoder_arm, :genome_fraction,
                :n_contigs],
            cross_host_metrics) |> unique)
    println("Track A cells: $(DataFrames.nrow(track_a)) from $(track_a_path)")
    # The row count above is NOT the denominator any figure uses — ok_cells removes
    # non-measurements first — so print the accounting immediately beside it. Without
    # this line the run reported "Track A cells: 288" and then computed on whatever
    # survived, with the drop count appearing nowhere: not in stdout, not in a warning,
    # not in a caption. Reproduced with QUAST failing on ~1/3 of cells: 105 rows
    # vanished in silence and the headline verdict moved 9/40 -> 5/30.
    exclusions = exclusion_counts(track_a)
    println("  measurements: $(exclusions.n_kept); excluded: " *
            "$(exclusions.n_non_ok) with status != ok, " *
            "$(exclusions.n_unscored) QUAST-unscored " *
            "(n_contigs > 0 and largest_contig == 0)")
    # Warn as well as print, matching write_power_analysis in
    # benchmarking/track_a_baseline_benchmark.jl — a stdout line scrolls past, and the
    # exclusion is one-directional on the verdict via the group-size denominator.
    exclusions.n_non_ok + exclusions.n_unscored > 0 &&
        @warn("figures exclude non-ok and QUAST-unscored cells; group sizes are of "*
            "SURVIVING seeds, so a dropped seed can make a group non-computable",
            n_non_ok=exclusions.n_non_ok, n_unscored=exclusions.n_unscored,
            n_kept=exclusions.n_kept, n_total=exclusions.n_total)

    replicate = nothing
    if replicate_path !== nothing
        isfile(replicate_path) || error("--track-a-replicate not found: $(replicate_path)")
        replicate = load_track_a(replicate_path)
        require_columns(replicate, replicate_path,
            vcat([:organism, :technology, :coverage, :seed, :decoder_arm],
                cross_host_metrics) |> unique)
    end

    rgv = nothing
    if rgv_path !== nothing
        isfile(rgv_path) || error("--rgv file not found: $(rgv_path)")
        rgv = CSV.read(rgv_path, DataFrames.DataFrame)
        require_columns(rgv, rgv_path,
            [:regime, :arm, :error_rate, :largest_contig, :runtime_s, :genome_len])
    end

    cv_table = coefficient_of_variation_table(track_a)
    summary = figure_cv_vs_threshold(cv_table, outdir; exclusions = exclusions)
    # State WHICH kind of undefined, for the same reason the caption does: "(NGA50 = 0)"
    # was an assertion about groups that may simply have too few surviving seeds.
    println("CV verdict: $(summary.exceed)/$(summary.total) computable groups exceed " *
            "$(CV_THRESHOLD); $(summary.undefined) groups undefined " *
            "($(summary.undefined_all_zero) all-zero NGA50, " *
            "$(summary.undefined_too_few_seeds) with < 2 surviving seeds).")

    figure_technology_envelope(track_a, outdir)

    if replicate !== nothing
        result = figure_cross_host(track_a, replicate, outdir)
        println("Cross-host: $(result.shared) shared cells; " *
                "differing by technology = $(result.differing)")
    else
        println("(no --track-a-replicate given; skipping cross-host figure)")
    end

    if rgv !== nothing
        figure_correction_sweep(rgv, outdir)
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
