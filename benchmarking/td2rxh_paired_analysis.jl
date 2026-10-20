# td-2rxh pre-registered paired analysis of `indel_reroute_evidence.jl --cells td2rxh`.
# ===================================================================================
#
# The decision rule below was recorded on bead td-2rxh BEFORE any td2rxh run. It is
# reproduced here verbatim in code so the analysis cannot drift from it:
#
#   H1' (decode cost no longer blocking): per cell, the median over seeds of the
#       PAIRED incremental cost (wall_nanopore - wall_illumina) / n_reads is
#       <= 0.30 s/read. Refuted if ANY cell's median exceeds 0.30.
#   H2  (accuracy): per cell, delta = nanopore - illumina on the PRIMARY metric
#       best_contig_reference_coverage, paired over seeds. Two-sided Wilcoxon
#       signed-rank, Holm across the cells, alpha 0.05. "nanopore beats illumina"
#       in a cell = median delta > 0 AND Holm-adjusted p < 0.05.
#   ACCEPTANCE: met iff H2 holds at BOTH read lengths 1000 and 2000 for at least
#       one error rate AND the median best_contig_fit_identity delta is >= 0 in
#       those cells.
#
# Unpaired rows (an arm that threw writes no row) are DROPPED AND REPORTED, never
# imputed. The fit-identity, n_contigs and counter deltas are descriptive only.
#
# Usage:
#   julia --project=. benchmarking/td2rxh_paired_analysis.jl \
#       --csv seed1.csv,seed2.csv,... --output-dir benchmarking/results/td2rxh-analysis

import CairoMakie
import CSV
import DataFrames
import Distributions
import Statistics

# `wilcoxon_signed_rank` (exact when valid, tie-corrected normal otherwise),
# cross-validated against scipy in test/4_assembly/rgv_paired_wilcoxon_test.jl.
include(joinpath(@__DIR__, "rgv_paired_wilcoxon.jl"))

const TD2RXH_H1_MAX_INCREMENTAL_SECONDS_PER_READ = 0.30
const TD2RXH_ALPHA = 0.05
const TD2RXH_READ_LENGTHS = (1000, 2000)
const TD2RXH_PRIMARY_METRIC = :best_contig_reference_coverage
const TD2RXH_SECONDARY_METRICS = (:best_contig_fit_identity, :n_contigs,
    :largest_contig, :substitution_length_divergences, :window_anchor_rejections,
    :indel_engaged)

"""
    holm_adjust(pvalues) -> Vector{Float64}

Holm step-down adjusted p-values, monotone and capped at 1. `NaN` inputs (an
undefined test) stay `NaN` and are excluded from the family size, so an
undefined cell can never make another cell look significant.
"""
function holm_adjust(pvalues::AbstractVector{<:Real})::Vector{Float64}
    adjusted = fill(NaN, length(pvalues))
    defined = [i for i in eachindex(pvalues) if !isnan(pvalues[i])]
    m = length(defined)
    order = sort(defined; by = i -> pvalues[i])
    running = 0.0
    for (rank, i) in enumerate(order)
        running = max(running, min(1.0, (m - rank + 1) * pvalues[i]))
        adjusted[i] = running
    end
    return adjusted
end

"""
    pair_arms(df) -> (paired::DataFrame, dropped::DataFrame)

One row per (cell, seed) holding both arms side by side (`<col>_illumina`,
`<col>_nanopore`). A (cell, seed) missing an arm, or with duplicate rows for an
arm, goes to `dropped` with a `reason`.
"""
function pair_arms(df::DataFrames.DataFrame)::Tuple{
        DataFrames.DataFrame, DataFrames.DataFrame}
    metrics = vcat([TD2RXH_PRIMARY_METRIC], collect(TD2RXH_SECONDARY_METRICS),
        [:wall_seconds, :n_reads])
    paired_rows = NamedTuple[]
    dropped_rows = NamedTuple[]
    for group in DataFrames.groupby(df, [:cell, :seed])
        cell = string(first(group.cell))
        seed = Int(first(group.seed))
        n_ill = count(==("illumina"), string.(group.arm))
        n_ont = count(==("nanopore"), string.(group.arm))
        if n_ill != 1 || n_ont != 1
            push!(dropped_rows,
                (cell = cell, seed = seed,
                    reason = "arms illumina=$(n_ill) nanopore=$(n_ont)"))
            continue
        end
        ill = group[findfirst(==("illumina"), string.(group.arm)), :]
        ont = group[findfirst(==("nanopore"), string.(group.arm)), :]
        values = Pair{Symbol, Any}[:cell => cell, :seed => seed]
        for metric in metrics
            push!(values, Symbol(metric, "_illumina") => ill[metric])
            push!(values, Symbol(metric, "_nanopore") => ont[metric])
        end
        push!(paired_rows, (; values...))
    end
    paired = DataFrames.DataFrame(paired_rows)
    dropped = isempty(dropped_rows) ?
              DataFrames.DataFrame(cell = String[], seed = Int[], reason = String[]) :
              DataFrames.DataFrame(dropped_rows)
    return paired, dropped
end

"""
Parse the read length and error rate out of a `cell_id` string
(`"<genome>/<readlen>/<cov>x/e<err>"`).
"""
function parse_cell_id(cell::AbstractString)::NamedTuple
    parts = split(cell, "/")
    length(parts) == 4 && startswith(parts[4], "e") ||
        throw(ArgumentError("unrecognised cell id \"$(cell)\""))
    return (genome_length = parse(Int, parts[1]), read_length = parse(Int, parts[2]),
        error_rate = parse(Float64, parts[4][2:end]))
end

"""
    cell_verdicts(paired) -> DataFrame

Per-cell H1' and H2 statistics and verdicts, with Holm applied across all cells.
"""
function cell_verdicts(paired::DataFrames.DataFrame)::DataFrames.DataFrame
    rows = NamedTuple[]
    for group in DataFrames.groupby(paired, :cell)
        cell = string(first(group.cell))
        parsed = parse_cell_id(cell)
        incremental = (group.wall_seconds_nanopore .- group.wall_seconds_illumina) ./
                      group.n_reads_illumina
        primary = Float64.(group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_nanopore")] .-
                           group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_illumina")])
        fit = Float64.(group.best_contig_fit_identity_nanopore .-
                       group.best_contig_fit_identity_illumina)
        test = wilcoxon_signed_rank(primary)
        push!(rows,
            (cell = cell, read_length = parsed.read_length,
                error_rate = parsed.error_rate, n_pairs = DataFrames.nrow(group),
                median_incremental_s_per_read = Statistics.median(incremental),
                h1_holds = Statistics.median(incremental) <=
                           TD2RXH_H1_MAX_INCREMENTAL_SECONDS_PER_READ,
                median_primary_delta = Statistics.median(primary),
                primary_illumina_median = Statistics.median(
                    group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_illumina")]),
                primary_nanopore_median = Statistics.median(
                    group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_nanopore")]),
                wilcoxon_n = test.n, wilcoxon_zero_dropped = test.n_zero_dropped,
                wilcoxon_method = test.method, p_raw = test.pvalue,
                median_fit_identity_delta = Statistics.median(fit)))
    end
    verdicts = DataFrames.sort!(DataFrames.DataFrame(rows), [:read_length, :error_rate])
    verdicts.p_holm = holm_adjust(verdicts.p_raw)
    verdicts.h2_holds = [(d > 0) && !isnan(p) && (p < TD2RXH_ALPHA)
                         for (d, p) in zip(verdicts.median_primary_delta, verdicts.p_holm)]
    return verdicts
end

"""
    acceptance(verdicts) -> NamedTuple

The pre-registered td-2rxh acceptance rule, evaluated per error rate. Requires
every pre-registered read length to be present for an error rate to qualify; a
missing cell can therefore only make acceptance harder, never easier.
"""
function acceptance(verdicts::DataFrames.DataFrame)::NamedTuple
    h1_all = !isempty(verdicts.h1_holds) && all(verdicts.h1_holds)
    qualifying = Float64[]
    for error_rate in sort(unique(verdicts.error_rate))
        at_rate = verdicts[verdicts.error_rate .== error_rate, :]
        lengths_present = all(rl -> rl in at_rate.read_length, TD2RXH_READ_LENGTHS)
        lengths_present || continue
        rows = at_rate[in.(at_rate.read_length, Ref(TD2RXH_READ_LENGTHS)), :]
        if all(rows.h2_holds) && all(rows.median_fit_identity_delta .>= 0)
            push!(qualifying, error_rate)
        end
    end
    return (h1_holds_all_cells = h1_all, qualifying_error_rates = qualifying,
        accepted = !isempty(qualifying))
end

function write_figure(paired::DataFrames.DataFrame, verdicts::DataFrames.DataFrame,
        output_dir::AbstractString)::Nothing
    cells = verdicts.cell
    fig = CairoMakie.Figure(size = (1100, 480))
    ax1 = CairoMakie.Axis(fig[1, 1]; xticks = (1:length(cells), cells),
        xticklabelrotation = pi / 5,
        ylabel = "nanopore − illumina\nbest-contig reference coverage",
        title = "H2: paired accuracy delta per seed (n = $(maximum(verdicts.n_pairs)))")
    ax2 = CairoMakie.Axis(fig[1, 2]; xticks = (1:length(cells), cells),
        xticklabelrotation = pi / 5,
        ylabel = "(wall_nanopore − wall_illumina) / read  [s]",
        title = "H1': incremental indel-profile cost per read")
    CairoMakie.hlines!(ax1, [0.0]; color = :gray40, linestyle = :dash)
    CairoMakie.hlines!(ax2, [TD2RXH_H1_MAX_INCREMENTAL_SECONDS_PER_READ];
        color = :firebrick, linestyle = :dash, label = "0.30 s/read threshold")
    for (x, cell) in enumerate(cells)
        group = paired[paired.cell .== cell, :]
        delta = group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_nanopore")] .-
                group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_illumina")]
        incremental = (group.wall_seconds_nanopore .- group.wall_seconds_illumina) ./
                      group.n_reads_illumina
        jitter = range(-0.15, 0.15; length = max(length(delta), 2))[1:length(delta)]
        CairoMakie.scatter!(ax1, x .+ jitter, Float64.(delta); color = :steelblue)
        CairoMakie.scatter!(ax2, x .+ jitter, Float64.(incremental); color = :darkorange)
        CairoMakie.text!(ax1, x, Float64(maximum(delta; init = 0.0));
            text = "p_Holm=$(round(verdicts.p_holm[x]; sigdigits = 2))",
            align = (:center, :bottom), fontsize = 10)
    end
    CairoMakie.axislegend(ax2; position = :rt)
    CairoMakie.Label(fig[2, 1:2],
        "6 kb synthetic genome, 20x, Mycelia.observe tech=:nanopore reads; both arms see " *
        "identical reads and differ only in assemble_genome sequencing_tech. " *
        "Two-sided Wilcoxon signed-rank, Holm across cells. PRELIMINARY tier.";
        fontsize = 10, tellwidth = false)
    for ext in ("png", "svg")
        CairoMakie.save(joinpath(output_dir, "td2rxh_paired_deltas.$(ext)"), fig)
    end
    return nothing
end

function main(args::Vector{String} = ARGS)::Int
    parsed = Dict{String, String}()
    index = 1
    while index < length(args)
        parsed[lstrip(args[index], '-')] = args[index + 1]
        index += 2
    end
    haskey(parsed, "csv") ||
        (println(stderr, "usage: --csv PATH --output-dir DIR"); return 2)
    output_dir = get(parsed, "output-dir", "td2rxh-analysis")
    mkpath(output_dir)
    # Comma-separated shards (one per array task / seed). A seed present in two
    # shards produces duplicate arm rows, which pair_arms drops and reports.
    csv_paths = String.(split(parsed["csv"], ","))
    df = reduce(vcat, [CSV.read(path, DataFrames.DataFrame) for path in csv_paths])
    paired, dropped = pair_arms(df)
    verdicts = cell_verdicts(paired)
    verdict = acceptance(verdicts)
    CSV.write(joinpath(output_dir, "td2rxh_cell_verdicts.csv"), verdicts)
    CSV.write(joinpath(output_dir, "td2rxh_paired.csv"), paired)
    CSV.write(joinpath(output_dir, "td2rxh_dropped.csv"), dropped)
    write_figure(paired, verdicts, output_dir)
    open(joinpath(output_dir, "report.md"), "w") do io
        println(io, "# td-2rxh pre-registered paired analysis\n")
        println(io, "Source CSV: `$(parsed["csv"])`\n")
        println(io,
            "Paired (cell, seed): $(DataFrames.nrow(paired)); dropped: $(DataFrames.nrow(dropped))\n")
        println(io,
            "| cell | n | median incr s/read | H1' | cov illumina | cov nanopore | median Δcov | p raw | p Holm | H2 | median Δfit-id |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
        for r in eachrow(verdicts)
            println(io,
                "| $(r.cell) | $(r.n_pairs) | $(round(r.median_incremental_s_per_read; digits = 4)) | ",
                "$(r.h1_holds) | $(round(r.primary_illumina_median; digits = 4)) | ",
                "$(round(r.primary_nanopore_median; digits = 4)) | $(round(r.median_primary_delta; digits = 4)) | ",
                "$(round(r.p_raw; sigdigits = 3)) | $(round(r.p_holm; sigdigits = 3)) | $(r.h2_holds) | ",
                "$(round(r.median_fit_identity_delta; digits = 5)) |")
        end
        println(io, "\nH1' holds in every cell: $(verdict.h1_holds_all_cells)")
        println(io, "Qualifying error rates for acceptance: $(verdict.qualifying_error_rates)")
        println(io, "td-2rxh ACCEPTANCE: $(verdict.accepted ? "MET" : "NOT MET")")
    end
    println(read(joinpath(output_dir, "report.md"), String))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
