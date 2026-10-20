# td-2rxh pre-registered paired analysis of `indel_reroute_evidence.jl --cells td2rxh`.
# ===================================================================================
#
# This file implements the decision rule recorded on bead td-2rxh BEFORE any td2rxh
# run: the original pre-registration plus AMENDMENTS 1-3 (also recorded before any
# data). The bead notes are the authoritative text; this is a condensed
# transcription, not a verbatim copy. Summary:
#
#   CELLS: 6000 bp genome x read length {1000, 2000} x 20x x err {0.01, 0.02, 0.05};
#       16 seeds (1..16); arms differ only in assemble_genome sequencing_tech.
#   H1' (decode cost): per cell, median over seeds of the PAIRED incremental cost
#       (wall_nanopore - wall_illumina) / n_reads. A cell is H1'-complete when it
#       has 16 pairs and >= 1 engaged nanopore run, whether or not its coverage
#       deltas are defined. HOLDS only if all 6 cells are H1'-complete with median
#       <= 0.30 s/read; REFUTED if any H1'-complete cell exceeds 0.30; otherwise
#       INCOMPLETE. The difference includes rerouting/window overhead, so it is an
#       upper bound on decode cost alone.
#   H2 (contiguity): per cell, delta = nanopore - illumina on the PRIMARY metric
#       best_contig_reference_coverage (a contiguity measure, not sequence
#       accuracy; see indel_bench_best_reference_alignment), rounded to 1e-6.
#       Exact two-sided sign-flip permutation test on the Wilcoxon signed-rank
#       statistic (average ranks, zeros dropped), exact under ties. Holm across the
#       6 pre-registered cells, fixed family: a cell that is MISSING, INCOMPLETE
#       (< 16 paired seeds), UNDEFINED (no non-zero delta) or TREATMENT-NOT-APPLIED
#       (no nanopore run with indel_engaged >= 1) enters at p = 1 and cannot
#       qualify. H2 holds in a cell iff it is COMPLETE, median delta > 0 and
#       Holm-adjusted p < 0.05.
#   ACCEPTANCE: MET iff, for at least one error rate, H2 holds at BOTH read lengths
#       and the median best_contig_fit_identity delta (over pairs where both arms
#       are finite) is >= 0 in both. NOT MET if none qualifies and every cell is
#       COMPLETE or UNDEFINED (a measured null: more data cannot change it).
#       INCOMPLETE whenever any cell is MISSING, INCOMPLETE or TREATMENT-NOT-APPLIED
#       (conservative: the amendment-2 list, applied whether or not that cell could
#       still change the outcome). Fit identity carries no p-value but gates
#       acceptance; n_contigs and the divergence/anchor counters are written per arm
#       to td2rxh_paired.csv for inspection and are not tested.
#
# Input integrity is enforced, not assumed: malformed or mixed run labels, cells
# outside the pre-registered six, seeds outside 1..16, duplicate (cell, seed, arm)
# rows, and (amendment 3) a pair with no contigs in either arm -- a pipeline
# failure that would otherwise score as a confident null -- are errors, and every
# expected (cell, seed) that lacks a pair is listed in the report with a reason.
#
# Usage:
#   julia --project=. benchmarking/td2rxh_paired_analysis.jl \
#       --csv seed1.csv,seed2.csv,... --output-dir benchmarking/results/td2rxh-analysis

import CairoMakie
import CSV
import DataFrames
import Distributions
import Statistics

# `average_ranks` (tie-averaged ranks). rgv_paired_wilcoxon.jl's own
# `wilcoxon_signed_rank` is deliberately NOT used: it is exact only without ties or
# zeros and otherwise switches to the normal approximation, and at n = 16 under
# Holm x6 that branch switch can decide a cell.
include(joinpath(@__DIR__, "rgv_paired_wilcoxon.jl"))

const TD2RXH_EXPECTED_CELLS = [
    "6000/1000/20x/e0.01", "6000/1000/20x/e0.02", "6000/1000/20x/e0.05",
    "6000/2000/20x/e0.01", "6000/2000/20x/e0.02", "6000/2000/20x/e0.05"]
const TD2RXH_N_SEEDS = 16
const TD2RXH_H1_MAX_INCREMENTAL_SECONDS_PER_READ = 0.30
const TD2RXH_ALPHA = 0.05
const TD2RXH_READ_LENGTHS = (1000, 2000)
const TD2RXH_DELTA_DIGITS = 6
const TD2RXH_MAX_EXACT_N = 20
const TD2RXH_PRIMARY_METRIC = :best_contig_reference_coverage
const TD2RXH_PAIRED_METRICS = (:best_contig_reference_coverage, :best_contig_fit_identity,
    :n_contigs, :largest_contig, :substitution_length_divergences,
    :window_anchor_rejections, :indel_engaged, :wall_seconds, :n_reads)
const TD2RXH_REQUIRED_COLUMNS = (:label, :cell, :seed, :arm, TD2RXH_PAIRED_METRICS...)
# The label run_td2rxh_evidence_lrc.sbatch writes: HEAD and Manifest.toml digests.
const TD2RXH_LABEL_PATTERN = r"^td2rxh@[0-9a-f]{12}\+manifest[0-9a-f]{12}$"

"""
    signflip_signed_rank_test(diffs) -> NamedTuple

Exact two-sided sign-flip permutation test on the Wilcoxon signed-rank statistic.
Differences are rounded to `TD2RXH_DELTA_DIGITS` decimals so floating-point noise
cannot create or break ties; zeros are dropped and counted; the statistic is
`T = sum(sign(d) * rank(|d|))` with tie-averaged ranks; the p-value is the exact
fraction of all `2^n` sign assignments with `|T| >= |T_obs|`. Exact under ties,
unlike the textbook exact-null DP, which assumes untied ranks. With no non-zero
difference the test is undefined and `pvalue` is `NaN`.

Throws `ArgumentError` on a non-finite difference (coverage is always finite --
an empty assembly scores 0.0 -- so this signals corrupt input, and the whole
analysis stops rather than scoring it) and when `n > TD2RXH_MAX_EXACT_N`.
"""
function signflip_signed_rank_test(diffs::AbstractVector{<:Real})::NamedTuple
    rounded = round.(Float64.(diffs); digits = TD2RXH_DELTA_DIGITS)
    any(!isfinite, rounded) &&
        throw(ArgumentError("non-finite paired difference; drop or fix the pair first"))
    nonzero = filter(!iszero, rounded)
    n = length(nonzero)
    n_zero = length(rounded) - n
    n == 0 && return (n = 0, n_zero_dropped = n_zero, statistic = NaN, pvalue = NaN,
        method = "undefined (no non-zero differences)")
    n > TD2RXH_MAX_EXACT_N &&
        throw(ArgumentError("n = $(n) exceeds the exact enumeration limit $(TD2RXH_MAX_EXACT_N)"))
    ranks = average_ranks(abs.(nonzero))
    observed = abs(sum(sign.(nonzero) .* ranks))
    # Ranks are multiples of 1/2, so a tolerance far below 1/2 compares exactly.
    tolerance = 1e-9
    extreme = 0
    for mask in 0:((1 << n) - 1)
        t = 0.0
        for i in 1:n
            t += ((mask >> (i - 1)) & 1) == 1 ? ranks[i] : -ranks[i]
        end
        abs(t) >= observed - tolerance && (extreme += 1)
    end
    return (n = n, n_zero_dropped = n_zero, statistic = observed,
        pvalue = extreme / (1 << n), method = "exact sign-flip ($(1 << n) assignments)")
end

"""
    holm_adjust(pvalues) -> Vector{Float64}

Holm step-down adjusted p-values over a FIXED family: every entry counts toward
the family size, and a `NaN` entry is treated as `p = 1`. Pass all six
pre-registered cells, never only the ones that produced a p-value: excluding a
cell shrinks the multiplier and makes every other cell easier to pass.
"""
function holm_adjust(pvalues::AbstractVector{<:Real})::Vector{Float64}
    p = [isnan(x) ? 1.0 : Float64(x) for x in pvalues]
    m = length(p)
    adjusted = similar(p)
    running = 0.0
    for (rank, i) in enumerate(sortperm(p))
        running = max(running, min(1.0, (m - rank + 1) * p[i]))
        adjusted[i] = running
    end
    return adjusted
end

"""
    validate_shards(df) -> Nothing

Refuse inputs the pre-registered rule cannot be applied to: missing columns, a
run label that is absent or not in the sbatch's `td2rxh@<head>+manifest<hash>`
form, more than one label (shards from different commits or dependency sets),
cells outside the six pre-registered ids, seeds outside `1:TD2RXH_N_SEEDS`, or
duplicate (cell, seed, arm) rows. Each is an error rather than a warning because
each would otherwise be silent: pooled labels measure different code, an extra
cell or seed would simply be ignored by the fixed loops in `pair_arms`, and a
duplicate would have to be resolved by picking one row.
"""
function validate_shards(df::DataFrames.DataFrame)::Nothing
    missing_cols = [c for c in TD2RXH_REQUIRED_COLUMNS if !(string(c) in names(df))]
    isempty(missing_cols) ||
        throw(ArgumentError("input lacks required columns: $(join(missing_cols, ", "))"))
    any(ismissing, df.label) && throw(ArgumentError("rows with an empty run label"))
    labels = unique(string.(df.label))
    length(labels) == 1 ||
        throw(ArgumentError("shards carry $(length(labels)) run labels ($(join(labels, ", "))); " *
                            "the analysis requires one commit"))
    occursin(TD2RXH_LABEL_PATTERN, only(labels)) ||
        throw(ArgumentError("run label \"$(only(labels))\" is not td2rxh@<12 hex>+manifest<12 hex>"))
    unknown = setdiff(unique(string.(df.cell)), TD2RXH_EXPECTED_CELLS)
    isempty(unknown) ||
        throw(ArgumentError("cells outside the pre-registered six: $(join(unknown, ", "))"))
    stray = setdiff(unique(Int.(df.seed)), 1:TD2RXH_N_SEEDS)
    isempty(stray) ||
        throw(ArgumentError("seeds outside 1:$(TD2RXH_N_SEEDS): $(join(sort(stray), ", "))"))
    keys = [(string(r.cell), Int(r.seed), string(r.arm)) for r in eachrow(df)]
    duplicated = unique(k for k in keys if count(==(k), keys) > 1)
    isempty(duplicated) ||
        throw(ArgumentError("duplicate (cell, seed, arm) rows: $(join(string.(duplicated), "; "))"))
    return nothing
end

"""
    pair_arms(df) -> (paired::DataFrame, unpaired::DataFrame)

One row per (cell, seed) with both arms side by side (`<col>_illumina`,
`<col>_nanopore`). Every EXPECTED (cell, seed) without both arms -- including ones
with no rows at all -- is listed in `unpaired` with a reason, so a timed-out task
or an arm that always throws is visible rather than silently absent.
"""
function pair_arms(df::DataFrames.DataFrame;
        expected_cells::Vector{String} = TD2RXH_EXPECTED_CELLS,
        n_seeds::Int = TD2RXH_N_SEEDS)::Tuple{DataFrames.DataFrame, DataFrames.DataFrame}
    paired_rows = NamedTuple[]
    unpaired_rows = NamedTuple[]
    for cell in expected_cells, seed in 1:n_seeds

        rows = df[(string.(df.cell) .== cell) .& (Int.(df.seed) .== seed), :]
        arms = string.(rows.arm)
        has_ill, has_ont = "illumina" in arms, "nanopore" in arms
        if !(has_ill && has_ont)
            reason = !has_ill && !has_ont ? "no rows" :
                     has_ill ? "illumina only" : "nanopore only"
            push!(unpaired_rows, (cell = cell, seed = seed, reason = reason))
            continue
        end
        ill = rows[findfirst(==("illumina"), arms), :]
        ont = rows[findfirst(==("nanopore"), arms), :]
        values = Pair{Symbol, Any}[:cell => cell, :seed => seed]
        for metric in TD2RXH_PAIRED_METRICS
            push!(values, Symbol(metric, "_illumina") => ill[metric])
            push!(values, Symbol(metric, "_nanopore") => ont[metric])
        end
        push!(paired_rows, (; values...))
    end
    template = DataFrames.DataFrame(
        [Symbol(m, suffix) => Any[] for m in TD2RXH_PAIRED_METRICS
         for suffix in ("_illumina", "_nanopore")])
    paired = isempty(paired_rows) ?
             DataFrames.hcat(DataFrames.DataFrame(cell = String[], seed = Int[]), template) :
             DataFrames.DataFrame(paired_rows)
    unpaired = isempty(unpaired_rows) ?
               DataFrames.DataFrame(cell = String[], seed = Int[], reason = String[]) :
               DataFrames.DataFrame(unpaired_rows)
    return paired, unpaired
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

"""An `indel_engaged` value counts as engaged only when it parses to an integer >= 1."""
function _engaged(value)::Bool
    value isa Integer && return value >= 1
    value isa Real && return isfinite(value) && value >= 1
    parsed = tryparse(Int, string(value))
    return parsed !== nothing && parsed >= 1
end

function _finite_median(x)
    (f = filter(isfinite, Float64.(x)); isempty(f) ? NaN : Statistics.median(f))
end

"""
    cell_verdicts(paired) -> DataFrame

One row per PRE-REGISTERED cell (always six, in pre-registered order), with its
completeness status, H1' and H2 statistics, and the Holm adjustment over the fixed
six-cell family.
"""
function cell_verdicts(paired::DataFrames.DataFrame;
        expected_cells::Vector{String} = TD2RXH_EXPECTED_CELLS,
        n_seeds::Int = TD2RXH_N_SEEDS)::DataFrames.DataFrame
    rows = NamedTuple[]
    for cell in expected_cells
        parsed = parse_cell_id(cell)
        group = paired[string.(paired.cell) .== cell, :]
        n_pairs = DataFrames.nrow(group)
        primary_ill = Float64.(group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_illumina")])
        primary_ont = Float64.(group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_nanopore")])
        delta = primary_ont .- primary_ill
        incremental = (Float64.(group.wall_seconds_nanopore) .-
                       Float64.(group.wall_seconds_illumina)) ./
                      Float64.(group.n_reads_illumina)
        fit = Float64.(group.best_contig_fit_identity_nanopore) .-
              Float64.(group.best_contig_fit_identity_illumina)
        n_engaged = count(_engaged, group.indel_engaged_nanopore)
        empty_ill = Int.(group.n_contigs_illumina) .== 0
        empty_ont = Int.(group.n_contigs_nanopore) .== 0
        # Amendment 3: a pair with no contigs in EITHER arm is a pipeline failure,
        # not a measured null; scored as 0.0 vs 0.0 it would make the cell
        # UNDEFINED and the verdict a confident NOT MET.
        any(empty_ill .& empty_ont) &&
            throw(ArgumentError("$(cell): seeds $(group.seed[empty_ill .& empty_ont]) " *
                                "produced no contigs in either arm"))
        no_reads = (Int.(group.n_reads_illumina) .<= 0) .|
                   (Int.(group.n_reads_nanopore) .<= 0)
        any(no_reads) &&
            throw(ArgumentError("$(cell): seeds $(group.seed[no_reads]) report n_reads == 0"))
        test = n_pairs == 0 ?
               (
            n = 0, n_zero_dropped = 0, statistic = NaN, pvalue = NaN, method = "no pairs") :
               signflip_signed_rank_test(delta)
        status = n_pairs == 0 ? "MISSING" :
                 n_pairs < n_seeds ? "INCOMPLETE" :
                 n_engaged == 0 ? "TREATMENT-NOT-APPLIED" :
                 test.n == 0 ? "UNDEFINED" : "COMPLETE"
        # Amendment 2: H1' completeness is about timing data, not about whether the
        # coverage deltas happen to be defined.
        h1_complete = n_pairs == n_seeds && n_engaged > 0
        median_incremental = n_pairs == 0 ? NaN : Statistics.median(incremental)
        if h1_complete
            Int.(group.n_reads_illumina) == Int.(group.n_reads_nanopore) ||
                throw(ArgumentError("$(cell): arms report different read counts"))
            isfinite(median_incremental) ||
                throw(ArgumentError("$(cell): non-finite incremental cost $(median_incremental)"))
        end
        push!(rows,
            (cell = cell, read_length = parsed.read_length,
                error_rate = parsed.error_rate, status = status, n_pairs = n_pairs,
                n_nanopore_engaged = n_engaged,
                n_empty_illumina = count(empty_ill), n_empty_nanopore = count(empty_ont),
                median_incremental_s_per_read = median_incremental,
                h1_cell = h1_complete ?
                          (median_incremental <=
                           TD2RXH_H1_MAX_INCREMENTAL_SECONDS_PER_READ ?
                           "holds" : "exceeds") : "n/a",
                primary_illumina_median = n_pairs == 0 ? NaN :
                                          Statistics.median(primary_ill),
                primary_nanopore_median = n_pairs == 0 ? NaN :
                                          Statistics.median(primary_ont),
                median_primary_delta = n_pairs == 0 ? NaN : Statistics.median(delta),
                test_n_nonzero = test.n, test_zero_dropped = test.n_zero_dropped,
                test_method = test.method, p_raw = test.pvalue,
                p_family = status == "COMPLETE" ? test.pvalue : 1.0,
                median_fit_identity_delta = _finite_median(fit),
                n_fit_pairs = count(isfinite, fit)))
    end
    verdicts = DataFrames.DataFrame(rows)
    verdicts.p_holm = holm_adjust(verdicts.p_family)
    verdicts.h2_holds = [(s == "COMPLETE") && (d > 0) && (p < TD2RXH_ALPHA)
                         for (s, d, p) in
                             zip(verdicts.status, verdicts.median_primary_delta,
        verdicts.p_holm)]
    return verdicts
end

"""
    h1_verdict(verdicts) -> String

"REFUTED" if any H1'-complete cell's median exceeds the threshold, "HOLDS" only
when all six cells are H1'-complete and within it, otherwise "INCOMPLETE". A
missing cell can therefore never be read as H1' holding, and an UNDEFINED cell's
timing still counts (amendment 2).
"""
function h1_verdict(verdicts::DataFrames.DataFrame)::String
    any(verdicts.h1_cell .== "exceeds") && return "REFUTED"
    all(verdicts.h1_cell .== "holds") && return "HOLDS"
    return "INCOMPLETE"
end

"""Cell states after which no further data can change the cell's contribution."""
const TD2RXH_FINAL_STATUSES = ("COMPLETE", "UNDEFINED")

"""
    acceptance(verdicts) -> NamedTuple

The pre-registered td-2rxh acceptance rule (with amendments 1-3). Returns the
qualifying error rates and a three-state verdict: "MET"; "NOT MET" (none
qualifies and every cell is COMPLETE or UNDEFINED); or "INCOMPLETE" (any cell
MISSING, INCOMPLETE or TREATMENT-NOT-APPLIED -- conservative, whether or not that
cell could still change the outcome).
"""
function acceptance(verdicts::DataFrames.DataFrame)::NamedTuple
    qualifying = Float64[]
    for error_rate in sort(unique(verdicts.error_rate))
        at_rate = verdicts[verdicts.error_rate .== error_rate, :]
        rows = [findfirst(==(rl), at_rate.read_length) for rl in TD2RXH_READ_LENGTHS]
        any(isnothing, rows) && continue
        cells = at_rate[collect(rows), :]
        if all(cells.h2_holds) && all(isfinite, cells.median_fit_identity_delta) &&
           all(cells.median_fit_identity_delta .>= 0)
            push!(qualifying, error_rate)
        end
    end
    verdict = !isempty(qualifying) ? "MET" :
              all(in(TD2RXH_FINAL_STATUSES), verdicts.status) ? "NOT MET" : "INCOMPLETE"
    return (qualifying_error_rates = qualifying, verdict = verdict)
end

function write_figure(paired::DataFrames.DataFrame, verdicts::DataFrames.DataFrame,
        output_dir::AbstractString)::Nothing
    cells = verdicts.cell
    ticks = ["$(c)\n$(s), n=$(n)"
             for (c, s, n) in zip(cells, verdicts.status, verdicts.n_pairs)]
    fig = CairoMakie.Figure(size = (1200, 520))
    ax1 = CairoMakie.Axis(fig[1, 1]; xticks = (1:length(cells), ticks),
        xticklabelrotation = pi / 5,
        ylabel = "nanopore − illumina\nbest-contig reference coverage",
        title = "H2: paired contiguity delta per seed")
    ax2 = CairoMakie.Axis(fig[1, 2]; xticks = (1:length(cells), ticks),
        xticklabelrotation = pi / 5,
        ylabel = "(wall_nanopore − wall_illumina) / read  [s]",
        title = "H1': incremental indel-profile cost per read")
    CairoMakie.hlines!(ax1, [0.0]; color = :gray40, linestyle = :dash)
    CairoMakie.hlines!(ax2, [TD2RXH_H1_MAX_INCREMENTAL_SECONDS_PER_READ];
        color = :firebrick, linestyle = :dash, label = "0.30 s/read threshold")
    for (x, cell) in enumerate(cells)
        group = paired[string.(paired.cell) .== cell, :]
        DataFrames.nrow(group) == 0 && continue
        delta = Float64.(group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_nanopore")]) .-
                Float64.(group[!, Symbol(TD2RXH_PRIMARY_METRIC, "_illumina")])
        incremental = (Float64.(group.wall_seconds_nanopore) .-
                       Float64.(group.wall_seconds_illumina)) ./
                      Float64.(group.n_reads_illumina)
        jitter = collect(range(-0.2, 0.2; length = max(length(delta), 2)))[1:length(delta)]
        CairoMakie.scatter!(ax1, x .+ jitter, delta; color = :steelblue)
        CairoMakie.scatter!(ax2, x .+ jitter, incremental; color = :darkorange)
        CairoMakie.text!(ax1, x, maximum(delta; init = 0.0);
            text = "p_Holm=$(round(verdicts.p_holm[x]; sigdigits = 2))",
            align = (:center, :bottom), fontsize = 10)
    end
    CairoMakie.axislegend(ax2; position = :rt)
    CairoMakie.Label(fig[2, 1:2],
        "6 kb synthetic genome, 20x, Mycelia.observe tech=:nanopore reads; both arms see " *
        "identical reads and differ only in assemble_genome sequencing_tech. Exact " *
        "sign-flip signed-rank test, Holm over the 6 pre-registered cells. PRELIMINARY tier.";
        fontsize = 10, tellwidth = false)
    for ext in ("png", "svg")
        CairoMakie.save(joinpath(output_dir, "td2rxh_paired_deltas.$(ext)"), fig)
    end
    return nothing
end

"""
Parse `--key value` pairs, rejecting unknown keys, a dangling flag, and a value
that looks like a flag, so a typo cannot silently fall back to a default.
"""
function parse_analysis_args(args::Vector{String})::Dict{String, String}
    allowed = ("csv", "output-dir")
    isodd(length(args)) && throw(ArgumentError("expected --key value pairs, got $(args)"))
    parsed = Dict{String, String}()
    for i in 1:2:length(args)
        startswith(args[i], "--") ||
            throw(ArgumentError("expected a --flag, got $(args[i])"))
        key = args[i][3:end]
        key in allowed || throw(ArgumentError("unknown flag --$(key)"))
        startswith(args[i + 1], "--") &&
            throw(ArgumentError("--$(key) is missing its value (got $(args[i + 1]))"))
        parsed[key] = args[i + 1]
    end
    haskey(parsed, "csv") || throw(ArgumentError("--csv is required"))
    # Required, not defaulted: the output directory's previous report is removed
    # before reading, so a shared default would let one job's mistaken run delete
    # another job's report.
    haskey(parsed, "output-dir") || throw(ArgumentError("--output-dir is required"))
    return parsed
end

_fmt(x::Real; digits = 4) = isnan(x) ? "NaN" : string(round(x; sigdigits = digits))
_fmt(x) = string(x)

function write_report(io::IO, csv_paths, label, paired, unpaired, verdicts, h1, verdict)
    println(io, "# td-2rxh pre-registered paired analysis\n")
    println(io, "Run label: `$(label)`; shards: $(length(csv_paths))\n")
    println(io, "Expected pairs: $(length(TD2RXH_EXPECTED_CELLS) * TD2RXH_N_SEEDS); ",
        "paired: $(DataFrames.nrow(paired)); unpaired: $(DataFrames.nrow(unpaired))\n")
    println(
        io, "| cell | status | pairs | nanopore engaged | empty asm ill/ont | ",
        "median incr s/read | H1' cell | ",
        "cov illumina | cov nanopore | median Δcov | nonzero | zeros | p raw | p Holm | H2 | ",
        "median Δfit-id (pairs) |")
    println(io, "|", repeat("---|", 16))
    for r in eachrow(verdicts)
        println(
            io, "| $(r.cell) | $(r.status) | $(r.n_pairs) | $(r.n_nanopore_engaged) | ",
            "$(r.n_empty_illumina)/$(r.n_empty_nanopore) | ",
            "$(_fmt(r.median_incremental_s_per_read)) | $(r.h1_cell) | ",
            "$(_fmt(r.primary_illumina_median)) | $(_fmt(r.primary_nanopore_median)) | ",
            "$(_fmt(r.median_primary_delta)) | $(r.test_n_nonzero) | $(r.test_zero_dropped) | ",
            "$(_fmt(r.p_raw; digits = 3)) | $(_fmt(r.p_holm; digits = 3)) | $(r.h2_holds) | ",
            "$(_fmt(r.median_fit_identity_delta)) ($(r.n_fit_pairs)) |")
    end
    if DataFrames.nrow(unpaired) > 0
        println(io, "\n## Unpaired expected (cell, seed)\n")
        for r in eachrow(unpaired)
            println(io, "- $(r.cell) seed $(r.seed): $(r.reason)")
        end
    end
    println(io, "\nH1': $(h1)")
    println(io, "Qualifying error rates: $(verdict.qualifying_error_rates)")
    println(io, "td-2rxh ACCEPTANCE: $(verdict.verdict)")
end

function td2rxh_analysis_main(args::Vector{String} = ARGS)::Int
    parsed = parse_analysis_args(args)
    output_dir = parsed["output-dir"]
    mkpath(output_dir)
    # Remove the previous generation BEFORE reading anything, so any failure --
    # validation, a statistic, or part-way through writing -- leaves no stale
    # report that could be read as the current result.
    for name in ("report.md", "td2rxh_cell_verdicts.csv", "td2rxh_paired.csv",
        "td2rxh_unpaired.csv", "td2rxh_paired_deltas.png", "td2rxh_paired_deltas.svg")
        rm(joinpath(output_dir, name); force = true)
    end
    csv_paths = String.(split(parsed["csv"], ","))
    df = reduce(vcat, [CSV.read(path, DataFrames.DataFrame) for path in csv_paths])
    validate_shards(df)
    paired, unpaired = pair_arms(df)
    verdicts = cell_verdicts(paired)
    h1 = h1_verdict(verdicts)
    verdict = acceptance(verdicts)
    CSV.write(joinpath(output_dir, "td2rxh_cell_verdicts.csv"), verdicts)
    CSV.write(joinpath(output_dir, "td2rxh_paired.csv"), paired)
    CSV.write(joinpath(output_dir, "td2rxh_unpaired.csv"), unpaired)
    write_figure(paired, verdicts, output_dir)
    open(joinpath(output_dir, "report.md"), "w") do io
        write_report(io, csv_paths, first(string.(df.label)), paired, unpaired,
            verdicts, h1, verdict)
    end
    println(read(joinpath(output_dir, "report.md"), String))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(td2rxh_analysis_main())
end
