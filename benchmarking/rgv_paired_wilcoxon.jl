# Pre-registered paired-Wilcoxon analysis of the RGV correction-validation sweep.
# ==============================================================================
#
# WHAT THIS IS, AND WHAT IT IS NOT
# --------------------------------
# It applies the pre-registration's STATISTICAL RULE — paired Wilcoxon signed-rank
# over replicate seeds, Benjamini-Hochberg across metrics, and the
# supported/partial/falsified/null decision thresholds
# (`rhizomorph-paper/planning/PLAN-2026-03-28-rhizomorph-preregistration.md`) — to
# the RGV correction-validation sweep.
#
# It is NOT a pre-registered hypothesis test, and must not be reported as one.
# The pre-registration's H1 is "Full Viterbi DP improves assembly quality over
# greedy Viterbi" (:344), and its procedure says "comparing Viterbi DP vs. greedy"
# (:52). This sweep compares `corrector=:none` against `corrector=:iterative`.
# Grepping the entire pre-registration for `iterativ|corrector|naive` returns ZERO
# hits: none of H1-H7 covers error correction. An earlier version of this header
# quoted H1's procedure with the arms elided to "[treatment] vs [control]", which
# read as though this comparison were pre-registered. It is not.
#
# Two further mismatches that must travel with any number this produces:
#   * The pre-registration's design is 4 coverage levels x 3 replicates per
#     organism/technology (~360 paired comparisons). The RGV sweep runs ONE
#     coverage and varies `error_rate`, an axis the pre-registration never names.
#   * `k` is pooled here; the pre-registration designates k=31 primary with
#     k=21/51 as sensitivity analyses.
# `write_paired_report` therefore prints the axes actually swept and the pair
# count next to every p-value, and labels the analysis exploratory.
#
# It exists for two reasons beyond producing the number.
#
# 1. IT IS THE EXECUTABILITY PROOF FOR BEAD td-59o7. The pairing key INCLUDES
#    `seed`. Run this against a sweep CSV that predates the `seed` column and it
#    raises with a pointer to `rgv_seed_backfill.jl`; run it against the current
#    schema and it produces a p-value. "Column added" is not the deliverable —
#    "the pre-registered test runs" is.
#
# 2. IT IS THE ANALYSIS PATH BEAD td-9p91 GUARDS. Before any aggregation it calls
#    `assert_single_metric_definition`, so a table mixing `metric_source="quast"`
#    with `"internal:quast-failed"` — or mixing two `quast_min_contig` thresholds —
#    RAISES instead of quietly comparing an alignment-validated NGA50 against a
#    size-ratio proxy. Selecting a definition is explicit (`--metric-source quast`)
#    or the run fails.
#
#    A definition column can also be SUPPLIED after the fact by
#    `rgv_seed_backfill.jl --metric-source`, on operator say-so, with nothing to
#    verify it against. The guard cannot see the difference — a backfilled label
#    agrees with itself by construction — so a genuinely mixed legacy table becomes
#    one the guard certifies as uniform. This script therefore reads the
#    `<column>_provenance` column the backfill writes and refuses to print the
#    "guard enforced it" assurance over operator-asserted rows; `results.json`
#    carries the same fact in `metric_definition_operator_asserted`.
#
# HOW PAIRS ARE FORMED
# --------------------
# One pair per (reference, error_rate, readlen, target_coverage, k, seed) cell:
# the treatment arm and the control arm ran on the SAME simulated reads, which is
# what makes the comparison paired. A cell missing either arm, holding a `missing`
# metric, or carrying duplicate rows for one arm is DROPPED AND REPORTED — never
# imputed, never silently averaged.
#
# STATISTICS
# ----------
# Wilcoxon signed-rank, implemented here rather than pulled from a new dependency:
#   * zero differences dropped (Wilcoxon's original treatment), count reported
#   * ties in |difference| get average ranks
#   * EXACT two-sided p from the null rank-sum distribution when there are no ties
#     and no zeros and n <= 50 (dynamic program over subset sums, not enumeration)
#   * otherwise the normal approximation with the standard tie correction
#     sigma^2 = [n(n+1)(2n+1) - sum(t^3 - t)/2] / 24
# Both branches are cross-validated against scipy.stats.wilcoxon in
# `test/4_assembly/rgv_paired_wilcoxon_test.jl`.
# Multiplicity across metrics uses Benjamini-Hochberg, as the pre-registration
# specifies.
#
# USAGE
# -----
#   julia --project=. benchmarking/rgv_paired_wilcoxon.jl \
#       --csv results/rhizomorph_correction_validation_sweep_20260725_*.csv \
#       --metric-source quast \
#       --output-dir benchmarking/results/rgv_paired_wilcoxon
#
#   Flags:
#     --csv PATH               (repeatable) sweep CSV; several are concatenated,
#                              which is the point of having `seed` — shards from
#                              different seeds merge into one pairable table
#     --treatment NAME         arm under test (default "iterative")
#     --control NAME           baseline arm (default "naive")
#     --metrics A,B            metrics to test (default quast_nga50,quast_num_misassemblies)
#     --metric-source VALUE    keep only rows with this metric_source (e.g. "quast")
#     --allow-mixed-src        proceed across mixed metric definitions (loud warning)
#     --keep-not-ok            keep rows with ok=false (default: drop and report)
#     --output-dir PATH        write report.md + results.json here
#
# The figure is produced separately by `rgv_paired_wilcoxon_figure.jl`, which
# reads this script's `results.json`, so the analysis stays fast enough to unit-test.

import CSV
import DataFrames
import Distributions
import JSON
import Dates
import Statistics

# Pure metric-definition guard (bead td-9p91) — the reason this analysis cannot
# silently aggregate across two metric definitions.
include(joinpath(@__DIR__, "metric_source_guard.jl"))

"""
Columns that identify one sweep CELL: the unit within which the treatment and
control arms are paired. `seed` is REQUIRED — it is what makes two replicate rows
distinguishable after shards are merged (bead td-59o7).
"""
const RGV_PAIR_KEYS = (:reference, :error_rate, :readlen, :target_coverage, :k, :seed)

"""
Default metrics, matching the pre-registration's primary metrics for H1: NGA50
(contiguity) and misassembly count (correctness).
"""
const RGV_DEFAULT_METRICS = (:quast_nga50, :quast_num_misassemblies)

"""
Direction in which each metric IMPROVES. `:higher` means larger is better,
`:lower` means smaller is better, `:none` means report the difference but make no
improvement claim (a duplication ratio is best near 1, not extremal).
"""
const RGV_METRIC_DIRECTION = Dict(
    :quast_nga50 => :higher,
    :quast_genome_fraction => :higher,
    :quast_num_misassemblies => :lower,
    :quast_duplication_ratio => :none,
    :n50 => :higher,
    :largest_contig => :higher,
    :genome_fraction => :higher
)

const RGV_ALPHA = 0.05
const RGV_IMPROVEMENT_THRESHOLD = 0.10  # pre-reg's "> 10% median improvement"

# === Wilcoxon signed-rank ===================================================

"""
    average_ranks(x) -> Vector{Float64}

Ranks of `x` with tied values receiving the average of the ranks they span (the
standard treatment for Wilcoxon when |differences| tie).
"""
function average_ranks(x::AbstractVector{<:Real})
    n = length(x)
    perm = sortperm(x)
    ranks = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && x[perm[j + 1]] == x[perm[i]]
            j += 1
        end
        r = (i + j) / 2
        for t in i:j
            ranks[perm[t]] = r
        end
        i = j + 1
    end
    return ranks
end

"""
    exact_signed_rank_pvalue(n, w_plus) -> Float64

Exact two-sided p-value for a signed-rank statistic `w_plus` under H0 with `n`
non-zero differences and no ties.

Under H0 each rank 1..n carries a `+` or `-` sign with probability 1/2
independently, so the null distribution of W+ is the distribution of a random
subset sum of `{1, ..., n}`. A dynamic program over that distribution is exact and
costs `O(n^3)` — far cheaper than the `2^n` enumeration, and it stays exact well
past the point enumeration becomes impractical.
"""
function exact_signed_rank_pvalue(n::Integer, w_plus::Real)
    maxw = n * (n + 1) ÷ 2
    dist = zeros(Float64, maxw + 1)
    dist[1] = 1.0                      # index w+1 holds P(W+ == w)
    for r in 1:n
        nxt = zeros(Float64, maxw + 1)
        for w in 0:maxw
            pw = dist[w + 1]
            pw == 0.0 && continue
            nxt[w + 1] += 0.5 * pw
            (w + r) <= maxw && (nxt[w + r + 1] += 0.5 * pw)
        end
        dist = nxt
    end
    w = Int(round(w_plus))
    w = clamp(w, 0, maxw)
    lower = sum(@view dist[1:(w + 1)])
    upper = sum(@view dist[(w + 1):end])
    return min(1.0, 2 * min(lower, upper))
end

"""
    wilcoxon_signed_rank(diffs; exact_max_n=50) -> NamedTuple

Two-sided Wilcoxon signed-rank test on paired differences `diffs`.

Returns `(; n, n_zero_dropped, w_plus, w_minus, statistic, pvalue, method)` where
`statistic = min(w_plus, w_minus)` (the conventional reported W) and `method`
names which branch produced the p-value.

Zero differences are dropped and counted. The exact branch is used only when it is
valid (no zeros, no ties in |difference|, `n <= exact_max_n`); otherwise the normal
approximation with the standard tie correction is used. With no non-zero
differences the test is undefined and `pvalue` is `NaN` — reported, not silently
treated as a null result.
"""
function wilcoxon_signed_rank(diffs::AbstractVector{<:Real}; exact_max_n::Integer = 50)
    nonzero = Float64[d for d in diffs if d != 0]
    n_zero = length(diffs) - length(nonzero)
    n = length(nonzero)
    if n == 0
        return (n = 0, n_zero_dropped = n_zero, w_plus = NaN, w_minus = NaN,
            statistic = NaN, pvalue = NaN,
            method = "undefined (no non-zero differences)")
    end
    absd = abs.(nonzero)
    ranks = average_ranks(absd)
    w_plus = sum(Float64[ranks[i] for i in 1:n if nonzero[i] > 0]; init = 0.0)
    w_minus = sum(Float64[ranks[i] for i in 1:n if nonzero[i] < 0]; init = 0.0)
    statistic = min(w_plus, w_minus)

    has_ties = length(unique(absd)) < n
    if !has_ties && n_zero == 0 && n <= exact_max_n
        return (n = n, n_zero_dropped = n_zero, w_plus = w_plus, w_minus = w_minus,
            statistic = statistic, pvalue = exact_signed_rank_pvalue(n, w_plus),
            method = "exact (null rank-sum distribution)")
    end

    mu = n * (n + 1) / 4
    tie_term = 0.0
    for t in values(_value_counts(absd))
        tie_term += t^3 - t
    end
    sigma2 = (n * (n + 1) * (2n + 1) - tie_term / 2) / 24
    if sigma2 <= 0
        return (n = n, n_zero_dropped = n_zero, w_plus = w_plus, w_minus = w_minus,
            statistic = statistic, pvalue = NaN,
            method = "undefined (zero variance under H0)")
    end
    z = (w_plus - mu) / sqrt(sigma2)
    p = 2 * Distributions.ccdf(Distributions.Normal(), abs(z))
    return (n = n, n_zero_dropped = n_zero, w_plus = w_plus, w_minus = w_minus,
        statistic = statistic, pvalue = min(1.0, p),
        method = "normal approximation with tie correction (z = $(round(z; digits = 4)))")
end

function _value_counts(x)
    counts = Dict{Float64, Int}()
    for v in x
        counts[Float64(v)] = get(counts, Float64(v), 0) + 1
    end
    return counts
end

"""
    benjamini_hochberg(pvalues) -> Vector{Float64}

Benjamini-Hochberg FDR-adjusted p-values, in the input order. `NaN` inputs are
excluded from the adjustment (an undefined test does not consume multiplicity
budget) and returned as `NaN`.
"""
function benjamini_hochberg(pvalues::AbstractVector{<:Real})
    adj = fill(NaN, length(pvalues))
    valid = [i for i in eachindex(pvalues) if !isnan(pvalues[i])]
    isempty(valid) && return adj
    m = length(valid)
    order = sort(valid; by = i -> pvalues[i])
    prev = 1.0
    for rank in m:-1:1
        idx = order[rank]
        val = min(prev, pvalues[idx] * m / rank)
        adj[idx] = val
        prev = val
    end
    return adj
end

# === Loading + pairing ======================================================

"""
    load_sweep_csvs(paths) -> DataFrames.DataFrame

Read and vertically concatenate sweep CSVs. `cols=:union` is used so schema drift
between files surfaces as `missing` values that the pairing step then rejects,
rather than as a silent read failure.
"""
function load_sweep_csvs(paths::AbstractVector{<:AbstractString})
    isempty(paths) && error("no CSV supplied")
    frames = DataFrames.DataFrame[]
    for p in paths
        isfile(p) || error("CSV not found: $p")
        push!(frames, CSV.read(p, DataFrames.DataFrame))
    end
    return reduce((a, b) -> vcat(a, b; cols = :union), frames)
end

const RGV_BACKFILL_HINT = "  julia --project=. benchmarking/rgv_seed_backfill.jl --csv <csv> " *
                          "(--seed N | --run-log PATH | --assume-default-seed) " *
                          "[--metric-source <src> --quast-min-contig <bp>]"

"""
    require_pairing_schema(df)

Assert the table can be paired at all: every pairing key must be PRESENT **and
fully populated**.

Checking presence alone is not enough, and the gap was a real defect. Shards are
concatenated with `cols = :union`, so a shard predating the `seed` column
contributes rows whose `seed` is `missing` — the column is then present and a
presence-only check passes. `build_pairs` stringifies the key, so `missing`
becomes the literal pseudo-seed `"missing"`, a distinct cell that pairs normally.
The same physical runs supplied twice then report n = 12 for 6 independent
replicates and deflate p by ~27x. Pseudo-replication is the single worst failure
available to this file, because it makes the result look BETTER, so the key
columns are checked for missing VALUES, not just for existence.
"""
function require_pairing_schema(df)
    cols = Set(Symbol.(DataFrames.names(df)))
    missing_keys = [k for k in RGV_PAIR_KEYS if !(k in cols)]
    if :seed in missing_keys
        error("this sweep CSV has no `seed` column, so replicate rows cannot be " *
              "paired and the pre-registration's paired-Wilcoxon rule over seeds " *
              "42/123/456 is NOT runnable against it (bead td-59o7). Backfill it " *
              "first:\n" * RGV_BACKFILL_HINT)
    end
    if !isempty(missing_keys)
        error("sweep CSV is missing pairing column(s): " *
              join(("`$k`" for k in missing_keys), ", "))
    end
    :arm in cols || error("sweep CSV is missing the `arm` column")

    # Fully-populated check. A `missing` in any pairing key means some shard did
    # not carry that column; `cols = :union` filled it in. Refuse rather than let
    # `missing` become a pseudo-level of that key.
    for k in RGV_PAIR_KEYS
        col = getproperty(df, k)
        n_missing = count(ismissing, col)
        n_missing == 0 && continue
        error("$(n_missing) of $(length(col)) rows have a `missing` value in the " *
              "pairing column `$k`. That happens when shards with different schemas " *
              "are concatenated (`cols = :union` fills the absent column with " *
              "`missing`), and it is NOT safe to proceed: `missing` would become a " *
              "distinct pseudo-level of `$k`, so the same physical runs supplied " *
              "twice would be counted as independent replicates — inflating n and " *
              "deflating p. Backfill the incomplete shard before merging it:\n" *
              RGV_BACKFILL_HINT)
    end

    # `arm` is not a pairing key but it selects which side of the pair a row is on;
    # a `missing` arm silently belongs to neither and would be dropped invisibly.
    n_missing_arm = count(ismissing, getproperty(df, :arm))
    n_missing_arm == 0 ||
        error("$(n_missing_arm) row(s) have a `missing` `arm` and cannot be " *
              "assigned to either side of a pair")
    return nothing
end

"""
    build_pairs(df, metric; treatment, control) -> (pairs, dropped)

Form one paired observation per cell for `metric`.

`pairs` is a vector of `(; key, control_value, treatment_value, difference)` where
`difference = treatment - control`. `dropped` is a vector of
`(; key, reason)` for every cell that could not contribute — missing arm,
duplicate rows for one arm, or a `missing`/non-finite metric on either side. Both
are returned so the report can state exactly how many cells the test actually used.
"""
function build_pairs(df, metric::Symbol; treatment::AbstractString, control::AbstractString)
    metric in Symbol.(DataFrames.names(df)) ||
        error("metric column `$metric` not present in the sweep table")
    keyvals = [getproperty(df, k) for k in RGV_PAIR_KEYS]
    arms = getproperty(df, :arm)
    metvals = getproperty(df, metric)
    nrows = length(arms)

    order = String[]
    cells = Dict{String, Dict{String, Vector{Any}}}()
    for i in 1:nrows
        key = join((string(v[i]) for v in keyvals), " | ")
        if !haskey(cells, key)
            cells[key] = Dict{String, Vector{Any}}()
            push!(order, key)
        end
        push!(get!(cells[key], string(arms[i]), Any[]), metvals[i])
    end

    pairs = Any[]
    dropped = Any[]
    for key in order
        byarm = cells[key]
        if !haskey(byarm, control) || !haskey(byarm, treatment)
            missing_arm = !haskey(byarm, control) ? control : treatment
            push!(dropped, (key = key, reason = "no `$missing_arm` row for this cell"))
            continue
        end
        if length(byarm[control]) != 1 || length(byarm[treatment]) != 1
            push!(dropped,
                (key = key,
                    reason = "duplicate rows for one arm ($(length(byarm[control])) " *
                             "$control, $(length(byarm[treatment])) $treatment) — the " *
                             "pairing key does not uniquely identify a cell"))
            continue
        end
        cv, tv = first(byarm[control]), first(byarm[treatment])
        if cv === missing || tv === missing
            push!(dropped, (key = key, reason = "`missing` $metric on at least one arm"))
            continue
        end
        cvf, tvf = Float64(cv), Float64(tv)
        if !isfinite(cvf) || !isfinite(tvf)
            push!(dropped, (key = key, reason = "non-finite $metric on at least one arm"))
            continue
        end
        push!(pairs,
            (key = key, control_value = cvf, treatment_value = tvf,
                difference = tvf - cvf))
    end
    return pairs, dropped
end

"""
    analyze_metric(df, metric; treatment, control) -> NamedTuple

Pair, test, and size the effect for one metric. Reports the median paired
difference (the pre-registration's effect size) AND the median per-pair relative
improvement (what its > 10% decision threshold is stated in), oriented by
`RGV_METRIC_DIRECTION` so "improvement" is positive for both higher-is-better and
lower-is-better metrics.

Pairs whose control value is zero are excluded from the RELATIVE improvement only
(the ratio is undefined) and counted, since dropping them can bias that statistic.
"""
function analyze_metric(df, metric::Symbol;
        treatment::AbstractString, control::AbstractString)
    pairs, dropped = build_pairs(df, metric; treatment = treatment, control = control)
    diffs = Float64[p.difference for p in pairs]
    test = wilcoxon_signed_rank(diffs)
    direction = get(RGV_METRIC_DIRECTION, metric, :none)
    sign_factor = direction == :lower ? -1.0 : 1.0

    median_diff = isempty(diffs) ? NaN : Statistics.median(diffs)
    rel = Float64[]
    n_zero_control = 0
    for p in pairs
        if p.control_value == 0
            n_zero_control += 1
            continue
        end
        push!(rel, sign_factor * (p.treatment_value - p.control_value) /
                   abs(p.control_value))
    end
    median_rel = isempty(rel) ? NaN : Statistics.median(rel)

    return (metric = metric, direction = direction,
        n_pairs = length(pairs), n_dropped = length(dropped), dropped = dropped,
        pairs = pairs,
        median_paired_difference = median_diff,
        median_relative_improvement = median_rel,
        n_zero_control_excluded_from_relative = n_zero_control,
        test = test)
end

"""
    decide(result, adjusted_p) -> String

The pre-registration's decision rule, applied to one metric.

Four outcomes, not three. The rule's falsification clause is
`"Falsified if p >= 0.05, THE DIRECTION IS NEGATIVE, or ..."`, so a
statistically significant result in the WRONG direction is falsification, not
weak support. An earlier version branched only on
`improved > RGV_IMPROVEMENT_THRESHOLD`, which made every significant result that
failed the 10% bar — including outright harm — read as "PARTIALLY SUPPORTED":
a treatment worse on all 9 pairs, median -45%, FDR p = 0.0039, was reported as
partial support. That is the one error class this whole PR exists to prevent, so
the sign is now tested before the magnitude.

`:none`-direction metrics get "reported (no directional claim)" rather than a
verdict — the pre-registration defines no improvement threshold for them.

An undefined effect size (`NaN`, e.g. every control value zero so no relative
improvement is computable) is reported as indeterminate rather than being allowed
to fall through the `> threshold` comparison, where `NaN > 0.1 === false` would
have silently rendered it as "<= 10%".
"""
function decide(result, adjusted_p)
    if result.direction == :none
        return "reported (no directional claim pre-registered for this metric)"
    end
    if result.n_pairs == 0 || isnan(adjusted_p)
        return "indeterminate (test undefined — see n_pairs / dropped cells)"
    end
    improved = result.median_relative_improvement
    pct = round(Int, RGV_IMPROVEMENT_THRESHOLD * 100)
    if adjusted_p >= RGV_ALPHA
        return "NOT SUPPORTED (FDR-adjusted p >= $RGV_ALPHA — report as a null with effect size)"
    end
    # Significant. The sign decides between falsified and (partially) supported.
    if isnan(improved)
        return "INDETERMINATE (significant, but the effect size is undefined — " *
               "every control value is zero, so relative improvement cannot be " *
               "computed; read the median paired difference instead)"
    end
    if improved < 0
        return "FALSIFIED (FDR-adjusted p < $RGV_ALPHA and the direction is " *
               "NEGATIVE — the treatment is significantly WORSE by " *
               "$(round(abs(improved) * 100; digits = 1))%)"
    end
    if improved > RGV_IMPROVEMENT_THRESHOLD
        return "SUPPORTED (FDR-adjusted p < $RGV_ALPHA and median improvement > $pct%)"
    end
    return "PARTIALLY SUPPORTED (significant, direction positive, but median " *
           "improvement <= $pct%)"
end

"""
Substring identifying a metric definition that was ASSERTED BY AN OPERATOR rather
than observed by the run that produced the metrics.

`rgv_seed_backfill.jl` writes `"operator-asserted (backfill)"` into
`<definition-column>_provenance` for every row it supplies a definition for. The
match is on the substring so the writer can extend the label (with a tool name, a
date) without breaking the reader; the coupling is pinned by a BEHAVIOURAL
round-trip test (backfill writes a frame, this reader classifies it), because a
grep of the writer's source passes unchanged while its call site drifts to a label
this reader does not understand.
"""
const RGV_ASSERTED_DEFINITION_MARKER = "operator-asserted"

"""
Label identifying a metric definition the run itself OBSERVED, written by
`rgv_seed_backfill.jl` as `DEFINITION_PROVENANCE_OBSERVED` for rows it did not
supply a value for.

Matched EXACTLY, not as a substring: this is the only value that earns the
unqualified assurance, so anything else — including a label a future writer
invents — must fall through to `undeterminable` rather than be read as evidence.
The coupling to the writer is pinned by a behavioural round-trip test, not by a
grep of the writer's source (a grep passes while the writer's call site drifts).
"""
const RGV_OBSERVED_DEFINITION_LABEL = "observed"

"""
    definition_provenance_axes(df) -> (asserted, undeterminable, n_asserted_rows)

Classify the provenance of every metric-definition axis PRESENT in `df` as
asserted, observed, or undeterminable, and count the rows carrying an assertion.

This is the fact `assert_single_metric_definition` structurally CANNOT establish.
The guard checks that the definition labels agree with each other; when those
labels were written by `rgv_seed_backfill.jl --metric-source`, agreement is
guaranteed by construction and says nothing about how the runs were scored. A
legacy table that genuinely mixed QUAST-scored and degraded rows, backfilled with
a single `--metric-source`, passes the guard perfectly. Only the provenance column
distinguishes that case, so the report and `results.json` must surface it.

# Why three states and not two

The two-state reader FAILED OPEN in exactly the situation it exists for. Absence
of the `"operator-asserted"` marker was read as evidence of measurement, so an
axis whose sibling `<axis>_provenance` column is missing — which is EVERY CSV the
documented backfill/analyse workflow produced before this axis existed — was
reported as `metric_definition_operator_asserted: false`, an affirmative claim the
definition was measured, and the report went on to print the unqualified
"the guard enforced it — no override was used".

An unrecognised label failed the same way and is worse: a column present but
holding, say, `"backfilled-by-tool"` is STRONGER evidence that something wrote
provenance metadata than no column at all, yet a blacklist reader treats it more
permissively. Both are the same epistemic state — the table carries no usable
provenance — and neither is "observed". This mirrors the doctrine applied twelve
lines away to `ok`, and the one `metric_source_guard.jl` applies when one
definition axis is present and another absent: unverifiable is not verified.

An axis can appear in BOTH returned lists (some rows asserted, others carrying an
unreadable label); each list answers its own question.
"""
function definition_provenance_axes(df)
    cols = Set(String.(DataFrames.names(df)))
    asserted = String[]
    undeterminable = String[]
    asserted_rows = falses(DataFrames.nrow(df))
    for col in METRIC_DEFINITION_COLUMNS
        # Only axes the table actually carries are classified: an axis that is not
        # in the table is not part of the definition being reported on.
        String(col) in cols || continue
        prov = String(col) * "_provenance"
        if !(prov in cols)
            push!(undeterminable, String(col))
            continue
        end
        values = getproperty(df, Symbol(prov))
        flags = Bool[!ismissing(v) &&
                     occursin(RGV_ASSERTED_DEFINITION_MARKER, string(v))
                     for v in values]
        unreadable = Bool[!flags[i] &&
                          (ismissing(values[i]) ||
                           string(values[i]) != RGV_OBSERVED_DEFINITION_LABEL)
                          for i in eachindex(values)]
        if any(flags)
            push!(asserted, String(col))
            asserted_rows .|= flags
        end
        any(unreadable) && push!(undeterminable, String(col))
    end
    return asserted, undeterminable, count(asserted_rows)
end

"""
    operator_asserted_definitions(df) -> (axes, n_rows)

Definition axes whose value was operator-asserted for at least one row of `df`,
and the number of rows carrying an assertion on any axis.

Thin view over [`definition_provenance_axes`](@ref) that answers only the
asserted question. Callers that need to distinguish "observed" from "cannot be
determined" must use the three-state function.
"""
function operator_asserted_definitions(df)
    asserted, _, n_rows = definition_provenance_axes(df)
    return asserted, n_rows
end

"""
    run_paired_analysis(df; treatment, control, metrics, metric_source=nothing,
                        allow_mixed_src=false, keep_not_ok=false) -> NamedTuple

Full pre-registered analysis over one already-loaded sweep table.

Order of operations is deliberate: schema check, then row filtering, then the
metric-definition guard, then pairing and testing. The guard runs on the FILTERED
rows — the ones that will actually be aggregated — so selecting
`metric_source="quast"` is what makes the run legal, not a bypass of the check.
"""
function run_paired_analysis(df;
        treatment::AbstractString = "iterative",
        control::AbstractString = "naive",
        # Any iterable of Symbol: a Tuple (the default) or a Vector (from --metrics).
        metrics = RGV_DEFAULT_METRICS,
        metric_source::Union{Nothing, AbstractString} = nothing,
        allow_mixed_src::Bool = false,
        keep_not_ok::Bool = false)
    require_pairing_schema(df)

    work = df
    n_input = DataFrames.nrow(work)
    # "The ok check ran and every row passed" and "the ok check COULD NOT RUN"
    # are different facts about the data, and they used to render identically:
    # `0 dropped for ok=false` was emitted verbatim whether the column was present
    # and all-true or absent entirely, and `results.json` carried neither count nor
    # a presence flag, so no consumer could recover the distinction from either
    # artifact. `metric_source_guard.jl` already applies the opposite doctrine to
    # the definition columns — a table with no provenance column is not a table
    # shown to be consistent, it is a table that cannot be checked — and `ok` is
    # the same kind of claim.
    ok_column_present = :ok in Symbol.(DataFrames.names(work))
    ok_filter_applied = !keep_not_ok && ok_column_present
    n_dropped_not_ok = 0
    n_dropped_ok_missing = 0
    if ok_filter_applied
        okcol = work.ok
        keep = coalesce.(okcol, false) .== true
        # `coalesce(missing, false)` drops an unknown-status row, which is the right
        # default, but counting it as ok=false conflates "the harness SAID this run
        # failed" with "we do not know whether it succeeded". Only the first is an
        # observed failure, so the two are counted separately.
        n_dropped_ok_missing = count(ismissing, okcol)
        n_dropped_not_ok = count(!, keep) - n_dropped_ok_missing
        work = work[keep, :]
    end
    n_dropped_source = 0
    if metric_source !== nothing
        before = DataFrames.nrow(work)
        work = work[string.(work.metric_source) .== String(metric_source), :]
        n_dropped_source = before - DataFrames.nrow(work)
    end
    DataFrames.nrow(work) == 0 &&
        error("no rows left after filtering (ok=$(!keep_not_ok ? "true" : "any"), " *
              "metric_source=$(something(metric_source, "any")))")

    # Bead td-9p91: refuse to aggregate rows scored under two definitions.
    # The override must leave a trace in the DELIVERABLE, not only in a stderr
    # @warn that lands in a SLURM .out nobody reads: `report.md` / `results.json`
    # are what reach a manuscript. Record both whether the override was PASSED and
    # whether it actually BOUND (i.e. the table really was mixed), because
    # "override available but unused" and "override load-bearing" are very
    # different provenance claims.
    definition = assert_single_metric_definition(work;
        context = "RGV paired-Wilcoxon analysis",
        allow_mixed_src = allow_mixed_src)
    mixed_axes = [String(col) for (col, vals) in definition if length(vals) > 1]
    override_bound = allow_mixed_src && !isempty(mixed_axes)

    # A definition the guard CHECKED is not the same as a definition anyone
    # OBSERVED. `rgv_seed_backfill.jl --metric-source` can write a uniform label
    # onto a table that really did mix definitions, which the guard then certifies
    # as consistent. That fact has to reach the deliverable, for the same reason
    # `override_bound` does.
    asserted_axes, undeterminable_axes,
    n_asserted_rows = definition_provenance_axes(work)
    definition_asserted = !isempty(asserted_axes)
    definition_undeterminable = !isempty(undeterminable_axes)

    results = [analyze_metric(work, m; treatment = treatment, control = control)
               for m in metrics]
    adjusted = benjamini_hochberg(Float64[r.test.pvalue for r in results])
    verdicts = [decide(results[i], adjusted[i]) for i in eachindex(results)]

    seeds = sort(unique(work.seed))
    # The axes actually swept, so the report can state them next to every p-value
    # (the pre-registration's thresholds were calibrated for a different design).
    axes = Dict{Symbol, Vector{Any}}()
    for col in (:reference, :error_rate, :readlen, :target_coverage, :k)
        hasproperty(work, col) || continue
        axes[col] = sort(unique(getproperty(work, col)))
    end
    return (results = results, adjusted_p = adjusted, verdicts = verdicts,
        definition = definition, seeds = seeds,
        axes = axes, allow_mixed_src = allow_mixed_src, mixed_axes = mixed_axes,
        override_bound = override_bound,
        operator_asserted_axes = asserted_axes,
        n_rows_operator_asserted = n_asserted_rows,
        definition_operator_asserted = definition_asserted,
        undeterminable_definition_axes = undeterminable_axes,
        definition_provenance_undeterminable = definition_undeterminable,
        treatment = String(treatment), control = String(control),
        n_input_rows = n_input, n_rows_analyzed = DataFrames.nrow(work),
        ok_column_present = ok_column_present, ok_filter_applied = ok_filter_applied,
        keep_not_ok = keep_not_ok,
        n_dropped_not_ok = n_dropped_not_ok,
        n_dropped_ok_missing = n_dropped_ok_missing,
        n_dropped_metric_source = n_dropped_source,
        pair_keys = RGV_PAIR_KEYS)
end

# === Reporting ==============================================================

_pw_fmt(x) = x isa Real ? (isnan(x) ? "n/a" : string(round(x; digits = 6))) : string(x)

"""
    write_paired_report(path, analysis; csv_paths) -> String

Markdown report: what was analyzed, the metric definition it was analyzed under,
the per-metric test and effect size, the pre-registration verdict, and every
dropped cell with its reason.
"""
function write_paired_report(path::AbstractString, analysis; csv_paths)
    open(path, "w") do io
        println(io, "# RGV correction sweep — paired Wilcoxon (EXPLORATORY)\n")
        println(io, "Generated: $(Dates.now())\n")
        println(io,
            "> **Status: exploratory, not a pre-registered test.** This applies the " *
            "statistical RULE from `PLAN-2026-03-28-rhizomorph-preregistration.md` " *
            "(paired Wilcoxon signed-rank over seeds, Benjamini-Hochberg across " *
            "metrics, and its decision thresholds) to a comparison the " *
            "pre-registration does not describe. Its H1 is *Viterbi DP vs greedy*; " *
            "this compares `corrector=:none` vs `corrector=:iterative`, and no " *
            "hypothesis H1-H7 covers error correction. Do not report these as " *
            "confirmatory results.\n")
        println(io, "## Inputs\n")
        for p in csv_paths
            println(io, "- `$p`")
        end
        println(io,
            "\n- Arms: **$(analysis.treatment)** (treatment) vs " *
            "**$(analysis.control)** (control)")
        println(io, "- Pairing key: " *
                    join(("`$k`" for k in analysis.pair_keys), " x "))
        println(io, "- Seeds present: $(join(string.(analysis.seeds), ", "))")
        # The axes ACTUALLY swept, next to the numbers. A p-value read without them
        # invites the reader to assume the pre-registered design produced it.
        for (label, col) in (("Error rates", :error_rate), ("Read lengths", :readlen),
            ("Coverages", :target_coverage), ("k", :k), ("References", :reference))
            haskey(analysis.axes, col) || continue
            println(io, "- $label swept: " * join(string.(analysis.axes[col]), ", "))
        end
        println(io,
            "\n_The pre-registration's design is 4 coverage levels x 3 replicates " *
            "per organism/technology (~360 paired comparisons), varies coverage " *
            "rather than error rate, and designates k=31 primary (k=21/51 " *
            "sensitivity). Compare the axes above against that before reading the " *
            "decision column — the thresholds were calibrated for that design._")
        # An absent `ok` column must never render as a check that ran and passed.
        ok_clause = if !analysis.ok_column_present
            "**`ok` column ABSENT — no row was verified; the ok=false check could " *
            "NOT run**"
        elseif !analysis.ok_filter_applied
            "`ok` filter DISABLED via --keep-not-ok — failed rows were KEPT"
        else
            "$(analysis.n_dropped_not_ok) dropped for ok=false, " *
            "$(analysis.n_dropped_ok_missing) dropped for ok=missing (status " *
            "unknown, not an observed failure)"
        end
        println(io,
            "- Rows: $(analysis.n_input_rows) read, " *
            "$(analysis.n_rows_analyzed) analyzed " *
            "($ok_clause; " *
            "$(analysis.n_dropped_metric_source) dropped by metric_source filter)")
        println(io, "\n### Metric definition analyzed under\n")
        for (col, vals) in analysis.definition
            println(io, "- `$col`: " * join(("`$v`" for v in vals), ", "))
        end
        # Same doctrine as the override block below, applied to the other way the
        # assurance can be false. A definition SUPPLIED by `rgv_seed_backfill.jl`
        # is uniform by construction, so the guard passing over it is not evidence
        # about how the runs were scored — a legacy table that genuinely mixed
        # QUAST-scored and degraded rows, backfilled with one `--metric-source`,
        # produces exactly this artifact. The reader must be told before reading
        # the numbers, and the unqualified assurance below is suppressed.
        if analysis.definition_operator_asserted
            println(io,
                "\n> ## :warning: METRIC DEFINITION OPERATOR-ASSERTED, NOT OBSERVED\n>\n" *
                "> $(analysis.n_rows_operator_asserted) of " *
                "$(analysis.n_rows_analyzed) analyzed rows carry a definition that " *
                "was **supplied by `rgv_seed_backfill.jl`** on operator say-so, on " *
                join(("`$a`" for a in analysis.operator_asserted_axes), ", ") *
                " (see the `*_provenance` column(s) in the input CSV).\n>\n" *
                "> Nothing verified that claim: unlike `seed`, which can be " *
                "recovered from a run log, there is no record of which definition " *
                "produced a historical metric. `metric_source_guard.jl` can only " *
                "confirm that the asserted labels agree with EACH OTHER, and a " *
                "backfilled label agrees with itself by construction — a table that " *
                "genuinely mixed an alignment-validated QUAST metric with an " *
                "internal size-ratio proxy passes this check once it is backfilled " *
                "with a single value. **Treat the definition below as an operator " *
                "assertion, not as measured provenance.**\n")
        end
        # The third state. "No contrary signal" is not evidence of measurement, and
        # reading it as such is how the two-state reader certified every pre-existing
        # sweep CSV — none of which has a `*_provenance` column — as measured.
        if analysis.definition_provenance_undeterminable
            println(io,
                "\n> ## :warning: METRIC DEFINITION PROVENANCE UNDETERMINABLE\n>\n" *
                "> This table carries no usable provenance for " *
                join(("`$a`" for a in analysis.undeterminable_definition_axes), ", ") *
                ": the sibling `<axis>_provenance` column is missing, empty, or " *
                "holds a label this reader does not recognise (only " *
                "`$(RGV_OBSERVED_DEFINITION_LABEL)` and labels containing " *
                "`$(RGV_ASSERTED_DEFINITION_MARKER)` are understood).\n>\n" *
                "> Whether those definitions were MEASURED by the runs or ASSERTED " *
                "by an operator therefore **cannot be determined from this table**. " *
                "That is not the same as their having been measured, and it is the " *
                "state every sweep CSV written before the provenance column existed " *
                "is in. Recover it with `rgv_seed_backfill.jl`, or read the numbers " *
                "below as resting on an unverified definition.\n")
        end
        # The artifact must never assert the guard was enforced when it was
        # overridden. Previously this line printed unconditionally, so a run with
        # `--allow-mixed-src` listed two definitions and then, one line later, told
        # the reader single-definition-ness had been enforced. That is strictly
        # worse than having no guard: it launders the override into an assurance.
        if analysis.override_bound
            println(io,
                "\n> ## :warning: MIXED METRIC DEFINITIONS — GUARD OVERRIDDEN\n>\n" *
                "> This analysis was run with `--allow-mixed-src`, and the override " *
                "was **load-bearing**: the rows really do span more than one " *
                "definition on " *
                join(("`$a`" for a in analysis.mixed_axes), ", ") * ".\n>\n" *
                "> Every number below therefore compares quantities produced under " *
                "different definitions — for `metric_source`, an alignment-validated " *
                "QUAST metric against an internal size-ratio proxy that is blind to " *
                "misassembly. **These results are not validation-grade and must not " *
                "be reported as a measured effect.**\n")
        elseif analysis.allow_mixed_src
            # Deliberately does NOT say "so the aggregate is well-defined": that is
            # a claim about the VALUES, and it is settled by the provenance block
            # below, not by whether a flag bound. This branch reports only what it
            # actually knows — that the override did not fire.
            println(io,
                "\n_`--allow-mixed-src` was passed but did NOT bind: the rows are " *
                "single-definition on every checked axis._\n")
        elseif !analysis.definition_operator_asserted &&
               !analysis.definition_provenance_undeterminable
            # The unqualified assurance requires AFFIRMATIVELY-OBSERVED provenance on
            # every present axis, not merely the absence of a contrary signal. Gating
            # it on `!definition_operator_asserted` alone handed it to every table
            # with no provenance column at all.
            println(io,
                "\n_A single value per definition column is what makes the aggregate " *
                "meaningful; `metric_source_guard.jl` (bead td-9p91) enforced it — " *
                "no override was used._\n")
        end

        # INDEPENDENT of the override state above. Whether the guard was overridden
        # and whether the values it compared were observed are two orthogonal facts,
        # and chaining them through one `if/elseif` let the weaker one shadow the
        # stronger: a harmless NON-BINDING `--allow-mixed-src` suppressed this
        # qualifier and restored an unqualified "the aggregate is well-defined" over
        # a table whose definition was never observed — the exact assurance this
        # machinery exists to withhold. Emitted last so it is the final word.
        if analysis.definition_operator_asserted
            println(io,
                "\n_`metric_source_guard.jl` (bead td-9p91) checked what it was given, " *
                "but the values it compared were operator-asserted for some rows " *
                "(see the warning above), so this run does NOT establish that the " *
                "aggregate spans a single MEASURED definition._\n")
        elseif analysis.definition_provenance_undeterminable
            println(io,
                "\n_`metric_source_guard.jl` (bead td-9p91) checked what it was given, " *
                "but the table carries no usable provenance for the axes named above " *
                "(see the warning above), so this run does NOT establish that the " *
                "aggregate spans a single MEASURED definition._\n")
        end

        println(io, "## Results\n")
        println(io,
            "| metric | better | n pairs | dropped | W | p | FDR p | median paired diff | median % improvement | verdict |")
        println(io, "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for (i, r) in enumerate(analysis.results)
            pct = isnan(r.median_relative_improvement) ? "n/a" :
                  string(round(r.median_relative_improvement * 100; digits = 2)) * "%"
            println(io,
                "| `$(r.metric)` | $(r.direction) | $(r.n_pairs) | $(r.n_dropped) | " *
                "$(_pw_fmt(r.test.statistic)) | $(_pw_fmt(r.test.pvalue)) | " *
                "$(_pw_fmt(analysis.adjusted_p[i])) | $(_pw_fmt(r.median_paired_difference)) | " *
                "$pct | $(analysis.verdicts[i]) |")
        end
        println(io)
        for r in analysis.results
            println(io, "### `$(r.metric)`\n")
            println(io, "- Test: $(r.test.method)")
            println(io,
                "- Non-zero differences used: $(r.test.n) " *
                "($(r.test.n_zero_dropped) zero difference(s) dropped)")
            println(io, "- W+ = $(_pw_fmt(r.test.w_plus)), W- = $(_pw_fmt(r.test.w_minus))")
            if r.n_zero_control_excluded_from_relative > 0
                println(io,
                    "- **$(r.n_zero_control_excluded_from_relative) pair(s) " *
                    "excluded from the % improvement statistic** (control value " *
                    "is 0, ratio undefined). The absolute median paired " *
                    "difference uses all pairs.")
            end
            if !isempty(r.dropped)
                println(io, "\n#### Cells dropped ($(length(r.dropped)))\n")
                for d in r.dropped
                    println(io, "- `$(d.key)` — $(d.reason)")
                end
            end
            println(io)
        end
        # The legend must enumerate EVERY outcome `decide` can emit, or a reader
        # meets a verdict string the report never defines. It listed three while
        # the function emitted five.
        pct = round(Int, RGV_IMPROVEMENT_THRESHOLD * 100)
        println(io, "## Decision rule applied\n")
        println(io,
            "- **SUPPORTED**: FDR-adjusted p < $RGV_ALPHA, direction positive, " *
            "median improvement > $pct%")
        println(io,
            "- **PARTIALLY SUPPORTED**: FDR-adjusted p < $RGV_ALPHA, direction " *
            "positive, median improvement <= $pct%")
        println(io,
            "- **FALSIFIED**: FDR-adjusted p < $RGV_ALPHA and the direction is " *
            "NEGATIVE — the treatment is significantly worse. This is the " *
            "pre-registration's falsification clause: a significant harm is not " *
            "weak support.")
        println(io,
            "- **NOT SUPPORTED**: FDR-adjusted p >= $RGV_ALPHA (report as a null " *
            "result with effect size)")
        println(io,
            "- **INDETERMINATE**: the test or the effect size is undefined (no " *
            "non-zero differences, or every control value is zero so relative " *
            "improvement cannot be computed) — read the median paired difference " *
            "instead")
    end
    return path
end

"""
    _json_num(x)

Render a number for JSON, mapping non-finite values to `nothing` (`null`).

`NaN` is an ORDINARY outcome here — an undefined test, or an undefined effect
size when every control value is zero — so it reaches the payload on normal runs.
JSON has no NaN literal. Julia's JSON.jl happens to emit `null` for it rather
than raising, but relying on that leaves the intent implicit and version-dependent;
doing it explicitly keeps the contract with `rgv_paired_wilcoxon_figure.jl`, which
already treats `null` as missing, visible in this file.
"""
_json_num(x::Real) = isfinite(x) ? x : nothing
_json_num(x) = x

"""
    paired_analysis_json(analysis; csv_paths) -> Dict

Machine-readable form of the analysis, also consumed by
`rgv_paired_wilcoxon_figure.jl` to draw the paired-difference figure.
"""
function paired_analysis_json(analysis; csv_paths)
    return Dict(
        "generated_at" => string(Dates.now()),
        "inputs" => String[String(p) for p in csv_paths],
        "treatment" => analysis.treatment,
        "control" => analysis.control,
        "pair_keys" => String[string(k) for k in analysis.pair_keys],
        "seeds" => collect(analysis.seeds),
        "n_input_rows" => analysis.n_input_rows,
        "n_rows_analyzed" => analysis.n_rows_analyzed,
        # Whether the ok=false check RAN, not only what it dropped. Without the
        # flag, "0 dropped" is indistinguishable from "could not be checked".
        "ok_column_present" => analysis.ok_column_present,
        "ok_filter_applied" => analysis.ok_filter_applied,
        "n_dropped_not_ok" => analysis.n_dropped_not_ok,
        "n_dropped_ok_missing" => analysis.n_dropped_ok_missing,
        "n_dropped_metric_source" => analysis.n_dropped_metric_source,
        "metric_definition" => Dict(string(col) => vals
        for (col, vals) in analysis.definition),
        "exploratory" => true,
        "preregistration_status" => "Applies the pre-registration's statistical rule to a comparison it " *
                                    "does not describe (H1 is Viterbi DP vs greedy; this is " *
                                    "corrector=:none vs :iterative). Not a confirmatory test.",
        "axes_swept" => Dict(string(k) => v for (k, v) in analysis.axes),
        "allow_mixed_src" => analysis.allow_mixed_src,
        "mixed_definition_axes" => analysis.mixed_axes,
        "metric_definition_override_bound" => analysis.override_bound,
        # Whether the definition the guard checked was MEASURED or merely ASSERTED.
        # A consumer that reads only `metric_definition` cannot tell the two apart,
        # and they are very different provenance claims.
        "metric_definition_operator_asserted" => analysis.definition_operator_asserted,
        "operator_asserted_definition_axes" => analysis.operator_asserted_axes,
        "n_rows_operator_asserted_definition" => analysis.n_rows_operator_asserted,
        # The THIRD state, kept separate from the boolean above rather than folded
        # into it: `metric_definition_operator_asserted` still means "at least one
        # row was asserted", and `false` on it no longer implies "observed".
        "metric_definition_provenance_undeterminable" => analysis.definition_provenance_undeterminable,
        "undeterminable_definition_provenance_axes" => analysis.undeterminable_definition_axes,
        "alpha" => RGV_ALPHA,
        "improvement_threshold" => RGV_IMPROVEMENT_THRESHOLD,
        "metrics" => [Dict(
                          "metric" => string(r.metric),
                          "direction" => string(r.direction),
                          "n_pairs" => r.n_pairs,
                          "n_dropped" => r.n_dropped,
                          "n_nonzero_differences" => r.test.n,
                          "n_zero_differences_dropped" => r.test.n_zero_dropped,
                          "w_plus" => _json_num(r.test.w_plus),
                          "w_minus" => _json_num(r.test.w_minus),
                          "statistic" => _json_num(r.test.statistic),
                          "pvalue" => _json_num(r.test.pvalue),
                          "pvalue_fdr" => _json_num(analysis.adjusted_p[i]),
                          "method" => r.test.method,
                          "median_paired_difference" => _json_num(r.median_paired_difference),
                          "median_relative_improvement" => _json_num(r.median_relative_improvement),
                          "n_zero_control_excluded_from_relative" => r.n_zero_control_excluded_from_relative,
                          "verdict" => analysis.verdicts[i],
                          "dropped" => [Dict("key" => d.key, "reason" => d.reason)
                                        for d in r.dropped],
                          # Per-pair observations, so rgv_paired_wilcoxon_figure.jl can
                          # draw the paired differences without re-reading the CSVs.
                          "pairs" => [Dict("key" => pr.key,
                                          "control" => _json_num(pr.control_value),
                                          "treatment" => _json_num(pr.treatment_value),
                                          "difference" => _json_num(pr.difference))
                                      for pr in r.pairs]
                      )
                      for (i, r) in enumerate(analysis.results)]
    )
end

# === CLI ====================================================================

function _pw_args(flag)
    # Consume EVERY token after the flag until the next `--flag`, not just one, so
    # both documented forms work:
    #   --csv a.csv --csv b.csv       (repeated flag)
    #   --csv results/sweep_*.csv     (shell glob -> --csv a.csv b.csv c.csv)
    # Taking only `ARGS[i + 1]` silently used the first value and discarded the
    # rest. Merging per-shard files is the entire purpose of the `seed` column, so
    # that failure mode quietly analysed a fraction of the data behind an
    # invocation this file's own header documents.
    out = String[]
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == flag
            j = i + 1
            while j <= length(ARGS) && !startswith(ARGS[j], "--")
                push!(out, ARGS[j])
                j += 1
            end
            i = j
        else
            i += 1
        end
    end
    return out
end

function _pw_arg(flag)
    i = findfirst(==(flag), ARGS)
    return (i !== nothing && i < length(ARGS)) ? ARGS[i + 1] : nothing
end

function main()
    csv_paths = _pw_args("--csv")
    if isempty(csv_paths)
        println(stderr,
            "ERROR: at least one --csv PATH is required.\n" *
            "  julia --project=. benchmarking/rgv_paired_wilcoxon.jl --csv <sweep.csv> " *
            "[--csv <sweep2.csv>] --metric-source quast")
        return 2
    end
    treatment = something(_pw_arg("--treatment"), "iterative")
    control = something(_pw_arg("--control"), "naive")
    metrics = let v = _pw_arg("--metrics")
        v === nothing ? RGV_DEFAULT_METRICS : Symbol.(strip.(split(v, ",")))
    end
    output_dir = something(_pw_arg("--output-dir"),
        joinpath(@__DIR__, "results", "rgv_paired_wilcoxon"))

    println("=== RGV correction sweep — paired Wilcoxon (EXPLORATORY) ===")
    println("    Applies the pre-registration's statistical RULE to a comparison it")
    println("    does not describe (H1 is Viterbi DP vs greedy). Not confirmatory.")
    df = load_sweep_csvs(csv_paths)
    analysis = run_paired_analysis(df;
        treatment = treatment, control = control, metrics = metrics,
        metric_source = _pw_arg("--metric-source"),
        allow_mixed_src = "--allow-mixed-src" in ARGS,
        keep_not_ok = "--keep-not-ok" in ARGS)

    mkpath(output_dir)
    report = write_paired_report(
        joinpath(output_dir, "report.md"), analysis; csv_paths = csv_paths)
    json_path = joinpath(output_dir, "results.json")
    open(json_path, "w") do io
        JSON.print(io, paired_analysis_json(analysis; csv_paths = csv_paths), 2)
    end

    println("Seeds present : $(join(string.(analysis.seeds), ", "))")
    println("Rows analyzed : $(analysis.n_rows_analyzed) / $(analysis.n_input_rows)")
    for (i, r) in enumerate(analysis.results)
        println("  $(r.metric): n=$(r.n_pairs) pairs, p=$(_pw_fmt(r.test.pvalue)) " *
                "(FDR $(_pw_fmt(analysis.adjusted_p[i]))) -> $(analysis.verdicts[i])")
    end
    println("\nWrote:\n  $report\n  $json_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
