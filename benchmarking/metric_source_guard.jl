# Metric-definition consistency guard for the benchmarking analysis path.
# =======================================================================
#
# Pure and dependency-free (like `rhizomorph_scale_guard.jl`,
# `quast_report_parsing.jl` and `quast_min_contig.jl`): no Mycelia, no
# DataFrames import, no network. It reaches columns through `hasproperty` /
# `getproperty`, so it works on a `DataFrames.DataFrame`, a `NamedTuple` of
# vectors, or anything else column-shaped — and its unit test runs in
# milliseconds.
#
# WHY THIS FILE EXISTS (bead td-9p91)
# -----------------------------------
# When QUAST fails, the harnesses degrade gracefully to internal size-ratio
# metrics and record `metric_source="internal:quast-failed"` while still
# reporting `ok=true`. Keeping the job alive is correct. The hazard is downstream:
# a single results table can then hold rows whose metrics come from TWO DIFFERENT
# DEFINITIONS —
#
#   * `metric_source="quast"`   — alignment-validated (QUAST aligns contigs to the
#                                 reference; NGA50 and misassemblies are real)
#   * `metric_source="internal*"` — a size ratio (`total_length / genome_len`) and a
#                                 length-sorted N50, BLIND to misassembly and to
#                                 mis-alignment
#
# A naive `groupby` across those rows compares incompatible quantities, and the
# resulting difference reads as a real effect.
#
# WHAT WAS ACTUALLY OBSERVED. An earlier version of this comment said the mix was
# "observed concretely on Lawrencium job 24247925: T4 rows MAY carry internal
# metrics while Lambda rows carry QUAST metrics" — but observed and "may" cannot
# both be true, and the T4 and Lambda rows came from different jobs and were never
# in one table. The real observation is stronger and needs no cross-job inference:
# the completed Lambda sweep produced a SINGLE 12-row CSV holding 4 rows with
# `metric_source="quast"` and 8 with `"internal:quast-failed"` (bead td-4e19d.1).
# And because `run_quast` is invoked PER ARM, the naive and iterative arms of one
# cell can carry different definitions — a mixed PAIR, which defeats the paired
# design directly rather than merely biasing a group mean.
#
# This guard is the DURABLE fix; `quast_min_contig.jl` (td-28o0) is the specific
# one. Fixing the `--min-contig` floor removes one trigger, but graceful
# degradation will always be able to produce mixed `metric_source`, so any
# aggregation must assert consistency rather than assume it.
#
# WHAT COUNTS AS "THE METRIC DEFINITION"
# --------------------------------------
# `metric_source` is the obvious axis, but it is not the only one. The threshold
# QUAST was run at is equally part of the definition: two rows both labelled
# `metric_source="quast"` but scored at `--min-contig` 4,850 vs 5,000 are also not
# comparable. `METRIC_DEFINITION_COLUMNS` therefore covers both, which is why the
# RGV sweep records `quast_min_contig` per row — it makes a future change to the
# threshold policy trip THIS guard instead of passing silently.
#
# FAIL-CLOSED ON A MISSING COLUMN — PER AXIS, NOT PER SET
# -------------------------------------------------------
# If a table carries none of the definition columns, the guard RAISES: a table
# with no provenance column is not a table that has been shown to be consistent,
# it is a table that cannot be checked.
#
# The same reasoning applies to each axis INDIVIDUALLY, and an earlier version got
# this wrong. It raised only when EVERY column was absent, so a table with a
# uniform `metric_source` and no `quast_min_contig` — the shape of every CSV
# produced before this policy existed — passed and reported a single consistent
# definition over an axis it had never examined. A table concatenating runs scored
# at different thresholds has exactly that shape. Any requested axis that is absent
# while others are present therefore also raises; pass `require_all_columns=false`
# to accept the narrower check (a warning then names which axes went unchecked), or
# pass an explicit `columns=` to record the narrower claim deliberately.

"""
Columns that together define "which metric definition produced this row".

  * `:metric_source`     — `"quast"` (alignment-validated) vs `"internal*"` (size-ratio proxy)
  * `:quast_min_contig`  — the `--min-contig` threshold QUAST was run at

Both must be constant across an aggregation for the aggregate to mean anything.
"""
const METRIC_DEFINITION_COLUMNS = (:metric_source, :quast_min_contig)

"""
    MixedMetricDefinitionError

Raised when an aggregation would span more than one metric definition. Carries
the offending column(s), the distinct values observed, and (for the grouped
check) the group keys at fault, so the message is actionable without re-running
the analysis.
"""
struct MixedMetricDefinitionError <: Exception
    context::String
    offenders::Vector{Tuple{Union{Nothing, String}, Symbol, Vector{String}}}
end

_definition_label(v) = v === missing ? "missing" : string(v)

function _render_offenders(offenders)
    lines = String[]
    for (group, col, vals) in offenders
        prefix = group === nothing ? "  " : "  [group $group] "
        push!(lines, prefix * "$col: " * join(("\"$v\"" for v in vals), ", "))
    end
    return join(lines, "\n")
end

function Base.showerror(io::IO, e::MixedMetricDefinitionError)
    print(io,
        """
        MixedMetricDefinitionError: refusing to aggregate in $(e.context) — the rows \
        were scored under more than one metric definition.

        $(_render_offenders(e.offenders))

        Rows labelled metric_source="quast" carry ALIGNMENT-VALIDATED metrics; rows \
        labelled "internal*" carry size-ratio proxies that are blind to misassembly. \
        Rows scored at different quast_min_contig thresholds were filtered differently. \
        Aggregating across either boundary compares incompatible quantities, and the \
        difference reads as a real effect.

        Fix the upstream cause (see beads td-28o0 / td-9p91), or filter to a single \
        definition, or pass allow_mixed_src=true to override deliberately (a loud \
        warning is emitted and the mixed values are named).""")
    return nothing
end

"""
    metric_definition_summary(df; columns=METRIC_DEFINITION_COLUMNS)
        -> Vector{Pair{Symbol, Vector{String}}}

Distinct values of each present definition column, sorted, with `missing`
rendered as the literal `"missing"` so an absent provenance value participates in
the consistency check rather than being skipped.

Columns in `columns` that the table does not carry are omitted from the summary
(use `assert_single_metric_definition` to fail closed when none are present).
"""
function metric_definition_summary(df; columns = METRIC_DEFINITION_COLUMNS)
    summary = Pair{Symbol, Vector{String}}[]
    for col in columns
        hasproperty(df, col) || continue
        vals = sort(unique(String[_definition_label(v) for v in getproperty(df, col)]))
        push!(summary, col => vals)
    end
    return summary
end

function _require_definition_columns(df, columns, context; require_all::Bool = true)
    present = Symbol[c for c in columns if hasproperty(df, c)]
    absent = Symbol[c for c in columns if !hasproperty(df, c)]
    if isempty(present)
        throw(ArgumentError(
            "no metric-definition column present in $context (looked for " *
            join((":$c" for c in columns), ", ") * "). A table with no provenance " *
            "column cannot be shown to be consistent — add one, or pass an explicit " *
            "`columns=` naming this table's provenance column."))
    end
    # PARTIAL absence is the subtler hole, and it was a real one. The doctrine in
    # this file's header is per-AXIS ("a table that cannot be checked"), but the
    # implementation was per-SET: it raised only when EVERY column was missing, and
    # `metric_definition_summary` silently skipped any individually-absent column.
    # So a table with a uniform `metric_source` and no `quast_min_contig` — the
    # shape of every CSV produced before this PR — passed and reported a single
    # consistent definition, over an axis that was never examined. A table
    # concatenating runs scored at different thresholds is exactly that shape.
    if require_all && !isempty(absent)
        throw(ArgumentError(
            "metric-definition axis not checkable in $context: " *
            join(("`$c`" for c in absent), ", ") *
            " absent while " * join(("`$c`" for c in present), ", ") *
            " present. Reporting 'consistent' would assert something about an axis " *
            "this table cannot answer for. Either add the column, or pass an " *
            "explicit `columns=` limited to the axes this table actually carries " *
            "(which records the narrower claim), or `require_all_columns=false` to " *
            "accept the narrower check with a warning."))
    end
    if !isempty(absent)
        @warn "metric-definition axis NOT checked (column absent) — the " *
              "consistency claim below covers only the axes listed as checked." context=String(context) absent=join(("$c" for c in absent), ", ") checked=join(("$c" for c in present), ", ")
    end
    return present
end

"""
    assert_single_metric_definition(df; columns, context, allow_mixed_src=false)
        -> Vector{Pair{Symbol, Vector{String}}}

Assert that every row of `df` was scored under ONE metric definition, i.e. that
each column in `columns` holds a single distinct value. Returns the definition
summary on success.

Throws `MixedMetricDefinitionError` when a definition column holds more than one
value, and `ArgumentError` when the table carries none of `columns` (fail-closed:
unverifiable is not the same as verified).

Set `allow_mixed_src=true` to proceed anyway; the mix is then reported through a
loud `@warn` that names every mixed value, so the override is visible in the run
log rather than implicit.

# Example

```julia
# One definition: passes, returns the summary.
assert_single_metric_definition(quast_only_rows; context = "H1 NGA50 aggregate")

# Mixed: raises, naming "internal:quast-failed" and "quast".
assert_single_metric_definition(all_rows; context = "H1 NGA50 aggregate")
```
"""
function assert_single_metric_definition(df;
        columns = METRIC_DEFINITION_COLUMNS,
        context::AbstractString = "aggregation",
        allow_mixed_src::Bool = false,
        require_all_columns::Bool = true)
    present = _require_definition_columns(df, columns, context;
        require_all = require_all_columns)
    summary = metric_definition_summary(df; columns = present)
    offenders = Tuple{Union{Nothing, String}, Symbol, Vector{String}}[(nothing, col, vals)
                                                                      for (col, vals) in
                                                                          summary
                                                                      if length(vals) > 1]
    isempty(offenders) && return summary
    if allow_mixed_src
        @warn "MIXED METRIC DEFINITION accepted via allow_mixed_src=true — the "*
              "aggregate below compares rows scored under different definitions and "*
              "is NOT validation-grade." context=String(context) mixed=_render_offenders(offenders)
        return summary
    end
    throw(MixedMetricDefinitionError(String(context), offenders))
end

"""
    assert_single_metric_definition_per_group(df, groupcols; columns, context,
                                             allow_mixed_src=false) -> Int

Grouped form of [`assert_single_metric_definition`](@ref): assert that EACH group
defined by `groupcols` is internally single-definition. This is the check a
`groupby`-then-aggregate needs, since each group becomes one aggregate.

Every offending group is collected before raising, so one run reports the full
extent of the problem instead of the first instance. Returns the number of groups
checked.
"""
function assert_single_metric_definition_per_group(df, groupcols;
        columns = METRIC_DEFINITION_COLUMNS,
        context::AbstractString = "grouped aggregation",
        allow_mixed_src::Bool = false,
        require_all_columns::Bool = true)
    present = _require_definition_columns(df, columns, context;
        require_all = require_all_columns)
    for col in groupcols
        hasproperty(df, col) ||
            throw(ArgumentError("group column :$col not present in $context"))
    end
    keyvecs = [getproperty(df, col) for col in groupcols]
    defvecs = [getproperty(df, col) for col in present]
    nrows = length(first(defvecs))

    group_order = String[]
    group_values = Dict{String, Vector{Set{String}}}()
    for i in 1:nrows
        key = join((_definition_label(v[i]) for v in keyvecs), " | ")
        if !haskey(group_values, key)
            group_values[key] = [Set{String}() for _ in present]
            push!(group_order, key)
        end
        for (j, v) in enumerate(defvecs)
            push!(group_values[key][j], _definition_label(v[i]))
        end
    end

    offenders = Tuple{Union{Nothing, String}, Symbol, Vector{String}}[]
    for key in group_order
        for (j, col) in enumerate(present)
            vals = group_values[key][j]
            length(vals) > 1 && push!(offenders, (key, col, sort(collect(vals))))
        end
    end
    isempty(offenders) && return length(group_order)
    if allow_mixed_src
        @warn "MIXED METRIC DEFINITION accepted via allow_mixed_src=true in "*
              "$(length(offenders)) group(s) — those aggregates compare rows scored "*
              "under different definitions and are NOT validation-grade." context=String(context) mixed=_render_offenders(offenders)
        return length(group_order)
    end
    throw(MixedMetricDefinitionError(String(context), offenders))
end

"""
    guarded_aggregate(f, df; columns, context, allow_mixed_src=false)

Run `f(df)` only after [`assert_single_metric_definition`](@ref) passes. The
one-line way to make an analysis step definition-safe:

```julia
median_nga50 = guarded_aggregate(rows; context = "median NGA50") do d
    Statistics.median(skipmissing(d.quast_nga50))
end
```
"""
function guarded_aggregate(f, df;
        columns = METRIC_DEFINITION_COLUMNS,
        context::AbstractString = "aggregation",
        allow_mixed_src::Bool = false,
        require_all_columns::Bool = true)
    assert_single_metric_definition(df; columns, context, allow_mixed_src,
        require_all_columns)
    return f(df)
end
