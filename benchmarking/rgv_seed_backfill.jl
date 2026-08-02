# Backfill the `seed` column into pre-existing RGV sweep CSVs (bead td-59o7).
# ===========================================================================
#
# `rhizomorph_correction_validation_sweep.jl` now writes `seed` as a first-class
# column. CSVs produced BEFORE that change do not carry it, so their replicate
# rows are indistinguishable once merged and the pre-registration's
# paired-Wilcoxon rule over seeds 42/123/456 cannot be run against them.
#
# This tool adds the column to an existing CSV, and is deliberately strict about
# provenance, because a seed value invented for a historical run is worse than a
# missing one: it makes an unpairable run look pairable.
#
#   * The seed must come from an explicit `--seed N`, or be RECOVERED from a run
#     log via `--run-log`.
#   * `--assume-default-seed` is the only way to fall back on the harness's
#     documented default (42), and it records that weaker provenance in the
#     sidecar rather than pretending the value was observed.
#   * The original CSV is NEVER modified in place. Output goes to a new path and a
#     sidecar JSON records the source file, the seed, and exactly how the seed was
#     determined, so the audit trail survives.
#
# WHAT IS RECOVERABLE
# -------------------
# `run_rgv_sweep_{lrc,nersc}.sbatch` export `MYCELIA_RGV_SEED` and now echo
# `SEED=<n>` into the run log; the harness banner now prints `Seed : <n>`. Any of
# those three forms is recognized. The completed Lovelace sweep predates all
# three, so its seed is recoverable only as the documented default 42 — which is
# exactly the case `--assume-default-seed` exists to label honestly.
#
# USAGE
# -----
#   # seed stated explicitly
#   julia --project=. benchmarking/rgv_seed_backfill.jl \
#       --csv benchmarking/results/rhizomorph_correction_validation_sweep_20260711_125013.csv \
#       --seed 42
#
#   # seed recovered from a run log
#   julia --project=. benchmarking/rgv_seed_backfill.jl \
#       --csv results/....csv --run-log ~/workspace/slurmlogs/lovelace-jobs/rgv_sweep_*.log
#
#   # no log available; accept the documented default, labelled as such
#   julia --project=. benchmarking/rgv_seed_backfill.jl --csv results/....csv \
#       --assume-default-seed
#
#   Flags:
#     --csv PATH                 (required) sweep CSV to backfill
#     --seed N                   seed to write; provenance "explicit"
#     --run-log PATH             recover the seed from a run log
#     --assume-default-seed      fall back to 42; provenance "documented-default"
#     --metric-source VALUE      also backfill `metric_source` (e.g. "quast")
#     --quast-min-contig BP      also backfill `quast_min_contig` (e.g. 500)
#     --output PATH              output CSV (default <csv-stem>_seedbackfill.csv)
#     --force                    overwrite an existing output file
#
# WHY THE DEFINITION COLUMNS ARE HERE TOO. A `seed`-only backfill DEAD-ENDS: the
# analysis accepts the seed and then refuses the same file for carrying no
# metric-definition column, and nothing else could supply one. Both are opt-in and
# operator-supplied — never inferred.
#
# THE ASSERTION TRAVELS IN THE TABLE, NOT ONLY IN THE SIDECAR. `--metric-source`
# states, on operator say-so, WHICH DEFINITION PRODUCED EVERY METRIC IN THE FILE.
# That is a much stronger claim than a seed label, and unlike `--seed` there is no
# `--run-log` that could recover it: nothing verifies it. Recording it only in
# `<out>.seed_backfill.json` is not enough, because NO consumer reads that sidecar
# — `rgv_paired_wilcoxon.jl` calls `load_sweep_csvs(csv_paths)` and nothing else.
# A genuinely mixed legacy table (QUAST-scored rows plus `internal:quast-failed`
# rows, no `metric_source` column) backfilled with `--metric-source quast` would
# therefore produce a file the guard certifies as uniform, and the analysis report
# would print "`metric_source_guard.jl` enforced it — no override was used" over
# numbers mixing an alignment-validated NGA50 with a size-ratio proxy. That is
# strictly worse than the dead end it replaced, because the dead end failed closed.
#
# So every row this tool SUPPLIES a definition for is labelled in a sibling
# `<column>_provenance` column carrying `DEFINITION_PROVENANCE_ASSERTED`. Rows that
# already carried the value keep `"observed"`. The column travels with the CSV into
# every downstream consumer, and `rgv_paired_wilcoxon.jl` reads it: an analysis over
# operator-asserted rows cannot print the unqualified guard assurance, and
# `results.json` carries the fact in machine-readable form.

import CSV
import DataFrames
import JSON
import Dates

"""
The harness default for `MYCELIA_RGV_SEED`, mirrored from
`rhizomorph_correction_validation_sweep.jl`. Only used when
`--assume-default-seed` is passed, and always recorded as the weaker
`"documented-default"` provenance.
"""
const RGV_DEFAULT_SEED = 42

# Three forms in which the seed can appear in a run log:
#   sbatch echo    "    ERR=...  K=21  SEED=42  QUAST=true"
#   sbatch export  'export MYCELIA_RGV_SEED="42"'  /  MYCELIA_RGV_SEED=42
#   harness banner "Seed           : 42"
const SEED_LOG_PATTERNS = (
    r"\bSEED=\"?([\d,]+)\"?"i,
    r"\bMYCELIA_RGV_SEED=\"?([\d,]+)\"?",
    r"^\s*Seed\s*:\s*(\d+)\s*$"
)

"""
    recover_seed_from_log(log_path) -> Union{Int, Nothing}

Scan a run log for the seed. Returns the seed, or `nothing` if no recognized form
appears.

Throws an `ErrorException` if the log contains MORE THAN ONE distinct seed value:
that means the log covers several runs, so attributing one seed to the CSV would
be a guess.
"""
function recover_seed_from_log(log_path::AbstractString)
    isfile(log_path) || error("run log not found: $log_path")
    found = Set{Int}()
    for line in eachline(log_path)
        for pat in SEED_LOG_PATTERNS
            m = match(pat, line)
            m === nothing && continue
            # The capture may be a COMMA LIST (`SEED=42,123,456`) since
            # MYCELIA_RGV_SEED became a list. Taking `parse(Int, ...)` of the whole
            # match would throw; taking only the first number would silently
            # attribute one seed to a multi-seed run.
            for piece in split(m.captures[1], ",")
                isempty(strip(piece)) && continue
                push!(found, parse(Int, strip(piece)))
            end
        end
    end
    isempty(found) && return nothing
    if length(found) > 1
        error("this run used MORE THAN ONE seed ($(join(sort(collect(found)), ", ")) " *
              "in $log_path), so a single value cannot be attributed to every row of " *
              "the CSV. Since MYCELIA_RGV_SEED became a list, this is the normal shape " *
              "of a replicate run: such a CSV must carry its own per-row `seed` column " *
              "from the harness. Backfilling one seed here would label every replicate " *
              "identically and make an unpairable table look pairable — the exact " *
              "failure this tool exists to prevent (see the pseudo-replication guard " *
              "in rgv_paired_wilcoxon.jl).")
    end
    return first(found)
end

"""
    insert_seed_column!(df, seed) -> DataFrames.DataFrame

Insert `seed` immediately after `k`, matching the live sweep schema's column
order. If `k` is absent the column is appended.

Throws if `df` already has a `seed` column holding a value other than `seed` —
silently rewriting an existing replicate identifier would corrupt provenance
rather than add it.
"""
function insert_seed_column!(df::DataFrames.DataFrame, seed::Integer)
    if "seed" in DataFrames.names(df)
        # A legacy CSV can carry the HEADER with empty cells, which CSV.jl reads as
        # `missing`. That is "present but unpopulated", not "present with a
        # conflicting value" — and it is exactly the historical shape this tool
        # targets, so treating it as a conflict made the documented workflow
        # impossible on its intended input. Fill it instead.
        existing = unique(skipmissing(df.seed))
        if isempty(existing)
            df[!, :seed] = fill(Int(seed), DataFrames.nrow(df))
            return df
        end
        if length(existing) == 1 && first(existing) == seed
            # PARTIALLY populated is the case between "all missing" and "conflicting":
            # some rows carry `seed`, others are `missing`. `skipmissing` above hides
            # those holes, so returning early here reported success while leaving the
            # table unpairable — and `require_pairing_schema` then rejects it with a
            # hint pointing back at THIS tool, so the operator ping-pongs forever.
            # The populated values already agree with `seed`, so filling the holes
            # asserts nothing new.
            if any(ismissing, df.seed)
                df[!, :seed] = fill(Int(seed), DataFrames.nrow(df))
            end
            return df
        end
        error("CSV already has a `seed` column with value(s) " *
              "$(join(string.(existing), ", ")); refusing to overwrite it with $seed.")
    end
    n = DataFrames.nrow(df)
    pos = findfirst(==("k"), DataFrames.names(df))
    idx = pos === nothing ? DataFrames.ncol(df) + 1 : pos + 1
    DataFrames.insertcols!(df, idx, :seed => fill(Int(seed), n))
    return df
end

"""
Value written into `<definition-column>_provenance` for every row whose metric
definition was SUPPLIED BY THIS TOOL rather than observed by the run that produced
the metrics.

`rgv_paired_wilcoxon.jl` matches on the `"operator-asserted"` substring, so the
prefix is part of the cross-file contract and is pinned by a test in
`test/4_assembly/rgv_paired_wilcoxon_test.jl`.
"""
const DEFINITION_PROVENANCE_ASSERTED = "operator-asserted (backfill)"

"""
Value kept for rows that already carried the definition value before the backfill
ran. Those rows are evidence; the asserted ones are a claim.
"""
const DEFINITION_PROVENANCE_OBSERVED = "observed"

"""
    definition_provenance_column(name) -> Symbol

Sibling provenance column for a metric-definition column, e.g.
`:metric_source` -> `:metric_source_provenance`.
"""
definition_provenance_column(name::Symbol) = Symbol(String(name) * "_provenance")

"""
    _mark_definition_provenance!(df, name, asserted) -> DataFrames.DataFrame

Record, per row, whether this tool SUPPLIED the value in definition column `name`.

`asserted[i]` is true for rows written by the backfill. Rows it did not write keep
whatever provenance they already carried (or `"observed"` if the column is new), so
re-running the tool can only ever ADD assertions — it can never launder a
previously asserted row back into an observed one.

No column is created when nothing was asserted and none exists: a table this tool
did not change should not gain a column claiming it did.
"""
function _mark_definition_provenance!(df::DataFrames.DataFrame, name::Symbol,
        asserted::AbstractVector{Bool})
    prov = definition_provenance_column(name)
    n = DataFrames.nrow(df)
    present = String(prov) in DataFrames.names(df)
    (!any(asserted) && !present) && return df
    prior = present ?
            String[ismissing(v) ? DEFINITION_PROVENANCE_OBSERVED : string(v)
                   for v in getproperty(df, prov)] :
            fill(DEFINITION_PROVENANCE_OBSERVED, n)
    labels = String[asserted[i] ? DEFINITION_PROVENANCE_ASSERTED : prior[i]
                    for i in 1:n]
    if present
        df[!, prov] = labels
    else
        DataFrames.insertcols!(df, DataFrames.ncol(df) + 1, prov => labels)
    end
    return df
end

"""
    insert_definition_column!(df, name, value) -> DataFrames.DataFrame

Append a metric-definition column (`:metric_source` / `:quast_min_contig`) with a
single operator-supplied `value`, and label every row it supplies in the sibling
`<name>_provenance` column.

Refuses to overwrite an existing column holding a different value, for the same
reason `insert_seed_column!` does: asserting a definition a run did not have makes
an unusable table look usable, and the metric-definition guard would then pass on
a claim nobody verified.

The provenance column exists because refusing to overwrite is not sufficient. When
the column is ABSENT there is nothing to conflict with, so the value is accepted
unconditionally — and a legacy table that really did mix definitions has exactly
that shape. The guard downstream can then only confirm that the asserted labels
agree with each other, which is not evidence about how the runs were scored. The
provenance column is what carries that distinction to the consumer.
"""
function insert_definition_column!(df::DataFrames.DataFrame, name::Symbol, value)
    col = String(name)
    if col in DataFrames.names(df)
        existing = unique(skipmissing(getproperty(df, name)))
        if isempty(existing)          # present but unpopulated — fill, do not refuse
            df[!, name] = fill(value, DataFrames.nrow(df))
            return _mark_definition_provenance!(df, name,
                trues(DataFrames.nrow(df)))
        end
        if length(existing) == 1 && first(existing) == value
            # PARTIALLY populated, exactly as in `insert_seed_column!`: `skipmissing`
            # above hides the holes, so returning early reported success while
            # leaving rows whose definition is still `missing`. The metric-definition
            # guard renders `missing` as the literal `"missing"` and counts it as a
            # DISTINCT value, so those holes make the guard raise later with a
            # message that points nowhere. The populated values already agree with
            # `value`, so filling the holes asserts nothing new about them.
            holes = BitVector(ismissing(v) for v in getproperty(df, name))
            if any(holes)
                df[!, name] = fill(value, DataFrames.nrow(df))
            end
            # Only the FILLED rows are asserted; the rows that already carried the
            # value are evidence and stay labelled as observed.
            return _mark_definition_provenance!(df, name, holes)
        end
        error("CSV already has a `$col` column with value(s) " *
              "$(join(string.(existing), ", ")); refusing to overwrite it with " *
              "$value.")
    end
    DataFrames.insertcols!(df, DataFrames.ncol(df) + 1,
        name => fill(value, DataFrames.nrow(df)))
    return _mark_definition_provenance!(df, name, trues(DataFrames.nrow(df)))
end

"""
    backfill_seed(csv_path; seed, provenance, output_path=nothing, force=false)
        -> (output_path, sidecar_path, nrows)

Write a copy of `csv_path` with a `seed` column, plus a sidecar
`<output>.seed_backfill.json` recording the source CSV, the seed, the provenance
label, and the timestamp. The input file is never modified.
"""
function backfill_seed(csv_path::AbstractString;
        seed::Integer,
        provenance::AbstractString,
        metric_source::Union{Nothing, AbstractString} = nothing,
        quast_min_contig::Union{Nothing, Integer} = nothing,
        output_path::Union{Nothing, AbstractString} = nothing,
        force::Bool = false)
    isfile(csv_path) || error("CSV not found: $csv_path")
    out = output_path === nothing ?
          replace(csv_path, r"\.csv$" => "") * "_seedbackfill.csv" :
          String(output_path)
    if isfile(out) && !force
        error("output already exists: $out (pass --force to overwrite)")
    end
    df = CSV.read(csv_path, DataFrames.DataFrame)
    insert_seed_column!(df, seed)
    # A `seed`-only backfill DEAD-ENDS: the analysis then refuses the same file for
    # having no metric-definition column, and nothing else could supply one. The
    # documented backfill -> analyse workflow has to be completable, so the same
    # provenance discipline is extended to the definition columns — each is opt-in,
    # explicit, labelled per row IN THE CSV, and recorded in the sidecar, never
    # inferred.
    if metric_source !== nothing
        insert_definition_column!(df, :metric_source, String(metric_source))
    end
    if quast_min_contig !== nothing
        insert_definition_column!(df, :quast_min_contig, Int(quast_min_contig))
    end
    CSV.write(out, df)
    # The sidecar is an audit record, NOT the disclosure mechanism: no consumer
    # reads it. The disclosure that has to reach the analysis rides in the CSV's
    # `<column>_provenance` columns, named here so the two cannot drift apart.
    provenance_columns = String[String(definition_provenance_column(c))
                                for c in (:metric_source, :quast_min_contig)
                                if String(definition_provenance_column(c)) in DataFrames.names(df)]
    sidecar = out * ".seed_backfill.json"
    open(sidecar, "w") do io
        JSON.print(io,
            Dict(
                "source_csv" => abspath(csv_path),
                "output_csv" => abspath(out),
                "seed" => Int(seed),
                "seed_provenance" => String(provenance),
                "metric_source_backfilled" => metric_source === nothing ? nothing :
                                              String(metric_source),
                "quast_min_contig_backfilled" => quast_min_contig === nothing ? nothing :
                                                 Int(quast_min_contig),
                "definition_provenance_columns" => provenance_columns,
                "definition_provenance_asserted_label" => DEFINITION_PROVENANCE_ASSERTED,
                "rows" => DataFrames.nrow(df),
                "backfilled_at" => string(Dates.now()),
                "tool" => "benchmarking/rgv_seed_backfill.jl",
                "bead" => "td-59o7"
            ), 2)
    end
    # `provenance_columns` is returned, not just written to the sidecar, so the CLI
    # can report what it ACTUALLY wrote. It previously announced "rows labelled
    # <asserted> in metric_source_provenance" unconditionally — including on the one
    # shape where `_mark_definition_provenance!` writes nothing at all (a definition
    # column that already carries the requested value, so no row is asserted and no
    # provenance column exists to update). Trailing positional, so existing
    # three-way destructuring keeps working.
    return out, sidecar, DataFrames.nrow(df), provenance_columns
end

# === CLI ====================================================================

function _bf_arg(flag)
    i = findfirst(==(flag), ARGS)
    return (i !== nothing && i < length(ARGS)) ? ARGS[i + 1] : nothing
end

function main()
    csv_path = _bf_arg("--csv")
    if csv_path === nothing
        println(stderr,
            "ERROR: --csv PATH is required.\n" *
            "  julia --project=. benchmarking/rgv_seed_backfill.jl --csv <sweep.csv> " *
            "(--seed N | --run-log PATH | --assume-default-seed)")
        return 2
    end

    explicit = _bf_arg("--seed")
    run_log = _bf_arg("--run-log")
    assume_default = "--assume-default-seed" in ARGS

    seed = nothing
    provenance = ""
    if explicit !== nothing
        seed = parse(Int, explicit)
        provenance = "explicit (--seed $seed)"
    elseif run_log !== nothing
        recovered = recover_seed_from_log(run_log)
        if recovered === nothing
            if assume_default
                seed = RGV_DEFAULT_SEED
                provenance = "documented-default ($RGV_DEFAULT_SEED); no seed found in " *
                             "$(abspath(run_log))"
            else
                println(stderr,
                    "ERROR: no seed found in $run_log. That log predates SEED being " *
                    "echoed. Pass --seed N if you know it, or --assume-default-seed to " *
                    "record the documented default ($RGV_DEFAULT_SEED) with weaker provenance.")
                return 1
            end
        else
            seed = recovered
            provenance = "recovered-from-run-log ($(abspath(run_log)))"
        end
    elseif assume_default
        seed = RGV_DEFAULT_SEED
        provenance = "documented-default ($RGV_DEFAULT_SEED); no run log supplied"
    else
        println(stderr,
            "ERROR: the seed must be justified. Pass one of --seed N, --run-log PATH, " *
            "or --assume-default-seed. Inventing a seed would make an unpairable run " *
            "look pairable.")
        return 2
    end

    ms = _bf_arg("--metric-source")
    qmc = _bf_arg("--quast-min-contig")
    out, sidecar,
    n,
    prov_cols = backfill_seed(csv_path;
        seed = seed, provenance = provenance,
        metric_source = ms,
        quast_min_contig = qmc === nothing ? nothing : parse(Int, qmc),
        output_path = _bf_arg("--output"), force = "--force" in ARGS)
    println("=== RGV seed backfill ===")
    println("Source     : $csv_path")
    println("Seed       : $seed")
    println("Provenance : $provenance")
    # Report what was WRITTEN, not what was requested. Naming a definition the CSV
    # already carries asserts nothing, so no provenance column is created — and
    # announcing one anyway is an assurance the sidecar contradicts in the same run.
    _bf_wrote(col) = String(definition_provenance_column(col)) in prov_cols
    ms === nothing || println("metric_source    : $ms (operator-supplied)" *
            (_bf_wrote(:metric_source) ?
             "; rows supplied are labelled " *
             "\"$DEFINITION_PROVENANCE_ASSERTED\" in " *
             "metric_source_provenance" :
             "; already present with this value — nothing asserted, " *
             "so NO metric_source_provenance column was written"))
    qmc === nothing || println("quast_min_contig : $qmc (operator-supplied)" *
            (_bf_wrote(:quast_min_contig) ?
             "; rows supplied are labelled " *
             "\"$DEFINITION_PROVENANCE_ASSERTED\" in " *
             "quast_min_contig_provenance" :
             "; already present with this value — nothing asserted, " *
             "so NO quast_min_contig_provenance column was written"))
    if !isempty(prov_cols)
        println("             NOTE: an operator-asserted definition is a CLAIM, not an " *
                "observation. rgv_paired_wilcoxon.jl reads the provenance column and " *
                "will refuse to report the guard as enforced over these rows.")
    end
    println("Rows       : $n")
    println("Wrote      : $out")
    println("Sidecar    : $sidecar")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
