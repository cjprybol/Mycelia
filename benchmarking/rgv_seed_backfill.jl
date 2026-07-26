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
#     --output PATH              output CSV (default <csv-stem>_seedbackfill.csv)
#     --force                    overwrite an existing output file

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
    r"\bSEED=\"?(\d+)\"?"i,
    r"\bMYCELIA_RGV_SEED=\"?(\d+)\"?",
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
            push!(found, parse(Int, m.captures[1]))
        end
    end
    isempty(found) && return nothing
    if length(found) > 1
        error("ambiguous seed in $log_path: found $(join(sort(collect(found)), ", ")). " *
              "The log appears to cover more than one run — split it or pass --seed explicitly.")
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
        existing = unique(df.seed)
        if length(existing) == 1 && first(existing) == seed
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
    backfill_seed(csv_path; seed, provenance, output_path=nothing, force=false)
        -> (output_path, sidecar_path, nrows)

Write a copy of `csv_path` with a `seed` column, plus a sidecar
`<output>.seed_backfill.json` recording the source CSV, the seed, the provenance
label, and the timestamp. The input file is never modified.
"""
function backfill_seed(csv_path::AbstractString;
        seed::Integer,
        provenance::AbstractString,
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
    CSV.write(out, df)
    sidecar = out * ".seed_backfill.json"
    open(sidecar, "w") do io
        JSON.print(io,
            Dict(
                "source_csv" => abspath(csv_path),
                "output_csv" => abspath(out),
                "seed" => Int(seed),
                "seed_provenance" => String(provenance),
                "rows" => DataFrames.nrow(df),
                "backfilled_at" => string(Dates.now()),
                "tool" => "benchmarking/rgv_seed_backfill.jl",
                "bead" => "td-59o7"
            ), 2)
    end
    return out, sidecar, DataFrames.nrow(df)
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

    out, sidecar,
    n = backfill_seed(csv_path;
        seed = seed, provenance = provenance,
        output_path = _bf_arg("--output"), force = "--force" in ARGS)
    println("=== RGV seed backfill ===")
    println("Source     : $csv_path")
    println("Seed       : $seed")
    println("Provenance : $provenance")
    println("Rows       : $n")
    println("Wrote      : $out")
    println("Sidecar    : $sidecar")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
