# Track A cross-host checkpoint merge (bead td-bblmi, epic td-yddpr).
# ===================================================================
#
# The Track A baseline 288-cell matrix is computed on MORE THAN ONE HOST at once,
# and its checkpoints are per-host and per-cell:
#
#   <results_root>/cells/<cell_id>/cell_result.json
#
# Concretely (2026-07-26):
#   * Lovelace   — Lambda complete (72/72), T4 in progress; phi29 / SARS-CoV-2 untouched
#   * Lawrencium — the `--organisms phi29,SARS-CoV-2` shard (144 cells)
#
# Those trees live on separate filesystems, so the matrix cannot be analyzed as
# one until they are merged. This tool does the merge, and it is deliberately
# PARANOID about the two ways a merge can silently produce a wrong 288-cell table.
#
# 1. COLLISION DETECTION, NOT LAST-WRITE-WINS
#    The shards were supposed to be disjoint. If the same `cell_id` exists on two
#    hosts with DIFFERENT content, that is a real anomaly — overlapping shards, a
#    stale checkpoint, or two different Mycelia SHAs — and picking a winner would
#    bury it. Byte-equivalent duplicates are benign and reported separately from
#    content collisions, which fail the merge unless `--allow-collisions`.
#
# 2. PER-CELL PROVENANCE, BECAUSE THE HOSTS DO NOT MEASURE THE SAME THING
#    QUAST works on Lawrencium but was failing on Lovelace (conda
#    `ProcessExited(4)`), where the harness degrades to internal metrics. So the
#    merged matrix MIXES METRIC SOURCES BY HOST. The merge therefore records which
#    host produced each cell and classifies whether QUAST actually scored it (see
#    `quast_evidence` below), and reports the breakdown by host — which is the
#    input the `metric_source_guard.jl` check (bead td-9p91) needs downstream.
#
# WHY `quast_evidence` IS AN INFERENCE, AND WHY IT IS A SOUND ONE
# --------------------------------------------------------------
# `track_a_baseline_benchmark.jl`'s row schema has no `metric_source` column: on
# QUAST failure it substitutes `empty_metrics()` (all five metrics zero) and still
# writes `status="ok"`. Adding a `metric_source` column to that schema would
# invalidate every existing checkpoint (`canonical()` would `KeyError`) and break
# resume for the runs currently in flight, so this tool does NOT touch the
# harness schema. It infers instead, from one sound signal:
#
#   n_contigs > 0 AND largest_contig == 0  =>  QUAST did not score this cell
#
# A contig has a length, so a QUAST run that scored a non-empty assembly cannot
# report `Largest contig = 0`. Only the `empty_metrics()` fallback can. This is an
# implication, not a heuristic — unlike "all metrics are zero", which a genuinely
# unalignable assembly could also produce (`genome_fraction = 0` with a real
# `largest_contig`). A follow-up should add a real `metric_source` column to the
# harness once the in-flight runs have finished.
#
# USAGE
# -----
#   julia --project=. benchmarking/track_a_merge_hosts.jl \
#       --host lovelace=/path/to/lovelace/track_a_baseline \
#       --host lrc=/path/to/lrc/track_a_baseline \
#       --output-dir benchmarking/results/track_a_baseline_merged
#
#   Flags:
#     --host LABEL=PATH        source tree (repeatable). PATH may be the results
#                              root or the `cells/` directory itself.
#     --output-dir PATH        merge destination (default results/track_a_baseline_merged)
#     --no-copy-cells          write tables + report only; do not materialize cells/
#     --allow-collisions       report content collisions but exit 0
#     --allow-incomplete       report missing cells but exit 0
#     --assert-single-metric-definition
#                              additionally require ONE quast_evidence class among
#                              non-empty assemblies (fails a host-mixed matrix)
#     --organisms / --technologies / --coverages / --seeds / --arms
#                              narrow the EXPECTED matrix, mirroring the harness's
#                              shard flags (comma-separated)
#
# IDEMPOTENT AND RE-RUNNABLE: both source runs are checkpointed and resumable, so
# cells appear on a host after a first merge. Every output is rebuilt from the
# sources on each run; a cell is copied only when absent or content-identical, and
# a destination cell that DIFFERS from its source is itself reported as a
# collision (host label `<merged-output>`) rather than overwritten.

import JSON
import SHA
import DataFrames
import CSV
import Dates

# Pure, dependency-free metric-definition guard (bead td-9p91). Included at top
# level, matching how the RGV sweep includes its own pure helpers, so
# `--assert-single-metric-definition` can reuse the one definition of the check
# rather than reimplementing it.
include(joinpath(@__DIR__, "metric_source_guard.jl"))

# === Expected matrix ========================================================
#
# MIRRORS `track_a_baseline_benchmark.jl` (which cannot be `include`d — it runs
# the whole benchmark at include time). `test/4_assembly/track_a_merge_hosts_test.jl`
# asserts these constants and the cell-id format still match the harness source,
# so drift fails a test rather than silently mis-computing "missing".

const TRACK_A_ORGANISMS = ("Lambda", "T4", "phi29", "SARS-CoV-2")
const TRACK_A_TECHNOLOGIES = ("illumina", "pacbio", "ont")
const TRACK_A_COVERAGES = (10, 30, 50, 100)
const TRACK_A_SEEDS = (42, 123, 456)
const TRACK_A_ARMS = ("qualmer", "kmer")

"""
    track_a_cell_id(organism, technology, coverage, seed, arm) -> String

The harness's per-cell directory name. Must stay byte-identical to
`track_a_baseline_benchmark.jl`'s
`cell_id = "\$(org)__\$(tech)__\$(cov)x__seed\$(seed)__\$(arm)"`.
"""
function track_a_cell_id(organism, technology, coverage, seed, arm)
    return "$(organism)__$(technology)__$(coverage)x__seed$(seed)__$(arm)"
end

"""
    expected_cell_ids(; organisms, technologies, coverages, seeds, arms) -> Vector{String}

Every cell id in the (possibly narrowed) expected matrix. The full matrix is
4 organisms x 3 technologies x 4 coverages x 3 seeds x 2 arms = 288 cells.
"""
function expected_cell_ids(;
        organisms = TRACK_A_ORGANISMS,
        technologies = TRACK_A_TECHNOLOGIES,
        coverages = TRACK_A_COVERAGES,
        seeds = TRACK_A_SEEDS,
        arms = TRACK_A_ARMS)
    ids = String[]
    for org in organisms, tech in technologies, cov in coverages, seed in seeds, arm in arms
        push!(ids, track_a_cell_id(org, tech, cov, seed, arm))
    end
    return ids
end

# === Canonicalization + digest ==============================================

# Normalize numbers so a pure representation difference (2 vs 2.0, which JSON
# round-tripping can introduce) never reads as a content collision, while any
# genuine value difference always does.
_canon_value(v::Integer) = string(v)
function _canon_value(v::AbstractFloat)
    if isfinite(v) && isinteger(v) && abs(v) < 2.0^53
        return string(Integer(v))
    end
    return repr(v)
end
_canon_value(v::AbstractString) = String(v)
_canon_value(v::Nothing) = "null"
_canon_value(v::Missing) = "missing"
_canon_value(v) = string(v)

"""
    canonical_cell_string(d) -> String

Stable, key-sorted textual rendering of one parsed `cell_result.json`. Two cells
with the same canonical string are the same result regardless of key order or
integer/float formatting.
"""
function canonical_cell_string(d::AbstractDict)
    ks = sort(String[string(k) for k in keys(d)])
    return join(("$(k)=$(_canon_value(d[k]))" for k in ks), "\n")
end

"""
    cell_digest(d) -> String

SHA-256 hex digest of [`canonical_cell_string`](@ref). Used for collision
detection and recorded per cell in the provenance table.
"""
cell_digest(d::AbstractDict) = bytes2hex(SHA.sha256(canonical_cell_string(d)))

"""
    differing_fields(a, b) -> Vector{String}

Field-level diff of two parsed cells, rendered for the collision report. Keys
present in only one side are reported as `<absent>` on the other.
"""
function differing_fields(a::AbstractDict, b::AbstractDict)
    ks = sort(collect(union(Set(string.(keys(a))), Set(string.(keys(b))))))
    diffs = String[]
    for k in ks
        va = haskey(a, k) ? _canon_value(a[k]) : "<absent>"
        vb = haskey(b, k) ? _canon_value(b[k]) : "<absent>"
        va == vb || push!(diffs, "$k: $va vs $vb")
    end
    return diffs
end

# === Metric provenance inference ============================================

"""
    quast_evidence(cell) -> String

Classify whether QUAST actually scored this cell. See the file header for why the
`largest_contig == 0` signal is an implication rather than a heuristic.

  * `"n/a:cell-error"`          — the harness recorded `status="error"`
  * `"n/a:empty-assembly"`      — no contigs to score
  * `"quast:scored"`            — QUAST reported a real `Largest contig`
  * `"unknown:quast-unscored"`  — contigs exist but `largest_contig == 0`, so the
                                  metrics are the `empty_metrics()` fallback, not
                                  alignment-validated values
"""
function quast_evidence(cell::AbstractDict)
    status = get(cell, "status", "")
    status == "error" && return "n/a:cell-error"
    n_contigs = _as_number(get(cell, "n_contigs", 0))
    n_contigs <= 0 && return "n/a:empty-assembly"
    largest = _as_number(get(cell, "largest_contig", 0))
    return largest > 0 ? "quast:scored" : "unknown:quast-unscored"
end

_as_number(v::Real) = Float64(v)
function _as_number(v)
    parsed = tryparse(Float64, string(v))
    return parsed === nothing ? 0.0 : parsed
end

# === Discovery ==============================================================

"""
    cells_root(path) -> String

Resolve a `--host` path to the directory that holds per-cell subdirectories.
Accepts either the results root (which contains `cells/`) or `cells/` itself.
"""
function cells_root(path::AbstractString)
    candidate = joinpath(path, "cells")
    isdir(candidate) && return candidate
    return String(path)
end

"""
    discover_cells(root) -> (cells, unreadable)

Scan `root` for `<cell_id>/cell_result.json`.

Returns `cells::Dict{String, Any}` mapping cell id to the parsed checkpoint, and
`unreadable::Vector{Tuple{String, String}}` of `(cell_id, reason)` for
directories whose checkpoint is missing (cell still running) or unparseable.
Both are reported; neither is silently dropped.
"""
function discover_cells(root::AbstractString)
    cells = Dict{String, Any}()
    unreadable = Tuple{String, String}[]
    isdir(root) || return cells, unreadable
    for entry in sort(readdir(root))
        cell_dir = joinpath(root, entry)
        isdir(cell_dir) || continue
        ckpt = joinpath(cell_dir, "cell_result.json")
        if !isfile(ckpt)
            push!(unreadable, (entry, "no cell_result.json (cell in progress?)"))
            continue
        end
        try
            parsed = JSON.parsefile(ckpt)
            parsed isa AbstractDict ||
                error("checkpoint is not a JSON object (got $(typeof(parsed)))")
            cells[entry] = parsed
        catch e
            push!(unreadable, (entry, "unparseable: $(sprint(showerror, e))"))
        end
    end
    return cells, unreadable
end

# === Merge core =============================================================

"""
    merge_hosts(sources; expected_ids) -> NamedTuple

Merge per-host cell checkpoints into one matrix.

`sources` is a vector of `(label, root_path)` pairs, where `root_path` is
anything [`cells_root`](@ref) accepts. `expected_ids` is the matrix the merge is
checked against.

Returns a NamedTuple with:

  * `merged`     — `Dict{String, Any}`: cell id -> parsed checkpoint (winner-free;
                   only cells with no content conflict are included)
  * `origin`     — `Dict{String, String}`: cell id -> producing host label
  * `digests`    — `Dict{String, String}`: cell id -> canonical digest
  * `collisions` — `Vector`: one entry per cell id whose content DIFFERS between
                   hosts, with the per-host digests and the field-level diff
  * `duplicates` — `Vector`: cell ids present on >1 host with IDENTICAL content
  * `missing_ids`, `unexpected_ids` — expected-matrix reconciliation
  * `per_host`   — per-source discovery counts and unreadable entries

A cell involved in a content collision is deliberately EXCLUDED from `merged`:
there is no defensible winner, so it must not enter the analyzed matrix at all.
"""
function merge_hosts(sources::AbstractVector; expected_ids::AbstractVector{<:AbstractString})
    merged = Dict{String, Any}()
    origin = Dict{String, String}()
    digests = Dict{String, String}()
    collisions = Any[]
    duplicates = Any[]
    per_host = Any[]

    # cell id -> Vector of (host, digest, parsed)
    seen = Dict{String, Vector{Tuple{String, String, Any}}}()

    for (label, path) in sources
        root = cells_root(path)
        cells, unreadable = discover_cells(root)
        push!(per_host,
            (host = String(label), root = root, discovered = length(cells),
                unreadable = unreadable))
        for (cell_id, parsed) in cells
            entry = (String(label), cell_digest(parsed), parsed)
            push!(get!(seen, cell_id, Tuple{String, String, Any}[]), entry)
        end
    end

    for cell_id in sort(collect(keys(seen)))
        entries = sort(seen[cell_id]; by = first)
        unique_digests = unique(e[2] for e in entries)
        if length(unique_digests) == 1
            first_host, digest, parsed = entries[1]
            merged[cell_id] = parsed
            origin[cell_id] = first_host
            digests[cell_id] = digest
            if length(entries) > 1
                push!(duplicates,
                    (cell_id = cell_id, hosts = String[e[1] for e in entries],
                        digest = digest))
            end
        else
            # Real anomaly: the shards were supposed to be disjoint.
            pairs = Any[]
            for i in 1:(length(entries) - 1)
                for j in (i + 1):length(entries)
                    entries[i][2] == entries[j][2] && continue
                    push!(pairs,
                        (a = entries[i][1], b = entries[j][1],
                            diffs = differing_fields(entries[i][3], entries[j][3])))
                end
            end
            push!(collisions,
                (cell_id = cell_id,
                    hosts = String[e[1] for e in entries],
                    host_digests = [e[1] => e[2] for e in entries],
                    pairs = pairs))
        end
    end

    expected = Set{String}(expected_ids)
    found = Set{String}(keys(merged))
    # A colliding cell is NOT merged, so it is genuinely absent from the matrix.
    missing_ids = sort(collect(setdiff(expected, found)))
    unexpected_ids = sort(collect(setdiff(Set{String}(keys(seen)), expected)))

    return (merged = merged, origin = origin, digests = digests,
        collisions = collisions, duplicates = duplicates,
        missing_ids = missing_ids, unexpected_ids = unexpected_ids,
        per_host = per_host, expected_n = length(expected))
end

# === Output tables ==========================================================

# Harness row order, so `track_a_results.tsv` written here is drop-in compatible
# with `track_a_baseline_benchmark.jl`'s own aggregate.
const TRACK_A_ROW_KEYS = String[
"organism", "accession", "technology", "coverage", "seed", "decoder_arm", "k",
"n_reads", "n_contigs", "NGA50", "misassemblies", "genome_fraction",
"duplication_ratio", "largest_contig", "wall_seconds", "peak_rss_bytes", "status"
]

"""
    merged_tables(result) -> (results_df, provenance_df)

Build the two output tables from a [`merge_hosts`](@ref) result:

  * `results_df`    — harness schema and column order (drop-in for downstream
                      readers of `track_a_results.tsv`)
  * `provenance_df` — the same rows plus `host`, `source_digest`, and
                      `quast_evidence`, which is what the metric-definition guard
                      consumes
"""
function merged_tables(result)
    ids = sort(collect(keys(result.merged)))
    rows = Any[]
    prov_rows = Any[]
    for cell_id in ids
        cell = result.merged[cell_id]
        base = Dict{String, Any}(k => get(cell, k, missing) for k in TRACK_A_ROW_KEYS)
        push!(rows, NamedTuple{Tuple(Symbol.(TRACK_A_ROW_KEYS))}(
            Tuple(base[k] for k in TRACK_A_ROW_KEYS)))
        prov_keys = vcat(["cell_id"], TRACK_A_ROW_KEYS,
            ["host", "source_digest", "quast_evidence"])
        prov_vals = Any[cell_id]
        append!(prov_vals, Any[base[k] for k in TRACK_A_ROW_KEYS])
        push!(prov_vals, result.origin[cell_id])
        push!(prov_vals, result.digests[cell_id])
        push!(prov_vals, quast_evidence(cell))
        push!(prov_rows,
            NamedTuple{Tuple(Symbol.(prov_keys))}(Tuple(prov_vals)))
    end
    results_df = isempty(rows) ?
                 DataFrames.DataFrame([Symbol(k) => Any[] for k in TRACK_A_ROW_KEYS]) :
                 DataFrames.DataFrame(rows)
    provenance_df = isempty(prov_rows) ?
                    DataFrames.DataFrame(
        [Symbol(k) => Any[]
         for k in vcat(["cell_id"], TRACK_A_ROW_KEYS,
        ["host", "source_digest", "quast_evidence"])]) :
                    DataFrames.DataFrame(prov_rows)
    return results_df, provenance_df
end

"""
    evidence_by_host(result) -> Dict{String, Dict{String, Int}}

Counts of each `quast_evidence` class per host — the "does the merged matrix mix
metric sources by host?" table. QUAST works on Lawrencium but was failing on
Lovelace, so a non-trivial `unknown:quast-unscored` count on one host and not the
other is the expected, and dangerous, shape.
"""
function evidence_by_host(result)
    out = Dict{String, Dict{String, Int}}()
    for (cell_id, cell) in result.merged
        host = result.origin[cell_id]
        per = get!(out, host, Dict{String, Int}())
        cls = quast_evidence(cell)
        per[cls] = get(per, cls, 0) + 1
    end
    return out
end

const EVIDENCE_CLASSES = ("quast:scored", "unknown:quast-unscored",
    "n/a:empty-assembly", "n/a:cell-error")

"""
    write_merge_report(path, result; sources, output_dir) -> String

Human-readable merge report: per-source counts, content collisions with
field-level diffs, benign duplicates, expected-matrix reconciliation (exactly
which cells are missing), and the `quast_evidence` x host breakdown.
"""
function write_merge_report(path::AbstractString, result; output_dir::AbstractString)
    open(path, "w") do io
        println(io, "# Track A cross-host merge report\n")
        println(io, "Generated: $(Dates.now())")
        println(io, "Output dir: `$output_dir`\n")

        println(io, "## Sources\n")
        println(io, "| host | cells root | checkpoints found | not readable |")
        println(io, "| --- | --- | --- | --- |")
        for h in result.per_host
            println(io,
                "| $(h.host) | `$(h.root)` | $(h.discovered) | $(length(h.unreadable)) |")
        end
        println(io)
        for h in result.per_host
            isempty(h.unreadable) && continue
            println(io, "### Not readable on $(h.host)\n")
            for (cell_id, reason) in h.unreadable
                println(io, "- `$cell_id` — $reason")
            end
            println(io)
        end

        n_merged = length(result.merged)
        println(io, "## Merge outcome\n")
        println(io, "- Expected matrix cells: **$(result.expected_n)**")
        println(io, "- Unique cells merged: **$n_merged**")
        println(io, "- Missing from the expected matrix: **$(length(result.missing_ids))**")
        println(io, "- Identical duplicates across hosts (benign): $(length(result.duplicates))")
        println(io, "- **CONTENT COLLISIONS: $(length(result.collisions))** (must be 0)")
        println(io, "- Cell ids outside the expected matrix: $(length(result.unexpected_ids))")
        complete = length(result.missing_ids) == 0 && isempty(result.collisions)
        println(io,
            "\n**Matrix complete and collision-free: $(complete ? "YES" : "NO")**\n")

        if !isempty(result.collisions)
            println(io, "## Content collisions\n")
            println(io,
                "The shards were supposed to be disjoint, so the same cell id with " *
                "DIFFERENT content on two hosts is a real anomaly (overlapping " *
                "shards, a stale checkpoint, or two different Mycelia SHAs). These " *
                "cells are EXCLUDED from the merged matrix — there is no defensible " *
                "winner to pick.\n")
            for c in result.collisions
                println(io, "### `$(c.cell_id)`\n")
                for (host, digest) in c.host_digests
                    println(io, "- $host: `$(first(digest, 16))…`")
                end
                for p in c.pairs
                    println(io, "\n$(p.a) vs $(p.b):\n")
                    if isempty(p.diffs)
                        println(io,
                            "- (digests differ but no field-level diff — " *
                            "check for key-set differences)")
                    else
                        for d in p.diffs
                            println(io, "- $d")
                        end
                    end
                end
                println(io)
            end
        end

        if !isempty(result.duplicates)
            println(io, "## Identical duplicates across hosts\n")
            println(io,
                "Same cell id, byte-equivalent content. Benign (the shard overlapped " *
                "but agreed); recorded so the overlap is visible.\n")
            for d in result.duplicates
                println(io, "- `$(d.cell_id)` on $(join(d.hosts, ", "))")
            end
            println(io)
        end

        if !isempty(result.missing_ids)
            println(io, "## Missing cells ($(length(result.missing_ids)))\n")
            by_org = Dict{String, Vector{String}}()
            for id in result.missing_ids
                org = first(split(id, "__"))
                push!(get!(by_org, org, String[]), id)
            end
            for org in sort(collect(keys(by_org)))
                ids = by_org[org]
                println(io, "### $org ($(length(ids)))\n")
                for id in ids
                    println(io, "- `$id`")
                end
                println(io)
            end
        end

        if !isempty(result.unexpected_ids)
            println(io, "## Cell ids outside the expected matrix\n")
            println(io,
                "Either the expected matrix was narrowed incorrectly, or the " *
                "harness's cell-id format drifted from this tool's mirror of it.\n")
            for id in result.unexpected_ids
                println(io, "- `$id`")
            end
            println(io)
        end

        println(io, "## Metric provenance: quast_evidence x host\n")
        println(io,
            "QUAST works on Lawrencium but was failing on Lovelace (conda " *
            "`ProcessExited(4)`), where the harness degrades to internal metrics. A " *
            "non-trivial `unknown:quast-unscored` count on one host and not another " *
            "means the merged matrix MIXES METRIC DEFINITIONS BY HOST — see bead " *
            "td-9p91 and `metric_source_guard.jl` before aggregating.\n")
        breakdown = evidence_by_host(result)
        println(io, "| host | " * join(EVIDENCE_CLASSES, " | ") * " |")
        println(io, "| --- |" * repeat(" --- |", length(EVIDENCE_CLASSES)))
        for host in sort(collect(keys(breakdown)))
            per = breakdown[host]
            println(io,
                "| $host | " *
                join((string(get(per, cls, 0)) for cls in EVIDENCE_CLASSES), " | ") * " |")
        end
        println(io)
        classes_present = unique(vcat(
            [collect(keys(v)) for v in values(breakdown)]...))
        scoring = filter(c -> !startswith(c, "n/a:"), classes_present)
        if length(scoring) > 1
            println(io,
                "> **WARNING:** more than one scoring class present " *
                "($(join(scoring, ", "))). Do not aggregate across them without " *
                "`assert_single_metric_definition`.\n")
        end
    end
    return path
end

"""
    copy_merged_cells(result, output_dir) -> (copied, skipped, conflicts)

Materialize `<output_dir>/cells/<cell_id>/cell_result.json` for every merged cell
so a single tree can be analyzed (or resumed) as one matrix.

Idempotent: a destination that is already content-identical is skipped. A
destination that DIFFERS from the source is never overwritten — it is returned as
a conflict, because the previously merged state disagreeing with the source is
the same class of anomaly as a cross-host collision.
"""
function copy_merged_cells(result, output_dir::AbstractString)
    dest_root = joinpath(output_dir, "cells")
    mkpath(dest_root)
    copied = String[]
    skipped = String[]
    conflicts = Any[]
    for (cell_id, cell) in result.merged
        cell_dir = joinpath(dest_root, cell_id)
        dest = joinpath(cell_dir, "cell_result.json")
        if isfile(dest)
            existing = try
                JSON.parsefile(dest)
            catch
                nothing
            end
            if existing isa AbstractDict && cell_digest(existing) == result.digests[cell_id]
                push!(skipped, cell_id)
                continue
            end
            push!(conflicts,
                (cell_id = cell_id, hosts = ["<merged-output>", result.origin[cell_id]],
                    host_digests = [
                        "<merged-output>" => (existing isa AbstractDict ?
                                              cell_digest(existing) : "unparseable"),
                        result.origin[cell_id] => result.digests[cell_id]],
                    pairs = existing isa AbstractDict ?
                            [(a = "<merged-output>", b = result.origin[cell_id],
                        diffs = differing_fields(existing, cell))] : Any[]))
            continue
        end
        mkpath(cell_dir)
        open(dest, "w") do io
            JSON.print(io, cell, 2)
        end
        push!(copied, cell_id)
    end
    return sort(copied), sort(skipped), conflicts
end

"""
    merge_exit_status(result; allow_collisions=false, allow_incomplete=false)
        -> (code, problems)

The merge's pass/fail decision, factored out of the CLI so it is unit-testable
(the exit code IS the contract for callers that chain a merge into an analysis).

Returns `code` (0 pass, 1 fail) and `problems`, a vector of
`(; kind, message, fatal)` describing each condition found. A condition that was
waived by its `--allow-*` flag still appears, with `fatal = false`, so a waiver is
never invisible.
"""
function merge_exit_status(result;
        allow_collisions::Bool = false, allow_incomplete::Bool = false,
        report_hint::AbstractString = "the merge report")
    problems = Any[]
    if !isempty(result.collisions)
        push!(problems,
            (kind = :collisions, fatal = !allow_collisions,
                message = "$(length(result.collisions)) content collision(s) — same " *
                          "cell id, different content across hosts. See $report_hint."))
    end
    if !isempty(result.missing_ids)
        push!(problems,
            (kind = :incomplete, fatal = !allow_incomplete,
                message = "matrix incomplete: $(length(result.missing_ids)) of " *
                          "$(result.expected_n) cells missing. See $report_hint."))
    end
    code = any(p -> p.fatal, problems) ? 1 : 0
    return code, problems
end

# === CLI ====================================================================

function _arg_all(flag)
    String[ARGS[i + 1]
           for i in eachindex(ARGS)
           if ARGS[i] == flag && i < length(ARGS)]
end

function _arg_value(flag)
    i = findfirst(==(flag), ARGS)
    return (i !== nothing && i < length(ARGS)) ? ARGS[i + 1] : nothing
end

function _arg_list(flag)
    v = _arg_value(flag)
    return v === nothing ? nothing : String.(split(v, ","))
end

function _parse_host_spec(spec::AbstractString)
    idx = findfirst('=', spec)
    idx === nothing && error(
        "--host expects LABEL=PATH, got \"$spec\" (e.g. --host lrc=/scratch/track_a_baseline)")
    label = strip(spec[1:(idx - 1)])
    path = strip(spec[(idx + 1):end])
    isempty(label) && error("--host label is empty in \"$spec\"")
    isempty(path) && error("--host path is empty in \"$spec\"")
    return String(label) => String(path)
end

function main()
    specs = _arg_all("--host")
    if isempty(specs)
        println(stderr,
            "ERROR: at least one --host LABEL=PATH is required.\n" *
            "  julia --project=. benchmarking/track_a_merge_hosts.jl \\\n" *
            "      --host lovelace=/path/to/track_a_baseline \\\n" *
            "      --host lrc=/path/to/track_a_baseline \\\n" *
            "      --output-dir benchmarking/results/track_a_baseline_merged")
        return 2
    end
    sources = [_parse_host_spec(s) for s in specs]
    labels = String[first(s) for s in sources]
    if length(unique(labels)) != length(labels)
        println(stderr, "ERROR: duplicate --host labels: $(join(labels, ", "))")
        return 2
    end

    output_dir = something(_arg_value("--output-dir"),
        joinpath(@__DIR__, "results", "track_a_baseline_merged"))
    copy_cells = !("--no-copy-cells" in ARGS)
    allow_collisions = "--allow-collisions" in ARGS
    allow_incomplete = "--allow-incomplete" in ARGS
    assert_single_def = "--assert-single-metric-definition" in ARGS

    organisms = something(_arg_list("--organisms"), TRACK_A_ORGANISMS)
    technologies = something(_arg_list("--technologies"), TRACK_A_TECHNOLOGIES)
    coverages = let v = _arg_list("--coverages")
        v === nothing ? TRACK_A_COVERAGES : parse.(Int, v)
    end
    seeds = let v = _arg_list("--seeds")
        v === nothing ? TRACK_A_SEEDS : parse.(Int, v)
    end
    arms = something(_arg_list("--arms"), TRACK_A_ARMS)

    ids = expected_cell_ids(; organisms, technologies, coverages, seeds, arms)

    println("=== Track A cross-host merge ===")
    println("Sources:")
    for (label, path) in sources
        println("  $label -> $path")
    end
    println("Expected matrix: $(length(ids)) cells")
    println("Output dir     : $output_dir")

    result = merge_hosts(sources; expected_ids = ids)

    mkpath(output_dir)
    results_df, provenance_df = merged_tables(result)
    results_path = joinpath(output_dir, "track_a_results.tsv")
    provenance_path = joinpath(output_dir, "track_a_results_provenance.tsv")
    CSV.write(results_path, results_df; delim = '\t')
    CSV.write(provenance_path, provenance_df; delim = '\t')

    if copy_cells
        copied, skipped, dest_conflicts = copy_merged_cells(result, output_dir)
        println("Cells materialized: $(length(copied)) new, $(length(skipped)) already identical")
        if !isempty(dest_conflicts)
            append!(result.collisions, dest_conflicts)
            println(stderr,
                "WARNING: $(length(dest_conflicts)) previously-merged cell(s) differ " *
                "from their source and were NOT overwritten.")
        end
    end

    report_path = write_merge_report(
        joinpath(output_dir, "merge_report.md"), result; output_dir = output_dir)

    println("\n--- Summary ---")
    for h in result.per_host
        println("  $(h.host): $(h.discovered) checkpoints, $(length(h.unreadable)) not readable")
    end
    println("  merged: $(length(result.merged)) / $(length(ids)) expected")
    println("  identical duplicates: $(length(result.duplicates))")
    println("  CONTENT COLLISIONS:   $(length(result.collisions))")
    println("  missing:              $(length(result.missing_ids))")
    println("  unexpected ids:       $(length(result.unexpected_ids))")
    println("\nWrote:\n  $results_path\n  $provenance_path\n  $report_path")

    rc,
    problems = merge_exit_status(result;
        allow_collisions = allow_collisions, allow_incomplete = allow_incomplete,
        report_hint = report_path)
    for p in problems
        if p.fatal
            println(stderr, "ERROR: $(p.message)")
        else
            @warn "WAIVED via --allow-$(p.kind == :collisions ? "collisions" : "incomplete"). $(p.message)"
        end
    end
    if assert_single_def
        try
            assert_single_metric_definition(provenance_df;
                columns = (:quast_evidence,),
                context = "merged Track A matrix")
            println("\nmetric-definition check: single quast_evidence class — OK")
        catch e
            println(stderr, "\nERROR: metric-definition check failed:")
            showerror(stderr, e)
            println(stderr)
            rc = 1
        end
    end
    return rc
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
