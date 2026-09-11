# Is ONT "nothing aligned" a threshold effect or a broken assembly? (td-4e19d.28)
#
# WHY THIS EXISTS
# ---------------
# QUAST scores contigs with a minimum alignment identity, and its default is
# 95.0% (visible in every quast.log as `min alignment IDY: 95.0`). The measured
# identity of these simulated ONT reads is 94.4% BLAST (see
# benchmarking/ont_read_identity.jl and results/ont_read_identity/).
#
# Those two numbers straddle. A naked k-mer assembler with no consensus or
# polishing step emits contigs by walking observed k-mers, so contig error rate
# is inherited from the reads rather than averaged away across them. Such
# contigs should land at roughly read identity — just BELOW QUAST's default
# threshold. If that is what is happening, then "NGA50 could not be computed"
# means "aligned at ~94%, which is under 95%", NOT "produced unrelated
# sequence", and those are very different claims about the assembler.
#
# This rescores already-assembled contigs at a ladder of identity thresholds.
# It runs NO new assembly — it reuses the contigs.fasta each sweep cell already
# wrote, so it cannot perturb the sweep's results and costs only QUAST time.
#
# INTERPRETATION
#
# FILTER ON `rescore_status` FIRST. A row whose status is not `ok` or
# `nonzero_with_report` carries all-missing metrics because QUAST produced no
# report — not because nothing aligned. Read unfiltered, such a row is
# indistinguishable from the strongest evidence this script can produce for an
# assembler defect, which is the most consequential hypothesis the sweep tests.
# Both readings below assume the uninterpretable rows have been dropped.
#
#   * Genome fraction rises sharply as the threshold is relaxed  -> the assembly
#     recovers the genome at an accuracy below QUAST's default cut. The
#     degeneracy is a real property of uncorrected k-mer assembly on long reads
#     (no consensus step), and the censored NGA50 is that property meeting a
#     scoring threshold.
#   * Genome fraction stays near zero even at 80%  -> the contigs are not
#     approximate reconstructions at all, which would point at a defect rather
#     than an accuracy limit.
#
# This script deliberately bypasses `Mycelia.run_quast`, which does not expose
# --min-identity. Bypassing it is the point: the parameter under test is exactly
# the one the wrapper pins to its default.
#
# Usage:
#   julia --project=. benchmarking/ont_alignment_threshold_diagnostic.jl
#   julia --project=. benchmarking/ont_alignment_threshold_diagnostic.jl --identities 95,90,85,80
#   julia --project=. benchmarking/ont_alignment_threshold_diagnostic.jl --cells Lambda__ont__k31__30x__seed42

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import CSV
import DataFrames
import Dates
import JSON

include(joinpath(@__DIR__, "ont_k_sweep.jl"))  # parse_quast_metrics, cell_id_for

const SWEEP_DIR = something(arg_value("--sweep-dir"),
    joinpath(@__DIR__, "results", "ont_k_sweep"))
const OUT_DIR = something(arg_value("--output-dir"),
    joinpath(@__DIR__, "results", "ont_alignment_threshold"))
const IDENTITIES = let v = arg_value("--identities")
    v === nothing ? [95.0, 90.0, 85.0, 80.5] : parse.(Float64, split(v, ","))
end
const CELL_FILTER = arg_list("--cells")

"""
    candidate_cells(sweep_dir) -> Vector{Dict}

Sweep cells worth rescoring: those that produced contigs long enough to score
but where QUAST could not compute NGA50. That is BOTH censoring causes:

  * `censored_no_alignment`  — nothing aligned at the default 95% identity.
    Relaxing the threshold asks whether the contigs are approximate
    reconstructions that merely fell under the cut.
  * `censored_partial_alignment`   — contigs aligned, but under the 50% genome
    fraction floor where NGA50 is undefined. Relaxing the threshold asks
    whether more of the genome comes into alignment range, which is the
    quantity that determines whether NGA50 becomes defined at all.

Cells with `no_contigs_ge_min` are deliberately EXCLUDED — they have nothing at
or above --min-contig for QUAST to align at any identity threshold, so relaxing
the threshold cannot tell us anything about them. Their failure is contig
LENGTH, which is a different mechanism and is already quantified by the sweep's
`asm_max_contig` column.
"""
const CENSORED_STATUSES = ("censored_no_alignment", "censored_partial_alignment")

function candidate_cells(sweep_dir)
    cells = Dict{String, Any}[]
    cells_dir = joinpath(sweep_dir, "cells")
    isdir(cells_dir) || return cells
    for entry in sort(readdir(cells_dir))
        checkpoint = joinpath(cells_dir, entry, "cell_result.json")
        isfile(checkpoint) || continue
        # Reclassify on read rather than trusting the stored label. The sweep's
        # `reclassify` exists precisely because stored labels go stale — an
        # earlier classifier wrote `censored_unaligned`, which matches neither
        # entry of CENSORED_STATUSES, so such a checkpoint would be silently
        # dropped and this script would report "nothing to do" over a tree full
        # of censored cells.
        # Guarded for the same reason `load_all_checkpoints` guards its read: a
        # checkpoint truncated by a crash mid-write makes JSON.parsefile throw,
        # and an unguarded throw here would abort the whole diagnostic over one
        # unreadable cell instead of skipping it.
        parsed = try
            JSON.parsefile(checkpoint)
        catch e
            @warn "unreadable checkpoint; skipping cell" cell=entry exception=e
            continue
        end
        row = Dict{String, Any}(parsed)
        try
            row["nga50_status"] = canonical(parsed).nga50_status
        catch e
            @warn "could not reclassify checkpoint; using stored label" cell=entry exception=e
        end
        row["cell_id"] = entry
        row["cell_dir"] = joinpath(cells_dir, entry)
        push!(cells, row)
    end
    return cells
end

"""
    rescore_status_for(exited_nonzero, has_report) -> String

How a QUAST invocation terminated, from its exit status and whether it left a
report. Split out from `rescore` so the classification is reachable without a
QUAST installation; the meaning of each value is documented on `rescore`.
"""
function rescore_status_for(exited_nonzero::Bool, has_report::Bool)
    if !exited_nonzero
        return has_report ? "ok" : "no_report"
    else
        return has_report ? "nonzero_with_report" : "nonzero_no_report"
    end
end

"""
    rescore(contigs, reference, outdir, min_identity) -> (; metrics, status)

Run QUAST over `contigs` at an explicit `--min-identity`, returning the parsed
metrics AND how the run terminated.

`--min-contig 500` matches the sweep and the Track-A pilot, so the only variable
between rows of the output is the identity threshold.

The status is returned, and persisted per row, because an all-missing metric row
has two completely different meanings and they are otherwise indistinguishable:

  * QUAST ran and nothing aligned at this threshold — a MEASUREMENT, and the
    input to this script's INTERPRETATION block.
  * QUAST could not run — missing conda env, OOM, an out-of-range
    `--min-identity` (QUAST rejects < 80.0), an unreadable reference. Not a
    measurement of anything.

Collapsing the second into the first is how a tool failure becomes evidence:
"genome fraction stays near zero even at 80%" is exactly the observation the
header reads as support for an assembler defect, the most consequential
hypothesis the sweep tests. A row that reached this table only because QUAST
was broken would supply that evidence while looking entirely normal.

Statuses:

  * `ok`                   — exited 0 and wrote report.tsv. Metrics are a
                             measurement; all-missing means nothing aligned.
  * `nonzero_with_report`  — exited non-zero but still emitted a report. Metrics
                             are parseable but the run complained.
  * `nonzero_no_report`    — exited non-zero and wrote nothing. THIS is the
                             ambiguous class: "nothing survived filtering" and
                             "QUAST is broken" look identical from here. Exclude
                             from threshold interpretation.
  * `no_report`            — exited 0 and wrote nothing. Anomalous; treat as
                             `nonzero_no_report`.
"""
function rescore(contigs, reference, outdir, min_identity)
    mkpath(outdir)
    Mycelia.add_bioconda_env("quast")
    exited_nonzero = false
    try
        run(pipeline(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n quast quast.py
             --output-dir $(outdir) --threads 1 --min-contig 500
             --min-identity $(min_identity) --reference $(reference) $(contigs)`,
            stdout = devnull, stderr = devnull))
    catch e
        exited_nonzero = true
        @warn "QUAST returned non-zero at this threshold; the row is recorded " *
              "with a non-ok rescore_status and must be excluded from threshold " *
              "interpretation. If this fires at EVERY threshold for a cell, " *
              "check it is not an environment or argument failure." min_identity exception=e
    end
    report = joinpath(outdir, "report.tsv")
    has_report = isfile(report)
    status = rescore_status_for(exited_nonzero, has_report)
    metrics = has_report ? parse_quast_metrics(report) : empty_metrics()
    return (; metrics, status)
end

# Rows whose metrics are a measurement rather than a report about the tool. Only
# these may be read as evidence about the assembly.
const INTERPRETABLE_RESCORE_STATUSES = ("ok", "nonzero_with_report")

if abspath(PROGRAM_FILE) == @__FILE__
    println("=== ONT alignment-threshold diagnostic (td-4e19d.28) ===")
    println("Start: $(Dates.now())")
    println("Identity ladder: $(join(IDENTITIES, ", "))%")

    # Reference PER ORGANISM. This was hardcoded to Lambda's NC_001416.fna, which
    # silently rescored every T4 cell against the wrong genome and produced a
    # full set of well-formed all-NA rows. Those NAs are not harmless: this
    # script's own INTERPRETATION block reads "genome fraction stays near zero
    # even at 80%" as evidence that the contigs are not approximate
    # reconstructions at all — i.e. as support for an assembler defect, the most
    # consequential hypothesis the sweep tests. A wrong-reference comparison
    # would have supplied that evidence while looking entirely normal.
    references = Dict{String, String}()
    for (org, acc, _size) in ORGANISMS
        path = joinpath(SWEEP_DIR, "refs", "$(acc).fna")
        isfile(path) && (references[org] = path)
    end
    isempty(references) &&
        error("no reference FASTA found under $(joinpath(SWEEP_DIR, "refs")). " *
              "Run the sweep first.")
    mkpath(OUT_DIR)

    all_cells = candidate_cells(SWEEP_DIR)
    selected = filter(all_cells) do cell
        CELL_FILTER !== nothing && return cell["cell_id"] in CELL_FILTER
        cell["nga50_status"] in CENSORED_STATUSES
    end
    println("Cells available: $(length(all_cells)); selected for rescoring: $(length(selected))")
    isempty(selected) &&
        println("  (nothing to do — no censored cells in this sweep tree)")

    rows = NamedTuple[]
    for cell in selected
        contigs = joinpath(cell["cell_dir"], "contigs.fasta")
        if !isfile(contigs)
            @warn "contigs missing; skipping" cell = cell["cell_id"]
            continue
        end
        organism = cell["organism"]
        if !haskey(references, organism)
            @warn "no reference for organism; skipping cell" organism cell=cell["cell_id"]
            continue
        end
        for min_identity in IDENTITIES
            outdir = joinpath(OUT_DIR, cell["cell_id"], "idy$(min_identity)")
            result = rescore(contigs, references[organism], outdir, min_identity)
            metrics = result.metrics
            push!(rows,
                (
                    cell_id = cell["cell_id"], organism = cell["organism"],
                    technology = cell["technology"],
                    k = cell["k"], coverage = cell["coverage"], seed = cell["seed"],
                    min_identity = min_identity,
                    asm_contigs_ge_min = cell["asm_contigs_ge_min"],
                    asm_max_contig = cell["asm_max_contig"],
                    rescore_status = result.status,
                    genome_fraction = metrics.genome_fraction,
                    NGA50 = metrics.NGA50,
                    NA50 = metrics.NA50,
                    largest_alignment = metrics.largest_alignment,
                    unaligned_length = metrics.unaligned_length,
                    misassemblies = metrics.misassemblies))
            println("  $(cell["cell_id"]) @ IDY $(min_identity)%: " *
                    "[$(result.status)] " *
                    "GF=$(ismissing(metrics.genome_fraction) ? "NA" : metrics.genome_fraction)% " *
                    "NGA50=$(ismissing(metrics.NGA50) ? "NA" : metrics.NGA50) " *
                    "largest_aln=$(ismissing(metrics.largest_alignment) ? "NA" : metrics.largest_alignment)")
            flush(stdout)
        end
    end

    if !isempty(rows)
        df = DataFrames.DataFrame(rows)
        sort!(df, [:technology, :k, :coverage, :seed, :min_identity])
        CSV.write(joinpath(OUT_DIR, "alignment_threshold_diagnostic.tsv"), df;
            delim = '\t', missingstring = "NA")
        println("\nWrote $(joinpath(OUT_DIR, "alignment_threshold_diagnostic.tsv"))")

        # State the interpretable/uninterpretable split at the point of use. An
        # uninterpretable row is all-missing and therefore reads, to the naked
        # eye and to any consumer that does not filter, as "nothing aligned even
        # at 80%" — this script's stated evidence for an assembler defect.
        n_bad = count(r -> !(r.rescore_status in INTERPRETABLE_RESCORE_STATUSES),
            eachrow(df))
        println("Rows: $(DataFrames.nrow(df)) total, " *
                "$(DataFrames.nrow(df) - n_bad) interpretable, $(n_bad) not.")
        if n_bad > 0
            @warn "rows with a non-interpretable rescore_status are present. " *
                  "They carry all-missing metrics because QUAST did not " *
                  "produce a report, NOT because nothing aligned. Filter on " *
                  "rescore_status before reading this table as evidence." n_bad
        end
    end
    println("End: $(Dates.now())")
end  # PROGRAM_FILE guard
