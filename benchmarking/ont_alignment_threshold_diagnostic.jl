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
  * `censored_gf_below_50`   — contigs aligned, but under the 50% genome
    fraction floor where NGA50 is undefined. Relaxing the threshold asks
    whether more of the genome comes into alignment range, which is the
    quantity that determines whether NGA50 becomes defined at all.

Cells with `no_contigs_ge_min` are deliberately EXCLUDED — they have nothing at
or above --min-contig for QUAST to align at any identity threshold, so relaxing
the threshold cannot tell us anything about them. Their failure is contig
LENGTH, which is a different mechanism and is already quantified by the sweep's
`asm_max_contig` column.
"""
const CENSORED_STATUSES = ("censored_no_alignment", "censored_gf_below_50")

function candidate_cells(sweep_dir)
    cells = Dict{String, Any}[]
    cells_dir = joinpath(sweep_dir, "cells")
    isdir(cells_dir) || return cells
    for entry in sort(readdir(cells_dir))
        checkpoint = joinpath(cells_dir, entry, "cell_result.json")
        isfile(checkpoint) || continue
        row = JSON.parsefile(checkpoint)
        row["cell_id"] = entry
        row["cell_dir"] = joinpath(cells_dir, entry)
        push!(cells, row)
    end
    return cells
end

"""
    rescore(contigs, reference, outdir, min_identity) -> NamedTuple

Run QUAST over `contigs` at an explicit `--min-identity`, returning the parsed
metrics (all-missing if QUAST declines or errors).

`--min-contig 500` matches the sweep and the Track-A pilot, so the only variable
between rows of the output is the identity threshold.
"""
function rescore(contigs, reference, outdir, min_identity)
    mkpath(outdir)
    Mycelia.add_bioconda_env("quast")
    try
        run(pipeline(
            `$(Mycelia.CONDA_RUNNER) run --live-stream -n quast quast.py
             --output-dir $(outdir) --threads 1 --min-contig 500
             --min-identity $(min_identity) --reference $(reference) $(contigs)`,
            stdout = devnull, stderr = devnull))
    catch
        # QUAST exits non-zero when nothing survives filtering. That is a
        # result, not an error: it means no contig aligned even at this
        # threshold, and the all-missing metrics below say exactly that.
        return empty_metrics()
    end
    report = joinpath(outdir, "report.tsv")
    return isfile(report) ? parse_quast_metrics(report) : empty_metrics()
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("=== ONT alignment-threshold diagnostic (td-4e19d.28) ===")
    println("Start: $(Dates.now())")
    println("Identity ladder: $(join(IDENTITIES, ", "))%")

    reference = joinpath(SWEEP_DIR, "refs", "NC_001416.fna")
    isfile(reference) || error("reference not found: $(reference). Run the sweep first.")
    mkpath(OUT_DIR)

    all_cells = candidate_cells(SWEEP_DIR)
    selected = filter(all_cells) do cell
        CELL_FILTER !== nothing && return cell["cell_id"] in CELL_FILTER
        cell["nga50_status"] in CENSORED_STATUSES
    end
    println("Cells available: $(length(all_cells)); selected for rescoring: $(length(selected))")
    isempty(selected) && println("  (nothing to do — no censored_unaligned cells yet)")

    rows = NamedTuple[]
    for cell in selected
        contigs = joinpath(cell["cell_dir"], "contigs.fasta")
        if !isfile(contigs)
            @warn "contigs missing; skipping" cell = cell["cell_id"]
            continue
        end
        for min_identity in IDENTITIES
            outdir = joinpath(OUT_DIR, cell["cell_id"], "idy$(min_identity)")
            metrics = rescore(contigs, reference, outdir, min_identity)
            push!(rows,
                (
                    cell_id = cell["cell_id"], technology = cell["technology"],
                    k = cell["k"], coverage = cell["coverage"], seed = cell["seed"],
                    min_identity = min_identity,
                    asm_contigs_ge_min = cell["asm_contigs_ge_min"],
                    asm_max_contig = cell["asm_max_contig"],
                    genome_fraction = metrics.genome_fraction,
                    NGA50 = metrics.NGA50,
                    NA50 = metrics.NA50,
                    largest_alignment = metrics.largest_alignment,
                    unaligned_length = metrics.unaligned_length,
                    misassemblies = metrics.misassemblies))
            println("  $(cell["cell_id"]) @ IDY $(min_identity)%: " *
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
    end
    println("End: $(Dates.now())")
end  # PROGRAM_FILE guard
