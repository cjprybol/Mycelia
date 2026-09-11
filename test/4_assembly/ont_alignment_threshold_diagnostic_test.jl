# Contract tests for the ONT alignment-threshold diagnostic (td-4e19d.28).
#
# The defect these pin is a DATA-INTEGRITY one, not a crash. The script maps
# every QUAST failure to all-missing metrics, and an all-missing row is exactly
# what its own INTERPRETATION block reads as "genome fraction stays near zero
# even at 80%" — i.e. as evidence for an assembler defect, the most consequential
# hypothesis the sweep tests. A broken conda env would have supplied that
# evidence while looking entirely normal in the output table.
#
# Run:
#   julia --project=. test/4_assembly/ont_alignment_threshold_diagnostic_test.jl

# Wrapped in a module: runtests.jl includes every test file into one shared
# `Main`, and this script (via ont_k_sweep.jl) defines top-level consts that
# collide with the other benchmarking harnesses already included there.
module OntAlignmentThresholdDiagnosticTest

import Test
import JSON

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "ont_alignment_threshold_diagnostic.jl"))

Test.@testset "ONT alignment-threshold diagnostic helpers" begin
    Test.@testset "a QUAST failure is distinguishable from a no-alignment result" begin
        # Exit 0 with a report is the only unambiguous measurement. Everything
        # else has to be excluded before the table is read as evidence, because
        # all four cases produce the same all-missing metrics.
        Test.@test rescore_status_for(false, true) == "ok"
        Test.@test rescore_status_for(true, true) == "nonzero_with_report"
        Test.@test rescore_status_for(true, false) == "nonzero_no_report"
        Test.@test rescore_status_for(false, false) == "no_report"

        # Only the two report-bearing cases carry metrics at all, so only those
        # can be read as a statement about the assembly.
        Test.@test "ok" in INTERPRETABLE_RESCORE_STATUSES
        Test.@test "nonzero_with_report" in INTERPRETABLE_RESCORE_STATUSES
        Test.@test !("nonzero_no_report" in INTERPRETABLE_RESCORE_STATUSES)
        Test.@test !("no_report" in INTERPRETABLE_RESCORE_STATUSES)

        # Every status must be classifiable — a value outside the set below would
        # silently fall through the interpretation filter.
        all_statuses = [rescore_status_for(a, b)
                        for a in (true, false), b in (true, false)]
        Test.@test length(unique(all_statuses)) == 4
    end

    Test.@testset "candidate_cells skips an unreadable checkpoint" begin
        # A checkpoint truncated by a crash mid-write used to abort the whole
        # diagnostic from JSON.parsefile, rather than costing one cell.
        mktempdir() do dir
            cells = joinpath(dir, "cells")
            mkpath(cells)

            good_id = cell_id_for("Lambda", "ont", 21, 30, 42)
            good = cell_row("Lambda", "NC_001416", "ont", 21, 30, 42;
                n_reads = 10, asm = contig_stats(["A"^600], MIN_CONTIG),
                metrics = merge(empty_metrics(),
                    (; genome_fraction = 31.663, quast_contigs = 1.0)),
                nga50_status = "censored_partial_alignment",
                outcome = "degenerate", wall_seconds = 1.0, status = "ok")
            mkpath(joinpath(cells, good_id))
            save_cell_json(joinpath(cells, good_id, "cell_result.json"), good)

            bad_id = cell_id_for("Lambda", "ont", 31, 30, 42)
            mkpath(joinpath(cells, bad_id))
            write(joinpath(cells, bad_id, "cell_result.json"), "{ truncated")

            found = candidate_cells(dir)
            Test.@test length(found) == 1
            Test.@test found[1]["cell_id"] == good_id
            # Reclassification still happens on the surviving cell.
            Test.@test found[1]["nga50_status"] == "censored_partial_alignment"
        end
    end
end

end  # module
