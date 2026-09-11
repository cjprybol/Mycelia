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

    Test.@testset "write_threshold_table refuses to shrink the committed table" begin
        # td-4blm in its worse form. This script had no union at all — it wrote
        # only the rows the current invocation computed — and --output-dir
        # defaults to the git-tracked results directory. Its OWN documented
        # usage line, `--cells Lambda__ont__k31__30x__seed42`, would therefore
        # have replaced the committed 152-row table with 4 rows. Rescoring is
        # expensive enough that narrowing is the normal way to run it, so the
        # truncating shape was the common one, not the exotic one.
        #
        # Conditional, not past tense: git history shows both committed tables
        # only ever GREW (the diagnostic table 57 -> 153 lines, the sweep table
        # 97 -> 241). Nothing establishes the loss actually occurred.
        row(cell,
            k,
            seed,
            idy) = (
            cell_id = cell, organism = "Lambda", technology = "ont",
            k = k, coverage = 30, seed = seed, min_identity = idy,
            asm_contigs_ge_min = 1, asm_max_contig = 600,
            rescore_status = "ok", genome_fraction = 31.663,
            NGA50 = missing, NA50 = missing, largest_alignment = 400,
            unaligned_length = 200, misassemblies = 0)
        # 80.5 is in the real default ladder and is the float most likely to
        # expose a key-normalisation bug, so it stays in the fixture.
        idys = (95.0, 90.0, 85.0, 80.5)
        a = cell_id_for("Lambda", "ont", 31, 30, 42)
        b = cell_id_for("Lambda", "ont", 21, 30, 42)
        full = vcat([row(a, 31, 42, i) for i in idys],
            [row(b, 21, 42, i) for i in idys])
        table_of(dir) = joinpath(dir, "alignment_threshold_diagnostic.tsv")
        nrows(dir) = DataFrames.nrow(CSV.read(table_of(dir),
            DataFrames.DataFrame; delim = '\t', missingstring = "NA"))

        # Seed with a plain CSV.write, NOT through write_threshold_table. A
        # fixture built by the symbol under test cannot run against pre-change
        # code at all — it fails to resolve the name — so every assertion below
        # it would be an artifact of symbol resolution rather than a statement
        # about behaviour. The sweep testset was corrected for exactly this and
        # this file was missed.
        seed_table(dir) = CSV.write(table_of(dir), DataFrames.DataFrame(full);
            delim = '\t', missingstring = "NA")

        mktempdir() do dir
            seed_table(dir)
            Test.@test nrows(dir) == 8

            # The documented --cells example: one cell, four thresholds.
            one_cell = filter(r -> r.cell_id == a, full)
            err = try
                write_threshold_table(dir, one_cell)
                nothing
            catch e
                e
            end
            Test.@test err isa ErrorException
            Test.@test occursin("refusing to shrink", err.msg)
            Test.@test occursin("--allow-shrink", err.msg)
            # The committed table survives the attempt intact.
            Test.@test nrows(dir) == 8

            # An empty run writes nothing at all rather than a headerless file
            # over a populated table.
            Test.@test write_threshold_table(dir, NamedTuple[]) === nothing
            Test.@test nrows(dir) == 8

            # Rewriting the same key set is not a shrink, so a legitimate
            # full re-run still lands.
            Test.@test write_threshold_table(dir, reverse(full)) !== nothing
            Test.@test nrows(dir) == 8
        end
    end

    Test.@testset "preflight refuses before any QUAST work is spent" begin
        # write_threshold_table is the backstop; this is the check that matters
        # operationally. The output key set is exactly selected x identities and
        # both are known before the rescoring loop, so an unpublishable run is
        # detectable up front. Without this the refusal lands AFTER every QUAST
        # invocation has completed — 152 of them for the committed cell set —
        # and the rows are then discarded in memory, with no per-cell
        # checkpoints to salvage the work.
        idys = (95.0, 90.0, 85.0, 80.5)
        a = cell_id_for("Lambda", "ont", 31, 30, 42)
        b = cell_id_for("Lambda", "ont", 21, 30, 42)
        committed = DataFrames.DataFrame(
            cell_id = [c for c in (a, b) for _ in idys],
            min_identity = [i for _ in (a, b) for i in idys],
            value = 1:8)
        cell_of(id) = Dict{String, Any}("cell_id" => id)

        mktempdir() do dir
            target = joinpath(dir, "alignment_threshold_diagnostic.tsv")
            CSV.write(target, committed; delim = '\t', missingstring = "NA")

            # The documented --cells invocation: one of the two cells.
            err = try
                preflight_threshold_table(dir, [cell_of(a)], idys)
                nothing
            catch e
                e
            end
            Test.@test err isa ErrorException
            Test.@test occursin("refusing to shrink", err.msg)

            # A narrowed --identities ladder is the other narrowing flag, and
            # reaches the same refusal by dropping half of every cell's rows.
            err2 = try
                preflight_threshold_table(dir, [cell_of(a), cell_of(b)],
                    (95.0, 90.0))
                nothing
            catch e
                e
            end
            Test.@test err2 isa ErrorException
            Test.@test occursin("refusing to shrink", err2.msg)

            # The full run passes preflight, so the check does not block the
            # invocation it is meant to permit.
            Test.@test preflight_threshold_table(
                dir, [cell_of(a), cell_of(b)], idys) === nothing

            # No table yet => nothing to protect => never refuse.
            mktempdir() do fresh
                Test.@test preflight_threshold_table(
                    fresh, [cell_of(a)], idys) === nothing
            end
        end
    end
end

end  # module
