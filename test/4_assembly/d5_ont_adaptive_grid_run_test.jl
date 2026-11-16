# D5 grid RUNNER (td-yddpr.5): resume / replan, Phase B reuse, the run manifest,
# and the merge aggregations that feed the admissibility rule and tau. Uses
# synthetic checkpoints in a temp dir; no simulator, assembler or QUAST runs.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/d5_ont_adaptive_grid_run_test.jl")'

module D5OntAdaptiveGridRunTest

import Test
import DataFrames
import JSON

include(joinpath(@__DIR__, "..", "..", "benchmarking", "d5_ont_adaptive_grid_run.jl"))

function error_message(f)
    try
        f()
    catch e
        return sprint(showerror, e)
    end
    return nothing
end

const CTX = D5OntAdaptiveGridDriver.EMPTY_CONTEXT
function measured_metrics(nga50)
    merge(D5OntAdaptiveGridDriver.empty_metrics(),
        (; NGA50 = nga50, genome_fraction = 95.0))
end
as_dict(nt) = Dict{String, Any}(String(k) => v for (k, v) in pairs(nt))

function row_for(cell; phase = :A, status = "ok", nga50 = 1000.0,
        nga50_status = "measured", kwargs...)
    metrics = nga50 === missing ? D5OntAdaptiveGridDriver.empty_metrics() :
              measured_metrics(nga50)
    return D5OntAdaptiveGridDriver.cell_row(cell, CTX; phase, status, metrics,
        nga50_status, kwargs...)
end

function save_row(root, row)
    D5OntAdaptiveGridDriver.save_json_atomically(
        joinpath(root, "cells", row.cell_id, "cell_result.json"),
        D5OntAdaptiveGridDriver.row_to_json(row))
end

Test.@testset "D5 grid runner" begin
    cell = first(D5OntAdaptiveGridDriver.phase_a_cells())

    Test.@testset "run_or_resume: cache, retry, failure row, stale refusal, replan" begin
        mktempdir() do root
            calls = Ref(0)
            make = d -> (calls[] += 1; row_for(cell))
            Test.@test run_or_resume(make, root, cell, CTX; phase = :A) == "ok"
            Test.@test calls[] == 1
            # Completed cell is served from its checkpoint without recomputing.
            Test.@test run_or_resume(make, root, cell, CTX; phase = :A) == "ok"
            Test.@test calls[] == 1
            # An exception becomes an `error` row, which the next run retries.
            mktempdir() do root2
                boom = d -> error("simulated assembler crash")
                Test.@test run_or_resume(boom, root2, cell, CTX; phase = :A) == "error"
                saved = D5OntAdaptiveGridDriver.load_row(checkpoint_path(root2, cell))
                Test.@test occursin("simulated assembler crash", saved["error_message"])
                Test.@test run_or_resume(make, root2, cell, CTX; phase = :A) == "ok"
                Test.@test calls[] == 2
            end
            # A cached row judged stale is refused unless replan, then recomputed.
            stale = prior -> "computed under another tau."
            msg = error_message(() -> run_or_resume(make, root, cell, CTX; phase = :A,
                stale))
            Test.@test msg !== nothing && occursin("Pass --replan", msg)
            Test.@test calls[] == 2
            Test.@test run_or_resume(make, root, cell, CTX; phase = :A, stale,
                replan = true) == "ok"
            Test.@test calls[] == 3
            Test.@test needs_retry(["ok", "error", "quast_failed", "absent"]) == 2
        end
    end

    Test.@testset "Phase B checkpoints are tied to tau and tau_hat" begin
        adaptive = first(D5OntAdaptiveGridDriver.phase_b_cells())
        row = as_dict(row_for(adaptive; phase = :B, tau = 0.8, tau_hat = 0.6))
        cli = (; tau = 0.8, tau_hat = 0.6)
        Test.@test phase_b_stale(cli)(row) === nothing
        changed = phase_b_stale((; tau = 0.75, tau_hat = 0.6))(row)
        Test.@test changed !== nothing && occursin("tau=0.8", changed)
        Test.@test phase_b_stale((; tau = 0.8, tau_hat = 0.5))(row) !== nothing
        # A row with no recorded tau is never accepted as current.
        Test.@test phase_b_stale(cli)(as_dict(row_for(adaptive; phase = :B))) !== nothing
    end

    Test.@testset "reuse_row copies the comparator cell and keeps censoring" begin
        mktempdir() do root
            adaptive = first(filter(c -> c.technology === :ont,
                D5OntAdaptiveGridDriver.phase_b_cells()))
            source_cell = merge(adaptive, (; k_policy = 21))
            save_row(root,
                row_for(source_cell; nga50 = missing,
                    nga50_status = "censored_partial_alignment", n_contigs = 40))
            plan = (; k_star = 21,
                reuse_cell_id = D5OntAdaptiveGridDriver.cell_id(source_cell))
            reused = reuse_row(root, adaptive, plan, CTX, (; tau = 0.8, tau_hat = 0.6))
            Test.@test reused.status == "ok" && reused.phase == "B"
            Test.@test ismissing(reused.NGA50)
            Test.@test reused.nga50_status == "censored_partial_alignment"
            Test.@test reused.reused_from_cell_id == plan.reuse_cell_id
            Test.@test reused.k_selected == 21 && reused.tau == 0.8
            Test.@test reused.policy_level == "adaptive"
            save_row(root, row_for(source_cell; status = "error", nga50 = missing,
                nga50_status = missing))
            msg = error_message(() -> reuse_row(root, adaptive, plan, CTX, (;)))
            Test.@test msg !== nothing && occursin("not reusable", msg)
            missing_plan = (; k_star = 13, reuse_cell_id = "no_such_cell")
            msg = error_message(() -> reuse_row(root, adaptive, missing_plan, CTX, (;)))
            Test.@test msg !== nothing && occursin("has no checkpoint", msg)
        end
    end

    Test.@testset "run manifest: one SHA, clean tree, same pins, toy kept apart" begin
        mktempdir() do root
            git = (; sha = "a" ^ 40, dirty = false)
            check_run_manifest!(root, git; allow_mismatch = false, toy = true)
            Test.@test isfile(joinpath(root, "run_manifest.json"))
            Test.@test check_run_manifest!(root, git; allow_mismatch = false,
                toy = true)["mycelia_git_sha"] == git.sha
            msg = error_message(() -> check_run_manifest!(
                root, (; sha = "b" ^ 40,
                    dirty = false);
                allow_mismatch = false, toy = true))
            Test.@test msg !== nothing && occursin("records Mycelia", msg)
            check_run_manifest!(root, (; sha = "b" ^ 40, dirty = false);
                allow_mismatch = true, toy = true)  # override warns, does not throw
            msg = error_message(() -> check_run_manifest!(root, git;
                allow_mismatch = false, toy = false))
            Test.@test msg !== nothing && occursin("refusing to mix toy and full", msg)
            msg = error_message(() -> check_run_manifest!(
                root, (; sha = "unknown",
                    dirty = false);
                allow_mismatch = false, toy = true))
            Test.@test msg !== nothing && occursin("cannot read the Mycelia git SHA", msg)
            msg = error_message(() -> check_run_manifest!(
                root, (; sha = git.sha,
                    dirty = true); allow_mismatch = false, toy = true))
            Test.@test msg !== nothing && occursin("uncommitted changes", msg)
            # A pin edited under an existing grid is refused.
            manifest = JSON.parsefile(joinpath(root, "run_manifest.json"))
            manifest["badread_error_model"] = "nanopore2020"
            D5OntAdaptiveGridDriver.save_json_atomically(
                joinpath(root, "run_manifest.json"), manifest)
            msg = error_message(() -> check_run_manifest!(root, git;
                allow_mismatch = false, toy = true))
            Test.@test msg !== nothing && occursin("badread_error_model", msg)
        end
    end

    Test.@testset "replicate evaluability counts failed, absent and adaptive seeds" begin
        seeds = collect(D5OntAdaptiveGridDriver.D5OntAdaptiveGridAxes.D5_SEEDS)
        Test.@test length(seeds) == 3
        same_group(c) = (c.technology, c.identity, c.coverage, c.corrector, c.k_policy,
            c.organism)
        trio = filter(c -> same_group(c) == same_group(cell),
            D5OntAdaptiveGridDriver.phase_a_cells())
        Test.@test length(trio) == 3
        cols = D5OntAdaptiveGridDriver.cell_row_columns()
        # Seed 3 crashed: the group must NOT read as fully measured.
        rows = [as_dict(row_for(trio[1])), as_dict(row_for(trio[2])),
            as_dict(row_for(trio[3]; status = "error", nga50 = missing,
                nga50_status = missing))]
        ev = replicate_evaluability(frame(rows, cols))
        Test.@test DataFrames.nrow(ev) == 1
        Test.@test ev.n_rows[1] == 3 && ev.n_measured[1] == 2 && ev.n_error[1] == 1
        Test.@test ev.all_replicates_measured[1] == false
        # Seed 3 never ran: an absent row still counts against the group.
        rows[3] = as_dict(D5OntAdaptiveGridDriver.absent_row(trio[3]; phase = :A))
        ev = replicate_evaluability(frame(rows, cols))
        Test.@test ev.n_absent[1] == 1 && ev.all_replicates_measured[1] == false
        # All three measured: true, with a CV.
        rows[3] = as_dict(row_for(trio[3]; nga50 = 1200.0))
        ev = replicate_evaluability(frame(rows, cols))
        Test.@test ev.all_replicates_measured[1] == true && !ismissing(ev.cv_nga50[1])
        # Adaptive seeds that chose different k* stay ONE replicate group.
        adaptive = filter(
            c -> same_group(merge(c, (; k_policy = cell.k_policy))) ==
                 same_group(cell),
            D5OntAdaptiveGridDriver.phase_b_cells())
        Test.@test length(adaptive) == 3
        brows = [as_dict(row_for(c; phase = :B, k_selected = k))
                 for (c, k) in zip(adaptive, (17, 21, 23))]
        ev = replicate_evaluability(frame(brows, cols))
        Test.@test DataFrames.nrow(ev) == 1 && ev.n_rows[1] == 3
        Test.@test ev.policy_level[1] == "adaptive"
    end

    Test.@testset "Illumina optimum: only fully measured rungs can win" begin
        calib = filter(
            c -> c.technology === :illumina && c.corrector === :none &&
                 c.organism === :Lambda && c.coverage == 30,
            D5OntAdaptiveGridDriver.phase_a_cells())
        cols = D5OntAdaptiveGridDriver.cell_row_columns()
        rows = Dict{String, Any}[]
        for c in calib
            if c.k_policy == 21
                push!(rows, as_dict(row_for(c; nga50 = 9000.0)))  # 3/3 measured
            elseif c.k_policy == 31
                # Higher NGA50 on the one seed that survived: must not win.
                push!(rows,
                    as_dict(c.seed ==
                            first(D5OntAdaptiveGridDriver.D5OntAdaptiveGridAxes.D5_SEEDS) ?
                            row_for(c; nga50 = 20000.0) :
                            row_for(c; nga50 = missing, nga50_status = "no_contigs_ge_min")))
            end
        end
        table = illumina_optimum(frame(rows, cols))
        per_cov = table[table.scope .== "per_coverage", :]
        best = per_cov[per_cov.is_optimum, :]
        Test.@test DataFrames.nrow(best) == 1 && best.k[1] == 21
        k31 = per_cov[per_cov.k .== 31, :]
        Test.@test k31.n_measured[1] == 1 && k31.all_measured[1] == false
        Test.@test k31.n_expected[1] == 3
    end

    Test.@testset "expected cells for the toy and the full grid" begin
        Test.@test length(expected_cells(true, false)) == 7
        Test.@test length(expected_cells(false, false)) == 1008 + 96 + 336
        Test.@test length(expected_cells(false, true)) == 1008 + 96 + 336 + 144 + 48
    end
end

end # module
