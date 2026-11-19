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

    toy_sets = D5OntAdaptiveGridDriver.toy_read_sets()
    toy_cells = reduce(vcat, D5OntAdaptiveGridDriver.toy_cells(rs) for rs in toy_sets)

    Test.@testset "merge: absent and unreadable cells, end to end through run_merge" begin
        mktempdir() do root
            save_row(root, row_for(toy_cells[1]))
            save_row(root, row_for(toy_cells[2]))
            bad = joinpath(root, "cells", D5OntAdaptiveGridDriver.cell_id(toy_cells[3]))
            mkpath(bad)
            write(joinpath(bad, "cell_result.json"), "{not json")
            collected = collect_merge_rows(root, true)
            Test.@test length(collected.rows) == 7
            Test.@test collected.n_absent == 4 && collected.unreadable == [basename(bad)]
            Test.@test collected.has_b == false
            p = merge_problems(collected.rows;
                n_unreadable_entries = length(collected.unreadable))
            Test.@test p.n_absent == 4 && p.n_unreadable == 1 && p.n_incomplete == 5
            # Through the real entry point: incomplete -> non-zero, tables still written.
            # (Rows built from EMPTY_CONTEXT have git_dirty = missing, i.e. not provably
            # clean; --allow-mixed-code waives exactly that. --allow-sha-mismatch is for the manifest.)
            cli = (; toy = true, allow_sha_mismatch = true, allow_mixed_code = true,
                allow_incomplete = false)
            Test.@test run_merge(cli, root) == 5
            tsv = read(joinpath(root, "d5_cells.tsv"), String)
            Test.@test count("\tabsent\t", tsv) == 4 && count("\tunreadable\t", tsv) == 1
            prov = JSON.parsefile(joinpath(root, "merge_provenance.json"))
            Test.@test prov["n_absent"] == 4 && prov["blocking_problems"] == 5
            Test.@test prov["overrides"]["allow_incomplete"] == false
            Test.@test length(prov["blocking_reasons"]) == 1
            Test.@test length(prov["waived_reasons"]) == 1  # unclean rows, waived
            # Without --allow-mixed-code the unclean (git_dirty = missing) rows block too,
            # and --allow-incomplete does not waive them.
            Test.@test run_merge(
                merge(cli, (; allow_incomplete = true,
                    allow_mixed_code = false)), root) == 2
            Test.@test run_merge(merge(cli, (; allow_incomplete = true)), root) == 0
            # A stray unreadable directory that is not a registered cell still counts.
            stray = joinpath(root, "cells", "not_a_registered_cell")
            mkpath(stray)
            write(joinpath(stray, "cell_result.json"), "{")
            again = collect_merge_rows(root, true)
            Test.@test length(again.unreadable) == 2
            Test.@test merge_problems(again.rows;
                n_unreadable_entries = length(again.unreadable)).n_unreadable == 2
            # One Phase B checkpoint makes the toy's 5 Phase B cells expected.
            pb = first(filter(c -> c.k_policy === :adaptive,
                reduce(vcat, D5OntAdaptiveGridDriver.cells_for_read_set(rs; phase = :B)
                for rs in toy_sets)))
            save_row(root, row_for(pb; phase = :B, tau = 0.8, tau_hat = 0.6))
            withb = collect_merge_rows(root, true)
            Test.@test withb.has_b && length(withb.rows) == 12
            Test.@test count(
                r -> isequal(r["phase"], "B") &&
                     isequal(r["status"], "absent"), withb.rows) == 4
        end
    end

    Test.@testset "merge verdict: completeness, code version and tau are waived differently" begin
        b = filter(c -> c.k_policy === :adaptive, D5OntAdaptiveGridDriver.phase_b_cells())
        clean(r) = merge(r, Dict{String, Any}("git_dirty" => false,
            "mycelia_git_sha" => "a" ^ 40))
        rows = [clean(as_dict(row_for(b[1]; phase = :B, tau = 0.8, tau_hat = 0.6))),
            clean(as_dict(row_for(b[2]; phase = :B, tau = 0.8, tau_hat = 0.6)))]
        p = merge_problems(rows)
        Test.@test p.n_incomplete == 0 && p.n_tau_conflicts == 0
        Test.@test p.n_sha_conflicts == 0 && p.n_dirty_rows == 0
        none = (; allow_incomplete = false, allow_mixed_code = false)
        Test.@test merge_verdict(p; none...).blocking == 0
        # A retryable row is incomplete; --allow-incomplete waives it (and says so).
        err = clean(as_dict(row_for(b[3]; phase = :B, status = "error", nga50 = missing,
            nga50_status = missing, tau = 0.8, tau_hat = 0.6)))
        pe = merge_problems(vcat(rows, [err]))
        Test.@test pe.n_retry == 1 && pe.n_incomplete == 1
        Test.@test merge_verdict(pe; none...).blocking == 1
        ve = merge_verdict(pe; allow_incomplete = true, allow_mixed_code = false)
        Test.@test isempty(ve.blocking_reasons) && length(ve.waived_reasons) == 1
        Test.@test ve.blocking == 0 && ve.waived == 1
        # Mixed tau is NEVER waivable, not even with both flags.
        mixed = vcat(rows,
            [clean(as_dict(row_for(b[3]; phase = :B, tau = 0.75, tau_hat = 0.6)))])
        pm = merge_problems(mixed)
        Test.@test pm.n_tau_conflicts == 1 && length(pm.tau_pairs) == 2
        Test.@test merge_verdict(pm; allow_incomplete = true,
            allow_mixed_code = true).blocking == 1
        # Mixed SHAs / unclean rows block unless --allow-mixed-code, never by
        # --allow-incomplete.
        shas = [merge(rows[1], Dict{String, Any}("mycelia_git_sha" => "b" ^ 40)), rows[2]]
        ps = merge_problems(shas)
        Test.@test ps.n_sha_conflicts == 1
        Test.@test merge_verdict(ps; allow_incomplete = true,
            allow_mixed_code = false).blocking == 1
        Test.@test merge_verdict(ps; allow_incomplete = false,
            allow_mixed_code = true).blocking == 0
        dirty = [merge(rows[1], Dict{String, Any}("git_dirty" => true)), rows[2]]
        Test.@test merge_problems(dirty).n_dirty_rows == 1
        # Absent rows carry no SHA / dirty state and are not "mixed code".
        absent = as_dict(D5OntAdaptiveGridDriver.absent_row(b[4]; phase = :B))
        pa = merge_problems(vcat(rows, [absent]))
        Test.@test pa.n_absent == 1 && pa.n_dirty_rows == 0 && pa.n_sha_conflicts == 0
        # Corrected-survival diagnostic outcomes are reported, never problems.
        diag = [merge(r, Dict{String, Any}("corrected_survival_status" => s))
                for (r, s) in
                    zip(vcat(rows, rows[1:1]), ("failed", "unavailable", "measured"))]
        pd = merge_problems(diag)
        Test.@test pd.n_corrected_survival_failed == 1
        Test.@test pd.n_corrected_survival_unavailable == 1
        Test.@test merge_verdict(pd; none...).blocking == 0
    end

    Test.@testset "run_or_resume carries the Phase B plan onto an error row" begin
        mktempdir() do root
            boom = d -> error("assembler crashed")
            adaptive = first(D5OntAdaptiveGridDriver.phase_b_cells())
            run_or_resume(boom, root, adaptive, CTX; phase = :B,
                error_kwargs = (; tau = 0.8, tau_hat = 0.6, k_selected = "NoFeasibleK"))
            saved = D5OntAdaptiveGridDriver.load_row(checkpoint_path(root, adaptive))
            Test.@test saved["status"] == "error" && saved["tau"] == 0.8
            Test.@test saved["tau_hat"] == 0.6 && saved["k_selected"] == "NoFeasibleK"
        end
    end

    Test.@testset "Phase B wiring: plan metadata and the toy reuse rule" begin
        adaptive = first(filter(c -> c.technology === :ont,
            D5OntAdaptiveGridDriver.phase_b_cells()))
        plan = (; action = :reuse, k_star = 21, reason = "admissible", floor_k = 12,
            k_star_hat = D5OntAdaptiveGridDriver.NO_FEASIBLE_K,
            reason_hat = "no_rung_meets_threshold", reuse_cell_id = "x")
        kw = phase_b_kwargs(plan, (; tau = 0.8, tau_hat = 0.6))
        Test.@test kw.tau == 0.8 && kw.tau_hat == 0.6 && kw.k_floor == 12
        Test.@test kw.k_star_hat == "NoFeasibleK" && kw.k_star_reason == "admissible"
        # Toy: a source outside the toy's Phase A cells is assembled at the same k*;
        # a source INSIDE the toy keeps :reuse (so a missing checkpoint fails loudly);
        # the full grid always keeps :reuse.
        rs95 = first(filter(rs -> rs.technology === :ont && rs.identity == 95, toy_sets))
        in_toy = first(D5OntAdaptiveGridDriver.toy_cells(rs95))
        inside = merge(plan, (; reuse_cell_id = D5OntAdaptiveGridDriver.cell_id(in_toy)))
        Test.@test phase_b_action(inside, rs95, true) === :reuse
        Test.@test phase_b_action(plan, rs95, true) === :assemble
        Test.@test phase_b_action(plan, rs95, false) === :reuse
        Test.@test phase_b_action(merge(plan, (; action = :infeasible)), rs95, true) ===
                   :infeasible
    end

    Test.@testset "reference length and read count guards" begin
        info = D5OntAdaptiveGridDriver.D5_ORGANISM_REFERENCES[:Lambda]
        Test.@test info.expected_bp == 48_502
        Test.@test D5OntAdaptiveGridDriver.D5_ORGANISM_REFERENCES[:T4].expected_bp ==
                   168_903
        Test.@test check_reference_length(48_502, info, "/r.fna") === nothing
        msg = error_message(() -> check_reference_length(48_501, info, "/r.fna"))
        Test.@test msg !== nothing && occursin("the RefSeq record is 48502 bp", msg)
        meta = Dict{String, Any}("read_set_id" => "rs", "n_reads" => 3)
        Test.@test length(check_read_count([1, 2, 3], meta, "/d")) == 3
        msg = error_message(() -> check_read_count([1, 2], meta, "/d"))
        Test.@test msg !== nothing && occursin("has 2 reads on disk", msg)
    end

    Test.@testset "manifest: pins are not waivable; unknown dirty state refuses" begin
        mktempdir() do root
            git = (; sha = "c" ^ 40, dirty = false)
            check_run_manifest!(root, git; allow_mismatch = false, toy = true)
            msg = error_message(() -> check_run_manifest!(
                root, (; sha = git.sha,
                    dirty = missing);
                allow_mismatch = false, toy = true))
            Test.@test msg !== nothing && occursin("could not be read", msg)
            path = joinpath(root, "run_manifest.json")
            manifest = JSON.parsefile(path)
            delete!(manifest, "badread_length_dist")  # a legacy manifest
            D5OntAdaptiveGridDriver.save_json_atomically(path, manifest)
            msg = error_message(() -> check_run_manifest!(root, git;
                allow_mismatch = true, toy = true))
            Test.@test msg !== nothing && occursin("badread_length_dist", msg)
            Test.@test occursin("= nothing", msg)
            manifest = JSON.parsefile(path)
            manifest["badread_length_dist"] = D5OntAdaptiveGridDriver.D5_BADREAD_LENGTH_DIST
            manifest["badread_qscore_model"] = "nanopore2020"
            D5OntAdaptiveGridDriver.save_json_atomically(path, manifest)
            msg = error_message(() -> check_run_manifest!(root, git;
                allow_mismatch = true, toy = true))
            Test.@test msg !== nothing && occursin("Use a fresh output root", msg)
            Test.@test occursin("badread_qscore_model", msg)
        end
    end

    Test.@testset "CLI: new switches, and unknown arguments are refused" begin
        cli = D5OntAdaptiveGridDriver.parse_cli(["--phase", "B", "--tau", "0.8",
            "--tau-hat", "0.6", "--replan"])
        Test.@test cli.replan && !cli.allow_incomplete
        Test.@test D5OntAdaptiveGridDriver.parse_cli(["--phase", "merge",
            "--allow-incomplete"]).allow_incomplete
        msg = error_message(() -> D5OntAdaptiveGridDriver.parse_cli(["--phase", "B",
            "--tau", "0.8", "--tau-hat", "0.6", "--replann"]))
        Test.@test msg !== nothing && occursin("unknown argument", msg)
        every = D5OntAdaptiveGridDriver.parse_cli(["--phase", "A", "--read-set-index",
            "0", "--output-dir", "/x", "--toy", "--allow-sha-mismatch", "--replan",
            "--allow-incomplete", "--allow-mixed-code", "--tmp-dir", "/t"])
        Test.@test every.phase === :A && every.read_set_index == 0
        Test.@test every.output_dir == "/x" && every.toy && every.allow_sha_mismatch
        Test.@test every.replan && every.allow_incomplete && every.allow_mixed_code
        Test.@test every.tmp_dir == "/t"
        Test.@test D5OntAdaptiveGridDriver.parse_cli(["--phase", "A"]).tmp_dir === nothing
        msg = error_message(() -> D5OntAdaptiveGridDriver.parse_cli(["--phase", "A",
            "--tmp-dir", " "]))
        Test.@test msg !== nothing && occursin("--tmp-dir must be a non-empty path", msg)
        msg = error_message(() -> D5OntAdaptiveGridDriver.parse_cli(["--phase", "A",
            "--tmp-dir", "/a", "--tmp-dir", "/b"]))
        Test.@test msg !== nothing && occursin("--tmp-dir given more than once", msg)
        msg = error_message(() -> D5OntAdaptiveGridDriver.parse_cli(["--phase", "merge",
            "--output-dir", ""]))
        Test.@test msg !== nothing && occursin("non-empty path", msg)
        both = D5OntAdaptiveGridDriver.parse_cli(["--phase", "B", "--tau", "0.8",
            "--tau-hat", "0.6", "--list"])
        Test.@test both.tau == 0.8 && both.tau_hat == 0.6
    end

    Test.@testset "temp root: resolved, pinned to TMPDIR, refused when rename fails" begin
        withenv("TMPDIR" => get(ENV, "TMPDIR", nothing)) do
            root = realpath(mktempdir())
            tmp = prepare_temp_root(nothing, root; rename_errno = _ -> 0)
            Test.@test tmp == joinpath(root, "tmp") && isdir(tmp)
            Test.@test ENV["TMPDIR"] == tmp
            # A requested root reached through a symlink is probed and pinned by its
            # real path.
            mkdir(joinpath(root, "real"))
            symlink(joinpath(root, "real"), joinpath(root, "link"))
            seen = Ref("")
            chosen = prepare_temp_root(joinpath(root, "link"), root;
                rename_errno = d -> (seen[] = d; 0))
            Test.@test chosen == joinpath(root, "real") && ENV["TMPDIR"] == chosen
            Test.@test seen[] == chosen
            # EINVAL (Lustre) refuses before TMPDIR moves, naming the errno.
            ENV["TMPDIR"] = "unchanged"
            msg = error_message(() -> prepare_temp_root(nothing, root;
                rename_errno = _ -> 22))
            Test.@test msg !== nothing && occursin("RENAME_NOREPLACE", msg) &&
                       occursin("--tmp-dir", msg) &&
                       occursin(Base.Libc.strerror(22), msg)
            Test.@test ENV["TMPDIR"] == "unchanged"
            msg = error_message(() -> prepare_temp_root(nothing, root;
                rename_errno = _ -> -1))
            Test.@test msg !== nothing && occursin("errno unreadable", msg)
            msg = error_message(() -> prepare_temp_root("relative/tmp", root;
                rename_errno = _ -> 0))
            Test.@test msg !== nothing && occursin("must be an absolute path", msg)
            Test.@test ENV["TMPDIR"] == "unchanged"
        end
    end

    Test.@testset "temp root: phases A and B get one, merge leaves TMPDIR alone" begin
        root = realpath(mktempdir())  # before TMPDIR points nowhere
        withenv("TMPDIR" => "unchanged") do
            merge_cli = D5OntAdaptiveGridDriver.parse_cli(["--phase", "merge"])
            Test.@test temp_root_for(merge_cli, root; rename_errno = _ -> 22) === nothing
            Test.@test ENV["TMPDIR"] == "unchanged" && !ispath(joinpath(root, "tmp"))
            for args in (["--phase", "A"], ["--phase", "B", "--tau", "0.5", "--tau-hat",
                "0.5"])
                cli = D5OntAdaptiveGridDriver.parse_cli(vcat(args,
                    ["--tmp-dir", joinpath(root, "node_local")]))
                ENV["TMPDIR"] = "unchanged"
                got = temp_root_for(cli, root; rename_errno = _ -> 0)
                Test.@test got == joinpath(root, "node_local") && ENV["TMPDIR"] == got
                # A and B are refused on a rejecting filesystem (the guard is not
                # skipped), and the refusal leaves TMPDIR alone.
                ENV["TMPDIR"] = "unchanged"
                msg = error_message(() -> temp_root_for(cli, root; rename_errno = _ -> 22))
                Test.@test msg !== nothing && occursin("RENAME_NOREPLACE", msg)
                Test.@test ENV["TMPDIR"] == "unchanged"
            end
        end
    end

    Test.@testset "temp root probe: real wrapper, errno and cleanup (Linux)" begin
        dir = realpath(mktempdir())
        if Sys.islinux()
            # The real wrapper with RENAME_NOREPLACE onto an existing directory gives
            # EEXIST (17): pins the flag value and the errno read through the wrapper.
            # Without the flag, renaming onto an empty directory would succeed.
            mkdir(joinpath(dir, "a"))
            mkdir(joinpath(dir, "b"))
            rc = Mycelia._metamdbg_linux_renameat2(Cint(-100), joinpath(dir, "a"),
                Cint(-100), joinpath(dir, "b"), Cuint(1))
            errno = Base.Libc.errno()
            Test.@test rc != 0 && errno == 17
            rm(joinpath(dir, "a"))
            rm(joinpath(dir, "b"))
            # The probe on this machine's temp filesystem succeeds and cleans up.
            Test.@test noreplace_rename_errno(dir) == 0
            Test.@test isempty(readdir(dir))
            # The probe passes AT_FDCWD and RENAME_NOREPLACE (a flag of 0 would let a
            # Lustre rename succeed here while the corrector still fails).
            calls = []
            Test.@test noreplace_rename_errno(dir;
                rename = (
                    sfd, s, dfd, d, flags) -> (push!(calls, (sfd, dfd, flags));
                    Cint(0))) == 0
            Test.@test calls == [(Cint(-100), Cint(-100), Cuint(1))]
            # A failure reports its errno; one whose errno reads 0 still fails closed.
            Test.@test noreplace_rename_errno(dir;
                rename = (_...) -> (Base.Libc.errno(22); Cint(-1))) == 22
            Test.@test noreplace_rename_errno(dir;
                rename = (_...) -> (Base.Libc.errno(0); Cint(-1))) == -1
            # A throwing rename still removes the probe directory.
            msg = error_message(() -> noreplace_rename_errno(dir;
                rename = (_...) -> error("probe boom")))
            Test.@test msg !== nothing && occursin("probe boom", msg)
            Test.@test isempty(readdir(dir))
        else
            # Elsewhere nothing is probed (macOS uses renameatx_np).
            Test.@test noreplace_rename_errno(dir;
                rename = (_...) -> error("not called")) == 0
        end
    end

    Test.@testset "sbatch temp-root cleanup: evidence kept on failure and on signals" begin
        helper = joinpath(@__DIR__, "..", "..", "benchmarking", "d5_temp_root.sh")
        Test.@test isfile(helper)
        if Sys.islinux() && Sys.which("bash") !== nothing
            # Runs the helper under bash: `script` gets $1 = helper, $2 = temp root,
            # $3 = evidence dir, $4 = exit code.
            # With via_symlink the temp root is passed through a symlinked parent
            # (a D5_TMPDIR that is a symlink); real_root is where it actually is.
            function run_helper(script; rc = 0, populate = true, keep_exists = false,
                    via_symlink = false)
                base = realpath(mktempdir())
                real_root = joinpath(base, "real", "root")
                mkpath(real_root)
                root = real_root
                if via_symlink
                    symlink(joinpath(base, "real"), joinpath(base, "link"))
                    root = joinpath(base, "link", "root")
                end
                populate && write(joinpath(real_root, "evidence.txt"), "kept")
                keep = joinpath(base, "tmp_evidence", "job1")
                keep_exists && mkpath(keep)
                err = IOBuffer()
                # No real cgroup is read here: the memory branches have their own
                # test on a fake tree below.
                cmd = addenv(
                    Cmd(["bash", "-c", script, "_", helper, root, keep, string(rc)]),
                    "D5_PROC_CGROUP" => "/nonexistent")
                proc = run(pipeline(ignorestatus(cmd); stderr = err))
                return (; root, real_root, keep, code = proc.exitcode,
                    log = String(take!(err)))
            end
            cleanup = raw"""source "$1" && d5_temp_root_cleanup "$4" "$2" "$3" """
            # Non-zero exit: the temp root is copied whole, then removed. Both the
            # path as given and its resolved form are recorded, so error rows that
            # cite resolved TMPDIR paths map into the evidence dir.
            r = run_helper(cleanup; rc = 3, via_symlink = true)
            Test.@test r.root != r.real_root
            Test.@test read(joinpath(r.keep, "evidence.txt"), String) == "kept"
            Test.@test strip(read(joinpath(r.keep, "ORIGINAL_TEMP_ROOT"), String)) == r.root
            Test.@test strip(read(joinpath(r.keep, "RESOLVED_TEMP_ROOT"), String)) ==
                       r.real_root
            Test.@test !ispath(r.real_root) && !ispath(r.keep * ".partial")
            Test.@test occursin("bytes (end state, not the peak)", r.log)
            r = run_helper(cleanup; rc = 3)
            Test.@test !ispath(r.root) && !ispath(r.keep * ".partial")
            Test.@test occursin("preserved at", r.log)
            Test.@test occursin("cgroup memory.peak unavailable (no cgroup v2 directory", r.log)
            # The copy is announced before it starts, so a copy cut off by SIGKILL
            # shows as a "copying" line with no "preserved" line after it.
            Test.@test occursin("copying $(r.root) -> $(r.keep).partial", r.log)
            # Clean exit, or nothing to keep: removed, nothing copied.
            r = run_helper(cleanup; rc = 0)
            Test.@test !ispath(r.root) && !ispath(r.keep)
            r = run_helper(cleanup; rc = 3, populate = false)
            Test.@test !ispath(r.root) && !ispath(r.keep)
            # An existing evidence dir is never overwritten; the WARN says so.
            r = run_helper(cleanup; rc = 3, keep_exists = true)
            Test.@test isempty(readdir(r.keep)) && !ispath(r.root)
            Test.@test occursin("WARN", r.log) && occursin("already exists", r.log)
            # SIGTERM (time limit, scancel) while the job runs: the traps turn it into
            # exit 143, so the evidence is kept rather than deleted with status 0.
            r = run_helper(raw"""source "$1" && d5_install_temp_root_traps "$2" "$3" &&
                true && kill -TERM $$ && sleep 30""")
            Test.@test r.code == 143
            Test.@test read(joinpath(r.keep, "evidence.txt"), String) == "kept"
            Test.@test !ispath(r.root) && occursin("rc=143", r.log)
            # A normal exit code passes through the EXIT trap unchanged.
            r = run_helper(raw"""source "$1" && d5_install_temp_root_traps "$2" "$3" &&
                exit 3""")
            Test.@test r.code == 3 && isdir(r.keep) && !ispath(r.root)
            # A second SIGTERM during the copy (a repeated scancel) is ignored: the
            # copy completes and is renamed, and the temp root is still removed.
            # `cp` is shadowed by a function that signals the shell first.
            r = run_helper(raw"""source "$1" && d5_install_temp_root_traps "$2" "$3" &&
                cp() { kill -TERM $$; command cp "$@"; } && exit 3""")
            Test.@test r.code == 3
            Test.@test read(joinpath(r.keep, "evidence.txt"), String) == "kept"
            Test.@test !ispath(r.keep * ".partial") && !ispath(r.root)
            # A failed copy (read-only evidence parent) warns, still removes the
            # temp root, and leaves the job's exit code unchanged. Root ignores
            # write permission, so this case is skipped when run as root.
            if ccall(:geteuid, Cuint, ()) == 0
                Test.@test_skip "read-only evidence parent is writable by root"
            else
                r = run_helper(raw"""source "$1" &&
                    d5_install_temp_root_traps "$2" "$3" && mkdir -p "$(dirname "$3")" &&
                    chmod a-w "$(dirname "$3")" && exit 3""")
                chmod(dirname(r.keep), 0o755)
                Test.@test r.code == 3 && !ispath(r.root) && !ispath(r.keep)
                Test.@test occursin("cannot create $(r.keep).partial", r.log)
                # A copy that fails part-way (an unreadable file) warns, leaves the
                # .partial copy for inspection, and still removes the temp root.
                r = run_helper(raw"""chmod 000 "$2/evidence.txt" && source "$1" &&
                    d5_install_temp_root_traps "$2" "$3" && exit 3""")
                Test.@test r.code == 3 && !ispath(r.root) && !ispath(r.keep)
                Test.@test isdir(r.keep * ".partial")
                Test.@test occursin("could not preserve", r.log)
                # An unreadable subdirectory: du reports a lower bound, and rm cannot
                # remove it, which is reported rather than passed over.
                r = run_helper(raw"""mkdir "$2/sub" && touch "$2/sub/f" &&
                    chmod 000 "$2/sub" && source "$1" &&
                    d5_temp_root_cleanup 0 "$2" "$3" """)
                Test.@test occursin(r"at least [0-9]+ bytes \(du incomplete", r.log)
                Test.@test occursin("not fully removed", r.log)
                chmod(joinpath(r.root, "sub"), 0o755)
            end
            # A temp root that is already gone: du fails and the size is reported
            # as unknown rather than as du's error text followed by "bytes".
            r = run_helper(raw"""rm -rf -- "$2" &&
                source "$1" && d5_temp_root_cleanup 0 "$2" "$3" """)
            Test.@test occursin("size unknown, du failed", r.log)
            Test.@test !occursin(r"cannot access.*bytes", r.log)
            # Memory is read at the first cgroup level with a numeric memory.max, but
            # never above job_<id>: higher levels hold every job on the node. A fake
            # tree with a numeric limit above the job checks the bound.
            cg = realpath(mktempdir())
            scope = joinpath(cg, "slurmstepd.scope")
            job = joinpath(scope, "job_1")
            leaf = joinpath(job, "step_batch", "user", "task_0")
            mkpath(leaf)
            write(joinpath(cg, "memory.max"), "99999\n")
            write(joinpath(scope, "memory.max"), "99999\n")
            write(joinpath(scope, "memory.peak"), "99999\n")
            write(joinpath(leaf, "memory.max"), "max\n")
            write(joinpath(leaf, "memory.peak"), "5\n")
            write(joinpath(job, "memory.max"), "1000\n")
            write(joinpath(job, "memory.peak"), "900\n")
            write(joinpath(job, "memory.current"), "800\n")
            write(joinpath(job, "memory.events"), "low 0\nhigh 0\nmax 3\noom 2\noom_kill 1\n")
            in_job = joinpath(cg, "proc_self_cgroup")
            write(in_job, "0::/slurmstepd.scope/job_1/step_batch/user/task_0\n")
            function memory_log(proc_cgroup = in_job)
                err = IOBuffer()
                script = raw"""source "$1" &&
                    D5_CGROUP_ROOT="$2" D5_PROC_CGROUP="$3" d5_log_memory_peak"""
                run(pipeline(Cmd(["bash", "-c", script, "_", helper, cg, proc_cgroup]);
                    stderr = err))
                return String(take!(err))
            end
            mlog = memory_log()
            Test.@test occursin("limiting level $(job):", mlog)
            Test.@test occursin("memory.max 1000 bytes", mlog)
            Test.@test occursin("memory.peak (RSS + tmpfs + page cache) 900 bytes", mlog)
            Test.@test occursin("oom_kill 1", mlog)
            # No memory.peak at the limiting level (kernel < 5.19): memory.current,
            # labelled as not a peak, and the by-hand advice.
            mv(joinpath(job, "memory.peak"), joinpath(job, "memory.peak.off"))
            mlog = memory_log()
            Test.@test occursin("memory.peak unreadable (memory.current now 800", mlog)
            Test.@test occursin("by hand", mlog)
            mv(joinpath(job, "memory.peak.off"), joinpath(job, "memory.peak"))
            # No numeric limit up to the job level: the walk stops there and does not
            # report the node-wide level above it.
            write(joinpath(job, "memory.max"), "max\n")
            mlog = memory_log()
            Test.@test occursin("no numeric memory.max", mlog)
            Test.@test occursin("leaf memory.peak: 5 bytes;", mlog)
            Test.@test !occursin("99999", mlog)
            # Not inside a job cgroup (a hand run): nothing is walked.
            mkpath(joinpath(cg, "user.slice", "session"))
            outside = joinpath(cg, "proc_outside")
            write(outside, "0::/user.slice/session\n")
            mlog = memory_log(outside)
            Test.@test occursin("not inside a SLURM job_<id> cgroup", mlog)
            Test.@test !occursin("99999", mlog)
            # An unreadable memory.max stops the walk instead of skipping the level.
            if ccall(:geteuid, Cuint, ()) != 0
                write(joinpath(job, "memory.max"), "1000\n")
                chmod(joinpath(job, "memory.max"), 0o000)
                mlog = memory_log()
                chmod(joinpath(job, "memory.max"), 0o644)
                Test.@test occursin("memory.max at $(job) unreadable", mlog)
                Test.@test !occursin("99999", mlog)
            end
        else
            Test.@test_skip "temp-root cleanup runs under bash with GNU coreutils on Linux"
        end
        # A requeued job keeps its SLURM_JOB_ID; the restart count keeps each
        # attempt's evidence in its own directory. A text check on the shared task
        # body only: it does not run a requeue.
        task_body = joinpath(@__DIR__, "..", "..", "benchmarking", "d5_grid_task.sh")
        Test.@test occursin(
            raw"tmp_evidence/${SLURM_JOB_ID:-$$}.r${SLURM_RESTART_COUNT:-0}",
            read(task_body, String))
    end

    Test.@testset "cell task: needs a prepared read set, never simulates" begin
        # A read set as its read-set task leaves it: readset.json written last,
        # reads on disk, and the reference staged where prepare_reference looks
        # (download_genome_by_accession reuses an existing <refs>/<accession>.fna).
        function plant_read_set(root, rs)
            refs = mkpath(joinpath(root, "refs"))
            info = D5OntAdaptiveGridDriver.D5_ORGANISM_REFERENCES[rs.organism]
            write(joinpath(refs, info.accession * ".fna"),
                ">$(info.accession)\n" * "ACGT"^(info.expected_bp ÷ 4) *
                "A"^(info.expected_bp % 4) * "\n")
            rs_dir = mkpath(joinpath(root, "readsets",
                D5OntAdaptiveGridDriver.read_set_id(rs)))
            reads = joinpath(rs_dir, "reads.fq")
            write(reads, "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTTGGCCAA\n+\nIIIIIIII\n")
            D5OntAdaptiveGridDriver.save_json_atomically(joinpath(rs_dir, "readset.json"),
                Dict("read_set_id" => D5OntAdaptiveGridDriver.read_set_id(rs),
                    "reads_paths" => [reads], "n_reads" => 2, "read_bases" => 16,
                    "identity_params" => nothing, "simulator_argv" => "planted",
                    "error_model" => "planted", "qscore_model" => "planted",
                    "error_rate_estimate_diagnostic" => 0.0))
            return rs_dir
        end
        phase_a(args,
            root) = run_phase_a(
            D5OntAdaptiveGridDriver.parse_cli(vcat(
                ["--phase", "A",
                    "--allow-sha-mismatch"], args)), root)
        toy_rs = D5OntAdaptiveGridDriver.toy_read_sets()

        mktempdir() do root
            # Through run_phase_a: a cell index dispatches to the cell task, which
            # refuses before staging anything.
            msg = error_message(() -> phase_a(["--cell-index", "2", "--toy"], root))
            Test.@test msg !== nothing && occursin("is not prepared", msg)
            # The toy's 99%/100x read set is read-set index 1.
            Test.@test occursin("--read-set-index 1 --skip-iterative --toy", msg)
            Test.@test occursin("resubmit --cell-index 2", msg)
            Test.@test !isdir(joinpath(root, "refs"))
            Test.@test !isdir(joinpath(root, "readsets"))
            Test.@test !isdir(joinpath(root, "cells"))
            msg = error_message(() -> phase_a(["--toy"], root))
            Test.@test msg !== nothing &&
                       occursin("requires --read-set-index or --cell-index", msg)
            for bad in ("4", "-1")
                out = error_message(() -> select_cell(D5OntAdaptiveGridDriver.parse_cli(
                    ["--phase", "A", "--cell-index", bad, "--toy"])))
                Test.@test out !== nothing &&
                           occursin("--cell-index $(bad) out of range 0:3", out)
            end

            # A prepared read set loads with its recorded read count checked.
            cell = select_cell(D5OntAdaptiveGridDriver.parse_cli(
                ["--phase", "A", "--cell-index", "2", "--toy"]))
            rs = D5OntAdaptiveGridDriver.read_set_of(cell)
            Test.@test load_prepared_read_set(rs, root) === nothing
            rs_dir = plant_read_set(root, rs)
            Test.@test length(load_prepared_read_set(rs, root).records) == 2
            meta = JSON.parsefile(joinpath(rs_dir, "readset.json"))
            D5OntAdaptiveGridDriver.save_json_atomically(joinpath(rs_dir, "readset.json"),
                merge(meta, Dict("n_reads" => 3)))
            msg = error_message(() -> load_prepared_read_set(rs, root))
            Test.@test msg !== nothing && occursin("has 2 reads on disk", msg)
        end

        mktempdir() do root
            # A cell task on a prepared read set runs exactly the cell its index
            # names, and nothing else: cell 2 of the toy is ONT 99%/100x k29
            # (checkpointed ok, so it is served from cache, no assembly).
            alone = D5OntAdaptiveGridDriver.cells_run_alone(; toy = true)
            target = alone[3]
            Test.@test D5OntAdaptiveGridDriver.cell_id(target) ==
                       "ont__id99__cov100x__Lambda__seed42__iterative__k29__validation__scalable"
            plant_read_set(root, D5OntAdaptiveGridDriver.read_set_of(target))
            save_row(root, row_for(target))
            Test.@test phase_a(["--cell-index", "2", "--toy"], root) == ["ok"]
            Test.@test readdir(joinpath(root, "cells")) ==
                       [D5OntAdaptiveGridDriver.cell_id(target)]
        end

        mktempdir() do root
            # Full grid: the remedy names the FULL-grid read-set index and no --toy.
            # ONT 95%/10x Lambda seed 42 is read set 24 (24 ONT 90% read sets first).
            alone = D5OntAdaptiveGridDriver.cells_run_alone(; toy = false)
            j = findfirst(
                c -> D5OntAdaptiveGridDriver.read_set_id(
                    D5OntAdaptiveGridDriver.read_set_of(c)) ==
                     "ont__id95__cov10x__Lambda__seed42",
                alone) - 1
            msg = error_message(() -> phase_a(["--cell-index", string(j)], root))
            Test.@test msg !== nothing && occursin("is not prepared", msg)
            Test.@test occursin("--read-set-index 24 --skip-iterative)", msg)
            Test.@test !occursin("--toy", msg)
        end

        mktempdir() do root
            # Toy read-set task with --skip-iterative, through run_phase_a: its two
            # :none cells are checkpointed ok (served from cache, no assembly), and
            # its two :iterative cells are left to cell tasks.
            rs = toy_rs[1]
            plant_read_set(root, rs)
            cells = D5OntAdaptiveGridDriver.toy_cells(rs)
            for c in cells
                c.corrector === :none && save_row(root, row_for(c))
            end
            statuses = phase_a(["--read-set-index", "0", "--skip-iterative", "--toy"],
                root)
            Test.@test statuses == ["ok", "ok"]
            for c in cells
                Test.@test isfile(checkpoint_path(root, c)) == (c.corrector === :none)
            end
            # Run by hand without --skip-iterative, the same read-set task runs all
            # four cells (the :iterative two now checkpointed ok, so cached).
            for c in cells
                c.corrector === :iterative && save_row(root, row_for(c))
            end
            Test.@test phase_a(["--read-set-index", "0", "--toy"], root) ==
                       ["ok", "ok", "ok", "ok"]
        end

        mktempdir() do root
            # Full-grid ONT 99%/100x Lambda seed 42 (read set 66): 7 :none cells
            # (cached ok) and 8 validation cells above the assembler's k limit run
            # in the read-set task; the 8 infeasible rows need no assembler. The 7
            # grid and 8 assemblable validation :iterative cells are deferred.
            rs = D5OntAdaptiveGridDriver.read_sets()[67]
            Test.@test D5OntAdaptiveGridDriver.read_set_id(rs) ==
                       "ont__id99__cov100x__Lambda__seed42"
            plant_read_set(root, rs)
            cells = D5OntAdaptiveGridDriver.cells_for_read_set(rs; phase = :A)
            for c in cells
                c.corrector === :none && save_row(root, row_for(c))
            end
            statuses = phase_a(["--read-set-index", "66", "--skip-iterative"], root)
            Test.@test length(statuses) == 15
            Test.@test count(==("ok"), statuses) == 7
            Test.@test count(==(D5OntAdaptiveGridDriver.STATUS_INFEASIBLE_K), statuses) == 8
            for c in cells
                ran_here = c.corrector === :none ||
                           !D5OntAdaptiveGridDriver.assembler_accepts_k(c.k_policy)
                Test.@test isfile(checkpoint_path(root, c)) == ran_here
            end
        end
    end

    Test.@testset "--list prints both array sizes" begin
        function printed(toy)
            path, io = mktemp()
            try
                redirect_stdout(() -> print_counts(toy), io)
            finally
                close(io)
            end
            return read(path, String)
        end
        full = printed(false)
        Test.@test occursin("read sets: 96 -> --array=0-95", full)
        Test.@test occursin("(:iterative, --cell-index): 720 -> --array=0-719", full)
        toy = printed(true)
        Test.@test occursin("read sets (toy): 3 -> --array=0-2", toy)
        Test.@test occursin("(:iterative, --cell-index) (toy): 4 -> --array=0-3", toy)
        for (j, id) in enumerate((
            "ont__id95__cov10x__Lambda__seed42__iterative__k21__grid__scalable",
            "ont__id95__cov10x__Lambda__seed42__iterative__k31__grid__scalable",
            "ont__id99__cov100x__Lambda__seed42__iterative__k29__validation__scalable",
            "illumina__cov30x__Lambda__seed42__iterative__k21__illumina_control__scalable"))
            Test.@test occursin("  cell $(j - 1): $(id)\n", toy)
        end
    end
    Test.@testset "expected cells for the toy and the full grid" begin
        Test.@test length(expected_cells(true, false)) == 7
        Test.@test length(expected_cells(false, false)) == 1008 + 96 + 336
        Test.@test length(expected_cells(false, true)) == 1008 + 96 + 336 + 144 + 48
    end

    # Both site wrappers source the same task body (d5_grid_task.sh); each is
    # exercised end to end, so a wrapper that sources the wrong file, or none,
    # fails here.
    for wrapper in ("run_d5_grid_nersc.sbatch", "run_d5_grid_lrc.sbatch")
        Test.@testset "sbatch full array requires D5_FULL=1; default is the toy ($(wrapper))" begin
            sbatch = joinpath(@__DIR__, "..", "..", "benchmarking", wrapper)
            lines = readlines(sbatch)
            directives = filter(line -> startswith(line, "#SBATCH"), lines)
            Test.@test count(d -> occursin("--array=0-2", d), directives) == 1
            Test.@test !any(d -> occursin("0-95", d), directives)
            if Sys.which("bash") === nothing
                Test.@test_skip "bash unavailable"
            else
                # No D5_OUTDIR: a refused array exits in the gate, and an allowed
                # array dies on the unbound-variable check. Neither path reaches
                # module load, mkdir, or mktemp.
                function run_gate(; extra = Dict{String, String}())
                    env = Dict{String, String}(
                        "PATH" => get(ENV, "PATH", "/usr/bin"),
                        "HOME" => homedir(),
                        "D5_FULL" => "0")
                    merge!(env, extra)
                    out = IOBuffer()
                    err = IOBuffer()
                    proc = run(pipeline(ignorestatus(Cmd(`bash $sbatch`; env = env));
                        stdout = out, stderr = err))
                    return (; code = proc.exitcode,
                        log = String(take!(out)) * String(take!(err)))
                end
                refused = run_gate(;
                    extra = Dict(
                    "SLURM_ARRAY_TASK_ID" => "0", "SLURM_ARRAY_TASK_MAX" => "95"))
                Test.@test refused.code == 1
                Test.@test occursin("D5_FULL=1", refused.log)
                Test.@test occursin("outside the toy", refused.log)
                Test.@test !occursin("set D5_OUTDIR", refused.log)
                # A task id past the toy is refused even when the max is unset.
                past = run_gate(; extra = Dict("SLURM_ARRAY_TASK_ID" => "5"))
                Test.@test past.code == 1 && occursin("outside the toy", past.log)
                # D5_FULL must be exactly 1.
                almost = run_gate(;
                    extra = Dict(
                    "D5_FULL" => "true", "SLURM_ARRAY_TASK_MAX" => "95"))
                Test.@test almost.code == 1 && occursin("outside the toy", almost.log)
                garbage = run_gate(; extra = Dict("SLURM_ARRAY_TASK_MAX" => "abc"))
                Test.@test garbage.code == 1 &&
                           occursin("non-negative integers", garbage.log)
                # Toy max, and an approved full array, both pass the gate and stop
                # on the missing output root.
                toy = run_gate(;
                    extra = Dict(
                    "SLURM_ARRAY_TASK_ID" => "2", "SLURM_ARRAY_TASK_MAX" => "2"))
                Test.@test toy.code == 1 && occursin("set D5_OUTDIR", toy.log)
                Test.@test !occursin("outside the toy", toy.log)
                approved = run_gate(;
                    extra = Dict(
                    "D5_FULL" => "1", "SLURM_ARRAY_TASK_ID" => "95",
                    "SLURM_ARRAY_TASK_MAX" => "95"))
                Test.@test approved.code == 1 && occursin("set D5_OUTDIR", approved.log)
                Test.@test !occursin("outside the toy", approved.log)
                # Cell tasks: the toy has 4 :iterative cells, so 0-3 pass and 4 does not.
                cells_toy = run_gate(;
                    extra = Dict("D5_PHASE" => "A-cells", "SLURM_ARRAY_TASK_ID" => "3",
                    "SLURM_ARRAY_TASK_MAX" => "3"))
                Test.@test cells_toy.code == 1 && occursin("set D5_OUTDIR", cells_toy.log)
                cells_past = run_gate(;
                    extra = Dict("D5_PHASE" => "A-cells", "SLURM_ARRAY_TASK_ID" => "0",
                    "SLURM_ARRAY_TASK_MAX" => "719"))
                Test.@test cells_past.code == 1 &&
                           occursin("outside the toy", cells_past.log)
                Test.@test occursin("D5_PHASE=A-cells", cells_past.log)
                # The boundary itself: cell task 4 is the first index past the toy.
                cells_four = run_gate(;
                    extra = Dict("D5_PHASE" => "A-cells", "SLURM_ARRAY_TASK_ID" => "4",
                    "SLURM_ARRAY_TASK_MAX" => "4"))
                Test.@test cells_four.code == 1 &&
                           occursin("outside the toy", cells_four.log)
                # Read-set tasks keep the 0-2 toy: index 3 is refused there.
                rs_past = run_gate(; extra = Dict("SLURM_ARRAY_TASK_ID" => "3"))
                Test.@test rs_past.code == 1 && occursin("outside the toy", rs_past.log)
                bad = run_gate(; extra = Dict("D5_PHASE" => "cells"))
                Test.@test bad.code == 1 && occursin("D5_PHASE must be A, A-cells", bad.log)
                # An unknown phase is refused even on an approved (D5_FULL=1) array.
                bad_full = run_gate(; extra = Dict("D5_PHASE" => "cells", "D5_FULL" => "1"))
                Test.@test bad_full.code == 1 &&
                           occursin("D5_PHASE must be A, A-cells", bad_full.log)
            end
            # The wrapper -> task body -> runner path, end to end through the real
            # sbatch: a fake `julia` records its arguments and environment, and a
            # fake `module` succeeds unless told to fail. Needs bash and GNU
            # coreutils (mktemp --directory, realpath, rm --recursive).
            gnu = Sys.which("bash") !== nothing &&
                  try
                success(pipeline(`mktemp --version`; stdout = devnull, stderr = devnull))
            catch
                false
            end
            if !gnu
                Test.@test_skip "sbatch argv check needs bash and GNU coreutils"
            else
                mktempdir() do dir
                    bin = mkpath(joinpath(dir, "bin"))
                    argv_file = joinpath(dir, "argv.txt")
                    write(joinpath(bin, "julia"),
                        "#!/bin/bash\n" *
                        "if [[ \"\${1:-}\" == \"--version\" ]]; then\n" *
                        # A julia that misreports unless probed with LD_LIBRARY_PATH
                        # cleared, and one that cannot start (stderr, exit 1).
                        "  [[ -z \"\${LD_LIBRARY_PATH:-}\" ]] || " *
                        "{ echo 'julia version ldpath-leak'; exit 0; }\n" *
                        "  [[ -z \"\${D5_FAKE_BROKEN:-}\" ]] || " *
                        "{ echo \"\${D5_FAKE_BROKEN}\" >&2; exit 1; }\n" *
                        "  echo \"\${D5_FAKE_VERSION:-julia version 1.10.10}\"; exit 0\n" *
                        "fi\n" *
                        "printf '%s\\n' \"\$@\" > \"\${D5_FAKE_ARGV}\"\n" *
                        "printf '%s\\n' \"\${JULIA_NUM_THREADS:-unset}\" " *
                        "\"\${JULIA_DEPOT_PATH:-unset}\" \"\${LD_LIBRARY_PATH-unset}\" " *
                        "\"\${JULIAUP_CHANNEL-unset}\" > \"\${D5_FAKE_ARGV}.env\"\n" *
                        "[[ \"\${D5_FAKE_EXIT:-0}\" == 0 ]] || " *
                        "echo evidence > \"\${TMPDIR}/evidence.txt\"\n" *
                        "exit \"\${D5_FAKE_EXIT:-0}\"\n")
                    write(joinpath(bin, "module"),
                        "#!/bin/bash\nexit \"\${D5_FAKE_MODULE_EXIT:-0}\"\n")
                    chmod(joinpath(bin, "julia"), 0o755)
                    chmod(joinpath(bin, "module"), 0o755)
                    function runner_argv(phase, task; d5_args = "--toy", tag = "",
                            extra = Dict{String, String}(), script = sbatch,
                            submit_dir = dirname(sbatch))
                        rm(argv_file; force = true)
                        rm(argv_file * ".env"; force = true)
                        out = mkpath(joinpath(dir, "out_$(phase)_$(task)$(tag)"))
                        env = Dict{String, String}(
                            "PATH" => bin * ":" * get(ENV, "PATH", "/usr/bin"),
                            "HOME" => homedir(), "D5_FULL" => "0", "D5_PHASE" => phase,
                            "D5_ARGS" => d5_args, "D5_OUTDIR" => out,
                            "D5_TMPDIR" => mkpath(joinpath(dir, "tmp")),
                            "D5_FAKE_ARGV" => argv_file, "D5_PROC_CGROUP" => "/nonexistent",
                            # A login shell's library path, which the body clears.
                            "LD_LIBRARY_PATH" => "/bogus/lib",
                            "SLURM_SUBMIT_DIR" => submit_dir,
                            "SLURM_CPUS_PER_TASK" => "12",
                            "SLURM_ARRAY_TASK_ID" => string(task),
                            "SLURM_ARRAY_TASK_MAX" => string(task))
                        merge!(env, extra)
                        log = IOBuffer()
                        proc = run(pipeline(ignorestatus(Cmd(`bash $script`; env = env));
                            stdout = log, stderr = log))
                        argv = isfile(argv_file) ? readlines(argv_file) : String[]
                        julia_env = isfile(argv_file * ".env") ?
                                    readlines(argv_file * ".env") : String[]
                        return (; code = proc.exitcode, argv, julia_env, out,
                            log = String(take!(log)))
                    end
                    flag_value(argv, flag) =
                        let i = findfirst(==(flag), argv)
                            i === nothing ? nothing : argv[i + 1]
                        end
                    site_depot = Dict(
                        "run_d5_grid_nersc.sbatch" => "/global/cfs/projectdirs/m4269/cjprybol/.julia_depot",
                        "run_d5_grid_lrc.sbatch" => "/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot_d5")[wrapper]

                    a = runner_argv("A", 1)
                    Test.@test a.code == 0
                    Test.@test a.argv[1:4] ==
                               [
                        "--project=..", "d5_ont_adaptive_grid_run.jl", "--phase", "A"]
                    Test.@test flag_value(a.argv, "--read-set-index") == "1"
                    Test.@test "--skip-iterative" in a.argv && "--toy" in a.argv
                    Test.@test !("--cell-index" in a.argv)
                    # 4 Julia threads whatever --cpus-per-task is (12 here), the
                    # site's depot, LD_LIBRARY_PATH cleared, and on Lawrencium the
                    # pinned juliaup channel.
                    channel = wrapper == "run_d5_grid_lrc.sbatch" ? "1.10.10" : "unset"
                    Test.@test a.julia_env == ["4", site_depot, "", channel]
                    Test.@test occursin("threads 4; depot $(site_depot) (site default)",
                        a.log)
                    # Nothing inherited, or the same depot inherited: nothing ignored.
                    Test.@test !occursin("Ignored the submit shell", a.log)
                    same = runner_argv("A", 1; tag = "_samedepot",
                        extra = Dict("JULIA_DEPOT_PATH" => site_depot))
                    Test.@test same.code == 0 &&
                               !occursin("Ignored the submit shell", same.log)
                    # A login shell's depot is ignored, as Lawrencium's dotfiles
                    # export the group depot (trailing colon included).
                    inherited = runner_argv("A", 1; tag = "_depot",
                        extra = Dict("JULIA_DEPOT_PATH" => "/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot:"))
                    Test.@test inherited.julia_env[1:2] == ["4", site_depot]
                    Test.@test occursin(
                        "Ignored the submit shell's JULIA_DEPOT_PATH=" *
                        "/global/home/groups-sw/pc_mfnanofabio/cjprybol/.julia_depot:",
                        inherited.log)
                    # Only D5_DEPOT chooses another depot.
                    override = runner_argv("A", 1; tag = "_d5depot",
                        extra = Dict("D5_DEPOT" => "/elsewhere/depot"))
                    Test.@test override.julia_env[1:2] == ["4", "/elsewhere/depot"]
                    Test.@test occursin("depot /elsewhere/depot (D5_DEPOT override)",
                        override.log)

                    cells = runner_argv("A-cells", 2)
                    Test.@test cells.code == 0
                    Test.@test cells.argv[3:4] == ["--phase", "A"]
                    Test.@test flag_value(cells.argv, "--cell-index") == "2"
                    Test.@test !("--read-set-index" in cells.argv)
                    Test.@test !("--skip-iterative" in cells.argv)
                    # An exact element: the output root's path contains "A-cells".
                    Test.@test !("A-cells" in cells.argv)
                    # D5_ARGS is word-split: each flag and value is its own argument.
                    b = runner_argv("B", 0; d5_args = "--toy --tau 0.8 --tau-hat 0.6")
                    Test.@test b.code == 0 && b.argv[3:4] == ["--phase", "B"]
                    Test.@test flag_value(b.argv, "--read-set-index") == "0"
                    Test.@test flag_value(b.argv, "--tau") == "0.8"
                    Test.@test flag_value(b.argv, "--tau-hat") == "0.6"
                    Test.@test !("--skip-iterative" in b.argv)
                    m = runner_argv("merge", 0)
                    Test.@test m.code == 0 && m.argv[3:4] == ["--phase", "merge"]
                    Test.@test !any(f -> f in m.argv,
                        ("--read-set-index", "--cell-index", "--skip-iterative"))
                    # Every task hands the runner a node-local temp root and the
                    # resolved output root.
                    Test.@test startswith(flag_value(cells.argv, "--tmp-dir"),
                        realpath(joinpath(dir, "tmp")))
                    Test.@test flag_value(cells.argv, "--output-dir") ==
                               realpath(cells.out)

                    # The runner's exit code is the task's, and a failed task keeps
                    # what it left in the temp root.
                    retry = runner_argv("A", 1; tag = "_exit3",
                        extra = Dict("D5_FAKE_EXIT" => "3"))
                    Test.@test retry.code == 3 && occursin(">>> RETRY phase=A", retry.log)
                    evidence = joinpath(retry.out, "tmp_evidence")
                    Test.@test isdir(evidence) && length(readdir(evidence)) == 1 &&
                               isfile(joinpath(evidence, only(readdir(evidence)),
                                   "evidence.txt"))
                    blocked = runner_argv("merge", 0; tag = "_exit3",
                        extra = Dict("D5_FAKE_EXIT" => "3"))
                    Test.@test blocked.code == 3 &&
                               occursin(">>> BLOCKED merge", blocked.log)
                    failed = runner_argv("A-cells", 2; tag = "_exit7",
                        extra = Dict("D5_FAKE_EXIT" => "7"))
                    Test.@test failed.code == 7 && occursin("(exit 7)", failed.log)

                    # Any julia but exactly 1.10.10 is refused before the runner
                    # starts, including near misses a prefix or substring match
                    # would let through.
                    for (n, v) in enumerate(("julia version 1.12.1",
                        "julia version 1.10.13", "julia version 1.10.9",
                        "julia version 1.10.10-DEV"))
                        wrong = runner_argv("A", 1; tag = "_version$(n)",
                            extra = Dict("D5_FAKE_VERSION" => v))
                        Test.@test wrong.code == 1 && isempty(wrong.argv)
                        Test.@test occursin("requires julia 1.10.10", wrong.log)
                    end
                    # A julia that cannot start is refused with its own message
                    # quoted (e.g. juliaup's channel-not-installed error).
                    broken = runner_argv("A", 1; tag = "_broken",
                        extra = Dict("D5_FAKE_BROKEN" => "ERROR: `1.10.10` from environment variable JULIAUP_CHANNEL is not installed."))
                    Test.@test broken.code == 1 && isempty(broken.argv)
                    Test.@test occursin("requires julia 1.10.10", broken.log)
                    Test.@test occursin("JULIAUP_CHANNEL is not installed", broken.log)
                    if wrapper == "run_d5_grid_nersc.sbatch"
                        nomod = runner_argv("A", 1; tag = "_nomodule",
                            extra = Dict("D5_FAKE_MODULE_EXIT" => "1"))
                        Test.@test nomod.code == 1 && isempty(nomod.argv)
                        Test.@test occursin("d5_load_julia failed", nomod.log)
                    end

                    # sbatch runs a spooled copy of the wrapper, so the body is
                    # found through SLURM_SUBMIT_DIR, not next to the script.
                    spool = mkpath(joinpath(dir, "spool"))
                    spooled = joinpath(spool, "slurm_script")
                    cp(sbatch, spooled; force = true)
                    via_spool = runner_argv("A-cells", 2; tag = "_spool", script = spooled)
                    Test.@test via_spool.code == 0
                    Test.@test flag_value(via_spool.argv, "--cell-index") == "2"
                    # Submitted from the wrong directory: refused, nothing run.
                    elsewhere = runner_argv("A-cells", 2; tag = "_wrongdir",
                        script = spooled, submit_dir = mkpath(joinpath(dir, "elsewhere")))
                    Test.@test elsewhere.code == 1 && isempty(elsewhere.argv)
                    Test.@test occursin("submit from the benchmarking/ directory",
                        elsewhere.log)
                end
            end
        end
    end
end

end # module
