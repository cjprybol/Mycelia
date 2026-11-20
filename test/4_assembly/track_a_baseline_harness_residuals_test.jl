# Track A harness residuals (td-p8e6).
#
# Pins the four contracts left in benchmarking/track_a_baseline_benchmark.jl after
# the atomic-publish split (#548):
#   (a) aggregation fails closed on a missing or unreadable cell instead of
#       publishing a smaller table;
#   (b) a cached status=error checkpoint is retried, not skipped as done;
#   (c) paired and grouped CLI flags reject partial or contradictory
#       combinations with exit code 2;
#   (d) parse_k is the single k-value helper.
#
# Run:
#   julia --project=. test/4_assembly/track_a_baseline_harness_residuals_test.jl

module TrackABaselineHarnessResidualsTest

import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "track_a_baseline_benchmark.jl"))

function _read_stale(root)
    return JSON.parsefile(aggregate_stale_path(root))
end

function _usage(args)
    try
        parse_track_a_invocation(args)
        return nothing
    catch e
        return e
    end
end

function _checkpoint(root, row)
    id = row_cell_id(row)
    cell_dir = joinpath(root, "cells", id)
    mkpath(cell_dir)
    save_cell_json(joinpath(cell_dir, "cell_result.json"), row)
    return id
end

function _sample(; org = "Lambda", tech = "ont", cov = 10, seed = 42, arm = "kmer",
        status = "ok", nga = 1000.0)
    return cell_row(org, "NC_001416", tech, cov, seed, arm;
        n_reads = 10, n_contigs = status == "empty_assembly" ? 0 : 1,
        wall_seconds = 1.0, peak_rss_bytes = 1, rss_baseline_bytes = -1,
        peak_rss_method = "unknown",
        metrics = merge(empty_metrics(),
            (NGA50 = nga, genome_fraction = 1.0, largest_contig = 1000.0)),
        status = status)
end

Test.@testset "track A harness residuals" begin

    Test.@testset "parse_k is the single odd-k helper" begin
        Test.@test parse_k("31") == 31
        Test.@test parse_k("19") == 19
        Test.@test parse_k(" 21 ") == 21
        for (raw, needle) in (
                ("20", "odd"), ("0", ">= 1"), ("-3", ">= 1"),
                ("abc", "integer"), ("19.0", "integer"))
            err = _usage(["--k", raw])
            # The invocation must go through parse_k, so the same needles fire
            # for a direct call and for `--k`.
            direct = try
                parse_k(raw)
                nothing
            catch e
                e
            end
            Test.@test direct isa UsageError
            Test.@test occursin(needle, direct.msg)
            Test.@test err isa UsageError
            Test.@test occursin(needle, err.msg)
            Test.@test track_a_usage_exit_code(["--k", raw]) == 2
        end
        # Default is the caller's job. parse_k does not invent 31.
        Test.@test parse_track_a_invocation(String[]).k == 31
        Test.@test parse_track_a_invocation(["--k=19"]).k == 19
    end

    Test.@testset "sibling flags exit 2" begin
        # Bare pairs. The message has to say why the default is refused: k and
        # the weighting are part of the cell id.
        bare_k = _usage(["--k"])
        Test.@test bare_k isa UsageError
        Test.@test occursin("bare trailing --k", bare_k.msg)
        Test.@test occursin("default to 31", bare_k.msg)
        Test.@test track_a_usage_exit_code(["--k"]) == 2

        bare_w = _usage(["--traversal-weighting"])
        Test.@test bare_w isa UsageError
        Test.@test occursin("bare trailing", bare_w.msg)
        Test.@test occursin("evidence", bare_w.msg)
        Test.@test track_a_usage_exit_code(["--traversal-weighting"]) == 2

        bare_out = _usage(["--output-dir"])
        Test.@test bare_out isa UsageError
        Test.@test occursin("--output-dir requires a value", bare_out.msg)

        # The next token is itself a flag: a partial pair, not a value.
        for args in (
                ["--k", "--traversal-weighting", "quality"],
                ["--output-dir", "--smoke"],
                ["--organisms", "--arms", "kmer"],
                ["--coverages", "--seeds", "42"])
            err = _usage(args)
            Test.@test err isa UsageError
            Test.@test occursin("next token", err.msg)
            Test.@test track_a_usage_exit_code(args) == 2
        end

        # Repeated flags used to be first-wins and drop the later value.
        twice = _usage(["--k", "19", "--k", "31"])
        Test.@test twice isa UsageError
        Test.@test occursin("more than once", twice.msg)
        Test.@test track_a_usage_exit_code(["--k=19", "--k", "31"]) == 2

        # `--flag=value` is a real spelling. Ignoring it used to run the default.
        Test.@test parse_track_a_invocation(["--k=19"]).k == 19
        Test.@test parse_track_a_invocation(
            ["--traversal-weighting=quality"]).traversal_weighting == "quality"
        empty_eq = _usage(["--k="])
        Test.@test empty_eq isa UsageError
        Test.@test occursin("non-empty", empty_eq.msg)
        even_eq = _usage(["--k=20"])
        Test.@test even_eq isa UsageError
        Test.@test occursin("odd", even_eq.msg)

        # `--smoke` replaces the shard group. Combining them used to ignore it.
        smoke_shard = _usage(["--smoke", "--organisms", "Lambda"])
        Test.@test smoke_shard isa UsageError
        Test.@test occursin("--smoke", smoke_shard.msg)
        Test.@test occursin("--organisms", smoke_shard.msg)
        Test.@test track_a_usage_exit_code(["--smoke", "--seeds", "42"]) == 2
        smoke_switch = _usage(["--smoke=true"])
        Test.@test smoke_switch isa UsageError
        Test.@test occursin("takes no value", smoke_switch.msg)

        # Unknown names used to be dropped, which ran a smaller matrix.
        unknown_org = _usage(["--organisms", "Lambda,Nope"])
        Test.@test unknown_org isa UsageError
        Test.@test occursin("Nope", unknown_org.msg)
        Test.@test occursin("smaller matrix", unknown_org.msg)
        Test.@test track_a_usage_exit_code(["--organisms", "Lambda,Nope"]) == 2
        unknown_tech = _usage(["--technologies", "illumina,ontt"])
        Test.@test unknown_tech isa UsageError
        Test.@test occursin("ontt", unknown_tech.msg)
        unknown_arm = _usage(["--arms", "qualmer,nope"])
        Test.@test unknown_arm isa UsageError
        Test.@test occursin("nope", unknown_arm.msg)

        bad_int = _usage(["--coverages", "30,10x"])
        Test.@test bad_int isa UsageError
        Test.@test occursin("integers", bad_int.msg)
        Test.@test occursin("10x", bad_int.msg)
        empty_entry = _usage(["--seeds", "42,,123"])
        Test.@test empty_entry isa UsageError
        Test.@test occursin("empty entry", empty_entry.msg)
        repeated = _usage(["--seeds", "42,42"])
        Test.@test repeated isa UsageError
        Test.@test occursin("repeats", repeated.msg)

        unknown_flag = _usage(["--aggregate-only"])
        Test.@test unknown_flag isa UsageError
        Test.@test occursin("unknown argument", unknown_flag.msg)
        positional = _usage(["Lambda"])
        Test.@test positional isa UsageError
        Test.@test occursin("unexpected argument", positional.msg)

        # The printed usage path returns nothing (exit 2) and does not throw.
        # Julia 1.10's redirect_stderr accepts a Pipe, not an IOBuffer.
        err_pipe = Pipe()
        parsed = redirect_stderr(err_pipe) do
            track_a_cli_or_usage(["--smoke", "--arms", "kmer"])
        end
        close(err_pipe.in)
        Test.@test parsed === nothing
        printed = read(err_pipe.out, String)
        Test.@test occursin("ERROR:", printed)
        Test.@test occursin("--smoke", printed)
        Test.@test occursin("--arms", printed)

        # Ctrl-C is not a bad flag.
        interrupted = try
            track_a_cli_or_usage(String[];
                parse_args = _ -> throw(InterruptException()))
            nothing
        catch e
            e
        end
        Test.@test interrupted isa InterruptException

        # Combinations the header actually advertises still parse, and smoke
        # still composes with --k / --traversal-weighting / --output-dir.
        header = parse_track_a_invocation(
            ["--coverages", "30,100", "--seeds", "42", "--technologies", "illumina,ont"])
        Test.@test header.coverages == [30, 100]
        Test.@test header.seeds == [42]
        Test.@test header.technologies == ["illumina", "ont"]
        Test.@test track_a_usage_exit_code(
            ["--organisms", "Lambda,T4", "--arms", "kmer"]) == 0
        orgs = parse_track_a_invocation(["--organisms", "T4,Lambda"])
        Test.@test [o[1] for o in orgs.organisms] == ["Lambda", "T4"]

        smoked = parse_track_a_invocation(
            ["--smoke", "--k", "19", "--traversal-weighting", "quality",
                "--output-dir", "/tmp/track_a_smoke"])
        Test.@test smoked.smoke
        Test.@test smoked.k == 19
        Test.@test smoked.traversal_weighting == "quality"
        Test.@test smoked.output_dir == "/tmp/track_a_smoke"
        Test.@test only(smoked.organisms)[1] == "Lambda"
        Test.@test smoked.technologies == ["illumina"]
        Test.@test smoked.coverages == [30]
        Test.@test smoked.seeds == [42]
        Test.@test smoked.arms == ["qualmer"]
        Test.@test track_a_usage_exit_code(String[
            "--smoke", "--k", "19", "--output-dir", "/tmp/track_a_smoke"]) == 0

        full = parse_track_a_invocation(String[])
        Test.@test full.k == 31
        Test.@test full.traversal_weighting == "evidence"
        Test.@test !full.smoke
        Test.@test length(full.organisms) == length(ORGANISMS)
        Test.@test full.technologies == TECHNOLOGIES
        Test.@test full.coverages == COVERAGES
        Test.@test full.seeds == SEEDS
        Test.@test full.arms == DECODER_ARMS
        Test.@test full.output_dir === nothing
        Test.@test track_a_usage_exit_code(String[]) == 0
    end

    Test.@testset "cached status=error is retried, not skipped" begin
        Test.@test cached_checkpoint_action("error") == :retry
        Test.@test cached_checkpoint_action("ok") == :skip
        Test.@test cached_checkpoint_action("empty_assembly") == :skip
        Test.@test "error" in RETRYABLE_STATUSES
        Test.@test !("empty_assembly" in RETRYABLE_STATUSES)
        Test.@test !("ok" in RETRYABLE_STATUSES)

        mktempdir() do dir
            ok = _sample()
            ok_id = _checkpoint(dir, ok)
            ok_decision = resume_decision(
                joinpath(dir, "cells", ok_id, "cell_result.json"))
            Test.@test ok_decision.action == :skip
            Test.@test ok_decision.row.status == "ok"

            err = error_row("Lambda", "NC_001416", "ont", 10, 123, "kmer")
            err_id = _checkpoint(dir, err)
            err_decision = resume_decision(
                joinpath(dir, "cells", err_id, "cell_result.json"))
            Test.@test err_decision.action == :retry
            Test.@test err_decision.row.status == "error"
            Test.@test err_id != ok_id

            empty = _sample(; tech = "illumina", cov = 30, seed = 456,
                arm = "qualmer", status = "empty_assembly", nga = 0.0)
            empty_id = _checkpoint(dir, empty)
            empty_decision = resume_decision(
                joinpath(dir, "cells", empty_id, "cell_result.json"))
            Test.@test empty_decision.action == :skip
            Test.@test empty_decision.row.status == "empty_assembly"

            Test.@test resume_decision(joinpath(dir, "absent.json")).action == :run

            bad = joinpath(dir, "truncated.json")
            write(bad, "{")
            bad_decision = resume_decision(bad)
            Test.@test bad_decision.action == :run
            Test.@test bad_decision.row === nothing

            interrupted = try
                resume_decision(bad; read = _ -> throw(InterruptException()))
                nothing
            catch e
                e
            end
            Test.@test interrupted isa InterruptException

            # Retry deletes reads/ (a truncated FASTQ would otherwise be reused)
            # and leaves the checkpoint in place.
            cell_dir = joinpath(dir, "cells", err_id)
            mkpath(joinpath(cell_dir, "reads"))
            write(joinpath(cell_dir, "reads", "ont_10x.fq.gz"), "truncated")
            prepare_cell_retry(cell_dir)
            Test.@test !isdir(joinpath(cell_dir, "reads"))
            Test.@test isfile(joinpath(cell_dir, "cell_result.json"))
        end
    end

    Test.@testset "aggregate refuses a smaller table" begin
        mktempdir() do dir
            ok = _sample()
            other = cell_row("T4", "NC_000866", "pacbio", 30, 123, "qualmer";
                n_reads = 10, n_contigs = 1, wall_seconds = 1.0, peak_rss_bytes = 1,
                rss_baseline_bytes = -1, peak_rss_method = "unknown",
                metrics = merge(empty_metrics(),
                    (NGA50 = 2000.0, genome_fraction = 1.0, largest_contig = 2000.0)),
                status = "ok")
            ok_id = _checkpoint(dir, ok)
            other_id = _checkpoint(dir, other)

            # A narrower in-memory slice must not replace the store.
            published = write_aggregate(dir, [ok])
            Test.@test DataFrames.nrow(published) == 2
            Test.@test Set(String.(published.organism)) == Set(["Lambda", "T4"])

            # A readable error row is a real checkpoint, not a hole.
            err = error_row("phi29", "NC_011048", "ont", 10, 42, "kmer")
            _checkpoint(dir, err)
            with_error = write_aggregate(dir, NamedTuple[])
            Test.@test DataFrames.nrow(with_error) == 3
            Test.@test "error" in String.(with_error.status)

            # A stray file in cells/ is not a cell.
            write(joinpath(dir, "cells", "notes.txt"), "ignore me")
            Test.@test DataFrames.nrow(write_aggregate(dir, NamedTuple[])) == 3

            sentinel = joinpath(dir, "track_a_results.tsv")
            write(sentinel, "SENTINEL\n")

            missing_id = "Lambda__illumina__50x__seed456__qualmer"
            mkpath(joinpath(dir, "cells", missing_id))
            broken_id = "T4__ont__100x__seed42__kmer"
            broken_dir = joinpath(dir, "cells", broken_id)
            mkpath(broken_dir)
            write(joinpath(broken_dir, "cell_result.json"), "{not json")

            warn_pipe = Pipe()
            thrown = redirect_stderr(warn_pipe) do
                try
                    write_aggregate(dir, NamedTuple[]; run_id = "residual-test-run")
                    nothing
                catch e
                    e
                end
            end
            close(warn_pipe.in)
            warning = read(warn_pipe.out, String)
            Test.@test thrown isa AggregateIncomplete
            Test.@test occursin("refusing to publish", thrown.msg)
            Test.@test occursin("smaller table", thrown.msg)
            Test.@test occursin(missing_id, thrown.msg)
            Test.@test occursin("missing", thrown.msg)
            Test.@test occursin(broken_id, thrown.msg)
            Test.@test occursin("incomplete", thrown.msg)
            Test.@test Set(g.cell_id for g in thrown.gaps) == Set([missing_id, broken_id])
            # The previous table is untouched, including the three good cells
            # that a skip-and-publish path would have written without the holes.
            Test.@test read(sentinel, String) == "SENTINEL\n"
            Test.@test !any(n -> startswith(n, "jl_"), readdir(dir))
            # Stderr names the kept table. The sidecar is what a later reader
            # of the results directory sees.
            Test.@test occursin("WARNING:", warning)
            Test.@test occursin(sentinel, warning)
            Test.@test occursin("STALE", warning)
            stale_path = aggregate_stale_path(dir)
            Test.@test isfile(stale_path)
            Test.@test basename(stale_path) == AGGREGATE_STALE_NAME
            stale = _read_stale(dir)
            Test.@test stale["stale_table"] == sentinel
            Test.@test stale["run_id"] == "residual-test-run"
            Test.@test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}$",
                stale["timestamp"])
            Test.@test occursin("refusing to publish", stale["reason"])
            Test.@test occursin(missing_id, stale["reason"])
            Test.@test occursin(broken_id, stale["reason"])
            kinds = Dict(c["cell_id"] => c["kind"] for c in stale["cells"])
            Test.@test kinds[missing_id] == "missing"
            Test.@test kinds[broken_id] == "incomplete"
            Test.@test occursin("no cell_result.json",
                only(c["reason"] for c in stale["cells"] if c["cell_id"] == missing_id))

            # An in-memory row whose id is the missing directory covers that
            # gap. The unreadable checkpoint is still unresolved, so the table
            # still does not publish.
            covering = _sample(; tech = "illumina", cov = 50, seed = 456, arm = "qualmer")
            Test.@test row_cell_id(covering) == missing_id
            still = redirect_stderr(devnull) do
                try
                    write_aggregate(dir, [covering]; run_id = "residual-test-run")
                    nothing
                catch e
                    e
                end
            end
            Test.@test still isa AggregateIncomplete
            Test.@test occursin(broken_id, still.msg)
            Test.@test !occursin(missing_id, still.msg)
            Test.@test read(sentinel, String) == "SENTINEL\n"
            # The sidecar is rewritten for the gaps that are still unresolved.
            rewritten = _read_stale(dir)
            Test.@test only(rewritten["cells"])["cell_id"] == broken_id
            Test.@test only(rewritten["cells"])["kind"] == "incomplete"
            Test.@test rewritten["run_id"] == "residual-test-run"
            Test.@test !occursin(missing_id, rewritten["reason"])

            rm(broken_dir; recursive = true)
            covered = write_aggregate(dir, [covering]; run_id = "residual-test-run")
            Test.@test !isfile(aggregate_stale_path(dir))
            Test.@test read(sentinel, String) != "SENTINEL\n"
            Test.@test occursin("organism", read(sentinel, String))
            Test.@test DataFrames.nrow(covered) == 4
            Test.@test any(i -> i == missing_id,
                (row_cell_id((;
                    organism = String(r.organism), technology = String(r.technology),
                    coverage = Int(r.coverage), seed = Int(r.seed),
                    decoder_arm = String(r.decoder_arm), k = Int(r.k),
                    traversal_weighting = String(r.traversal_weighting)))
                    for r in DataFrames.eachrow(covered)))
        end

        # A checkpoint that parses as JSON but fails the row schema is incomplete.
        mktempdir() do dir
            cell = joinpath(dir, "cells", "Lambda__ont__10x__seed42__kmer")
            mkpath(cell)
            write(joinpath(cell, "cell_result.json"), "{\"organism\":\"Lambda\"}")
            write(joinpath(dir, "track_a_results.tsv"), "SENTINEL\n")
            thrown = redirect_stderr(devnull) do
                try
                    write_aggregate(dir, NamedTuple[]; run_id = "schema-gap-run")
                    nothing
                catch e
                    e
                end
            end
            Test.@test thrown isa AggregateIncomplete
            Test.@test occursin("incomplete", thrown.msg)
            Test.@test occursin("missing required key", thrown.msg)
            Test.@test read(joinpath(dir, "track_a_results.tsv"), String) == "SENTINEL\n"
            schema_stale = _read_stale(dir)
            Test.@test schema_stale["run_id"] == "schema-gap-run"
            Test.@test schema_stale["stale_table"] == joinpath(dir, "track_a_results.tsv")
            cell = only(schema_stale["cells"])
            Test.@test cell["cell_id"] == "Lambda__ont__10x__seed42__kmer"
            Test.@test cell["kind"] == "incomplete"
            Test.@test occursin("missing required key", cell["reason"])
            Test.@test occursin("missing required key", schema_stale["reason"])
        end

        # Interrupt during the scan is not an incomplete cell.
        mktempdir() do dir
            cell = joinpath(dir, "cells", "Lambda__ont__10x__seed42__kmer")
            mkpath(cell)
            write(joinpath(cell, "cell_result.json"), "{}")
            thrown = try
                scan_cell_store(dir; read = _ -> throw(InterruptException()))
                nothing
            catch e
                e
            end
            Test.@test thrown isa InterruptException
        end
    end
end

end  # module TrackABaselineHarnessResidualsTest
