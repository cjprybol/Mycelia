# Unit test for the RGV `seed` column backfill (bead td-59o7).
#
# The tool's job is narrow and its risk is specific: writing a seed value that was
# never actually used would make an unpairable historical run LOOK pairable, which
# is worse than leaving it unpairable. So the tests concentrate on provenance
# discipline:
#
#   * a seed must be justified (explicit, recovered from a log, or an explicitly
#     labelled fall back to the documented default)
#   * an ambiguous log raises rather than picking one value
#   * the input file is never modified, and the sidecar records how the seed was
#     determined
#   * an existing `seed` column is never silently overwritten
#
# Dependency-free apart from CSV / DataFrames / JSON (direct Mycelia deps).

import Test
import CSV
import DataFrames
import JSON

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rgv_seed_backfill.jl"))

function _sbf_throws_with(f, needles::Vector{String})
    try
        f()
    catch e
        msg = sprint(showerror, e)
        for needle in needles
            Test.@test occursin(needle, msg)
        end
        return msg
    end
    Test.@test false
    return ""
end

# A pre-`seed` sweep CSV, matching the real one committed at
# benchmarking/results/rhizomorph_correction_validation_sweep_20260711_125013.csv.
const _SBF_LEGACY_HEADER = "reference,genome_len,error_rate,regime,readlen,tech," *
                           "target_coverage,effective_coverage,k,arm,ok,n_contigs,total_length," *
                           "largest_contig,n50,genome_fraction,runtime_s,mode,scale_metric_bases," *
                           "scale_floor_bases"

function _sbf_legacy_csv(path)
    open(path, "w") do io
        println(io, _SBF_LEGACY_HEADER)
        println(io,
            "synthetic_2000bp,2000,0.05,short-low-error,150,illumina,10.0,10.05," *
            "21,naive,true,935,49291,152,60,2464.5,9.606,SMOKE-ONLY,20100.0,1.0e6")
        println(io,
            "synthetic_2000bp,2000,0.05,short-low-error,150,illumina,10.0,10.05," *
            "21,iterative,true,244,14097,153,62,704.8,17.107,SMOKE-ONLY,20100.0,1.0e6")
    end
    return path
end

Test.@testset "RGV seed backfill (td-59o7)" begin
    Test.@testset "seed recovery from each run-log form" begin
        mktempdir() do dir
            # The sbatch echo form, as patched into run_rgv_sweep_{lrc,nersc}.sbatch.
            sbatch_log = joinpath(dir, "sbatch.log")
            write(sbatch_log,
                ">>> RUN rhizomorph_correction_validation_sweep.jl  Fri Jul 25\n" *
                "    ERR=0.01,0.05,0.10  READLEN=150,5000  COVERAGE=30x  K=21  " *
                "SEED=123  QUAST=true\n")
            Test.@test recover_seed_from_log(sbatch_log) == 123

            # The export form.
            export_log = joinpath(dir, "export.log")
            write(export_log, "export MYCELIA_RGV_SEED=\"456\"\n")
            Test.@test recover_seed_from_log(export_log) == 456

            # The harness banner form, as added to the sweep's config printout.
            banner_log = joinpath(dir, "banner.log")
            write(banner_log, "k              : 21\nSeed           : 42\nArms: naive\n")
            Test.@test recover_seed_from_log(banner_log) == 42

            # All three forms agreeing in one log is fine.
            combined = joinpath(dir, "combined.log")
            write(combined,
                "export MYCELIA_RGV_SEED=\"42\"\n  SEED=42  QUAST=true\nSeed           : 42\n")
            Test.@test recover_seed_from_log(combined) == 42

            # A log with no seed at all: nothing, not a guess.
            silent = joinpath(dir, "silent.log")
            write(silent, "ERR=0.01  READLEN=150  COVERAGE=30x  K=21  QUAST=true\n")
            Test.@test recover_seed_from_log(silent) === nothing

            _sbf_throws_with(() -> recover_seed_from_log(joinpath(dir, "nope.log")),
                ["run log not found"])
        end
    end

    Test.@testset "the log form the launchers NOW emit (a seed LIST)" begin
        # MYCELIA_RGV_SEED became a comma-list, so the wrappers echo
        # `SEED=42,123,456`. The old single-digit regex captured only `42`, which
        # would have labelled every row of a 3-seed replicate run with one seed —
        # making an unpairable table look pairable, the exact failure this tool is
        # supposed to prevent.
        mktempdir() do dir
            log = joinpath(dir, "multi.log")
            write(log,
                "    ERR=0.01,0.05,0.10  READLEN=150,5000  COVERAGE=30x  K=21  " *
                "SEED=42,123,456  QUAST=true\n")
            msg = _sbf_throws_with(() -> recover_seed_from_log(log),
                ["MORE THAN ONE seed", "42, 123, 456", "pairable"])
            Test.@test occursin("per-row `seed` column", msg)

            # A single-seed run in the new format still recovers cleanly.
            single = joinpath(dir, "single.log")
            write(single, "  ERR=0.01  K=21  SEED=123  QUAST=true\n")
            Test.@test recover_seed_from_log(single) == 123

            # And the export form with a list is caught too.
            exp = joinpath(dir, "export.log")
            write(exp, "export MYCELIA_RGV_SEED=\"42,123,456\"\n")
            _sbf_throws_with(() -> recover_seed_from_log(exp), ["MORE THAN ONE seed"])
        end
    end

    Test.@testset "a column present but entirely empty is filled, not refused" begin
        # Legacy CSVs can carry the header with empty cells; CSV.jl reads those as
        # `missing`. Treating that as a conflict made the documented workflow
        # impossible on exactly the historical tables this tool targets.
        mktempdir() do dir
            src = joinpath(dir, "empty_col.csv")
            write(src, "arm,k,seed,metric_source\nnaive,21,,\niterative,21,,\n")
            out, _,
            _ = backfill_seed(src;
                seed = 42, provenance = "explicit (test)", metric_source = "quast")
            df = CSV.read(out, DataFrames.DataFrame)
            Test.@test all(df.seed .== 42)
            Test.@test all(df.metric_source .== "quast")
        end
    end

    Test.@testset "a PARTIALLY populated column is completed, not reported as done" begin
        # The hole between "entirely empty" (filled above) and "conflicting"
        # (refused below): SOME rows carry the value and others are `missing`.
        # `unique(skipmissing(...))` hides the holes, so the equal-value branch
        # returned the frame UNCHANGED and reported success while leaving the table
        # unpairable — and `require_pairing_schema` in rgv_paired_wilcoxon.jl then
        # rejects it with a hint pointing back at THIS tool, so the operator
        # ping-pongs between the two forever.
        mktempdir() do dir
            src = joinpath(dir, "partial.csv")
            write(src,
                "arm,k,seed,metric_source\n" *
                "naive,21,42,quast\n" *
                "iterative,21,,\n")

            out, _,
            _ = backfill_seed(src;
                seed = 42, provenance = "explicit (test)", metric_source = "quast")
            df = CSV.read(out, DataFrames.DataFrame)
            # The contract: after a successful backfill NO row is left unlabelled.
            Test.@test count(ismissing, df.seed) == 0
            Test.@test all(df.seed .== 42)
            Test.@test count(ismissing, df.metric_source) == 0
            Test.@test all(df.metric_source .== "quast")

            # The same holds at the function level, for both columns.
            partial_seed = DataFrames.DataFrame(
                arm = ["naive", "iterative"],
                k = [21, 21],
                seed = [42, missing])
            filled = insert_seed_column!(partial_seed, 42)
            Test.@test collect(filled.seed) == [42, 42]

            partial_def = DataFrames.DataFrame(
                arm = ["naive", "iterative"],
                metric_source = ["quast", missing])
            filled_def = insert_definition_column!(partial_def, :metric_source, "quast")
            Test.@test collect(filled_def.metric_source) == ["quast", "quast"]

            # Filling holes is only legitimate because the populated rows AGREE with
            # the supplied value. A partially populated column whose populated value
            # DISAGREES is a conflict and must still raise — completing it would
            # rewrite an observed value, which is the failure this tool exists to
            # prevent.
            _sbf_throws_with(
                () -> insert_seed_column!(
                    DataFrames.DataFrame(k = [21, 21], seed = [42, missing]), 99),
                ["already has a `seed` column", "refusing to overwrite", "99"])
            _sbf_throws_with(
                () -> insert_definition_column!(
                    DataFrames.DataFrame(metric_source = ["quast", missing]),
                    :metric_source, "internal:quast-failed"),
                ["already has a `metric_source` column", "refusing to overwrite"])
            # And a fully populated CONFLICTING column is refused as before.
            _sbf_throws_with(
                () -> insert_seed_column!(
                    DataFrames.DataFrame(k = [21, 21], seed = [42, 99]), 42),
                ["already has a `seed` column", "refusing to overwrite"])
        end
    end

    Test.@testset "an ambiguous log raises instead of choosing" begin
        mktempdir() do dir
            log = joinpath(dir, "two-runs.log")
            write(log, "SEED=42\n...\nSEED=123\n")
            _sbf_throws_with(() -> recover_seed_from_log(log),
                ["MORE THAN ONE seed", "42, 123"])
        end
    end

    Test.@testset "backfill writes a new file and leaves the original alone" begin
        mktempdir() do dir
            src = _sbf_legacy_csv(joinpath(dir, "sweep.csv"))
            before = read(src, String)

            out, sidecar, n = backfill_seed(src; seed = 42, provenance = "explicit (test)")
            Test.@test n == 2
            Test.@test out == joinpath(dir, "sweep_seedbackfill.csv")
            Test.@test read(src, String) == before      # input untouched

            df = CSV.read(out, DataFrames.DataFrame)
            Test.@test "seed" in DataFrames.names(df)
            Test.@test all(df.seed .== 42)
            # Placed right after `k`, matching the live sweep's column order.
            cols = DataFrames.names(df)
            Test.@test cols[findfirst(==("k"), cols) + 1] == "seed"
            # Nothing else changed.
            Test.@test DataFrames.nrow(df) == 2
            Test.@test df.arm == ["naive", "iterative"]

            meta = JSON.parsefile(sidecar)
            Test.@test meta["seed"] == 42
            Test.@test occursin("explicit", meta["seed_provenance"])
            Test.@test meta["rows"] == 2
            Test.@test meta["bead"] == "td-59o7"
            Test.@test abspath(src) == meta["source_csv"]
        end
    end

    Test.@testset "existing seed column: idempotent if equal, refused if not" begin
        mktempdir() do dir
            src = _sbf_legacy_csv(joinpath(dir, "sweep.csv"))
            out, _, _ = backfill_seed(src; seed = 42, provenance = "explicit (test)")

            # Re-running with the SAME seed is a no-op, so a re-run is safe.
            df = CSV.read(out, DataFrames.DataFrame)
            Test.@test insert_seed_column!(copy(df), 42).seed == [42, 42]

            # A DIFFERENT seed would rewrite provenance — refuse.
            _sbf_throws_with(() -> insert_seed_column!(copy(df), 123),
                ["already has a `seed` column", "refusing to overwrite", "123"])
        end
    end

    Test.@testset "appends when there is no `k` column to anchor to" begin
        mktempdir() do dir
            src = joinpath(dir, "nok.csv")
            write(src, "arm,n50\nnaive,60\niterative,62\n")
            out, _, _ = backfill_seed(src; seed = 42, provenance = "explicit (test)")
            cols = DataFrames.names(CSV.read(out, DataFrames.DataFrame))
            Test.@test cols == ["arm", "n50", "seed"]
        end
    end

    Test.@testset "existing output is not clobbered without --force" begin
        mktempdir() do dir
            src = _sbf_legacy_csv(joinpath(dir, "sweep.csv"))
            out, _, _ = backfill_seed(src; seed = 42, provenance = "explicit (test)")
            _sbf_throws_with(
                () -> backfill_seed(src; seed = 42, provenance = "explicit (test)"),
                ["output already exists", "--force"])
            # With force it succeeds.
            out2, _,
            _ = backfill_seed(src;
                seed = 42, provenance = "explicit (test)", force = true)
            Test.@test out2 == out
        end
    end

    Test.@testset "missing input raises" begin
        mktempdir() do dir
            _sbf_throws_with(
                () -> backfill_seed(joinpath(dir, "nope.csv");
                    seed = 42, provenance = "explicit (test)"),
                ["CSV not found"])
        end
    end

    Test.@testset "the documented default matches the harness" begin
        # RGV_DEFAULT_SEED is mirrored from the sweep's MYCELIA_RGV_SEED default.
        # If the harness default changes, --assume-default-seed would attribute the
        # wrong seed to historical runs.
        sweep = read(
            joinpath(@__DIR__, "..", "..", "benchmarking",
                "rhizomorph_correction_validation_sweep.jl"), String)
        Test.@test occursin("MYCELIA_RGV_SEED", sweep)
        Test.@test occursin("[$(RGV_DEFAULT_SEED)]", sweep)   # the list default
        Test.@test RGV_DEFAULT_SEED == 42
    end

    Test.@testset "the sbatch wrappers echo SEED so future runs are recoverable" begin
        for wrapper in ("run_rgv_sweep_lrc.sbatch", "run_rgv_sweep_nersc.sbatch")
            text = read(
                joinpath(@__DIR__, "..", "..", "benchmarking", wrapper), String)
            Test.@test occursin("SEED=\${MYCELIA_RGV_SEED}", text)
        end
    end

    Test.@testset "C3: the documented backfill -> analyse workflow completes" begin
        # A seed-only backfill DEAD-ENDED: the analysis accepted the seed and then
        # refused the same file for carrying no metric-definition column, and
        # nothing could supply one. The workflow the tools document has to finish.
        mktempdir() do dir
            src = _sbf_legacy_csv(joinpath(dir, "sweep.csv"))
            out, sidecar,
            _ = backfill_seed(src;
                seed = 42, provenance = "explicit (test)",
                metric_source = "internal:quast-disabled", quast_min_contig = 200)
            df = CSV.read(out, DataFrames.DataFrame)
            for col in ("seed", "metric_source", "quast_min_contig")
                Test.@test col in DataFrames.names(df)
            end
            Test.@test all(df.metric_source .== "internal:quast-disabled")
            Test.@test all(df.quast_min_contig .== 200)

            # Backfilled definitions are DISTINGUISHABLE from observed ones.
            meta = JSON.parsefile(sidecar)
            Test.@test meta["metric_source_backfilled"] == "internal:quast-disabled"
            Test.@test meta["quast_min_contig_backfilled"] == 200

            # Never inferred: with no flags the columns stay absent, so the
            # analysis still fails closed rather than getting a fabricated value.
            out2, _,
            _ = backfill_seed(src;
                seed = 42, provenance = "explicit (test)",
                output_path = joinpath(dir, "seed_only.csv"))
            Test.@test !("metric_source" in DataFrames.names(CSV.read(out2, DataFrames.DataFrame)))
        end
    end

    Test.@testset "an existing definition column is not silently overwritten" begin
        mktempdir() do dir
            src = _sbf_legacy_csv(joinpath(dir, "sweep.csv"))
            out, _,
            _ = backfill_seed(src;
                seed = 42, provenance = "t", metric_source = "quast")
            df = CSV.read(out, DataFrames.DataFrame)
            # Same value is a no-op, so a re-run is safe: the column keeps its
            # value and no duplicate column is appended.
            again = insert_definition_column!(copy(df), :metric_source, "quast")
            Test.@test all(again.metric_source .== "quast")
            Test.@test DataFrames.ncol(again) == DataFrames.ncol(df)
            # A different value would assert a definition the run did not have.
            _sbf_throws_with(
                () -> insert_definition_column!(copy(df), :metric_source, "internal"),
                ["already has a `metric_source` column", "refusing to overwrite"])
        end
    end

    Test.@testset "B: a supplied definition is labelled IN THE CSV, not just the sidecar" begin
        # `--metric-source` states, on operator say-so, which definition produced
        # every metric in the file — a claim nothing verifies, and one no
        # `--run-log` equivalent could recover. Recording it only in the sidecar is
        # not a disclosure, because NO consumer reads the sidecar: the analysis
        # calls `load_sweep_csvs(csv_paths)` and nothing else. So a legacy table
        # that genuinely mixed QUAST-scored and degraded rows, backfilled with one
        # `--metric-source`, produced a file the guard certifies as uniform and a
        # report asserting the guard was enforced. The disclosure has to ride in
        # the CSV.
        mktempdir() do dir
            src = _sbf_legacy_csv(joinpath(dir, "sweep.csv"))
            out, sidecar,
            _ = backfill_seed(src;
                seed = 42, provenance = "explicit (test)",
                metric_source = "quast", quast_min_contig = 500)
            df = CSV.read(out, DataFrames.DataFrame)

            for col in ("metric_source_provenance", "quast_min_contig_provenance")
                Test.@test col in DataFrames.names(df)
            end
            # Every row of a column this tool SUPPLIED is an assertion.
            Test.@test all(occursin.("operator-asserted", df.metric_source_provenance))
            Test.@test all(occursin.("operator-asserted", df.quast_min_contig_provenance))
            # The sidecar names the columns that carry the disclosure, so the two
            # records cannot drift apart.
            meta = JSON.parsefile(sidecar)
            Test.@test "metric_source_provenance" in meta["definition_provenance_columns"]
            Test.@test occursin("operator-asserted",
                meta["definition_provenance_asserted_label"])

            # A backfill that supplies no definition asserts nothing, so it must not
            # gain a provenance column claiming otherwise.
            out2, _,
            _ = backfill_seed(src;
                seed = 42, provenance = "explicit (test)",
                output_path = joinpath(dir, "seed_only.csv"))
            names2 = DataFrames.names(CSV.read(out2, DataFrames.DataFrame))
            Test.@test !any(endswith.(names2, "_provenance"))
        end

        # Rows that ALREADY carried the value are evidence, not assertions: only
        # the rows the tool actually writes are labelled asserted.
        partial = DataFrames.DataFrame(
            arm = ["naive", "iterative"],
            metric_source = ["quast", missing])
        marked = insert_definition_column!(partial, :metric_source, "quast")
        Test.@test marked.metric_source_provenance ==
                   [DEFINITION_PROVENANCE_OBSERVED, DEFINITION_PROVENANCE_ASSERTED]

        # Re-running can only ADD assertions. Downgrading an already-asserted row to
        # "observed" would launder the claim away on a second pass.
        again = insert_definition_column!(copy(marked), :metric_source, "quast")
        Test.@test again.metric_source_provenance ==
                   [DEFINITION_PROVENANCE_OBSERVED, DEFINITION_PROVENANCE_ASSERTED]
    end

    Test.@testset "I1: the launchers actually request the replicate seeds" begin
        # The seed axis was asserted by the schema but not produced by anything: a
        # scalar MYCELIA_RGV_SEED and no loop meant a pairable multi-seed table
        # needed three manual submissions and an undocumented merge.
        sweep = read(
            joinpath(@__DIR__, "..", "..", "benchmarking",
                "rhizomorph_correction_validation_sweep.jl"), String)
        Test.@test occursin("_parse_int_list(get(ENV, \"MYCELIA_RGV_SEED\"", sweep)
        Test.@test occursin("for seed in seeds", sweep)
        for wrapper in ("run_rgv_sweep_lrc.sbatch", "run_rgv_sweep_nersc.sbatch")
            text = read(
                joinpath(@__DIR__, "..", "..", "benchmarking", wrapper), String)
            Test.@test occursin("MYCELIA_RGV_SEED:-42,123,456", text)
        end
    end
end
