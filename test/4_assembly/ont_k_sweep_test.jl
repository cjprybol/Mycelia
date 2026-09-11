# ONT k-selection sweep contract tests (td-4e19d.28).
#
# These pin the pure helpers in benchmarking/ont_k_sweep.jl inside the Pkg.test()
# sweep, since benchmarking/ is not otherwise run in CI — the same pattern as
# track_a_baseline_benchmark_test.jl.
#
# What is under test is the CENSORING LOGIC, and it is worth testing because a
# defect there is silent: every function here returns a well-formed, plausible
# value on the failure path, and the failure path is exactly the regime the
# sweep exists to characterise. The Track-A pilot's ONT rows read "NGA50 = 0"
# not because anything measured zero but because its QUAST parser coerced
# QUAST's "-" to 0.0; that coercion is a one-character difference from correct
# and it made a censored floor indistinguishable from a measurement.
#
# Run:
#   julia --project=. test/4_assembly/ont_k_sweep_test.jl

# Wrapped in a module: runtests.jl includes every test file into one shared
# `Main`, and this driver defines top-level consts (COVERAGES, SEEDS, KS,
# ORGANISM, ...) that collide with benchmarking/track_a_baseline_benchmark.jl,
# which another test already includes there.
module OntKSweepTest

import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking", "ont_k_sweep.jl"))

const LAMBDA = genome_size_for("Lambda")
const T4 = genome_size_for("T4")

# Read a table back the way the harness wrote it, so an assertion about row
# count is an assertion about the FILE rather than about the DataFrame the
# writer happened to return.
read_tsv(path) = CSV.read(path, DataFrames.DataFrame;
    delim = '\t', missingstring = "NA")

Test.@testset "ONT k-sweep helpers" begin
    Test.@testset "contig_stats is independent of QUAST" begin
        # The low-k regime: many contigs, none long enough for QUAST to score.
        # These numbers must still be real, because they are the only
        # quantitative description of that regime.
        contigs = ["A"^12, "C"^338, "G"^14, "T"^499]
        stats = contig_stats(contigs, 500)
        Test.@test stats.n == 4
        Test.@test stats.n_ge_min == 0
        Test.@test stats.max_length == 499
        Test.@test stats.total_bp == 12 + 338 + 14 + 499

        # Boundary: min_contig is inclusive.
        Test.@test contig_stats(["A"^500], 500).n_ge_min == 1
        Test.@test contig_stats(["A"^499], 500).n_ge_min == 0

        empty_stats = contig_stats(String[], 500)
        Test.@test empty_stats.n == 0
        Test.@test empty_stats.max_length == 0
        Test.@test empty_stats.total_bp == 0
    end

    Test.@testset "nga50_status_for separates the censoring causes" begin
        measured = (; NGA50 = 1234.0, genome_fraction = 99.9)
        unaligned = (; NGA50 = missing, genome_fraction = missing)
        aligned_low = (; NGA50 = missing, genome_fraction = 31.663)

        # Nothing assembled at all.
        Test.@test nga50_status_for(unaligned, 0, 0, false) == "no_contigs"

        # Contigs exist but none reaches min_contig. This is the signature of
        # the low-k regime and MUST NOT be reported as a tool failure — QUAST
        # declining to score an assembly with no scorable contig is a correct
        # refusal, and conflating it with `quast_failed` would make the most
        # informative stratum of the sweep look like broken infrastructure.
        Test.@test nga50_status_for(unaligned, 149_115, 0, false) == "no_contigs_ge_min"

        # Scorable contigs exist, QUAST ran, and NOTHING aligned.
        Test.@test nga50_status_for(unaligned, 13_656, 42, true) == "censored_no_alignment"

        # Scorable contigs exist, QUAST ran, contigs DID align — but under the
        # 50% genome-fraction floor below which NGA50 has no definition. This
        # is the real ONT/k=21/30x cell: 31.663% of the genome was recovered.
        # Labelling it "nothing aligned" would assert the opposite of the
        # measurement sitting in the same row.
        Test.@test nga50_status_for(aligned_low, 13_656, 24, true) ==
                   "censored_partial_alignment"

        # Scorable contigs exist and QUAST genuinely failed.
        Test.@test nga50_status_for(unaligned, 13_656, 42, false) == "quast_failed"

        # The measured case.
        Test.@test nga50_status_for(measured, 5_357, 900, true) == "measured"
    end

    Test.@testset "reclassify corrects a stale label from recorded values" begin
        # A checkpoint written by the FIRST version of the classifier, which
        # called every absent NGA50 "censored_unaligned". The raw measurements
        # in the row are correct; only the derived label is wrong.
        stale = (; organism = "Lambda", accession = "NC_001416", technology = "ont",
            k = 21, coverage = 30, seed = 42, decoder_arm = "kmer",
            n_reads = 139, n_contigs = 44_000, asm_contigs_ge_min = 24,
            asm_max_contig = 4682, asm_total_bp = 1_000_000,
            quast_contigs = 24.0, total_length = 60_000.0, N50 = 900.0,
            largest_contig = 4682.0, NGA50 = missing,
            nga50_status = "censored_unaligned", NA50 = missing,
            largest_alignment = 762.0, genome_fraction = 31.663,
            duplication_ratio = 1.0, misassemblies = 0.0,
            unaligned_contigs = 10.0, unaligned_length = 5000.0,
            outcome = "degenerate", wall_seconds = 1.0, status = "ok")
        fixed = reclassify(stale)
        Test.@test fixed.nga50_status == "censored_partial_alignment"
        # Relabelling must be LOSSLESS — no measurement may be altered.
        Test.@test fixed.genome_fraction == 31.663
        Test.@test ismissing(fixed.NGA50)
        Test.@test fixed.asm_max_contig == 4682
        # Idempotent: reclassifying an already-correct row is a no-op.
        Test.@test reclassify(fixed).nga50_status == fixed.nga50_status
    end

    Test.@testset "classify_outcome never treats a censored floor as a measurement" begin
        # A censored NGA50 is degenerate regardless of what genome fraction
        # says. If QUAST could not compute NGA50 then nothing aligned, so a
        # nonzero genome fraction alongside it would be internally
        # inconsistent rather than evidence of partial success.
        for status in ("no_contigs", "no_contigs_ge_min", "censored_no_alignment",
            "censored_partial_alignment", "quast_failed")
            Test.@test classify_outcome(0.0, 0.0, status, LAMBDA) == "degenerate"
            Test.@test classify_outcome(48_000.0, 99.9, status, LAMBDA) == "degenerate"
        end

        # Measured cells fall through the tier ladder.
        Test.@test classify_outcome(0.0, 9.7, "measured", LAMBDA) == "degenerate"     # pilot ONT 30x
        Test.@test classify_outcome(552.0, 51.9, "measured", LAMBDA) == "partial"     # pilot ONT 50x
        Test.@test classify_outcome(2356.0, 96.8, "measured", LAMBDA) == "partial"    # pilot ONT 100x
        Test.@test classify_outcome(48_058.0, 99.9, "measured", LAMBDA) == "near_complete"  # Illumina 30x

        # Boundary behaviour of the two thresholds that define "substantial":
        # genome fraction >= 90 AND NGA50 >= 10% of the genome.
        Test.@test classify_outcome(0.10 * LAMBDA, 90.0, "measured", LAMBDA) ==
                   "substantial"
        Test.@test classify_outcome(0.10 * LAMBDA - 1, 90.0, "measured", LAMBDA) ==
                   "partial"
        Test.@test classify_outcome(0.10 * LAMBDA, 89.9, "measured", LAMBDA) == "partial"
        # High genome fraction with fragmented contigs is NOT substantial — the
        # genome being present is not the same as it being assembled.
        Test.@test classify_outcome(500.0, 99.0, "measured", LAMBDA) == "partial"
    end

    Test.@testset "outcome tiers scale with the genome, not a constant" begin
        # The tiers are FRACTIONS of the genome, so the same NGA50 must classify
        # differently on a 48.5 kb and a 168.9 kb reference. Getting this wrong
        # does not throw — it silently rescales every T4 row against Lambda's
        # size and produces a well-formed table of misclassified cells, which is
        # exactly the failure a second organism was added to avoid.
        Test.@test T4 > 3 * LAMBDA   # sanity: the two scales really do differ

        # NGA50 = 24,251 is 50% of Lambda but only ~14% of T4.
        half_lambda = 0.50 * LAMBDA
        Test.@test classify_outcome(half_lambda, 99.0, "measured", LAMBDA) ==
                   "near_complete"
        Test.@test classify_outcome(half_lambda, 99.0, "measured", T4) ==
                   "substantial"

        # NGA50 = 4,850 is 10% of Lambda but under 3% of T4.
        tenth_lambda = 0.10 * LAMBDA
        Test.@test classify_outcome(tenth_lambda, 92.0, "measured", LAMBDA) ==
                   "substantial"
        Test.@test classify_outcome(tenth_lambda, 92.0, "measured", T4) == "partial"
    end

    Test.@testset "write_aggregate cannot shrink the table (partial-run case)" begin
        # THE bug this guards: write_aggregate used to emit only the rows the
        # current process held, and OUTPUT_DIR defaults to the git-tracked
        # results directory — so `--smoke` (one cell) replaced a 240-row
        # deliverable with one row, silently, recoverable only via git.
        #
        # The happy path ("all cells in memory") and the honest-empty path
        # ("no cells at all") both looked fine. The defect lived in PARTIAL:
        # some on disk, some in memory. That is the case tested here.
        mktempdir() do dir
            cells = joinpath(dir, "cells")
            mkpath(cells)
            function put(org, tech, k, cov, seed; gf = 99.0, nga = 4000.0)
                row = cell_row(org, "ACC", tech, k, cov, seed;
                    n_reads = 10, asm = contig_stats(["A"^600], MIN_CONTIG),
                    metrics = merge(empty_metrics(),
                        (; NGA50 = nga, genome_fraction = gf, quast_contigs = 1.0)),
                    nga50_status = "measured", outcome = "partial",
                    wall_seconds = 1.0, status = "ok")
                id = cell_id_for(org, tech, k, cov, seed)
                mkpath(joinpath(cells, id))
                save_cell_json(joinpath(cells, id, "cell_result.json"), row)
                return row
            end
            a = put("Lambda", "ont", 15, 30, 42)
            put("Lambda", "ont", 15, 30, 123)
            put("Lambda", "ont", 15, 30, 456)

            # PARTIAL: one row in memory, three on disk. Must emit three.
            df = write_aggregate(dir, [a])
            Test.@test DataFrames.nrow(df) == 3

            # ABSENT: zero rows in memory. Must still emit three, not zero.
            Test.@test DataFrames.nrow(write_aggregate(dir, NamedTuple[])) == 3

            # The in-memory row must WIN over its on-disk twin, so a
            # just-recomputed cell supersedes a stale checkpoint.
            superseding = merge(a, (; wall_seconds = 999.0))
            df3 = write_aggregate(dir, [superseding])
            hit = df3[(df3.seed .== 42) .& (df3.k .== 15), :]
            Test.@test DataFrames.nrow(hit) == 1
            Test.@test hit.wall_seconds[1] == 999.0

            # An unreadable checkpoint must not take the whole aggregation down
            # — one truncated kilobyte would otherwise destroy a long run's
            # authoritative table. The READER skips it and keeps going...
            write(
                joinpath(cells, "Lambda__ont__k15__30x__seed42",
                    "cell_result.json"), "{ truncated")
            Test.@test length(load_all_checkpoints(dir)) == 2

            # ...and the WRITER refuses to publish the resulting 2-row view over
            # the 3-row table already on disk (td-4blm). Seed 42 was genuinely
            # measured; losing its checkpoint is not grounds to delete its row,
            # and cells/ is gitignored so the TSV is the only surviving copy.
            # The table is left intact, which is the property the comment above
            # is actually about — what the refusal costs is that the run stops.
            err = try
                write_aggregate(dir, NamedTuple[])
                nothing
            catch e
                e
            end
            Test.@test err isa ErrorException
            Test.@test occursin("refusing to shrink", err.msg)
            Test.@test occursin("Lambda / ont / 15 / 30 / 42", err.msg)
            Test.@test DataFrames.nrow(read_tsv(
                joinpath(dir, "ont_k_sweep_results.tsv"))) == 3
        end
    end

    Test.@testset "a fresh-clone invocation cannot shrink the committed table" begin
        # THE td-4blm bug. `write_aggregate`'s union is not sufficient on its
        # own: it unions against OUTPUT_DIR/cells/, which this harness's
        # .gitignore deliberately excludes, while the aggregate TSV is tracked.
        # So on a FRESH CLONE — table present, cells/ absent — the union
        # degenerates to "whatever this invocation computed" and the default
        # 96-cell grid would overwrite the committed 240-row deliverable,
        # silently,
        # recoverable only via git checkout.
        #
        # Every case below is the fresh-clone shape: a populated table and NO
        # cells/ directory at all.
        row(org,
            tech,
            k,
            cov,
            seed;
            wall = 1.0) = cell_row(
            org, "ACC", tech, k, cov, seed;
            n_reads = 10, asm = contig_stats(["A"^600], MIN_CONTIG),
            metrics = merge(empty_metrics(),
                (; NGA50 = 4000.0, genome_fraction = 99.0, quast_contigs = 1.0)),
            nga50_status = "measured", outcome = "partial",
            wall_seconds = wall, status = "ok")
        full = [row("Lambda", "ont", 15, 30, s) for s in (42, 123, 456)]
        append!(full, [row("T4", "ont", 15, 30, s) for s in (42, 123, 456)])
        results_name = "ont_k_sweep_results.tsv"

        # Seed the fixture with a PLAIN CSV.write, deliberately not through the
        # guard. The guard is the thing under test, and a fixture that depends
        # on it cannot run against the pre-change code at all — it fails to
        # resolve the symbol, which proves nothing about what the old code DID.
        # Seeding this way keeps "a partial grid over a full table is refused"
        # a real positive control: against the pre-change write_aggregate it
        # fails because the table is silently shrunk, not because it won't load.
        seed_table(dir) = CSV.write(joinpath(dir, results_name),
            DataFrames.DataFrame(full); delim = '\t', missingstring = "NA")

        Test.@testset "a partial grid over a full table is refused" begin
            mktempdir() do dir
                seed_table(dir)
                Test.@test !isdir(joinpath(dir, "cells"))   # fresh-clone shape

                # The default grid's shape: Lambda only, T4 absent — PLUS one
                # new T4 seed, so 4 rows are written while 3 keys are lost.
                #
                # The asymmetry is load-bearing. Writing plain `lambda_only`
                # makes lost == 3 and nrow == 3, and an oracle asserting
                # "has 3 key(s)" then cannot tell the two apart: swapping the
                # message to interpolate nrow(df) leaves it green. That is
                # precisely the confusion the message exists to avoid, so the
                # test must be able to see it.
                partial = vcat(filter(r -> r.organism == "Lambda", full),
                    [row("T4", "ont", 15, 30, 789)])
                err = try
                    write_aggregate(dir, partial)
                    nothing
                catch e
                    e
                end
                Test.@test err isa ErrorException
                Test.@test occursin("refusing to shrink", err.msg)
                # The count of DROPPED keys...
                Test.@test occursin("has 3 key(s)", err.msg)
                # ...and separately the count WRITTEN, which differs from it.
                Test.@test occursin("writing 4 rows", err.msg)
                Test.@test occursin("on disk 6", err.msg)
                Test.@test occursin("--allow-shrink", err.msg)
                # And the deliverable is untouched.
                Test.@test DataFrames.nrow(
                    read_tsv(joinpath(dir, results_name))) == 6
            end
        end

        Test.@testset "--allow-shrink is a real escape hatch" begin
            mktempdir() do dir
                seed_table(dir)
                lambda_only = DataFrames.DataFrame(
                    filter(r -> r.organism == "Lambda", full))
                write_table_guarded(joinpath(dir, results_name), lambda_only,
                    RESULTS_KEYCOLS; allow_shrink = true)
                Test.@test DataFrames.nrow(
                    read_tsv(joinpath(dir, results_name))) == 3
            end
        end

        Test.@testset "same keys with changed values is not a shrink" begin
            mktempdir() do dir
                seed_table(dir)
                touched = DataFrames.DataFrame(
                    [merge(r, (; wall_seconds = 999.0)) for r in full])
                write_table_guarded(joinpath(dir, results_name), touched,
                    RESULTS_KEYCOLS)
                back = read_tsv(joinpath(dir, results_name))
                Test.@test DataFrames.nrow(back) == 6
                Test.@test all(back.wall_seconds .== 999.0)
            end
        end

        Test.@testset "a superset write is not a shrink" begin
            mktempdir() do dir
                seed_table(dir)
                grown = DataFrames.DataFrame(
                    vcat(full, [row("T4", "ont", 15, 30, 789)]))
                write_table_guarded(joinpath(dir, results_name), grown,
                    RESULTS_KEYCOLS)
                Test.@test DataFrames.nrow(
                    read_tsv(joinpath(dir, results_name))) == 7
            end
        end

        Test.@testset "the guard does not fire on the case the union handles" begin
            # Regression fence around the EXISTING protection: with cells/
            # populated, a one-row invocation still emits the full table, and
            # the guard must stay silent. A guard that fired here would have
            # broken every ordinary resume.
            mktempdir() do dir
                cells = joinpath(dir, "cells")
                mkpath(cells)
                for r in full
                    id = cell_id_for(r.organism, r.technology, r.k, r.coverage,
                        r.seed)
                    mkpath(joinpath(cells, id))
                    save_cell_json(joinpath(cells, id, "cell_result.json"), r)
                end
                seed_table(dir)
                Test.@test DataFrames.nrow(write_aggregate(dir, [full[1]])) == 6
            end
        end

        Test.@testset "an unprovable write is refused, not waved through" begin
            # A write that cannot be SHOWN to be non-shrinking is refused.
            # Failing OPEN here would let the guard be defeated by exactly the
            # corruption it should catch. Three causes, three messages — and
            # they are asserted separately, because an earlier version of this
            # testset built one fixture and asserted two phrases against the
            # single error it produced, which is one assertion wearing two.
            guarded_error(f) =
                try
                    f()
                    nothing
                catch e
                    e
                end

            Test.@testset "on-disk table lacks a key column (schema change)" begin
                mktempdir() do dir
                    target = joinpath(dir, results_name)
                    CSV.write(target,
                        DataFrames.DataFrame(organism = ["Lambda"], k = [15]);
                        delim = '\t')
                    err = guarded_error(() -> write_table_guarded(
                        target, DataFrames.DataFrame(full), RESULTS_KEYCOLS))
                    Test.@test err isa ErrorException
                    Test.@test occursin("missing key column", err.msg)
                    Test.@test occursin("schema change", err.msg)
                    # Distinct from the corruption message, not a synonym.
                    Test.@test !occursin("could not be read", err.msg)
                    # Overridable, like every other refusal here.
                    write_table_guarded(target, DataFrames.DataFrame(full),
                        RESULTS_KEYCOLS; allow_shrink = true)
                    Test.@test DataFrames.nrow(read_tsv(target)) == 6
                end
            end

            Test.@testset "the frame being WRITTEN lacks a key column" begin
                # The symmetric case, which used to escape as a bare
                # ArgumentError with no path and no remedy.
                mktempdir() do dir
                    target = joinpath(dir, results_name)
                    seed_table(dir)
                    err = guarded_error(() -> write_table_guarded(
                        target,
                        DataFrames.DataFrame(organism = ["Lambda"], k = [15]),
                        RESULTS_KEYCOLS))
                    Test.@test err isa ErrorException
                    Test.@test occursin("the table being written", err.msg)
                    Test.@test occursin("--allow-shrink", err.msg)
                    Test.@test DataFrames.nrow(read_tsv(target)) == 6
                end
            end

            Test.@testset "an on-disk table that is genuinely unreadable" begin
                # CSV.jl is lenient, so this establishes what it actually does
                # with a binary file rather than assuming it throws. Either way
                # the write must not go through unproven: it either fails to
                # parse (corruption message) or parses into something without
                # the key columns (schema message). What must NOT happen is a
                # silent pass.
                mktempdir() do dir
                    target = joinpath(dir, results_name)
                    write(target, UInt8[0x00, 0xff, 0xfe, 0x00, 0x01, 0x02])
                    err = guarded_error(() -> write_table_guarded(
                        target, DataFrames.DataFrame(full), RESULTS_KEYCOLS))
                    Test.@test err isa ErrorException
                    Test.@test occursin("refusing to overwrite", err.msg)
                end
            end
        end

        Test.@testset "a Float64 key round-trips rather than reading as lost" begin
            # 80.5 is the value in the real identity ladder most likely to
            # expose a normalisation bug, so an identical rewrite of a table
            # keyed on it must not read as a loss.
            #
            # Scope note, because an earlier version of this comment overstated
            # it: this does NOT discriminate string keys from typed-tuple keys.
            # Julia's Set compares by isequal/hash, which is value-based, and a
            # typed setdiff across this same round trip was measured empty too.
            # The test pins the round trip; it is not evidence for the
            # representation choice. See table_row_keys' docstring for the real
            # argument.
            mktempdir() do dir
                target = joinpath(dir, "threshold.tsv")
                df = DataFrames.DataFrame(
                    cell_id = ["a", "a", "b", "b"],
                    min_identity = [95.0, 80.5, 95.0, 80.5],
                    value = [1, 2, 3, 4])
                CSV.write(target, df; delim = '\t', missingstring = "NA")
                write_table_guarded(target, df, (:cell_id, :min_identity))
                Test.@test DataFrames.nrow(read_tsv(target)) == 4
            end
        end

        Test.@testset "the partial cases, where set bugs actually live" begin
            # For an invariant over a SET, "all present" and "none present" are
            # the easy paths. The bugs live in SOME present, and in particular
            # in shapes where the ROW COUNT does not fall — a guard comparing
            # sizes rather than membership passes all of these.
            Test.@testset "simultaneous add and drop at equal row count" begin
                mktempdir() do dir
                    seed_table(dir)
                    swapped = vcat(
                        filter(r -> r.organism == "Lambda", full),
                        [row("T4", "ont", 15, 30, s) for s in (777, 888, 999)])
                    Test.@test length(swapped) == length(full)   # 6 vs 6
                    err = try
                        write_table_guarded(joinpath(dir, results_name),
                            DataFrames.DataFrame(swapped), RESULTS_KEYCOLS)
                        nothing
                    catch e
                        e
                    end
                    Test.@test err isa ErrorException
                    Test.@test occursin("has 3 key(s)", err.msg)
                    Test.@test occursin("writing 6 rows", err.msg)
                end
            end

            Test.@testset "a wholly disjoint key set is refused" begin
                mktempdir() do dir
                    seed_table(dir)
                    disjoint = DataFrames.DataFrame(
                        [row("T4", "illumina", 31, 100, s) for s in (1, 2, 3)])
                    err = try
                        write_table_guarded(joinpath(dir, results_name),
                            disjoint, RESULTS_KEYCOLS)
                        nothing
                    catch e
                        e
                    end
                    Test.@test err isa ErrorException
                    Test.@test occursin("has 6 key(s)", err.msg)
                end
            end

            Test.@testset "a 0-row frame over a populated table is refused" begin
                # The absent case on the WRITE side — an empty set is a subset
                # of everything, so a containment check written the wrong way
                # round waves this through.
                mktempdir() do dir
                    seed_table(dir)
                    empty_df = DataFrames.DataFrame(full)[1:0, :]
                    err = try
                        write_table_guarded(joinpath(dir, results_name),
                            empty_df, RESULTS_KEYCOLS)
                        nothing
                    catch e
                        e
                    end
                    Test.@test err isa ErrorException
                    Test.@test occursin("writing 0 rows", err.msg)
                    Test.@test DataFrames.nrow(
                        read_tsv(joinpath(dir, results_name))) == 6
                end
            end

            Test.@testset "a header-only table on disk blocks nothing" begin
                # The absent case on the DISK side. Nothing is lost by writing
                # over an empty table, so this must NOT refuse — a guard that
                # fired here would block every first real write.
                mktempdir() do dir
                    target = joinpath(dir, results_name)
                    CSV.write(target, DataFrames.DataFrame(full)[1:0, :];
                        delim = '\t', missingstring = "NA")
                    write_table_guarded(target, DataFrames.DataFrame(full),
                        RESULTS_KEYCOLS)
                    Test.@test DataFrames.nrow(read_tsv(target)) == 6
                end
            end
        end

        Test.@testset "a missing results table with siblings present is refused" begin
            # The hole created by unguarding the derived tables. They are
            # unguarded because the results table refuses first — which stops
            # being true the moment the results table is the one that is gone.
            # Reachable through this harness's own former advice ("move or
            # delete the file"), and the siblings CANNOT be rebuilt without
            # cells/, so the loss is permanent short of git checkout.
            mktempdir() do dir
                cells = joinpath(dir, "cells")
                mkpath(cells)
                for r in full
                    id = cell_id_for(r.organism, r.technology, r.k, r.coverage,
                        r.seed)
                    mkpath(joinpath(cells, id))
                    save_cell_json(joinpath(cells, id, "cell_result.json"), r)
                end
                # Siblings present, results table absent.
                CSV.write(joinpath(dir, "ont_k_sweep_summary.tsv"),
                    DataFrames.DataFrame(organism = ["Lambda"],
                        technology = ["ont"], k = [15], coverage = [30]);
                    delim = '\t')
                Test.@test !isfile(joinpath(dir, results_name))

                err = try
                    write_aggregate(dir, NamedTuple[])
                    nothing
                catch e
                    e
                end
                Test.@test err isa ErrorException
                Test.@test occursin("MISSING while its derived siblings", err.msg)
                Test.@test occursin("ont_k_sweep_summary.tsv", err.msg)
                # The sibling is untouched.
                Test.@test isfile(joinpath(dir, "ont_k_sweep_summary.tsv"))
            end

            # A genuinely fresh tree has neither, so this must NOT fire there —
            # otherwise it would block every first run.
            mktempdir() do dir
                cells = joinpath(dir, "cells")
                mkpath(cells)
                for r in full
                    id = cell_id_for(r.organism, r.technology, r.k, r.coverage,
                        r.seed)
                    mkpath(joinpath(cells, id))
                    save_cell_json(joinpath(cells, id, "cell_result.json"), r)
                end
                Test.@test DataFrames.nrow(
                    write_aggregate(dir, NamedTuple[])) == 6
            end
        end

        Test.@testset "the key separator cannot be forged from key values" begin
            # Pins the choice of _KEY_SEP, which was otherwise unpinned prose:
            # no key column contains an ordinary punctuation character today,
            # so swapping the separator for one changed nothing observable.
            #
            # A collision fixture only catches the separator it was built
            # against — an earlier version of this testset used ("a b","c") vs
            # ("a","b c"), which collides under a space and NOT under "/", so a
            # mutation to "/" survived it. So carry one colliding pair per
            # plausible separator: for candidate c, ("x"*c*"y", "z") and
            # ("x", "y"*c*"z") fold to the same string when c is the separator.
            #
            # Why a collision is fatal: two distinct keys become one, the
            # on-disk key set shrinks to match, and a genuinely shrinking write
            # then reads as non-shrinking — the guard fails OPEN.
            candidates = (' ', '/', '-', '_', ',', '|', ':')
            rows = NamedTuple[]
            for c in candidates
                push!(rows, (left = "x$(c)y", right = "z", v = 1))
                push!(rows, (left = "x", right = "y$(c)z", v = 2))
            end
            df = DataFrames.DataFrame(rows)
            n = DataFrames.nrow(df)

            # Every row must stay distinct. Under any candidate separator the
            # matching pair collapses, so this drops below n.
            Test.@test length(table_row_keys(df, (:left, :right))) == n

            mktempdir() do dir
                target = joinpath(dir, "sep.tsv")
                # And the write must SUCCEED: check_keycols_are_a_key turns a
                # separator collision into a loud "is not a key for this table"
                # refusal, so a bad separator cannot pass quietly here either.
                write_table_guarded(target, df, (:left, :right))
                Test.@test DataFrames.nrow(read_tsv(target)) == n

                # Dropping one row must still be refused.
                err = try
                    write_table_guarded(target, df[1:(n - 1), :],
                        (:left, :right))
                    nothing
                catch e
                    e
                end
                Test.@test err isa ErrorException
                Test.@test occursin("refusing to shrink", err.msg)
                Test.@test DataFrames.nrow(read_tsv(target)) == n
            end
        end

        Test.@testset "declared key columns must actually be a key" begin
            # A Set discards multiplicity, so without this check a write that
            # collapses distinct rows into duplicates has an unchanged key set
            # and sails through while silently dropping rows.
            mktempdir() do dir
                target = joinpath(dir, "dup.tsv")
                dup = DataFrames.DataFrame(
                    organism = ["Lambda", "Lambda"], k = [15, 15], v = [1, 2])
                err = try
                    write_table_guarded(target, dup, (:organism, :k))
                    nothing
                catch e
                    e
                end
                Test.@test err isa ErrorException
                Test.@test occursin("is not a key for this table", err.msg)
                Test.@test occursin("2 rows collapse to 1", err.msg)
            end
        end

        Test.@testset "the published write is atomic" begin
            # The guard proves a write non-shrinking and then replaces the
            # committed file. Truncating in place would mean a crash mid-write
            # produced exactly the loss the guard exists to prevent, so the
            # write lands via a temp file and an atomic rename.
            good = DataFrames.DataFrame(organism = ["Lambda", "T4"],
                k = [15, 21], v = [1, 2])

            Test.@testset "a failure inside write! leaves the original intact" begin
                mktempdir() do dir
                    target = joinpath(dir, "atomic.tsv")
                    CSV.write(target, good; delim = '\t', missingstring = "NA")
                    Test.@test_throws ErrorException publish_atomically(target) do tmp
                        write(tmp, "partial")
                        error("simulated crash mid-write")
                    end
                    Test.@test DataFrames.nrow(read_tsv(target)) == 2
                    Test.@test isempty(filter(f -> occursin(".tmp.", f),
                        readdir(dir)))
                end
            end

            Test.@testset "the publish step never deletes before it replaces" begin
                # THIS is the assertion that discriminates the implementation.
                # The version above passes under `mv(tmp, path; force = true)`
                # too, because the failure happens before the publish step is
                # reached — so it could not tell the two apart.
                #
                # Julia's `mv` does rm(dst) THEN rename, so a failure in the
                # publish step itself leaves the committed table DELETED.
                # Measured: with an unrenameable source, mv leaves
                # isfile(dst) == false. Base.Filesystem.rename leaves it
                # byte-for-byte intact. Simulate a publish-step failure by
                # removing the temp out from under it, which is what a crash
                # between write and rename looks like from the target's side.
                mktempdir() do dir
                    target = joinpath(dir, "atomic.tsv")
                    CSV.write(target, good; delim = '\t', missingstring = "NA")
                    before = read(target, String)
                    threw = false
                    try
                        publish_atomically(target) do tmp
                            write(tmp, "x")
                            rm(tmp; force = true)   # publish step now must fail
                        end
                    catch
                        threw = true
                    end
                    Test.@test threw
                    # Under `mv(...; force = true)` the file is gone here.
                    Test.@test isfile(target)
                    Test.@test read(target, String) == before
                end
            end

            Test.@testset "a successful publish replaces the contents" begin
                mktempdir() do dir
                    target = joinpath(dir, "atomic.tsv")
                    CSV.write(target, good; delim = '\t', missingstring = "NA")
                    bigger = DataFrames.DataFrame(
                        organism = ["Lambda", "T4", "T4"], k = [15, 21, 31],
                        v = [1, 2, 3])
                    publish_atomically(target) do tmp
                        CSV.write(tmp, bigger; delim = '\t',
                            missingstring = "NA")
                    end
                    Test.@test DataFrames.nrow(read_tsv(target)) == 3
                    Test.@test isempty(filter(f -> occursin(".tmp.", f),
                        readdir(dir)))
                end
            end
        end
    end

    Test.@testset "write_summary survives a tree in which nothing succeeded" begin
        # The reachable case: `--aggregate-only` over a tree where every
        # checkpoint carries error / quast_failed / empty_assembly. That is the
        # RECOVERY path — the one an operator reaches for precisely when the run
        # went badly — and it used to raise a column-not-found error from `sort!`
        # on a 0x0 DataFrame, aborting before verdict_stats.tsv was written. The
        # operator got a DataFrames internal error instead of the diagnosis
        # "nothing succeeded", and lost the second output too.
        function row_with(status, nga50_status)
            return cell_row("Lambda", "NC_001416", "ont", 21, 30, 42;
                n_reads = 10, asm = contig_stats(["A"^600], MIN_CONTIG),
                metrics = empty_metrics(), nga50_status = nga50_status,
                outcome = "degenerate", wall_seconds = 1.0, status = status)
        end

        mktempdir() do dir
            failed = DataFrames.DataFrame([
                row_with("quast_failed", "quast_failed"),
                row_with("error", "quast_failed"),
                row_with("empty_assembly", "no_contigs")
            ])
            summary = write_summary(dir, failed)
            Test.@test DataFrames.nrow(summary) == 0
            # An empty summary must still carry the SCHEMA. A 0x0 frame writes a
            # headerless file, which no downstream reader can distinguish from a
            # truncated write.
            Test.@test Tuple(Symbol.(DataFrames.names(summary))) == SUMMARY_KEYS
            out = joinpath(dir, "ont_k_sweep_summary.tsv")
            Test.@test isfile(out)
            Test.@test split(first(eachline(out)), '\t') == collect(String.(SUMMARY_KEYS))

            # The aggregation must reach its second output rather than aborting
            # in the first — losing verdict_stats.tsv was the actual damage.
            write_verdict_stats(dir, failed)
            Test.@test isfile(joinpath(dir, "verdict_stats.tsv"))
        end

        # SCHEMA SYNC: the hardcoded empty-case column list must equal the
        # columns the populated path actually produces, or the empty file
        # silently documents a table that no longer exists.
        mktempdir() do dir
            ok = DataFrames.DataFrame([
                cell_row("Lambda", "NC_001416", "ont", 21, 30, 42;
                n_reads = 10, asm = contig_stats(["A"^600], MIN_CONTIG),
                metrics = merge(empty_metrics(),
                    (; NGA50 = 4000.0, genome_fraction = 99.0,
                        quast_contigs = 1.0)),
                nga50_status = "measured", outcome = "partial",
                wall_seconds = 1.0, status = "ok")
            ])
            populated = write_summary(dir, ok)
            Test.@test DataFrames.nrow(populated) == 1
            Test.@test Tuple(Symbol.(DataFrames.names(populated))) == SUMMARY_KEYS
        end
    end

    Test.@testset "classify_outcome rejects a censored value instead of zeroing it" begin
        # Both call sites used to coerce `missing` to 0.0 before calling this,
        # reintroducing the collapse the harness exists to prevent. It was live
        # for genome_fraction, because nga50_status_for returns "measured"
        # without ever inspecting genome_fraction.
        threw = false
        try
            classify_outcome(4000.0, missing, "measured", LAMBDA)
        catch err
            threw = true
            Test.@test occursin("censored", sprint(showerror, err))
        end
        Test.@test threw

        threw2 = false
        try
            classify_outcome(missing, 99.0, "measured", LAMBDA)
        catch
            threw2 = true
        end
        Test.@test threw2

        # A censored row that is correctly LABELLED censored still short-circuits
        # to degenerate without reaching the guard.
        Test.@test classify_outcome(missing, missing, "censored_no_alignment",
            LAMBDA) == "degenerate"
    end

    Test.@testset "quast_failed is retryable and error_row matches the schema" begin
        # A QUAST infrastructure failure used to be recorded as status="ok",
        # which let it pass the summary filter as evidence for degeneracy, print
        # "Grid is COMPLETE", and be cached forever.
        Test.@test "quast_failed" in RETRYABLE_STATUSES
        Test.@test "error" in RETRYABLE_STATUSES
        # empty_assembly is a real measurement — re-running only reproduces it.
        Test.@test !("empty_assembly" in RETRYABLE_STATUSES)

        # error_row is built by hand and is called from inside a catch handler,
        # so schema drift there would abort the sweep from its own recovery path.
        er = error_row("Lambda", "NC_001416", "ont", 21, 30, 42)
        Test.@test keys(er) === ROW_KEYS
        Test.@test er.status == "error"
    end

    Test.@testset "genome_size_for refuses an unknown organism" begin
        Test.@test genome_size_for("Lambda") == 48_502
        Test.@test genome_size_for("T4") == 168_903
        # Defaulting here would rescale the tiers silently, so it must throw.
        threw = false
        try
            genome_size_for("NotAGenome")
        catch err
            threw = true
            Test.@test occursin("NotAGenome", sprint(showerror, err))
        end
        Test.@test threw
    end

    Test.@testset "parse_quast_metrics maps QUAST's \"-\" to missing, never 0.0" begin
        # This is the exact defect the sweep was written around: the Track-A
        # pilot's parser did `something(tryparse(Float64, "-"), 0.0)`, turning
        # "QUAST could not compute this" into a measured zero.
        mktempdir() do dir
            report = joinpath(dir, "report.tsv")
            write(report,
                "Assembly\tcontigs\n" *
                "# contigs\t13656\n" *
                "Total length\t6800000\n" *
                "Largest contig\t16434\n" *
                "N50\t612\n" *
                "NGA50\t-\n" *
                "Genome fraction (%)\t-\n" *
                "# misassemblies\t0\n")
            metrics = parse_quast_metrics(report)
            Test.@test metrics.quast_contigs == 13656.0
            Test.@test metrics.largest_contig == 16434.0
            Test.@test metrics.N50 == 612.0
            Test.@test ismissing(metrics.NGA50)
            Test.@test ismissing(metrics.genome_fraction)
            # A metric absent from the report is also missing, not zero.
            Test.@test ismissing(metrics.NA50)
            # "# misassemblies" = 0 is REAL here, but it is only interpretable
            # where NGA50 was measured; QUAST reports 0 on unaligned contigs.
            Test.@test metrics.misassemblies == 0.0
        end

        # An absent report yields all-missing rather than all-zero.
        Test.@test ismissing(parse_quast_metrics("/nonexistent/report.tsv").NGA50)
    end

    Test.@testset "canonical refuses to default a missing measurement" begin
        # `organism` must be a real one: canonical() reclassifies on read, and
        # reclassification looks the genome size up by name rather than
        # defaulting, so a placeholder organism is (correctly) rejected.
        complete = Dict(String(key) => (key in STR_KEYS ? "x" : 1) for key in ROW_KEYS)
        complete["organism"] = "Lambda"
        Test.@test canonical(complete) isa NamedTuple

        # Every column here is a measurement or a grouping key. Silently
        # defaulting one would let a truncated checkpoint enter the analysis as
        # though it were a result — the same class of failure as the coercion
        # above, one layer down.
        truncated = copy(complete)
        delete!(truncated, "genome_fraction")
        threw = false
        try
            canonical(truncated)
        catch err
            threw = true
            Test.@test occursin("genome_fraction", sprint(showerror, err))
        end
        Test.@test threw

        # JSON round-trips absent values as `nothing`; those must land as
        # `missing`, not as 0.0.
        with_null = copy(complete)
        with_null["NGA50"] = nothing
        Test.@test ismissing(canonical(with_null).NGA50)
    end
end

end  # module
