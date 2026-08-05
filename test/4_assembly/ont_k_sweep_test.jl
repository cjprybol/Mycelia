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
            # authoritative table.
            write(
                joinpath(cells, "Lambda__ont__k15__30x__seed42",
                    "cell_result.json"), "{ truncated")
            Test.@test DataFrames.nrow(write_aggregate(dir, NamedTuple[])) == 2
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
