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
        Test.@test nga50_status_for(aligned_low, 13_656, 24, true) == "censored_gf_below_50"

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
        Test.@test fixed.nga50_status == "censored_gf_below_50"
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
            "censored_gf_below_50", "quast_failed")
            Test.@test classify_outcome(0.0, 0.0, status) == "degenerate"
            Test.@test classify_outcome(48_000.0, 99.9, status) == "degenerate"
        end

        # Measured cells fall through the tier ladder.
        Test.@test classify_outcome(0.0, 9.7, "measured") == "degenerate"     # pilot ONT 30x
        Test.@test classify_outcome(552.0, 51.9, "measured") == "partial"     # pilot ONT 50x
        Test.@test classify_outcome(2356.0, 96.8, "measured") == "partial"    # pilot ONT 100x
        Test.@test classify_outcome(48_058.0, 99.9, "measured") == "near_complete"  # Illumina 30x

        # Boundary behaviour of the two thresholds that define "substantial":
        # genome fraction >= 90 AND NGA50 >= 10% of the genome.
        Test.@test classify_outcome(0.10 * GENOME_SIZE, 90.0, "measured") == "substantial"
        Test.@test classify_outcome(0.10 * GENOME_SIZE - 1, 90.0, "measured") == "partial"
        Test.@test classify_outcome(0.10 * GENOME_SIZE, 89.9, "measured") == "partial"
        # High genome fraction with fragmented contigs is NOT substantial — the
        # genome being present is not the same as it being assembled.
        Test.@test classify_outcome(500.0, 99.0, "measured") == "partial"
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
        complete = Dict(String(key) => (key in STR_KEYS ? "x" : 1) for key in ROW_KEYS)
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
