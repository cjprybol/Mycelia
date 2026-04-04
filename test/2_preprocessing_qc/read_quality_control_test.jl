import Test
import FASTX
import BioSequences
import Mycelia

Test.@testset "Read Quality Control" begin
    make_record(id, seq, quality_scores = fill(UInt8(40), length(seq))) = FASTX.FASTQ.Record(
        id,
        seq,
        quality_scores
    )

    Test.@testset "gc_content_per_read — type safety regression (td-zyf99)" begin
        # Regression guard: the function must use FASTX.sequence(LongDNA{4}, record)
        # so BioSequences symbol comparisons work.  If it reverts to plain
        # FASTX.sequence(record) (returns String), all GC fractions would be 0.0.

        all_gc = make_record("all_gc", "GCGCGCGC")   # 100% GC
        all_at = make_record("all_at", "ATATATATAT")  # 0% GC
        half_gc = make_record("half_gc", "ATATGCGC")   # 50% GC
        empty_r = make_record("empty", "")

        records = [all_gc, all_at, half_gc, empty_r]
        fracs = Mycelia.gc_content_per_read(records)

        Test.@test length(fracs) == 4
        Test.@test fracs[1] ≈ 1.0
        Test.@test fracs[2] ≈ 0.0
        Test.@test fracs[3] ≈ 0.5
        Test.@test fracs[4] == 0.0   # empty record → 0.0, not NaN
    end

    Test.@testset "per-base and per-read quality summaries" begin
        records = [
            make_record("r1", "ATGC", UInt8[40, 40, 40, 40]),
            make_record("r2", "AT", UInt8[20, 20]),
            make_record("r3", "ATG", UInt8[30, 30, 30])
        ]

        Test.@test Mycelia.per_base_quality_scores(FASTX.FASTQ.Record[]) == Float64[]
        Test.@test Mycelia.per_base_quality_scores(records) == [30.0, 30.0, 35.0, 40.0]
        Test.@test Mycelia.per_read_quality_scores(records) == [40.0, 20.0, 30.0]
        Test.@test Mycelia.read_length_distribution(records) == [4, 2, 3]
    end

    Test.@testset "duplication statistics" begin
        records = [
            make_record("dup_a1", "AAAA"),
            make_record("dup_a2", "AAAA"),
            make_record("dup_g1", "GGGG"),
            make_record("dup_g2", "GGGG"),
            make_record("unique", "TTTT")
        ]

        stats = Mycelia.duplication_stats(records; min_fraction = 0.4, top_n = 1)
        Test.@test stats.unique_sequences == 3
        Test.@test stats.duplicate_reads == 2
        Test.@test stats.duplicate_fraction ≈ 0.4
        Test.@test length(stats.overrepresented) == 1
        Test.@test stats.overrepresented[1].count == 2
        Test.@test stats.overrepresented[1].fraction ≈ 0.4
        Test.@test stats.overrepresented[1].sequence in ("AAAA", "GGGG")

        empty_stats = Mycelia.duplication_stats(FASTX.FASTQ.Record[])
        Test.@test empty_stats.unique_sequences == 0
        Test.@test empty_stats.duplicate_reads == 0
        Test.@test empty_stats.duplicate_fraction == 0.0
        Test.@test isempty(empty_stats.overrepresented)
    end

    Test.@testset "analyze_fastq_quality and summarize_fastq" begin
        temp_dir = mktempdir()
        fastq_path = joinpath(temp_dir, "reads.fastq")
        records = [
            make_record("gc_high", "GCGC", UInt8[40, 40, 40, 40]),
            make_record("gc_low", "ATATAT", UInt8[30, 30, 30, 30, 30, 30])
        ]
        Mycelia.write_fastq(records = records, filename = fastq_path)

        quality_stats = Mycelia.analyze_fastq_quality(fastq_path)
        Test.@test quality_stats.n_reads == 2
        Test.@test quality_stats.mean_quality ≈ 35.0
        Test.@test quality_stats.mean_length ≈ 5.0
        Test.@test quality_stats.gc_content ≈ 40.0
        Test.@test quality_stats.quality_distribution.q20_percent ≈ 100.0
        Test.@test quality_stats.quality_distribution.q30_percent ≈ 100.0
        Test.@test quality_stats.quality_distribution.q40_percent ≈ 50.0

        summary_output = mktemp() do path, io
            summary = redirect_stdout(io) do
                Mycelia.summarize_fastq("test sample", fastq_path, 20)
            end
            flush(io)
            seekstart(io)
            text = read(io, String)
            (summary, text)
        end
        summary, summary_output = summary_output

        Test.@test summary.read_lengths == [4, 6]
        Test.@test summary.coverage ≈ 0.5
        Test.@test summary.error_rate_est ≈ Mycelia.q_value_to_error_rate(35.0)
        Test.@test occursin("test sample", summary_output)
        Test.@test occursin("reads: 2", summary_output)
        Test.@test occursin("estimated coverage: 0.5x", summary_output)

        missing_result = mktemp() do path, io
            missing_summary = redirect_stdout(io) do
                Mycelia.summarize_fastq("missing sample", joinpath(temp_dir, "missing.fastq"), 20)
            end
            flush(io)
            seekstart(io)
            text = read(io, String)
            (missing_summary, text)
        end
        missing_summary, missing_output = missing_result

        Test.@test missing_summary === nothing
        Test.@test occursin("missing file, skipping", missing_output)
    end
end
