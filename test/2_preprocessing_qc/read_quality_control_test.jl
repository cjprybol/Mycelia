import Test
import FASTX
import Mycelia

Test.@testset "Read Quality Control" begin
    q40(n) = fill(40, n)

    Test.@testset "gc_content_per_read — type safety regression (td-zyf99)" begin
        # Regression guard: the function must use FASTX.sequence(LongDNA{4}, record)
        # so BioSequences symbol comparisons work.  If it reverts to plain
        # FASTX.sequence(record) (returns String), all GC fractions would be 0.0.

        all_gc = Mycelia.fastq_record(identifier = "all_gc", sequence = "GCGCGCGC", quality_scores = q40(8))
        all_at = Mycelia.fastq_record(identifier = "all_at", sequence = "ATATATATAT", quality_scores = q40(10))
        half_gc = Mycelia.fastq_record(identifier = "half_gc", sequence = "ATATGCGC", quality_scores = q40(8))
        empty_r = Mycelia.fastq_record(identifier = "empty", sequence = "", quality_scores = Int[])

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
            Mycelia.fastq_record(identifier = "r1", sequence = "ATGC", quality_scores = [
                40, 40, 40, 40]),
            Mycelia.fastq_record(identifier = "r2", sequence = "AT", quality_scores = [
                20, 20]),
            Mycelia.fastq_record(identifier = "r3", sequence = "ATG", quality_scores = [
                30, 30, 30])
        ]

        Test.@test Mycelia.per_base_quality_scores(FASTX.FASTQ.Record[]) == Float64[]
        Test.@test Mycelia.per_base_quality_scores(records) == [30.0, 30.0, 35.0, 40.0]
        Test.@test Mycelia.per_read_quality_scores(records) == [40.0, 20.0, 30.0]
        Test.@test Mycelia.read_length_distribution(records) == [4, 2, 3]
    end

    Test.@testset "duplication statistics" begin
        records = [
            Mycelia.fastq_record(identifier = "dup_a1", sequence = "AAAA", quality_scores = q40(4)),
            Mycelia.fastq_record(identifier = "dup_a2", sequence = "AAAA", quality_scores = q40(4)),
            Mycelia.fastq_record(identifier = "dup_g1", sequence = "GGGG", quality_scores = q40(4)),
            Mycelia.fastq_record(identifier = "dup_g2", sequence = "GGGG", quality_scores = q40(4)),
            Mycelia.fastq_record(identifier = "unique", sequence = "TTTT", quality_scores = q40(4))
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
        mktempdir() do temp_dir
            fastq_path = joinpath(temp_dir, "reads.fastq")
            records = [
                Mycelia.fastq_record(identifier = "gc_high", sequence = "GCGC", quality_scores = q40(4)),
                Mycelia.fastq_record(identifier = "gc_low", sequence = "ATATAT", quality_scores = fill(30, 6))
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

            summary_result = mktemp() do _, io
                summary = redirect_stdout(io) do
                    Mycelia.summarize_fastq("test sample", fastq_path, 20)
                end
                flush(io)
                seekstart(io)
                text = read(io, String)
                (summary, text)
            end
            summary, summary_text = summary_result

            Test.@test summary.read_lengths == [4, 6]
            Test.@test summary.coverage ≈ 0.5
            Test.@test summary.error_rate_est ≈ Mycelia.q_value_to_error_rate(35.0)
            Test.@test occursin("test sample", summary_text)
            Test.@test occursin("reads: 2", summary_text)
            Test.@test occursin("estimated coverage: 0.5x", summary_text)

            missing_result = mktemp() do _, io
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
end
