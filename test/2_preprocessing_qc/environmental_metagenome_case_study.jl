import Test
import Mycelia
import DataFrames
import FASTX
import LinearAlgebra

function _write_environmental_demo_fastq(path::String, records::Vector{FASTX.FASTQ.Record})
    Mycelia.write_fastq(records = records, filename = path)
    return path
end

Test.@testset "Environmental metagenome helpers" begin
    case_study = Mycelia.environmental_metagenome_case_study(max_runs = 3)
    Test.@test case_study.case_study_id == :yukon_river_pilot_station
    Test.@test case_study.habitat == "freshwater river water"
    Test.@test DataFrames.nrow(case_study.samples) == 3
    Test.@test collect(case_study.samples.run) == ["SRR12197245", "SRR12197243", "SRR12197244"]

    mktempdir() do dir
        fastq_a = _write_environmental_demo_fastq(joinpath(dir, "sample_a.fastq"), [
            FASTX.FASTQ.Record("a1", "ACGTACGTACGT", "IIIIIIIIIIII"),
            FASTX.FASTQ.Record("a2", "ACGTACGTACGA", "HHHHHHHHHHHH"),
            FASTX.FASTQ.Record("a3", "ACGTACGTACGG", "GGGGGGGGGGGG")
        ])
        fastq_b = _write_environmental_demo_fastq(joinpath(dir, "sample_b.fastq"), [
            FASTX.FASTQ.Record("b1", "GGGGCCCCAAAA", "IIIIIIIIIIII"),
            FASTX.FASTQ.Record("b2", "GGGGCCCCAAAT", "HHHHHHHHHHHH"),
            FASTX.FASTQ.Record("b3", "GGGGCCCCAAAG", "GGGGGGGGGGGG")
        ])
        fastq_c = _write_environmental_demo_fastq(joinpath(dir, "sample_c.fastq"), [
            FASTX.FASTQ.Record("c1", "TATATATATATA", "IIIIIIIIIIII"),
            FASTX.FASTQ.Record("c2", "TATATATATATC", "HHHHHHHHHHHH"),
            FASTX.FASTQ.Record("c3", "TATATATATATG", "GGGGGGGGGGGG")
        ])

        fastqs = [fastq_a, fastq_b, fastq_c]
        sample_names = ["sample_a", "sample_b", "sample_c"]

        quality_summary = Mycelia.environmental_fastq_quality_summary(
            fastqs;
            sample_names = sample_names
        )
        Test.@test DataFrames.nrow(quality_summary) == 3
        Test.@test quality_summary.sample == sample_names
        Test.@test all(quality_summary.n_reads .== 3)

        kmer_profiles = Mycelia.analyze_environmental_kmer_profiles(
            fastqs;
            sample_names = sample_names,
            k = 3,
            max_reads = 2
        )
        Test.@test size(kmer_profiles.distance_matrix) == (3, 3)
        Test.@test all(LinearAlgebra.diag(kmer_profiles.distance_matrix) .== 0)
        Test.@test LinearAlgebra.issymmetric(kmer_profiles.distance_matrix)
        Test.@test DataFrames.nrow(kmer_profiles.spectrum_summary) == 3
        Test.@test DataFrames.nrow(kmer_profiles.pcoa_df) == 3
        Test.@test Set(kmer_profiles.pcoa_df.sample) == Set(sample_names)
        Test.@test all(kmer_profiles.spectrum_summary.reads_used .== 2)
    end
end
