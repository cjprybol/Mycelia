import Test
import Mycelia

Test.@testset "Alignment Parsing Helpers" begin
    Test.@testset "dnadiff report parsing" begin
        mktempdir() do dir
            report_path = joinpath(dir, "sample.report")
            open(report_path, "w") do io
                println(io, "1-to-1")
                println(io, "TotalBases 1000 900")
                println(io, "AlignedBases 800 700 80.0 77.8")
                println(io, "UnalignedBases 200 200 20.0 22.2")
                println(io, "AvgIdentity 98.5")
                println(io, "AvgIdentity(Aligned) 99.0")
                println(io, "SNPs 5")
                println(io, "Indels 2 3")
            end

            parsed = Mycelia.parse_dnadiff_report(report_path)
            Test.@test parsed.summary.avg_identity == 98.5
            Test.@test parsed.summary.snps == 5
            Test.@test parsed.distance ≈ 0.015
            Test.@test parsed.distance_aligned ≈ 0.01
        end
    end

    Test.@testset "MUMmer coords parsing and summary" begin
        mktempdir() do dir
            coords_path = joinpath(dir, "sample.coords")
            open(coords_path, "w") do io
                println(io, "1 100 5 104 100 100 99.5 90.0 90.0 ref1 qry1")
            end

            coords = Mycelia.parse_mummer_coords_table(coords_path)
            Test.@test Mycelia.DataFrames.nrow(coords) == 1
            Test.@test coords[1, "ref_start"] == 1
            Test.@test coords[1, "identity"] == 99.5

            summary = Mycelia.summarize_mummer_coords(coords; reference_length = 1000, query_length = 900)
            Test.@test summary.num_alignments == 1
            Test.@test summary.aligned_bases_ref == 100
            Test.@test summary.distance ≈ 0.005
            Test.@test summary.aligned_pct_ref ≈ 10.0
        end
    end

    Test.@testset "Alignment scoring helpers" begin
        alignment = Mycelia.assess_alignment("ATCG", "ATCG")
        Test.@test alignment.total_matches == 4
        Test.@test alignment.total_edits == 0
        Test.@test Mycelia.assess_alignment_accuracy(alignment) == 1.0

        seq = Mycelia.BioSequences.LongDNA{4}("ATG")
        alignment_result, orientation = Mycelia.assess_optimal_kmer_alignment(seq, seq)
        Test.@test alignment_result.total_matches == 3
        Test.@test orientation == true

        rc_seq = Mycelia.BioSequences.reverse_complement(seq)
        Test.@test Mycelia.is_equivalent(seq, rc_seq)
    end

    Test.@testset "Minimap and label helpers" begin
        size = Mycelia.system_mem_to_minimap_index_size(system_mem_gb = 8.0, denominator = 4.0)
        Test.@test size == "2G"

        tag = Mycelia.build_sample_tag("foo bar.fastq.gz")
        Test.@test tag == "foo_bar"

        label = Mycelia.build_output_label(["/path/a.fastq.gz", "/path/b.fastq.gz"])
        Test.@test label == "a__b"
    end
end
