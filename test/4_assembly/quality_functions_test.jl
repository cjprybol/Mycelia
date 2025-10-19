# Quality Functions Test
#
# Tests for joint quality calculation and quality-aware operations.
#
# Run with: julia --project=. test/4_assembly/quality_functions_test.jl

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "Quality Functions" begin

    Test.@testset "combine_phred_scores - Basic" begin
        # Two Q10 observations combine to Q20
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[10, 10])
        Test.@test combined == UInt8(20)

        # Three Q20 observations combine to Q60
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[20, 20, 20])
        Test.@test combined == UInt8(60)

        # Single observation returns same score
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[30])
        Test.@test combined == UInt8(30)
    end

    Test.@testset "combine_phred_scores - Can Exceed 60" begin
        # Three Q30 observations combine to Q90
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[30, 30, 30])
        Test.@test combined == UInt8(90)

        # Ten Q20 observations combine to Q200
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[20, 20, 20, 20, 20, 20, 20, 20, 20, 20])
        Test.@test combined == UInt8(200)
    end

    Test.@testset "combine_phred_scores - Clamps at 255" begin
        # Many high-quality observations clamp at 255
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[60, 60, 60, 60, 60])
        Test.@test combined == UInt8(255)
    end

    Test.@testset "combine_phred_scores - Empty" begin
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[])
        Test.@test combined == UInt8(0)
    end

    Test.@testset "get_vertex_joint_quality - Single Observation" begin
        seq = "ATGC"
        qual = [30, 30, 30, 30]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record], 3; dataset_id="test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        joint_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "test")

        # Single observation: joint quality = original quality
        Test.@test joint_qual == UInt8[30, 30, 30]
    end

    Test.@testset "get_vertex_joint_quality - Multiple Observations" begin
        # Two observations with Q10
        seq = "ATGC"
        qual1 = [10, 10, 10, 10]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [10, 10, 10, 10]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id="test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        joint_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "test")

        # Two Q10 observations combine to Q20
        Test.@test joint_qual == UInt8[20, 20, 20]
    end

    Test.@testset "get_vertex_joint_quality - Different Quality Profiles" begin
        # Two observations with different qualities
        seq = "ATGC"
        qual1 = [20, 30, 25, 15]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [30, 20, 25, 35]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id="test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        joint_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "test")

        # Position-wise combination: [20+30, 30+20, 25+25]
        Test.@test joint_qual == UInt8[50, 50, 50]
    end

    Test.@testset "get_vertex_mean_quality" begin
        seq = "ATGC"
        qual1 = [20, 30, 40, 50]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [30, 30, 30, 30]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id="test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        mean_qual = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, "test")

        # Mean of [20,30] [30,30] [40,30]
        Test.@test mean_qual ≈ [25.0, 30.0, 35.0]
    end

    Test.@testset "get_vertex_min_quality" begin
        seq = "ATGC"
        qual1 = [20, 30, 40, 50]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [30, 20, 30, 40]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id="test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        min_qual = Mycelia.Rhizomorph.get_vertex_min_quality(vertex_data, "test")

        # Min of [20,30] [30,20] [40,30]
        Test.@test min_qual == UInt8[20, 20, 30]
    end

    Test.@testset "filter_by_quality" begin
        # Create graph with k-mers of varying quality
        seq1 = "ATGC"  # High quality
        qual1 = [40, 40, 40, 40]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq1, qual_str1)

        seq2 = "CGTA"  # Low quality
        qual2 = [10, 10, 10, 10]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq2, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id="test")

        # Filter for Q30+
        high_qual = Mycelia.Rhizomorph.filter_by_quality(graph, 30, "test")

        # Should include ATG and TGC (from high-quality read)
        Test.@test Kmers.DNAKmer{3}("ATG") in high_qual
        Test.@test Kmers.DNAKmer{3}("TGC") in high_qual

        # Should NOT include CGT and GTA (from low-quality read)
        Test.@test !(Kmers.DNAKmer{3}("CGT") in high_qual)
        Test.@test !(Kmers.DNAKmer{3}("GTA") in high_qual)
    end

    Test.@testset "Quality Functions - Nonexistent Dataset" begin
        seq = "ATGC"
        qual = [30, 30, 30, 30]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record], 3; dataset_id="test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        # Query nonexistent dataset
        Test.@test isnothing(Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "nonexistent"))
        Test.@test isnothing(Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, "nonexistent"))
        Test.@test isnothing(Mycelia.Rhizomorph.get_vertex_min_quality(vertex_data, "nonexistent"))
    end
end

println("✓ Quality functions tests completed")
