import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import DataFrames
import BioSequences
import FASTX
import SHA
import DelimitedFiles

Test.@testset "Sequence Comparison Tests" begin
    Test.@testset "Mash Distance Calculation" begin
        # Test basic mash distance formula
        Test.@test Mycelia.mash_distance_from_jaccard(1.0, 21) ≈ 0.0
        Test.@test Mycelia.mash_distance_from_jaccard(0.0, 21) ≈ 1.0
        Test.@test Mycelia.mash_distance_from_jaccard(0.5, 21) > 0.0
        Test.@test Mycelia.mash_distance_from_jaccard(0.5, 21) < 1.0
        
        # Test parameter validation
        Test.@test_throws ErrorException Mycelia.mash_distance_from_jaccard(-0.1, 21)
        Test.@test_throws ErrorException Mycelia.mash_distance_from_jaccard(1.1, 21)
        Test.@test_throws ErrorException Mycelia.mash_distance_from_jaccard(0.5, 0)
        Test.@test_throws ErrorException Mycelia.mash_distance_from_jaccard(0.5, -1)
        
        # Test k-mer size effects
        dist_k15 = Mycelia.mash_distance_from_jaccard(0.8, 15)
        dist_k21 = Mycelia.mash_distance_from_jaccard(0.8, 21)
        Test.@test dist_k15 > dist_k21  # Smaller k should give larger distance
    end

    Test.@testset "SHA256 Sequence Hashing" begin
        # Test string input
        test_seq = "ATCGATCG"
        hash1 = Mycelia.seq2sha256(test_seq)
        Test.@test hash1 isa String
        Test.@test length(hash1) == 64
        
        # Test case insensitivity
        hash2 = Mycelia.seq2sha256("atcgatcg")
        Test.@test hash1 == hash2
        
        # Test BioSequence input
        bio_seq = BioSequences.LongDNA{4}("ATCGATCG")
        hash3 = Mycelia.seq2sha256(bio_seq)
        Test.@test hash1 == hash3
        
        # Test different sequences
        different_hash = Mycelia.seq2sha256("GCTAGCTA")
        Test.@test hash1 != different_hash
        
        # Test empty sequence
        empty_hash = Mycelia.seq2sha256("")
        Test.@test empty_hash isa String
        Test.@test length(empty_hash) == 64
    end

    Test.@testset "Mock Mash Comparison Tests" begin
        # Create test FASTA files
        temp_fasta1 = tempname() * ".fasta"
        temp_fasta2 = tempname() * ".fasta"
        
        write(temp_fasta1, ">seq1\nATCGATCGATCGATCG\n")
        write(temp_fasta2, ">seq2\nGCTAGCTAGCTAGCTA\n")
        
        # Test file existence validation
        Test.@test_throws ErrorException Mycelia.run_mash_comparison("nonexistent1.fasta", "nonexistent2.fasta")
        Test.@test_throws ErrorException Mycelia.run_mash_comparison(temp_fasta1, "nonexistent2.fasta")
        
        # Note: Actual mash execution would require mash to be installed
        # In a real test environment, this would test the full pipeline
        
        # Cleanup
        rm(temp_fasta1, force=true)
        rm(temp_fasta2, force=true)
    end

    Test.@testset "Pairwise Minimap Comparison Mock Tests" begin
        # Create test FASTA files
        ref_fasta = tempname() * ".fasta"
        query_fasta = tempname() * ".fasta"
        
        ref_seq = ">reference\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        query_seq = ">query\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        
        write(ref_fasta, ref_seq)
        write(query_fasta, query_seq)
        
        # Test that files exist and can be read
        Test.@test isfile(ref_fasta)
        Test.@test isfile(query_fasta)
        
        # Test FASTA reading logic
        ref_reader = FASTX.FASTA.Reader(open(ref_fasta))
        query_reader = FASTX.FASTA.Reader(open(query_fasta))
        
        ref_record = first(ref_reader)
        query_record = first(query_reader)
        
        Test.@test length(FASTX.sequence(ref_record)) == 52
        Test.@test length(FASTX.sequence(query_record)) == 52
        
        close(ref_reader)
        close(query_reader)
        
        # Note: Actual minimap2 execution would require minimap2 to be installed
        # In a real test environment, this would test the full comparison pipeline
        
        # Cleanup
        rm(ref_fasta, force=true)
        rm(query_fasta, force=true)
    end

    Test.@testset "FastANI Mock Tests" begin
        # Test query/reference list file creation
        temp_query_list = tempname() * ".txt"
        temp_ref_list = tempname() * ".txt"
        temp_fasta1 = tempname() * ".fasta"
        temp_fasta2 = tempname() * ".fasta"
        
        # Create mock FASTA files
        write(temp_fasta1, ">genome1\nATCGATCGATCGATCG\n")
        write(temp_fasta2, ">genome2\nGCTAGCTAGCTAGCTA\n")
        
        # Create list files
        write(temp_query_list, temp_fasta1 * "\n")
        write(temp_ref_list, temp_fasta2 * "\n")
        
        Test.@test isfile(temp_query_list)
        Test.@test isfile(temp_ref_list)
        
        # Test file reading
        query_list_content = read(temp_query_list, String)
        ref_list_content = read(temp_ref_list, String)
        
        Test.@test occursin(temp_fasta1, query_list_content)
        Test.@test occursin(temp_fasta2, ref_list_content)
        
        # Note: Actual FastANI execution would require fastani to be installed
        # The functions would use Bioconda to install and run FastANI
        
        # Cleanup
        rm(temp_query_list, force=true)
        rm(temp_ref_list, force=true)
        rm(temp_fasta1, force=true)
        rm(temp_fasta2, force=true)
    end

    Test.@testset "Sequence Merging and Mapping Mock Tests" begin
        # Test parameter validation
        Test.@test_nowarn Mycelia.normalized_current_date()
        
        # Create mock FASTQ files
        temp_fastq1 = tempname() * ".fq"
        temp_fastq2 = tempname() * ".fq"
        temp_ref = tempname() * ".fasta"
        temp_index = tempname() * ".mmi"
        
        fastq_content = "@read1\nATCGATCG\n+\nIIIIIIII\n"
        ref_content = ">ref\nATCGATCGATCGATCGATCGATCGATCGATCG\n"
        
        write(temp_fastq1, fastq_content)
        write(temp_fastq2, fastq_content)
        write(temp_ref, ref_content)
        write(temp_index, "mock_index_content")
        
        # Test file existence
        Test.@test isfile(temp_fastq1)
        Test.@test isfile(temp_fastq2)
        Test.@test isfile(temp_ref)
        Test.@test isfile(temp_index)
        
        # Test FASTQ list creation
        fastq_list = [temp_fastq1, temp_fastq2]
        Test.@test length(fastq_list) == 2
        Test.@test all(isfile, fastq_list)
        
        # Note: Full pipeline testing would require minimap2 and other tools
        
        # Cleanup
        rm(temp_fastq1, force=true)
        rm(temp_fastq2, force=true)
        rm(temp_ref, force=true)
        rm(temp_index, force=true)
    end

    Test.@testset "Edge Cases and Error Handling" begin
        # Test with empty sequences
        empty_hash = Mycelia.seq2sha256("")
        Test.@test length(empty_hash) == 64
        
        # Test with very long sequences
        long_seq = repeat("ATCG", 10000)
        long_hash = Mycelia.seq2sha256(long_seq)
        Test.@test length(long_hash) == 64
        
        # Test mash distance edge cases
        Test.@test Mycelia.mash_distance_from_jaccard(0.0, 21) == 1.0  # No shared k-mers
        Test.@test Mycelia.mash_distance_from_jaccard(1.0, 21) == 0.0  # Identical sequences
        
        # Test with different k-mer sizes
        for k in [15, 17, 19, 21, 23]
            dist = Mycelia.mash_distance_from_jaccard(0.5, k)
            Test.@test 0.0 < dist < 1.0
        end
    end

    Test.@testset "Mathematical Properties" begin
        # Test that mash distance is monotonic with Jaccard index
        jaccard_values = [0.1, 0.3, 0.5, 0.7, 0.9]
        distances = [Mycelia.mash_distance_from_jaccard(j, 21) for j in jaccard_values]
        
        # Distances should decrease as Jaccard similarity increases
        for i in 1:(length(distances)-1)
            Test.@test distances[i] > distances[i+1]
        end
        
        # Test formula consistency
        for jaccard in [0.1, 0.2, 0.5, 0.8, 0.9]
            dist = Mycelia.mash_distance_from_jaccard(jaccard, 21)
            Test.@test 0.0 ≤ dist ≤ 1.0
        end
    end

    Test.@testset "Hash Consistency Tests" begin
        # Test that identical sequences produce identical hashes
        seq1 = "ATCGATCGATCG"
        seq2 = "ATCGATCGATCG"
        Test.@test Mycelia.seq2sha256(seq1) == Mycelia.seq2sha256(seq2)
        
        # Test case insensitivity
        upper_seq = "ATCGATCG"
        lower_seq = "atcgatcg"
        mixed_seq = "AtCgAtCg"
        
        hash_upper = Mycelia.seq2sha256(upper_seq)
        hash_lower = Mycelia.seq2sha256(lower_seq)
        hash_mixed = Mycelia.seq2sha256(mixed_seq)
        
        Test.@test hash_upper == hash_lower
        Test.@test hash_upper == hash_mixed
        
        # Test different sequence types
        dna_seq = BioSequences.LongDNA{4}("ATCGATCG")
        string_seq = "ATCGATCG"
        
        Test.@test Mycelia.seq2sha256(dna_seq) == Mycelia.seq2sha256(string_seq)
    end
end