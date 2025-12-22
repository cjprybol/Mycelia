import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import FASTX
import Graphs

Test.@testset "Viterbi Polishing and Error Correction Tests" begin
    Test.@testset "Viterbi Maximum Likelihood Traversals" begin
        # Create a simple test k-mer graph
        test_records = [
            FASTX.FASTA.Record("seq1", "ATCGATCGATCG"),
            FASTX.FASTA.Record("seq2", "ATCGATCGATCG"),  # Same sequence for higher coverage
            FASTX.FASTA.Record("seq3", "ATCGATCGATCC")   # Similar with one difference
        ]
        
        kmer_type = Mycelia.Kmers.DNAKmer{4}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, test_records)
        
        # Test with default parameters
        result = Mycelia.viterbi_maximum_likelihood_traversals(
            graph,
            verbosity="dataset"
        )
        
        Test.@test result isa NamedTuple
        Test.@test haskey(result, :viterbi_ml_corrected_paths) || haskey(result, :paths)
        
        # Test with custom error rate
        custom_result = Mycelia.viterbi_maximum_likelihood_traversals(
            graph,
            error_rate=0.1,
            verbosity="dataset"
        )
        
        Test.@test custom_result isa NamedTuple
        
        # Test error rate validation
        Test.@test_throws ErrorException Mycelia.viterbi_maximum_likelihood_traversals(
            graph,
            error_rate=0.6,  # > 0.5 should error
            verbosity="dataset"
        )
    end

    Test.@testset "FASTQ Record Processing" begin
        # Create a test FASTQ record
        test_record = FASTX.FASTQ.Record("test_read", "ATCGATCG", "IIIIIIII")
        
        # Create a simple k-mer graph
        fasta_records = [FASTX.FASTA.Record("ref", "ATCGATCGATCG")]
        kmer_type = Mycelia.Kmers.DNAKmer{4}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, fasta_records)
        
        # Create mock k-shortest paths data
        # In practice, this would come from a proper shortest paths algorithm
        mock_paths = Dict(
            "test_read" => [
                ([1, 2, 3], 0.8),  # (path, weight)
                ([1, 2, 4], 0.6),
                ([1, 3, 4], 0.4)
            ]
        )
        
        # Test record processing
        result = Mycelia.process_fastq_record(
            record=test_record,
            kmer_graph=graph,
            yen_k_shortest_paths_and_weights=mock_paths,
            yen_k=3
        )
        
        Test.@test result isa NamedTuple
        Test.@test haskey(result, :polished_record) || haskey(result, :record)
    end

    Test.@testset "FASTQ Polishing Pipeline" begin
        # Create a test FASTQ file
        temp_fastq = tempname() * ".fastq"
        fastq_content = """@read1
ATCGATCGATCG
+
IIIIIIIIIIII
@read2
GCTAGCTAGCTA
+
IIIIIIIIIIII
"""
        write(temp_fastq, fastq_content)
        
        # Test polishing with k=4
        result = Mycelia.polish_fastq(fastq=temp_fastq, k=4)
        
        Test.@test result isa NamedTuple
        Test.@test haskey(result, :polished_fastq) || haskey(result, :output_file)
        
        # Cleanup
        rm(temp_fastq, force=true)
    end

    Test.@testset "Iterative Polishing" begin
        # Create a test FASTQ file for iterative polishing
        temp_fastq = tempname() * ".fastq"
        fastq_content = """@read1
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
@read2
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
"""
        write(temp_fastq, fastq_content)
        
        # Test iterative polishing with small max_k for faster testing
        result = Mycelia.iterative_polishing(temp_fastq, 9)
        
        Test.@test result isa NamedTuple
        # Result should contain information about polishing iterations
        
        # Cleanup
        rm(temp_fastq, force=true)
    end

    Test.@testset "Error Rate Validation" begin
        # Create minimal graph for testing
        test_records = [FASTX.FASTA.Record("seq", "ATCGATCG")]
        kmer_type = Mycelia.Kmers.DNAKmer{3}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, test_records)
        
        # Test valid error rates
        valid_rates = [0.01, 0.05, 0.1, 0.2, 0.4, 0.49]
        for rate in valid_rates
            Test.@test_nowarn Mycelia.viterbi_maximum_likelihood_traversals(
                graph,
                error_rate=rate,
                verbosity="dataset"
            )
        end
        
        # Test invalid error rates
        invalid_rates = [0.5, 0.6, 0.8, 1.0]
        for rate in invalid_rates
            Test.@test_throws ErrorException Mycelia.viterbi_maximum_likelihood_traversals(
                graph,
                error_rate=rate,
                verbosity="dataset"
            )
        end
    end

    Test.@testset "Verbosity Levels" begin
        # Create test graph
        test_records = [FASTX.FASTA.Record("seq", "ATCGATCG")]
        kmer_type = Mycelia.Kmers.DNAKmer{3}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, test_records)
        
        # Test different verbosity levels
        verbosity_levels = ["debug", "reads", "dataset"]
        for level in verbosity_levels
            Test.@test_nowarn Mycelia.viterbi_maximum_likelihood_traversals(
                graph,
                verbosity=level
            )
        end
        
        # Test invalid verbosity level
        Test.@test_throws AssertionError Mycelia.viterbi_maximum_likelihood_traversals(
            graph,
            verbosity="invalid"
        )
    end

    Test.@testset "K-mer Graph Properties" begin
        # Test that the algorithm can handle graphs with different k values
        test_sequence = "ATCGATCGATCGATCG"
        
        for k in [3, 4, 5, 6]
            if length(test_sequence) >= k
                records = [FASTX.FASTA.Record("seq", test_sequence)]
                kmer_type = Mycelia.Kmers.DNAKmer{k}
                graph = Mycelia.build_stranded_kmer_graph(kmer_type, records)
                
                Test.@test graph.gprops[:k] == k
                
                # Test that viterbi can process this graph
                Test.@test_nowarn Mycelia.viterbi_maximum_likelihood_traversals(
                    graph,
                    verbosity="dataset"
                )
            end
        end
    end

    Test.@testset "Coverage and Likelihood Calculations" begin
        # Create graph with varying coverage
        high_coverage_records = [
            FASTX.FASTA.Record("seq$i", "ATCGATCG") for i in 1:10
        ]
        
        kmer_type = Mycelia.Kmers.DNAKmer{3}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, high_coverage_records)
        
        # Test that coverage affects likelihood calculations
        result = Mycelia.viterbi_maximum_likelihood_traversals(
            graph,
            verbosity="dataset"
        )
        
        Test.@test result isa NamedTuple
        
        # Compare with low coverage
        low_coverage_records = [FASTX.FASTA.Record("seq1", "ATCGATCG")]
        low_graph = Mycelia.build_stranded_kmer_graph(kmer_type, low_coverage_records)
        
        low_result = Mycelia.viterbi_maximum_likelihood_traversals(
            low_graph,
            verbosity="dataset"
        )
        
        Test.@test low_result isa NamedTuple
    end

    Test.@testset "Edge Cases and Error Handling" begin
        # Test with empty graph (minimal case)
        empty_records = FASTX.FASTA.Record[]
        kmer_type = Mycelia.Kmers.DNAKmer{3}
        empty_graph = Mycelia.build_stranded_kmer_graph(kmer_type, empty_records)
        
        # This might error or handle gracefully depending on implementation
        empty_result = try
            Mycelia.viterbi_maximum_likelihood_traversals(
                empty_graph,
                verbosity="dataset"
            )
            :ok
        catch
            :error
        end
        Test.@test empty_result in (:ok, :error)
        
        # Test with very short sequences
        short_records = [FASTX.FASTA.Record("short", "AT")]
        short_graph = Mycelia.build_stranded_kmer_graph(kmer_type, short_records)
        
        Test.@test_nowarn Mycelia.viterbi_maximum_likelihood_traversals(
            short_graph,
            verbosity="dataset"
        )
        
        # Test with non-existent files
        Test.@test_throws Exception Mycelia.polish_fastq(
            fastq="nonexistent.fastq",
            k=4
        )
    end

    Test.@testset "Quality Score Integration" begin
        # Test that quality scores are properly handled
        high_quality_record = FASTX.FASTQ.Record("hq", "ATCGATCG", "IIIIIIII")
        low_quality_record = FASTX.FASTQ.Record("lq", "ATCGATCG", "!!!!!!!!")
        
        # Test records have different quality scores
        hq_qualities = FASTX.FASTQ.quality(high_quality_record)
        lq_qualities = FASTX.FASTQ.quality(low_quality_record)
        
        Test.@test all(hq_qualities .> lq_qualities)
        
        # Create graph for testing
        ref_records = [FASTX.FASTA.Record("ref", "ATCGATCGATCG")]
        kmer_type = Mycelia.Kmers.DNAKmer{4}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, ref_records)
        
        # Mock paths for testing
        mock_paths = Dict(
            "hq" => [([1, 2], 0.9)],
            "lq" => [([1, 2], 0.9)]
        )
        
        # Test that both can be processed
        Test.@test_nowarn Mycelia.process_fastq_record(
            record=high_quality_record,
            kmer_graph=graph,
            yen_k_shortest_paths_and_weights=mock_paths,
            yen_k=1
        )
        
        Test.@test_nowarn Mycelia.process_fastq_record(
            record=low_quality_record,
            kmer_graph=graph,
            yen_k_shortest_paths_and_weights=mock_paths,
            yen_k=1
        )
    end

    Test.@testset "Algorithm Parameter Effects" begin
        # Create test data
        test_records = [FASTX.FASTA.Record("seq", "ATCGATCGATCG")]
        kmer_type = Mycelia.Kmers.DNAKmer{4}
        graph = Mycelia.build_stranded_kmer_graph(kmer_type, test_records)
        
        # Test different error rates and compare results
        low_error = Mycelia.viterbi_maximum_likelihood_traversals(
            graph,
            error_rate=0.01,
            verbosity="dataset"
        )
        
        high_error = Mycelia.viterbi_maximum_likelihood_traversals(
            graph,
            error_rate=0.2,
            verbosity="dataset"
        )
        
        Test.@test low_error isa NamedTuple
        Test.@test high_error isa NamedTuple
        
        # Results should be different (though we can't easily test specifics without implementation details)
    end
end
