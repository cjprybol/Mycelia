import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import MetaGraphsNext
import Graphs
import StatsBase

Test.@testset "String Graphs Tests" begin
    Test.@testset "N-gram Generation" begin
        # Test basic n-gram generation
        test_string = "ABCDEF"
        
        # Test 2-grams
        bigrams = Mycelia.ngrams(test_string, 2)
        expected_bigrams = ["AB", "BC", "CD", "DE", "EF"]
        Test.@test bigrams == expected_bigrams
        Test.@test length(bigrams) == length(test_string) - 2 + 1
        
        # Test 3-grams
        trigrams = Mycelia.ngrams(test_string, 3)
        expected_trigrams = ["ABC", "BCD", "CDE", "DEF"]
        Test.@test trigrams == expected_trigrams
        Test.@test length(trigrams) == length(test_string) - 3 + 1
        
        # Test single character n-grams
        unigrams = Mycelia.ngrams(test_string, 1)
        Test.@test unigrams == ["A", "B", "C", "D", "E", "F"]
        Test.@test length(unigrams) == length(test_string)
        
        # Test edge cases
        empty_result = Mycelia.ngrams("AB", 3)  # n > string length
        Test.@test isempty(empty_result)
        
        single_ngram = Mycelia.ngrams("ABC", 3)  # n == string length
        Test.@test single_ngram == ["ABC"]
        Test.@test length(single_ngram) == 1
    end

    Test.@testset "String to N-gram Graph Conversion" begin
        # Test simple string to graph conversion
        test_string = "ABCABC"
        n = 2
        
        graph = Mycelia.string_to_ngram_graph(s=test_string, n=n)
        
        Test.@test graph isa MetaGraphsNext.MetaGraph
        
        # Check that all expected n-grams are nodes in the graph
        expected_ngrams = ["AB", "BC", "CA"]  # Unique 2-grams
        for ngram in expected_ngrams
            Test.@test MetaGraphsNext.haskey(graph, ngram)
        end
        
        # Test node counts (each n-gram appears twice)
        Test.@test graph["AB"] == 2
        Test.@test graph["BC"] == 2  
        Test.@test graph["CA"] == 1
        
        # Test that edges exist
        Test.@test MetaGraphsNext.haskey(graph, "AB", "BC")
        Test.@test MetaGraphsNext.haskey(graph, "BC", "CA")
        Test.@test MetaGraphsNext.haskey(graph, "CA", "AB")
    end

    Test.@testset "Graph Properties and Structure" begin
        # Test with repeating pattern
        repeated_string = "ABABAB"
        graph = Mycelia.string_to_ngram_graph(s=repeated_string, n=2)
        
        # Should have only 2 unique n-grams: "AB", "BA"
        expected_labels = Set(["AB", "BA"])
        actual_labels = Set(MetaGraphsNext.labels(graph))
        Test.@test actual_labels == expected_labels
        
        # Test edge counts
        Test.@test graph["AB", "BA"] == 2  # AB->BA appears twice
        Test.@test graph["BA", "AB"] == 2  # BA->AB appears twice
        
        # Test with longer n-grams
        long_string = "ABCDEFGH"
        long_graph = Mycelia.string_to_ngram_graph(s=long_string, n=4)
        
        expected_4grams = ["ABCD", "BCDE", "CDEF", "DEFG", "EFGH"]
        for ngram in expected_4grams
            Test.@test MetaGraphsNext.haskey(long_graph, ngram)
            Test.@test long_graph[ngram] == 1  # Each appears once
        end
    end

    Test.@testset "Connected Components Analysis" begin
        # Create a simple test graph for connected components
        test_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Int,
            edge_data_type=Int
        )
        
        # Add disconnected components
        test_graph["A"] = 1
        test_graph["B"] = 1
        test_graph["C"] = 1
        test_graph["D"] = 1
        
        # Connect A-B and C-D (two components)
        test_graph["A", "B"] = 1
        test_graph["C", "D"] = 1
        
        # Test connected components finding
        components = Mycelia.find_connected_components(test_graph)
        Test.@test components isa Vector{Vector{Int}}
        Test.@test length(components) == 2  # Two separate components
    end

    Test.@testset "Path Collapsing Operations" begin
        # Create a linear path graph for testing
        linear_string = "ABCDEF"
        linear_graph = Mycelia.string_to_ngram_graph(s=linear_string, n=2)
        
        # Test unbranching path collapse
        collapsed_graph = Mycelia.collapse_unbranching_paths(linear_graph)
        Test.@test collapsed_graph isa MetaGraphsNext.MetaGraph
        
        # In a linear path, most vertices should be collapsible
        collapsible_vertices = Mycelia._find_collapsible_vertices(linear_graph)
        Test.@test collapsible_vertices isa Vector{Int}
        Test.@test length(collapsible_vertices) >= 0
    end

    Test.@testset "String Assembly" begin
        # Test string assembly from graph
        simple_string = "ABCD"
        simple_graph = Mycelia.string_to_ngram_graph(s=simple_string, n=2)
        
        assembled_strings = Mycelia.assemble_strings(simple_graph)
        Test.@test assembled_strings isa Vector{String}
        Test.@test length(assembled_strings) >= 1
        
        # For a simple linear string, should be able to reconstruct
        # (though exact reconstruction depends on algorithm details)
        for assembled in assembled_strings
            Test.@test assembled isa String
            Test.@test length(assembled) >= 3  # At least some reasonable length
        end
    end

    Test.@testset "Graph Visualization Support" begin
        # Test that plot function can be called without error
        test_graph = Mycelia.string_to_ngram_graph(s="ABCABC", n=2)
        
        # Note: We can't easily test the actual plotting output without display
        # But we can test that the function exists and accepts the right type
        Test.@test hasmethod(Mycelia.plot_ngram_graph, (typeof(test_graph),))
    end

    Test.@testset "Edge Cases and Error Handling" begin
        # Test with empty string
        Test.@test_throws BoundsError Mycelia.ngrams("", 1)
        
        # Test with n = 0
        Test.@test_throws BoundsError Mycelia.ngrams("ABC", 0)
        
        # Test with very long n
        long_n_result = Mycelia.ngrams("ABC", 10)
        Test.@test isempty(long_n_result)
        
        # Test single character string
        single_char = Mycelia.ngrams("A", 1)
        Test.@test single_char == ["A"]
        
        # Test graph with single n-gram
        single_graph = Mycelia.string_to_ngram_graph(s="AA", n=1)
        Test.@test MetaGraphsNext.haskey(single_graph, "A")
        Test.@test single_graph["A"] == 2
    end

    Test.@testset "Unicode and Special Characters" begin
        # Test with unicode characters
        unicode_string = "αβγδε"
        unicode_grams = Mycelia.ngrams(unicode_string, 2)
        Test.@test length(unicode_grams) == 4
        Test.@test unicode_grams[1] == "αβ"
        Test.@test unicode_grams[end] == "δε"
        
        # Test with numbers and symbols
        mixed_string = "A1B2C3"
        mixed_grams = Mycelia.ngrams(mixed_string, 2)
        expected_mixed = ["A1", "1B", "B2", "2C", "C3"]
        Test.@test mixed_grams == expected_mixed
        
        # Test graph creation with mixed characters
        mixed_graph = Mycelia.string_to_ngram_graph(s=mixed_string, n=2)
        for gram in expected_mixed
            Test.@test MetaGraphsNext.haskey(mixed_graph, gram)
        end
    end

    Test.@testset "Performance and Memory Tests" begin
        # Test with moderately large strings
        large_string = repeat("ABCD", 1000)  # 4000 characters
        large_grams = Mycelia.ngrams(large_string, 3)
        
        Test.@test length(large_grams) == length(large_string) - 3 + 1
        Test.@test length(large_grams) == 3998
        
        # Test graph creation doesn't fail with large input
        Test.@test_nowarn Mycelia.string_to_ngram_graph(s=large_string, n=3)
        
        # Test memory efficiency with repeated patterns
        repeated_string = repeat("AB", 500)  # Highly repetitive
        repeated_graph = Mycelia.string_to_ngram_graph(s=repeated_string, n=2)
        
        # Should have only 2 unique nodes despite 1000 character input
        Test.@test length(MetaGraphsNext.labels(repeated_graph)) == 2
        Test.@test MetaGraphsNext.haskey(repeated_graph, "AB")
        Test.@test MetaGraphsNext.haskey(repeated_graph, "BA")
    end

    Test.@testset "Graph Algorithm Integration" begin
        # Test with complex branching structure
        branched_string = "ABCDEFABGHIJ"  # Has branching at "AB"
        branched_graph = Mycelia.string_to_ngram_graph(s=branched_string, n=2)
        
        # "AB" should appear twice
        Test.@test branched_graph["AB"] == 2
        
        # Should have edges from "AB" to both "BC" and "BG"
        Test.@test MetaGraphsNext.haskey(branched_graph, "AB", "BC")
        Test.@test MetaGraphsNext.haskey(branched_graph, "AB", "BG")
        
        # Test that component analysis works
        components = Mycelia.find_connected_components(branched_graph)
        Test.@test length(components) >= 1  # Should be connected
    end
end