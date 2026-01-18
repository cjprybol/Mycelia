# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/in_development/viterbi_next.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/in_development/viterbi_next.jl", "test/in_development", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX

Test.@testset "Viterbi Next-Generation Algorithm Tests" begin
    
    # Helper function to create test graph
    function create_viterbi_test_graph()
        graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph(),
            label_type=String,
            vertex_data_type=Mycelia.KmerVertexData,
            edge_data_type=Mycelia.KmerEdgeData,
            weight_function=Mycelia.Rhizomorph.edge_data_weight,
            default_weight=0.0
        )
        
        # Add vertices for 3-mers: ATC -> TCG -> CGA -> GAT
        for kmer in ["ATC", "TCG", "CGA", "GAT"]
            graph[kmer] = Mycelia.KmerVertexData(kmer)
        end
        
        # Add edges forming a path
        # High coverage edge (3 observations)
        high_coverage = [
            ((1, 1, Mycelia.Rhizomorph.Forward), (1, 2, Mycelia.Rhizomorph.Forward)),
            ((2, 1, Mycelia.Rhizomorph.Forward), (2, 2, Mycelia.Rhizomorph.Forward)),
            ((3, 1, Mycelia.Rhizomorph.Forward), (3, 2, Mycelia.Rhizomorph.Forward))
        ]
        graph["ATC", "TCG"] = Mycelia.KmerEdgeData(
            high_coverage, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward
        )
        # Medium coverage edge (2 observations)
        medium_coverage = [
            ((1, 2, Mycelia.Rhizomorph.Forward), (1, 3, Mycelia.Rhizomorph.Forward)),
            ((2, 2, Mycelia.Rhizomorph.Forward), (2, 3, Mycelia.Rhizomorph.Forward))
        ]
        graph["TCG", "CGA"] = Mycelia.KmerEdgeData(
            medium_coverage, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward
        )
        # Low coverage edge (1 observation)
        low_coverage = [
            ((1, 3, Mycelia.Rhizomorph.Forward), (1, 4, Mycelia.Rhizomorph.Forward))
        ]
        graph["CGA", "GAT"] = Mycelia.KmerEdgeData(
            low_coverage, Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Forward
        )
        
        return graph
    end
    
    Test.@testset "ViterbiState Construction" begin
        state = Mycelia.ViterbiState("ATC", Mycelia.Rhizomorph.Forward, 0.95, 1)
        Test.@test state.vertex_label == "ATC"
        Test.@test state.strand == Mycelia.Rhizomorph.Forward
        Test.@test state.emission_prob == 0.95
        Test.@test state.position == 1
        
        # Test invalid probability
        Test.@test_throws AssertionError Mycelia.ViterbiState("ATC", Mycelia.Rhizomorph.Forward, 1.5, 1)
        Test.@test_throws AssertionError Mycelia.ViterbiState("ATC", Mycelia.Rhizomorph.Forward, -0.1, 1)
        
        # Test invalid position
        Test.@test_throws AssertionError Mycelia.ViterbiState("ATC", Mycelia.Rhizomorph.Forward, 0.95, -1)
    end
    
    Test.@testset "ViterbiConfig Construction and Validation" begin
        # Default config
        config = Mycelia.ViterbiConfig()
        Test.@test config.match_prob == 0.95
        Test.@test config.mismatch_prob == 0.04
        Test.@test config.insertion_prob == 0.005
        Test.@test config.deletion_prob == 0.005
        Test.@test config.use_log_space == true
        
        # Custom config
        custom_config = Mycelia.ViterbiConfig(
            match_prob=0.9,
            mismatch_prob=0.08,
            insertion_prob=0.01,
            deletion_prob=0.01,
            batch_size=500
        )
        Test.@test custom_config.match_prob == 0.9
        Test.@test custom_config.batch_size == 500
        
        # Test invalid probabilities (should not sum to 1)
        Test.@test_throws AssertionError Mycelia.ViterbiConfig(
            match_prob=0.9, mismatch_prob=0.2, insertion_prob=0.005, deletion_prob=0.005
        )
    end
    
    Test.@testset "HMM Creation from Graph" begin
        graph = create_viterbi_test_graph()
        config = Mycelia.ViterbiConfig()
        
        states, transitions, emissions = Mycelia.create_hmm_from_graph(graph, config)
        
        # Should have states for each vertex (and reverse complements if enabled)
        expected_n_states = 4 * (config.consider_reverse_complement ? 2 : 1)
        Test.@test length(states) == expected_n_states
        Test.@test size(transitions) == (expected_n_states, expected_n_states)
        Test.@test length(emissions) == expected_n_states
        
        # Check that transition matrix rows sum to 1 (approximately)
        for i in 1:size(transitions, 1)
            row_sum = sum(transitions[i, :])
            if row_sum > 0  # Only check non-zero rows
                Test.@test abs(row_sum - 1.0) < 1e-10
            end
        end
        
        # States should have correct structure
        for state in states
            Test.@test state isa Mycelia.ViterbiState
            Test.@test state.strand in [Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse]
            Test.@test !isempty(state.vertex_label)
        end
    end
    
    Test.@testset "Emission Probability Calculation" begin
        config = Mycelia.ViterbiConfig()
        state = Mycelia.ViterbiState("ATC", Mycelia.Rhizomorph.Forward, 1.0, 1)
        
        # Perfect match
        match_prob = Mycelia.calculate_emission_probability(state, "ATC", config)
        Test.@test match_prob == config.match_prob
        
        # Single mismatch
        mismatch_prob = Mycelia.calculate_emission_probability(state, "ATG", config)
        Test.@test mismatch_prob == config.mismatch_prob
        
        # Multiple mismatches (should be lower probability)
        multi_mismatch_prob = Mycelia.calculate_emission_probability(state, "GGG", config)
        Test.@test multi_mismatch_prob < config.mismatch_prob
        Test.@test multi_mismatch_prob > 0
        
        # Test reverse complement state
        reverse_state = Mycelia.ViterbiState("ATC", Mycelia.Rhizomorph.Reverse, 1.0, 1)
        # GAT is reverse complement of ATC
        rc_match_prob = Mycelia.calculate_emission_probability(reverse_state, "GAT", config)
        Test.@test rc_match_prob == config.match_prob
    end
    
    Test.@testset "Simple Edit Distance" begin
        Test.@test Mycelia.simple_edit_distance("ATC", "ATC") == 0
        Test.@test Mycelia.simple_edit_distance("ATC", "ATG") == 1
        Test.@test Mycelia.simple_edit_distance("ATC", "GGG") == 3
        Test.@test Mycelia.simple_edit_distance("ATC", "ATCG") == 1  # Length difference
        Test.@test Mycelia.simple_edit_distance("ATCG", "ATC") == 1
    end
    
    Test.@testset "Viterbi Decoding" begin
        graph = create_viterbi_test_graph()
        config = Mycelia.ViterbiConfig()
        
        # Perfect sequence that matches the graph path
        observations = ["ATC", "TCG", "CGA", "GAT"]
        result = Mycelia.viterbi_decode_next(graph, observations, config)
        
        Test.@test result isa Mycelia.ViterbiPath
        Test.@test length(result.states) == length(observations)
        Test.@test result.log_probability > -Inf
        Test.@test !isempty(result.polished_sequence)
        
        # Check that states follow the observations
        for (i, state) in enumerate(result.states)
            Test.@test state.position == i
            Test.@test state isa Mycelia.ViterbiState
        end
        
        # Test with sequence containing errors
        noisy_observations = ["ATC", "TTG", "CGA", "GAT"]  # TTG instead of TCG
        noisy_result = Mycelia.viterbi_decode_next(graph, noisy_observations, config)
        
        Test.@test noisy_result isa Mycelia.ViterbiPath
        Test.@test length(noisy_result.corrections_made) > 0  # Should detect and correct error
        
        # Test empty input
        empty_result = Mycelia.viterbi_decode_next(graph, String[], config)
        Test.@test isempty(empty_result.states)
        Test.@test empty_result.log_probability == -Inf
    end
    
    Test.@testset "Sequence Polishing" begin
        graph = create_viterbi_test_graph()
        config = Mycelia.ViterbiConfig()
        
        # Perfect sequence
        perfect_seq = "ATCGAT"  # Should generate k-mers: ATC, TCG, CGA, GAT
        result = Mycelia.polish_sequence_next(graph, perfect_seq, config)
        
        Test.@test result isa Mycelia.ViterbiPath
        Test.@test !isempty(result.polished_sequence)
        Test.@test length(result.states) == 4  # Four 3-mers
        
        # Sequence with errors
        noisy_seq = "ATTGAT"  # Should generate: ATT, TTG, TGA, GAT
        noisy_result = Mycelia.polish_sequence_next(graph, noisy_seq, config)
        
        Test.@test noisy_result isa Mycelia.ViterbiPath
        # May have corrections depending on graph structure
        
        # Short sequence (shorter than k)
        short_seq = "AT"
        short_result = Mycelia.polish_sequence_next(graph, short_seq, config)
        Test.@test isempty(short_result.states)
        Test.@test short_result.polished_sequence == short_seq
    end
    
    Test.@testset "Error Correction" begin
        graph = create_viterbi_test_graph()
        config = Mycelia.ViterbiConfig()
        
        # Create test FASTA records
        perfect_record = FASTX.FASTA.Record("seq1", "ATCGAT")
        noisy_record = FASTX.FASTA.Record("seq2", "ATTGAT")
        sequences = [perfect_record, noisy_record]
        
        corrected = Mycelia.correct_errors_next(graph, sequences, config)
        
        Test.@test length(corrected) == length(sequences)
        Test.@test all(rec -> rec isa FASTX.FASTA.Record, corrected)
        
        # Check that IDs are modified
        for (i, rec) in enumerate(corrected)
            Test.@test occursin("polished", FASTX.identifier(rec))
        end
        
        # Test with empty input
        empty_corrected = Mycelia.correct_errors_next(graph, FASTX.FASTA.Record[], config)
        Test.@test isempty(empty_corrected)
    end
    
    Test.@testset "Batch Processing" begin
        graph = create_viterbi_test_graph()
        config = Mycelia.ViterbiConfig(batch_size=2)
        
        # Create multiple test sequences
        sequences = [
            FASTX.FASTA.Record("seq1", "ATCGAT"),
            FASTX.FASTA.Record("seq2", "ATTGAT"),
            FASTX.FASTA.Record("seq3", "ATCGAT"),
            FASTX.FASTA.Record("seq4", "ATCGTT")
        ]
        
        results = Mycelia.viterbi_batch_process(graph, sequences, config)
        
        Test.@test length(results) == length(sequences)
        Test.@test all(result -> result isa Mycelia.ViterbiPath, results)
        
        # Test with empty input
        empty_results = Mycelia.viterbi_batch_process(graph, FASTX.FASTA.Record[], config)
        Test.@test isempty(empty_results)
    end
    
    Test.@testset "Transition Probability Estimation" begin
        graph = create_viterbi_test_graph()
        
        # Create sequences that follow the graph structure
        sequences = [
            FASTX.FASTA.Record("seq1", "ATCGAT"),
            FASTX.FASTA.Record("seq2", "ATCGAT"),
            FASTX.FASTA.Record("seq3", "TCGAT")
        ]
        
        transition_probs = Mycelia.estimate_transition_probabilities(graph, sequences)
        
        # Should be a square matrix
        n_vertices = length(MetaGraphsNext.labels(graph))
        Test.@test size(transition_probs) == (n_vertices, n_vertices)
        
        # Rows should sum to 1 (or 0 for vertices with no outgoing transitions)
        for i in 1:size(transition_probs, 1)
            row_sum = sum(transition_probs[i, :])
            Test.@test row_sum ≈ 1.0 || row_sum ≈ 0.0
        end
        
        # Test with empty sequences
        empty_probs = Mycelia.estimate_transition_probabilities(graph, FASTX.FASTA.Record[])
        Test.@test size(empty_probs) == (n_vertices, n_vertices)
        Test.@test all(empty_probs .== 0.0)
    end
    
    Test.@testset "Reverse Complement Handling" begin
        # Test reverse complement function
        Test.@test Mycelia.reverse_complement("ATC") == "GAT"
        Test.@test Mycelia.reverse_complement("ATCG") == "CGAT"
        Test.@test Mycelia.reverse_complement("AAAA") == "TTTT"
        Test.@test Mycelia.reverse_complement("GCGC") == "GCGC"
        
        # Test with lowercase
        Test.@test Mycelia.reverse_complement("atc") == "gat"
        
        # Test strand-aware state
        reverse_state = Mycelia.ViterbiState("ATC", Mycelia.Rhizomorph.Reverse, 1.0, 1)
        config = Mycelia.ViterbiConfig()
        
        # Should match reverse complement
        prob = Mycelia.calculate_emission_probability(reverse_state, "GAT", config)
        Test.@test prob == config.match_prob
    end
    
    Test.@testset "Memory and Performance Configuration" begin
        config = Mycelia.ViterbiConfig(
            memory_limit=512,  # 512 MB
            use_log_space=false,
            batch_size=100
        )
        
        Test.@test config.memory_limit == 512
        Test.@test config.use_log_space == false
        Test.@test config.batch_size == 100
        
        graph = create_viterbi_test_graph()
        observations = ["ATC", "TCG", "CGA"]
        
        # Should work with log_space disabled
        result = Mycelia.viterbi_decode_next(graph, observations, config)
        Test.@test result isa Mycelia.ViterbiPath
        Test.@test !isempty(result.states)
    end
    
    Test.@testset "Integration with Real K-mer Graph" begin
        # Create a real k-mer graph and test Viterbi on it
        seq1 = FASTX.FASTA.Record("test1", "ATCGATCG")
        seq2 = FASTX.FASTA.Record("test2", "TCGATCGA")
        observations = [seq1, seq2]
        
        graph = Mycelia.Rhizomorph.build_kmer_graph(observations, 3; dataset_id="test", mode=:singlestrand)
        
        if !isempty(MetaGraphsNext.labels(graph))
            config = Mycelia.ViterbiConfig()
            test_sequence = "ATCGATC"
            
            result = Mycelia.polish_sequence_next(graph, test_sequence, config)
            Test.@test result isa Mycelia.ViterbiPath
            Test.@test !isempty(result.polished_sequence)
            
            # Test error correction on the original sequences
            corrected = Mycelia.correct_errors_next(graph, observations, config)
            Test.@test length(corrected) == length(observations)
        end
    end
end
