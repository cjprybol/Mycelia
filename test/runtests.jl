import Mycelia
import Test
import Documenter
import BioSequences
import Random
import Primes
import Graphs

Test.@testset "Don't support even kmers" begin
    k = 2
    seqlen = 2
    msg = "Even kmers are not supported"
    Test.@test_throws ErrorException(msg) Mycelia.KmerGraph(BioSequences.DNAMer{k}, [BioSequences.randdnaseq(seqlen)])
end

Test.@testset "viterbi, error-free" begin
    Test.@testset "k=1, error_rate = 0.0, # observations = 1" begin
        k = 1
        for seqlen in k:k+2
            error_rate = 0.0
            for nucleotides in Iterators.product([Mycelia.DNA_ALPHABET for i in 1:seqlen]...)
                # convert tuple into actual longdnaseq
                observation = BioSequences.LongDNASeq([nucleotides...])
                graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
            end
        end
    end

    Test.@testset "k=3, error_rate = 0.0, # observations = 1" begin
        k = 3
        for seqlen in k:k+2
            error_rate = 0.0
            for nucleotides in Iterators.product([Mycelia.DNA_ALPHABET for i in 1:seqlen]...)
                # convert tuple into actual longdnaseq
                observation = BioSequences.LongDNASeq([nucleotides...])
                graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
            end
        end
    end
end

Test.@testset "viterbi, error-free, both orientations" begin
    Test.@testset "k=1, error_rate = 0.0, # observations = 1" begin
        k = 1
        for seqlen in k:k+2
            error_rate = 0.0
            for nucleotides in Iterators.product([Mycelia.DNA_ALPHABET for i in 1:seqlen]...)
                # convert tuple into actual longdnaseq
                observation = BioSequences.LongDNASeq([nucleotides...])
                rc_observation = BioSequences.reverse_complement(observation)
                
                graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
                
                graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation, rc_observation])
                optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
            end
        end
    end

    Test.@testset "k=3, error_rate = 0.0, # observations = 1" begin
        k = 3
        for seqlen in k:k+2
            error_rate = 0.0
            for nucleotides in Iterators.product([Mycelia.DNA_ALPHABET for i in 1:seqlen]...)
                # convert tuple into actual longdnaseq
                observation = BioSequences.LongDNASeq([nucleotides...])
                rc_observation = BioSequences.reverse_complement(observation)
                
                graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
                
                graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation, rc_observation])
                optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
            end
        end
    end
end

Test.@testset "viterbi, variable-error-rates" begin
    Test.@testset "k=1" begin
        k = 1
        for seqlen in k:k+2
            for error_rate in [0.01, 0.05, 0.10, 0.15]
                for nucleotides in Iterators.product([Mycelia.DNA_ALPHABET for i in 1:seqlen]...)
                    # convert tuple into actual longdnaseq
                    observation = BioSequences.LongDNASeq([nucleotides...])
                    graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                    optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                    Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                end
            end
        end
    end

    Test.@testset "k=3" begin
        k = 3
        for seqlen in k:k+2
            for error_rate in [0.01, 0.05, 0.10, 0.15]
                for nucleotides in Iterators.product([Mycelia.DNA_ALPHABET for i in 1:seqlen]...)
                    # convert tuple into actual longdnaseq
                    observation = BioSequences.LongDNASeq([nucleotides...])
                    graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                    optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                    Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                end
            end
        end
    end
end

Test.@testset "viterbi, variable-error-rates, both orientations" begin
    Test.@testset "k=1" begin
        k = 1
        for seqlen in k:k+2
            for error_rate in [0.01, 0.05, 0.10, 0.15]
                for nucleotides in Iterators.product([Mycelia.DNA_ALPHABET for i in 1:seqlen]...)
                    # convert tuple into actual longdnaseq
                    observation = BioSequences.LongDNASeq([nucleotides...])
                    rc_observation = BioSequences.reverse_complement(observation)

                    graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                    optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                    Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation

                    graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation, rc_observation])
                    optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                    Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                    optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                    Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
                end
            end
        end
    end

    Test.@testset "k=3" begin
        k = 3
        for seqlen in k:k+2
            for error_rate in [0.01, 0.05, 0.10, 0.15]
                for nucleotides in Iterators.product([Mycelia.DNA_ALPHABET for i in 1:seqlen]...)
                    # convert tuple into actual longdnaseq
                    observation = BioSequences.LongDNASeq([nucleotides...])
                    rc_observation = BioSequences.reverse_complement(observation)

                    graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                    optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                    Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation

                    graph = Mycelia.KmerGraph(BioSequences.DNAMer{k}, [observation, rc_observation])
                    optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                    Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                    optimal_path = Mycelia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                    Test.@test Mycelia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
                end
            end
        end
    end
end

Test.@testset "graph reconstruction with tip clipping" begin
    n_sequences = 1
    seqlen = 10
    n_observations = 100
    error_rate = 0.05

    sequences = [BioSequences.randdnaseq(Random.seed!(i), seqlen) for i in 1:n_sequences]    
    Random.seed!(1)
    observations = [
        Mycelia.observe(rand(sequences), error_rate = error_rate) 
            for i in 1:n_observations
    ]
    
    graph, corrected_observations = Mycelia.iterate_until_convergence(Primes.primes(3, 7), observations, error_rate);
    pruned_graph = Mycelia.clip_low_coverage_tips(graph, corrected_observations)
    
    for connected_component in Graphs.connected_components(pruned_graph.graph)
        primary_path = Mycelia.maximum_likelihood_walk(pruned_graph, connected_component)
        reconstructed_sequence = Mycelia.path_to_sequence(pruned_graph, primary_path)
        Test.@test Mycelia.is_equivalent(reconstructed_sequence, first(sequences))
    end
end

Test.@testset "documentation" begin
    # https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Doctesting-as-Part-of-Testing
    # doctest(Mycelia; manual = false)
    Documenter.doctest(Mycelia)
end

# # Benchmarking
# using PkgBenchmark
# # compares this current state against master branch
# results = judge(Mycelia, "master")
# export_markdown(stdout, results, export_invariants=true)
