import Eisenia
import Test
import Documenter
import BioSequences

Test.@testset "documentation" begin
    # https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Doctesting-as-Part-of-Testing
    # doctest(Eisenia; manual = false)
    Documenter.doctest(Eisenia)
end


Test.@testset "Don't support even kmers" begin
    k = 2
    seqlen = 2
    msg = "Even kmers are not supported"
    Test.@test_throws ErrorException(msg) Eisenia.KmerGraph(BioSequences.DNAMer{k}, [BioSequences.randdnaseq(seqlen)])
end

Test.@testset "viterbi, error-free" begin
    Test.@testset "k=1, error_rate = 0.0, # observations = 1" begin
        k = 1
        for seqlen in k:k+2
            error_rate = 0.0
            for nucleotides in Iterators.product([Eisenia.DNA_ALPHABET for i in 1:seqlen]...)
                # convert tuple into actual longdnaseq
                observation = BioSequences.LongDNASeq([nucleotides...])
                graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
            end
        end
    end

    Test.@testset "k=3, error_rate = 0.0, # observations = 1" begin
        k = 3
        for seqlen in k:k+2
            error_rate = 0.0
            for nucleotides in Iterators.product([Eisenia.DNA_ALPHABET for i in 1:seqlen]...)
                # convert tuple into actual longdnaseq
                observation = BioSequences.LongDNASeq([nucleotides...])
                graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
            end
        end
    end
end

Test.@testset "viterbi, error-free, both orientations" begin
    Test.@testset "k=1, error_rate = 0.0, # observations = 1" begin
        k = 1
        for seqlen in k:k+2
            error_rate = 0.0
            for nucleotides in Iterators.product([Eisenia.DNA_ALPHABET for i in 1:seqlen]...)
                # convert tuple into actual longdnaseq
                observation = BioSequences.LongDNASeq([nucleotides...])
                rc_observation = BioSequences.reverse_complement(observation)
                
                graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
                
                graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation, rc_observation])
                optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
            end
        end
    end

    Test.@testset "k=3, error_rate = 0.0, # observations = 1" begin
        k = 3
        for seqlen in k:k+2
            error_rate = 0.0
            for nucleotides in Iterators.product([Eisenia.DNA_ALPHABET for i in 1:seqlen]...)
                # convert tuple into actual longdnaseq
                observation = BioSequences.LongDNASeq([nucleotides...])
                rc_observation = BioSequences.reverse_complement(observation)
                
                graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
                
                graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation, rc_observation])
                optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
            end
        end
    end
end

Test.@testset "viterbi, variable-error-rates" begin
    Test.@testset "k=1" begin
        k = 1
        for seqlen in k:k+2
            for error_rate in [0.01, 0.05, 0.10, 0.15]
                for nucleotides in Iterators.product([Eisenia.DNA_ALPHABET for i in 1:seqlen]...)
                    # convert tuple into actual longdnaseq
                    observation = BioSequences.LongDNASeq([nucleotides...])
                    graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                    optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                    Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                end
            end
        end
    end

    Test.@testset "k=3" begin
        k = 3
        for seqlen in k:k+2
            for error_rate in [0.01, 0.05, 0.10, 0.15]
                for nucleotides in Iterators.product([Eisenia.DNA_ALPHABET for i in 1:seqlen]...)
                    # convert tuple into actual longdnaseq
                    observation = BioSequences.LongDNASeq([nucleotides...])
                    graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                    optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                    Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
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
                for nucleotides in Iterators.product([Eisenia.DNA_ALPHABET for i in 1:seqlen]...)
                    # convert tuple into actual longdnaseq
                    observation = BioSequences.LongDNASeq([nucleotides...])
                    rc_observation = BioSequences.reverse_complement(observation)

                    graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                    optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                    Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation

                    graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation, rc_observation])
                    optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                    Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                    optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                    Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
                end
            end
        end
    end

    Test.@testset "k=3" begin
        k = 3
        for seqlen in k:k+2
            for error_rate in [0.01, 0.05, 0.10, 0.15]
                for nucleotides in Iterators.product([Eisenia.DNA_ALPHABET for i in 1:seqlen]...)
                    # convert tuple into actual longdnaseq
                    observation = BioSequences.LongDNASeq([nucleotides...])
                    rc_observation = BioSequences.reverse_complement(observation)

                    graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation])
                    optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                    Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation

                    graph = Eisenia.KmerGraph(BioSequences.DNAMer{k}, [observation, rc_observation])
                    optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, observation, error_rate)
                    Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == observation
                    optimal_path = Eisenia.viterbi_maximum_likelihood_path(graph, rc_observation, error_rate)
                    Test.@test Eisenia.oriented_path_to_sequence(optimal_path[1], graph.kmers) == rc_observation
                end
            end
        end
    end
end

# # Benchmarking
# using PkgBenchmark
# # compares this current state against master branch
# results = judge(Eisenia, "master")
# export_markdown(stdout, results, export_invariants=true)
