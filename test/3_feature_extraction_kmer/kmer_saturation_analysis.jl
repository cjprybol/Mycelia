import Test
import Mycelia
import BioSequences
import Kmers

Test.@testset "K-mer saturation analysis" begin
    sequences = [
        BioSequences.LongDNA{4}("ATGCA"),
        BioSequences.LongDNA{4}("ATGGA"),
        BioSequences.LongDNA{4}("ATGTA"),
        BioSequences.LongDNA{4}("ATGAA"),
    ]

    result = Mycelia.analyze_kmer_saturation(
        sequences,
        [3, 4];
        max_sampling_points=3,
        min_sequences=2,
    )

    Test.@test result.optimal_k in [3, 4]
    Test.@test length(result.saturation_levels) == 2
    Test.@test length(result.sampling_points) == length(result.unique_kmer_counts[3])
    Test.@test length(result.sampling_points) == length(result.unique_kmer_counts[4])
    Test.@test all(k -> haskey(result.unique_kmer_counts, k), [3, 4])

    Test.@test_throws ArgumentError Mycelia.analyze_kmer_saturation(sequences, Int[]; min_sequences=2)
    Test.@test_throws ArgumentError Mycelia.analyze_kmer_saturation(sequences, [3]; min_sequences=10)

    kmer_counts = Dict(
        Kmers.DNAKmer{3}(BioSequences.LongDNA{4}("ATG")) => 5,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{4}("TGC")) => 4,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{4}("GCA")) => 1,
    )

    threshold_result = Mycelia.find_optimal_assembly_threshold(
        kmer_counts,
        3;
        min_threshold=1,
        max_threshold=5,
    )

    Test.@test threshold_result.optimal_threshold in threshold_result.thresholds_tested
    Test.@test length(threshold_result.thresholds_tested) == length(threshold_result.connectivity_scores)
    Test.@test length(threshold_result.thresholds_tested) == length(threshold_result.component_counts)
    Test.@test length(threshold_result.thresholds_tested) == length(threshold_result.average_connectivity)
end
