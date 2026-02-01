import Test
import Mycelia

Test.@testset "Coverage Clustering" begin
    kmer_counts = Dict(
        "AAA" => 1,
        "AAT" => 2,
        "ATG" => 15,
        "TGC" => 18
    )

    result = Mycelia.k_medoids_coverage_clustering(kmer_counts, 2; max_iterations = 20)
    Test.@test length(result.assignments) == length(kmer_counts)
    Test.@test length(result.coverage_thresholds) == 2
    Test.@test result.cost >= 0.0
    Test.@test length(result.silhouettes) == length(kmer_counts)

    stats = Mycelia.analyze_coverage_clustering(result, kmer_counts)
    Test.@test stats.error_cluster_size + stats.signal_cluster_size == length(kmer_counts)
    Test.@test 0.0 <= stats.separation_quality <= 1.0
    Test.@test stats.optimal_threshold >= 0.0

    filtered, filter_stats,
    removed = Mycelia.automatic_error_filtering(
        kmer_counts;
        separation_threshold = 0.0,
        max_error_rate = 1.0
    )
    Test.@test length(filtered) + length(removed) == length(kmer_counts)
    Test.@test filter_stats.error_cluster_size + filter_stats.signal_cluster_size ==
               length(kmer_counts)

    Test.@test_throws ArgumentError Mycelia.k_medoids_coverage_clustering(kmer_counts, 10)
    Test.@test_throws ArgumentError Mycelia._calculate_distance([1.0], [1.0], :invalid)

    Test.@test Mycelia._calculate_distance([1.0, 2.0], [1.0, 2.0], :euclidean) == 0.0
    Test.@test Mycelia._calculate_distance([1.0, 2.0], [1.0, 2.0], :manhattan) == 0.0
    Test.@test Mycelia._calculate_distance([1.0, 0.0], [1.0, 0.0], :cosine) == 0.0
end
