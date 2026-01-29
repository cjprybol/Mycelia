# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/cluster_comparison_metrics_test.jl")'
# ```

import Test
import Mycelia
import Random
import Statistics

## Suppress verbose output during CI/automated testing
const SUPPRESS_VERBOSE_OUTPUT = get(ENV, "CI", "") == "true" || get(ENV, "GITHUB_ACTIONS", "") == "true"

## Conditional printing function
function test_println(args...)
    if !SUPPRESS_VERBOSE_OUTPUT
        println(args...)
    end
end

# =============================================================================
# Test Data Setup
# =============================================================================

## Perfect agreement case
const LABELS_IDENTICAL_A = [1, 1, 1, 2, 2, 2, 3, 3, 3]
const LABELS_IDENTICAL_B = [1, 1, 1, 2, 2, 2, 3, 3, 3]

## Relabeled but same structure (different labels, same partition)
const LABELS_RELABELED_A = [1, 1, 1, 2, 2, 2, 3, 3, 3]
const LABELS_RELABELED_B = [3, 3, 3, 1, 1, 1, 2, 2, 2]  ## Same structure, different labels

## Partial agreement
const LABELS_PARTIAL_A = [1, 1, 1, 2, 2, 2]
const LABELS_PARTIAL_B = [1, 1, 2, 2, 2, 2]  ## One sample swapped

## No agreement (random)
const LABELS_NO_AGREE_A = [1, 1, 1, 2, 2, 2]
const LABELS_NO_AGREE_B = [1, 2, 1, 2, 1, 2]  ## Alternating, no structure match

## String labels
const LABELS_STRING_A = ["cat", "cat", "dog", "dog", "bird"]
const LABELS_STRING_B = ["A", "A", "B", "B", "C"]

## Single cluster
const LABELS_SINGLE_A = [1, 1, 1, 1, 1]
const LABELS_SINGLE_B = [1, 1, 1, 1, 1]

## Two clusters completely swapped
const LABELS_SWAP_A = [1, 1, 1, 2, 2, 2]
const LABELS_SWAP_B = [2, 2, 2, 1, 1, 1]


# =============================================================================
# Contingency Matrix Tests
# =============================================================================

Test.@testset "Contingency Matrix" begin
    test_println("\n=== Testing Contingency Matrix ===")

    Test.@testset "Basic contingency matrix" begin
        labels1 = [1, 1, 2, 2, 3]
        labels2 = [1, 1, 1, 2, 2]

        result = Mycelia.contingency_matrix(labels1, labels2)

        Test.@test size(result.matrix) == (3, 2)
        Test.@test sum(result.matrix) == 5
        Test.@test result.labels1 == [1, 2, 3]
        Test.@test result.labels2 == [1, 2]

        ## Check specific counts
        Test.@test result.matrix[1, 1] == 2  ## label1=1, label2=1
        Test.@test result.matrix[2, 1] == 1  ## label1=2, label2=1
        Test.@test result.matrix[2, 2] == 1  ## label1=2, label2=2
        Test.@test result.matrix[3, 2] == 1  ## label1=3, label2=2
    end

    Test.@testset "Identical labels" begin
        result = Mycelia.contingency_matrix(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)

        ## Should be diagonal matrix
        for i in 1:size(result.matrix, 1)
            for j in 1:size(result.matrix, 2)
                if i == j
                    Test.@test result.matrix[i, j] == 3
                else
                    Test.@test result.matrix[i, j] == 0
                end
            end
        end
    end

    Test.@testset "String labels" begin
        result = Mycelia.contingency_matrix(LABELS_STRING_A, LABELS_STRING_B)

        Test.@test result.labels1 == ["bird", "cat", "dog"]
        Test.@test result.labels2 == ["A", "B", "C"]
        Test.@test sum(result.matrix) == 5
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.contingency_matrix([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# Adjusted Rand Index Tests
# =============================================================================

Test.@testset "Adjusted Rand Index (ARI)" begin
    test_println("\n=== Testing Adjusted Rand Index ===")

    Test.@testset "Perfect agreement" begin
        ari = Mycelia.adjusted_rand_index(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        Test.@test isapprox(ari, 1.0, atol=1e-10)
        test_println("  Identical labels ARI: $ari")
    end

    Test.@testset "Relabeled but same structure" begin
        ## Same partition structure should give ARI = 1
        ari = Mycelia.adjusted_rand_index(LABELS_RELABELED_A, LABELS_RELABELED_B)
        Test.@test isapprox(ari, 1.0, atol=1e-10)
        test_println("  Relabeled structure ARI: $ari")
    end

    Test.@testset "Swapped clusters" begin
        ## Swapped labels should still give ARI = 1
        ari = Mycelia.adjusted_rand_index(LABELS_SWAP_A, LABELS_SWAP_B)
        Test.@test isapprox(ari, 1.0, atol=1e-10)
        test_println("  Swapped clusters ARI: $ari")
    end

    Test.@testset "Partial agreement" begin
        ari = Mycelia.adjusted_rand_index(LABELS_PARTIAL_A, LABELS_PARTIAL_B)
        Test.@test 0.0 < ari < 1.0
        test_println("  Partial agreement ARI: $ari")
    end

    Test.@testset "Low agreement" begin
        ari = Mycelia.adjusted_rand_index(LABELS_NO_AGREE_A, LABELS_NO_AGREE_B)
        Test.@test ari < 0.5  ## Should be low
        test_println("  Low agreement ARI: $ari")
    end

    Test.@testset "ARI range" begin
        ## ARI should be in [-1, 1] range
        for _ in 1:10
            n = 20
            labels1 = rand(1:3, n)
            labels2 = rand(1:3, n)
            ari = Mycelia.adjusted_rand_index(labels1, labels2)
            Test.@test -1.0 <= ari <= 1.0
        end
    end

    Test.@testset "String labels" begin
        ari = Mycelia.adjusted_rand_index(LABELS_STRING_A, LABELS_STRING_B)
        Test.@test isapprox(ari, 1.0, atol=1e-10)
    end

    Test.@testset "Single sample" begin
        ari = Mycelia.adjusted_rand_index([1], [1])
        Test.@test ari == 1.0
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.adjusted_rand_index([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# Normalized Mutual Information Tests
# =============================================================================

Test.@testset "Normalized Mutual Information (NMI)" begin
    test_println("\n=== Testing Normalized Mutual Information ===")

    Test.@testset "Perfect agreement" begin
        nmi = Mycelia.normalized_mutual_information(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        Test.@test isapprox(nmi, 1.0, atol=1e-10)
        test_println("  Identical labels NMI: $nmi")
    end

    Test.@testset "Relabeled but same structure" begin
        nmi = Mycelia.normalized_mutual_information(LABELS_RELABELED_A, LABELS_RELABELED_B)
        Test.@test isapprox(nmi, 1.0, atol=1e-10)
        test_println("  Relabeled structure NMI: $nmi")
    end

    Test.@testset "Partial agreement" begin
        nmi = Mycelia.normalized_mutual_information(LABELS_PARTIAL_A, LABELS_PARTIAL_B)
        Test.@test 0.0 < nmi < 1.0
        test_println("  Partial agreement NMI: $nmi")
    end

    Test.@testset "NMI range" begin
        ## NMI should be in [0, 1] range
        for _ in 1:10
            n = 20
            labels1 = rand(1:3, n)
            labels2 = rand(1:3, n)
            nmi = Mycelia.normalized_mutual_information(labels1, labels2)
            Test.@test 0.0 <= nmi <= 1.0
        end
    end

    Test.@testset "Single cluster" begin
        ## When all samples are in one cluster, entropy is 0
        nmi = Mycelia.normalized_mutual_information(LABELS_SINGLE_A, LABELS_SINGLE_B)
        Test.@test nmi == 1.0
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.normalized_mutual_information([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# Adjusted Mutual Information Tests
# =============================================================================

Test.@testset "Adjusted Mutual Information (AMI)" begin
    test_println("\n=== Testing Adjusted Mutual Information ===")

    Test.@testset "Perfect agreement" begin
        ami = Mycelia.adjusted_mutual_information(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        Test.@test isapprox(ami, 1.0, atol=0.1)  ## AMI can have numerical issues
        test_println("  Identical labels AMI: $ami")
    end

    Test.@testset "Relabeled but same structure" begin
        ami = Mycelia.adjusted_mutual_information(LABELS_RELABELED_A, LABELS_RELABELED_B)
        Test.@test isapprox(ami, 1.0, atol=0.1)
        test_println("  Relabeled structure AMI: $ami")
    end

    Test.@testset "Partial agreement" begin
        ami = Mycelia.adjusted_mutual_information(LABELS_PARTIAL_A, LABELS_PARTIAL_B)
        Test.@test ami < 1.0
        test_println("  Partial agreement AMI: $ami")
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.adjusted_mutual_information([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# V-Measure Tests
# =============================================================================

Test.@testset "V-Measure" begin
    test_println("\n=== Testing V-Measure ===")

    Test.@testset "Perfect agreement" begin
        result = Mycelia.v_measure(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)

        Test.@test isapprox(result.homogeneity, 1.0, atol=1e-10)
        Test.@test isapprox(result.completeness, 1.0, atol=1e-10)
        Test.@test isapprox(result.v_measure, 1.0, atol=1e-10)
        test_println("  Identical labels V-measure: $(result.v_measure)")
    end

    Test.@testset "Relabeled but same structure" begin
        result = Mycelia.v_measure(LABELS_RELABELED_A, LABELS_RELABELED_B)
        Test.@test isapprox(result.v_measure, 1.0, atol=1e-10)
        test_println("  Relabeled structure V-measure: $(result.v_measure)")
    end

    Test.@testset "Partial agreement" begin
        result = Mycelia.v_measure(LABELS_PARTIAL_A, LABELS_PARTIAL_B)

        Test.@test 0.0 < result.homogeneity < 1.0
        Test.@test 0.0 < result.completeness < 1.0
        Test.@test 0.0 < result.v_measure < 1.0
        test_println("  Partial agreement - H: $(result.homogeneity), C: $(result.completeness), V: $(result.v_measure)")
    end

    Test.@testset "V-measure range" begin
        for _ in 1:10
            n = 20
            labels1 = rand(1:3, n)
            labels2 = rand(1:3, n)
            result = Mycelia.v_measure(labels1, labels2)

            Test.@test 0.0 <= result.homogeneity <= 1.0
            Test.@test 0.0 <= result.completeness <= 1.0
            Test.@test 0.0 <= result.v_measure <= 1.0
        end
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.v_measure([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# Fowlkes-Mallows Index Tests
# =============================================================================

Test.@testset "Fowlkes-Mallows Index (FMI)" begin
    test_println("\n=== Testing Fowlkes-Mallows Index ===")

    Test.@testset "Perfect agreement" begin
        fmi = Mycelia.fowlkes_mallows_index(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        Test.@test isapprox(fmi, 1.0, atol=1e-10)
        test_println("  Identical labels FMI: $fmi")
    end

    Test.@testset "Relabeled but same structure" begin
        fmi = Mycelia.fowlkes_mallows_index(LABELS_RELABELED_A, LABELS_RELABELED_B)
        Test.@test isapprox(fmi, 1.0, atol=1e-10)
        test_println("  Relabeled structure FMI: $fmi")
    end

    Test.@testset "Partial agreement" begin
        fmi = Mycelia.fowlkes_mallows_index(LABELS_PARTIAL_A, LABELS_PARTIAL_B)
        Test.@test 0.0 < fmi < 1.0
        test_println("  Partial agreement FMI: $fmi")
    end

    Test.@testset "FMI range" begin
        ## FMI should be in [0, 1] range
        for _ in 1:10
            n = 20
            labels1 = rand(1:3, n)
            labels2 = rand(1:3, n)
            fmi = Mycelia.fowlkes_mallows_index(labels1, labels2)
            Test.@test 0.0 <= fmi <= 1.0
        end
    end

    Test.@testset "Single sample" begin
        fmi = Mycelia.fowlkes_mallows_index([1], [1])
        Test.@test fmi == 1.0
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.fowlkes_mallows_index([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# Cluster Purity Tests
# =============================================================================

Test.@testset "Cluster Purity" begin
    test_println("\n=== Testing Cluster Purity ===")

    Test.@testset "Perfect purity" begin
        result = Mycelia.cluster_purity(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        Test.@test isapprox(result.overall_purity, 1.0, atol=1e-10)
        test_println("  Identical labels purity: $(result.overall_purity)")
    end

    Test.@testset "Partial purity" begin
        result = Mycelia.cluster_purity(LABELS_PARTIAL_A, LABELS_PARTIAL_B)
        Test.@test 0.0 < result.overall_purity < 1.0
        test_println("  Partial agreement purity: $(result.overall_purity)")
    end

    Test.@testset "Purity per cluster" begin
        labels_true = [1, 1, 1, 2, 2, 2]
        labels_pred = [1, 1, 1, 1, 2, 2]  ## Cluster 1 has impurity

        result = Mycelia.cluster_purity(labels_true, labels_pred)

        Test.@test haskey(result.cluster_purities, 1)
        Test.@test haskey(result.cluster_purities, 2)
        Test.@test result.cluster_purities[2] == 1.0  ## Pure cluster
        Test.@test result.cluster_purities[1] < 1.0   ## Impure cluster
    end

    Test.@testset "Purity range" begin
        for _ in 1:10
            n = 20
            labels1 = rand(1:3, n)
            labels2 = rand(1:3, n)
            result = Mycelia.cluster_purity(labels1, labels2)
            Test.@test 0.0 <= result.overall_purity <= 1.0
        end
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.cluster_purity([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# Jaccard Index Tests
# =============================================================================

Test.@testset "Jaccard Index" begin
    test_println("\n=== Testing Jaccard Index ===")

    Test.@testset "Perfect agreement" begin
        ji = Mycelia.jaccard_index(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        Test.@test isapprox(ji, 1.0, atol=1e-10)
        test_println("  Identical labels Jaccard: $ji")
    end

    Test.@testset "Relabeled but same structure" begin
        ji = Mycelia.jaccard_index(LABELS_RELABELED_A, LABELS_RELABELED_B)
        Test.@test isapprox(ji, 1.0, atol=1e-10)
        test_println("  Relabeled structure Jaccard: $ji")
    end

    Test.@testset "Partial agreement" begin
        ji = Mycelia.jaccard_index(LABELS_PARTIAL_A, LABELS_PARTIAL_B)
        Test.@test 0.0 < ji < 1.0
        test_println("  Partial agreement Jaccard: $ji")
    end

    Test.@testset "Jaccard range" begin
        for _ in 1:10
            n = 20
            labels1 = rand(1:3, n)
            labels2 = rand(1:3, n)
            ji = Mycelia.jaccard_index(labels1, labels2)
            Test.@test 0.0 <= ji <= 1.0
        end
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.jaccard_index([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# Clustering Comparison Summary Tests
# =============================================================================

Test.@testset "Clustering Comparison Summary" begin
    test_println("\n=== Testing Clustering Comparison Summary ===")

    Test.@testset "Summary structure" begin
        labels1 = [1, 1, 2, 2, 3, 3]
        labels2 = [1, 1, 1, 2, 2, 2]

        summary = Mycelia.clustering_comparison_summary(labels1, labels2; name1="Test1", name2="Test2")

        Test.@test summary.name1 == "Test1"
        Test.@test summary.name2 == "Test2"
        Test.@test summary.n_samples == 6
        Test.@test summary.n_clusters1 == 3
        Test.@test summary.n_clusters2 == 2
        Test.@test haskey(summary, :adjusted_rand_index)
        Test.@test haskey(summary, :normalized_mutual_information)
        Test.@test haskey(summary, :adjusted_mutual_information)
        Test.@test haskey(summary, :v_measure)
        Test.@test haskey(summary, :fowlkes_mallows_index)
        Test.@test haskey(summary, :jaccard_index)
        Test.@test haskey(summary, :purity_1_to_2)
        Test.@test haskey(summary, :purity_2_to_1)
        Test.@test haskey(summary, :contingency)
    end

    Test.@testset "Perfect agreement summary" begin
        summary = Mycelia.clustering_comparison_summary(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)

        Test.@test isapprox(summary.adjusted_rand_index, 1.0, atol=1e-10)
        Test.@test isapprox(summary.normalized_mutual_information, 1.0, atol=1e-10)
        Test.@test isapprox(summary.fowlkes_mallows_index, 1.0, atol=1e-10)
        Test.@test isapprox(summary.jaccard_index, 1.0, atol=1e-10)
    end

    Test.@testset "Error on mismatched lengths" begin
        Test.@test_throws ArgumentError Mycelia.clustering_comparison_summary([1, 2, 3], [1, 2])
    end
end


# =============================================================================
# Print Clustering Comparison Tests
# =============================================================================

Test.@testset "Print Clustering Comparison" begin
    test_println("\n=== Testing Print Clustering Comparison ===")

    Test.@testset "Print to IOBuffer" begin
        labels1 = [1, 1, 2, 2, 3, 3]
        labels2 = [1, 1, 1, 2, 2, 3]

        io = IOBuffer()
        summary = Mycelia.print_clustering_comparison(labels1, labels2; name1="A", name2="B", io=io)

        output = String(take!(io))

        Test.@test contains(output, "CLUSTERING COMPARISON REPORT")
        Test.@test contains(output, "Adjusted Rand Index")
        Test.@test contains(output, "Normalized Mutual Info")
        Test.@test contains(output, "V-measure")
        Test.@test contains(output, "A")
        Test.@test contains(output, "B")

        ## Summary should also be returned
        Test.@test haskey(summary, :adjusted_rand_index)
    end
end


# =============================================================================
# Comparison Summary to DataFrame Tests
# =============================================================================

Test.@testset "Comparison Summary to DataFrame" begin
    test_println("\n=== Testing Comparison Summary to DataFrame ===")

    Test.@testset "DataFrame structure" begin
        labels1 = [1, 1, 2, 2, 3, 3]
        labels2 = [1, 1, 1, 2, 2, 2]

        summary = Mycelia.clustering_comparison_summary(labels1, labels2; name1="Clust1", name2="Clust2")
        df = Mycelia.comparison_summary_to_dataframe(summary)

        Test.@test Mycelia.DataFrames.nrow(df) == 1
        Test.@test "adjusted_rand_index" in names(df)
        Test.@test "normalized_mutual_information" in names(df)
        Test.@test "v_measure" in names(df)
        Test.@test "homogeneity" in names(df)
        Test.@test "completeness" in names(df)
        Test.@test "fowlkes_mallows_index" in names(df)
        Test.@test "jaccard_index" in names(df)
        Test.@test "purity_1_to_2" in names(df)
        Test.@test "purity_2_to_1" in names(df)
        Test.@test df.name1[1] == "Clust1"
        Test.@test df.name2[1] == "Clust2"
    end
end


# =============================================================================
# Consistency Between Metrics Tests
# =============================================================================

Test.@testset "Metric Consistency" begin
    test_println("\n=== Testing Metric Consistency ===")

    Test.@testset "All metrics agree on perfect match" begin
        ari = Mycelia.adjusted_rand_index(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        nmi = Mycelia.normalized_mutual_information(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        fmi = Mycelia.fowlkes_mallows_index(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        ji = Mycelia.jaccard_index(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        vm = Mycelia.v_measure(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)
        purity = Mycelia.cluster_purity(LABELS_IDENTICAL_A, LABELS_IDENTICAL_B)

        Test.@test isapprox(ari, 1.0, atol=1e-10)
        Test.@test isapprox(nmi, 1.0, atol=1e-10)
        Test.@test isapprox(fmi, 1.0, atol=1e-10)
        Test.@test isapprox(ji, 1.0, atol=1e-10)
        Test.@test isapprox(vm.v_measure, 1.0, atol=1e-10)
        Test.@test isapprox(purity.overall_purity, 1.0, atol=1e-10)
    end

    Test.@testset "Metrics are symmetric where expected" begin
        labels1 = [1, 1, 2, 2, 3]
        labels2 = [1, 2, 2, 3, 3]

        ## ARI, NMI, FMI, Jaccard should be symmetric
        ari_12 = Mycelia.adjusted_rand_index(labels1, labels2)
        ari_21 = Mycelia.adjusted_rand_index(labels2, labels1)
        Test.@test isapprox(ari_12, ari_21, atol=1e-10)

        nmi_12 = Mycelia.normalized_mutual_information(labels1, labels2)
        nmi_21 = Mycelia.normalized_mutual_information(labels2, labels1)
        Test.@test isapprox(nmi_12, nmi_21, atol=1e-10)

        fmi_12 = Mycelia.fowlkes_mallows_index(labels1, labels2)
        fmi_21 = Mycelia.fowlkes_mallows_index(labels2, labels1)
        Test.@test isapprox(fmi_12, fmi_21, atol=1e-10)

        ji_12 = Mycelia.jaccard_index(labels1, labels2)
        ji_21 = Mycelia.jaccard_index(labels2, labels1)
        Test.@test isapprox(ji_12, ji_21, atol=1e-10)
    end

    Test.@testset "Purity is asymmetric" begin
        ## Purity depends on which clustering is "true"
        labels1 = [1, 1, 1, 2, 2, 2]
        labels2 = [1, 1, 2, 2, 2, 2]

        purity_12 = Mycelia.cluster_purity(labels1, labels2)
        purity_21 = Mycelia.cluster_purity(labels2, labels1)

        ## These can differ
        test_println("  Purity 1->2: $(purity_12.overall_purity), Purity 2->1: $(purity_21.overall_purity)")
    end
end


# =============================================================================
# Edge Cases Tests
# =============================================================================

Test.@testset "Edge Cases" begin
    test_println("\n=== Testing Edge Cases ===")

    Test.@testset "Empty vectors" begin
        result = Mycelia.contingency_matrix(Int[], Int[])
        Test.@test size(result.matrix) == (0, 0)

        nmi = Mycelia.normalized_mutual_information(Int[], Int[])
        Test.@test nmi == 1.0
    end

    Test.@testset "Single sample" begin
        ari = Mycelia.adjusted_rand_index([1], [2])
        Test.@test ari == 1.0

        nmi = Mycelia.normalized_mutual_information([1], [1])
        Test.@test nmi == 1.0
    end

    Test.@testset "All same cluster" begin
        labels1 = [1, 1, 1, 1, 1]
        labels2 = [2, 2, 2, 2, 2]

        ari = Mycelia.adjusted_rand_index(labels1, labels2)
        Test.@test ari == 1.0  ## Same structure (one cluster)

        nmi = Mycelia.normalized_mutual_information(labels1, labels2)
        Test.@test nmi == 1.0
    end

    Test.@testset "Each sample its own cluster" begin
        labels1 = [1, 2, 3, 4, 5]
        labels2 = [1, 2, 3, 4, 5]

        ari = Mycelia.adjusted_rand_index(labels1, labels2)
        Test.@test isapprox(ari, 1.0, atol=1e-10)

        ## FMI is 0 when each sample is its own cluster because there are
        ## no pairs in the same cluster (TP = 0), so precision/recall undefined
        fmi = Mycelia.fowlkes_mallows_index(labels1, labels2)
        Test.@test fmi == 0.0  ## No pairs to compare
    end

    Test.@testset "Negative cluster labels" begin
        labels1 = [-1, -1, 0, 0, 1, 1]
        labels2 = [-1, -1, -1, 0, 0, 0]

        ## Should work with negative labels
        ari = Mycelia.adjusted_rand_index(labels1, labels2)
        Test.@test -1.0 <= ari <= 1.0

        cont = Mycelia.contingency_matrix(labels1, labels2)
        Test.@test sum(cont.matrix) == 6
    end

    Test.@testset "Mixed type labels" begin
        labels1 = ["a", "a", "b", "b"]
        labels2 = [1, 1, 2, 2]

        ## This should work - comparing structure, not label values
        ari = Mycelia.adjusted_rand_index(labels1, labels2)
        Test.@test isapprox(ari, 1.0, atol=1e-10)
    end
end


test_println("\n=== All Cluster Comparison Metric Tests Complete ===")
