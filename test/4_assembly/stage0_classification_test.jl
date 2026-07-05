# Unit tests for Stage 0 k-mer classification (td-4z2p).
# Deterministic tests of the classifier math independent of the stochastic
# benchmark, including regression guards for review findings C1 and C2.

import Test
import Mycelia
import Statistics

const R = Mycelia.Rhizomorph

Test.@testset "Stage 0 k-mer classification" begin

    Test.@testset "C1 guard: _quality_correct_prob never returns exactly 1.0" begin
        # Accumulated/extreme Phred must not collapse the (1 - q) error term.
        for q in (0.0, 20.0, 40.0, 60.0, 255.0, 1e6)
            p = R._quality_correct_prob(Float64(q))
            Test.@test 0.0 <= p <= 1.0 - 1e-6
            Test.@test 1.0 - p >= 1e-6            # error channel stays strictly positive
        end
        Test.@test R._quality_correct_prob(0.0) == 0.0
    end

    Test.@testset "C2 guard: _mean_phred handles nothing and empty" begin
        Test.@test R._mean_phred(nothing) == 0.0           # accessor returns nothing on no evidence
        Test.@test R._mean_phred(Float64[]) == 0.0
        Test.@test R._mean_phred([30.0, 40.0]) == 35.0
        Test.@test R._mean_phred(UInt8[20, 40]) == 30.0    # still accepts integer vectors
    end

    Test.@testset "_histogram_valley finds the bimodal minimum" begin
        # error peak at coverage 1, valley near 3, genomic peak near 8.
        hist = [50, 20, 5, 3, 8, 20, 40, 60, 40, 15]
        v = R._histogram_valley(hist, 2)
        Test.@test 2 <= v <= 5                              # valley sits between the two peaks
        Test.@test hist[v] <= hist[1]                       # below the error peak
        Test.@test hist[v] <= hist[8]                       # below the genomic peak
    end

    Test.@testset "fit_logistic_fusion recovers a separable boundary" begin
        # coverage cleanly separates the classes → learned coverage weight > 0,
        # and the fitted classifier reproduces the labels.
        covs = vcat(fill(1, 20), fill(30, 20))
        phreds = vcat(fill(15.0, 20), fill(35.0, 20))
        ys = vcat(fill(false, 20), fill(true, 20))
        clf = R.fit_logistic_fusion(covs, phreds, ys; iters = 2000)
        Test.@test clf.w_cov > 0
        correct = 0
        for i in eachindex(ys)
            post = R._posterior(clf, covs[i], phreds[i])
            (post >= 0.5) == ys[i] && (correct += 1)
        end
        Test.@test correct == length(ys)                    # perfectly separable → 100%
    end

    Test.@testset "_posterior(BayesianMixture) increases with coverage" begin
        clf = R.BayesianMixtureClassifier(genomic_mean_coverage = 30.0)
        p_low = R._posterior(clf, 1, 35.0)
        p_high = R._posterior(clf, 30, 35.0)
        Test.@test 0.0 <= p_low <= 1.0
        Test.@test 0.0 <= p_high <= 1.0
        Test.@test p_high > p_low                           # more coverage ⇒ more likely genomic
    end

    Test.@testset "BloomFilter has no false negatives" begin
        bf = R._BloomFilter(1 << 12, 4)
        inserted = ["ACGTA", "TTGCA", "GGGGG", "CATCA"]
        for x in inserted
            R._bloom_insert!(bf, x)
        end
        for x in inserted
            Test.@test R._bloom_contains(bf, x)             # membership is never missed
        end
    end
end
