# Test for the EXHAUSTIVE hard-biological-repeat prime(coprime) vs
# composite(factor-sharing) k ablation (bead td-tjym).
#
# Covers tandem / microsatellite (STR) / satellite / nested / interspersed /
# immune-VDJ / MHC-polymorphic / antigenic-cassette / SINE-Alu / LINE-truncated /
# mixed. Encodes the HONEST, POWERED negative. It asserts:
#   (1) every hard repeat class sits BELOW the resolution ceiling (powered per
#       class, not a saturation artifact);
#   (2) the PRIMARY resolution metric has dynamic range on tandems while the
#       w-mer-recovery metric saturates (documenting why the latter is blind);
#   (3) the size-controlled ISOLATING COMPARISON — no factor-sharing k
#       (gcd(k,period)>1) robustly dips below the coprime-k envelope (interpolated
#       between coprime k below & above) AND below the lower-size coprime k, across
#       all classes incl. the immune-VDJ and SINE internal-period wildcards; and
#   (4) the HARNESS-SENSITIVITY collision control fires — the metric CAN detect an
#       aliasing deficit, so the factor-sharing null is real.

import DataFrames
import Test

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_prime_k_ablation_benchmark.jl"))

Test.@testset "rhizomorph prime(coprime)-vs-composite(factor-sharing) degree+class ablation" begin
    output_dir = mktempdir()
    artifacts = run_prime_k_ablation_benchmark(output_dir; write_plots = true)

    fixtures = ablation_fixtures()
    n_fixtures = length(fixtures)
    per_run = artifacts.per_run_table

    Test.@testset "structure" begin
        Test.@test artifacts.per_run_rows == n_fixtures * length(ABLATION_K)
        Test.@test isfile(artifacts.per_run_csv)
        Test.@test isfile(artifacts.factor_sharing_csv)
        Test.@test isfile(artifacts.class_summary_csv)
        Test.@test isfile(artifacts.collision_control_csv)
        Test.@test isfile(artifacts.figure_png)
        Test.@test isfile(artifacts.biological_figure_png)
    end

    Test.@testset "exhaustive hard-repeat coverage present" begin
        # Every named hard-repeat class the manuscript claims to sweep is present.
        classes = Set(per_run.category)
        for expected in ("tandem", "microsatellite", "satellite", "nested", "interspersed",
            "immune_vdj", "immune_mhc", "immune_cassette", "sine_alu", "line", "mixed")
            Test.@test expected in classes
        end
        # The immune-VDJ (period 27) and SINE-Alu internal period-5 wildcard both
        # enter the factor-sharing comparison (a composite k aligning with each
        # period is in the grid).
        fs = artifacts.factor_sharing_table
        Test.@test "immune_vdj_segment_array" in Set(fs.repeat_class)
        Test.@test "sine_alu_interspersed" in Set(fs.repeat_class)
        Test.@test any((fs.repeat_class .== "sine_alu_interspersed") .&
                       (fs.composite_k .== 25))
    end

    Test.@testset "every hard class is below the resolution ceiling (powered)" begin
        cs = artifacts.class_summary_table
        # Non-control, non-collision classes must sit below the ceiling so the null
        # is a powered negative for EACH class (recovery is visibly below 1.0).
        hard = cs[(cs.category .!= "control") .& (cs.category .!= "collision_shared") .& (cs.category .!= "collision_distinct"), :]
        Test.@test all(hard.below_resolution_ceiling)
        Test.@test all(hard.max_resolution .< 0.95)
    end

    Test.@testset "metric ranges" begin
        Test.@test all(0.0 .<= per_run.mean_wmer_recovery .<= 1.0)
        Test.@test all(0.0 .<= per_run.largest_correct_contig_fraction .<= 1.0)
        Test.@test Set(unique(per_run.k)) == Set(ABLATION_K)
        for group in DataFrames.groupby(per_run, :dataset_id)
            Test.@test DataFrames.nrow(group) == length(ABLATION_K)
        end
    end

    Test.@testset "resolution metric is powered (NOT saturated on tandems)" begin
        tandem = per_run[per_run.category .== "tandem", :]
        resolution = tandem.largest_correct_contig_fraction
        # Dynamic range: tandems never resolve into one contig -> stays well below
        # the ceiling, so a prime benefit (if any) COULD appear. This is the fix
        # for the v3 saturation flaw.
        Test.@test maximum(resolution) < 0.90
        Test.@test (maximum(resolution) - minimum(resolution)) > 0.15
    end

    Test.@testset "w-mer recovery saturates above threshold (why it is blind)" begin
        # Documents the contrast: the presence-only metric hits the ceiling, so it
        # cannot detect the resolution deficits the primary metric can.
        top = per_run[
        (per_run.category .== "tandem") .& (per_run.k .== maximum(ABLATION_K)), :]
        Test.@test all(top.mean_wmer_recovery .>= 0.95)
    end

    Test.@testset "no isolated factor-sharing deficit (size-controlled, powered null)" begin
        fs = artifacts.factor_sharing_table
        Test.@test DataFrames.nrow(fs) > 0
        Test.@test all(fs.composite_gcd .> 1)                       # every row IS factor-sharing
        Test.@test all(fs.coprime_k_below .< fs.composite_k .< fs.coprime_k_above)
        # robust_striking requires: penalty>=margin AND >seed-noise AND a genuine
        # dip below the lower-size coprime k (rejects the SINE-k=15 interpolation
        # artifact where res sits on a steep rising ramp yet above its lower
        # coprime neighbour). No cell clears all three.
        Test.@test !any(fs.robust_striking)
        Test.@test artifacts.striking_prime_regime_found == false
    end

    Test.@testset "harness-sensitivity collision control fires" begin
        cc = artifacts.collision_control_table
        penalty_9 = only(cc[cc.k .== 9, :collision_penalty])
        penalty_11 = only(cc[cc.k .== 11, :collision_penalty])
        Test.@test penalty_9 >= COLLISION_MIN_GAP
        Test.@test (penalty_9 - penalty_11) >= COLLISION_MIN_GAP
        Test.@test artifacts.collision_control_fires
    end

    # Report the largest factor-sharing penalty and the verdict (informational).
    fs = artifacts.factor_sharing_table
    worst = fs[argmax(fs.factor_sharing_penalty), :]
    println("\nLargest factor-sharing penalty (below coprime envelope): " *
            "$(worst.repeat_class) p=$(worst.period) k=$(worst.composite_k) " *
            "penalty=$(round(worst.factor_sharing_penalty; digits = 3)) " *
            "(margin $(STRIKING_PRIME_ADVANTAGE), seed_noise $(round(worst.seed_noise; digits = 3)))")
    println("Verdict — striking prime-k regime found: $(artifacts.striking_prime_regime_found) | " *
            "collision control fires: $(artifacts.collision_control_fires) | " *
            "resolution range on tandems: [" *
            "$(round(minimum(per_run[per_run.category .== "tandem", :largest_correct_contig_fraction]); digits = 3)), " *
            "$(round(maximum(per_run[per_run.category .== "tandem", :largest_correct_contig_fraction]); digits = 3))]")
end
