# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/momentum_fork_resolution_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/momentum_fork_resolution_test.jl", "test/4_assembly", execute=false)'
# ```

import Test
import Mycelia

Test.@testset "Momentum Fork Resolution" begin
    Test.@testset "Weight functions" begin
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.uniform_weight(0) == 0.0
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.uniform_weight(3) == 1.0
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.uniform_weight(
            3; unit_weight = 2.5) == 2.5

        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.linear_weight(3) == 3.0
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.linear_weight(
            4; slope = 1.5, intercept = 0.5) == 6.5
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.linear_weight(
            10; cap = 4.0) == 4.0

        saturating = Mycelia.Rhizomorph.MomentumForkResolver.saturating_weight(
            9; max_weight = 3.0, half_saturation = 3.0)
        Test.@test saturating ≈ 2.25
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.saturating_weight(0) == 0.0

        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.llr_weight(8, 2) > 0
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.llr_weight(2, 8) < 0
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.llr_weight(5, 5) ≈ 0.0
    end

    Test.@testset "SPRT thresholds" begin
        thresholds = Mycelia.Rhizomorph.MomentumForkResolver.SPRTThresholds(
            alpha = 0.05, beta = 0.1)
        Test.@test thresholds.accept_alternative > 0
        Test.@test thresholds.accept_reference < 0

        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.sprt_decision(
            thresholds.accept_alternative + 0.01; alpha = 0.05, beta = 0.1) == :accept_alternative
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.sprt_decision(
            thresholds.accept_reference - 0.01; alpha = 0.05, beta = 0.1) == :accept_reference
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.sprt_decision(
            0.0; alpha = 0.05, beta = 0.1) == :continue
    end

    Test.@testset "Single-branch weighted updates" begin
        state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
            "read-1", [:major, :minor]; reference_branch = :major)

        decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state, :minor, 2; weight_fn = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight)
        Test.@test decision == :continue
        Test.@test state.current_branch == :minor
        Test.@test state.branch_support[:minor] == 2.0
        Test.@test state.log_likelihood_ratio == 2.0

        decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state,
            :minor,
            5;
            weight_fn = Mycelia.Rhizomorph.MomentumForkResolver.saturating_weight,
            max_weight = 4.0,
            half_saturation = 1.0
        )
        Test.@test decision == :accept_alternative
        Test.@test state.resolved_branch == :minor
        Test.@test state.observations == 2
    end

    Test.@testset "Reference rescue" begin
        state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
            "read-2", [:major, :minor]; reference_branch = :major)

        Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state, :minor, 2; weight_fn = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight)
        decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state, :major, 7; weight_fn = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight)

        Test.@test decision == :accept_reference
        Test.@test state.current_branch == :major
        Test.@test state.resolved_branch == :major
        Test.@test state.log_likelihood_ratio < 0
    end

    Test.@testset "Direct LLR updates" begin
        state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
            "read-3", [:path_a, :path_b]; reference_branch = :path_a)

        decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state, :path_b, :path_a, 40, 1; alpha = 0.05, beta = 0.05)

        Test.@test decision == :accept_alternative
        Test.@test state.current_branch == :path_b
        Test.@test state.resolved_branch == :path_b
        Test.@test state.branch_support[:path_b] == 40.0
        Test.@test state.branch_support[:path_a] == 1.0

        resolution = Mycelia.Rhizomorph.MomentumForkResolver.fork_resolution(state)
        Test.@test resolution.decision == :accept_alternative
        Test.@test resolution.resolved_branch == :path_b
    end

    Test.@testset "Paired weighted updates" begin
        state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
            "read-3b", [:path_a, :path_b]; reference_branch = :path_a)

        decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_weighted_fork_event!(
            state,
            :path_a,
            1,
            :path_b,
            6;
            weight_fn = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight
        )

        Test.@test decision == :accept_alternative
        Test.@test state.observations == 1
        Test.@test state.branch_support[:path_a] == 1.0
        Test.@test state.branch_support[:path_b] == 6.0
        Test.@test state.current_branch == :path_b
        Test.@test state.resolved_branch == :path_b
    end

    Test.@testset "Direct LLR updates require reference branch" begin
        state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
            "read-3c", [:path_a, :path_b, :path_c]; reference_branch = :path_a)

        Test.@test_throws ErrorException Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state, :path_b, :path_c, 8, 2; alpha = 0.05, beta = 0.05)
    end

    Test.@testset "Best-branch helpers" begin
        state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
            "read-4", ["left", "right", "center"]; reference_branch = "center")
        Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state, "left", 3; weight_fn = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight)
        Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
            state, "right", 2; weight_fn = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight)

        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.best_branch(state) == "left"
        Test.@test Mycelia.Rhizomorph.MomentumForkResolver.strongest_alternative(state) == "left"
    end
end
