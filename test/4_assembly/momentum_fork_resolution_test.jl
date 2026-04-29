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
    function simulate_diploid_case(
            case_id::AbstractString,
            support_pairs::Vector{Tuple{Int, Int}};
            reference_branch::Symbol = :hap_a
    )
        state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
            case_id, [:hap_a, :hap_b]; reference_branch)
        decisions = Symbol[]

        for (reference_support, alternative_support) in support_pairs
            decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
                state,
                :hap_b,
                :hap_a,
                alternative_support,
                reference_support
            )
            push!(decisions, decision)
        end

        resolution = Mycelia.Rhizomorph.MomentumForkResolver.fork_resolution(state)
        return (
            case_id = case_id,
            decisions = decisions,
            resolution = resolution,
            preserved = resolution.resolved_branch === nothing,
            collapsed = resolution.resolved_branch !== nothing
        )
    end

    function summarize_diploid_cases(cases)
        results = [simulate_diploid_case(case_id, support_pairs) for (case_id, support_pairs) in cases]
        total_cases = length(results)
        collapsed_cases = count(result -> result.collapsed, results)
        preserved_cases = total_cases - collapsed_cases

        return (
            results = results,
            total_cases = total_cases,
            collapsed_cases = collapsed_cases,
            preserved_cases = preserved_cases,
            collapse_rate = total_cases == 0 ? 0.0 : collapsed_cases / total_cases,
            preservation_rate = total_cases == 0 ? 0.0 : preserved_cases / total_cases
        )
    end

    Test.@testset "Weight functions" begin
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

    Test.@testset "Diploid heterozygosity stress" begin
        preserved_cases = [
            ("perfect_balance", [(12, 12), (11, 11), (13, 13), (10, 10)]),
            ("alternating_microbias", [(13, 12), (12, 13), (13, 12), (12, 13), (13, 12), (12, 13)]),
            ("mild_reference_advantage", [(20, 15), (20, 15), (20, 15), (20, 15)]),
            ("mild_alternative_advantage", [(15, 20), (15, 20), (15, 20), (15, 20)])
        ]
        preserved_summary = summarize_diploid_cases(preserved_cases)

        Test.@test preserved_summary.total_cases == length(preserved_cases)
        Test.@test preserved_summary.collapse_rate == 0.0
        Test.@test preserved_summary.preservation_rate == 1.0
        for result in preserved_summary.results
            Test.@test result.preserved
            Test.@test !result.collapsed
            Test.@test result.resolution.decision == :continue
            Test.@test result.resolution.resolved_branch === nothing
        end

        collapsing_cases = [
            ("decisive_reference", [(40, 1)]),
            ("decisive_alternative", [(1, 40)]),
            ("compound_reference", [(20, 2), (20, 2)]),
            ("compound_alternative", [(2, 20), (2, 20)])
        ]
        collapsing_summary = summarize_diploid_cases(collapsing_cases)
        resolved_branches = Set(result.resolution.resolved_branch for result in collapsing_summary.results)

        Test.@test collapsing_summary.total_cases == length(collapsing_cases)
        Test.@test collapsing_summary.collapse_rate == 1.0
        Test.@test collapsing_summary.preservation_rate == 0.0
        Test.@test resolved_branches == Set([:hap_a, :hap_b])
        for result in collapsing_summary.results
            Test.@test !result.preserved
            Test.@test result.collapsed
            Test.@test result.resolution.decision in (:accept_reference, :accept_alternative)
            Test.@test result.resolution.resolved_branch in (:hap_a, :hap_b)
        end
    end
end
