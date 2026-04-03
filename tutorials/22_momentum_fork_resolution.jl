# # Tutorial 22: Momentum Fork Resolution
#
# This tutorial demonstrates how to use the `MomentumForkResolver` API to
# accumulate read evidence at assembly forks and make an SPRT-based branch
# decision.

# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/22_momentum_fork_resolution.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia

println("=== Momentum Fork Resolution Tutorial ===")

# ## 1. Define a fork and initialize the active-read state
#
# The reference branch is the null hypothesis (`H0`). Competing branches become
# the alternative hypothesis (`H1`) whenever they accumulate more evidence than
# the reference.

state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
    "read_001",
    [:reference_path, :variant_path];
    reference_branch = :reference_path
)

println("Initial state:")
println("  reference branch: $(state.reference_branch)")
println("  current branch: $(state.current_branch)")
println("  initial LLR: $(state.log_likelihood_ratio)")

# ## 2. Compare the available weighting functions

uniform = Mycelia.Rhizomorph.MomentumForkResolver.uniform_weight(6)
linear = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight(6; slope = 1.0)
saturating = Mycelia.Rhizomorph.MomentumForkResolver.saturating_weight(
    6; max_weight = 3.0, half_saturation = 2.0)
llr = Mycelia.Rhizomorph.MomentumForkResolver.llr_weight(12, 3)

println("\nExample weights:")
println("  uniform weight for support=6: $uniform")
println("  linear weight for support=6: $linear")
println("  saturating weight for support=6: $(round(saturating, digits=3))")
println("  direct LLR for support 12 vs 3: $(round(llr, digits=3))")

# ## 3. Accumulate sequential evidence with a saturating response
#
# Saturating weights are useful when repeated support should still help, but
# each extra read should contribute a little less than the previous one.

for support in (2, 3, 4)
    decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
        state,
        :variant_path,
        support;
        weight_fn = Mycelia.Rhizomorph.MomentumForkResolver.saturating_weight,
        max_weight = 2.5,
        half_saturation = 1.5,
        alpha = 0.05,
        beta = 0.05
    )
    println(
        "  update with support=$(support): branch=$(state.current_branch), " *
        "llr=$(round(state.log_likelihood_ratio, digits=3)), decision=$(decision)"
    )
end

# ## 4. Inject a direct log-likelihood-ratio update
#
# When evidence already arrives as a branch-vs-branch comparison, apply it
# directly instead of forcing it through a univariate weight curve.

decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
    state,
    :variant_path,
    :reference_path,
    50,
    2;
    alpha = 0.05,
    beta = 0.05
)

resolution = Mycelia.Rhizomorph.MomentumForkResolver.fork_resolution(state)

println("\nDirect LLR update:")
println("  decision: $(decision)")
println("  resolved branch: $(resolution.resolved_branch)")
println("  current branch: $(resolution.current_branch)")
println("  cumulative LLR: $(round(resolution.log_likelihood_ratio, digits=3))")
println("  accept alternative at: $(round(resolution.thresholds.accept_alternative, digits=3))")
println("  accept reference at: $(round(resolution.thresholds.accept_reference, digits=3))")

# ## 5. Interpretation
#
# A positive LLR favors the strongest competing branch over the reference.
# Once the upper SPRT threshold is crossed, the fork can be resolved in favor of
# the competing branch. Strong negative evidence would instead drive the state
# below the lower threshold and resolve back to the reference branch.
