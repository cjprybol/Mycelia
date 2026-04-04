"""
Sequential fork-resolution helpers for read-evidence assembly graphs.

The resolver models fork decisions as a sequential hypothesis test between a
reference branch (H0) and whichever competing branch currently has the highest
weighted evidence (H1). Evidence can be accumulated with linear or saturating
weighting, or injected directly as a log-likelihood ratio.
"""
module MomentumForkResolver

"""
    SPRTThresholds

Decision thresholds for Wald's sequential probability ratio test.

# Fields
- `alpha::Float64`: False-positive rate for accepting the alternative branch
- `beta::Float64`: False-negative rate for accepting the reference branch
- `accept_alternative::Float64`: Upper LLR threshold
- `accept_reference::Float64`: Lower LLR threshold

# Example
```julia
thresholds = Mycelia.Rhizomorph.MomentumForkResolver.SPRTThresholds(; alpha=0.05, beta=0.1)
```
"""
struct SPRTThresholds
    alpha::Float64
    beta::Float64
    accept_alternative::Float64
    accept_reference::Float64
end

"""
    SPRTThresholds(; alpha=0.05, beta=0.05)

Construct sequential probability ratio test thresholds from false-positive and
false-negative error budgets.

# Arguments
- `alpha::Real=0.05`: False-positive rate for accepting the alternative branch
- `beta::Real=0.05`: False-negative rate for accepting the reference branch

# Returns
- `SPRTThresholds`: Threshold container with precomputed acceptance boundaries

# Example
```julia
thresholds = Mycelia.Rhizomorph.MomentumForkResolver.SPRTThresholds(; alpha=0.05, beta=0.1)
```
"""
function SPRTThresholds(; alpha::Real = 0.05, beta::Real = 0.05)
    0.0 < alpha < 1.0 || error("alpha must satisfy 0 < alpha < 1, got $(alpha)")
    0.0 < beta < 1.0 || error("beta must satisfy 0 < beta < 1, got $(beta)")
    alpha64 = Float64(alpha)
    beta64 = Float64(beta)
    return SPRTThresholds(
        alpha64,
        beta64,
        log((1.0 - beta64) / alpha64),
        log(beta64 / (1.0 - alpha64))
    )
end

"""
    ActiveReadState(read_id, branches; reference_branch=first(branches))

Track sequential branch evidence for one read traversing a fork.

# Fields
- `read_id::String`: Read identifier
- `reference_branch`: Null-hypothesis branch used by the SPRT
- `current_branch`: Branch with the highest cumulative weighted evidence
- `branch_scores::Dict`: Per-branch weighted evidence totals
- `branch_support::Dict`: Per-branch raw support totals
- `log_likelihood_ratio::Float64`: Signed LLR relative to `reference_branch`
- `observations::Int`: Number of sequential updates applied
- `resolved_branch`: Final resolved branch when a threshold is crossed
- `decision::Symbol`: `:continue`, `:accept_reference`, or `:accept_alternative`

# Example
```julia
state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
    "read-1",
    [:major, :minor];
    reference_branch=:major,
)
```
"""
mutable struct ActiveReadState{B}
    read_id::String
    reference_branch::B
    current_branch::B
    branch_scores::Dict{B, Float64}
    branch_support::Dict{B, Float64}
    log_likelihood_ratio::Float64
    observations::Int
    resolved_branch::Union{Nothing, B}
    decision::Symbol
end

"""
    ActiveReadState(read_id, branches; reference_branch=first(branches))

Initialize per-read fork-resolution state with zeroed evidence totals for each
candidate branch.

# Arguments
- `read_id::AbstractString`: Read identifier used in diagnostics
- `branches::AbstractVector{B}`: Candidate branch labels tracked for the read
- `reference_branch::Union{Nothing, B}=nothing`: Branch treated as the null
  hypothesis. Defaults to the first unique entry in `branches`

# Returns
- `ActiveReadState{B}`: Mutable state object ready for sequential evidence updates

# Example
```julia
state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
    "read-2",
    [:hap_a, :hap_b];
    reference_branch=:hap_a,
)
```
"""
function ActiveReadState(
        read_id::AbstractString,
        branches::AbstractVector{B};
        reference_branch::Union{Nothing, B} = nothing
) where {B}
    isempty(branches) && error("ActiveReadState requires at least one branch label")
    unique_branches = unique(collect(branches))
    chosen_reference = isnothing(reference_branch) ? first(unique_branches) : reference_branch
    chosen_reference in unique_branches ||
        error("reference_branch must be one of the supplied branches")

    branch_scores = Dict{B, Float64}(branch => 0.0 for branch in unique_branches)
    branch_support = Dict{B, Float64}(branch => 0.0 for branch in unique_branches)
    return ActiveReadState{B}(
        String(read_id),
        chosen_reference,
        chosen_reference,
        branch_scores,
        branch_support,
        0.0,
        0,
        nothing,
        :continue
    )
end

"""
    linear_weight(support; slope=1.0, intercept=0.0, cap=nothing)

Convert raw support into an additive linear evidence weight.

# Arguments
- `support::Real`: Non-negative raw support for the observed branch
- `slope::Real=1.0`: Multiplicative scaling applied to `support`
- `intercept::Real=0.0`: Constant offset added before clamping
- `cap::Union{Nothing, Real}=nothing`: Optional upper bound on the resulting weight

# Returns
- `Float64`: Non-negative additive evidence weight

# Example
```julia
weight = Mycelia.Rhizomorph.MomentumForkResolver.linear_weight(
    4;
    slope=1.5,
    intercept=0.5,
)
```
"""
function linear_weight(
        support::Real;
        slope::Real = 1.0,
        intercept::Real = 0.0,
        cap::Union{Nothing, Real} = nothing
)
    support >= 0 || error("support must be non-negative, got $(support)")
    slope >= 0 || error("slope must be non-negative, got $(slope)")
    weight = Float64(slope) * Float64(support) + Float64(intercept)
    if !isnothing(cap)
        weight = min(weight, Float64(cap))
    end
    return max(0.0, weight)
end

"""
    saturating_weight(support; max_weight=1.0, half_saturation=1.0)

Convert raw support into an additive saturating evidence weight using a
Michaelis-Menten style response curve.

# Arguments
- `support::Real`: Non-negative raw support for the observed branch
- `max_weight::Real=1.0`: Asymptotic maximum additive evidence
- `half_saturation::Real=1.0`: Support value that yields half of `max_weight`

# Returns
- `Float64`: Saturating evidence contribution bounded by `max_weight`

# Example
```julia
weight = Mycelia.Rhizomorph.MomentumForkResolver.saturating_weight(
    9;
    max_weight=3.0,
    half_saturation=3.0,
)
```
"""
function saturating_weight(
        support::Real;
        max_weight::Real = 1.0,
        half_saturation::Real = 1.0
)
    support >= 0 || error("support must be non-negative, got $(support)")
    max_weight > 0 || error("max_weight must be positive, got $(max_weight)")
    half_saturation > 0 || error("half_saturation must be positive, got $(half_saturation)")
    support64 = Float64(support)
    if support64 == 0.0
        return 0.0
    end
    return Float64(max_weight) * support64 / (Float64(half_saturation) + support64)
end

"""
    llr_weight(favored_support, competing_support; pseudocount=0.5)

Compute a direct log-likelihood-ratio contribution for two competing branches.
Positive values favor the first branch, negative values favor the second.

# Arguments
- `favored_support::Real`: Non-negative support for the favored branch
- `competing_support::Real`: Non-negative support for the competing branch
- `pseudocount::Real=0.5`: Positive smoothing value added to both supports

# Returns
- `Float64`: Signed log-likelihood-ratio contribution

# Example
```julia
llr = Mycelia.Rhizomorph.MomentumForkResolver.llr_weight(8, 2; pseudocount=0.5)
```
"""
function llr_weight(
        favored_support::Real,
        competing_support::Real;
        pseudocount::Real = 0.5
)
    favored_support >= 0 || error("favored_support must be non-negative, got $(favored_support)")
    competing_support >= 0 || error("competing_support must be non-negative, got $(competing_support)")
    pseudocount > 0 || error("pseudocount must be positive, got $(pseudocount)")
    return log(
        (Float64(favored_support) + Float64(pseudocount)) /
        (Float64(competing_support) + Float64(pseudocount))
    )
end

"""
    sprt_decision(log_likelihood_ratio; alpha=0.05, beta=0.05)

Apply SPRT threshold logic to a signed log-likelihood ratio.

# Arguments
- `log_likelihood_ratio::Real`: Signed evidence balance between the strongest
  alternative branch and the reference branch
- `alpha::Real=0.05`: False-positive rate for accepting the alternative branch
- `beta::Real=0.05`: False-negative rate for accepting the reference branch

# Returns
- `Symbol`: One of `:accept_alternative`, `:accept_reference`, or `:continue`

# Example
```julia
decision = Mycelia.Rhizomorph.MomentumForkResolver.sprt_decision(
    0.0;
    alpha=0.05,
    beta=0.1,
)
```
"""
function sprt_decision(log_likelihood_ratio::Real; alpha::Real = 0.05, beta::Real = 0.05)
    thresholds = SPRTThresholds(; alpha, beta)
    llr = Float64(log_likelihood_ratio)
    if llr >= thresholds.accept_alternative
        return :accept_alternative
    elseif llr <= thresholds.accept_reference
        return :accept_reference
    end
    return :continue
end

"""
    best_branch(state)

Return the branch with the highest cumulative weighted evidence.

# Arguments
- `state::ActiveReadState`: Sequential fork-resolution state to inspect

# Returns
- Branch label: Branch with the largest weighted evidence total after tie-breaking

# Example
```julia
state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState("read-3", ["left", "right"])
branch = Mycelia.Rhizomorph.MomentumForkResolver.best_branch(state)
```
"""
function best_branch(state::ActiveReadState)
    ranked = _ranked_branches(state)
    return first(ranked)
end

"""
    strongest_alternative(state)

Return the strongest branch that is not the configured reference branch.

# Arguments
- `state::ActiveReadState`: Sequential fork-resolution state to inspect

# Returns
- Branch label: Highest-ranked non-reference branch, or the reference branch
  itself when no alternative exists

# Example
```julia
state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
    "read-4",
    [:center, :left, :right];
    reference_branch=:center,
)
branch = Mycelia.Rhizomorph.MomentumForkResolver.strongest_alternative(state)
```
"""
function strongest_alternative(state::ActiveReadState)
    for branch in _ranked_branches(state)
        if branch != state.reference_branch
            return branch
        end
    end
    return state.reference_branch
end

"""
    observe_fork_event!(state, branch, support=1; weight_fn=linear_weight, alpha=0.05, beta=0.05, weight_kwargs...)

Accumulate evidence for a single branch observation, then update the SPRT
decision state.

# Arguments
- `state::ActiveReadState{B}`: Mutable fork-resolution state
- `branch::B`: Observed branch to update
- `support::Real=1`: Non-negative support contributed by this observation
- `weight_fn=linear_weight`: Weighting callback applied to `support`
- `alpha::Real=0.05`: False-positive rate for accepting the alternative branch
- `beta::Real=0.05`: False-negative rate for accepting the reference branch
- `weight_kwargs...`: Additional keyword arguments forwarded to `weight_fn`

# Returns
- `Symbol`: Updated SPRT decision after incorporating the observation

# Example
```julia
state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
    "read-5",
    [:major, :minor];
    reference_branch=:major,
)
decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
    state,
    :minor,
    2;
    weight_fn=Mycelia.Rhizomorph.MomentumForkResolver.linear_weight,
)
```
"""
function observe_fork_event!(
        state::ActiveReadState{B},
        branch::B,
        support::Real = 1;
        weight_fn = linear_weight,
        alpha::Real = 0.05,
        beta::Real = 0.05,
        weight_kwargs...
) where {B}
    _require_known_branch(state, branch)
    support >= 0 || error("support must be non-negative, got $(support)")

    state.branch_support[branch] += Float64(support)
    state.branch_scores[branch] += Float64(weight_fn(support; weight_kwargs...))
    state.observations += 1
    _refresh_state!(state; alpha, beta)
    return state.decision
end

"""
    observe_fork_event!(state, branch_a, branch_b, support_a, support_b; alpha=0.05, beta=0.05, pseudocount=0.5)

Accumulate a competitive update between two branches using a direct
log-likelihood-ratio contribution.

# Arguments
- `state::ActiveReadState{B}`: Mutable fork-resolution state
- `branch_a::B`: Branch credited when the direct LLR is positive
- `branch_b::B`: Competing branch credited when the direct LLR is negative
- `support_a::Real`: Non-negative support for `branch_a`
- `support_b::Real`: Non-negative support for `branch_b`
- `alpha::Real=0.05`: False-positive rate for accepting the alternative branch
- `beta::Real=0.05`: False-negative rate for accepting the reference branch
- `pseudocount::Real=0.5`: Positive smoothing value for the direct LLR

# Returns
- `Symbol`: Updated SPRT decision after the competitive observation

# Example
```julia
state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState(
    "read-6",
    [:path_a, :path_b];
    reference_branch=:path_a,
)
decision = Mycelia.Rhizomorph.MomentumForkResolver.observe_fork_event!(
    state,
    :path_b,
    :path_a,
    40,
    1,
)
```
"""
function observe_fork_event!(
        state::ActiveReadState{B},
        branch_a::B,
        branch_b::B,
        support_a::Real,
        support_b::Real;
        alpha::Real = 0.05,
        beta::Real = 0.05,
        pseudocount::Real = 0.5
) where {B}
    _require_known_branch(state, branch_a)
    _require_known_branch(state, branch_b)
    support_a >= 0 || error("support_a must be non-negative, got $(support_a)")
    support_b >= 0 || error("support_b must be non-negative, got $(support_b)")

    direct_llr = llr_weight(support_a, support_b; pseudocount)
    state.branch_support[branch_a] += Float64(support_a)
    state.branch_support[branch_b] += Float64(support_b)
    state.branch_scores[branch_a] += direct_llr / 2
    state.branch_scores[branch_b] -= direct_llr / 2
    state.observations += 1
    _refresh_state!(state; alpha, beta)
    return state.decision
end

"""
    fork_resolution(state; alpha=0.05, beta=0.05)

Return the current resolution state as a named tuple.

# Arguments
- `state::ActiveReadState`: Sequential fork-resolution state to summarize
- `alpha::Real=0.05`: False-positive rate used to reconstruct the thresholds
- `beta::Real=0.05`: False-negative rate used to reconstruct the thresholds

# Returns
- Named tuple: Decision summary containing the resolved branch, current branch,
  reference branch, signed log-likelihood ratio, and `SPRTThresholds`

# Example
```julia
state = Mycelia.Rhizomorph.MomentumForkResolver.ActiveReadState("read-7", [:a, :b])
resolution = Mycelia.Rhizomorph.MomentumForkResolver.fork_resolution(state)
```
"""
function fork_resolution(state::ActiveReadState; alpha::Real = 0.05, beta::Real = 0.05)
    thresholds = SPRTThresholds(; alpha, beta)
    return (
        decision = state.decision,
        resolved_branch = state.resolved_branch,
        current_branch = state.current_branch,
        reference_branch = state.reference_branch,
        log_likelihood_ratio = state.log_likelihood_ratio,
        thresholds = thresholds
    )
end

function _ranked_branches(state::ActiveReadState)
    branches = collect(keys(state.branch_scores))
    sort!(
        branches;
        by = branch -> (
            state.branch_scores[branch],
            state.branch_support[branch],
            branch == state.reference_branch ? 1 : 0,
            string(branch)
        ),
        rev = true
    )
    return branches
end

function _refresh_state!(state::ActiveReadState; alpha::Real, beta::Real)
    state.current_branch = best_branch(state)
    reference_score = state.branch_scores[state.reference_branch]

    if state.current_branch == state.reference_branch
        alternative_branch = strongest_alternative(state)
        alternative_score = alternative_branch == state.reference_branch ? reference_score :
                            state.branch_scores[alternative_branch]
        state.log_likelihood_ratio = alternative_score - reference_score
    else
        state.log_likelihood_ratio = state.branch_scores[state.current_branch] - reference_score
    end

    state.decision = sprt_decision(state.log_likelihood_ratio; alpha, beta)
    if state.decision == :accept_alternative
        state.resolved_branch = state.current_branch
    elseif state.decision == :accept_reference
        state.resolved_branch = state.reference_branch
    else
        state.resolved_branch = nothing
    end

    return state
end

function _require_known_branch(state::ActiveReadState{B}, branch::B) where {B}
    haskey(state.branch_scores, branch) || error("Unknown branch $(branch) for read $(state.read_id)")
    return nothing
end

end
