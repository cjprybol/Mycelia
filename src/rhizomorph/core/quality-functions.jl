# Quality Score Functions for Qualmer Graphs
#
# Functions for computing joint quality scores and statistics
# from multiple observations of the same k-mer.
#
# Uses existing Mycelia quality functions from qualmer-analysis.jl

"""
    combine_phred_scores(scores::Vector{UInt8})

Combine multiple Phred quality scores from independent observations.

Phred scores are additive in log space:
- Phred10 + Phred10 = Phred20
- Error probabilities multiply: p_combined = p1 * p2 * ... * pn
- In Phred space: Q_combined = -10 * log10(∏ p_i) = ∑ Q_i

Combined scores can exceed 60 (the typical single-observation limit) and go up to 255.
This allows differentiation between:
- Low confidence: 10-30 (few low-quality observations)
- Moderate: 30-60 (typical single high-quality observation)
- High: 60-100 (multiple high-quality observations)
- Very high: 100+ (many high-quality observations)

# Arguments
- `scores::Vector{UInt8}`: Vector of Phred scores (NOT Phred+33 encoded, raw 0-60)

# Returns
- `UInt8`: Combined Phred score (clamped to 0-255, the UInt8 range)

# Examples
```julia
# Two Q10 observations combine to Q20
combined = combine_phred_scores(UInt8[10, 10])  # UInt8(20)

# Three Q30 observations combine to Q90
combined = combine_phred_scores(UInt8[30, 30, 30])  # UInt8(90)

# Ten Q20 observations combine to Q200
combined = combine_phred_scores(UInt8[20, 20, 20, 20, 20, 20, 20, 20, 20, 20])  # UInt8(200)
```
"""
function combine_phred_scores(scores::Vector{UInt8})
    if isempty(scores)
        return UInt8(0)
    end

    # Phred scores add in log space
    total = sum(Int(s) for s in scores)

    # Clamp to UInt8 range [0, 255]
    return UInt8(clamp(total, 0, 255))
end

# ASSUMPTION NOTICE (audit td-8zle; open question td-2tg8).
#
# Summing Phred is `aggregate_quality_scores_independence` from
# `planning-docs/rhizomorph-graph-ecosystem-plan.md:498-550`. Summing in log space is
# multiplying error probabilities, which is valid ONLY if the observations are
# independent.
#
# That assumption is load-bearing and unvalidated here. Correlated error sources —
# systematic/technical artifacts, strand, tile/flowcell, PCR duplicates, homopolymer
# context — make repeated observations less informative than independent ones, so this
# OVERSTATES confidence, and it does so fastest exactly where artifacts recur, which is
# the case quality-aware assembly exists to catch. The 0-255 clamp hides how quickly
# that saturates.
#
# The design anticipates this and specifies an alternative,
# `aggregate_quality_scores_conservative` — P(all correct) = prod(P_i), then
# Q = -10*log10(1 - P) — for when observations may share systematic errors. That
# alternative has NO implementation anywhere in `src/`; it exists only in the design
# document. Choosing between the two (or a middle model: capped independence,
# effective-sample-size shrinkage, or a per-position beta-binomial posterior) is
# td-2tg8, and it needs empirical calibration evidence, not a code change.
#
# Until then: this is the documented default, and anything reporting a joint quality
# should say which model produced it.

"""
    decode_quality_scores(scores::AbstractVector{UInt8})

Decode Phred+33 encoded quality scores into raw Phred values.

# Arguments
- `scores`: Vector of Phred+33 encoded quality scores

# Returns
- `Vector{Float64}`: Raw Phred scores as floating-point values
"""
function decode_quality_scores(scores::AbstractVector{UInt8})
    return Float64.(scores) .- 33.0
end

"""
    mean_joint_quality(graph, vertex_label, dataset_id::String)

Compute the mean joint quality score for a vertex in a qualmer graph.
Returns 0.0 when no quality evidence is available.
"""
function mean_joint_quality(graph::MetaGraphsNext.MetaGraph, vertex_label, dataset_id::String)
    joint = get_vertex_joint_quality(graph[vertex_label], dataset_id)
    return isnothing(joint) ? 0.0 : Statistics.mean(Float64.(joint))
end

"""
    mean_path_quality(graph, path, dataset_id::String)

Compute the mean joint quality score across a path of vertex labels.

!!! note "Unused by the core assembly paths (audit td-8zle, 2026-08)"
    No callers in `src/` and no test coverage. Kept as API rather than deleted, but
    it is NOT evidence that path quality informs any assembly decision — nothing
    calls it.

    Note also the double average: this means the joint quality across positions,
    then across path vertices, so a single low-quality vertex is diluted by the
    path's length. A path-level gate would want the WEAKEST vertex, mirroring
    `find_high_quality_sequences` (variable-length/fastq-graphs.jl:405) and
    `_qualmer_joint_confidence` (rhizomorph/assembly.jl), both of which require
    every position to clear the floor.
"""
function mean_path_quality(graph::MetaGraphsNext.MetaGraph, path, dataset_id::String)
    scores = Float64[]
    for vertex_label in path
        joint = get_vertex_joint_quality(graph[vertex_label], dataset_id)
        if !isnothing(joint)
            push!(scores, Statistics.mean(Float64.(joint)))
        end
    end
    return isempty(scores) ? 0.0 : Statistics.mean(scores)
end

"""
    get_vertex_joint_quality(vertex_data, dataset_id::String)

Compute joint quality scores for a k-mer vertex from all observations in a dataset.

For each position in the k-mer, combines quality scores across all observations
by adding Phred scores (which is equivalent to multiplying error probabilities).

# Arguments
- `vertex_data`: QualmerVertexData containing quality evidence
- `dataset_id::String`: Dataset identifier

# Returns
- `Vector{UInt8}`: Joint quality scores (raw Phred 0-255), or `nothing` if no quality data

# Examples
```julia
# Get joint quality for a k-mer across all reads in a dataset
joint_qual = get_vertex_joint_quality(vertex_data, "dataset_01")
```
"""
function get_vertex_joint_quality(vertex_data, dataset_id::String)
    # Get dataset evidence
    dataset_evidence = get_dataset_evidence(vertex_data, dataset_id)
    if isnothing(dataset_evidence)
        return nothing
    end

    # Collect all quality scores from all observations
    quality_vectors = Vector{Vector{UInt8}}()

    for observation_evidence in values(dataset_evidence)
        for evidence_entry in observation_evidence
            if evidence_entry isa QualityEvidenceEntry
                push!(quality_vectors, evidence_entry.quality_scores)
            end
        end
    end

    if isempty(quality_vectors)
        return nothing
    end

    # All observations should have same k-mer length
    k = length(quality_vectors[1])

    # Compute joint quality for each position
    joint_scores = Vector{UInt8}(undef, k)

    for pos in 1:k
        # Get all quality scores at this position (Phred+33 encoded)
        position_scores_encoded = [obs[pos] for obs in quality_vectors]

        # Decode from Phred+33 to raw Phred
        position_scores = [s - UInt8(33) for s in position_scores_encoded]

        # Combine (add in log space)
        joint_scores[pos] = combine_phred_scores(position_scores)
    end

    return joint_scores
end

"""
    get_vertex_mean_quality(vertex_data, dataset_id::String)

Compute mean quality score for each position from all observations in a dataset.

# Arguments
- `vertex_data`: QualmerVertexData containing quality evidence
- `dataset_id::String`: Dataset identifier

# Returns
- `Vector{Float64}`: Mean quality scores (raw Phred scale 0-60), or `nothing` if no quality data
"""
function get_vertex_mean_quality(vertex_data, dataset_id::String)
    dataset_evidence = get_dataset_evidence(vertex_data, dataset_id)
    if isnothing(dataset_evidence)
        return nothing
    end

    # Collect all quality scores
    quality_vectors = Vector{Vector{UInt8}}()
    for observation_evidence in values(dataset_evidence)
        for evidence_entry in observation_evidence
            if evidence_entry isa QualityEvidenceEntry
                push!(quality_vectors, evidence_entry.quality_scores)
            end
        end
    end

    if isempty(quality_vectors)
        return nothing
    end

    k = length(quality_vectors[1])
    mean_scores = Vector{Float64}(undef, k)

    for pos in 1:k
        # Decode from Phred+33, compute mean
        phred_scores = [Float64(obs[pos] - 33) for obs in quality_vectors]
        mean_scores[pos] = sum(phred_scores) / length(phred_scores)
    end

    return mean_scores
end

"""
    get_vertex_min_quality(vertex_data, dataset_id::String)

Get the minimum quality score at each position from all observations.

Useful for conservative quality thresholding.

!!! note "Unused by the core assembly paths (audit td-8zle, 2026-08)"
    The only caller in `src/` is `filter_by_quality` below, which is itself unused,
    so this is transitively dead in production.

    Read the semantics carefully before wiring it up: this takes the minimum ACROSS
    OBSERVATIONS per position, which DISCARDS the compounding that
    `get_vertex_joint_quality` provides. Under the design's model
    (`planning-docs/rhizomorph-graph-ecosystem-plan.md:498-550`) repeated
    high-quality observation of a position should INCREASE confidence; a per-position
    minimum instead lets one bad read veto a hundred good ones. That is the
    conservative posture, and it may be the right one if observations are correlated
    (see td-2tg8) — but it is a different model from the documented default, not an
    implementation of it.

# Arguments
- `vertex_data`: QualmerVertexData containing quality evidence
- `dataset_id::String`: Dataset identifier

# Returns
- `Vector{UInt8}`: Minimum quality scores (raw Phred 0-60), or `nothing` if no quality data
"""
function get_vertex_min_quality(vertex_data, dataset_id::String)
    dataset_evidence = get_dataset_evidence(vertex_data, dataset_id)
    if isnothing(dataset_evidence)
        return nothing
    end

    quality_vectors = Vector{Vector{UInt8}}()
    for observation_evidence in values(dataset_evidence)
        for evidence_entry in observation_evidence
            if evidence_entry isa QualityEvidenceEntry
                push!(quality_vectors, evidence_entry.quality_scores)
            end
        end
    end

    if isempty(quality_vectors)
        return nothing
    end

    k = length(quality_vectors[1])
    min_scores = Vector{UInt8}(undef, k)

    for pos in 1:k
        # Decode from Phred+33, find min
        phred_scores = [obs[pos] - UInt8(33) for obs in quality_vectors]
        min_scores[pos] = minimum(phred_scores)
    end

    return min_scores
end

"""
    filter_by_quality(graph, min_phred::Int, dataset_id::String)

Filter graph vertices by minimum quality score.

Returns k-mers where ALL positions have quality >= `min_phred` in at
least one observation from the specified dataset.

# Arguments
- `graph`: Qualmer graph
- `min_phred::Int`: Minimum Phred quality score (0-60, before +33 encoding)
- `dataset_id::String`: Dataset to filter on

# Returns
- `Vector`: K-mers meeting quality threshold

# Examples
```julia
# Get high-quality k-mers (Q30+)
high_qual_kmers = filter_by_quality(graph, 30, "dataset_01")
```

!!! note "Unused by the core assembly paths (audit td-8zle, 2026-08)"
    No callers in `src/`; exercised only by
    `test/4_assembly/quality_functions_test.jl` and
    `test/3_feature_extraction_kmer/qualmer-analysis.jl`.

    This is the natural quality gate for graph traversal, and its absence is
    precisely why the greedy qualmer arm was quality-blind (td-4e19d): a
    quality-filtering primitive existed, was tested, and was never called, so the
    arm's traversal stayed count-based while the framework was described as
    Phred-weighted. The opt-in `traversal_weighting = :quality` gate in
    `rhizomorph/assembly.jl` now does this job, but via
    `get_vertex_joint_quality` rather than this function, because it wants the
    compounding model rather than the per-position minimum this inherits from
    `get_vertex_min_quality`.
"""
function filter_by_quality(graph, min_phred::Int, dataset_id::String)
    filtered_kmers = []

    for kmer in MetaGraphsNext.labels(graph)
        vertex_data = graph[kmer]
        min_qual = get_vertex_min_quality(vertex_data, dataset_id)

        if !isnothing(min_qual) && all(q >= min_phred for q in min_qual)
            push!(filtered_kmers, kmer)
        end
    end

    return filtered_kmers
end

# -----------------------------------------------------------------------------
# H-11 P1: molecule collapse + Float64 posterior log-odds (additive APIs)
#
# These do NOT change `get_vertex_joint_quality` / `combine_phred_scores` (the
# uncalibrated `:quality` UInt8 path). Collapse keys by fragment id and reduces
# within-molecule entries by the mean of per-entry error probabilities (E5).
# -----------------------------------------------------------------------------

"""
Trailing mate / read-pair suffixes stripped to form a fragment observation key.

Matches (at end of string only): `/1`, `/2`, `.1`, `.2`, `_R1`, `_R2`.
"""
const FRAGMENT_OBSERVATION_ID_SUFFIX = r"(?:/[12]|\.[12]|_R[12])$"

"""
Default chemistry prior error probability for posterior log-odds (`≈ Q30`).

P2/P3 may specialize this per chemistry; P1 exposes it as a keyword default.
"""
const DEFAULT_POSTERIOR_PRIOR_ERROR = 1e-3

"""
Floor / ceiling for molecule error probabilities inside log-odds numerics.

Keeps `log(p)` / `log1p(-p)` finite so `Λ` is never ±Inf or NaN for finite Phred.
"""
const POSTERIOR_ERROR_PROB_EPS = 1e-300

"""
    fragment_observation_id(observation_id::AbstractString) -> String

Map a read / mate identifier to its fragment key by stripping a trailing
`/1`, `/2`, `.1`, `.2`, `_R1`, or `_R2` suffix.

Paired-end mates that share a fragment therefore collapse to one molecule under
H-11 E5 (see `collapse_molecule_qualities` / `get_vertex_posterior_log_odds`).
"""
function fragment_observation_id(observation_id::AbstractString)::String
    return replace(String(observation_id), FRAGMENT_OBSERVATION_ID_SUFFIX => "")
end

"""
    phred_error_probability(raw_phred::Real) -> Float64

Convert a raw (not Phred+33) Phred score to an error probability `10^(-Q/10)`.
"""
function phred_error_probability(raw_phred::Real)::Float64
    return 10.0^(-Float64(raw_phred) / 10.0)
end

"""
    collapse_molecule_qualities(quality_vectors) -> Vector{Float64}

Collapse multiple Phred+33-encoded quality vectors from **one molecule** to a
per-position mean error probability:

```
p_j = mean_e( 10^(-Q_{e,j}/10) )
```

where `Q` is the decoded (raw) Phred at position `j` for entry `e`.

# Arguments
- `quality_vectors`: Non-empty collection of equal-length Phred+33 `Vector{UInt8}`
  (as stored on `QualityEvidenceEntry.quality_scores`)

# Returns
- `Vector{Float64}`: length-k mean error probabilities (not Phred, not UInt8)

# Notes
- Empty input returns an empty `Float64` vector.
- This is the within-molecule reduction for H-11 E5. Across-molecule combination
  is done by `get_vertex_posterior_log_odds` / `_posterior_log_odds_from_molecule_errors`.
"""
function collapse_molecule_qualities(
        quality_vectors::AbstractVector{<:AbstractVector{UInt8}}
)::Vector{Float64}
    if isempty(quality_vectors)
        return Float64[]
    end

    k = length(quality_vectors[1])
    n = length(quality_vectors)
    mean_p = Vector{Float64}(undef, k)

    for pos in 1:k
        s = 0.0
        for obs in quality_vectors
            raw_phred = Float64(obs[pos]) - 33.0
            s += phred_error_probability(raw_phred)
        end
        mean_p[pos] = s / n
    end

    return mean_p
end

"""
    collapse_molecule_qualities(vertex_data, dataset_id) -> Union{Nothing, Vector{Vector{Float64}}}

Fragment-key all observation ids in `dataset_id`, collapse each molecule's
entries with mean-of-`p`, and return one length-k error-probability vector per
molecule (order is not meaningful; treat as a bag).

Returns `nothing` when the dataset is missing or has no quality evidence.
"""
function collapse_molecule_qualities(vertex_data, dataset_id::String)
    groups = _molecule_quality_groups(vertex_data, dataset_id)
    groups === nothing && return nothing

    molecule_ps = Vector{Vector{Float64}}()
    sizehint!(molecule_ps, length(groups))
    for quality_vectors in values(groups)
        push!(molecule_ps, collapse_molecule_qualities(quality_vectors))
    end
    return molecule_ps
end

"""
Group Phred+33 quality vectors by `fragment_observation_id(observation_id)`.

Returns `nothing` if the dataset is absent or no `QualityEvidenceEntry` exists.
"""
function _molecule_quality_groups(vertex_data, dataset_id::String)
    dataset_evidence = get_dataset_evidence(vertex_data, dataset_id)
    dataset_evidence === nothing && return nothing

    groups = Dict{String, Vector{Vector{UInt8}}}()
    for (observation_id, observation_evidence) in dataset_evidence
        frag = fragment_observation_id(observation_id)
        for evidence_entry in observation_evidence
            if evidence_entry isa QualityEvidenceEntry
                bucket = get!(groups, frag) do
                    Vector{Vector{UInt8}}()
                end
                push!(bucket, evidence_entry.quality_scores)
            end
        end
    end

    return isempty(groups) ? nothing : groups
end

"""
    _posterior_log_odds_from_molecule_errors(molecule_error_vectors; prior_error)

Per-position posterior log-odds that the consensus base is true, after molecule
collapse.

For each position `j` with molecule error probabilities `{p_m}`:

```
Λ_j = log((1-π)/π) + Σ_m [ log1p(-p_m) - log(p_m) ]
```

with each `p_m` clamped into `(ε, 1-ε)` so the result is finite. Accumulation is
log-space only (ADR §2.3). No UInt8 clamp.
"""
function _posterior_log_odds_from_molecule_errors(
        molecule_error_vectors::AbstractVector{<:AbstractVector{<:Real}};
        prior_error::Float64 = DEFAULT_POSTERIOR_PRIOR_ERROR
)::Vector{Float64}
    if isempty(molecule_error_vectors)
        return Float64[]
    end
    if !(0.0 < prior_error < 1.0)
        error("prior_error must be in (0, 1), got prior_error=$(prior_error)")
    end

    k = length(molecule_error_vectors[1])
    prior_log_odds = log1p(-prior_error) - log(prior_error)
    scores = fill(prior_log_odds, k)
    eps = POSTERIOR_ERROR_PROB_EPS

    for molecule_p in molecule_error_vectors
        length(molecule_p) == k || error(
            "molecule error vectors must share length k=$(k); got $(length(molecule_p))"
        )
        for pos in 1:k
            p = clamp(Float64(molecule_p[pos]), eps, 1.0 - eps)
            scores[pos] += log1p(-p) - log(p)
        end
    end

    return scores
end

"""
    get_vertex_posterior_log_odds_scores(vertex_data, dataset_id; prior_error) -> Union{Nothing, Vector{Float64}}

Per-position Float64 posterior log-odds `Λ_j` after fragment-keyed molecule
collapse (H-11 E5). Additive API — does not alter `get_vertex_joint_quality`.

# Arguments
- `vertex_data`: Qualmer (or other quality-evidence) vertex
- `dataset_id::String`: Dataset key
- `prior_error::Float64`: Chemistry prior P(error); default `DEFAULT_POSTERIOR_PRIOR_ERROR`

# Returns
- `Vector{Float64}` of length k, or `nothing` if no quality evidence
"""
function get_vertex_posterior_log_odds_scores(
        vertex_data,
        dataset_id::String;
        prior_error::Float64 = DEFAULT_POSTERIOR_PRIOR_ERROR
)
    molecule_ps = collapse_molecule_qualities(vertex_data, dataset_id)
    molecule_ps === nothing && return nothing
    return _posterior_log_odds_from_molecule_errors(molecule_ps; prior_error = prior_error)
end

"""
    get_vertex_posterior_log_odds(vertex_data, dataset_id; prior_error) -> Union{Nothing, Float64}

Scalar k-mer posterior log-odds `Λ(v)`: the **weakest** per-position score from
`get_vertex_posterior_log_odds_scores` (conjunction semantics matching
`_qualmer_joint_confidence` in `assembly.jl`).

Returns `nothing` when the vertex has no recoverable quality evidence.
"""
function get_vertex_posterior_log_odds(
        vertex_data,
        dataset_id::String;
        prior_error::Float64 = DEFAULT_POSTERIOR_PRIOR_ERROR
)
    scores = get_vertex_posterior_log_odds_scores(
        vertex_data, dataset_id; prior_error = prior_error
    )
    scores === nothing && return nothing
    isempty(scores) && return nothing
    return minimum(scores)
end
