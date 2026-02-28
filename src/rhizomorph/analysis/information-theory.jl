# Information-theoretic metrics for graph analysis.
#
# Provides entropy, divergence, and frequency distribution functions
# for analyzing Rhizomorph graphs and their transition structures.

"""
    shannon_entropy(probs::Vector{Float64})

Compute Shannon entropy H(X) = -Σ p(x) log₂ p(x) of a probability distribution.

Zero-probability entries are safely skipped (0 log 0 = 0 by convention).

# Returns
- `Float64`: Entropy in bits
"""
function shannon_entropy(probs::Vector{Float64})
    return -sum(p > 0 ? p * log2(p) : 0.0 for p in probs)
end

"""
    vertex_transition_entropy(graph::MetaGraphsNext.MetaGraph, vertex_label)

Compute Shannon entropy of the transition probability distribution from a vertex.
Higher entropy means more uniform (uncertain) transitions.

# Returns
- `Float64`: Entropy in bits. Returns 0.0 for vertices with no outgoing edges.
"""
function vertex_transition_entropy(graph::MetaGraphsNext.MetaGraph, vertex_label)
    transitions = _get_valid_transitions(graph, vertex_label, Forward)
    if isempty(transitions)
        return 0.0
    end
    probs = _calculate_transition_probabilities(transitions)
    return shannon_entropy(probs)
end

"""
    graph_transition_entropy(graph::MetaGraphsNext.MetaGraph)

Compute mean transition entropy across all vertices in the graph.
A measure of the overall predictability of random walks.

Vertices with no outgoing edges are excluded from the mean.

# Returns
- `Float64`: Mean entropy in bits. Returns 0.0 for empty graphs.
"""
function graph_transition_entropy(graph::MetaGraphsNext.MetaGraph)
    labels = collect(MetaGraphsNext.labels(graph))
    if isempty(labels)
        return 0.0
    end

    entropies = Float64[]
    for label in labels
        transitions = _get_valid_transitions(graph, label, Forward)
        if !isempty(transitions)
            push!(entropies, vertex_transition_entropy(graph, label))
        end
    end

    return isempty(entropies) ? 0.0 : Statistics.mean(entropies)
end

"""
    jensen_shannon_divergence(p::Dict, q::Dict)

Compute Jensen-Shannon divergence between two distributions.
JSD(P||Q) = 0.5 * KL(P||M) + 0.5 * KL(Q||M), where M = 0.5*(P+Q).

Both dictionaries map keys to probabilities. Keys present in one but not
the other are treated as having zero probability in the missing distribution.

# Returns
- `Float64`: JSD value in [0, 1] (in bits)
"""
function jensen_shannon_divergence(p::Dict, q::Dict)
    all_keys = union(keys(p), keys(q))
    if isempty(all_keys)
        return 0.0
    end

    jsd = 0.0
    for k in all_keys
        pk = get(p, k, 0.0)
        qk = get(q, k, 0.0)
        mk = 0.5 * (pk + qk)
        if mk > 0
            if pk > 0
                jsd += 0.5 * pk * log2(pk / mk)
            end
            if qk > 0
                jsd += 0.5 * qk * log2(qk / mk)
            end
        end
    end

    return max(0.0, jsd)
end

"""
    kl_divergence(p::Dict, q::Dict)

Compute KL divergence D_KL(P || Q) = Σ p(x) log₂(p(x) / q(x)).

Returns `Inf` if any key with nonzero probability in P has zero probability in Q.

# Returns
- `Float64`: KL divergence in bits (non-negative, possibly Inf)
"""
function kl_divergence(p::Dict, q::Dict)
    result = 0.0
    for (k, pk) in p
        if pk > 0
            qk = get(q, k, 0.0)
            if qk == 0
                return Inf
            end
            result += pk * log2(pk / qk)
        end
    end
    return max(0.0, result)
end

"""
    estimate_zipf_exponent(frequencies::Vector{<:Number})

Estimate the Zipf exponent α from a frequency vector using log-log linear regression.
Zipf's law: f(r) ∝ r^(-α), so log(f) = -α * log(r) + c.

# Returns
- `Float64`: Estimated Zipf exponent (positive value; ~1.0 for natural language)
- Returns 0.0 if fewer than 2 nonzero frequencies
"""
function estimate_zipf_exponent(frequencies::Vector{<:Number})
    sorted = sort(filter(x -> x > 0, frequencies); rev = true)
    if length(sorted) < 2
        return 0.0
    end

    ranks = collect(1:length(sorted))
    log_ranks = log.(Float64.(ranks))
    log_freqs = log.(Float64.(sorted))

    slope = Statistics.cov(log_ranks, log_freqs) / Statistics.var(log_ranks)
    return -slope
end

"""
    vertex_label_frequencies(graph::MetaGraphsNext.MetaGraph)

Compute observation counts per vertex in the graph. For evidence-based graphs,
returns the evidence count per vertex. For weighted graphs (no evidence data),
returns 1 per vertex.

# Returns
- `Dict{Any, Int}`: Map from vertex label to observation count
"""
function vertex_label_frequencies(graph::MetaGraphsNext.MetaGraph)
    freqs = Dict{Any, Int}()
    for label in MetaGraphsNext.labels(graph)
        vdata = graph[label]
        count = try
            count_evidence(vdata)
        catch
            1
        end
        freqs[label] = count
    end
    return freqs
end

"""
    edge_weight_distribution(graph::MetaGraphsNext.MetaGraph)

Collect all edge weights from the graph.

For weighted graphs (StrandWeightedEdgeData), extracts the `.weight` field.
For evidence graphs, computes weight via `count_evidence`.

# Returns
- `Vector{Float64}`: Edge weights
"""
function edge_weight_distribution(graph::MetaGraphsNext.MetaGraph)
    weights = Float64[]
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        edata = graph[src, dst]
        w = if edata isa StrandWeightedEdgeData
            edata.weight
        else
            try
                Float64(count_evidence(edata))
            catch
                1.0
            end
        end
        push!(weights, w)
    end
    return weights
end
