# Graph metrics for Rhizomorph graphs.
#
# Migrated from panlingual-metagraphs notebooks 13-14.
# Provides eigenvector centrality, modularity, closeness centrality,
# and Betti number computation.

import LinearAlgebra

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute eigenvector centrality using power iteration.

Works on both `Graphs.AbstractGraph` and `MetaGraphsNext.MetaGraph` objects.
For MetaGraphs, operates on the underlying `.graph` field.

# Arguments
- `g`: Graph (Graphs.jl AbstractGraph or MetaGraphsNext.MetaGraph)
- `max_iter::Int=100`: Maximum iterations
- `tol::Float64=1e-6`: Convergence tolerance

# Returns
- `Vector{Float64}`: Centrality values per vertex (normalized)

# Example
```julia
centrality = Mycelia.Rhizomorph.compute_eigenvector_centrality(graph; max_iter = 50)
```
"""
function compute_eigenvector_centrality(g; max_iter::Int = 100, tol::Float64 = 1e-6)
    graph = g isa MetaGraphsNext.MetaGraph ? g.graph : g
    n = Graphs.nv(graph)
    if n == 0
        return Float64[]
    end

    x = fill(1.0 / sqrt(n), n)
    adj = Graphs.adjacency_matrix(graph)

    for _ in 1:max_iter
        x_new = adj * x
        norm_val = LinearAlgebra.norm(x_new)
        if norm_val > 0
            x_new ./= norm_val
        end
        if LinearAlgebra.norm(x_new - x) < tol
            break
        end
        x = x_new
    end

    return x
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute modularity Q using connected components as communities.
Q = (1/m) Σ_edges [δ(c_i, c_j) * (1 - k_i*k_j/(2m))]

Works on both `Graphs.AbstractGraph` and `MetaGraphsNext.MetaGraph` objects.

# Arguments
- `g`: Graph (Graphs.jl `AbstractGraph` or `MetaGraphsNext.MetaGraph`).

# Returns
- `Tuple{Float64, Vector{Int}}`: (modularity Q, community labels per vertex)

# Example
```julia
modularity, labels = Mycelia.Rhizomorph.compute_modularity(graph)
```
"""
function compute_modularity(g)
    graph = g isa MetaGraphsNext.MetaGraph ? g.graph : g
    n = Graphs.nv(graph)
    m = Graphs.ne(graph)

    if n == 0 || m == 0
        return 0.0, Int[]
    end

    components = Graphs.connected_components(graph)
    labels = zeros(Int, n)
    for (i, comp) in enumerate(components)
        for v in comp
            labels[v] = i
        end
    end

    degrees = Graphs.degree(graph)
    Q = 0.0

    for e in Graphs.edges(graph)
        u, v = Graphs.src(e), Graphs.dst(e)
        if labels[u] == labels[v]
            Q += 1.0 - (degrees[u] * degrees[v]) / (2.0 * m)
        end
    end

    Q /= m
    return Q, labels
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute closeness centrality: C(v) = (n-1) / Σ d(v, u).
For large graphs, samples `sample_size` random vertices.

Works on both `Graphs.AbstractGraph` and `MetaGraphsNext.MetaGraph` objects.

# Arguments
- `g`: Graph (Graphs.jl `AbstractGraph` or `MetaGraphsNext.MetaGraph`).
- `sample_size::Int=1000`: Maximum number of vertices to sample when the graph is large.

# Returns
- `Vector{Float64}`: Closeness centrality per sampled vertex

# Example
```julia
closeness = Mycelia.Rhizomorph.compute_closeness_centrality(graph; sample_size = 100)
```
"""
function compute_closeness_centrality(g; sample_size::Int = 1000)
    graph = g isa MetaGraphsNext.MetaGraph ? g.graph : g
    n = Graphs.nv(graph)
    if n == 0
        return Float64[]
    end

    nodes = n > sample_size ? Mycelia.Random.rand(1:n, sample_size) : collect(1:n)
    closeness = zeros(length(nodes))

    for (i, v) in enumerate(nodes)
        dists = Graphs.gdistances(graph, v)
        finite_dists = filter(d -> d < typemax(Int) && d > 0, dists)
        if !isempty(finite_dists)
            closeness[i] = length(finite_dists) / sum(finite_dists)
        end
    end

    return closeness
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute Betti numbers for a graph:
- β₀ = number of connected components
- β₁ = number of independent cycles = |E| - |V| + β₀ (Euler characteristic)

Works on both `Graphs.AbstractGraph` and `MetaGraphsNext.MetaGraph` objects.

# Arguments
- `g`: Graph (Graphs.jl `AbstractGraph` or `MetaGraphsNext.MetaGraph`).

# Returns
- `Tuple{Int, Int}`: (β₀, β₁)

# Example
```julia
beta_0, beta_1 = Mycelia.Rhizomorph.compute_betti_numbers(graph)
```
"""
function compute_betti_numbers(g)
    graph = g isa MetaGraphsNext.MetaGraph ? g.graph : g
    nv_count = Graphs.nv(graph)
    ne_count = Graphs.ne(graph)
    components = Graphs.connected_components(graph)
    β0 = length(components)
    β1 = ne_count - nv_count + β0
    return β0, max(0, β1)
end
