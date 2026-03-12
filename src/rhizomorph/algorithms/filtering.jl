# Graph Filtering Functions
#
# Functions for extracting subgraphs based on connected components,
# coverage, and other criteria. Useful for producing readable
# visualizations of large graphs.

"""
    filter_largest_components(graph::MetaGraphsNext.MetaGraph;
        max_components::Union{Int, Nothing}=nothing,
        min_fraction::Float64=0.0,
        min_nodes::Int=1)

Extract a subgraph containing only the largest connected components.

# Arguments
- `graph`: The input MetaGraphsNext graph
- `max_components`: Maximum number of components to keep (by descending size).
  If `nothing`, uses `min_fraction` or `min_nodes` to filter.
- `min_fraction`: Keep components with at least this fraction of total nodes (0.0-1.0).
  Ignored if `max_components` is set.
- `min_nodes`: Keep components with at least this many nodes. Ignored if
  `max_components` is set.

# Returns
- `MetaGraphsNext.MetaGraph`: A new graph containing only the selected components

# Examples
```julia
# Keep only the largest connected component
subgraph = filter_largest_components(graph; max_components=1)

# Keep components representing the top 10% of nodes
subgraph = filter_largest_components(graph; min_fraction=0.10)

# Keep components with at least 5 nodes
subgraph = filter_largest_components(graph; min_nodes=5)
```
"""
function filter_largest_components(graph::MetaGraphsNext.MetaGraph;
        max_components::Union{Int, Nothing} = nothing,
        min_fraction::Float64 = 0.0,
        min_nodes::Int = 1)
    components = Graphs.connected_components(graph.graph)

    if isempty(components)
        return _empty_graph_like(graph)
    end

    # Sort components by size (descending)
    sort!(components; by = length, rev = true)

    total_nodes = Graphs.nv(graph)

    # Determine which components to keep
    keep_components = if max_components !== nothing
        components[1:min(max_components, length(components))]
    else
        threshold = max(min_nodes, ceil(Int, min_fraction * total_nodes))
        filter(c -> length(c) >= threshold, components)
    end

    if isempty(keep_components)
        keep_components = [components[1]]
    end

    # Collect vertex indices to keep
    keep_vertices = Set{Int}()
    for comp in keep_components
        union!(keep_vertices, comp)
    end

    return _induced_subgraph(graph, keep_vertices)
end

"""
    filter_top_coverage(graph::MetaGraphsNext.MetaGraph; top_fraction::Float64=0.1)

Extract a subgraph containing only nodes in the top fraction by depth/coverage.

# Arguments
- `graph`: The input MetaGraphsNext graph
- `top_fraction`: Keep nodes with coverage in the top this fraction (0.0-1.0).
  Default 0.1 (top 10%).

# Returns
- `MetaGraphsNext.MetaGraph`: A new graph with only high-coverage nodes and their edges
"""
function filter_top_coverage(graph::MetaGraphsNext.MetaGraph;
        top_fraction::Float64 = 0.1)
    labels_vec = collect(MetaGraphsNext.labels(graph))
    if isempty(labels_vec)
        return _empty_graph_like(graph)
    end

    # Extract coverage/depth for each vertex
    coverages = Float64[]
    for label in labels_vec
        vertex_data = graph[label]
        depth = if hasfield(typeof(vertex_data), :evidence)
            Float64(count_evidence(vertex_data))
        elseif hasfield(typeof(vertex_data), :coverage)
            Float64(length(vertex_data.coverage))
        else
            1.0
        end
        push!(coverages, depth)
    end

    # Find coverage threshold
    sorted_cov = sort(coverages; rev = true)
    cutoff_idx = max(1, ceil(Int, top_fraction * length(sorted_cov)))
    threshold = sorted_cov[min(cutoff_idx, length(sorted_cov))]

    # Collect vertex indices to keep
    keep_vertices = Set{Int}()
    for (i, cov) in enumerate(coverages)
        if cov >= threshold
            push!(keep_vertices, MetaGraphsNext.code_for(graph, labels_vec[i]))
        end
    end

    return _induced_subgraph(graph, keep_vertices)
end

# Internal: create an empty graph with the same type parameters
function _empty_graph_like(
        graph::MetaGraphsNext.MetaGraph{
        Code, G, Label, VData, EData}) where {
        Code, G, Label, VData, EData}
    return MetaGraphsNext.MetaGraph(
        G();
        label_type = Label,
        vertex_data_type = VData,
        edge_data_type = EData,
        weight_function = graph.weight_function
    )
end

# Internal: build induced subgraph from a set of vertex codes
function _induced_subgraph(graph::MetaGraphsNext.MetaGraph, keep_vertices::Set{Int})
    new_graph = _empty_graph_like(graph)

    labels_vec = collect(MetaGraphsNext.labels(graph))

    # Add kept vertices
    for vi in keep_vertices
        label = labels_vec[vi]
        new_graph[label] = graph[label]
    end

    # Add edges between kept vertices
    for edge_labels in MetaGraphsNext.edge_labels(graph)
        if length(edge_labels) == 2
            src_label, dst_label = edge_labels
            src_code = MetaGraphsNext.code_for(graph, src_label)
            dst_code = MetaGraphsNext.code_for(graph, dst_label)
            if src_code in keep_vertices && dst_code in keep_vertices
                new_graph[src_label, dst_label] = graph[src_label, dst_label]
            end
        end
    end

    return new_graph
end
