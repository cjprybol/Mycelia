import Graphs
import MetaGraphsNext

function is_topologically_equal(graph_a::MetaGraphsNext.MetaGraph, graph_b::MetaGraphsNext.MetaGraph)
    Graphs.nv(graph_a.graph) == Graphs.nv(graph_b.graph) || return false
    Graphs.ne(graph_a.graph) == Graphs.ne(graph_b.graph) || return false

    vertices_a = collect(MetaGraphsNext.labels(graph_a))
    vertices_b = collect(MetaGraphsNext.labels(graph_b))
    isempty(vertices_a) && return true

    adjacency_a = _adjacency_matrix(graph_a, vertices_a)
    adjacency_b = _adjacency_matrix(graph_b, vertices_b)

    initial_colors_a = _compress_colors(_initial_vertex_colors(adjacency_a))
    initial_colors_b = _compress_colors(_initial_vertex_colors(adjacency_b))
    sort(initial_colors_a) == sort(initial_colors_b) || return false

    refined_colors_a = _refine_vertex_colors(adjacency_a, initial_colors_a)
    refined_colors_b = _refine_vertex_colors(adjacency_b, initial_colors_b)
    sort(refined_colors_a) == sort(refined_colors_b) || return false

    classes_a = _group_vertices_by_color(refined_colors_a)
    classes_b = _group_vertices_by_color(refined_colors_b)
    sort(length.(values(classes_a))) == sort(length.(values(classes_b))) || return false

    ordered_classes = sort!(collect(keys(classes_a)); by = color -> length(classes_a[color]))
    mapping = fill(0, length(vertices_a))
    used = fill(false, length(vertices_b))

    return _search_isomorphism!(
        adjacency_a,
        adjacency_b,
        classes_a,
        classes_b,
        ordered_classes,
        mapping,
        used,
        1
    )
end

function is_semantically_equal(graph_a::MetaGraphsNext.MetaGraph, graph_b::MetaGraphsNext.MetaGraph)
    return _assembled_sequences(graph_a) == _assembled_sequences(graph_b)
end

function _adjacency_matrix(graph::MetaGraphsNext.MetaGraph, vertices)
    n_vertices = length(vertices)
    adjacency = falses(n_vertices, n_vertices)
    vertex_to_index = Dict(vertex => index for (index, vertex) in enumerate(vertices))

    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        adjacency[vertex_to_index[src], vertex_to_index[dst]] = true
    end

    return adjacency
end

function _initial_vertex_colors(adjacency::BitMatrix)
    n_vertices = size(adjacency, 1)
    colors = Vector{Tuple{Int, Int, Bool}}(undef, n_vertices)

    for index in 1:n_vertices
        indegree = count(@view adjacency[:, index])
        outdegree = count(@view adjacency[index, :])
        colors[index] = (indegree, outdegree, adjacency[index, index])
    end

    return colors
end

function _refine_vertex_colors(adjacency::BitMatrix, colors)
    refined = colors

    while true
        signatures = map(eachindex(refined)) do index
            incoming = Tuple(sort!([refined[neighbor] for neighbor in findall(@view adjacency[:, index])]))
            outgoing = Tuple(sort!([refined[neighbor] for neighbor in findall(@view adjacency[index, :])]))
            return (refined[index], incoming, outgoing)
        end
        next_colors = _compress_colors(signatures)

        next_colors == refined && return refined
        refined = next_colors
    end
end

function _compress_colors(signatures)
    signature_to_color = Dict{Any, Int}()
    colors = Vector{Int}(undef, length(signatures))

    for (index, signature) in enumerate(signatures)
        colors[index] = get!(signature_to_color, signature) do
            length(signature_to_color) + 1
        end
    end

    return colors
end

function _group_vertices_by_color(colors)
    groups = Dict{Any, Vector{Int}}()

    for (index, color) in enumerate(colors)
        push!(get!(groups, color, Int[]), index)
    end

    return groups
end

function _search_isomorphism!(
        adjacency_a::BitMatrix,
        adjacency_b::BitMatrix,
        classes_a::Dict,
        classes_b::Dict,
        ordered_classes,
        mapping::Vector{Int},
        used::Vector{Bool},
        class_index::Int
)
    class_index > length(ordered_classes) && return true

    color = ordered_classes[class_index]
    vertices_a = classes_a[color]
    vertices_b = classes_b[color]
    length(vertices_a) == length(vertices_b) || return false

    available_targets = [vertex for vertex in vertices_b if !used[vertex]]
    return _search_class_assignments!(
        adjacency_a,
        adjacency_b,
        classes_a,
        classes_b,
        ordered_classes,
        mapping,
        used,
        class_index,
        vertices_a,
        available_targets,
        1
    )
end

function _search_class_assignments!(
        adjacency_a::BitMatrix,
        adjacency_b::BitMatrix,
        classes_a::Dict,
        classes_b::Dict,
        ordered_classes,
        mapping::Vector{Int},
        used::Vector{Bool},
        class_index::Int,
        vertices_a::Vector{Int},
        available_targets::Vector{Int},
        vertex_index::Int
)
    if vertex_index > length(vertices_a)
        return _search_isomorphism!(
            adjacency_a,
            adjacency_b,
            classes_a,
            classes_b,
            ordered_classes,
            mapping,
            used,
            class_index + 1
        )
    end

    vertex_a = vertices_a[vertex_index]

    for target_position in eachindex(available_targets)
        vertex_b = available_targets[target_position]
        if !_mapping_is_consistent(adjacency_a, adjacency_b, mapping, vertex_a, vertex_b)
            continue
        end

        mapping[vertex_a] = vertex_b
        used[vertex_b] = true
        deleteat!(available_targets, target_position)

        if _search_class_assignments!(
            adjacency_a,
            adjacency_b,
            classes_a,
            classes_b,
            ordered_classes,
            mapping,
            used,
            class_index,
            vertices_a,
            available_targets,
            vertex_index + 1
        )
            return true
        end

        insert!(available_targets, target_position, vertex_b)
        used[vertex_b] = false
        mapping[vertex_a] = 0
    end

    return false
end

function _mapping_is_consistent(
        adjacency_a::BitMatrix,
        adjacency_b::BitMatrix,
        mapping::Vector{Int},
        candidate_a::Int,
        candidate_b::Int
)
    for other_a in eachindex(mapping)
        other_b = mapping[other_a]
        other_b == 0 && continue

        adjacency_a[candidate_a, other_a] == adjacency_b[candidate_b, other_b] || return false
        adjacency_a[other_a, candidate_a] == adjacency_b[other_b, candidate_b] || return false
    end

    return adjacency_a[candidate_a, candidate_a] == adjacency_b[candidate_b, candidate_b]
end

function _assembled_sequences(graph::MetaGraphsNext.MetaGraph)
    Graphs.nv(graph.graph) == 0 && return String[]

    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
    isempty(paths) && return sort!(string.(MetaGraphsNext.labels(graph)))

    sequences = String[]
    for path in paths
        push!(sequences, string(Mycelia.Rhizomorph.path_to_sequence(path, graph)))
    end

    return sort!(sequences)
end
