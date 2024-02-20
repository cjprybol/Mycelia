# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function graph_to_kmers(g)
#     kmers = [g.vprops[v][:kmer] for v in Graphs.vertices(g)]
#     return kmers
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function graph_to_edge_sequences(g)
#     edges = Set{BioSequences.BigDNAMer{g.gprops[:k]+1}}()
#     for edge in Graphs.edges(g)
#         src_kmer = g.vprops[edge.src][:kmer]
#         dst_kmer = g.vprops[edge.dst][:kmer]
#         for orientation in g.eprops[edge][:orientations]
#             if orientation.source_orientation
#                 oriented_src_kmer = src_kmer
#             else
#                 oriented_src_kmer = BioSequences.reverse_complement(src_kmer)
#             end
#             if orientation.destination_orientation
#                 oriented_dst_kmer = dst_kmer
#             else
#                 oriented_dst_kmer = BioSequences.reverse_complement(dst_kmer)
#             end
#             for i in 1:g.gprops[:k]-1
#                 @assert oriented_src_kmer[i+1] == oriented_dst_kmer[i]
#             end
#             edge_mer = BioSequences.BigDNAMer((nuc for nuc in oriented_src_kmer)..., last(oriented_dst_kmer))
#             push!(edges, BioSequences.canonical(edge_mer))
#         end
#     end
#     return edges
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function kmer_graph_distances(g1, g2)
#     g1_kmers = Set(graph_to_kmers(g1))
#     g1_edges = graph_to_edge_sequences(g1)
    
#     g2_kmers = Set(graph_to_kmers(g2))
#     g2_edges = graph_to_edge_sequences(g2)
    
#     kmer_distance = 1 - LSHFunctions.jaccard(g1_kmers, g2_kmers)
#     edge_distance = 1 - LSHFunctions.jaccard(g1_edges, g2_edges)
    
#     result = (
#         kmer_distance = kmer_distance,
#         edge_distance = edge_distance
#     )
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function set_metadata!(kmer_graph, vertex::V, key, value) where V <: Integer
#     if MetaGraphs.has_prop(kmer_graph, vertex, key)
#         push!(kmer_graph.vprops[vertex][key], value)
#     else
#         MetaGraphs.set_prop!(kmer_graph, vertex, key, Set([value]))
#     end
#     return true
# end

# function set_metadata!(kmer_graph, edge::E, key, value) where E <: Graphs.Edge
#     if MetaGraphs.has_prop(kmer_graph, edge, key)
#         current_value = MetaGraphs.get_prop(kmer_graph, edge, key)
#         updated_value = push!(current_value, value)
#         MetaGraphs.set_prop!(kmer_graph, edge, key, updated_value)
#     else
#         MetaGraphs.set_prop!(kmer_graph, edge, key, Set([value]))
#     end
#     return true
# end
