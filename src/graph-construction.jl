"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function construct(KMER_TYPE, fastx, out)
    mkpath(dirname(out))
    if !occursin(r"\.jld2$", out)
        out *= ".jld2"
    end
    if !isfile(out)
        graph = fastx_to_kmer_graph(KMER_TYPE, fastx)
        @info "saving graph"
        FileIO.save(out, Dict("graph" => graph))
        return graph
    else
        @info "graph $out already exists, loading existing"
        return load_graph(out)
    end
end

function construct(args)
    @show args
    @assert (0 < args["k"] < 64) && isodd(args["k"]) 
    KMER_TYPE = BioSequences.BigDNAMer{args["k"]}
    construct(KMER_TYPE, args["fastx"], args["out"])
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function load_graph(file)
    return FileIO.load(file)["graph"]
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_edge_to_graph(graph, edge_mer, kmers)
#     edge = BioSequences.LongDNASeq(edge_mer.fw)
#     k = length(first(kmers))
# #     canonical_src = BioSequences.DNAMer{k}(BioSequences.canonical!(edge[1:end-1]))
# #     canonical_dst = BioSequences.DNAMer{k}(BioSequences.canonical!(edge[2:end]))

#     canonical_src = BioSequences.canonical(BioSequences.DNAMer{k}(edge[1:end-1]))
    
#     src_index_range = searchsorted(kmers, canonical_src)
#     if isempty(src_index_range)
#         return
#     else
#         @assert length(src_index_range) == 1
#     end
#     src_index = first(src_index_range)

#     canonical_dst = BioSequences.canonical(BioSequences.DNAMer{k}(edge[2:end]))
#     dst_index_range = searchsorted(kmers, canonical_dst)
#     if isempty(dst_index_range)
#         return
#     else
#         @assert length(dst_index_range) == 1
#     end
#     dst_index = first(dst_index_range)
#     graph_edge = Graphs.Edge(src_index, dst_index)
#     Graphs.add_edge!(graph, graph_edge)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create an in-memory kmer-graph that records:
- all kmers
- counts
- all *observed* edges between kmers
- edge orientations
- edge counts

```jldoctest
julia> 1 + 1
2
```
"""

# note: indexing will break if we mix BigDNAMer and normal DNAMer types, so just force requirement of BigDNAMer
# can't do the restriction because the type here is `DataType` rather than the actual kmer type
function fastx_to_kmer_graph(KMER_TYPE, fastxs::AbstractVector{<:AbstractString})
    
    if !(KMER_TYPE <: BioSequences.BigDNAMer)
        error()
    end
    
    @info "counting kmers"
    @time kmer_counts = count_canonical_kmers(KMER_TYPE, fastxs)
    
    @info "initializing graph"
    K = length(keys(kmer_counts))
    k = length(first(keys(kmer_counts)))
    # create an undirected kmer graph from the sequence
    graph = MetaGraphs.MetaGraph(K)
    # graph = Graphs.SimpleGraph(K)

#     MetaGraphs.set_prop!(graph, :kmer_counts, kmer_counts)
    MetaGraphs.set_prop!(graph, :k, k)

    @info "adding node metadata"
    ProgressMeter.@showprogress for (i, (kmer, count)) in enumerate(kmer_counts)
    #     @show i, kmer, count
        MetaGraphs.set_prop!(graph, i, :kmer, kmer)
        MetaGraphs.set_prop!(graph, i, :count, count)
    end
    # allow graph[kmer, :kmer] to dict-lookup the index of a kmer
    MetaGraphs.set_indexing_prop!(graph, :kmer)

#     kmers = collect(keys(kmer_counts))

    # p = ProgressMeter.Progress(8452, 1)
    # 50 minutes
    # 40 minutes
    # 0:09:51!
    @info "adding edges"
    ProgressMeter.@showprogress for (i, fastx) in enumerate(fastxs)
#         @show i, fastx
#         n_records = count_records(fastx)
#         p = ProgressMeter.Progress(n_records, 1)
        for record in open_fastx(fastx)
            for edge_mer in BioSequences.each(BioSequences.BigDNAMer{k+1}, FASTX.sequence(record))
#                 add_edge_to_graph(graph, edge_mer, kmers)
                add_edge_to_simple_kmer_graph!(graph, edge_mer)
            end
#             ProgressMeter.next!(p)
        end
    end
    return graph
end

function fastx_to_kmer_graph(KMER_TYPE, fastx::AbstractString)
    fastx_to_kmer_graph(KMER_TYPE, [fastx])
end


function edgemer_to_vertex_kmers(edgemer)
    a = BioSequences.BigDNAMer(edgemer[i] for i in 1:length(edgemer)-1)
    b = BioSequences.BigDNAMer(edgemer[i] for i in 2:length(edgemer))
    return a, b
end

# @inline function add_edge_to_simple_kmer_graph!(simple_kmer_graph, kmers, sequence_edge)
@inline function add_edge_to_simple_kmer_graph!(simple_kmer_graph, sequence_edge)
    observed_source_kmer, observed_destination_kmer = edgemer_to_vertex_kmers(sequence_edge.fw)
    canonical_source_kmer = BioSequences.canonical(observed_source_kmer)
    canonical_destination_kmer = BioSequences.canonical(observed_destination_kmer)
    source_kmer_index = simple_kmer_graph[canonical_source_kmer, :kmer]
    desination_kmer_index = simple_kmer_graph[canonical_destination_kmer, :kmer]
    
    if source_kmer_index > desination_kmer_index
        observed_source_kmer, observed_destination_kmer = edgemer_to_vertex_kmers(sequence_edge.bw)
        canonical_source_kmer = BioSequences.canonical(observed_source_kmer)
        canonical_destination_kmer = BioSequences.canonical(observed_destination_kmer)
        source_kmer_index = simple_kmer_graph[canonical_source_kmer, :kmer]
        desination_kmer_index = simple_kmer_graph[canonical_destination_kmer, :kmer]
    end
    
    @assert source_kmer_index <= desination_kmer_index

    oriented_source_kmer = 
        (canonical_kmer = canonical_source_kmer,
         orientation = BioSequences.iscanonical(observed_source_kmer))

    oriented_destination_kmer = 
        (canonical_kmer = canonical_destination_kmer,
         orientation = BioSequences.iscanonical(observed_destination_kmer))

    oriented_source_vertex = 
#         (vertex = searchsortedfirst(kmers, oriented_source_kmer.canonical_kmer),
        (vertex = simple_kmer_graph[oriented_source_kmer.canonical_kmer, :kmer],
         orientation = oriented_source_kmer.orientation)

    oriented_destination_vertex = 
#         (vertex = searchsortedfirst(kmers, oriented_destination_kmer.canonical_kmer),
        (vertex = simple_kmer_graph[oriented_destination_kmer.canonical_kmer, :kmer],
         orientation = oriented_destination_kmer.orientation)

#     forward_edge = Graphs.Edge(oriented_source_vertex.vertex, oriented_destination_vertex.vertex)
#     forward_edge_orientations = 
#         (source_orientation = oriented_source_vertex.orientation,
#          destination_orientation = oriented_destination_vertex.orientation)
    
#     reverse_edge = Graphs.Edge(oriented_destination_vertex.vertex, oriented_source_vertex.vertex)
#     reverse_edge_orientations = 
#         (source_orientation = !oriented_destination_vertex.orientation,
#          destination_orientation = !oriented_source_vertex.orientation)
    
    edge = Graphs.Edge(oriented_source_vertex.vertex, oriented_destination_vertex.vertex)
    edge_orientations = 
        (source_orientation = oriented_source_vertex.orientation,
         destination_orientation = oriented_destination_vertex.orientation)
    
#     orientations = Set([forward_edge_orientations, reverse_edge_orientations])
    orientations = Set([edge_orientations])
    if Graphs.has_edge(simple_kmer_graph, edge)
        edge_weight = MetaGraphs.get_prop(simple_kmer_graph, edge, :weight) + 1
        orientations = union(MetaGraphs.get_prop(simple_kmer_graph, edge, :orientations), orientations)
    else
        Graphs.add_edge!(simple_kmer_graph, edge)
        edge_weight = 1
#         @show Graphs.ne(simple_kmer_graph)
#         @show forward_edge
    end
    MetaGraphs.set_prop!(simple_kmer_graph, edge, :weight, edge_weight)
    MetaGraphs.set_prop!(simple_kmer_graph, edge, :orientations, orientations)
end

# function fastx_to_simple_kmer_graph(KMER_TYPE, fastx::AbstractString; minimum_coverage::Int=1)
#     fastx_to_simple_kmer_graph(KMER_TYPE, [fastx], minimum_coverage=minimum_coverage)
# end

# function fastx_to_simple_kmer_graph(KMER_TYPE, fastxs::AbstractVector{<:AbstractString}; minimum_coverage::Int=1)
#     @info "counting kmers"
#     canonical_kmer_counts = count_canonical_kmers(KMER_TYPE, fastxs)
#     # hard filter any nodes that are less frequent than minimum coverage threshold
#     canonical_kmer_counts = filter(canonical_kmer_count -> last(canonical_kmer_count) >= minimum_coverage, canonical_kmer_counts)
#     simple_kmer_graph = MetaGraphs.MetaDiGraph(length(canonical_kmer_counts))
    
#     k = length(first(keys(canonical_kmer_counts)))

#     MetaGraphs.set_prop!(simple_kmer_graph, :k, k)

#     @info "setting metadata on vertices"
#     ProgressMeter.@showprogress for (vertex, (kmer, count)) in enumerate(canonical_kmer_counts)
#         MetaGraphs.set_prop!(simple_kmer_graph, vertex, :kmer, kmer)
#         MetaGraphs.set_prop!(simple_kmer_graph, vertex, :weight, count)
#     end

#     kmers = collect(keys(canonical_kmer_counts))

#     EDGE_MER = BioSequences.BigDNAMer{k+1}
#     @info "loading fastx files into graph"
#     ProgressMeter.@showprogress for fastx in fastxs
#         n_records = 0
#         for record in (open_fastx(fastx))
#             n_records += 1
#         end
#         p = ProgressMeter.Progress(n_records, 1)   # minimum update interval: 1 second
#         for record in (open_fastx(fastx))
#             sequence = FASTX.sequence(record)
#             edge_iterator = BioSequences.each(EDGE_MER, sequence)
#             for sequence_edge in edge_iterator
#                 add_edge_to_simple_kmer_graph!(simple_kmer_graph, kmers, sequence_edge)
#             end
#             ProgressMeter.next!(p)
#         end
#     end
#     return simple_kmer_graph
# end



# @inline function add_edge_to_kmer_graph!(kmer_graph, kmers, sequence_edge, record_identifier)
#     observed_source_kmer, observed_destination_kmer = edgemer_to_vertex_kmers(sequence_edge.fw)

#     oriented_source_kmer = 
#         (canonical_kmer = BioSequences.canonical(observed_source_kmer),
#          orientation = BioSequences.iscanonical(observed_source_kmer))

#     oriented_destination_kmer = 
#         (canonical_kmer = BioSequences.canonical(observed_destination_kmer),
#          orientation = BioSequences.iscanonical(observed_destination_kmer))

#     oriented_source_vertex = 
#         (vertex = searchsortedfirst(kmers, oriented_source_kmer.canonical_kmer),
#          orientation = oriented_source_kmer.orientation)

#     oriented_destination_vertex = 
#         (vertex = searchsortedfirst(kmers, oriented_destination_kmer.canonical_kmer),
#          orientation = oriented_destination_kmer.orientation)

#     source_evidence = 
#         (record = record_identifier,
#          index = sequence_edge.position,
#          orientation = oriented_source_vertex.orientation)

#     destination_evidence = 
#         (record = record_identifier,
#          index = sequence_edge.position + 1,
#          orientation = oriented_destination_vertex.orientation)

#     set_metadata!(kmer_graph, oriented_source_vertex.vertex, :evidence, source_evidence)
#     new_weight = length(kmer_graph.vprops[oriented_source_vertex.vertex][:evidence])
#     MetaGraphs.set_prop!(kmer_graph, oriented_source_vertex.vertex, :weight, new_weight)

#     set_metadata!(kmer_graph, oriented_destination_vertex.vertex, :evidence, destination_evidence)
#     new_weight = length(kmer_graph.vprops[oriented_destination_vertex.vertex][:evidence])
#     MetaGraphs.set_prop!(kmer_graph, oriented_destination_vertex.vertex, :weight, new_weight)
    

#     forward_edge = Graphs.Edge(oriented_source_vertex.vertex, oriented_destination_vertex.vertex)

#     Graphs.add_edge!(kmer_graph, forward_edge)

#     forward_edge_orientations = 
#         (source_orientation = oriented_source_vertex.orientation,
#          destination_orientation = oriented_destination_vertex.orientation)

#     set_metadata!(kmer_graph, forward_edge, :orientations, forward_edge_orientations)

#     forward_edge_evidence = (
#         record = record_identifier,
#         index = sequence_edge.position,
#         orientation = true
#     )

#     set_metadata!(kmer_graph, forward_edge, :evidence, forward_edge_evidence)
#     new_weight = length(kmer_graph.eprops[forward_edge][:evidence])
#     MetaGraphs.set_prop!(kmer_graph, forward_edge, :weight, new_weight)

#     reverse_edge = Graphs.Edge(oriented_destination_vertex.vertex, oriented_source_vertex.vertex)

#     Graphs.add_edge!(kmer_graph, reverse_edge)

#     reverse_edge_orientations = 
#         (source_orientation = !oriented_destination_vertex.orientation,
#          destination_orientation = !oriented_source_vertex.orientation)

#     set_metadata!(kmer_graph, reverse_edge, :orientations, reverse_edge_orientations)

#     reverse_edge_evidence = (
#         record = record_identifier,
#         index = sequence_edge.position,
#         orientation = false
#     )

#     set_metadata!(kmer_graph, reverse_edge, :evidence, reverse_edge_evidence)
#     new_weight = length(kmer_graph.eprops[reverse_edge][:evidence])
#     MetaGraphs.set_prop!(kmer_graph, reverse_edge, :weight, new_weight)
# end

# # function fastx_to_kmer_graph(::Type{KMER_TYPE}, fastxs) where {KMER_TYPE <: BioSequences.AbstractMer{A, K}} where {A, K}
# function fastx_to_kmer_graph(KMER_TYPE, fastx::AbstractString)
#     fastx_to_kmer_graph(KMER_TYPE, [fastx])
# end

# function fastx_to_kmer_graph(KMER_TYPE, fastxs::AbstractVector{<:AbstractString})
#     @info "assessing kmers"
#     kmer_counts = count_canonical_kmers(KMER_TYPE, fastxs)
# #     kmer_set = Set{KMER_TYPE}()
# #     for fastxs in fastxs
# #         kmer_set = union!(kmer_set, collect(keys(count_canonical_kmers(KMER_TYPE, fastxs))))
# #     end
# #     kmers = unique(sort(collect(kmer_set)))
    
#     kmer_graph = MetaGraphs.MetaDiGraph(length(kmer_counts))
#     k = length(first(keys(kmer_counts)))
#     kmers = collect(keys(kmer_counts))
#     MetaGraphs.set_prop!(kmer_graph, :k, k)
#     # don't set this since when we filter an induced subgraph, these don't update
# #     MetaGraphs.set_prop!(kmer_graph, :kmers, kmers)
#     for (vertex, (kmer, count)) in enumerate(kmer_counts)
#         MetaGraphs.set_prop!(kmer_graph, vertex, :kmer, kmer)
#         MetaGraphs.set_prop!(kmer_graph, vertex, :weight, count)
#     end
#     EDGE_MER = BioSequences.BigDNAMer{k+1}
#     @info "creating graph"
#     ProgressMeter.@showprogress for fastx in fastxs
#         fastx_io = open_fastx(fastx)
#         for record in fastx_io
#             sequence = FASTX.sequence(record)
#             record_identifier = FASTX.identifier(record) 
#             edge_iterator = BioSequences.each(EDGE_MER, sequence)
#             for sequence_edge in edge_iterator
#                 add_edge_to_kmer_graph!(kmer_graph, kmers, sequence_edge, record_identifier)
#             end
#         end
#         close(fastx_io)
#     end
#     return kmer_graph
# end