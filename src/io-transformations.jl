function graph_to_gfa(graph, kmer_size, outfile="$(kmer_size).gfa")
    kmer_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size})))
    # add fastq here too???
    record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
    edgemer_edges = collect(MetaGraphs.filter_edges(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size+1})))
    open(outfile, "w") do io
        println(io, "H\tVN:Z:1.0")
        for vertex in kmer_vertices
            if haskey(graph.vprops[vertex], :kmer)
                sequence = graph.vprops[vertex][:kmer]
            else
                sequence = graph.vprops[vertex][:sequence]
            end
            # depth = graph.vprops[vertex][:count]
            # depth = length(graph.vprops[vertex][:evidence])
            total_count = 0
            vertex_outneighbors = Graphs.outneighbors(graph, vertex)
            for connected_record in intersect(vertex_outneighbors, record_vertices)
                edge = Graphs.Edge(vertex, connected_record)
                total_count += MetaGraphs.get_prop(graph, edge, :count)
            end

            fields = ["S", "$vertex", sequence, "DP:f:$(total_count)"]
            line = join(fields, '\t')
            println(io, line)
        end
        overlap = kmer_size - 1
        for edgemer_edge in edgemer_edges
            source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
            link = ["L",
                        edgemer_edge.src,
                        BioSequences.iscanonical(source_kmer) ? '+' : '-',
                        edgemer_edge.dst,
                        BioSequences.iscanonical(dest_kmer) ? '+' : '-',
                        "$(overlap)M"]
            line = join(link, '\t')
            println(io, line)
        end
    end
    return
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function open_fastx(path::AbstractString)
    if isfile(path)
        io = open(path)
    elseif occursin(r"^ftp", path) || occursin(r"^http", path)
        path = replace(path, r"^ftp:" => "http:")
        io = IOBuffer(HTTP.get(path).body)
    else
        error("unable to locate file $path")
    end
    path_base = basename(path)
    if occursin(r"\.gz$", path_base)
        io = CodecZlib.GzipDecompressorStream(io)
        path_base = replace(path_base, ".gz" => "")
    end
    if occursin(r"\.(fasta|fna|faa|fa)$", path_base)
        fastx_io = FASTX.FASTA.Reader(io)
    elseif occursin(r"\.(fastq|fq)$", path_base)
        fastx_io = FASTX.FASTQ.Reader(io)
    else
        @show path_base
        error()
    end
    return fastx_io
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function open_gff(path::String)
    if isfile(path)
        io = open(path)
    elseif occursin(r"^ftp", path) || occursin(r"^http", path)
        path = replace(path, r"^ftp:" => "http:")
        io = IOBuffer(HTTP.get(path).body)
    else
        error("unable to locate file $path")
    end
    path_base = basename(path)
    if occursin(r"\.gz$", path_base)
        io = CodecZlib.GzipDecompressorStream(io)
    end
    gff_io = GFF3.Reader(io)
    return gff_io
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function plot_kmer_frequency_spectra(counts; log_scale = log2, kwargs...)
    kmer_counts_hist = StatsBase.countmap(c for c in counts)
    xs = collect(keys(kmer_counts_hist))
    ys = collect(values(kmer_counts_hist))
    if isa(log_scale, Function)
        xs = log_scale.(xs)
        ys = log_scale.(ys)
    end
    
    p = StatsPlots.plot(
        xs,
        ys,
        xlims = (0, maximum(xs) + max(1, ceil(0.1 * maximum(xs)))),
        ylims = (0, maximum(ys) + max(1, ceil(0.1 * maximum(ys)))),
        seriestype = :scatter,
        legend = false,
        xlabel = isa(log_scale, Function) ? "$(log_scale)(observed frequency)" : "observed frequency",
        ylabel = isa(log_scale, Function) ? "$(log_scale)(# of kmers)" : "observed frequency",
        ;kwargs...
    )
    return p
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
function convert(args)
    @show args
    in_type = missing
    out_type = missing
    if occursin(r"\.jld2$", args["in"])
        in_type = :jld2
    elseif occursin(r"\.gfa$", args["in"])
        in_type = :gfa
    elseif occursin(r"\.neo4j$", args["in"])
        in_type = :neo4j
    end

    if occursin(r"\.jld2$", args["out"])
        out_type = :jld2
    elseif occursin(r"\.gfa$", args["out"])
        out_type = :gfa
    elseif occursin(r"\.neo4j$", args["out"])
        out_type = :neo4j
    end
    
    if ismissing(in_type) || ismissing(out_type)
        error("unable to determine in and out types")
    end
    
    if (in_type == :jld2) && (out_type == :jld2)
        # done
    elseif (in_type == :jld2) && (out_type != :jld2)
        # convert out of jld2
        if out_type == :gfa
            loaded = FileIO.load(args["in"])
            @assert haskey(loaded, "graph")
            graph = loaded["graph"]
            graph_to_gfa(graph, args["out"])
        end
    elseif (in_type != :jld2) && (out_type == :jld2)
        # convert into jld2
    else
        # need to convert from input to jld2 as an intermediate first
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function plot_graph(graph)
    
#     kmer_counts = MetaGraphs.get_prop(graph, :kmer_counts)
    kmers = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph)]
    counts = [MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)]
    scale = 150
    
    n = Graphs.nv(graph)
    p = GraphRecipes.graphplot(
        graph,
#         markersize = 1/log2(n),
        markersize = 1/2^2,
        size = (2 * scale * log(n), scale * log(n)),
        node_weights = counts,
        names = kmers)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function save_graph(graph::Graphs.AbstractGraph, outfile::String)
    if !occursin(r"\.jld2$", outfile)
        outfile *= ".jld2"
    end
    FileIO.save(outfile, Dict("graph" => graph))
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function load_graph(file::String)
    return FileIO.load(file)["graph"]
end

# OLD FOR SIMPLE KMER GRAPHS
# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # Description

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# function graph_to_gfa(graph, outfile)
#     open(outfile, "w") do io
#         println(io, "H\tVN:Z:1.0")
#         for vertex in Graphs.vertices(graph)
#             if haskey(graph.vprops[vertex], :kmer)
#                 sequence = graph.vprops[vertex][:kmer]
#             else
#                 sequence = graph.vprops[vertex][:sequence]
#             end
# #             depth = graph.vprops[vertex][:weight]
#             depth = graph.vprops[vertex][:count]
#             fields = ["S", "$vertex", sequence, "DP:f:$(depth)"]
#             line = join(fields, '\t')
#             println(io, line)
#         end
#         for edge in Graphs.edges(graph)
#             overlap = graph.gprops[:k] - 1
#             for o in graph.eprops[edge][:orientations]
# #                 if !(!o.source_orientation && !o.destination_orientation)
#                 link = ["L",
#                             edge.src,
#                             o.source_orientation ? '+' : '-',
#                             edge.dst,
#                             o.destination_orientation ? '+' : '-',
#                             "$(overlap)M"]
#                 line = join(link, '\t')
#                 println(io, line)
# #                 end
#             end
#         end
#     end
#     return outfile
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function my_plot(graph::KmerGraph)
#     graph_hash = hash(sort(graph.graph.fadjlist), hash(graph.graph.ne))
#     filename = "/assets/images/$(graph_hash).svg"
#     p = plot_graph(graph)
#     Plots.savefig(p, dirname(pwd()) * filename)
#     display(p)
#     display("text/markdown", "![]($filename)")
# end