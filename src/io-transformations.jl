# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
function graph_to_gfa(graph, outfile)
    open(outfile, "w") do io
        println(io, "H\tVN:Z:1.0")
        for vertex in Graphs.vertices(graph)
            if haskey(graph.vprops[vertex], :kmer)
                sequence = graph.vprops[vertex][:kmer]
            else
                sequence = graph.vprops[vertex][:sequence]
            end
#             depth = graph.vprops[vertex][:weight]
            depth = graph.vprops[vertex][:count]
            fields = ["S", "$vertex", sequence, "DP:f:$(depth)"]
            line = join(fields, '\t')
            println(io, line)
        end
        for edge in Graphs.edges(graph)
            overlap = graph.gprops[:k] - 1
            for o in graph.eprops[edge][:orientations]
#                 if !(!o.source_orientation && !o.destination_orientation)
                link = ["L",
                            edge.src,
                            o.source_orientation ? '+' : '-',
                            edge.dst,
                            o.destination_orientation ? '+' : '-',
                            "$(overlap)M"]
                line = join(link, '\t')
                println(io, line)
#                 end
            end
        end
    end
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
function open_fastx(path::String)
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
function list_databases(;address, username, password)
    cmd = "show databases"
    database = "system"
    cmd = cypher(;address, username, password, database, cmd)
    return DataFrames.DataFrame(uCSV.read(open(cmd), header=1, quotes='"', encodings=Dict("FALSE" => false, "TRUE" => true))...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
function create_database(;database, address, username, password)
    current_databases = list_databases(;address, username, password)
    if database in current_databases[!, "name"]
        return
    else
        f = run
        cmd = "create database $(database)"
        # switch database to system, so that we can create the user-specific database in the system
        database = "system"
        run(cypher(;address, username, password, database, cmd, f))
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
function cypher(;address, username, password, database, cmd)
    return `cypher-shell --address $address --username $username --password $password --database $(database) --format auto $(cmd)`
end

# TODO: uncomment if PR gets merged or switch to GenomicAnnotations
# https://github.com/BioJulia/GFF3.jl/pull/12
# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function open_gff(path::String)
#     if isfile(path)
#         io = open(path)
#     elseif occursin(r"^ftp", path) || occursin(r"^http", path)
#         path = replace(path, r"^ftp:" => "http:")
#         io = IOBuffer(HTTP.get(path).body)
#     else
#         error("unable to locate file $path")
#     end
#     path_base = basename(path)
#     if occursin(r"\.gz$", path_base)
#         io = CodecZlib.GzipDecompressorStream(io)
# #         path_base = replace(path_base, ".gz" => "")
#     end
#     gff_io = GFF3.Reader(io)
#     return gff_io
# end

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

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function my_plot(graph::Mycelia.KmerGraph)
#     graph_hash = hash(sort(graph.graph.fadjlist), hash(graph.graph.ne))
#     filename = "/assets/images/$(graph_hash).svg"
#     p = Mycelia.plot_graph(graph)
#     Plots.savefig(p, dirname(pwd()) * filename)
#     display(p)
#     display("text/markdown", "![]($filename)")
# end

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