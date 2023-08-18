function update_gff_with_mmseqs(gff_file, mmseqs_file)
    mmseqs_results = DataFrames.DataFrame(uCSV.read(mmseqs_file, header=1, delim='\t')...)

    gdf = DataFrames.groupby(mmseqs_results, "query")
    for g in gdf
        @assert issorted(g[!, "evalue"])
    end
    top_hits = DataFrames.combine(gdf, first)

    id_to_product = Dict()
    for row in DataFrames.eachrow(top_hits)
        id = Dict(a => b for (a, b) in split.(split(last(split(row["qheader"], " # ")), ';'), '='))["ID"]
        product = row["theader"]
        id_to_product[id] = product
    end

    gff_table = Mycelia.GFF_to_table(gff_file)
    for (i, row) in enumerate(DataFrames.eachrow(gff_table))
        id = Dict(a => b for (a,b) in split.(split(row["attributes"], " "), "="))["ID"]
        product = get(id_to_product, id, "")
        gff_table[i, "attributes"] = "product=\"$(product)\" " * row["attributes"]
        # gff_table[i, "attributes"] = "product=" * replace(id_to_product[id], " "=>"__") *  " $(row["attributes"])"
    end
    return gff_table
end

function GFF_to_table(gff_file)
    gff_columns = [
        :seqid,
        :source,
        :featuretype,
        :seqstart,
        :seqend,
        :score,
        :strand,
        :phase,
        :attributes
    ]
    table = DataFrames.DataFrame([[] for x in gff_columns], gff_columns)
    for record in collect(Mycelia.open_gff(gff_file))
        row = Dict(
            :seqid => GFF3.seqid(record),
            :source => GFF3.source(record),
            :featuretype => GFF3.featuretype(record),
            :seqstart => GFF3.seqstart(record),
            :seqend => GFF3.seqend(record),
            :score => GFF3.score(record),
            :strand => GFF3.strand(record),
            :phase => GFF3.phase(record),
            :attributes => GFF3.attributes(record),
        )
        row[:attributes] = join(["$(attribute)=$(join(values, ','))" for (attribute, values) in row[:attributes]], " ")
        push!(table, row)
    end
    return table
end

function read_kraken_report(kraken_report)
    kraken_report_header = [
        "percentage_of_fragments_at_or_below_taxon",
        "number_of_fragments_at_or_below_taxon",
        "number_of_fragments_assigned_directly_to_taxon",
        "rank",
        "ncbi_taxonid",
        "scientific_name"
    ]

    data, header = uCSV.read(kraken_report, delim='\t')
    kraken_report_table = DataFrames.DataFrame(data, kraken_report_header)
    return kraken_report_table
end

function diamond_line_to_named_tuple(diamond_line)
    sline = split(line)
    values_named_tuple = (
        qseqid = sline[1],
        sseqid = sline[2],
        pident = parse(Float64, sline[3]),
        length = parse(Int, sline[4]),
        mismatch = parse(Int, sline[5]),
        gapopen = parse(Int, sline[6]),
        qlen = parse(Int, sline[7]),
        qstart = parse(Int, sline[8]),
        qend = parse(Int, sline[9]),
        slen = parse(Int, sline[10]),
        sstart = parse(Int, sline[11]),
        send = parse(Int, sline[12]),
        evalue = parse(Float64, sline[13]),
        bitscore = parse(Float64, sline[14])
        )
    return values_named_tuple
end

function read_diamond_alignments_file(diamond_file)
    column_names_to_types = [
        "qseqid" => String,
        "sseqid" => String,
        "pident" => Float64,
        "length" => Int,
        "mismatch" => Int,
        "gapopen" => Int,
        "qlen" => Int,
        "qstart" => Int,
        "qend" => Int,
        "slen" => Int,
        "sstart" => Int,
        "send" => Int,
        "evalue" => Float64,
        "bitscore" => Float64,
    ]
    types = Dict(i => t for (i, t) in enumerate(last.(column_names_to_types)))
    
    data, header = uCSV.read(diamond_file, header=0, delim='\t', types = types)
    header = first.(column_names_to_types)    
    
    # data, header = uCSV.read(diamond_file, header=1, delim='\t', types = types)
    # @assert header == first.(column_names_to_types)
    
    table = DataFrames.DataFrame(data, header)
    return table
end

function add_header_to_diamond_file(infile, outfile=replace(infile, ".tsv" => ".with-header.tsv"))
    column_names = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qlen",
        "qstart",
        "qend",
        "slen",
        "sstart",
        "send",
        "evalue",
        "bitscore"
    ]
    # dangerous but fast
    # try
    #     inserted_text = join(columns_names, '\t') * '\n'
    #     sed_cmd = "1s/^/$(inserted_text)/"
    #     full_cmd = `sed -i $sed_cmd $infile`
    # catch
    open(outfile, "w") do io
        println(io, join(column_names, "\t"))
        for line in eachline(infile)
            println(io, line)
        end
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse the contig coverage information from qualimap bamqc text report, which looks like the following:

```
# this is spades
>>>>>>> Coverage per contig

	NODE_1_length_107478_cov_9.051896	107478	21606903	201.0355886786133	60.39424208607496
	NODE_2_length_5444_cov_1.351945	5444	153263	28.152645113886848	5.954250612823136
	NODE_3_length_1062_cov_0.154390	1062	4294	4.043314500941619	1.6655384692688975
	NODE_4_length_776_cov_0.191489	776	3210	4.13659793814433	2.252009588980858

# below is megahit
>>>>>>> Coverage per contig

	k79_175	235	3862	16.43404255319149	8.437436249612457
	k79_89	303	3803	12.551155115511552	5.709975376279777
	k79_262	394	6671	16.931472081218274	7.579217802849293
	k79_90	379	1539	4.060686015831134	1.2929729111266581
	k79_91	211	3749	17.767772511848342	11.899185693011933
	k79_0	2042	90867	44.49902056807052	18.356525483516613
```

To make this more robust, consider reading in the names of the contigs from the assembled fasta
"""
function parse_qualimap_contig_coverage(qualimap_report_txt)
    coverage_line_regex = r"\t.*?\t\d+\t\d+\t[\d\.]+\t[\d\.]+$"
    lines = filter(x -> occursin(coverage_line_regex, x), readlines("$(qualimap_report_txt)"))
    io = IOBuffer(join(map(x -> join(split(x, '\t')[2:end], '\t'), lines), '\n'))
    header = ["Contig", "Length", "Mapped bases", "Mean coverage", "Standard Deviation"]
    types = [String, Int, Int, Float64, Float64]
    data, _ = uCSV.read(io, delim='\t', types=types)
    qualimap_results = DataFrames.DataFrame(data, header)
    qualimap_results[!, "% Mapped bases"] = qualimap_results[!, "Mapped bases"] ./ sum(qualimap_results[!, "Mapped bases"]) .* 100
    return qualimap_results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Imports results of fastani

```jldoctest
julia> 1 + 1
2
```
"""
function read_fastani(path::String)
    data, header = uCSV.read(path, delim='\t', typedetectrows=100)
    header = [
        "query_identifier",
        "reference_identifier",
        "%_identity",
        "fragments_mapped",
        "total_query_fragments"
    ]

    ani_table = DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Imports results of Diamond (or blast) in outfmt 6 as a DataFrame

```jldoctest
julia> 1 + 1
2
```
"""
function read_diamond(path::String)
  diamond_colnames = [ "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore" ]
  
  diamond_coltypes = Dict(
     1 => String, 
     2 => String, 
     3 => Float64,
     4 => Int,
     5 => Int,
     6 => Int,
     7 => Int,
     8 => Int,
     9 => Int,
     10 => Int,
     11 => String,
     12 => String
  )
    return DataFrames.DataFrame(uCSV.read(open(path), delim = '\t', header = diamond_colnames, types = diamond_coltypes)...)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Expects output type 7 from BLAST, default output type 6 doesn't have the header comments and won't auto-parse

```jldoctest
julia> 1 + 1
2
```
"""
function parse_blast_report(blast_report)
    # example header line 
    # "# Fields: query id, subject id, subject acc., subject acc.ver, subject title, query length, subject length, q. start, q. end, s. start, s. end, evalue, bit score, score, alignment length, % identity, identical, mismatches, subject tax id"
    header_lines = collect(Iterators.filter(x -> occursin(r"# Fields:", x), eachline(blast_report)))
    if isempty(header_lines)
        @info "not hits found, returning empty table"
        return DataFrames.DataFrame()
    end
    header_line = first(header_lines)
    header = split(last(split(header_line, ": ")), ", ")
    blast_col_types = Dict(
        "query id" => String,
        "query title" => String,
        "subject id" => String,
        "subject gi" => String,
        "subject acc." => String,
        "subject acc.ver" => String,
        "subject title" => String,
        "query length" => Int,
        "subject length" => Int,
        "q. start" => Int,
        "q. end" => Int,
        "s. start" => Int,
        "s. end" => Int,
        "evalue" => Float64,
        "bit score" => Float64,
        "score" => Float64,
        "alignment length" => Int,
        "% identity" => Float64,
        "identical" => Int,
        "mismatches" => Int,
        "subject tax id" => Int,
        "subject sci name" => String,
        "subject com names" => String,
        "subject blast name" => String,
        "subject super kingdom" => String,
        "subject tax ids" => String,
        "subject sci names" => String,
        "subject com names" => String,
        "subject blast names" => String,
        "subject super kingdoms" => String,
        "subject title" => String,
        "subject titles" => String
    )
    data, _ = uCSV.read(
        blast_report,
        delim='\t',
        comment='#',
        skipmalformed=true,
        allowmissing=true,
        encodings=Dict("N/A" => missing),
        types=[blast_col_types[h] for h in header])
    return DataFrames.DataFrame(data, header, makeunique=true)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a GFA file into a genome graph - need to finish implementation and assert contig normalization (i.e. is canonical) before using with my code

```jldoctest
julia> 1 + 1
2
```
"""
function parse_gfa(gfa)
    
    gfa_record_types = Dict(
        '#' => "Comment",
        'H' => "Header",
        'S' => "Segment",
        'L' => "Link",
        'J' => "Jump",
        'C' => "Containment",
        'P' => "Path",
        'W' => "Walk"
    )

    gfa_graph = MetaGraphs.MetaDiGraph()
    MetaGraphs.set_prop!(gfa_graph, :paths, Dict{String, Any}())
    for line in eachline(gfa)
        record_type = gfa_record_types[line[1]]
        if record_type == "Header"
            # metadata
            sline = split(line)
            # add me later
        elseif record_type == "Comment"
            # metadata
            # add me later
        elseif record_type == "Segment"
            # node
            record_type, record_name, sequence = split(line, '\t')
            Graphs.add_vertex!(gfa_graph)
            node_index = Graphs.nv(gfa_graph)
            MetaGraphs.set_prop!(gfa_graph, node_index, :identifier, record_name)
            MetaGraphs.set_indexing_prop!(gfa_graph, :identifier)
            MetaGraphs.set_prop!(gfa_graph, node_index, :sequence, sequence)
        elseif record_type == "Link"
            record_type, source_identifier, source_orientation, destination_identifier, destination_orientation, overlap_CIGAR = split(line, '\t')
            source_index = gfa_graph[source_identifier, :identifier]
            destination_index = gfa_graph[destination_identifier, :identifier]
            edge = Graphs.Edge(source_index, destination_index)
            Graphs.add_edge!(gfa_graph, edge)
            MetaGraphs.set_prop!(gfa_graph, edge, :source_identifier, source_identifier)
            MetaGraphs.set_prop!(gfa_graph, edge, :source_orientation, source_orientation)
            MetaGraphs.set_prop!(gfa_graph, edge, :destination_identifier, destination_identifier)
            MetaGraphs.set_prop!(gfa_graph, edge, :destination_orientation, destination_orientation)
            MetaGraphs.set_prop!(gfa_graph, edge, :overlap_CIGAR, overlap_CIGAR)
        elseif record_type == "Path"
            record_type, path_identifier, segments, overlaps = split(line, '\t')
            gfa_graph.gprops[:paths][path_identifier] = Dict("segments" => segments, "overlaps" => overlaps)
        else
            @warn "GFA line type $(record_type) not currently handled by the import - please add"
        end
    end
    return gfa_graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

A short description of the function

```jldoctest
julia> 1 + 1
2
```
"""
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