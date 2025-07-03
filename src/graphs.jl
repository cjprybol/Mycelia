"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create a weighted, strand-specific kmer (de bruijn) graph from a set of kmers
and a series of sequence observations in FASTA format.
"""
function build_stranded_kmer_graph(kmer_type, observations::AbstractVector{<:Union{FASTX.FASTA.Record, FASTX.FASTQ.Record}})
    
    canonical_kmer_counts = Mycelia.count_canonical_kmers(kmer_type, observations)
    canonical_kmers = collect(keys(canonical_kmer_counts))
    # @show canonical_kmers
    
    # if isempty(canonical_kmers)
    #     @error "isempty(canonical_kmers) = $(isempty(canonical_kmers))"
    # elseif isempty(observations)
    #     @error "isempty(observations) = $(isempty(observations))"
    # end
    stranded_kmers = sort!(vcat(canonical_kmers, [BioSequences.reverse_complement(kmer) for kmer in canonical_kmers]))
    stranded_kmer_to_reverse_complement_map = [
        findfirst(stranded_kmer -> BioSequences.reverse_complement(stranded_kmer) == kmer, stranded_kmers) for kmer in stranded_kmers
    ]
    stranded_kmer_graph = MetaGraphs.MetaDiGraph(length(stranded_kmers))
    stranded_kmer_graph.gprops[:stranded_kmers] = stranded_kmers
    stranded_kmer_graph.gprops[:reverse_complement_map] = stranded_kmer_to_reverse_complement_map
    stranded_kmer_graph.gprops[:k] = length(first(stranded_kmers))
    stranded_kmer_graph.gprops[:K] = length(stranded_kmers)
    stranded_kmer_graph.gprops[:observation_color_map] = Vector{Int}()
    stranded_kmer_graph.gprops[:observation_ids] = Vector{String}()
    stranded_kmer_graph.gprops[:observed_paths] = Vector{Vector{Pair{Int, Bool}}}()
    for vertex in 1:Graphs.nv(stranded_kmer_graph)
        stranded_kmer_graph.vprops[vertex] = Dict(:coverage => Vector{Pair{Int, Pair{Int, Bool}}}())
    end
    for (observation_index, observation) in enumerate(observations)
        observation_id = FASTX.FASTA.identifier(observation)
        observed_sequence = FASTX.FASTA.sequence(observation)
        if length(observed_sequence) < stranded_kmer_graph.gprops[:k]
            @error "skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        else
            observed_path = sequence_to_stranded_path(stranded_kmer_graph.gprops[:stranded_kmers], observed_sequence)
            i = 1
            ui, ui_orientation = observed_path[i]
            ui_coverage = (observation_index => (i => ui_orientation ))
            push!(stranded_kmer_graph.vprops[ui][:coverage], ui_coverage)
            for i in 2:length(observed_path)
                vi, vi_orientation = observed_path[i]
                vi_coverage = (observation_index => (i => vi_orientation))
                push!(stranded_kmer_graph.vprops[vi][:coverage], vi_coverage)
                edge_coverage = ui_coverage => vi_coverage
                if Graphs.has_edge(stranded_kmer_graph, ui, vi)
                    push!(stranded_kmer_graph.eprops[Graphs.Edge(ui, vi)][:coverage], edge_coverage)
                else
                    Graphs.add_edge!(stranded_kmer_graph, ui, vi, Dict(:coverage => [edge_coverage]))
                end
                # not sure this is necessary
#                 ui′ = stranded_kmer_graph.gprops[:reverse_complement_map][ui]
#                 vi′ = stranded_kmer_graph.gprops[:reverse_complement_map][vi]
#                 if !LightGraphs.has_edge(stranded_kmer_graph, vi′, ui′)
#                     LightGraphs.add_edge!(stranded_kmer_graph, vi′, ui′, Dict(:coverage => Vector{typeof(edge_coverage)}()))
#                 end
                ui, ui_orientation = vi, vi_orientation
                ui_coverage = vi_coverage
            end
            push!(stranded_kmer_graph.gprops[:observed_paths], observed_path)
            push!(stranded_kmer_graph.gprops[:observation_ids], observation_id)
            push!(stranded_kmer_graph.gprops[:observation_color_map], observation_index)
        end
    end
    @info Graphs.nv(stranded_kmer_graph)
    @info Graphs.ne(stranded_kmer_graph)
    return stranded_kmer_graph
end