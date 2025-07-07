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

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function path_to_sequence(kmer_graph, path)
#     sequence = BioSequences.LongDNASeq(oriented_kmer_to_sequence(kmer_graph, first(path)))
#     for oriented_kmer in path[2:end]
#         nucleotide = last(oriented_kmer_to_sequence(kmer_graph, oriented_kmer))
#         push!(sequence, nucleotide)
#     end
#     return sequence
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function oriented_kmer_to_sequence(kmer_graph, oriented_kmer)
#     kmer_sequence = kmer_graph.kmers[oriented_kmer.index]
#     if !oriented_kmer.orientation
#         kmer_sequence = BioSequences.reverse_complement(kmer_sequence)
#     end
#     return kmer_sequence
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a DNA sequence into a path through a collection of stranded k-mers.

# Arguments
- `stranded_kmers`: Collection of unique k-mers representing possible path vertices
- `sequence`: Input DNA sequence to convert to a path

# Returns
Vector of `Pair{Int,Bool}` where:
- First element (Int) is the index of the k-mer in `stranded_kmers`
- Second element (Bool) indicates orientation (true=forward, false=reverse)
"""
function sequence_to_stranded_path(stranded_kmers, sequence)
    KMER_TYPE = typeof(first(stranded_kmers))
    path = Vector{Pair{Int, Bool}}()
    for (i, kmer) in Kmers.EveryKmer{KMER_TYPE}(BioSequences.LongDNA{4}(sequence))
        kmer_index = findfirst(stranded_kmer -> kmer == stranded_kmer, stranded_kmers)
        orientation = true
        push!(path, kmer_index => orientation)
    end
    return path
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a path through k-mers into a single DNA sequence.

Takes a vector of k-mers and a path representing the order to traverse them,
reconstructs the original sequence by joining the k-mers according to the path.
The first k-mer is used in full, then only the last nucleotide from each subsequent k-mer is added.

# Arguments
- `kmers`: Vector of DNA k-mers (as LongDNA{4})
- `path`: Vector of tuples representing the path through the k-mers

# Returns
- `LongDNA{4}`: The reconstructed DNA sequence
"""
function path_to_sequence(kmers, path)
    # @show path
    sequence = BioSequences.LongDNA{4}(kmers[first(first(path))])
    for i in 2:length(path)
        push!(sequence, kmers[first(path[i])][end])
    end
    return sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the probability of traversing a specific edge in a stranded k-mer graph.

The probability is computed as the ratio of this edge's coverage weight to the sum
of all outgoing edge weights from the source vertex.

$(DocStringExtensions.TYPEDSIGNATURES)

# Arguments
- `stranded_kmer_graph`: A directed graph where edges represent k-mer connections
- `edge`: The edge for which to calculate the probability

# Returns
- `Float64`: Probability in range [0,1] representing likelihood of traversing this edge
  Returns 0.0 if sum of all outgoing edge weights is zero

# Note
Probability is based on the :coverage property of edges, using their length as weights
"""
function edge_probability(stranded_kmer_graph, edge)
    neighbors = Graphs.outneighbors(stranded_kmer_graph, edge.src)
    neighbor_under_consideration = findfirst(neighbor -> neighbor == edge.dst, neighbors)
    edge_weights = [length(stranded_kmer_graph.eprops[Graphs.Edge(edge.src, neighbor)][:coverage]) for neighbor in neighbors]
    if sum(edge_weights) == 0
        p = 0.0
    else
        edge_probabilities = edge_weights ./ sum(edge_weights)
        p = edge_probabilities[neighbor_under_consideration]
    end
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Converts a path of edges in a kmer graph into a DNA sequence by concatenating overlapping kmers.

# Arguments
- `kmer_graph`: A directed graph where vertices represent kmers and edges represent overlaps
- `edge_path`: Vector of edges representing a path through the graph

# Returns
A `BioSequences.LongDNASeq` containing the merged sequence represented by the path

# Details
The function:
1. Takes the first kmer from the source vertex of first edge
2. For each edge, handles orientation (forward/reverse complement)
3. Verifies overlaps between consecutive kmers
4. Concatenates unique bases to build final sequence
"""
function edge_path_to_sequence(kmer_graph, edge_path)
    edge = first(edge_path)
    sequence = BioSequences.LongDNASeq(kmers[edge.src])
    if !kmer_graph.eprops[edge][:orientations].source_orientation
        sequence = BioSequences.reverse_complement(sequence)
    end
    for edge in edge_path
        destination = BioSequences.LongDNASeq(kmers[edge.dst])
        if !kmer_graph.eprops[edge][:orientations].destination_orientation
            destination = BioSequences.reverse_complement(destination)
        end
        sequence_suffix = sequence[end-length(destination)+2:end]
        destination_prefix = destination[1:end-1]
        @assert sequence_suffix == destination_prefix
        push!(sequence, destination[end])
    end
    sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns the index position of a given k-mer in a sorted list of k-mers.

# Arguments
- `kmers`: A sorted vector of k-mers to search within
- `kmer`: The k-mer sequence to find

# Returns
Integer index position where `kmer` is found in `kmers`

# Throws
- `AssertionError`: If the k-mer is not found in the list
"""
function get_kmer_index(kmers, kmer)
    index = searchsortedfirst(kmers, kmer)
    @assert kmers[index] == kmer "$kmer not found in kmer list"
    return index
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a path of overlapping k-mers into a single DNA sequence.

# Arguments
- `kmer_path`: Vector of k-mers (DNA sequences) where each consecutive pair overlaps by k-1 bases

# Returns
- `BioSequences.LongDNA{2}`: Assembled DNA sequence from the k-mer path

# Description
Reconstructs the original DNA sequence by joining k-mers, validating that consecutive k-mers 
overlap correctly. The first k-mer is used in full, then each subsequent k-mer contributes 
its last base.
"""
function kmer_path_to_sequence(kmer_path)
    sequence = BioSequences.LongDNA{2}(first(kmer_path))
    for kmer in kmer_path[2:end]
        for i in 1:length(kmer)-1
            a = kmer[i]
            b = sequence[end-(length(kmer)-1)+i]
            @assert a == b
        end
        push!(sequence, kmer[end])
    end
    return sequence
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a GFA (Graphical Fragment Assembly) file to FASTA format.

# Arguments
- `gfa::String`: Path to input GFA file
- `fasta::String=gfa * ".fna"`: Path for output FASTA file. Defaults to input filename with ".fna" extension

# Returns
- `String`: Path to the generated FASTA file

# Details
Uses gfatools (via Conda) to perform the conversion. The function will:
1. Ensure gfatools is available in the Conda environment
2. Execute the conversion using gfatools gfa2fa
3. Write sequences to the specified FASTA file
"""
function gfa_to_fasta(;gfa, fasta=gfa * ".fna")
    Mycelia.add_bioconda_env("gfatools")
    p = pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n gfatools gfatools gfa2fa $(gfa)`, fasta)
    run(p)
    return fasta
    # collect(GraphicalFragmentAssembly.Reader(open(primary_contig_gfa)))
    # open(fasta, "w") do io
    #     fastx_io = FASTX.FASTA.Writer(io)
    #     gfa_graph = Mycelia.parse_gfa(gfa)
    #     for v in Graphs.vertices(gfa_graph)
    #         record = FASTX.FASTA.Record(gfa_graph.vprops[v][:identifier], gfa_graph.vprops[v][:sequence])
    #         write(fastx_io, record)
    #     end
    #     close(fastx_io)
    # end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Identifies sequence regions that require resampling based on kmer solidity patterns.

# Arguments
- `record_kmer_solidity::BitVector`: Boolean array where `true` indicates solid kmers
- `solid_branching_kmer_indices::Vector{Int}`: Indices of solid branching kmers

# Returns
- `Vector{UnitRange{Int64}}`: Array of ranges (start:stop) indicating stretches that need resampling

# Details
Finds continuous stretches of non-solid kmers and extends them to the nearest solid branching
kmers on either side. These stretches represent regions that need resampling.

If a stretch doesn't have solid branching kmers on both sides, it is excluded from the result.
Duplicate ranges are removed from the final output.
"""
function find_resampling_stretches(;record_kmer_solidity, solid_branching_kmer_indices)
    indices = findall(.!record_kmer_solidity)  # Find the indices of false values
    if isempty(indices)
        return UnitRange{Int64}[]
    end
    
    diffs = diff(indices)  # Calculate the differences between consecutive indices
    # @show diffs
    range_starts = [indices[1]]  # Start with the first false index
    range_ends = Int[]
    
    for (i, d) in enumerate(diffs)
        if d > 1
            push!(range_ends, indices[i])
            push!(range_starts, indices[i+1])
        end
    end
    
    push!(range_ends, indices[end])  # Add the last false index as a range end
    
    low_quality_runs = [(start, stop) for (start, stop) in zip(range_starts, range_ends)]
    
    resampling_stretches = UnitRange{Int64}[]
    
    for low_quality_run in low_quality_runs
        unders = filter(solid_branching_kmer -> solid_branching_kmer < first(low_quality_run), solid_branching_kmer_indices)
        overs = filter(solid_branching_kmer -> solid_branching_kmer > last(low_quality_run), solid_branching_kmer_indices)
        if isempty(overs) || isempty(unders)
            continue
        else
            nearest_under = maximum(unders)
            nearest_over = minimum(overs)
            push!(resampling_stretches, nearest_under:nearest_over)
        end
    end
    if !allunique(resampling_stretches)
        resampling_stretches = unique!(resampling_stretches)
    end
    return resampling_stretches
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Constructs a directed graph representation of k-mer transitions from FASTQ sequencing data.

# Arguments
- `fastq`: Path to input FASTQ file
- `k`: K-mer size (default: 1). Must be odd and prime. If k=1, optimal size is auto-determined
- `plot`: Boolean to display quality distribution plot (default: false)

# Returns
MetaDiGraph with properties:
- assembly_k: k-mer size used
- kmer_counts: frequency of each k-mer
- transition_likelihoods: edge weights between k-mers
- kmer_mean_quality, kmer_total_quality: quality metrics
- branching_nodes, unbranching_nodes: topological classification
- likely_valid_kmer_indices: k-mers above mean quality threshold
- likely_sequencing_artifact_indices: potential erroneous k-mers

# Note
For DNA assembly, quality scores are normalized across both strands.
"""
function build_directed_kmer_graph(;fastq, k=1, plot=false)
    if k == 1
        assembly_k = Mycelia.assess_dnamer_saturation([fastq])
    else
        @assert isodd(k)
        @assert Primes.isprime(k)
        assembly_k = k
    end
    kmer_type = Kmers.DNAKmer{assembly_k}

    # initializing the graph with kmer counts
    kmer_counts = Mycelia.count_kmers(kmer_type, fastq)
    ordered_kmers = collect(keys(kmer_counts))
    total_states = length(ordered_kmers)
    graph = MetaGraphs.MetaDiGraph(total_states)
    MetaGraphs.set_prop!(graph, :assembly_k, assembly_k)
    MetaGraphs.set_prop!(graph, :kmer_counts, kmer_counts)
    MetaGraphs.set_prop!(graph, :total_states, total_states)
    MetaGraphs.set_prop!(graph, :ordered_kmers, ordered_kmers)
    kmer_indices = sort(Dict(kmer => i for (i, kmer) in enumerate(keys(kmer_counts))))
    MetaGraphs.set_prop!(graph, :kmer_indices, kmer_indices)
    canonical_kmer_counts = Mycelia.count_canonical_kmers(kmer_type, fastq)
    MetaGraphs.set_prop!(graph, :canonical_kmer_counts, canonical_kmer_counts)
    canonical_kmer_indices = sort(Dict(kmer => i for (i, kmer) in enumerate(keys(canonical_kmer_counts))))
    MetaGraphs.set_prop!(graph, :canonical_kmer_indices, canonical_kmer_indices)
    
    
    # kmer quality and likelihoods
    
    records = collect(Mycelia.open_fastx(fastq))
    read_quality_scores = [collect(FASTX.quality_scores(record)) for record in records]
    all_kmer_quality_support = Dict{kmer_type, Vector{Float64}}()
    for record in records
        record_quality_scores = collect(FASTX.quality_scores(record))
        record_quality_score_slices = [record_quality_scores[i:i+assembly_k-1] for i in 1:length(record_quality_scores)-assembly_k+1]
        sequence = BioSequences.LongDNA{2}(FASTX.sequence(record))
        for ((i, kmer), kmer_base_qualities) in zip(Kmers.EveryKmer{kmer_type}(sequence), record_quality_score_slices)
            if haskey(all_kmer_quality_support, kmer)
                all_kmer_quality_support[kmer] = all_kmer_quality_support[kmer] .+ kmer_base_qualities
            else
                all_kmer_quality_support[kmer] = kmer_base_qualities
            end
        end
    end
    
    # strand normalization shares observational quality across strands - only relevant for non-stranded DNA genome assembly
    strand_normalized_quality_support = Dict{kmer_type, Vector{Float64}}()
    for (kmer, support) in all_kmer_quality_support
        strand_normalized_quality_support[kmer] = support
        if haskey(all_kmer_quality_support, BioSequences.reverse_complement(kmer))
            strand_normalized_quality_support[kmer] .+= all_kmer_quality_support[BioSequences.reverse_complement(kmer)]
        end
    end
    strand_normalized_quality_support
    kmer_mean_quality = sort(Dict(kmer => strand_normalized_quality_support[kmer] ./ canonical_kmer_counts[BioSequences.canonical(kmer)] for kmer in ordered_kmers))
    MetaGraphs.set_prop!(graph, :kmer_mean_quality, kmer_mean_quality)
    kmer_total_quality = sort(Dict(kmer => sum(quality_values) for (kmer, quality_values) in strand_normalized_quality_support))
    MetaGraphs.set_prop!(graph, :kmer_total_quality, kmer_total_quality)
    state_likelihoods = sort(Dict(kmer => total_quality / sum(values(kmer_total_quality)) for (kmer, total_quality) in kmer_total_quality))
    MetaGraphs.set_prop!(graph, :state_likelihoods, state_likelihoods)


    # all transition likelihood calculation
    transition_likelihoods = SparseArrays.spzeros(total_states, total_states)
    for record in records
        sequence = BioSequences.LongDNA{4}(FASTX.sequence(record))
        sources = Kmers.EveryKmer{kmer_type}(sequence[1:end-1])
        destinations = Kmers.EveryKmer{kmer_type}(sequence[2:end])
        for ((source_i, source), (destination_i, destination)) in zip(sources, destinations)
            source_index = kmer_indices[source]
            destination_index = kmer_indices[destination]
            transition_likelihoods[source_index, destination_index] += 1
        end
    end
    for source in 1:total_states
        outgoing_transition_counts = transition_likelihoods[source, :]
        if sum(outgoing_transition_counts) > 0
            transition_likelihoods[source, :] .= transition_likelihoods[source, :] ./ sum(transition_likelihoods[source, :]) 
        end
    end
    row_indices, column_indices, cell_values = SparseArrays.findnz(transition_likelihoods)
    for (row, col, value) in zip(row_indices, column_indices, cell_values)
        Graphs.add_edge!(graph, row, col)
        MetaGraphs.set_prop!(graph, row, col, :transition_likelihood, value)
    end
    MetaGraphs.set_prop!(graph, :transition_likelihoods, transition_likelihoods)

    # helpful for downstream processing
    unbranching_nodes = Set(Int[])
    for node in Graphs.vertices(graph)
        if (Graphs.indegree(graph, node) <= 1) && (Graphs.outdegree(graph, node) <= 1)
            push!(unbranching_nodes, node)
        end
    end
    branching_nodes = Set(setdiff(Graphs.vertices(graph), unbranching_nodes))
    MetaGraphs.set_prop!(graph, :unbranching_nodes, unbranching_nodes)
    MetaGraphs.set_prop!(graph, :branching_nodes, branching_nodes)
    
    
    # total_strand_normalized_quality_support = sum.(collect(values(strand_normalized_quality_support)))
    mean_total_support = Statistics.mean(collect(values(kmer_total_quality)))
    sorted_kmer_total_quality_values = collect(values(kmer_total_quality))
    mean_quality_value = Statistics.mean(sorted_kmer_total_quality_values)
    threshold = mean_quality_value

    xs = [
        [i for (i, y) in enumerate(sorted_kmer_total_quality_values) if y > threshold],
        [i for (i, y) in enumerate(sorted_kmer_total_quality_values) if y <= threshold]
        ]
    
    likely_valid_kmer_indices = xs[1]
    MetaGraphs.set_prop!(graph, :likely_valid_kmer_indices, likely_valid_kmer_indices)
    likely_sequencing_artifact_indices = xs[2]
    MetaGraphs.set_prop!(graph, :likely_sequencing_artifact_indices, likely_sequencing_artifact_indices)
    # likely_sequencing_artifact_kmers = Set(ordered_kmers[likely_sequencing_artifact_indices])
    # likely_valid_kmers = Set(ordered_kmers[likely_valid_kmer_indices])
    # kmer_to_index_map = Dict(kmer => i for (i, kmer) in enumerate(ordered_kmers))
    
    
    if plot
        ys = [
            [y for y in sorted_kmer_total_quality_values if y > threshold],
            [y for y in sorted_kmer_total_quality_values if y <= threshold]
        ]

        p = StatsPlots.scatter(
            xs,
            ys,
            title = "kmer qualities",
            ylabel = "canonical kmer cumulative QUAL value",
            label = ["above" "below"],
            legend = :outertopright,
            # size = (900, 500),
            margins=10StatsPlots.Plots.PlotMeasures.mm,
            xticks = false
        )
        p = StatsPlots.hline!(p, [mean_quality_value], label="mean")
        display(p)
    end
    return graph
end

# not a very good function yet, but good enough for the pinches I need it for
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse a GFA (Graphical Fragment Assembly) file into a MetaGraph representation.

# Arguments
- `gfa`: Path to GFA format file

# Returns
A `MetaGraph` where:
- Vertices represent segments (contigs)
- Edges represent links between segments
- Vertex properties include `:id` with segment identifiers
- Graph property `:records` contains the original FASTA records

# Format Support
Handles standard GFA v1 lines:
- `H`: Header lines (skipped)
- `S`: Segments (stored as nodes with FASTA records)
- `L`: Links (stored as edges)
- `P`: Paths (stored in paths dictionary)
- `A`: HiFiAsm specific lines (skipped)
"""
function parse_gfa(gfa)
    segments = Vector{FASTX.FASTA.Record}()
    links = Vector{Pair{String, String}}()
    paths = Dict{String, Vector{String}}()
    for l in eachline(open(gfa))
        s = split(l, '\t')
        if first(s) == "H"
            # header line
            continue
        elseif first(s) == "S"
            # segment
            # push!(segments, string(s[2]))
            identifier = string(s[2])
            description = string(s[4])
            sequence = string(s[3])
            push!(segments, FASTX.FASTA.Record("$(identifier) $(description)", sequence))
        elseif first(s) == "L"
            # link
            push!(links, string(s[2]) => string(s[4]))
        elseif first(s) == "P"
            # path
            paths[string(s[2])] = string.(split(replace(s[3], r"[+-]" => ""), ','))
        elseif first(s) == "A" # hifiasm https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#output-file-formats
            continue
        else
            println(l)
            error("unexpected line encountered while parsing GFA")
        end
    end

    g = MetaGraphs.MetaGraph(length(segments))

    for link in links
        (u, v) = link
        ui = findfirst(FASTX.identifier.(segments) .== u)
        vi = findfirst(FASTX.identifier.(segments) .== v)
        Graphs.add_edge!(g, ui => vi)
    end
    for (i, segment) in enumerate(segments)
        MetaGraphs.set_prop!(g, i, :id, FASTX.identifier(segment))
    end
    MetaGraphs.set_prop!(g, :records, segments)
    return g
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a GFA (Graphical Fragment Assembly) file into a structured representation.

# Arguments
- `gfa`: Path to GFA file or GFA content as string

# Returns
Named tuple containing:
- `contig_table`: DataFrame with columns:
  - `connected_component`: Integer ID for each component
  - `contigs`: Comma-separated list of contig IDs
  - `is_circular`: Boolean indicating if component forms a cycle
  - `is_closed`: Boolean indicating if single contig forms a cycle
  - `lengths`: Comma-separated list of contig lengths
- `records`: FASTA records from the GFA
"""
function gfa_to_structure_table(gfa)
    gfa_metagraph = parse_gfa(gfa)
    contig_table = DataFrames.DataFrame()
    records = MetaGraphs.get_prop(gfa_metagraph, :records)
    contig_lengths = Dict(FASTX.identifier(record) => length(FASTX.sequence(record)) for record in records)
    # @show String.(FASTX.description.(records))
    # try
    # contig_depths = Dict(FASTX.identifier(record) => first(match(r"^.*?dp:i:(\d+).*$", String(FASTX.description(record))).captures) for record in records)
    for (i, connected_component) in enumerate(Graphs.connected_components(gfa_metagraph))
        subgraph, node_map = Graphs.induced_subgraph(gfa_metagraph, connected_component)
        # display(subgraph)
        contigs = [MetaGraphs.get_prop(subgraph, v, :id) for v in Graphs.vertices(subgraph)]
        row = (
            connected_component = i,
            contigs = join(contigs, ","),
            is_circular = Graphs.is_cyclic(subgraph),
            is_closed = (length(contigs) == 1) && Graphs.is_cyclic(subgraph),
            lengths = join([contig_lengths[contig] for contig in contigs], ","),
            # depths = join([contig_depths[contig] for contig in contigs], ","),
            )
        push!(contig_table, row)
    end
    
    return (;contig_table, records)
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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Saves the given graph to a file in JLD2 format.

# Arguments
- `graph::Graphs.AbstractGraph`: The graph to be saved.
- `outfile::String`: The name of the output file. If the file extension is not `.jld2`, it will be appended automatically.

# Returns
- `String`: The name of the output file with the `.jld2` extension.
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

Loads a graph object from a serialized file.

# Arguments
- `file::String`: Path to the file containing the serialized graph data. The file should have been created using `save_graph`.

# Returns
- The deserialized graph object stored under the "graph" key.
"""
function load_graph(file::String)
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
# function graph_to_gfa(graph, kmer_size, outfile="$(kmer_size).gfa")
#     kmer_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size})))
#     # add fastq here too???
#     record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
#     edgemer_edges = collect(MetaGraphs.filter_edges(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size+1})))
#     open(outfile, "w") do io
#         println(io, "H\tVN:Z:1.0")
#         for vertex in kmer_vertices
#             if haskey(graph.vprops[vertex], :kmer)
#                 sequence = graph.vprops[vertex][:kmer]
#             else
#                 sequence = graph.vprops[vertex][:sequence]
#             end
#             # depth = graph.vprops[vertex][:count]
#             # depth = length(graph.vprops[vertex][:evidence])
#             total_count = 0
#             vertex_outneighbors = Graphs.outneighbors(graph, vertex)
#             for connected_record in intersect(vertex_outneighbors, record_vertices)
#                 edge = Graphs.Edge(vertex, connected_record)
#                 total_count += MetaGraphs.get_prop(graph, edge, :count)
#             end

#             fields = ["S", "$vertex", sequence, "DP:f:$(total_count)"]
#             line = join(fields, '\t')
#             println(io, line)
#         end
#         overlap = kmer_size - 1
#         for edgemer_edge in edgemer_edges
#             source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
#             link = ["L",
#                         edgemer_edge.src,
#                         BioSequences.iscanonical(source_kmer) ? '+' : '-',
#                         edgemer_edge.dst,
#                         BioSequences.iscanonical(dest_kmer) ? '+' : '-',
#                         "$(overlap)M"]
#             line = join(link, '\t')
#             println(io, line)
#         end
#     end
#     return
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Parse a GFA file into a genome graph - need to finish implementation and assert contig normalization (i.e. is canonical) before using with my code
# """
# function parse_gfa(gfa)
    
#     # collect(GraphicalFragmentAssembly.Reader(open(primary_contig_gfa)))
    
#     gfa_record_types = Dict(
#         '#' => "Comment",
#         'H' => "Header",
#         'S' => "Segment",
#         'L' => "Link",
#         'J' => "Jump",
#         'C' => "Containment",
#         'P' => "Path",
#         'W' => "Walk"
#     )

#     gfa_graph = MetaGraphs.MetaDiGraph()
#     MetaGraphs.set_prop!(gfa_graph, :paths, Dict{String, Any}())
#     for line in eachline(gfa)
#         record_type = gfa_record_types[line[1]]
#         if record_type == "Header"
#             # metadata
#             sline = split(line)
#             # add me later
#         elseif record_type == "Comment"
#             # metadata
#             # add me later
#         elseif record_type == "Segment"
#             # node
#             record_type, record_name, sequence = split(line, '\t')
#             Graphs.add_vertex!(gfa_graph)
#             node_index = Graphs.nv(gfa_graph)
#             MetaGraphs.set_prop!(gfa_graph, node_index, :identifier, record_name)
#             MetaGraphs.set_indexing_prop!(gfa_graph, :identifier)
#             MetaGraphs.set_prop!(gfa_graph, node_index, :sequence, sequence)
#         elseif record_type == "Link"
#             record_type, source_identifier, source_orientation, destination_identifier, destination_orientation, overlap_CIGAR = split(line, '\t')
#             source_index = gfa_graph[source_identifier, :identifier]
#             destination_index = gfa_graph[destination_identifier, :identifier]
#             edge = Graphs.Edge(source_index, destination_index)
#             Graphs.add_edge!(gfa_graph, edge)
#             MetaGraphs.set_prop!(gfa_graph, edge, :source_identifier, source_identifier)
#             MetaGraphs.set_prop!(gfa_graph, edge, :source_orientation, source_orientation)
#             MetaGraphs.set_prop!(gfa_graph, edge, :destination_identifier, destination_identifier)
#             MetaGraphs.set_prop!(gfa_graph, edge, :destination_orientation, destination_orientation)
#             MetaGraphs.set_prop!(gfa_graph, edge, :overlap_CIGAR, overlap_CIGAR)
#         elseif record_type == "Path"
#             record_type, path_identifier, segments, overlaps = split(line, '\t')
#             gfa_graph.gprops[:paths][path_identifier] = Dict("segments" => segments, "overlaps" => overlaps)
#         else
#             @warn "GFA line type $(record_type) not currently handled by the import - please add"
#         end
#     end
#     return gfa_graph
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a Mycelia graph to GFA (Graphical Fragment Assembly) format.

Writes a graph to GFA format, including:
- Header (H) line with GFA version
- Segment (S) lines for each vertex with sequence and depth
- Link (L) lines for edges with overlap size and orientations

# Arguments
- `graph`: MetaGraph containing sequence vertices and their relationships
- `outfile`: Path where the GFA file should be written

# Returns
- Path to the written GFA file
"""
function graph_to_gfa(;graph, outfile)
    # kmer_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size})))
    # add fastq here too???
    # record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
    # edgemer_edges = collect(MetaGraphs.filter_edges(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size+1})))
    open(outfile, "w") do io
        println(io, "H\tVN:Z:1.0")
        for vertex in Graphs.vertices(graph)
            if haskey(graph.vprops[vertex], :kmer)
                sequence = graph.vprops[vertex][:kmer]
            else
                sequence = graph.vprops[vertex][:sequence]
            end
            depth = graph.vprops[vertex][:count]
            # depth = length(graph.vprops[vertex][:evidence])
            # total_count = 0
            # vertex_outneighbors = Graphs.outneighbors(graph, vertex)
            # for connected_record in vertex_outneighbors
            #     edge = Graphs.Edge(vertex, connected_record)
            #     total_count += MetaGraphs.get_prop(graph, edge, :count)
            # end

            fields = ["S", "$vertex", sequence, "DP:f:$(depth)"]
            line = join(fields, '\t')
            println(io, line)
        end
        overlap = MetaGraphs.get_prop(graph, :k) - 1
        for edge in collect(Graphs.edges(graph))
            # @show edge
            # @show MetaGraphs.props(graph, edge)
            unique_orientations = unique(observation.orientation for observation in MetaGraphs.get_prop(graph, edge, :observations))
            # @show unique_orientations
            source_kmer = edge.src
            destination_kmer = edge.dst
            for unique_orientation in unique_orientations
                # source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
                link = ["L",
                            edge.src,
                            first(unique_orientation) ? '+' : '-',
                            edge.dst,
                            last(unique_orientation) ? '+' : '-',
                            "$(overlap)M"]
                line = join(link, '\t')
                println(io, line)
            end
                
            # source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
            # link = ["L",
            #             edgemer_edge.src,
            #             BioSequences.iscanonical(source_kmer) ? '+' : '-',
            #             edgemer_edge.dst,
            #             BioSequences.iscanonical(dest_kmer) ? '+' : '-',
            #             "$(overlap)M"]
            # line = join(link, '\t')
            # println(io, line)
        end
    end
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Load a graph structure from a serialized file.

# Arguments
- `file::AbstractString`: Path to the file containing the serialized graph data

# Returns
- The deserialized graph object
"""
function load_graph(file)
    return FileIO.load(file)["graph"]
end

# function bidirectional_dijkstra(graph, a::T, b::T) where T <: Kmers.Kmer
#     forward_distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     forward_distances[a] = 0

#     forward_arrival_paths = Dict{T, Vector{T}}()
#     forward_arrival_paths[a] = Int[]

#     forward_queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(forward_queue, a, 0)

#     reverse_b = BioSequences.reverse_complement(b)
#     reverse_distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     reverse_distances[reverse_b] = 0

#     reverse_arrival_paths = Dict{T, Vector{T}}()
#     reverse_arrival_paths[reverse_b] = Int[]

#     reverse_queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(reverse_queue, reverse_b, 0)
    
#     hits = Set{T}()
    
#     while isempty(hits)
#         dijkstra_step!(graph, forward_distances, forward_arrival_paths, forward_queue)
#         dijkstra_step!(graph, reverse_distances, reverse_arrival_paths, reverse_queue)
#         hits = intersect(
#                     keys(forward_arrival_paths), 
#                     BioSequences.reverse_complement.(keys(reverse_arrival_paths)))
#     end
#     lowest_cost, optimal_path = evaluate_hits(hits, forward_distances, reverse_distances, forward_arrival_paths, reverse_arrival_paths)
#     return lowest_cost, optimal_path
# end

# # # simple point to point
# # function dijkstra(graph, a::T, b::T) where {T <: Kmers.Kmer}
# #     distances = DataStructures.DefaultDict{T, Float64}(Inf)
# #     distances[a] = 0

# #     arrival_paths = Dict{T, Vector{T}}()
# #     arrival_paths[a] = Int[]

# #     queue = DataStructures.PriorityQueue{T, Float64}()
# #     DataStructures.enqueue!(queue, a, 0)

# #     current_kmer, cost = a, Inf
# #     while current_kmer != b
# #         distances, arrival_paths, queue = dijkstra_step!(graph, distances, arrival_paths, queue)
# #     end
# #     shortest_path = vcat(arrival_paths[b], b)
# #     return shortest_path, distances[b]    
# # end

# function dijkstra(graph, a::T, targets::Set{T}; search_strategy::Union{Symbol, Missing}=missing) where {T <: Kmers.Kmer}
#     if ismissing(search_strategy)
#         # local breadth first search to find the single target
#         # a-b shortest path
#         if length(targets) == 1
#             search_strategy = :BFS
#         # a - no targets - go find the ends
#         # a - series of any target
#         else
#             search_strategy = :DFS
#         end
#     end 
#     if search_strategy == :BFS
#         default_cost = 1
#     elseif search_strategy == :DFS
#         default_cost = 0
#     else
#         error("please specify either :BFS (breadth-first search) or :DFS (depth-first search)")
#     end
    
#     distances = DataStructures.DefaultDict{T, Float64}(Inf)
#     distances[a] = 0

#     arrival_paths = Dict{T, Vector{T}}()
#     arrival_paths[a] = Int[]

#     queue = DataStructures.PriorityQueue{T, Float64}()
#     DataStructures.enqueue!(queue, a, 0)
    
#     while !isempty(queue) && !(first(DataStructures.peek(queue)) in targets)
#         distances, arrival_paths, queue = 
#             dijkstra_step!(graph, distances, arrival_paths, queue, default_cost = default_cost)
#     end
#     if !isempty(queue)
#         discovered_destination = DataStructures.dequeue!(queue)
#         @assert discovered_destination in targets
#         shortest_path = vcat(arrival_paths[discovered_destination], discovered_destination)
#         return shortest_path, distances[discovered_destination]
#     else
#         longest_path = Int[]
#         distance = Inf
#         for (destination, arrival_path) in arrival_paths
#             this_distance = distances[destination]
#             this_path = vcat(arrival_path, destination)
#             if (length(this_path) > length(longest_path)) && (this_distance <= distance)
#                 longest_path = this_path
#                 distance = this_distance
#             end
            
#         end
#         return longest_path, distance
#     end
# end

# function update_remaining_targets(current_walk::AbstractVector{T}, remaining_targets::AbstractSet{T}) where T <: Kmers.Kmer
#     # assess whether targets have been hit in the canonical space
#     remaining_targets = setdiff(BioSequences.canonical.(remaining_targets), BioSequences.canonical.(current_walk))
#     # blow back out into forward and reverse_complement space
#     remaining_targets = Set{T}(vcat(remaining_targets, BioSequences.reverse_complement.(remaining_targets)))
#     return remaining_targets
# end

# function assess_downstream_weight(graph, kmer)
#     # here we look to see if walking forward or backward from the initial node gets us to heavier weight options
#     score = 0
#     for neighbor in BioSequences.neighbors(kmer)
#         try
#             score += MetaGraphs.get_prop(graph, graph[BioSequences.canonical(neighbor), :kmer], :count)
#         catch
#             continue
#         end
#     end
#     return score
# end

# # vertices should either be entire graph (by default) or a connected component
# # if people want to work on just the connected component, let them induce a subgraph
# function find_graph_core(graph; seed=rand(Int))
    
#     Random.seed!(seed)
    
#     T = typeof(MetaGraphs.get_prop(graph, 1, :kmer))
    
#     # targets = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph)]
#     # sortperm
    
#     # targets = [
#     #         MetaGraphs.get_prop(graph, v_index, :count) for v_index in 
#     #         sortperm([MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)], rev=true)[1:Int(ceil(Graphs.nv(graph)/10))]]
    
#     targets = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph) if MetaGraphs.get_prop(graph, v, :count) > 1]
#             # sortperm([MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)], rev=true)[1:Int(ceil(Graphs.nv(graph)/10))]]
    
#     @show length(targets)
#     starting_kmer = first(targets)
#     max_degree = 0
#     for node in targets
#         node_degree = Graphs.degree(graph, graph[node, :kmer])
#         if node_degree > max_degree
#             max_degree = node_degree
#             starting_kmer = node
#         end
#     end
        
#     current_walk = [starting_kmer]
#     prior_walk_length = length(current_walk)
#     remaining_targets = update_remaining_targets(current_walk, Set(targets))
#     done = isempty(remaining_targets)
    
#     while !done
#         # here we look to see if walking forward or backward from the current ends gets us to heavier weight options
#         # we want to prioritize walks toward higher coverage nodes
#         forward_score = assess_downstream_weight(graph, last(current_walk))
#         reverse_score = assess_downstream_weight(graph, BioSequences.reverse_complement(first(current_walk)))
#         if reverse_score > forward_score
#             current_walk = reverse(BioSequences.reverse_complement.(current_walk))
#         end
        
#         forward_source = last(current_walk)
#         forward_walk, forward_distance = dijkstra(graph, forward_source, remaining_targets, search_strategy=:DFS)
#         current_walk = vcat(current_walk, forward_walk[2:end])
#         remaining_targets = update_remaining_targets(current_walk, remaining_targets)
#         if isempty(remaining_targets)
#             done = true
#         else
#             reverse_source = BioSequences.reverse_complement(first(current_walk))
#             reverse_walk, reverse_distance = dijkstra(graph, reverse_source, remaining_targets, search_strategy=:DFS)
#             current_walk = vcat(reverse(BioSequences.reverse_complement.(reverse_walk))[1:end-1], current_walk)
#             remaining_targets = update_remaining_targets(current_walk, remaining_targets)
#         end
#         @show length(current_walk)
#         failed_this_expansion = length(current_walk) == prior_walk_length
#         prior_walk_length = length(current_walk)
#         if isempty(remaining_targets)
#             done = true
#         elseif failed_this_expansion
#             done = true
#         end
#     end
#     return current_walk
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_edge_probabilities(graph, strand)
#     kmers = graph.gprops[:kmers]
#     outgoing_edge_probabilities = SparseArrays.spzeros(length(kmers), length(kmers))
    
#     for (kmer_index, kmer) in enumerate(kmers)
#         if !strand
#             kmer = BioSequences.reverse_complement(kmer)
#         end
        
#         downstream_neighbor_indices = Int[]
#         for neighbor in BioSequences.neighbors(kmer)
#             index = get_kmer_index(kmers, BioSequences.canonical(neighbor))
#             # kmer must be in our dataset and there must be a connecting edge
#             if !isnothing(index) && Graphs.has_edge(graph, ordered_edge(kmer_index, index))
#                 push!(downstream_neighbor_indices, index)
#             end
#         end
#         sort!(unique!(downstream_neighbor_indices))
        
#         downstream_edge_weights = Int[
#             length(get(graph.edge_evidence, ordered_edge(kmer_index, neighbor_index), EdgeEvidence[])) for neighbor_index in downstream_neighbor_indices
#         ]
        
#         non_zero_indices = downstream_edge_weights .> 0
#         downstream_neighbor_indices = downstream_neighbor_indices[non_zero_indices]
#         downstream_edge_weights = downstream_edge_weights[non_zero_indices]
        
#         downstream_edge_likelihoods = downstream_edge_weights ./ sum(downstream_edge_weights)
        
#         for (neighbor_index, likelihood) in zip(downstream_neighbor_indices, downstream_edge_likelihoods)
#             outgoing_edge_probabilities[kmer_index, neighbor_index] = likelihood
#         end
#     end
#     return outgoing_edge_probabilities
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_edge_probabilities(graph)
#     outgoing_edge_probabilities = determine_edge_probabilities(graph, true)
#     incoming_edge_probabilities = determine_edge_probabilities(graph, false)
#     return outgoing_edge_probabilities, incoming_edge_probabilities
# end





# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function reverse_oriented_path(oriented_path)
# #     reversed_path = copy(oriented_path)
# #     for (index, state) in enumerate(oriented_path)
# #         reversed_path[index] = OrientedKmer(index = state.index, orientation = !state.orientation)
# #     end
# #     return reverse!(reversed_path)
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function maximum_likelihood_walk(graph, connected_component)
# #     max_count = maximum(graph.counts[connected_component])
# #     max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
# #     initial_node_index = rand(max_count_indices)
# #     initial_node = connected_component[initial_node_index]
# #     outgoing_edge_probabilities, incoming_edge_probabilities = determine_edge_probabilities(graph)
# #     forward_walk = maximum_likelihood_walk(graph, [OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
# #     reverse_walk = maximum_likelihood_walk(graph, [OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
# #     reversed_reverse_walk = reverse!(
# #         [
# #             OrientedKmer(index = oriented_kmer.index, orientation = oriented_kmer.orientation)
# #             for oriented_kmer in reverse_walk[2:end]
# #         ]
# #         )
# #     full_path = [reversed_reverse_walk..., forward_walk...]
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function maximum_likelihood_walk(graph, path::Vector{OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
# #     done = false
# #     while !done
# #         maximum_path_likelihood = 0.0
# #         maximum_likelihood_path = Vector{OrientedKmer}()
# #         for neighbor in Graphs.neighbors(graph.graph, last(path).index)
# #             this_path = [last(path).index, neighbor]
# #             this_oriented_path, this_path_likelihood = 
# #                 assess_path(this_path,
# #                     graph.kmers,
# #                     graph.counts,
# #                     last(path).orientation,
# #                     outgoing_edge_probabilities,
# #                     incoming_edge_probabilities)
# #             if this_path_likelihood > maximum_path_likelihood
# #                 maximum_path_likelihood = this_path_likelihood
# #                 maximum_likelihood_path = this_oriented_path
# #             end
# #         end
# #         if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
# #             done = true
# #         else
# #             append!(path, maximum_likelihood_path[2:end])
# #         end
# #     end
# #     return path
# # end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function take_a_walk(graph, connected_component)
# #     max_count = maximum(graph.counts[connected_component])
# #     max_count_indices = findall(count -> count == max_count, graph.counts[connected_component])
# #     initial_node_index = rand(max_count_indices)
# #     initial_node = connected_component[initial_node_index]
# #     outgoing_edge_probabilities, incoming_edge_probabilities = determine_edge_probabilities(graph)
    
# #     # walk forwards from the initial starting node
# #     forward_walk = take_a_walk(graph, [OrientedKmer(index = initial_node, orientation = true)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
# #     # walk backwards from the initial starting node
# #     reverse_walk = take_a_walk(graph, [OrientedKmer(index = initial_node, orientation = false)], outgoing_edge_probabilities, incoming_edge_probabilities)
    
# #     # we need to reverse everything to re-orient against the forward walk
# #     reverse_walk = reverse_oriented_path(reverse_walk)
    
# #     # also need to drop the last node, which is equivalent to the first node of the 
# #     @assert last(reverse_walk) == first(forward_walk)
# #     full_path = [reverse_walk[1:end-1]..., forward_walk...]
# # #     @show full_path
# # end


# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # A short description of the function

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# # function take_a_walk(graph, path::Vector{OrientedKmer}, outgoing_edge_probabilities, incoming_edge_probabilities)
# #     done = false
# #     while !done
# #         maximum_path_likelihood = 0.0
# #         maximum_likelihood_path = Vector{OrientedKmer}()
# #         for neighbor in Graphs.neighbors(graph.graph, last(path).index)
# #             this_path = [last(path).index, neighbor]
# #             this_oriented_path, this_path_likelihood = 
# #                 assess_path(this_path,
# #                     graph.kmers,
# #                     graph.counts,
# #                     last(path).orientation,
# #                     outgoing_edge_probabilities,
# #                     incoming_edge_probabilities)
# #             if this_path_likelihood > maximum_path_likelihood
# #                 maximum_path_likelihood = this_path_likelihood
# #                 maximum_likelihood_path = this_oriented_path
# #             end
# #         end
# #         if isempty(maximum_likelihood_path) && (maximum_path_likelihood == 0.0)
# #             done = true
# #         else
# #             append!(path, maximum_likelihood_path[2:end])
# #         end
# #     end
# #     return path
# # end


# # use default cost 1 for breadth first seach
# # use default cost 0 for depth first search
# function dijkstra_step!(graph, distances, arrival_paths, queue; default_cost = 1)
#     current_kmer, cost = DataStructures.dequeue_pair!(queue)
#     current_orientation = BioSequences.iscanonical(current_kmer)
#     current_canonical_kmer = BioSequences.canonical(current_kmer)
#     current_index = graph[current_canonical_kmer, :kmer]

#     present_neighbors = Vector{typeof(current_kmer)}()
#     for neighbor in BioSequences.neighbors(current_kmer)
#         try
#             graph[BioSequences.canonical(neighbor), :kmer]
#             push!(present_neighbors, neighbor)
#         catch
#             continue
#         end
#     end
    
#     true_neighbors =
#         filter(neighbor -> 
#             Graphs.has_edge(
#                 graph,
#                 graph[BioSequences.canonical(current_kmer), :kmer],
#                 graph[BioSequences.canonical(neighbor), :kmer]), 
#             present_neighbors)
    
#     if !isempty(true_neighbors)
#         canonical_true_neighbors = BioSequences.canonical.(true_neighbors)
#         true_neighbor_is_canonical = BioSequences.iscanonical.(true_neighbors)
#         neighbor_indices = map(canonical_neighbor -> graph[canonical_neighbor, :kmer], canonical_true_neighbors)
#         neighbor_weights = map(neighbor_i -> MetaGraphs.get_prop(graph, current_index, neighbor_i, :weight), neighbor_indices)


#         # inverse probability!!
#         neighbor_costs = 1 .- (neighbor_weights ./ sum(neighbor_weights))
#         neighbor_costs .+= default_cost
#         # afterwards, 100% likelihood edges cost 1
#         # 50% likelihood costs 1.5
#         # 0% likelihood is impossible and not included

#         for (neighbor, cost) in zip(true_neighbors, neighbor_costs)    
#             candidate_path = vcat(arrival_paths[current_kmer], current_kmer)
#             candidate_path_cost = distances[current_kmer] + cost  
#             if distances[neighbor] > candidate_path_cost
#                 arrival_paths[neighbor] = candidate_path
#                 DataStructures.enqueue!(queue, neighbor => cost)
#                 distances[neighbor] = candidate_path_cost
#             end  
#         end
#     end
#     return distances, arrival_paths, queue
# end

# function apply_kmedoids_treshold(graph)
#     kmer_counts = [MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)]

#     kmer_counts_histogram = sort(collect(StatsBase.countmap(values(kmer_counts))), by=x->x[1])

# #     scale = 250
# #     p = plot_kmer_frequency_spectra(values(kmer_counts), size=(2scale,scale), log_scale=log2, title="kmer frequencies")
# #     display(p)

# #     p = StatsPlots.scatter(log2.(first.(kmer_counts_histogram)))
# #     display(p)

#     kmer_depth_of_coverage_bins = log2.(first.(kmer_counts_histogram))

#     distance_matrix = zeros((length(kmer_depth_of_coverage_bins), length(kmer_depth_of_coverage_bins)))
#     for (row, depth_of_coverage_bin_1) in enumerate(kmer_depth_of_coverage_bins)
#         for (col, depth_of_coverage_bin_2) in enumerate(kmer_depth_of_coverage_bins)
#             distance = abs(depth_of_coverage_bin_1 - depth_of_coverage_bin_2)
#             distance_matrix[row, col] = distance
#         end
#     end
#     distance_matrix

#     # max out k at the same max k we use for DNAMers
#     max_k = min(length(kmer_depth_of_coverage_bins), 63)
#     ks = Primes.primes(2, max_k)
#     ys = map(k ->
#                 Statistics.mean(Statistics.mean(Clustering.silhouettes(Clustering.kmedoids(distance_matrix, k), distance_matrix)) for i in 1:100),
#                 ks)

#     p = StatsPlots.plot(ks, ys, label="silhouette score", ylabel = "silhouette score", xlabel = "number of clusters")
#     display(p)

#     ymax, ymax_index = findmax(ys)
#     optimal_k = ks[ymax_index]
#     clusterings = [Clustering.kmedoids(distance_matrix, optimal_k) for i in 1:10]
#     max_value, max_value_index = findmax(clustering -> Statistics.mean(Clustering.silhouettes(clustering, distance_matrix)), clusterings)
#     optimal_clustering = clusterings[max_value_index]
#     # optimal_clustering.assignments
#     min_medoid_value, min_medoid_index = findmin(optimal_clustering.medoids)
#     indices_to_include = map(assignment -> assignment .!= min_medoid_index, optimal_clustering.assignments)
#     # kmer_depth_of_coverage_bins
#     threshold = Int(ceil(2^maximum(kmer_depth_of_coverage_bins[.!indices_to_include]))) + 1

#     scale = 250
#     p = plot_kmer_frequency_spectra(values(kmer_counts), log_scale = log2, size=(2scale,scale), title="kmer frequencies")
#     StatsPlots.vline!(p, log2.([threshold]))
#     display(p)

#     # find all vertices with count > threshold
#     vertices_to_keep = [v for v in Graphs.vertices(graph) if (MetaGraphs.get_prop(graph, v, :count) > threshold)]
#     # induce subgraph
#     induced_subgraph, vertex_map = Graphs.induced_subgraph(graph, vertices_to_keep)

#     # set kmer as indexing prop
#     MetaGraphs.set_indexing_prop!(induced_subgraph, :kmer)
#     return induced_subgraph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function simplify_kmer_graph(kmer_graph)
#     @info "simplifying kmer graph"
#     @info "resolving untigs..."
#     @time untigs = resolve_untigs(kmer_graph)
#     @info "determining untig orientations..."
#     oriented_untigs = determine_oriented_untigs(kmer_graph, untigs)
#     simplified_graph = MetaGraphs.MetaDiGraph(length(oriented_untigs))
#     MetaGraphs.set_prop!(simplified_graph, :k, kmer_graph.gprops[:k])
#     @info "initializing graph node metadata"
#     for (vertex, untig) in enumerate(oriented_untigs)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :sequence, untig.sequence)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :path, untig.path)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :orientations, untig.orientations)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :evidence, untig.evidence)
#     end
    
#     # determine oriented edges of simplified graph
#     simplified_untigs = []
#     @info "creating simplified untigs to resolve connections"
#     for vertex in Graphs.vertices(simplified_graph)
#         in_kmer = simplified_graph.vprops[vertex][:path][1] => simplified_graph.vprops[vertex][:orientations][1]
#         out_kmer = simplified_graph.vprops[vertex][:path][end] => simplified_graph.vprops[vertex][:orientations][end]
#     #     @show vertex, in_kmer, out_kmer
#         push!(simplified_untigs, in_kmer => out_kmer)
#     end

#     @info "resolving connections"
#     ProgressMeter.@showprogress for (ui, u) in enumerate(simplified_untigs)
#         for (vi, v) in enumerate(simplified_untigs)
#     #         + => +
#             source_kmer_index, source_orientation = last(u)
#             destination_kmer_index, destination_orientation = first(v)
#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! + +"

#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         + => -
#             source_kmer_index, source_orientation = last(u)
#             destination_kmer_index, destination_orientation = last(v)
#             destination_orientation = !destination_orientation

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! + -"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         - => +
#             source_kmer_index, source_orientation = first(u)
#             source_orientation = !source_orientation
#             destination_kmer_index, destination_orientation = first(v)

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! - +"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#     #         - => -
#             source_kmer_index, source_orientation = first(u)
#             source_orientation = !source_orientation
#             destination_kmer_index, destination_orientation = last(v)
#             destination_orientation = !destination_orientation

#             edge = Graphs.Edge(source_kmer_index, destination_kmer_index)
#             if Graphs.has_edge(kmer_graph, edge)
#                 source_orientation_matches = (kmer_graph.eprops[edge][:orientations].source_orientation == source_orientation)
#                 destination_orientation_matches = (kmer_graph.eprops[edge][:orientations].destination_orientation == destination_orientation)
#                 if source_orientation_matches && destination_orientation_matches
# #                     @show "right orientation!! - -"
#                     simplified_graph_edge = Graphs.Edge(ui, vi)

#                     Graphs.add_edge!(simplified_graph, simplified_graph_edge)
#                     edge_orientations = (
#                         source_orientation = source_orientation,
#                         destination_orientation = destination_orientation
#                     )
#                     set_metadata!(simplified_graph, simplified_graph_edge, :orientations, edge_orientations)
#                 end
#             end
#         end
#     end
#     return simplified_graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function simplify_kmer_graph(kmer_graph)
#     @info "simplifying kmer graph"
#     @info "resolving untigs..."
#     @time untigs = resolve_untigs(kmer_graph)

#     @info "determining untig orientations..."
#     @time oriented_untigs = determine_oriented_untigs(kmer_graph, untigs)

#     simplified_graph = MetaGraphs.MetaDiGraph(length(oriented_untigs))
#     MetaGraphs.set_prop!(simplified_graph, :k, kmer_graph.gprops[:k])
#     @info "initializing graph node metadata"
#     ProgressMeter.@showprogress for (vertex, untig) in enumerate(oriented_untigs)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :sequence, untig.sequence)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :path, untig.path)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :orientations, untig.orientations)
#         MetaGraphs.set_prop!(simplified_graph, vertex, :weight, untig.weight)
#     end

#     # determine oriented edges of simplified graph
#     simplified_untigs = Vector{Pair{Pair{Int64,Bool},Pair{Int64,Bool}}}(undef, length(Graphs.vertices(simplified_graph)))
#     @info "creating simplified unitgs to help resolve connections"
#     # use a pre-allocated array here to speed up
#     ProgressMeter.@showprogress for vertex in Graphs.vertices(simplified_graph)
#         in_kmer = simplified_graph.vprops[vertex][:path][1] => simplified_graph.vprops[vertex][:orientations][1]
#         out_kmer = simplified_graph.vprops[vertex][:path][end] => simplified_graph.vprops[vertex][:orientations][end]
#     #     @show vertex, in_kmer, out_kmer
#         simplified_untigs[vertex] = in_kmer => out_kmer
#     #     push!(simplified_untigs, )
#     end

#     # make a dictionary mapping endcap to oriented_untig index

#     end_mer_map = Dict()
#     ProgressMeter.@showprogress for (i, oriented_untig) in enumerate(oriented_untigs)
#         end_mer_map[first(oriented_untig.path)] = i
#         end_mer_map[last(oriented_untig.path)] = i
#     end

#     ProgressMeter.@showprogress for (untig_index, oriented_untig) in enumerate(oriented_untigs)
#     #     @show untig_index
#         true_in_overlap = oriented_untig.sequence[1:simplified_graph.gprops[:k]-1]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[1])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[2])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_out_overlap = neighboring_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_true_out_overlap == true_in_overlap
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = true => true
#                 o = (source_orientation = true, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_out_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_false_out_overlap == true_in_overlap        
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = false => true
#                 o = (source_orientation = false, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         true_out_overlap = oriented_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[end])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[end-1])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_in_overlap = neighboring_untig.sequence[1:simplified_graph.gprops[:k]-1]
#             if true_out_overlap == neighbor_true_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = true => true
#                 o = (source_orientation = true, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_in_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[1:simplified_graph.gprops[:k]-1]
#             if true_out_overlap == neighbor_false_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = true => false
#                 o = (source_orientation = true, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         false_in_overlap = BioSequences.reverse_complement(oriented_untig.sequence)[1:simplified_graph.gprops[:k]-1]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[end])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[end-1])
#         end
#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_out_overlap = neighboring_untig.sequence[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_true_out_overlap == false_in_overlap
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = true => false
#                 o = (source_orientation = true, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_out_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]
#             if neighbor_false_out_overlap == false_in_overlap        
#                 e = Graphs.Edge(neighboring_untig_index, untig_index)
#     #             o = false => false
#                 o = (source_orientation = false, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end

#         false_out_overlap = BioSequences.reverse_complement(oriented_untig.sequence)[end-simplified_graph.gprops[:k]+2:end]

#         non_backtracking_neighbors = Graphs.neighbors(kmer_graph, oriented_untig.path[1])
#         if length(oriented_untig.path) > 1
#             non_backtracking_neighbors = setdiff(non_backtracking_neighbors, oriented_untig.path[2])
#         end

#         for non_backtracking_neighbor in non_backtracking_neighbors
#             neighboring_untig_index = end_mer_map[non_backtracking_neighbor]
#             neighboring_untig = oriented_untigs[neighboring_untig_index]

#             neighbor_true_in_overlap = neighboring_untig.sequence[1:simplified_graph.gprops[:k]-1]
#             if false_out_overlap == neighbor_true_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = false => true
#                 o = (source_orientation = false, destination_orientation = true)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end

#             neighbor_false_in_overlap = BioSequences.reverse_complement(neighboring_untig.sequence)[1:simplified_graph.gprops[:k]-1]
#             if false_out_overlap == neighbor_false_in_overlap
#                 e = Graphs.Edge(untig_index, neighboring_untig_index)
#     #             o = false => false
#                 o = (source_orientation = false, destination_orientation = false)
#                 Graphs.add_edge!(simplified_graph, e)
#                 set_metadata!(simplified_graph, e, :orientations, o)    
#             end
#         end
#     end
#     return simplified_graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Description

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function determine_oriented_untigs(kmer_graph, untigs)
#     oriented_untigs = []
#     for path in untigs
#         sequence = BioSequences.LongDNASeq(kmer_graph.vprops[first(path)][:kmer])
#         if length(path) == 1
#             orientations = [true]
#         elseif length(path) > 1
#             initial_edge = Graphs.Edge(path[1], path[2])
#             initial_orientation = first(kmer_graph.eprops[initial_edge][:orientations]).source_orientation
#             orientations = [initial_orientation]
#             if !initial_orientation
#                 sequence = BioSequences.reverse_complement(sequence)
#             end

#             for (src, dst) in zip(path[1:end-1], path[2:end])
#                 edge = Graphs.Edge(src, dst)
#                 destination = BioSequences.LongDNASeq(kmer_graph.vprops[edge.dst][:kmer])
#                 destination_orientation = first(kmer_graph.eprops[edge][:orientations]).destination_orientation
#                 push!(orientations, destination_orientation)
#                 if !destination_orientation
#                     destination = BioSequences.reverse_complement(destination)
#                 end
#                 sequence_suffix = sequence[end-length(destination)+2:end]
#                 destination_prefix = destination[1:end-1]
#                 @assert sequence_suffix == destination_prefix
#                 push!(sequence, destination[end])
#             end
#         end

#         oriented_untig = 
#         (
#             sequence = BioSequences.canonical(sequence),
#             path = BioSequences.iscanonical(sequence) ? path : reverse(path),
#             orientations = BioSequences.iscanonical(sequence) ? orientations : reverse(.!orientations),
#             weight = Statistics.median([kmer_graph.vprops[v][:weight] for v in path])
#         )

#         push!(oriented_untigs, oriented_untig)
#     end
#     return oriented_untigs
# end

# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # Description

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# function resolve_untigs(kmer_graph)
#     untigs = Vector{Int}[]
#     visited = falses(Graphs.nv(kmer_graph))
#     first_unvisited = findfirst(!, visited)
#     while first_unvisited != nothing
#         forward_walk = oriented_unbranching_walk(kmer_graph, first_unvisited, true)
#         reverse_walk = oriented_unbranching_walk(kmer_graph, first_unvisited, false)
#         inverted_reverse_walk = [Graphs.Edge(e.dst, e.src) for e in reverse(reverse_walk)]
#         edges = vcat(inverted_reverse_walk, forward_walk)
#         if isempty(edges)
#             untig = [first_unvisited]
#         else
#             untig = vcat([first(edges).src], [edge.dst for edge in edges])
#         end
#         push!(untigs, untig)
#         for vertex in untig
#             visited[vertex] = true
#         end
#         first_unvisited = findfirst(!, visited)
#     end
#     return untigs
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function oriented_unbranching_walk(kmer_graph, vertex, orientation)
#     walk = []
#     viable_neighbors = find_unbranched_neighbors(kmer_graph, vertex, orientation)
#     while length(viable_neighbors) == 1
# #         @show "found a viable neighbor!!"
#         viable_neighbor = first(viable_neighbors)
#         edge = Graphs.Edge(vertex, viable_neighbor)
#         push!(walk, edge)
#         vertex = edge.dst
#         viable_neighbors = Set{Int}()
#         destination_orientations = [o.destination_orientation for o in kmer_graph.eprops[edge][:orientations]]
#         for destination_orientation in destination_orientations
#             union!(viable_neighbors, find_unbranched_neighbors(kmer_graph, vertex, destination_orientation))
#         end
#     end
#     return walk
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
# or get fasta directly from FTP site

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function find_unbranched_neighbors(kmer_graph, vertex, orientation)
#     downstream_vertices = find_downstream_vertices(kmer_graph, vertex, orientation)
# #     backtrack_vertices
# #     @show downstream_vertices
#     if length(downstream_vertices) == 1
#         downstream_vertex = first(downstream_vertices)
# #         @show downstream_vertex
#         edge = Graphs.Edge(vertex, downstream_vertex)
# #         @show edge
#         destination_orientations = [o.destination_orientation for o in kmer_graph.eprops[edge][:orientations]]
# #         @show destination_orientations
#         for destination_orientation in destination_orientations
#             backtrack_vertices = find_downstream_vertices(kmer_graph, downstream_vertex, !destination_orientation)
# #             @show backtrack_vertices
#             # if the only backtrack is the vertex we're on, then we can simplify
#             if backtrack_vertices == Set([vertex])
#                 return downstream_vertices
#             end
#         end
#     end
#     return Int[]
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
# or get fasta directly from FTP site

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function find_downstream_vertices(kmer_graph, vertex, orientation)
#     viable_neighbors = Set{Int}()
#     for neighbor in Graphs.neighbors(kmer_graph, vertex)
#         not_same_vertex = vertex != neighbor
#         candidate_edge = Graphs.Edge(vertex, neighbor)
#         # palindromes can have multiple viable orientations
#         # check each viable orientation individually
#         edge_src_orientations = [e.source_orientation for e in kmer_graph.eprops[candidate_edge][:orientations]]
#         for edge_src_orientation in edge_src_orientations
# #             edge_src_orientation = kmer_graph.eprops[candidate_edge][:orientations].source_orientation
#             viable_orientation = edge_src_orientation == orientation
#             if not_same_vertex && viable_orientation
#                 push!(viable_neighbors, neighbor)
#             end
#         end
#     end
#     return viable_neighbors
# end

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


################################################################################
# defunct bcalm usage
# run(`conda install -c bioconda bcalm`)

# fasta_file_list = joinpath(working_directory, repr(hash(fastx_files)) * ".fasta_list.txt")
# open(fasta_file_list, "w") do io
#     for f in fastx_files
#         @show f
#         println(io, f)
#     end
# end

# k = 3
# outfile = fasta_file_list * ".bcalm.$(k).fna"
# cmd = `bcalm -in $(fastx_files[1]) -abundance-min 1 -kmer-size 11 -all-abundance-counts -out $(outfile)`
# run(cmd)

# cmds = [
#     `bcalm`,
#     `-in $(fasta_list_file)`,
#     `-abundance-min 1`,
#     `-kmer-size 3`,
#     `-all-abundance-counts`,
#     `-abundance-max $(typemax(UInt64))`
# ]
# run(cmds)

# ls -1 *.fastq > list_reads
# ./bcalm -in list_reads [..]
# ./bcalm -in [reads.fa] -kmer-size [kmer_size] -abundance-min [abundance_threshold]

# scripts/convertToGFA.py
##################################################################################

# # @inline function add_edge_to_simple_kmer_graph!(simple_kmer_graph, kmers, sequence_edge)
# @inline function add_edge_to_simple_kmer_graph!(simple_kmer_graph, sequence_edge)
#     observed_source_kmer, observed_destination_kmer = edgemer_to_vertex_kmers(sequence_edge.fw)
#     canonical_source_kmer = BioSequences.canonical(observed_source_kmer)
#     canonical_destination_kmer = BioSequences.canonical(observed_destination_kmer)
#     source_kmer_index = simple_kmer_graph[canonical_source_kmer, :kmer]
#     desination_kmer_index = simple_kmer_graph[canonical_destination_kmer, :kmer]
    
#     if source_kmer_index > desination_kmer_index
#         observed_source_kmer, observed_destination_kmer = edgemer_to_vertex_kmers(sequence_edge.bw)
#         canonical_source_kmer = BioSequences.canonical(observed_source_kmer)
#         canonical_destination_kmer = BioSequences.canonical(observed_destination_kmer)
#         source_kmer_index = simple_kmer_graph[canonical_source_kmer, :kmer]
#         desination_kmer_index = simple_kmer_graph[canonical_destination_kmer, :kmer]
#     end
    
#     @assert source_kmer_index <= desination_kmer_index

#     oriented_source_kmer = 
#         (canonical_kmer = canonical_source_kmer,
#          orientation = BioSequences.iscanonical(observed_source_kmer))

#     oriented_destination_kmer = 
#         (canonical_kmer = canonical_destination_kmer,
#          orientation = BioSequences.iscanonical(observed_destination_kmer))

#     oriented_source_vertex = 
# #         (vertex = searchsortedfirst(kmers, oriented_source_kmer.canonical_kmer),
#         (vertex = simple_kmer_graph[oriented_source_kmer.canonical_kmer, :kmer],
#          orientation = oriented_source_kmer.orientation)

#     oriented_destination_vertex = 
# #         (vertex = searchsortedfirst(kmers, oriented_destination_kmer.canonical_kmer),
#         (vertex = simple_kmer_graph[oriented_destination_kmer.canonical_kmer, :kmer],
#          orientation = oriented_destination_kmer.orientation)

# #     forward_edge = Graphs.Edge(oriented_source_vertex.vertex, oriented_destination_vertex.vertex)
# #     forward_edge_orientations = 
# #         (source_orientation = oriented_source_vertex.orientation,
# #          destination_orientation = oriented_destination_vertex.orientation)
    
# #     reverse_edge = Graphs.Edge(oriented_destination_vertex.vertex, oriented_source_vertex.vertex)
# #     reverse_edge_orientations = 
# #         (source_orientation = !oriented_destination_vertex.orientation,
# #          destination_orientation = !oriented_source_vertex.orientation)
    
#     edge = Graphs.Edge(oriented_source_vertex.vertex, oriented_destination_vertex.vertex)
#     edge_orientations = 
#         (source_orientation = oriented_source_vertex.orientation,
#          destination_orientation = oriented_destination_vertex.orientation)
    
# #     orientations = Set([forward_edge_orientations, reverse_edge_orientations])
#     orientations = Set([edge_orientations])
#     if Graphs.has_edge(simple_kmer_graph, edge)
#         edge_weight = MetaGraphs.get_prop(simple_kmer_graph, edge, :weight) + 1
#         orientations = union(MetaGraphs.get_prop(simple_kmer_graph, edge, :orientations), orientations)
#     else
#         Graphs.add_edge!(simple_kmer_graph, edge)
#         edge_weight = 1
# #         @show Graphs.ne(simple_kmer_graph)
# #         @show forward_edge
#     end
#     MetaGraphs.set_prop!(simple_kmer_graph, edge, :weight, edge_weight)
#     MetaGraphs.set_prop!(simple_kmer_graph, edge, :orientations, orientations)
# end

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

Determines the contig with the greatest number of total bases mapping to it

Identify the primary contig based on mapping coverage from Qualimap results.

# Arguments
- `qualimap_results::DataFrame`: DataFrame containing Qualimap alignment statistics with 
  columns "Contig" and "Mapped bases"

# Returns
- `String`: Name of the contig with the highest number of mapped bases

# Description
Takes Qualimap alignment results and determines which contig has the most total bases 
mapped to it, which often indicates the main chromosomal assembly.
"""
function determine_primary_contig(qualimap_results)
    primary_contig_index = last(findmax(qualimap_results[!, "Mapped bases"]))
    primary_contig = qualimap_results[primary_contig_index, "Contig"]
    return primary_contig
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Primary contig is defined as the contig with the most bases mapped to it

In the context of picking out phage from metagenomic assemblies
the longest contig is often bacteria whereas the highest coverage contigs are often primer-dimers or other PCR amplification artifacts.

Taking the contig that has the most bases mapped to it as a product of length * depth is cherry picked as our phage

Isolates and exports the primary contig from an assembly based on coverage depth × length.

The primary contig is defined as the contig with the highest total mapped bases 
(coverage depth × length). This method helps identify potential phage contigs in 
metagenomic assemblies, avoiding both long bacterial contigs and short high-coverage 
PCR artifacts.

# Arguments
- `assembled_fasta`: Path to the assembled contigs in FASTA format
- `assembled_gfa`: Path to the assembly graph in GFA format
- `qualimap_report_txt`: Path to Qualimap coverage report
- `identifier`: String identifier for the output file
- `k`: Integer representing k-mer size used in assembly
- `primary_contig_fasta`: Optional output filepath (default: "\${identifier}.primary_contig.fna")

# Returns
- Path to the output FASTA file containing the primary contig

# Notes
- For circular contigs, removes the k-mer closure scar if detected
- Trims k bases from the end if they match the first k bases
- Uses coverage × length to avoid both long bacterial contigs and short PCR artifacts
"""
function isolate_normalized_primary_contig(assembled_fasta, assembled_gfa, qualimap_report_txt, identifier, k::Int; primary_contig_fasta = "$(identifier).primary_contig.fna")
    
    qualimap_results = parse_qualimap_contig_coverage(qualimap_report_txt)
    primary_contig = determine_primary_contig(qualimap_results)

    # Find primary contig from scaffolds, then export as primary_contig.fasta
    for record in FASTX.FASTA.Reader(open(assembled_fasta))
        record_id = FASTX.identifier(record)
        if record_id == primary_contig
            primary_contig_sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)

            # If the primary contig is circular, need to trim to remove closure scar
            if Mycelia.contig_is_circular(assembled_gfa, primary_contig)
                # trim k-length from end before writing if it matches the first k of the contig
                if primary_contig_sequence[1:k] == primary_contig_sequence[end-k+1:end]
                    for i in 1:k pop!(primary_contig_sequence) end
                end
            end

        w = FASTX.FASTA.Writer(open(primary_contig_fasta, "w")) 
            write(w, FASTX.FASTA.Record(identifier, primary_contig_sequence))
            close(w)
        end
    end
    return primary_contig_fasta
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns bool indicating whether the contig is cleanly assembled.

By cleanly assembled we mean that the contig does not have other contigs attached in the same connected component.

graph_file = path to assembly graph.gfa file
contig_name = name of the contig

Check if a contig exists in isolation within its connected component in an assembly graph.

# Arguments
- `graph_file::String`: Path to the assembly graph file in GFA format
- `contig_name::String`: Name/identifier of the contig to check

# Returns
- `Bool`: `true` if the contig exists alone in its connected component, `false` otherwise

# Details
A contig is considered "cleanly assembled" if it appears as a single entry in the 
assembly graph's connected components. This function parses the GFA file and checks
the contig's isolation status using the graph structure.
"""
function contig_is_cleanly_assembled(graph_file::String, contig_name::String)
    contig_table, records = Mycelia.gfa_to_structure_table(graph_file)
    matching_connected_components = findfirst(contig_table[!, "contigs"] .== contig_name)
    # contigs should be comma-seperated identifiers as strings, so if we get a hit then
    # it must be cleanly assembled.
    return !isnothing(matching_connected_components)
    # gfa_graph = Mycelia.parse_gfa(graph_file)
    # if !isempty(gfa_graph.gprops[:paths])
    #     # probably spades
    #     # gfa segment identifiers have _1 appended to fasta sequence identifiers for spades
    #     segment_identifier = contig_name * "_1"
    #     segment_node_string = gfa_graph.gprops[:paths][segment_identifier]["segments"]
    #     nodes_in_segment = replace.(split(segment_node_string, ','), r"[^\d]" => "")
    #     node_indices = [gfa_graph[n, :identifier] for n in nodes_in_segment]
    # else
    #     # megahit
    #     node_indices = [gfa_graph[contig_name, :identifier]]
    # end
    # connected_components = Graphs.connected_components(gfa_graph)
    # component_of_interest = first(filter(cc -> all(n -> n in cc, node_indices), connected_components))
    # if (1 <= length(component_of_interest) <= 2)
    #     return true
    # else
    #     return false
    # end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Returns bool indicating whether the contig is a circle

graph_file = path to assembly graph.gfa file
contig_name = name of the contig

Determine if a contig represents a circular structure in the assembly graph.

A circular contig is one where the sequence forms a complete loop in the assembly graph,
typically representing structures like plasmids, circular chromosomes, or other circular DNA elements.

# Arguments
- `graph_file::String`: Path to the assembly graph in GFA format
- `contig_name::String`: Name/identifier of the contig to check

# Returns
- `Bool`: `true` if the contig forms a circular structure, `false` otherwise
"""
function contig_is_circular(graph_file::String, contig_name::String)
    contig_table, records = Mycelia.gfa_to_structure_table(graph_file)
    # display(contig_name)
    # display(contig_table)
    matching_connected_components = findfirst(contig_table[!, "contigs"] .== contig_name)
    if !isnothing(matching_connected_components)
        return contig_table[matching_connected_components, "is_circular"] && contig_table[matching_connected_components, "is_closed"]
    else
        return false
    end
    # gfa_graph = Mycelia.parse_gfa(graph_file)
    # if !isempty(gfa_graph.gprops[:paths])
    #     # probably spades
    #     # gfa segment identifiers have _1 appended to fasta sequence identifiers for spades
    #     segment_identifier = contig_name * "_1"
    #     segment_node_string = gfa_graph.gprops[:paths][segment_identifier]["segments"]
    #     nodes_in_segment = replace.(split(segment_node_string, ','), r"[^\d]" => "")
    #     node_indices = [gfa_graph[n, :identifier] for n in nodes_in_segment]
    # else
    #     # megahit
    #     node_indices = [gfa_graph[contig_name, :identifier]]
    # end
    # connected_components = Graphs.connected_components(gfa_graph)
    # component_of_interest = first(filter(cc -> all(n -> n in cc, node_indices), connected_components))
    # subgraph, vertex_map = Graphs.induced_subgraph(gfa_graph, component_of_interest)
    # if Graphs.is_cyclic(subgraph)
    #     return true
    # else
    #     return false
    # end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create an in-memory kmer-graph that records:
- all kmers
- counts
- all *observed* edges between kmers
- edge orientations
- edge counts

Construct a kmer-graph from one or more FASTX files (FASTA/FASTQ).

# Arguments
- `KMER_TYPE`: Type for kmer representation (e.g., `DNAKmer{K}`)
- `fastxs`: Vector of paths to FASTX files

# Returns
A `MetaGraph` where:
- Vertices represent unique kmers with properties:
  - `:kmer` => The kmer sequence
  - `:count` => Number of occurrences
- Edges represent observed kmer adjacencies with properties:
  - `:orientation` => Relative orientation of connected kmers
  - `:count` => Number of observed transitions
"""
function fastx_to_kmer_graph(KMER_TYPE, fastxs::AbstractVector{<:AbstractString})
    
    @info "counting kmers"
    @time kmer_counts = count_canonical_kmers(KMER_TYPE, fastxs)
    K = length(keys(kmer_counts))
    k = length(first(keys(kmer_counts)))
    
    @info "initializing graph with $(K) $(k)mer states"
    graph = MetaGraphs.MetaGraph(K)
    MetaGraphs.set_prop!(graph, :k, k)
    
    @info "adding kmers and kmer counts"
    ProgressMeter.@showprogress for (i, (kmer, count)) in enumerate(kmer_counts)
    #     @show i, kmer, count
        MetaGraphs.set_prop!(graph, i, :kmer, kmer)
        MetaGraphs.set_prop!(graph, i, :count, count)
    end
    
    @info "indexing kmers"
    # allow graph[kmer, :kmer] to dict-lookup the index of a kmer
    MetaGraphs.set_indexing_prop!(graph, :kmer)
    
    @info "adding records to graph"
    graph = add_fastx_records_to_graph!(graph, fastxs)

    @info "adding edges and edge counts"
    graph = add_record_edgemers_to_graph!(graph)
    
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Constructs a k-mer graph from a single FASTX format string.

# Arguments
- `KMER_TYPE`: The k-mer type specification (e.g., DNAKmer{K} where K is k-mer length)
- `fastx::AbstractString`: Input sequence in FASTX format (FASTA or FASTQ)

# Returns
- `KmerGraph`: A directed graph where vertices are k-mers and edges represent overlaps
"""
function fastx_to_kmer_graph(KMER_TYPE, fastx::AbstractString)
    return fastx_to_kmer_graph(KMER_TYPE, [fastx])
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add FASTX records from multiple files as a graph property.

# Arguments
- `graph`: A MetaGraph that will store the FASTX records
- `fastxs`: Collection of FASTA/FASTQ file paths to process

# Details
Creates a dictionary mapping sequence descriptions to their corresponding FASTX records,
then stores this dictionary as a graph property under the key `:records`.
Multiple input files are merged, with later files overwriting records with duplicate descriptions.

# Returns
The modified graph with added records property.
"""
function add_fastx_records_to_graph!(graph, fastxs)
    record_dict = Dict(String(FASTX.description(record)) => record for record in Mycelia.open_fastx(first(fastxs)))
    for fastx in fastxs[2:end]
        for record in Mycelia.open_fastx(fastx)
            record_dict[String(FASTX.description(record))] = record
        end
    end
    MetaGraphs.set_prop!(graph, :records, record_dict)
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Processes DNA sequence records stored in the graph and adds their edgemers (k+1 length subsequences) 
to build the graph structure.

# Arguments
- `graph`: A Mycelia graph object containing DNA sequence records and graph properties

# Details
- Uses the k-mer size specified in `graph.gprops[:k]` to generate k+1 length edgemers
- Iterates through each record in `graph.gprops[:records]`
- For each record, generates all possible overlapping edgemers
- Adds each edgemer to the graph with its position and record information

# Returns
- The modified graph object with added edgemer information
"""
function add_record_edgemers_to_graph!(graph)
    edgemer_size = graph.gprops[:k] + 1
    ProgressMeter.@showprogress for (record_id, record) in graph.gprops[:records]
        for (index, observed_edgemer) in Kmers.EveryKmer{Kmers.DNAKmer{edgemer_size}}(BioSequences.LongDNA{2}(FASTX.sequence(record)))
            add_edgemer_to_graph!(graph, record_id, index, observed_edgemer)
        end
    end
    return graph
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert an edgemer to two vertex kmers.

This function takes an edgemer (a sequence of DNA nucleotides) and converts it into two vertex kmers. 
A kmer is a substring of length k from a DNA sequence. The first kmer is created from the first 
n-1 elements of the edgemer, and the second kmer is created from the last n-1 elements of the edgemer.

# Arguments
- `edgemer::AbstractVector{T}`: A vector of DNA nucleotides where T is a subtype of `BioSequences.DNAAlphabet{2}`.

# Returns
- `Tuple{Kmers.Kmer{BioSequences.DNAAlphabet{2}}, Kmers.Kmer{BioSequences.DNAAlphabet{2}}}`: A tuple containing two kmers.
"""
function edgemer_to_vertex_kmers(edgemer)
    a = Kmers.Kmer{BioSequences.DNAAlphabet{2}}(collect(edgemer[i] for i in 1:length(edgemer)-1))
    b = Kmers.Kmer{BioSequences.DNAAlphabet{2}}(collect(edgemer[i] for i in 2:length(edgemer)))
    return a, b
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Add an observed edgemer to a graph with its associated metadata.

# Arguments
- `graph::MetaGraph`: The graph to modify
- `record_identifier`: Identifier for the source record
- `index`: Position where edgemer was observed
- `observed_edgemer`: The biological sequence representing the edgemer

# Details
Processes the edgemer by:
1. Splitting it into source and destination kmers
2. Converting kmers to their canonical forms
3. Creating or updating an edge with orientation metadata
4. Storing observation details (record, position, orientation)

# Returns
Modified graph with the new edge and metadata

# Note
If the edge already exists, the observation is added to the existing metadata.
"""
function add_edgemer_to_graph!(graph, record_identifier, index, observed_edgemer)
    # observed_orientation = BioSequences.iscanonical(observed_edgemer)
    # canonical_edgemer = BioSequences.canonical(observed_edgemer)
    observed_source_kmer, observed_destination_kmer = Mycelia.edgemer_to_vertex_kmers(observed_edgemer)
    
    canonical_source_kmer = BioSequences.canonical(observed_source_kmer)
    canonical_destination_kmer = BioSequences.canonical(observed_destination_kmer)
    
    source_kmer_orientation = BioSequences.iscanonical(observed_source_kmer)
    destination_kmer_orientation = BioSequences.iscanonical(observed_destination_kmer)
    
    source_kmer_index = graph[canonical_source_kmer, :kmer]
    destination_kmer_index = graph[canonical_destination_kmer, :kmer]
    
    edgemer = Graphs.Edge(source_kmer_index, destination_kmer_index)
    orientation = (source_kmer_orientation => destination_kmer_orientation)
    observation = (;record_identifier, index, orientation)
    
    if !Graphs.has_edge(graph, edgemer)
        Graphs.add_edge!(graph, edgemer)
        MetaGraphs.set_prop!(graph, edgemer, :observations, Set([observation]))
    else
        observations = push!(MetaGraphs.get_prop(graph, edgemer, :observations), observation)
        MetaGraphs.set_prop!(graph, edgemer, :observations, observations)
    end
    return graph
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_fastx_to_graph!(graph, fastx_file::AbstractString)
#     for record in Mycelia.open_fastx(fastx_file)
#         add_fastx_record_to_graph!(graph, record)
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_fastx_record_to_graph!(graph, record::FASTX.FASTA.Record)
#     try
#         graph[FASTX.identifier(record), :identifier]
#         @info "node $(FASTX.identifier(record)) already present"
#     catch
#         Graphs.add_vertex!(graph)
#         vertex_id = Graphs.nv(graph)

#         MetaGraphs.set_prop!(graph, vertex_id, :TYPE, typeof(record))

#         MetaGraphs.set_prop!(graph, vertex_id, :identifier, FASTX.identifier(record))

#         MetaGraphs.set_prop!(graph, vertex_id, :description, FASTX.description(record))

#         sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)
#         MetaGraphs.set_prop!(graph, vertex_id, :sequence, sequence)
#     end
#     return graph
# end
    
# function add_fastx_record_to_graph!(graph, record::FASTX.FASTQ.Record)
#     Graphs.add_vertex!(graph)
#     vertex_id = Graphs.nv(graph)
    
#     MetaGraphs.set_prop!(graph, vertex_id, :TYPE, typeof(record))
    
#     MetaGraphs.set_prop!(graph, vertex_id, :identifier, FASTX.identifier(record))
    
#     MetaGraphs.set_prop!(graph, vertex_id, :description, FASTX.description(record))
    
#     sequence = FASTX.sequence(BioSequences.LongDNA{4}, record)
#     MetaGraphs.set_prop!(graph, vertex_id, :sequence, sequence)
    
#     MetaGraphs.set_prop!(graph, vertex_id, :quality, FASTX.quality_scores(record))
    
#     return graph
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# # for kmer_size in kmer_sizes
# function add_fasta_record_kmers_to_graph!(graph, kmer_size)
#     record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
#     for vertex in record_vertices
#         record_identifier = graph.vprops[vertex][:identifier]
#         record_sequence = graph.vprops[vertex][:sequence]
#         kmer_counts = Mycelia.count_canonical_kmers(Kmers.Kmer{BioSequences.DNAAlphabet{4},kmer_size}, record_sequence)
#         Mycelia.add_kmers_to_graph!(graph, keys(kmer_counts))
#         Mycelia.add_record_kmer_counts_to_graph!(graph, kmer_counts, record_identifier)
#         Mycelia.add_record_edgemers_to_graph!(graph, record_identifier, kmer_size)
#     end
#     return graph
# end



# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_record_kmer_counts_to_graph!(graph, kmer_counts, record_identifier)
#     for (kmer, count) in kmer_counts
#         kmer_vertex = graph[kmer, :identifier]
#         record_vertex = graph[record_identifier, :identifier]
#         edge = Graphs.Edge(kmer_vertex, record_vertex)
#         if !Graphs.has_edge(graph, edge)
#             Graphs.add_edge!(graph, edge)
#             MetaGraphs.set_prop!(graph, edge, :TYPE, "RECORD_KMER_COUNT")
#             MetaGraphs.set_prop!(graph, edge, :count, count)
#         else
#             graph_count = MetaGraphs.get_prop(graph, edge, :count)
#             if graph_count != count
#                 @warn "edge found but this count $(count) != current count $(graph_count)"
#             # else
#                 # @info "edge exists and matches current data"
#             end
#         end
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_kmers_to_graph!(graph, kmers)
#     for kmer in kmers
#         if !haskey(graph.metaindex[:identifier], kmer)
#             Graphs.add_vertex!(graph)
#             v = Graphs.nv(graph)
#             MetaGraphs.set_prop!(graph, v, :identifier, kmer)
#             MetaGraphs.set_prop!(graph, v, :sequence, kmer)
#             MetaGraphs.set_prop!(graph, v, :TYPE, typeof(kmer))
#         end
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_metadata_from_table!(
#         graph::MetaGraphs.AbstractMetaGraph,
#         table::DataFrames.AbstractDataFrame;
#         identifier_column::Union{Symbol, AbstractString} = :identifier)
#     for row in DataFrames.eachrow(table)
#         add_metadata_from_table_row!(graph, row, identifier_column)
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_metadata_from_table_row!(graph, row, identifier_column)
#     other_columns = filter(n -> n != :identifier, Symbol.(names(row)))
#     row_metadata_dict = Dict(column => row[column] for column in other_columns)
#     metadata_dict = Dict(row[identifier_column] => row_metadata_dict)
#     Mycelia.add_metadata_to_graph!(graph, metadata_dict::AbstractDict)
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_key_value_pair_to_node!(graph, identifier, key, value)
#     node_id = graph[identifier, :identifier]
#     MetaGraphs.set_prop!(graph, node_id, Symbol(key), value)
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_metadata_to_node!(graph, identifier, metadata::AbstractVector{<:Pair})
#     for (key, value) in metadata
#         add_key_value_pair_to_node!(graph, identifier, key, value)
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function add_metadata_to_graph!(graph, metadata_dict::AbstractDict)
#     for (identifier, metadata) in metadata_dict
#         add_metadata_to_node!(graph, identifier, collect(metadata))
#     end
#     return graph
# end

# need to get identifier column and then all non-identifier columns and then pass that to above
# function add_metadata_to_graph!(graph, metadata_table::DataFrames.AbstractDataFrame)
#     for (identifier, metadata) in metadata_dict
#         add_metadata_to_node!(graph, identifier, collect(metadata))
#     end
#     return graph
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function initialize_graph()
#     graph = MetaGraphs.MetaDiGraph()
#     MetaGraphs.set_indexing_prop!(graph, :identifier)
#     return graph
# end


# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function construct(KMER_TYPE, fastx, out)
#     mkpath(dirname(out))
#     if !occursin(r"\.jld2$", out)
#         out *= ".jld2"
#     end
#     if !isfile(out)
#         graph = fastx_to_kmer_graph(KMER_TYPE, fastx)
#         @info "saving graph"
#         FileIO.save(out, Dict("graph" => graph))
#         return graph
#     else
#         @info "graph $out already exists, loading existing"
#         return load_graph(out)
#     end
# end