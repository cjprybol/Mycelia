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