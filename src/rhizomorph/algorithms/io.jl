# Graph I/O Functions
#
# This file contains functions for reading and writing graphs in various formats.
#
# Supported formats:
# - GFA v1.0 (Graphical Fragment Assembly format)
# - JLD2 (Julia serialization for lossless graph storage)
#
# Based on legacy graph utilities now ported to Rhizomorph

# ============================================================================
# GFA Export
# ============================================================================

"""
    write_gfa_next(graph::MetaGraphsNext.MetaGraph, outfile::AbstractString; segment_naming::Symbol=:numeric)

Export a MetaGraphsNext graph to GFA v1.0 format.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: The graph to export
- `outfile::AbstractString`: Path to output GFA file
- `segment_naming::Symbol`: How to assign GFA segment IDs.
  - `:numeric` (default): Sequential integers (1, 2, 3, ...). Standard GFA convention.
  - `:sequence`: Use the sequence text as the segment ID. Useful for Bandage
    visualization with `--names`, which displays the segment ID as the node label.

# Returns
- `String`: Path to the output file

# GFA Format
Writes:
- Header (H) line with version and creator information
- Segment (S) lines for each vertex with sequence and depth
- Link (L) lines for each edge with strand orientations and overlap

# Details
- Depth (DP:f:) is calculated from coverage data
- Overlap is inferred from k-mer length (k-1) or edge data
- Strand orientations from edge data (+ for Forward, - for Reverse)

# Example
```julia
graph = Mycelia.Rhizomorph.build_kmer_graph(records, 31; dataset_id="dataset_01", mode=:singlestrand)
Mycelia.Rhizomorph.write_gfa_next(graph, "assembly.gfa")

# For Bandage visualization with node labels
Mycelia.Rhizomorph.write_gfa_next(graph, "assembly_labeled.gfa"; segment_naming=:sequence)
```
"""
function write_gfa_next(graph::MetaGraphsNext.MetaGraph, outfile::AbstractString;
        segment_naming::Symbol = :numeric)
    segment_naming in (:numeric, :sequence) ||
        throw(ArgumentError("segment_naming must be :numeric or :sequence, got :$segment_naming"))

    open(outfile, "w") do io
        println(io, "H\tVN:Z:1.0\tMY:Z:Mycelia-Rhizomorph")

        vertex_id_map = Dict()
        for (i, label) in enumerate(MetaGraphsNext.labels(graph))
            vertex_data = graph[label]

            # Calculate depth from evidence
            depth = if hasfield(typeof(vertex_data), :evidence)
                count_evidence(vertex_data)
            elseif hasfield(typeof(vertex_data), :coverage)
                length(vertex_data.coverage)
            else
                1
            end

            # Extract sequence from vertex data
            sequence = if label isa Kmers.Kmer || label isa BioSequences.BioSequence ||
                          label isa AbstractString
                string(label)
            elseif hasfield(typeof(vertex_data), :Kmer)
                string(vertex_data.Kmer)
            elseif hasfield(typeof(vertex_data), :canonical_kmer)
                string(vertex_data.canonical_kmer)
            elseif hasfield(typeof(vertex_data), :canonical_qualmer)
                string(vertex_data.canonical_qualmer.kmer)
            elseif hasfield(typeof(vertex_data), :sequence)
                string(vertex_data.sequence)
            elseif hasfield(typeof(vertex_data), :string_value)
                vertex_data.string_value
            else
                string(label)  # Fallback: use label as sequence
            end

            seg_id = segment_naming == :sequence ? sequence : string(i)
            vertex_id_map[label] = seg_id

            println(io, "S\t$(seg_id)\t$(sequence)\tDP:f:$depth")
        end

        # Determine overlap length
        if !isempty(MetaGraphsNext.labels(graph))
            first_label = first(MetaGraphsNext.labels(graph))
            first_vertex_data = graph[first_label]
            overlap = if hasfield(typeof(first_vertex_data), :Kmer)
                length(first_vertex_data.Kmer) - 1
            elseif hasfield(typeof(first_vertex_data), :canonical_kmer)
                length(first_vertex_data.canonical_kmer) - 1
            elseif hasfield(typeof(first_vertex_data), :canonical_qualmer)
                length(first_vertex_data.canonical_qualmer.kmer) - 1
            else
                0  # Variable-length graphs need edge-specific overlap
            end
        else
            overlap = 0
        end

        # Write edges
        for edge_labels in MetaGraphsNext.edge_labels(graph)
            if length(edge_labels) == 2
                src_label, dst_label = edge_labels
                edge_data = graph[src_label, dst_label]
                src_id = vertex_id_map[src_label]
                dst_id = vertex_id_map[dst_label]

                # Get strand orientations
                src_orientation = if hasfield(typeof(edge_data), :src_strand)
                    edge_data.src_strand == Forward ? '+' : '-'
                else
                    '+'  # Default forward
                end

                dst_orientation = if hasfield(typeof(edge_data), :dst_strand)
                    edge_data.dst_strand == Forward ? '+' : '-'
                else
                    '+'  # Default forward
                end

                # Use edge-specific overlap if available
                edge_overlap = if hasfield(typeof(edge_data), :overlap_length)
                    edge_data.overlap_length
                else
                    overlap
                end

                println(io, "L\t$src_id\t$src_orientation\t$dst_id\t$dst_orientation\t$(edge_overlap)M")
            end
        end
    end

    return outfile
end

"""
    write_gfa_canonical(graph::MetaGraphsNext.MetaGraph, outfile::AbstractString; segment_naming::Symbol=:numeric)

Export a canonical (undirected) k-mer graph to GFA v1.0 with link orientations
derived from sequence overlap.

Intended for graphs produced by `convert_to_canonical` (for example via
`build_kmer_graph_canonical` or `build_qualmer_graph_canonical`). In such a
graph each vertex is the lexicographically smaller of a k-mer and its reverse
complement, and edges are undirected and carry no strand fields, so
`write_gfa_next` writes every link as `+ +`. That is only correct when both
k-mers were observed in their canonical orientation; an adjacency such as
`A -> reverse_complement(B)` must instead be written `L A + B - (k-1)M`.

# Orientation semantics
For each undirected edge `{u, v}`, all four orientation combinations
`(u±, v±)` are tested. Orientation `+` means the canonical k-mer as written in
the S line; `-` means its reverse complement. A link `L u ou v ov (k-1)M` is
written for every combination where the last `k-1` bases of the oriented `u`
equal the first `k-1` bases of the oriented `v`, so every emitted link is
literally spelled by the segment sequences. GFA-equivalent links are written
once: `L u + v -` is the same link as `L v + u -`, and `L u + v +` is the same
link as `L v - u -`. Self-loops (`u == v`) are handled the same way, e.g. a
k-mer whose `(k-1)`-suffix equals the `(k-1)`-prefix of its own reverse
complement yields `L u + u -`.

Orientation is derived from sequence rather than from edge evidence strand
flags, so it is correct regardless of how evidence was recorded.

# Segment lines
S lines are identical to those written by `write_gfa_next` with the same
`segment_naming` (numeric IDs follow enumeration order of
`MetaGraphsNext.labels(graph)`), so tables keyed by segment ID join across
both writers.

# Known limitation
`convert_to_canonical` collapses edges between the same canonical pair that
differ only in orientation into a single undirected edge. The links written
here therefore reflect every valid `(k-1)` overlap between adjacent canonical
k-mers, not the multiplicity or the specific orientations that were observed.
An edge between two canonical k-mers can, in principle, yield a link for an
orientation combination that no read traversed.

# Arguments
- `graph::MetaGraphsNext.MetaGraph`: Undirected graph with nucleotide k-mer labels
  (`Kmers.DNAKmer` or `Kmers.RNAKmer`)
- `outfile::AbstractString`: Path to output GFA file
- `segment_naming::Symbol`: `:numeric` (default) or `:sequence`, as in `write_gfa_next`

# Returns
- `String`: Path to the output file

# Throws
- `ArgumentError` if the graph is directed, `segment_naming` is invalid, or
  labels are not nucleotide k-mers.
- `ErrorException` if an edge has no valid `(k-1)` overlap in any orientation,
  which indicates the graph is not a k-mer adjacency graph.

# Example
```julia
graph = Mycelia.Rhizomorph.build_kmer_graph_canonical(records, 31)
Mycelia.Rhizomorph.write_gfa_canonical(graph, "assembly_canonical.gfa")
```
"""
function write_gfa_canonical(graph::MetaGraphsNext.MetaGraph, outfile::AbstractString;
        segment_naming::Symbol = :numeric)
    segment_naming in (:numeric, :sequence) ||
        throw(ArgumentError("segment_naming must be :numeric or :sequence, got :$segment_naming"))
    Graphs.is_directed(graph) &&
        throw(ArgumentError("write_gfa_canonical expects an undirected canonical graph; use write_gfa_next for directed graphs"))

    labels = collect(MetaGraphsNext.labels(graph))
    for label in labels
        if !(label isa Kmers.DNAKmer || label isa Kmers.RNAKmer)
            throw(ArgumentError("write_gfa_canonical requires nucleotide k-mer labels, got $(typeof(label))"))
        end
    end

    open(outfile, "w") do io
        println(io, "H\tVN:Z:1.0\tMY:Z:Mycelia-Rhizomorph")

        # Segment lines: same content and ID assignment as write_gfa_next.
        vertex_id_map = Dict{eltype(labels), String}()
        vertex_index = Dict{eltype(labels), Int}()
        for (i, label) in enumerate(labels)
            vertex_data = graph[label]
            depth = if hasfield(typeof(vertex_data), :evidence)
                count_evidence(vertex_data)
            elseif hasfield(typeof(vertex_data), :coverage)
                length(vertex_data.coverage)
            else
                1
            end
            sequence = string(label)
            seg_id = segment_naming == :sequence ? sequence : string(i)
            vertex_id_map[label] = seg_id
            vertex_index[label] = i
            println(io, "S\t$(seg_id)\t$(sequence)\tDP:f:$depth")
        end

        # Link lines: one per valid, GFA-distinct orientation combination.
        written = Set{Tuple{Int, Bool, Int, Bool}}()
        for edge_labels in MetaGraphsNext.edge_labels(graph)
            length(edge_labels) == 2 || continue
            u, v = edge_labels
            k = length(u)
            overlap = k - 1
            u_seqs = (string(u), string(BioSequences.reverse_complement(u)))
            v_seqs = (string(v), string(BioSequences.reverse_complement(v)))
            n_valid = 0
            for (u_fwd, u_seq) in zip((true, false), u_seqs),
                (v_fwd, v_seq) in zip((true, false), v_seqs)

                u_seq[(end - overlap + 1):end] == v_seq[1:overlap] || continue
                n_valid += 1
                key = _canonical_gfa_link_key(
                    vertex_index[u], u_fwd, vertex_index[v], v_fwd)
                key in written && continue
                push!(written, key)
                (a, a_fwd, b, b_fwd) = key
                a_id = vertex_id_map[labels[a]]
                b_id = vertex_id_map[labels[b]]
                println(io,
                    "L\t$a_id\t$(a_fwd ? '+' : '-')\t$b_id\t$(b_fwd ? '+' : '-')\t$(overlap)M")
            end
            n_valid > 0 ||
                error("Edge $(u) -- $(v) has no valid $(overlap)-base overlap in any orientation")
        end
    end

    return outfile
end

"""
    _canonical_gfa_link_key(a, a_fwd, b, b_fwd)

Return a representative of the GFA link equivalence class
`(a, a_fwd, b, b_fwd) ≡ (b, !b_fwd, a, !a_fwd)`, so that the same link is
written once regardless of the direction in which it was discovered.

The representative with more `+` orientations is preferred (so `+ +` is written
rather than the equivalent `- -`); ties are broken by segment index order.
"""
function _canonical_gfa_link_key(a::Int, a_fwd::Bool, b::Int, b_fwd::Bool)
    forward = (a, a_fwd, b, b_fwd)
    flipped = (b, !b_fwd, a, !a_fwd)
    n_plus(key) = key[2] + key[4]
    n_plus(forward) == n_plus(flipped) && return min(forward, flipped)
    return n_plus(forward) > n_plus(flipped) ? forward : flipped
end

# ============================================================================
# GFA Import
# ============================================================================

"""
    read_gfa_next(gfa_file::AbstractString, kmer_type::Type, graph_mode::GraphMode=DoubleStrand)

Read a GFA file and create a k-mer graph with the specified k-mer type.

# Arguments
- `gfa_file::AbstractString`: Path to input GFA file
- `kmer_type::Type`: K-mer type to use (e.g., `Kmers.DNAKmer{31}`)
- `graph_mode::GraphMode`: SingleStrand or DoubleStrand (default: DoubleStrand)

# Returns
- `MetaGraphsNext.MetaGraph`: K-mer graph constructed from GFA

# GFA Format Support
Parses GFA v1.0:
- Header (H) lines: Metadata (currently ignored)
- Segment (S) lines: Vertices with sequence data
- Link (L) lines: Edges with strand orientations
- Path (P) lines: Path information (stored for future use)

# Example
```julia
# Read as DNA k-mer graph in double-strand mode
graph = Mycelia.Rhizomorph.read_gfa_next("assembly.gfa", Kmers.DNAKmer{31})

# Read as single-strand
graph = Mycelia.Rhizomorph.read_gfa_next(
    "assembly.gfa",
    Kmers.DNAKmer{31},
    Mycelia.Rhizomorph.SingleStrand,
)
```
"""
function read_gfa_next(gfa_file::AbstractString, kmer_type::Type, graph_mode::GraphMode = DoubleStrand)
    lines = readlines(gfa_file)

    function parse_overlap_str(s::AbstractString)
        m = match(Regex("\\d+"), s)
        isnothing(m) && return 0
        try
            return parse(Int, m.captures[1])
        catch
            return 0
        end
    end

    # Parse GFA file content
    segments = Dict{String, String}()  # id -> sequence
    links = Vector{Tuple{String, Bool, String, Bool, String}}()  # (src_id, src_forward, dst_id, dst_forward, overlap_str)
    paths = Dict{String, Vector{String}}()  # path_name -> vertex_ids

    for line in lines
        fields = split(line, '\t')
        if isempty(fields)
            continue
        end

        line_type = first(fields)

        if line_type == "H"
            continue
        elseif line_type == "S"
            if length(fields) >= 3
                seg_id = fields[2]
                sequence = fields[3]
                segments[seg_id] = sequence
            end
        elseif line_type == "L"
            # Link line: L<tab>src<tab>src_orient<tab>dst<tab>dst_orient<tab>overlap
            if length(fields) >= 6
                src_id = fields[2]
                src_orient = fields[3] == "+"
                dst_id = fields[4]
                dst_orient = fields[5] == "+"
                overlap = fields[6]
                push!(links, (src_id, src_orient, dst_id, dst_orient, overlap))
            end
        elseif line_type == "P"
            if length(fields) >= 3
                path_name = fields[2]
                path_vertices = split(replace(fields[3], r"[+-]" => ""), ',')
                paths[path_name] = string.(path_vertices)
            end
        elseif line_type == "A"
            continue
        else
            @warn "Unknown GFA line type: $line_type in line: $line"
        end
    end

    if kmer_type == String
        segment_lengths = Set(length(sequence) for sequence in values(segments))
        all_zero_link_overlaps = all(parse_overlap_str(overlap_str) == 0 for (_, _, _, _, overlap_str) in links)
        default_string_overlap = if length(segment_lengths) == 1 && all_zero_link_overlaps
            max(first(segment_lengths) - 1, 0)
        else
            0
        end

        graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph();
            label_type = String,
            vertex_data_type = StringVertexData,
            edge_data_type = StringEdgeData,
            weight_function = compute_edge_weight
        )

        id_to_string = Dict{String, String}()
        for (seg_id, sequence) in segments
            graph[sequence] = StringVertexData(sequence)
            id_to_string[seg_id] = sequence
        end

        for (src_id, _, dst_id, _, overlap_str) in links
            if haskey(id_to_string, src_id) && haskey(id_to_string, dst_id)
                src_string = id_to_string[src_id]
                dst_string = id_to_string[dst_id]
                overlap = parse_overlap_str(overlap_str)
                if overlap == 0 && default_string_overlap > 0
                    overlap = default_string_overlap
                end
                graph[src_string, dst_string] = StringEdgeData(overlap)
            else
                @warn "Link references unknown segment: $src_id -> $dst_id"
            end
        end

        return graph
    end

    # Determine vertex data type from kmer_type
    if kmer_type <: Kmers.DNAKmer || kmer_type <: Kmers.RNAKmer || kmer_type <: Kmers.AAKmer
        sample_kmer_type = if isempty(segments)
            kmer_type
        else
            typeof(kmer_type(first(values(segments))))
        end
        vertex_data_type = KmerVertexData{sample_kmer_type}
        edge_data_type = KmerEdgeData
    else
        error("Unsupported k-mer type: $kmer_type")
    end

    # Create MetaGraphsNext graph
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = sample_kmer_type,
        vertex_data_type = vertex_data_type,
        edge_data_type = edge_data_type,
        weight_function = compute_edge_weight
    )

    # Add vertices (segments) - convert sequences to k-mer types without canonicalizing
    for (seg_id, sequence) in segments
        kmer_seq = kmer_type(sequence)
        graph[kmer_seq] = vertex_data_type(kmer_seq)
    end

    # Create mapping from segment IDs to observed k-mers (strand as provided in the GFA)
    id_to_kmer = Dict{String, sample_kmer_type}()
    for (seg_id, sequence) in segments
        kmer_seq = kmer_type(sequence)
        id_to_kmer[seg_id] = kmer_seq
    end

    # Add edges (links)
    # <<<<<<< HEAD
    #     for (src_id, src_forward, dst_id, dst_forward) in links
    #         if haskey(id_to_kmer, src_id) && haskey(id_to_kmer, dst_id)
    #             src_kmer = id_to_kmer[src_id]
    #             dst_kmer = id_to_kmer[dst_id]

    #             # Convert GFA orientations to StrandOrientation
    #             src_strand = src_forward ? Forward : Reverse
    #             dst_strand = dst_forward ? Forward : Reverse

    #             # Create strand-aware edge
    #             graph[src_kmer, dst_kmer] = edge_data_type(Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}())
    # =======
    for (src_id, src_forward, dst_id, dst_forward, _) in links
        if haskey(id_to_kmer, src_id) && haskey(id_to_kmer, dst_id)
            src_kmer = id_to_kmer[src_id]
            dst_kmer = id_to_kmer[dst_id]
            graph[src_kmer, dst_kmer] = edge_data_type()
        else
            @warn "Link references unknown segment: $src_id -> $dst_id"
        end
    end

    return graph
end

"""
    read_gfa_next(gfa_file::AbstractString, graph_mode::GraphMode=DoubleStrand; force_biosequence_graph::Bool=false)

Read a GFA file and auto-detect whether to create a k-mer graph or BioSequence graph.

# Arguments
- `gfa_file::AbstractString`: Path to input GFA file
- `graph_mode::GraphMode`: SingleStrand or DoubleStrand (default: DoubleStrand)
- `force_biosequence_graph::Bool`: Force creation of variable-length BioSequence graph (default: false)

# Returns
- `MetaGraphsNext.MetaGraph`: Graph with auto-detected vertex type

# Auto-Detection Logic
- **Fixed-length**: If all segments are the same length k, creates `DNAKmer{k}` graph
- **Variable-length**: If segments have different lengths, creates BioSequence graph
- **Override**: Use `force_biosequence_graph=true` to force variable-length graph

# Example
```julia
# Auto-detect graph type
graph = Mycelia.Rhizomorph.read_gfa_next("assembly.gfa")

# Force BioSequence graph
graph = Mycelia.Rhizomorph.read_gfa_next("assembly.gfa", force_biosequence_graph=true)
```
"""
function read_gfa_next(gfa_file::AbstractString, graph_mode::GraphMode = DoubleStrand;
        force_biosequence_graph::Bool = false)
    function parse_overlap_str(s::AbstractString)
        m = match(Regex("\\d+"), s)
        isnothing(m) && return 0
        try
            return parse(Int, m.captures[1])
        catch
            return 0
        end
    end

    # Parse GFA to detect segment lengths
    segment_lengths = Set{Int}()
    segments = Dict{String, String}()

    for line in eachline(gfa_file)
        fields = split(line, '\t')
        if isempty(fields) || first(fields) != "S"
            continue
        end

        if length(fields) >= 3
            seg_id = fields[2]
            sequence = fields[3]
            segments[seg_id] = sequence
            push!(segment_lengths, length(sequence))
        end
    end

    # Detect if fixed-length or variable-length
    if !force_biosequence_graph && length(segment_lengths) == 1
        # All segments same length - use k-mer graph
        k = first(segment_lengths)

        # Detect sequence type from first segment
        Mycelia_module = parentmodule(Rhizomorph)
        first_seq = first(values(segments))
        detected_alphabet = try
            Mycelia_module.detect_alphabet(first_seq)
        catch
            nothing
        end
        if detected_alphabet == :DNA
            kmer_type = Kmers.DNAKmer{k}
        elseif detected_alphabet == :RNA
            kmer_type = Kmers.RNAKmer{k}
        elseif detected_alphabet == :AA
            kmer_type = Kmers.AAKmer{k}
        else
            kmer_type = String
        end

        return read_gfa_next(gfa_file, kmer_type, graph_mode)
    else
        # Variable lengths or forced - use BioSequence graph
        # Use existing sequence detection utilities for type safety
        Mycelia_module = parentmodule(Rhizomorph)
        seqs = Dict{String, Any}()
        seq_types = Set{DataType}()
        for (seg_id, seq_str) in segments
            seq = try
                Mycelia_module.convert_sequence(seq_str)
            catch
                seq_str
            end
            seqs[seg_id] = seq
            push!(seq_types, typeof(seq))
        end

        if isempty(seq_types)
            return MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = StringVertexData,
                edge_data_type = StringEdgeData,
                weight_function = compute_edge_weight
            )
        end
        if length(seq_types) != 1
            error("Mixed sequence types in GFA segments; unable to construct a homogeneous BioSequence graph.")
        end

        SeqType = first(seq_types)

        vertex_data_type = SeqType == String ? StringVertexData :
                           BioSequenceVertexData{SeqType}
        edge_data_type = SeqType == String ? StringEdgeData : BioSequenceEdgeData

        graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph();
            label_type = SeqType,
            vertex_data_type = vertex_data_type,
            edge_data_type = edge_data_type,
            weight_function = compute_edge_weight
        )

        # Add vertices
        for seq in values(seqs)
            if SeqType == String
                graph[seq] = StringVertexData(seq)
            else
                graph[seq] = BioSequenceVertexData(seq)
            end
        end

        # Add edges
        lines = readlines(gfa_file)
        for line in lines
            fields = split(line, '\t')
            if isempty(fields) || first(fields) != "L" || length(fields) < 6
                continue
            end
            src_id = fields[2]
            dst_id = fields[4]
            overlap = parse_overlap_str(fields[6])

            if haskey(seqs, src_id) && haskey(seqs, dst_id)
                src_seq = seqs[src_id]
                dst_seq = seqs[dst_id]
                if SeqType == String
                    graph[src_seq, dst_seq] = StringEdgeData(overlap)
                else
                    graph[src_seq, dst_seq] = BioSequenceEdgeData(overlap)
                end
            end
        end

        return graph
    end
end
