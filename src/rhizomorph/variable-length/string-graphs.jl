# Variable-length String Graph Construction (Singlestrand only)
#
# Overlap-Layout-Consensus (OLC) graphs for variable-length strings.
#
# This module provides graph construction functions for variable-length string
# assembly using the OLC approach where complete strings are vertices and edges
# represent suffix-prefix overlaps.
#
# Key Features:
# - Stores complete strings as vertices (not fixed-length substrings)
# - Odd-length overlaps only (unambiguous resolution)
# - Evidence tracking with nested Dict structure
# - Suitable for assembly from reads, not pattern analysis
#
# Use Cases:
# - Read assembly from sequencing data (after string conversion)
# - Text fragment reconstruction
# - Sequence overlap analysis
#
# Note: Variable-length string graphs are singlestrand only; doublestrand/canonical
# conversions are not defined for generic strings.
#
# Based on rhizomorph-graph-ecosystem-plan.md Section 2.2

# ============================================================================
# Public API - Variable-length String Graph Construction
# ============================================================================

"""
    build_string_graph(strings; dataset_id="dataset_01", min_overlap=3)

Build variable-length string graph using Overlap-Layout-Consensus (OLC) approach.

This is the main user-facing function for variable-length string graph construction.
Unlike n-gram graphs, this stores COMPLETE strings as vertices and connects them
based on suffix-prefix overlaps.

# Arguments
- `strings::Vector{String}`: Input strings (complete sequences)
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `min_overlap::Int=3`: Minimum overlap length (odd-length overlaps only; even values are rounded up)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length string graph with overlap evidence

# Key Concepts
- **Vertices**: Complete strings (not fixed-length substrings)
- **Edges**: Suffix-prefix overlaps between strings
- **Overlap Requirement**: Only odd-length overlaps allowed
- **Evidence**: Tracks which strings overlap and where

# Examples
```julia
# Assembly-like scenario
reads = ["ATCGATCG", "TCGATCGA", "GATCGATT", "CGATTGCC"]
graph = build_string_graph(reads; min_overlap=5)

# Get statistics
stats = get_string_graph_statistics(graph)
println("Found \$(stats[:num_vertices]) strings with \$(stats[:num_edges]) overlaps")

# Find overlap paths
for (src, dst) in MetaGraphsNext.edge_labels(graph)
    edge_data = graph[src, dst]
    println("Overlap: \$src → \$dst (\$(edge_data.overlap_length) chars)")
end
```

# Use Cases

**Read Assembly:**
```julia
# Assemble overlapping reads
reads = ["ATCGATCG", "TCGATCGA", "GATCGATT"]
graph = build_string_graph(reads; min_overlap=5)
# Find assembly paths through the graph
```

**Text Reconstruction:**
```julia
# Reconstruct text from fragments
fragments = ["hello wor", "orld is", "is beautiful"]
graph = build_string_graph(fragments; min_overlap=3)
```

# See Also
- `build_string_graph_from_file`: Build from text file
- `get_string_graph_statistics`: Get graph statistics
"""
function build_string_graph(
    strings::Vector{String};
    dataset_id::String="dataset_01",
    min_overlap::Int=3
)
    return build_string_graph_olc(strings; dataset_id=dataset_id, min_overlap=min_overlap)
end

# Note: The core implementation (build_string_graph_olc, find_overlap_length)
# is defined in core/graph-construction.jl and is available through the parent module.

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    build_token_graph(token_sequences; dataset_id="dataset_01")

Build a token adjacency graph from tokenized sequences.

Vertices are unique token strings; edges connect adjacent tokens observed in each
sequence. This is analogous to n-gram graphs, but uses external tokenization
(e.g., SentencePiece) instead of fixed-length character windows.

# Arguments
- `token_sequences::Vector{Vector{String}}`: Tokenized sequences (one vector per observation)
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking

# Returns
- `MetaGraphsNext.MetaGraph`: Directed graph with token vertices and adjacency edges

# Examples
```julia
tokens = [
    ["The", "cat", "sat"],
    ["The", "dog", "sat"]
]
graph = build_token_graph(tokens; dataset_id="tokens_en")
stats = get_graph_statistics(graph)
println(stats)
```
"""
function build_token_graph(
    token_sequences::Vector{Vector{String}};
    dataset_id::String="dataset_01"
)
    if isempty(token_sequences)
        error("Cannot build token graph from empty token sequence vector")
    end

    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type=String,
        vertex_data_type=StringVertexData,
        edge_data_type=StringEdgeData,
        weight_function=compute_edge_weight
    )

    for (sequence_idx, tokens) in enumerate(token_sequences)
        if isempty(tokens)
            continue
        end
        observation_id = "sequence_$(lpad(sequence_idx, 6, '0'))"

        for (position, token) in enumerate(tokens)
            if !haskey(graph, token)
                graph[token] = StringVertexData(token)
            end
            vertex_data = graph[token]
            add_evidence!(vertex_data, dataset_id, observation_id,
                          EvidenceEntry(position, Forward))
        end

        for i in 1:(length(tokens) - 1)
            src_token = tokens[i]
            dst_token = tokens[i + 1]
            if !haskey(graph, src_token, dst_token)
                graph[src_token, dst_token] = StringEdgeData(0)
            end
            edge_data = graph[src_token, dst_token]
            add_evidence!(edge_data, dataset_id, observation_id,
                          EdgeEvidenceEntry(i, i + 1, Forward))
        end
    end

    return graph
end

"""
    build_string_graph_from_file(filepath; dataset_id=nothing, min_overlap=3)

Build variable-length string graph from a text file.

Reads all lines from the file and builds a string graph.
Each line is treated as a separate string.

# Arguments
- `filepath::String`: Path to text file
- `dataset_id::String=nothing`: Dataset identifier (defaults to filename)
- `min_overlap::Int=3`: Minimum overlap length (odd-length overlaps only; even values are rounded up)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length string graph

# Examples
```julia
# Analyze text file with overlapping lines
graph = build_string_graph_from_file("fragments.txt"; min_overlap=5)

# Read assembly data
graph = build_string_graph_from_file("reads.txt"; min_overlap=7)
```
"""
function build_string_graph_from_file(
    filepath::String;
    dataset_id::Union{String,Nothing}=nothing,
    min_overlap::Int=3
)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    # Use filename as dataset_id if not specified
    if isnothing(dataset_id)
        dataset_id = get_dataset_id_from_file(filepath)
    end

    # Read all lines from file
    strings = readlines(filepath)

    # Filter out empty lines
    strings = filter(!isempty, strings)

    if isempty(strings)
        error("No non-empty lines found in file: $filepath")
    end

    return build_string_graph(strings; dataset_id=dataset_id, min_overlap=min_overlap)
end

"""
    build_string_graph_from_files(filepaths; min_overlap=3)

Build variable-length string graph from multiple text files.

Each file is treated as a separate dataset, using the filename as dataset_id.

# Arguments
- `filepaths::Vector{String}`: List of text files
- `min_overlap::Int=3`: Minimum overlap length (odd-length overlaps only; even values are rounded up)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length string graph with evidence from all files

# Examples
```julia
# Compare multiple fragment sets
files = ["set1.txt", "set2.txt", "set3.txt"]
graph = build_string_graph_from_files(files; min_overlap=5)

# See which strings appear in which datasets
for string_label in MetaGraphsNext.labels(graph)
    vertex = graph[string_label]
    datasets = get_all_dataset_ids(vertex)
    if length(datasets) > 1
        println("String '\$string_label' appears in: \$(join(datasets, ", "))")
    end
end
```
"""
function build_string_graph_from_files(
    filepaths::Vector{String};
    min_overlap::Int=3
)
    if isempty(filepaths)
        error("No files provided")
    end

    # Build graph from first file
    graph = build_string_graph_from_file(filepaths[1]; min_overlap=min_overlap)

    # Add observations from remaining files
    # Note: We rebuild the graph with all strings to find ALL overlaps
    # This is different from incremental addition in fixed-length graphs
    all_strings = String[]
    all_dataset_ids = String[]

    for filepath in filepaths
        dataset_id = get_dataset_id_from_file(filepath)
        strings = filter(!isempty, readlines(filepath))

        for str in strings
            push!(all_strings, str)
            push!(all_dataset_ids, dataset_id)
        end
    end

    # Rebuild with all strings and proper dataset tracking
    return build_string_graph_olc_multi_dataset(all_strings, all_dataset_ids; min_overlap=min_overlap)
end

"""
Helper function for building string graph from multiple datasets.
"""
function build_string_graph_olc_multi_dataset(
    strings::Vector{String},
    dataset_ids::Vector{String};
    min_overlap::Int=3
)
    if isempty(strings)
        error("Cannot build graph from empty string vector")
    end

    if length(strings) != length(dataset_ids)
        error("Number of strings must match number of dataset IDs")
    end

    min_overlap = _normalize_min_overlap(min_overlap)

    # Create empty directed graph
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type=String,
        vertex_data_type=StringVertexData,
        edge_data_type=StringEdgeData,
        weight_function=compute_edge_weight
    )

    # Add all strings as vertices with their dataset IDs
    for (string_idx, (input_string, dataset_id)) in enumerate(zip(strings, dataset_ids))
        observation_id = "string_$(lpad(string_idx, 6, '0'))"

        # Create vertex if it doesn't exist
        if !haskey(graph, input_string)
            vertex_data = StringVertexData(input_string)
            graph[input_string] = vertex_data
        end

        # Add evidence
        vertex_data = graph[input_string]
        add_evidence!(vertex_data, dataset_id, observation_id,
                     EvidenceEntry(1, Forward))
    end

    # Find overlaps between all unique strings
    unique_strings = unique(strings)
    for i in 1:length(unique_strings)
        for j in 1:length(unique_strings)
            if i == j
                continue
            end

            str1 = unique_strings[i]
            str2 = unique_strings[j]

            overlap_len = find_overlap_length(str1, str2, min_overlap)

            if overlap_len >= min_overlap
                if !haskey(graph, str1, str2)
                    edge_data = StringEdgeData(overlap_len)
                    graph[str1, str2] = edge_data
                end

                # Add evidence from all occurrences
                edge_data = graph[str1, str2]
                for (idx, (s, ds_id)) in enumerate(zip(strings, dataset_ids))
                    if s == str1
                        observation_id = "string_$(lpad(idx, 6, '0'))"
                        src_pos = length(str1) - overlap_len + 1
                        dst_pos = 1
                        add_evidence!(edge_data, ds_id, observation_id,
                                    EdgeEvidenceEntry(src_pos, dst_pos, Forward))
                    end
                end
            end
        end
    end

    return graph
end

# ============================================================================
# String Graph Analysis Functions
# ============================================================================

"""
    get_string_graph_statistics(graph)

Get basic statistics about a variable-length string graph.

# Returns
Dictionary with:
- `:num_vertices`: Number of unique strings
- `:num_edges`: Number of overlaps
- `:num_datasets`: Number of datasets with evidence
- `:total_observations`: Total number of string observations
- `:mean_overlap_length`: Mean overlap length across all edges
- `:min_overlap_length`: Minimum overlap length
- `:max_overlap_length`: Maximum overlap length
- `:mean_string_length`: Mean string length across all vertices

# Examples
```julia
stats = get_string_graph_statistics(graph)
println("Found \$(stats[:num_vertices]) strings")
println("Mean overlap: \$(stats[:mean_overlap_length]) chars")
```
"""
function get_string_graph_statistics(graph::MetaGraphsNext.MetaGraph)
    if isempty(MetaGraphsNext.labels(graph))
        return Dict(
            :num_vertices => 0,
            :num_edges => 0,
            :num_datasets => 0,
            :total_observations => 0,
            :mean_overlap_length => nothing,
            :min_overlap_length => nothing,
            :max_overlap_length => nothing,
            :mean_string_length => nothing
        )
    end

    # Collect dataset IDs and count observations
    all_dataset_ids = Set{String}()
    total_observations = 0
    string_lengths = Int[]

    for string_label in MetaGraphsNext.labels(graph)
        vertex_data = graph[string_label]
        union!(all_dataset_ids, keys(vertex_data.evidence))
        total_observations += count_total_observations(vertex_data)
        push!(string_lengths, length(string_label))
    end

    # Collect overlap statistics
    overlap_lengths = Int[]
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        edge_data = graph[src, dst]
        push!(overlap_lengths, edge_data.overlap_length)
    end

    mean_overlap = isempty(overlap_lengths) ? nothing : Statistics.mean(overlap_lengths)
    min_overlap = isempty(overlap_lengths) ? nothing : minimum(overlap_lengths)
    max_overlap = isempty(overlap_lengths) ? nothing : maximum(overlap_lengths)
    mean_str_len = isempty(string_lengths) ? nothing : Statistics.mean(string_lengths)

    return Dict(
        :num_vertices => Graphs.nv(graph.graph),
        :num_edges => Graphs.ne(graph.graph),
        :num_datasets => length(all_dataset_ids),
        :total_observations => total_observations,
        :mean_overlap_length => mean_overlap,
        :min_overlap_length => min_overlap,
        :max_overlap_length => max_overlap,
        :mean_string_length => mean_str_len
    )
end

"""
    find_source_strings(graph)

Find strings that have no incoming edges (potential assembly starting points).

# Returns
- `Vector{String}`: Strings with no predecessors

# Examples
```julia
sources = find_source_strings(graph)
println("Found \$(length(sources)) source strings")
```
"""
function find_source_strings(graph::MetaGraphsNext.MetaGraph)
    sources = String[]

    for string_label in MetaGraphsNext.labels(graph)
        # Check if this vertex has any incoming edges
        has_incoming = false
        for (src, dst) in MetaGraphsNext.edge_labels(graph)
            if dst == string_label
                has_incoming = true
                break
            end
        end

        if !has_incoming
            push!(sources, string_label)
        end
    end

    return sources
end

"""
    find_sink_strings(graph)

Find strings that have no outgoing edges (potential assembly ending points).

# Returns
- `Vector{String}`: Strings with no successors

# Examples
```julia
sinks = find_sink_strings(graph)
println("Found \$(length(sinks)) sink strings")
```
"""
function find_sink_strings(graph::MetaGraphsNext.MetaGraph)
    sinks = String[]

    for string_label in MetaGraphsNext.labels(graph)
        # Check if this vertex has any outgoing edges
        has_outgoing = false
        for (src, dst) in MetaGraphsNext.edge_labels(graph)
            if src == string_label
                has_outgoing = true
                break
            end
        end

        if !has_outgoing
            push!(sinks, string_label)
        end
    end

    return sinks
end
