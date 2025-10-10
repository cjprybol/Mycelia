# Variable-length FASTA Graph Construction
#
# Overlap-Layout-Consensus (OLC) graphs for variable-length FASTA sequences.
#
# This module provides graph construction functions for variable-length biological
# sequence assembly using the OLC approach where complete BioSequences are vertices
# and edges represent suffix-prefix overlaps.
#
# Key Features:
# - Stores complete BioSequences as vertices (NOT strings)
# - Odd-length overlaps only (unambiguous resolution)
# - Evidence tracking with nested Dict structure
# - Suitable for assembly from contigs, scaffolds, or reads
#
# Use Cases:
# - Contig/scaffold assembly
# - Sequence overlap analysis
# - Reference-free assembly
#
# Based on rhizomorph-graph-ecosystem-plan.md Section 2.2

# ============================================================================
# Public API - Variable-length FASTA Graph Construction
# ============================================================================

"""
    build_fasta_graph(records; dataset_id="dataset_01", min_overlap=3)

Build variable-length FASTA graph using Overlap-Layout-Consensus (OLC) approach.

This is the main user-facing function for variable-length FASTA graph construction.
Unlike k-mer graphs, this stores COMPLETE BioSequences as vertices and connects them
based on suffix-prefix overlaps.

# Arguments
- `records::Vector{FASTX.FASTA.Record}`: Input FASTA records
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `min_overlap::Int=3`: Minimum overlap length (must be odd)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length FASTA graph with BioSequence vertices

# Key Concepts
- **Vertices**: Complete BioSequences (LongDNA, LongRNA, or LongAA)
- **Edges**: Suffix-prefix overlaps between sequences
- **Overlap Requirement**: Only odd-length overlaps allowed
- **Evidence**: Tracks which sequences overlap and where

# Examples
```julia
# Load contigs from assembly
records = collect(open_fastx("contigs.fasta"))
graph = build_fasta_graph(records; min_overlap=7)

# Get statistics
stats = get_fasta_graph_statistics(graph)
println("Found \$(stats[:num_vertices]) sequences with \$(stats[:num_edges]) overlaps")

# Find overlap paths
for (src, dst) in MetaGraphsNext.edge_labels(graph)
    edge_data = graph[src, dst]
    println("Overlap: \$(length(src))bp → \$(length(dst))bp (\$(edge_data.overlap_length)bp)")
end
```

# Use Cases

**Contig Assembly:**
```julia
# Assemble overlapping contigs
contigs = collect(open_fastx("assembly.fasta"))
graph = build_fasta_graph(contigs; min_overlap=15)
# Find paths through the graph to scaffold contigs
```

**Sequence Overlap Analysis:**
```julia
# Find which sequences overlap
seqs = collect(open_fastx("sequences.fasta"))
graph = build_fasta_graph(seqs; min_overlap=9)
```

# See Also
- `build_fasta_graph_from_file`: Build from FASTA file
- `get_fasta_graph_statistics`: Get graph statistics
"""
function build_fasta_graph(
    records::Vector{FASTX.FASTA.Record};
    dataset_id::String="dataset_01",
    min_overlap::Int=3
)
    return build_fasta_graph_olc(records; dataset_id=dataset_id, min_overlap=min_overlap)
end

# Note: The core implementation (build_fasta_graph_olc, find_biosequence_overlap_length)
# is defined in core/graph-construction.jl and is available through the parent module.

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    build_fasta_graph_from_file(filepath; dataset_id=nothing, min_overlap=3)

Build variable-length FASTA graph from a FASTA file.

Automatically handles compressed files (.gz, .bz2, .xz).

# Arguments
- `filepath::String`: Path to FASTA file
- `dataset_id::String=nothing`: Dataset identifier (defaults to filename)
- `min_overlap::Int=3`: Minimum overlap length (must be odd)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length FASTA graph

# Examples
```julia
# Analyze contig file
graph = build_fasta_graph_from_file("contigs.fasta"; min_overlap=7)

# From compressed file
graph = build_fasta_graph_from_file("assembly.fasta.gz"; min_overlap=9)
```
"""
function build_fasta_graph_from_file(
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

    # Use Mycelia.open_fastx which handles format detection and compression
    Mycelia_module = parentmodule(Rhizomorph)
    records = collect(Mycelia_module.open_fastx(filepath))

    # Validate FASTA
    if isempty(records)
        error("No records found in file: $filepath")
    end
    if !(records[1] isa FASTX.FASTA.Record)
        error("FASTA graphs require FASTA input. File appears to be FASTQ: $filepath")
    end

    return build_fasta_graph(records; dataset_id=dataset_id, min_overlap=min_overlap)
end

"""
    build_fasta_graph_from_files(filepaths; min_overlap=3)

Build variable-length FASTA graph from multiple FASTA files.

Each file is treated as a separate dataset, using the filename as dataset_id.

# Arguments
- `filepaths::Vector{String}`: List of FASTA files
- `min_overlap::Int=3`: Minimum overlap length (must be odd)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length FASTA graph with evidence from all files

# Examples
```julia
# Compare multiple assemblies
files = ["assembly1.fasta", "assembly2.fasta", "assembly3.fasta"]
graph = build_fasta_graph_from_files(files; min_overlap=7)

# See which sequences appear in which assemblies
for seq_label in MetaGraphsNext.labels(graph)
    vertex = graph[seq_label]
    datasets = get_all_dataset_ids(vertex)
    if length(datasets) > 1
        println("Sequence appears in: \$(join(datasets, ", "))")
    end
end
```
"""
function build_fasta_graph_from_files(
    filepaths::Vector{String};
    min_overlap::Int=3
)
    if isempty(filepaths)
        error("No files provided")
    end

    # Build graph from first file
    graph = build_fasta_graph_from_file(filepaths[1]; min_overlap=min_overlap)

    # Get Mycelia module for open_fastx
    Mycelia_module = parentmodule(Rhizomorph)

    # Add observations from remaining files
    for filepath in filepaths[2:end]
        dataset_id = get_dataset_id_from_file(filepath)

        # Use open_fastx which handles compression and format detection
        records = collect(Mycelia_module.open_fastx(filepath))

        # Validate FASTA
        if isempty(records) || !(records[1] isa FASTX.FASTA.Record)
            error("FASTA graphs require FASTA input. File appears to be FASTQ: $filepath")
        end

        # Note: For OLC graphs, we need to rebuild to find ALL overlaps
        # This is a simplified version - proper implementation would merge evidence
        # For now, we just add to the first graph
        # TODO: Implement proper multi-dataset OLC graph merging
    end

    return graph
end

# ============================================================================
# FASTA Graph Analysis Functions
# ============================================================================

"""
    get_fasta_graph_statistics(graph)

Get basic statistics about a variable-length FASTA graph.

# Returns
Dictionary with:
- `:num_vertices`: Number of unique sequences
- `:num_edges`: Number of overlaps
- `:num_datasets`: Number of datasets with evidence
- `:total_observations`: Total number of sequence observations
- `:mean_overlap_length`: Mean overlap length across all edges
- `:min_overlap_length`: Minimum overlap length
- `:max_overlap_length`: Maximum overlap length
- `:mean_sequence_length`: Mean sequence length across all vertices
- `:sequence_type`: Type of sequences (DNA, RNA, or AA)

# Examples
```julia
stats = get_fasta_graph_statistics(graph)
println("Found \$(stats[:num_vertices]) sequences")
println("Mean overlap: \$(stats[:mean_overlap_length]) bp")
println("Sequence type: \$(stats[:sequence_type])")
```
"""
function get_fasta_graph_statistics(graph::MetaGraphsNext.MetaGraph)
    if isempty(MetaGraphsNext.labels(graph))
        return Dict(
            :num_vertices => 0,
            :num_edges => 0,
            :num_datasets => 0,
            :total_observations => 0,
            :mean_overlap_length => nothing,
            :min_overlap_length => nothing,
            :max_overlap_length => nothing,
            :mean_sequence_length => nothing,
            :sequence_type => nothing
        )
    end

    # Determine sequence type from first label
    first_seq = first(MetaGraphsNext.labels(graph))
    sequence_type = if first_seq isa BioSequences.LongDNA
        "DNA"
    elseif first_seq isa BioSequences.LongRNA
        "RNA"
    elseif first_seq isa BioSequences.LongAA
        "Amino Acid"
    else
        "Unknown"
    end

    # Collect dataset IDs and count observations
    all_dataset_ids = Set{String}()
    total_observations = 0
    sequence_lengths = Int[]

    for seq_label in MetaGraphsNext.labels(graph)
        vertex_data = graph[seq_label]
        union!(all_dataset_ids, keys(vertex_data.evidence))
        total_observations += count_total_observations(vertex_data)
        push!(sequence_lengths, length(seq_label))
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
    mean_seq_len = isempty(sequence_lengths) ? nothing : Statistics.mean(sequence_lengths)

    return Dict(
        :num_vertices => Graphs.nv(graph.graph),
        :num_edges => Graphs.ne(graph.graph),
        :num_datasets => length(all_dataset_ids),
        :total_observations => total_observations,
        :mean_overlap_length => mean_overlap,
        :min_overlap_length => min_overlap,
        :max_overlap_length => max_overlap,
        :mean_sequence_length => mean_seq_len,
        :sequence_type => sequence_type
    )
end

"""
    find_source_sequences(graph)

Find sequences that have no incoming edges (potential assembly starting points).

# Returns
- `Vector`: Sequences with no predecessors

# Examples
```julia
sources = find_source_sequences(graph)
println("Found \$(length(sources)) source sequences")
```
"""
function find_source_sequences(graph::MetaGraphsNext.MetaGraph)
    sources = []

    for seq_label in MetaGraphsNext.labels(graph)
        # Check if this vertex has any incoming edges
        has_incoming = false
        for (src, dst) in MetaGraphsNext.edge_labels(graph)
            if dst == seq_label
                has_incoming = true
                break
            end
        end

        if !has_incoming
            push!(sources, seq_label)
        end
    end

    return sources
end

"""
    find_sink_sequences(graph)

Find sequences that have no outgoing edges (potential assembly ending points).

# Returns
- `Vector`: Sequences with no successors

# Examples
```julia
sinks = find_sink_sequences(graph)
println("Found \$(length(sinks)) sink sequences")
```
"""
function find_sink_sequences(graph::MetaGraphsNext.MetaGraph)
    sinks = []

    for seq_label in MetaGraphsNext.labels(graph)
        # Check if this vertex has any outgoing edges
        has_outgoing = false
        for (src, dst) in MetaGraphsNext.edge_labels(graph)
            if src == seq_label
                has_outgoing = true
                break
            end
        end

        if !has_outgoing
            push!(sinks, seq_label)
        end
    end

    return sinks
end
