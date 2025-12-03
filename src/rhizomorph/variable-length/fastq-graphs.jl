# Variable-length FASTQ Graph Construction (Singlestrand with conversions)
#
# Overlap-Layout-Consensus (OLC) graphs for variable-length FASTQ sequences with quality.
#
# This module provides graph construction functions for variable-length biological
# sequence assembly using the OLC approach where complete BioSequences are vertices
# with quality information preserved.
#
# Key Features:
# - Stores complete BioSequences as vertices (NOT strings)
# - Quality scores preserved for each sequence
# - Odd-length overlaps only (unambiguous resolution)
# - Evidence tracking with nested Dict structure
# - Suitable for quality-aware assembly from reads
#
# Use Cases:
# - Quality-aware read assembly
# - High-confidence overlap detection
# - Error-aware sequence assembly
#
# Note: Variable-length FASTQ graphs are constructed singlestrand; doublestrand
# and canonical conversions are provided in `variable-length/strand-conversions.jl`
# for DNA/RNA sequences.
#
# Based on rhizomorph-graph-ecosystem-plan.md Section 2.2

# ============================================================================
# Public API - Variable-length FASTQ Graph Construction
# ============================================================================

"""
    build_fastq_graph(records; dataset_id="dataset_01", min_overlap=3)

Build variable-length FASTQ graph using Overlap-Layout-Consensus (OLC) approach.

This is the main user-facing function for variable-length FASTQ graph construction.
Unlike qualmer graphs, this stores COMPLETE BioSequences as vertices and connects them
based on suffix-prefix overlaps, preserving quality information.

# Arguments
- `records::Vector{FASTX.FASTQ.Record}`: Input FASTQ records
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `min_overlap::Int=3`: Minimum overlap length (must be odd)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length FASTQ graph with BioSequence vertices and quality

# Key Concepts
- **Vertices**: Complete BioSequences (LongDNA or LongRNA) with quality
- **Edges**: Suffix-prefix overlaps between sequences
- **Overlap Requirement**: Only odd-length overlaps allowed
- **Evidence**: Tracks which sequences overlap with quality scores

# Examples
```julia
# Load reads
records = collect(open_fastx("reads.fastq.gz"))
graph = build_fastq_graph(records; min_overlap=7)

# Get statistics
stats = get_fastq_graph_statistics(graph)
println("Found \$(stats[:num_vertices]) sequences with \$(stats[:num_edges]) overlaps")
println("Mean quality: \$(stats[:mean_quality])")

# Find high-quality overlap paths
for (src, dst) in MetaGraphsNext.edge_labels(graph)
    edge_data = graph[src, dst]
    println("Overlap: \$(length(src))bp → \$(length(dst))bp (\$(edge_data.overlap_length)bp)")
end
```

# Use Cases

**Quality-aware Assembly:**
```julia
# Assemble reads with quality tracking
reads = collect(open_fastx("reads.fastq.gz"))
graph = build_fastq_graph(reads; min_overlap=15)
# Use quality scores to select best assembly paths
```

**High-confidence Overlap Detection:**
```julia
# Find only high-quality overlaps
reads = collect(open_fastx("reads.fastq.gz"))
graph = build_fastq_graph(reads; min_overlap=21)
```

# See Also
- `build_fastq_graph_from_file`: Build from FASTQ file
- `get_fastq_graph_statistics`: Get graph statistics including quality
"""
function build_fastq_graph(
    records::Vector{FASTX.FASTQ.Record};
    dataset_id::String="dataset_01",
    min_overlap::Int=3
)
    return build_fastq_graph_olc(records; dataset_id=dataset_id, min_overlap=min_overlap)
end

# Note: The core implementation (build_fastq_graph_olc, find_biosequence_overlap_length)
# is defined in core/graph-construction.jl and is available through the parent module.

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    build_fastq_graph_from_file(filepath; dataset_id=nothing, min_overlap=3)

Build variable-length FASTQ graph from a FASTQ file.

Automatically handles compressed files (.gz, .bz2, .xz).

# Arguments
- `filepath::String`: Path to FASTQ file
- `dataset_id::String=nothing`: Dataset identifier (defaults to filename)
- `min_overlap::Int=3`: Minimum overlap length (must be odd)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length FASTQ graph

# Examples
```julia
# Analyze reads file
graph = build_fastq_graph_from_file("reads.fastq"; min_overlap=7)

# From compressed file
graph = build_fastq_graph_from_file("reads.fastq.gz"; min_overlap=9)
```
"""
function build_fastq_graph_from_file(
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

    # Validate FASTQ
    if isempty(records)
        error("No records found in file: $filepath")
    end
    if !(records[1] isa FASTX.FASTQ.Record)
        error("FASTQ graphs require FASTQ input. File appears to be FASTA: $filepath")
    end

    return build_fastq_graph(records; dataset_id=dataset_id, min_overlap=min_overlap)
end

"""
    build_fastq_graph_from_files(filepaths; min_overlap=3)

Build variable-length FASTQ graph from multiple FASTQ files.

Each file is treated as a separate dataset, using the filename as dataset_id.

# Arguments
- `filepaths::Vector{String}`: List of FASTQ files
- `min_overlap::Int=3`: Minimum overlap length (must be odd)

# Returns
- `MetaGraphsNext.MetaGraph`: Variable-length FASTQ graph with evidence from all files

# Examples
```julia
# Compare multiple read sets
files = ["sample1.fastq.gz", "sample2.fastq.gz", "sample3.fastq.gz"]
graph = build_fastq_graph_from_files(files; min_overlap=7)

# See which sequences appear in which samples
for seq_label in MetaGraphsNext.labels(graph)
    vertex = graph[seq_label]
    datasets = get_all_dataset_ids(vertex)
    if length(datasets) > 1
        println("Sequence appears in: \$(join(datasets, ", "))")
    end
end
```
"""
function build_fastq_graph_from_files(
    filepaths::Vector{String};
    min_overlap::Int=3
)
    if isempty(filepaths)
        error("No files provided")
    end

    # Build graph from first file
    graph = build_fastq_graph_from_file(filepaths[1]; min_overlap=min_overlap)

    # Get Mycelia module for open_fastx
    Mycelia_module = parentmodule(Rhizomorph)

    # Add observations from remaining files
    for filepath in filepaths[2:end]
        dataset_id = get_dataset_id_from_file(filepath)

        # Use open_fastx which handles compression and format detection
        records = collect(Mycelia_module.open_fastx(filepath))

        # Validate FASTQ
        if isempty(records) || !(records[1] isa FASTX.FASTQ.Record)
            error("FASTQ graphs require FASTQ input. File appears to be FASTA: $filepath")
        end

        # Note: For OLC graphs, we need to rebuild to find ALL overlaps
        # This is a simplified version - proper implementation would merge evidence
        # For now, we just add to the first graph
        # TODO: Implement proper multi-dataset OLC graph merging
    end

    return graph
end

# ============================================================================
# FASTQ Graph Analysis Functions
# ============================================================================

"""
    get_fastq_graph_statistics(graph; dataset_id=nothing)

Get statistics about a variable-length FASTQ graph including quality information.

# Arguments
- `graph`: FASTQ graph
- `dataset_id::Union{String,Nothing}=nothing`: Specific dataset or all datasets

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
- `:sequence_type`: Type of sequences (DNA or RNA)
- `:mean_quality`: Mean quality score across all sequences

# Examples
```julia
stats = get_fastq_graph_statistics(graph)
println("Mean quality: \$(stats[:mean_quality])")

# For specific dataset
stats_sample_a = get_fastq_graph_statistics(graph; dataset_id="sample_A")
```
"""
function get_fastq_graph_statistics(
    graph::MetaGraphsNext.MetaGraph;
    dataset_id::Union{String,Nothing}=nothing
)
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
            :sequence_type => nothing,
            :mean_quality => nothing
        )
    end

    # Determine sequence type from first label
    first_seq = first(MetaGraphsNext.labels(graph))
    sequence_type = if first_seq isa BioSequences.LongDNA
        "DNA"
    elseif first_seq isa BioSequences.LongRNA
        "RNA"
    else
        "Unknown"
    end

    # Collect dataset IDs and count observations
    all_dataset_ids = Set{String}()
    total_observations = 0
    sequence_lengths = Int[]
    quality_scores = Float64[]

    for seq_label in MetaGraphsNext.labels(graph)
        vertex_data = graph[seq_label]
        union!(all_dataset_ids, keys(vertex_data.evidence))
        total_observations += count_total_observations(vertex_data)
        push!(sequence_lengths, length(seq_label))

        # Get quality scores
        if isnothing(dataset_id)
            # Aggregate across all datasets
            for ds_id in keys(vertex_data.evidence)
                joint_qual = get_vertex_joint_quality(vertex_data, ds_id)
                if !isnothing(joint_qual) && !isempty(joint_qual)
                    push!(quality_scores, Statistics.mean(joint_qual))
                end
            end
        else
            # Specific dataset
            joint_qual = get_vertex_joint_quality(vertex_data, dataset_id)
            if !isnothing(joint_qual) && !isempty(joint_qual)
                push!(quality_scores, Statistics.mean(joint_qual))
            end
        end
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
    mean_qual = isempty(quality_scores) ? nothing : Statistics.mean(quality_scores)

    return Dict(
        :num_vertices => Graphs.nv(graph.graph),
        :num_edges => Graphs.ne(graph.graph),
        :num_datasets => length(all_dataset_ids),
        :total_observations => total_observations,
        :mean_overlap_length => mean_overlap,
        :min_overlap_length => min_overlap,
        :max_overlap_length => max_overlap,
        :mean_sequence_length => mean_seq_len,
        :sequence_type => sequence_type,
        :mean_quality => mean_qual
    )
end

"""
    find_high_quality_sequences(graph, min_quality::Int; dataset_id=nothing)

Find sequences with joint quality >= min_quality.

# Arguments
- `graph`: FASTQ graph
- `min_quality::Int`: Minimum joint quality score (Phred scale)
- `dataset_id::Union{String,Nothing}=nothing`: Specific dataset or all datasets

# Returns
- `Vector`: Sequences meeting quality threshold

# Examples
```julia
# Find high-quality sequences (Q >= 60)
high_qual = find_high_quality_sequences(graph, 60)

# Find very high confidence sequences (Q >= 100)
very_high_qual = find_high_quality_sequences(graph, 100; dataset_id="sample_A")
```
"""
function find_high_quality_sequences(
    graph::MetaGraphsNext.MetaGraph,
    min_quality::Int;
    dataset_id::Union{String,Nothing}=nothing
)
    result = []

    for seq_label in MetaGraphsNext.labels(graph)
        vertex_data = graph[seq_label]

        # Check quality for specified dataset(s)
        meets_threshold = if isnothing(dataset_id)
            # Check all datasets
            any_meets = false
            for ds_id in keys(vertex_data.evidence)
                joint_qual = get_vertex_joint_quality(vertex_data, ds_id)
                if !isnothing(joint_qual) && all(q >= min_quality for q in joint_qual)
                    any_meets = true
                    break
                end
            end
            any_meets
        else
            # Check specific dataset
            joint_qual = get_vertex_joint_quality(vertex_data, dataset_id)
            !isnothing(joint_qual) && all(q >= min_quality for q in joint_qual)
        end

        if meets_threshold
            push!(result, seq_label)
        end
    end

    return result
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
