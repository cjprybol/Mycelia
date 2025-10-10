# Qualmer Graph Construction
#
# Fixed-length quality-aware k-mer de Bruijn graphs for DNA and RNA sequences.
#
# This module provides graph construction functions that preserve per-base quality
# scores throughout the assembly process, enabling quality-aware algorithms.
#
# Key Features:
# - Preserves Phred quality scores for each k-mer observation
# - Strand-specific by default (stores k-mers as observed)
# - Evidence tracking with nested Dict structure including quality data
# - Type-stable core functions for performance
# - Support for both SingleStrand and DoubleStrand modes
#
# Based on rhizomorph-graph-ecosystem-plan.md Section 2.1

# ============================================================================
# Public API - Qualmer Graph Construction
# ============================================================================

"""
    build_qualmer_graph(records, k; dataset_id="dataset_01", mode=:singlestrand)

Build a quality-aware k-mer de Bruijn graph from FASTQ records.

This is the main user-facing function for qualmer (quality k-mer) graph construction.
Preserves per-base quality scores for quality-aware assembly algorithms.

# Arguments
- `records::Vector{FASTX.FASTQ.Record}`: Input FASTQ records (MUST be FASTQ for quality)
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand or :doublestrand)

# Returns
- `MetaGraphsNext.MetaGraph`: Quality-aware k-mer de Bruijn graph with Phred scores

# Graph Modes
- `:singlestrand` - Records k-mers exactly as observed (strand-specific)
- `:doublestrand` - Merges forward and reverse-complement k-mers (canonical)

# Examples
```julia
# Load FASTQ records
records = collect(open_fastx("reads.fastq.gz"))

# Build strand-specific qualmer graph (default)
graph = build_qualmer_graph(records, 31)

# Build canonical (doublestrand) qualmer graph
graph_canonical = build_qualmer_graph(records, 31; mode=:doublestrand)

# Multiple datasets
graph = build_qualmer_graph(records1, 31; dataset_id="sample_A")
add_observations_to_graph!(graph, records2, 31; dataset_id="sample_B")

# Get joint quality for a k-mer
kmer = first(labels(graph))
vertex_data = graph[kmer]
joint_quality = get_vertex_joint_quality(vertex_data, "sample_A")
```

# Quality Score Handling
- Phred scores are stored as raw UInt8 values (0-60 for single observations)
- Joint quality across multiple observations can exceed 60 (up to 255)
- Quality scores add in log space: Q_combined = Q1 + Q2 + ... + Qn
- This allows differentiation between low (10-30), moderate (30-60),
  high (60-100), and very high (100+) confidence

# Notes
- REQUIRES FASTQ input (will error on FASTA)
- Automatically detects sequence type (DNA or RNA only, not amino acids)
- Uses type-stable core functions for performance
- All k-mers stored as observed, not canonical (unless mode=:doublestrand)

# See Also
- `build_kmer_graph`: Non-quality version for FASTA data
- `get_vertex_joint_quality`: Compute joint quality from multiple observations
- `combine_phred_scores`: Combine quality scores from independent observations
"""
function build_qualmer_graph(
    records::Vector{FASTX.FASTQ.Record},
    k::Int;
    dataset_id::String="dataset_01",
    mode::Symbol=:singlestrand
)
    if mode == :singlestrand
        return build_qualmer_graph_singlestrand(records, k; dataset_id=dataset_id)
    elseif mode == :doublestrand
        return build_qualmer_graph_doublestrand(records, k; dataset_id=dataset_id)
    else
        error("Invalid mode: $mode. Must be :singlestrand or :doublestrand")
    end
end

# Note: The core implementation functions (build_qualmer_graph_singlestrand,
# build_qualmer_graph_doublestrand, _build_qualmer_graph_core, etc.) are defined
# in core/graph-construction.jl and are available through the parent module.

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    build_qualmer_graph_from_file(filepath, k; dataset_id=nothing, mode=:singlestrand)

Build qualmer graph directly from a FASTQ file.

Automatically handles compressed files (.gz, .bz2, .xz) and uses filename as
dataset_id if not specified.

# Arguments
- `filepath::String`: Path to FASTQ file (compressed or uncompressed)
- `k::Int`: K-mer size
- `dataset_id::String=nothing`: Dataset identifier (defaults to filename without extension)
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand or :doublestrand)

# Returns
- `MetaGraphsNext.MetaGraph`: Quality-aware k-mer de Bruijn graph

# Examples
```julia
# From FASTQ file
graph = build_qualmer_graph_from_file("reads.fastq", 31)

# From compressed FASTQ file
graph = build_qualmer_graph_from_file("reads.fastq.gz", 31)

# With explicit dataset ID
graph = build_qualmer_graph_from_file("patient1_tumor.fastq.gz", 31;
                                      dataset_id="patient1_tumor")
```
"""
function build_qualmer_graph_from_file(
    filepath::String,
    k::Int;
    dataset_id::Union{String,Nothing}=nothing,
    mode::Symbol=:singlestrand
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

    # Validate that we got FASTQ records
    if isempty(records)
        error("No records found in file: $filepath")
    end
    if !(records[1] isa FASTX.FASTQ.Record)
        error("Qualmer graphs require FASTQ input. File appears to be FASTA: $filepath")
    end

    return build_qualmer_graph(records, k; dataset_id=dataset_id, mode=mode)
end

"""
    build_qualmer_graph_from_files(filepaths, k; mode=:singlestrand)

Build qualmer graph from multiple FASTQ files.

Each file is treated as a separate dataset, using the filename as dataset_id.
Handles compressed files automatically.

# Arguments
- `filepaths::Vector{String}`: List of FASTQ files (compressed or uncompressed)
- `k::Int`: K-mer size
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand or :doublestrand)

# Returns
- `MetaGraphsNext.MetaGraph`: Quality-aware k-mer de Bruijn graph with evidence from all files

# Examples
```julia
# Build graph from multiple samples
files = ["sample_A.fastq.gz", "sample_B.fastq.gz", "sample_C.fastq.gz"]
graph = build_qualmer_graph_from_files(files, 31)

# Each file's evidence is tracked separately by filename
dataset_ids = get_all_dataset_ids(graph[first(labels(graph))])
# ["sample_A", "sample_B", "sample_C"]
```
"""
function build_qualmer_graph_from_files(
    filepaths::Vector{String},
    k::Int;
    mode::Symbol=:singlestrand
)
    if isempty(filepaths)
        error("No files provided")
    end

    # Build graph from first file
    graph = build_qualmer_graph_from_file(filepaths[1], k; mode=:singlestrand)

    # Get Mycelia module for open_fastx
    Mycelia_module = parentmodule(Rhizomorph)

    # Add observations from remaining files
    for filepath in filepaths[2:end]
        dataset_id = get_dataset_id_from_file(filepath)

        # Use open_fastx which handles compression and format detection
        records = collect(Mycelia_module.open_fastx(filepath))

        # Validate FASTQ
        if isempty(records) || !(records[1] isa FASTX.FASTQ.Record)
            error("Qualmer graphs require FASTQ input. File appears to be FASTA: $filepath")
        end

        # Add observations
        add_observations_to_graph!(graph, records, k; dataset_id=dataset_id)
    end

    # Convert to doublestrand if requested
    if mode == :doublestrand
        graph = convert_to_doublestrand(graph)
    end

    return graph
end

# ============================================================================
# Quality-Specific Analysis Functions
# ============================================================================

"""
    get_qualmer_statistics(graph; dataset_id=nothing)

Get statistics about a qualmer graph including quality score information.

# Arguments
- `graph`: Qualmer graph
- `dataset_id::Union{String,Nothing}=nothing`: Specific dataset or all datasets

# Returns
Dictionary with:
- `:num_vertices`: Number of k-mers in graph
- `:num_edges`: Number of k-mer transitions
- `:num_datasets`: Number of datasets with evidence
- `:total_observations`: Total number of unique observations
- `:k`: K-mer size
- `:sequence_type`: Type of sequences (DNA or RNA)
- `:mean_joint_quality`: Mean of joint quality scores across all k-mers
- `:min_joint_quality`: Minimum joint quality score
- `:max_joint_quality`: Maximum joint quality score

# Examples
```julia
stats = get_qualmer_statistics(graph)
println("Mean quality: \$(stats[:mean_joint_quality])")

# For specific dataset
stats_sample_a = get_qualmer_statistics(graph; dataset_id="sample_A")
```
"""
function get_qualmer_statistics(
    graph::MetaGraphsNext.MetaGraph;
    dataset_id::Union{String,Nothing}=nothing
)
    if isempty(MetaGraphsNext.labels(graph))
        return Dict(
            :num_vertices => 0,
            :num_edges => 0,
            :num_datasets => 0,
            :total_observations => 0,
            :k => nothing,
            :sequence_type => nothing,
            :mean_joint_quality => nothing,
            :min_joint_quality => nothing,
            :max_joint_quality => nothing
        )
    end

    # Get first k-mer to determine type and size
    first_kmer = first(MetaGraphsNext.labels(graph))
    k = Kmers.ksize(typeof(first_kmer))

    # Determine sequence type
    sequence_type = if first_kmer isa Kmers.DNAKmer
        "DNA"
    elseif first_kmer isa Kmers.RNAKmer
        "RNA"
    else
        "Unknown"
    end

    # Collect dataset IDs and count observations
    all_dataset_ids = Set{String}()
    total_observations = 0

    # Collect quality statistics
    quality_scores = Float64[]

    for kmer in MetaGraphsNext.labels(graph)
        vertex_data = graph[kmer]
        union!(all_dataset_ids, keys(vertex_data.evidence))
        total_observations += count_total_observations(vertex_data)

        # Get joint quality for this k-mer
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

    mean_qual = isempty(quality_scores) ? nothing : Statistics.mean(quality_scores)
    min_qual = isempty(quality_scores) ? nothing : minimum(quality_scores)
    max_qual = isempty(quality_scores) ? nothing : maximum(quality_scores)

    return Dict(
        :num_vertices => Graphs.nv(graph.graph),
        :num_edges => Graphs.ne(graph.graph),
        :num_datasets => length(all_dataset_ids),
        :total_observations => total_observations,
        :k => k,
        :sequence_type => sequence_type,
        :mean_joint_quality => mean_qual,
        :min_joint_quality => min_qual,
        :max_joint_quality => max_qual
    )
end

"""
    find_high_quality_kmers(graph, min_quality::Int; dataset_id=nothing)

Find k-mers with joint quality >= min_quality.

# Arguments
- `graph`: Qualmer graph
- `min_quality::Int`: Minimum joint quality score (Phred scale)
- `dataset_id::Union{String,Nothing}=nothing`: Specific dataset or all datasets if nothing

# Returns
- `Vector`: K-mers meeting quality threshold

# Examples
```julia
# Find high-quality k-mers (joint Q >= 60)
high_qual = find_high_quality_kmers(graph, 60)

# Find k-mers with very high confidence in specific dataset (Q >= 100)
very_high_qual = find_high_quality_kmers(graph, 100; dataset_id="sample_A")
```
"""
function find_high_quality_kmers(
    graph::MetaGraphsNext.MetaGraph,
    min_quality::Int;
    dataset_id::Union{String,Nothing}=nothing
)
    result = []

    for kmer in MetaGraphsNext.labels(graph)
        vertex_data = graph[kmer]

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
            push!(result, kmer)
        end
    end

    return result
end
