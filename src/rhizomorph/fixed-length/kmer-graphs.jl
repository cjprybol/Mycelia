# K-mer Graph Construction
#
# Fixed-length k-mer de Bruijn graphs for DNA, RNA, and amino acid sequences.
#
# This module provides graph construction functions that use proper BioSequences
# k-mer types (DNAKmer, RNAKmer, AAKmer) with NO string conversions.
#
# Key Features:
# - Strand-specific by default (stores k-mers as observed)
# - Evidence tracking with nested Dict structure
# - Type-stable core functions for performance
# - Support for SingleStrand (default), DoubleStrand (DNA/RNA), and Canonical (DNA/RNA) modes
#
# Based on rhizomorph-graph-ecosystem-plan.md Section 2.1

# ============================================================================
# Public API - K-mer Graph Construction
# ============================================================================

"""
    build_kmer_graph(records, k; dataset_id="dataset_01", mode=:singlestrand)

Build a k-mer de Bruijn graph from FASTA or FASTQ records.

This is the main user-facing function for k-mer graph construction.

# Arguments
- `records::Vector{<:FASTX.Record}`: Input sequences (FASTA or FASTQ)
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand, :doublestrand, or :canonical)

# Returns
- `MetaGraphsNext.MetaGraph`: K-mer de Bruijn graph with evidence

# Graph Modes
- `:singlestrand` - Records k-mers exactly as observed (strand-specific)
- `:doublestrand` - Adds reverse-complement strand with directed edges (DNA/RNA only)
- `:canonical` - Merges reverse-complement pairs into canonical vertices with undirected edges (DNA/RNA only)

# Examples
```julia
# Load FASTA records
records = collect(open_fastx("reads.fasta"))

# Build strand-specific graph (default)
graph = build_kmer_graph(records, 31)

# Build doublestrand graph
graph_doublestrand = build_kmer_graph(records, 31; mode=:doublestrand)

# Build canonical graph
graph_canonical = build_kmer_graph(records, 31; mode=:canonical)

# Multiple datasets
graph = build_kmer_graph(records1, 31; dataset_id="sample_A")
add_observations_to_graph!(graph, records2, 31; dataset_id="sample_B")
```

# Notes
- Automatically detects sequence type (DNA, RNA, or amino acid)
- Uses type-stable core functions for performance
- FASTQ quality scores are ignored (use build_qualmer_graph for quality-aware graphs)
- All k-mers stored as observed, not canonical (unless mode=:doublestrand)

# See Also
- `build_qualmer_graph`: Quality-aware version for FASTQ data
- `add_observations_to_graph!`: Add more observations to existing graph
"""
function build_kmer_graph(
    records::Vector{<:FASTX.Record},
    k::Int;
    dataset_id::String="dataset_01",
    mode::Symbol=:singlestrand
)
    if mode == :singlestrand
        return build_kmer_graph_singlestrand(records, k; dataset_id=dataset_id)
    elseif mode == :doublestrand
        return build_kmer_graph_doublestrand(records, k; dataset_id=dataset_id)
    elseif mode == :canonical
        singlestrand_graph = build_kmer_graph_singlestrand(records, k; dataset_id=dataset_id)
        return convert_to_canonical(singlestrand_graph)
    else
        error("Invalid mode: $mode. Must be :singlestrand, :doublestrand, or :canonical")
    end
end

# Note: The core implementation functions (build_kmer_graph_singlestrand,
# build_kmer_graph_doublestrand, _build_kmer_graph_core, etc.) are defined
# in core/graph-construction.jl and are available through the parent module.

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    build_kmer_graph_from_file(filepath, k; dataset_id=nothing, mode=:singlestrand)

Build k-mer graph directly from a FASTA or FASTQ file.

Automatically handles compressed files (.gz, .bz2, .xz) and uses filename as
dataset_id if not specified.

# Arguments
- `filepath::String`: Path to FASTA or FASTQ file (compressed or uncompressed)
- `k::Int`: K-mer size
- `dataset_id::String=nothing`: Dataset identifier (defaults to filename without extension)
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand, :doublestrand, or :canonical)

# Returns
- `MetaGraphsNext.MetaGraph`: K-mer de Bruijn graph

# Examples
```julia
# From FASTA file
graph = build_kmer_graph_from_file("reads.fasta", 31)

# From compressed FASTQ file
graph = build_kmer_graph_from_file("reads.fastq.gz", 31)

# With explicit dataset ID
graph = build_kmer_graph_from_file("patient1_tumor.fastq.gz", 31;
                                   dataset_id="patient1_tumor")
```
"""
function build_kmer_graph_from_file(
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

    return build_kmer_graph(records, k; dataset_id=dataset_id, mode=mode)
end

"""
    build_kmer_graph_from_files(filepaths, k; mode=:singlestrand)

Build k-mer graph from multiple FASTA/FASTQ files.

Each file is treated as a separate dataset, using the filename as dataset_id.
Handles compressed files automatically.

# Arguments
- `filepaths::Vector{String}`: List of FASTA/FASTQ files (compressed or uncompressed)
- `k::Int`: K-mer size
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand, :doublestrand, or :canonical)

# Returns
- `MetaGraphsNext.MetaGraph`: K-mer de Bruijn graph with evidence from all files

# Examples
```julia
# Build graph from multiple samples
files = ["sample_A.fasta.gz", "sample_B.fastq.gz", "sample_C.fasta"]
graph = build_kmer_graph_from_files(files, 31)

# Each file's evidence is tracked separately by filename
dataset_ids = get_all_dataset_ids(graph[first(labels(graph))])
# ["sample_A", "sample_B", "sample_C"]
```
"""
function build_kmer_graph_from_files(
    filepaths::Vector{String},
    k::Int;
    mode::Symbol=:singlestrand
)
    if isempty(filepaths)
        error("No files provided")
    end

    if !(mode in (:singlestrand, :doublestrand, :canonical))
        error("Invalid mode: $mode. Must be :singlestrand, :doublestrand, or :canonical")
    end

    # Build graph from first file (strand-specific) and convert later if needed
    graph = build_kmer_graph_from_file(filepaths[1], k; mode=:singlestrand)

    # Get Mycelia module for open_fastx
    Mycelia_module = parentmodule(Rhizomorph)

    # Add observations from remaining files
    for filepath in filepaths[2:end]
        dataset_id = get_dataset_id_from_file(filepath)

        # Use open_fastx which handles compression and format detection
        records = collect(Mycelia_module.open_fastx(filepath))

        # Add observations
        add_observations_to_graph!(graph, records, k; dataset_id=dataset_id)
    end

    # Determine if canonical/doublestrand conversion is valid for the k-mer type
    if !isempty(MetaGraphsNext.labels(graph))
        first_label = first(MetaGraphsNext.labels(graph))
        if mode != :singlestrand && !(first_label isa Kmers.DNAKmer || first_label isa Kmers.RNAKmer)
            error("Mode $mode is only supported for DNA/RNA k-mers")
        end
    end

    # Convert to doublestrand or canonical if requested
    if mode == :doublestrand
        graph = convert_to_doublestrand(graph)
    elseif mode == :canonical
        graph = convert_to_canonical(graph)
    end

    return graph
end

# ============================================================================
# Graph Statistics and Analysis
# ============================================================================

"""
    get_kmer_statistics(graph)

Get basic statistics about a k-mer graph.

# Returns
Dictionary with:
- `:num_vertices`: Number of k-mers in graph
- `:num_edges`: Number of k-mer transitions
- `:num_datasets`: Number of datasets with evidence
- `:total_observations`: Total number of unique observations
- `:k`: K-mer size
- `:sequence_type`: Type of sequences (DNA, RNA, or AA)
- `:graph_mode`: SingleStrand or DoubleStrand

# Examples
```julia
stats = get_kmer_statistics(graph)
println("Graph has \$(stats[:num_vertices]) k-mers from \$(stats[:num_datasets]) datasets")
```
"""
function get_kmer_statistics(graph::MetaGraphsNext.MetaGraph)
    if isempty(MetaGraphsNext.labels(graph))
        return Dict(
            :num_vertices => 0,
            :num_edges => 0,
            :num_datasets => 0,
            :total_observations => 0,
            :k => nothing,
            :sequence_type => nothing,
            :graph_mode => nothing
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
    elseif first_kmer isa Kmers.AAKmer
        "Amino Acid"
    else
        "Unknown"
    end

    # Collect dataset IDs and count observations
    all_dataset_ids = Set{String}()
    total_observations = 0

    for kmer in MetaGraphsNext.labels(graph)
        vertex_data = graph[kmer]
        union!(all_dataset_ids, keys(vertex_data.evidence))
        total_observations += count_total_observations(vertex_data)
    end

    return Dict(
        :num_vertices => Graphs.nv(graph.graph),
        :num_edges => Graphs.ne(graph.graph),
        :num_datasets => length(all_dataset_ids),
        :total_observations => total_observations,
        :k => k,
        :sequence_type => sequence_type,
        :graph_mode => "Unknown"  # Would need to analyze to determine
    )
end

"""
    find_high_coverage_kmers(graph, min_coverage::Int; dataset_id=nothing)

Find k-mers with coverage >= min_coverage.

# Arguments
- `graph`: K-mer graph
- `min_coverage::Int`: Minimum number of observations
- `dataset_id::Union{String,Nothing}=nothing`: Specific dataset or all datasets if nothing

# Returns
- `Vector`: K-mers meeting coverage threshold

# Examples
```julia
# Find k-mers observed at least 10 times across all datasets
high_cov = find_high_coverage_kmers(graph, 10)

# Find k-mers with high coverage in specific dataset
high_cov_sample_a = find_high_coverage_kmers(graph, 5; dataset_id="sample_A")
```
"""
function find_high_coverage_kmers(
    graph::MetaGraphsNext.MetaGraph,
    min_coverage::Int;
    dataset_id::Union{String,Nothing}=nothing
)
    result = []

    for kmer in MetaGraphsNext.labels(graph)
        vertex_data = graph[kmer]

        coverage = if isnothing(dataset_id)
            count_total_observations(vertex_data)
        else
            count_dataset_observations(vertex_data, dataset_id)
        end

        if coverage >= min_coverage
            push!(result, kmer)
        end
    end

    return result
end
