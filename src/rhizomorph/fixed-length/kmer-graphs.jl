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
    build_kmer_graph(records, k; dataset_id="dataset_01", mode=:singlestrand,
                     type_hint=nothing, ambiguous_action=:dna,
                     memory_profile=:full, lightweight=nothing)

Build a k-mer de Bruijn graph from FASTA or FASTQ records.

This is the main user-facing function for k-mer graph construction.

# Arguments
- `records::Vector{<:FASTX.Record}`: Input sequences (FASTA or FASTQ)
- `k::Int`: K-mer size
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand, :doublestrand, or :canonical)
- `type_hint::Union{Nothing,Symbol}=nothing`: Optional alphabet hint (:DNA, :RNA, :AA)
- `ambiguous_action::Symbol=:dna`: Resolution for ambiguous alphabets
- `memory_profile::Symbol=:full`: Memory profile controlling information retention:
  - `:full` — Full evidence (positions, strands, per-observation)
  - `:lightweight` — Counts + observation ID tracking (no positions)
  - `:ultralight` — Counts only (no observation IDs, no positions)
  - `:ultralight_quality` — Counts + aggregated quality scores (FASTQ only)
  - `:lightweight_quality` — Counts + obs IDs + aggregated quality (FASTQ only)
- `lightweight::Union{Nothing,Bool}=nothing`: Deprecated. Use `memory_profile` instead.
  `lightweight=true` maps to `memory_profile=:lightweight`.

# Returns
- `MetaGraphsNext.MetaGraph`: K-mer de Bruijn graph with evidence

# Graph Modes
- `:singlestrand` - Records k-mers exactly as observed (strand-specific)
- `:doublestrand` - Adds reverse-complement strand with directed edges (DNA/RNA only)
- `:canonical` - Merges reverse-complement pairs into canonical vertices (DNA/RNA only)

# Memory Profile Selection Guide
```
Need individual read positions/strands?
  YES → :full
  NO  → Need to know WHICH reads observed each k-mer?
          YES → Need quality scores?
                  YES → :lightweight_quality
                  NO  → :lightweight
          NO  → Need quality scores?
                  YES → :ultralight_quality
                  NO  → :ultralight (maximum scalability)
```

# Examples
```julia
# Build full graph (default)
graph = build_kmer_graph(records, 31)

# Build ultralight graph (maximum memory savings)
graph = build_kmer_graph(records, 31; memory_profile=:ultralight)

# Build quality-aware ultralight graph from FASTQ
graph = build_kmer_graph(fastq_records, 31; memory_profile=:ultralight_quality)

# Backward compatible
graph = build_kmer_graph(records, 31; lightweight=true)
```

# See Also
- `build_qualmer_graph`: Quality-aware version for FASTQ data
- `add_observations_to_graph!`: Add more observations to existing graph
"""
function build_kmer_graph(
        records::Vector{<:FASTX.Record},
        k::Int;
        dataset_id::String = "dataset_01",
        mode::Symbol = :singlestrand,
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna,
        memory_profile::Symbol = :full,
        lightweight::Union{Nothing, Bool} = nothing
)
    # Backward compatibility: lightweight kwarg maps to memory_profile
    if !isnothing(lightweight)
        if lightweight
            memory_profile = :lightweight
        else
            memory_profile = :full
        end
    end

    # Quality profiles require FASTQ input
    if memory_profile in (:ultralight_quality, :lightweight_quality)
        if isempty(records)
            throw(ArgumentError("Cannot build graph from empty record set"))
        end
        if !(records[1] isa FASTX.FASTQ.Record)
            error("Memory profile :$memory_profile requires FASTQ input (quality scores needed)")
        end
    end

    # Dispatch based on memory_profile
    if memory_profile != :full
        # All reduced profiles: build singlestrand, then convert if needed
        singlestrand_graph = if memory_profile == :ultralight
            build_kmer_graph_singlestrand_ultralight(
                records, k;
                dataset_id = dataset_id,
                type_hint = type_hint,
                ambiguous_action = ambiguous_action
            )
        elseif memory_profile == :lightweight
            build_kmer_graph_singlestrand_lightweight(
                records, k;
                dataset_id = dataset_id,
                type_hint = type_hint,
                ambiguous_action = ambiguous_action
            )
        elseif memory_profile == :ultralight_quality
            build_kmer_graph_singlestrand_ultralight_quality(
                records, k;
                dataset_id = dataset_id,
                type_hint = type_hint,
                ambiguous_action = ambiguous_action
            )
        elseif memory_profile == :lightweight_quality
            build_kmer_graph_singlestrand_lightweight_quality(
                records, k;
                dataset_id = dataset_id,
                type_hint = type_hint,
                ambiguous_action = ambiguous_action
            )
        else
            error("Invalid memory_profile: :$memory_profile")
        end

        if mode == :singlestrand
            return singlestrand_graph
        elseif mode == :doublestrand
            return convert_to_doublestrand(singlestrand_graph)
        elseif mode == :canonical
            return convert_to_canonical(singlestrand_graph)
        else
            error("Invalid mode: $mode. Must be :singlestrand, :doublestrand, or :canonical")
        end
    elseif memory_profile == :full
        if mode == :singlestrand
            return build_kmer_graph_singlestrand(
                records, k;
                dataset_id = dataset_id,
                type_hint = type_hint,
                ambiguous_action = ambiguous_action
            )
        elseif mode == :doublestrand
            return build_kmer_graph_doublestrand(
                records, k;
                dataset_id = dataset_id,
                type_hint = type_hint,
                ambiguous_action = ambiguous_action
            )
        elseif mode == :canonical
            singlestrand_graph = build_kmer_graph_singlestrand(
                records, k;
                dataset_id = dataset_id,
                type_hint = type_hint,
                ambiguous_action = ambiguous_action
            )
            return convert_to_canonical(singlestrand_graph)
        else
            error("Invalid mode: $mode. Must be :singlestrand, :doublestrand, or :canonical")
        end
    else
        error("Invalid memory_profile: :$memory_profile. Must be :full, :lightweight, :ultralight, :ultralight_quality, or :lightweight_quality")
    end
end

# Note: The core implementation functions (build_kmer_graph_singlestrand,
# build_kmer_graph_doublestrand, _build_kmer_graph_core, etc.) are defined
# in core/graph-construction.jl and are available through the parent module.

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    build_kmer_graph_from_file(filepath, k; dataset_id=nothing, mode=:singlestrand,
                               type_hint=nothing, ambiguous_action=:dna)

Build k-mer graph directly from a FASTA or FASTQ file.

Automatically handles compressed files (.gz, .bz2, .xz) and uses filename as
dataset_id if not specified. If the filepath ends with `.fna`, `.frn`, or `.faa`
and `type_hint` is not provided, the extension is used as the alphabet hint.

# Arguments
- `filepath::String`: Path to FASTA or FASTQ file (compressed or uncompressed)
- `k::Int`: K-mer size
- `dataset_id::String=nothing`: Dataset identifier (defaults to filename without extension)
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand, :doublestrand, or :canonical)
- `type_hint::Union{Nothing,Symbol}=nothing`: Optional alphabet hint (:DNA, :RNA, :AA)
- `ambiguous_action::Symbol=:dna`: Resolution for ambiguous alphabets (:dna, :rna, :aa, :error)

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
        dataset_id::Union{String, Nothing} = nothing,
        mode::Symbol = :singlestrand,
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna,
        memory_profile::Symbol = :full,
        lightweight::Union{Nothing, Bool} = nothing
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

    file_hint = _alphabet_hint_from_path(filepath)
    final_hint = isnothing(type_hint) ? file_hint : type_hint

    return build_kmer_graph(
        records,
        k;
        dataset_id = dataset_id,
        mode = mode,
        type_hint = final_hint,
        ambiguous_action = ambiguous_action,
        memory_profile = memory_profile,
        lightweight = lightweight
    )
end

"""
    build_kmer_graph_from_files(filepaths, k; mode=:singlestrand,
                                type_hint=nothing, ambiguous_action=:dna)

Build k-mer graph from multiple FASTA/FASTQ files.

Each file is treated as a separate dataset, using the filename as dataset_id.
Handles compressed files automatically. If any filepath ends with `.fna`, `.frn`,
or `.faa` and `type_hint` is not provided, the extension hint is applied for all
files; conflicting hints raise an error.

# Arguments
- `filepaths::Vector{String}`: List of FASTA/FASTQ files (compressed or uncompressed)
- `k::Int`: K-mer size
- `mode::Symbol=:singlestrand`: Graph mode (:singlestrand, :doublestrand, or :canonical)
- `type_hint::Union{Nothing,Symbol}=nothing`: Optional alphabet hint (:DNA, :RNA, :AA)
- `ambiguous_action::Symbol=:dna`: Resolution for ambiguous alphabets (:dna, :rna, :aa, :error)

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
        mode::Symbol = :singlestrand,
        type_hint::Union{Nothing, Symbol} = nothing,
        ambiguous_action::Symbol = :dna,
        memory_profile::Symbol = :full,
        lightweight::Union{Nothing, Bool} = nothing
)
    if isempty(filepaths)
        error("No files provided")
    end

    # Backward compatibility
    if !isnothing(lightweight)
        memory_profile = lightweight ? :lightweight : :full
    end

    if !(mode in (:singlestrand, :doublestrand, :canonical))
        error("Invalid mode: $mode. Must be :singlestrand, :doublestrand, or :canonical")
    end

    file_hints = Set{Symbol}()
    for filepath in filepaths
        hint = _alphabet_hint_from_path(filepath)
        if !isnothing(hint)
            push!(file_hints, hint)
        end
    end

    if length(file_hints) > 1
        throw(ArgumentError("Conflicting alphabet hints from file extensions: $(sort!(collect(file_hints)))."))
    end

    inferred_hint = isempty(file_hints) ? nothing : first(file_hints)
    final_hint = isnothing(type_hint) ? inferred_hint : type_hint

    if !isnothing(type_hint) && !isnothing(inferred_hint)
        normalized_hint = Symbol(uppercase(String(type_hint)))
        if normalized_hint != inferred_hint
            throw(ArgumentError("Provided type_hint $type_hint conflicts with file extension hint $inferred_hint."))
        end
    end

    # Build graph from first file (strand-specific) and convert later if needed
    graph = build_kmer_graph_from_file(
        filepaths[1],
        k;
        mode = :singlestrand,
        type_hint = final_hint,
        ambiguous_action = ambiguous_action,
        memory_profile = memory_profile
    )

    # Get Mycelia module for open_fastx
    Mycelia_module = parentmodule(Rhizomorph)

    # Add observations from remaining files
    for filepath in filepaths[2:end]
        dataset_id = get_dataset_id_from_file(filepath)

        # Use open_fastx which handles compression and format detection
        records = collect(Mycelia_module.open_fastx(filepath))

        # Reduced mode graphs don't support add_observations_to_graph!
        # so we build a temporary graph and merge.
        if memory_profile != :full
            file_hint = _alphabet_hint_from_path(filepath)
            file_final_hint = isnothing(type_hint) ?
                              (isnothing(file_hint) ? final_hint : file_hint) : type_hint
            temp_graph = build_kmer_graph(
                records, k;
                dataset_id = dataset_id,
                type_hint = file_final_hint,
                ambiguous_action = ambiguous_action,
                memory_profile = memory_profile
            )
            # Merge temp_graph into main graph
            _merge_reduced_graphs!(graph, temp_graph, memory_profile)
        else
            add_observations_to_graph!(graph, records, k; dataset_id = dataset_id)
        end
    end

    # Determine if canonical/doublestrand conversion is valid for the k-mer type
    if !isempty(MetaGraphsNext.labels(graph))
        first_label = first(MetaGraphsNext.labels(graph))
        if mode != :singlestrand &&
           !(first_label isa Kmers.DNAKmer || first_label isa Kmers.RNAKmer)
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

"""
    _merge_reduced_graphs!(main_graph, temp_graph, memory_profile)

Merge a temporary reduced-mode graph into a main graph. Handles vertex/edge
count merging for all reduced memory profiles.
"""
function _merge_reduced_graphs!(main_graph, temp_graph, memory_profile::Symbol)
    for label in MetaGraphsNext.labels(temp_graph)
        temp_vdata = temp_graph[label]
        if !haskey(main_graph, label)
            # Create a new vertex with the same type
            main_graph[label] = _create_reduced_vertex(label, memory_profile)
        end
        vdata = main_graph[label]
        vdata.total_count += temp_vdata.total_count
        for (ds_id, ds_count) in temp_vdata.dataset_counts
            vdata.dataset_counts[ds_id] = get(vdata.dataset_counts, ds_id, 0) + ds_count
        end
        # Merge observation IDs if the type has them
        if hasfield(typeof(vdata), :dataset_observations)
            for (ds_id, obs_set) in temp_vdata.dataset_observations
                if !haskey(vdata.dataset_observations, ds_id)
                    vdata.dataset_observations[ds_id] = Set{String}()
                end
                union!(vdata.dataset_observations[ds_id], obs_set)
            end
        end
        # Merge quality scores if the type has them
        if hasfield(typeof(vdata), :joint_quality)
            if !isempty(temp_vdata.joint_quality)
                if isempty(vdata.joint_quality)
                    append!(vdata.joint_quality, temp_vdata.joint_quality)
                else
                    for i in eachindex(temp_vdata.joint_quality)
                        vdata.joint_quality[i] = UInt8(clamp(
                            Int(vdata.joint_quality[i]) + Int(temp_vdata.joint_quality[i]),
                            0, 255))
                    end
                end
            end
            if hasfield(typeof(vdata), :dataset_joint_quality)
                for (ds_id, ds_qual) in temp_vdata.dataset_joint_quality
                    existing = get(vdata.dataset_joint_quality, ds_id, nothing)
                    if isnothing(existing)
                        vdata.dataset_joint_quality[ds_id] = copy(ds_qual)
                    else
                        for i in eachindex(ds_qual)
                            existing[i] = UInt8(clamp(
                                Int(existing[i]) + Int(ds_qual[i]), 0, 255))
                        end
                    end
                end
            end
        end
    end
    for (src, dst) in MetaGraphsNext.edge_labels(temp_graph)
        temp_edata = temp_graph[src, dst]
        if !haskey(main_graph, src, dst)
            main_graph[src, dst] = _create_reduced_edge(
                temp_edata.overlap_length, memory_profile)
        end
        edata = main_graph[src, dst]
        edata.total_count += temp_edata.total_count
        for (ds_id, ds_count) in temp_edata.dataset_counts
            edata.dataset_counts[ds_id] = get(edata.dataset_counts, ds_id, 0) + ds_count
        end
        if hasfield(typeof(edata), :dataset_observations)
            for (ds_id, obs_set) in temp_edata.dataset_observations
                if !haskey(edata.dataset_observations, ds_id)
                    edata.dataset_observations[ds_id] = Set{String}()
                end
                union!(edata.dataset_observations[ds_id], obs_set)
            end
        end
        if hasfield(typeof(edata), :from_joint_quality)
            _merge_quality_vector!(edata.from_joint_quality, temp_edata.from_joint_quality)
            _merge_quality_vector!(edata.to_joint_quality, temp_edata.to_joint_quality)
        end
    end
end

function _merge_quality_vector!(dest::Vector{UInt8}, src::Vector{UInt8})
    if !isempty(src)
        if isempty(dest)
            append!(dest, src)
        else
            for i in eachindex(src)
                dest[i] = UInt8(clamp(Int(dest[i]) + Int(src[i]), 0, 255))
            end
        end
    end
end

function _create_reduced_vertex(kmer, memory_profile::Symbol)
    if memory_profile == :ultralight
        return UltralightKmerVertexData(kmer)
    elseif memory_profile == :lightweight
        return LightweightKmerVertexData(kmer)
    elseif memory_profile == :ultralight_quality
        return UltralightQualityKmerVertexData(kmer)
    elseif memory_profile == :lightweight_quality
        return LightweightQualityKmerVertexData(kmer)
    else
        error("Unknown memory_profile: $memory_profile")
    end
end

function _create_reduced_edge(overlap_length::Int, memory_profile::Symbol)
    if memory_profile == :ultralight
        return UltralightEdgeData(overlap_length)
    elseif memory_profile == :lightweight
        return LightweightEdgeData(overlap_length)
    elseif memory_profile == :ultralight_quality
        return UltralightQualityEdgeData(overlap_length)
    elseif memory_profile == :lightweight_quality
        return LightweightQualityEdgeData(overlap_length)
    else
        error("Unknown memory_profile: $memory_profile")
    end
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
        # Use get_all_dataset_ids for compatibility with both full and lightweight types
        for ds_id in get_all_dataset_ids(vertex_data)
            push!(all_dataset_ids, ds_id)
        end
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
        dataset_id::Union{String, Nothing} = nothing
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

"""
    find_low_coverage_kmers(graph, max_coverage::Int; dataset_id=nothing)

Find k-mers with coverage <= max_coverage.

# Arguments
- `graph`: K-mer graph
- `max_coverage::Int`: Maximum number of observations
- `dataset_id::Union{String,Nothing}=nothing`: Specific dataset or all datasets if nothing

# Returns
- `Vector`: K-mers meeting coverage threshold
"""
function find_low_coverage_kmers(
        graph::MetaGraphsNext.MetaGraph,
        max_coverage::Int;
        dataset_id::Union{String, Nothing} = nothing
)
    result = []

    for kmer in MetaGraphsNext.labels(graph)
        vertex_data = graph[kmer]

        coverage = if isnothing(dataset_id)
            count_total_observations(vertex_data)
        else
            count_dataset_observations(vertex_data, dataset_id)
        end

        if coverage <= max_coverage
            push!(result, kmer)
        end
    end

    return result
end
