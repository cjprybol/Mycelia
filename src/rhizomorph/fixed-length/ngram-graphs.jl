# N-gram Graph Construction
#
# Fixed-length n-gram de Bruijn graphs for Unicode text/string analysis.
#
# This module provides graph construction functions for analyzing fixed-length
# character patterns (n-grams) in text data.
#
# Key Features:
# - Works with full Unicode strings
# - No reverse complement (strings don't have biological RC)
# - Evidence tracking with nested Dict structure
# - Support for text analysis, NLP, and pattern detection
#
# Use Cases:
# - Natural language processing (word/character pattern analysis)
# - Text compression and pattern detection
# - Linguistic analysis
# - Code analysis and pattern matching
#
# Based on rhizomorph-graph-ecosystem-plan.md Section 2.1

# ============================================================================
# Public API - N-gram Graph Construction
# ============================================================================

"""
    build_ngram_graph(strings, n; dataset_id="dataset_01", memory_profile=:full,
                      lightweight=nothing)

Build an n-gram de Bruijn graph from Unicode strings.

This is the main user-facing function for n-gram graph construction.
N-grams are fixed-length character sequences with (n-1) character overlap.

# Arguments
- `strings::Vector{String}`: Input strings to analyze
- `n::Int`: N-gram size (number of characters)
- `dataset_id::String="dataset_01"`: Dataset identifier for evidence tracking
- `memory_profile::Symbol=:full`: Memory profile (:full, :lightweight, :ultralight).
  Quality profiles are not supported for n-gram graphs (text has no quality scores).
- `lightweight::Union{Nothing,Bool}=nothing`: Deprecated. Use `memory_profile` instead.

# Returns
- `MetaGraphsNext.MetaGraph`: N-gram de Bruijn graph with evidence

# Examples
```julia
# Analyze text patterns
texts = ["hello world", "world peace", "hello there"]
graph = build_ngram_graph(texts, 3)

# Ultralight mode for maximum scalability
graph = build_ngram_graph(texts, 3; memory_profile=:ultralight)

# Backward compatible
graph = build_ngram_graph(texts, 3; lightweight=true)
```

# See Also
- `build_ngram_graph_from_file`: Build from text file
- `add_observations_to_ngram_graph!`: Add more strings to existing graph
- `get_ngram_statistics`: Get graph statistics
"""
function build_ngram_graph(
        strings::Vector{String},
        n::Int;
        dataset_id::String = "dataset_01",
        memory_profile::Symbol = :full,
        lightweight::Union{Nothing, Bool} = nothing
)
    # Backward compatibility
    if !isnothing(lightweight)
        memory_profile = lightweight ? :lightweight : :full
    end

    # Quality profiles don't apply to text n-grams
    if memory_profile in (:ultralight_quality, :lightweight_quality)
        error("Quality memory profiles are not supported for n-gram graphs (text has no quality scores)")
    end

    if memory_profile == :ultralight
        return build_ngram_graph_singlestrand_ultralight(strings, n; dataset_id = dataset_id)
    elseif memory_profile == :lightweight
        return build_ngram_graph_singlestrand_lightweight(strings, n; dataset_id = dataset_id)
    elseif memory_profile == :full
        return build_ngram_graph_singlestrand(strings, n; dataset_id = dataset_id)
    else
        error("Invalid memory_profile: :$memory_profile. Must be :full, :lightweight, or :ultralight")
    end
end

# Note: The core implementation (build_ngram_graph_singlestrand,
# add_observations_to_ngram_graph!) is defined in core/graph-construction.jl
# and is available through the parent module.

# ============================================================================
# Convenience Functions
# ============================================================================

"""
    build_ngram_graph_from_file(filepath, n; dataset_id=nothing)

Build n-gram graph from a text file.

Reads all lines from the file and builds an n-gram graph.
Each line is treated as a separate string observation.

# Arguments
- `filepath::String`: Path to text file
- `n::Int`: N-gram size
- `dataset_id::String=nothing`: Dataset identifier (defaults to filename)

# Returns
- `MetaGraphsNext.MetaGraph`: N-gram de Bruijn graph

# Examples
```julia
# Analyze text file
graph = build_ngram_graph_from_file("novel.txt", 5)

# Analyze code file
graph = build_ngram_graph_from_file("source.jl", 10)
```
"""
function build_ngram_graph_from_file(
        filepath::String,
        n::Int;
        dataset_id::Union{String, Nothing} = nothing,
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

    # Read all lines from file
    strings = readlines(filepath)

    # Filter out empty lines
    strings = filter(!isempty, strings)

    if isempty(strings)
        error("No non-empty lines found in file: $filepath")
    end

    return build_ngram_graph(strings, n; dataset_id = dataset_id,
        memory_profile = memory_profile, lightweight = lightweight)
end

"""
    build_ngram_graph_from_files(filepaths, n)

Build n-gram graph from multiple text files.

Each file is treated as a separate dataset, using the filename as dataset_id.

# Arguments
- `filepaths::Vector{String}`: List of text files
- `n::Int`: N-gram size

# Returns
- `MetaGraphsNext.MetaGraph`: N-gram de Bruijn graph with evidence from all files

# Examples
```julia
# Compare multiple documents
files = ["doc1.txt", "doc2.txt", "doc3.txt"]
graph = build_ngram_graph_from_files(files, 4)

# See which n-grams appear in which documents
for ngram_label in MetaGraphsNext.labels(graph)
    vertex = graph[ngram_label]
    datasets = get_all_dataset_ids(vertex)
    if length(datasets) > 1
        println("N-gram '\$ngram_label' appears in: \$(join(datasets, ", "))")
    end
end
```
"""
function build_ngram_graph_from_files(
        filepaths::Vector{String},
        n::Int;
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

    # Build graph from first file
    graph = build_ngram_graph_from_file(filepaths[1], n; memory_profile = memory_profile)

    # Add observations from remaining files
    for filepath in filepaths[2:end]
        dataset_id = get_dataset_id_from_file(filepath)

        # Read lines from file
        strings = readlines(filepath)
        strings = filter(!isempty, strings)

        # Add observations using appropriate function for graph type
        if memory_profile == :lightweight
            add_observations_to_ngram_graph_lightweight!(
                graph, strings, n; dataset_id = dataset_id)
        elseif memory_profile == :ultralight
            add_observations_to_ngram_graph_ultralight!(
                graph, strings, n; dataset_id = dataset_id)
        else
            add_observations_to_ngram_graph!(graph, strings, n; dataset_id = dataset_id)
        end
    end

    return graph
end

# ============================================================================
# N-gram Specific Analysis Functions
# ============================================================================

"""
    get_ngram_statistics(graph)

Get basic statistics about an n-gram graph.

# Returns
Dictionary with:
- `:num_vertices`: Number of unique n-grams
- `:num_edges`: Number of n-gram transitions
- `:num_datasets`: Number of datasets with evidence
- `:total_observations`: Total number of string observations
- `:n`: N-gram size
- `:most_common_ngram`: Most frequently observed n-gram

# Examples
```julia
stats = get_ngram_statistics(graph)
println("Found \$(stats[:num_vertices]) unique n-grams")
println("Most common: \$(stats[:most_common_ngram])")
```
"""
function get_ngram_statistics(graph::MetaGraphsNext.MetaGraph)
    if isempty(MetaGraphsNext.labels(graph))
        return Dict(
            :num_vertices => 0,
            :num_edges => 0,
            :num_datasets => 0,
            :total_observations => 0,
            :n => nothing,
            :most_common_ngram => nothing
        )
    end

    # Get n-gram size from first label
    first_ngram = first(MetaGraphsNext.labels(graph))
    n = length(first_ngram)

    # Collect dataset IDs and count observations
    all_dataset_ids = Set{String}()
    total_observations = 0
    max_obs_count = 0
    most_common = nothing

    for ngram in MetaGraphsNext.labels(graph)
        vertex_data = graph[ngram]
        # Use get_all_dataset_ids for compatibility with both full and lightweight types
        for ds_id in get_all_dataset_ids(vertex_data)
            push!(all_dataset_ids, ds_id)
        end
        obs_count = count_total_observations(vertex_data)
        total_observations += obs_count

        if obs_count > max_obs_count
            max_obs_count = obs_count
            most_common = ngram
        end
    end

    return Dict(
        :num_vertices => Graphs.nv(graph.graph),
        :num_edges => Graphs.ne(graph.graph),
        :num_datasets => length(all_dataset_ids),
        :total_observations => total_observations,
        :n => n,
        :most_common_ngram => most_common
    )
end

"""
    find_high_coverage_ngrams(graph, min_coverage::Int; dataset_id=nothing)

Find n-grams that appear at least `min_coverage` times.

# Arguments
- `graph`: N-gram graph
- `min_coverage::Int`: Minimum number of observations
- `dataset_id::Union{String,Nothing}=nothing`: Specific dataset or all datasets if nothing

# Returns
- `Vector{String}`: N-grams meeting coverage threshold

# Examples
```julia
# Find n-grams that appear at least 5 times
common = find_high_coverage_ngrams(graph, 5)

# Find common patterns in specific dataset
common_doc1 = find_high_coverage_ngrams(graph, 3; dataset_id="doc1")
```
"""
function find_high_coverage_ngrams(
        graph::MetaGraphsNext.MetaGraph,
        min_coverage::Int;
        dataset_id::Union{String, Nothing} = nothing
)
    result = String[]

    for ngram in MetaGraphsNext.labels(graph)
        vertex_data = graph[ngram]

        coverage = if isnothing(dataset_id)
            count_total_observations(vertex_data)
        else
            count_dataset_observations(vertex_data, dataset_id)
        end

        if coverage >= min_coverage
            push!(result, ngram)
        end
    end

    return result
end

"""
    find_unique_ngrams(graph, dataset_id::String)

Find n-grams that appear only in a specific dataset.

Useful for identifying dataset-specific patterns.

# Arguments
- `graph`: N-gram graph with multiple datasets
- `dataset_id::String`: Dataset to find unique patterns for

# Returns
- `Vector{String}`: N-grams unique to the specified dataset

# Examples
```julia
# Find patterns unique to doc1
unique_to_doc1 = find_unique_ngrams(graph, "doc1")
```
"""
function find_unique_ngrams(
        graph::MetaGraphsNext.MetaGraph,
        dataset_id::String
)
    unique_ngrams = String[]

    for ngram in MetaGraphsNext.labels(graph)
        vertex_data = graph[ngram]
        dataset_ids = get_all_dataset_ids(vertex_data)

        # Check if this ngram only appears in the specified dataset
        if length(dataset_ids) == 1 && dataset_ids[1] == dataset_id
            push!(unique_ngrams, ngram)
        end
    end

    return unique_ngrams
end

"""
    find_shared_ngrams(graph, dataset_ids::Vector{String})

Find n-grams that appear in all specified datasets.

Useful for finding common patterns across multiple texts.

# Arguments
- `graph`: N-gram graph with multiple datasets
- `dataset_ids::Vector{String}`: Datasets to check for shared patterns

# Returns
- `Vector{String}`: N-grams that appear in ALL specified datasets

# Examples
```julia
# Find patterns common to all three documents
shared = find_shared_ngrams(graph, ["doc1", "doc2", "doc3"])
```
"""
function find_shared_ngrams(
        graph::MetaGraphsNext.MetaGraph,
        dataset_ids::Vector{String}
)
    shared_ngrams = String[]
    target_set = Set(dataset_ids)

    for ngram in MetaGraphsNext.labels(graph)
        vertex_data = graph[ngram]
        ngram_datasets = Set(get_all_dataset_ids(vertex_data))

        # Check if all target datasets are present
        if target_set ⊆ ngram_datasets
            push!(shared_ngrams, ngram)
        end
    end

    return shared_ngrams
end
