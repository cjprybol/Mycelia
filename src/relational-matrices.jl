# src/relational-matrices.jl
#
# Relational matrix data structures and transformations for entity-entity measurements.
# Supports phenotypic clustering workflows where entities are characterized by their
# relationships with other entities (e.g., host range profiles, drug susceptibilities).

"""
    RelationalMatrix{T}

Container for a matrix derived from pairwise measurements between two entity types.

Generalizes entity-entity relationship data (e.g., phage-host interactions,
drug-pathogen susceptibilities, gene-condition expression profiles).

# Fields
- `matrix::Matrix{T}`: The measurement matrix (entity_a rows × entity_b columns)
- `entity_a_ids::Vector{String}`: Row identifiers (entity_a)
- `entity_b_ids::Vector{String}`: Column identifiers (entity_b)
- `entity_a_name::String`: Label for entity_a type (e.g., "sample", "drug")
- `entity_b_name::String`: Label for entity_b type (e.g., "feature", "target")
- `measurement_name::String`: What the values represent (e.g., "activity", "expression")
- `metadata_a::Union{DataFrames.DataFrame, Nothing}`: Optional entity_a metadata
- `metadata_b::Union{DataFrames.DataFrame, Nothing}`: Optional entity_b metadata

# Example
```julia
import Mycelia
import DataFrames

df = DataFrames.DataFrame(
    sample = ["S1", "S1", "S2", "S2"],
    feature = ["F1", "F2", "F1", "F2"],
    value = [0.8, 0.2, 0.1, 0.9]
)

rm = Mycelia.long_to_relational_matrix(
    df, :sample, :feature, :value;
    entity_a_name="sample",
    entity_b_name="feature",
    measurement_name="activity"
)
```
"""
struct RelationalMatrix{T}
    matrix::Matrix{T}
    entity_a_ids::Vector{String}
    entity_b_ids::Vector{String}
    entity_a_name::String
    entity_b_name::String
    measurement_name::String
    metadata_a::Union{DataFrames.DataFrame, Nothing}
    metadata_b::Union{DataFrames.DataFrame, Nothing}
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Construct a RelationalMatrix from a DataFrame in long format.

# Arguments
- `df::DataFrames.DataFrame`: Long-format DataFrame with entity_a, entity_b, value columns
- `entity_a_col::Symbol`: Column name for entity_a identifiers
- `entity_b_col::Symbol`: Column name for entity_b identifiers
- `value_col::Symbol`: Column name for measurement values
- `entity_a_name::String`: Human-readable name for entity_a type (default: "entity_a")
- `entity_b_name::String`: Human-readable name for entity_b type (default: "entity_b")
- `measurement_name::String`: Human-readable name for measurement (default: "value")
- `missing_value`: Value to use for missing combinations (default: NaN)
- `agg_func::Function`: Function to aggregate duplicate entries (default: Statistics.median)
- `entity_a_order::Union{Nothing, Vector{String}}`: Optional explicit ordering for entity_a
- `entity_b_order::Union{Nothing, Vector{String}}`: Optional explicit ordering for entity_b

# Returns
- `RelationalMatrix{T}`: Constructed relational matrix

# Example
```julia
import Mycelia
import DataFrames

df = DataFrames.DataFrame(
    sample = ["S1", "S1", "S2", "S2"],
    feature = ["F1", "F2", "F1", "F2"],
    value = [0.8, 0.2, 0.1, 0.9]
)

rm = Mycelia.long_to_relational_matrix(df, :sample, :feature, :value)
```
"""
function long_to_relational_matrix(
        df::DataFrames.DataFrame,
        entity_a_col::Symbol,
        entity_b_col::Symbol,
        value_col::Symbol;
        entity_a_name::String = "entity_a",
        entity_b_name::String = "entity_b",
        measurement_name::String = "value",
        missing_value = NaN,
        agg_func::Function = Statistics.median,
        entity_a_order::Union{Nothing, Vector{String}} = nothing,
        entity_b_order::Union{Nothing, Vector{String}} = nothing
)
    ## Validate columns exist
    required_cols = [entity_a_col, entity_b_col, value_col]
    for col in required_cols
        if !(col in DataFrames.propertynames(df))
            throw(ArgumentError("Column $col not found in DataFrame. Available: $(DataFrames.propertynames(df))"))
        end
    end

    ## Aggregate duplicates if present
    grouped = DataFrames.combine(
        DataFrames.groupby(df, [entity_a_col, entity_b_col]),
        value_col => agg_func => :_agg_value
    )

    ## Determine entity ordering
    if isnothing(entity_a_order)
        entity_a_ids = sort(unique(string.(grouped[!, entity_a_col])))
    else
        entity_a_ids = entity_a_order
        ## Warn about entities in data but not in order
        data_a_ids = Set(string.(grouped[!, entity_a_col]))
        missing_from_order = setdiff(data_a_ids, Set(entity_a_order))
        if !isempty(missing_from_order)
            @warn "$(length(missing_from_order)) $(entity_a_col) values in data not in entity_a_order (will be dropped)"
        end
    end

    if isnothing(entity_b_order)
        entity_b_ids = sort(unique(string.(grouped[!, entity_b_col])))
    else
        entity_b_ids = entity_b_order
        data_b_ids = Set(string.(grouped[!, entity_b_col]))
        missing_from_order = setdiff(data_b_ids, Set(entity_b_order))
        if !isempty(missing_from_order)
            @warn "$(length(missing_from_order)) $(entity_b_col) values in data not in entity_b_order (will be dropped)"
        end
    end

    n_a = length(entity_a_ids)
    n_b = length(entity_b_ids)

    ## Create index lookups
    a_to_idx = Dict(a => i for (i, a) in enumerate(entity_a_ids))
    b_to_idx = Dict(b => i for (i, b) in enumerate(entity_b_ids))

    ## Initialize matrix with missing values
    T = promote_type(eltype(grouped._agg_value), typeof(missing_value))
    matrix = fill(convert(T, missing_value), n_a, n_b)

    ## Fill matrix
    n_filled = 0
    n_skipped = 0
    for row in DataFrames.eachrow(grouped)
        a_id = string(row[entity_a_col])
        b_id = string(row[entity_b_col])
        val = row._agg_value

        if !haskey(a_to_idx, a_id) || !haskey(b_to_idx, b_id)
            n_skipped += 1
            continue
        end

        i = a_to_idx[a_id]
        j = b_to_idx[b_id]
        matrix[i, j] = convert(T, val)
        n_filled += 1
    end

    ## Report statistics
    total_cells = n_a * n_b
    pct_filled = round(100 * n_filled / total_cells, digits = 1)
    @info "RelationalMatrix: $n_a $(entity_a_name) × $n_b $(entity_b_name), $pct_filled% filled"

    return RelationalMatrix{T}(
        matrix,
        entity_a_ids,
        entity_b_ids,
        entity_a_name,
        entity_b_name,
        measurement_name,
        nothing,
        nothing
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Transpose a RelationalMatrix, swapping entity_a and entity_b.

Useful for switching clustering perspective (e.g., from sample-centric to feature-centric).

# Arguments
- `rm::RelationalMatrix{T}`: Input relational matrix

# Returns
- `RelationalMatrix{T}`: Transposed matrix with swapped entity types
"""
function transpose_relational(rm::RelationalMatrix{T}) where {T}
    return RelationalMatrix{T}(
        permutedims(rm.matrix),
        rm.entity_b_ids,
        rm.entity_a_ids,
        rm.entity_b_name,
        rm.entity_a_name,
        rm.measurement_name,
        rm.metadata_b,
        rm.metadata_a
    )
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate summary statistics for a RelationalMatrix.

# Arguments
- `rm::RelationalMatrix{T}`: Input relational matrix
- `verbose::Bool`: Print summary to stdout (default: true)

# Returns
- `NamedTuple`: Summary statistics including shape, sparsity, and value distribution
"""
function summarize_relational_matrix(rm::RelationalMatrix{T}; verbose::Bool = true) where {T}
    n_a, n_b = size(rm.matrix)
    total_cells = n_a * n_b

    ## Count valid (non-missing, non-NaN) values
    valid_values = filter(x -> !ismissing(x) && !isnan(x), vec(rm.matrix))
    n_valid = length(valid_values)
    n_missing = total_cells - n_valid

    ## Per-row and per-column missing counts
    row_missing = [count(x -> ismissing(x) || isnan(x), rm.matrix[i, :]) for i in 1:n_a]
    col_missing = [count(x -> ismissing(x) || isnan(x), rm.matrix[:, j]) for j in 1:n_b]

    stats = (
        entity_a_name = rm.entity_a_name,
        entity_b_name = rm.entity_b_name,
        measurement_name = rm.measurement_name,
        n_entity_a = n_a,
        n_entity_b = n_b,
        n_total_cells = total_cells,
        n_valid_values = n_valid,
        n_missing_values = n_missing,
        pct_filled = round(100 * n_valid / total_cells, digits = 2),
        sparsity = round(1.0 - (n_valid / total_cells), digits = 4),
        rows_all_missing = count(==(n_b), row_missing),
        cols_all_missing = count(==(n_a), col_missing),
        value_min = isempty(valid_values) ? NaN : minimum(valid_values),
        value_max = isempty(valid_values) ? NaN : maximum(valid_values),
        value_mean = isempty(valid_values) ? NaN : Statistics.mean(valid_values),
        value_median = isempty(valid_values) ? NaN : Statistics.median(valid_values),
        value_std = isempty(valid_values) ? NaN : Statistics.std(valid_values)
    )

    if verbose
        println("=" ^ 50)
        println("RelationalMatrix Summary")
        println("=" ^ 50)
        println("  Rows ($(rm.entity_a_name)): $(stats.n_entity_a)")
        println("  Cols ($(rm.entity_b_name)): $(stats.n_entity_b)")
        println("  Values ($(rm.measurement_name)):")
        println("    Filled: $(stats.n_valid_values) / $(stats.n_total_cells) ($(stats.pct_filled)%)")
        println("    Range: [$(stats.value_min), $(stats.value_max)]")
        println("    Mean: $(round(stats.value_mean, digits=4))")
        println("    Std: $(round(stats.value_std, digits=4))")
        if stats.rows_all_missing > 0
            println("  WARNING: $(stats.rows_all_missing) rows with ALL missing values")
        end
        if stats.cols_all_missing > 0
            println("  WARNING: $(stats.cols_all_missing) columns with ALL missing values")
        end
        println("=" ^ 50)
    end

    return stats
end

## Accessor functions

"""Get entity_a identifiers."""
entity_a_ids(rm::RelationalMatrix) = rm.entity_a_ids

"""Get entity_b identifiers."""
entity_b_ids(rm::RelationalMatrix) = rm.entity_b_ids

"""Get the underlying measurement matrix."""
Base.Matrix(rm::RelationalMatrix) = rm.matrix

"""Get matrix dimensions (n_entity_a, n_entity_b)."""
Base.size(rm::RelationalMatrix) = size(rm.matrix)

"""Get matrix dimensions along a specific axis."""
Base.size(rm::RelationalMatrix, dim::Int) = size(rm.matrix, dim)

"""Number of entity_a (rows)."""
n_entity_a(rm::RelationalMatrix) = length(rm.entity_a_ids)

"""Number of entity_b (columns)."""
n_entity_b(rm::RelationalMatrix) = length(rm.entity_b_ids)

"""Count of missing/NaN values."""
function n_missing(rm::RelationalMatrix)
    count(x -> ismissing(x) || (x isa Number && isnan(x)), rm.matrix)
end

"""Count of valid (non-missing, non-NaN) values."""
n_filled(rm::RelationalMatrix) = length(rm.matrix) - n_missing(rm)

"""Fraction of cells that are filled (valid values)."""
coverage(rm::RelationalMatrix) = n_filled(rm) / length(rm.matrix)

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Index into the relational matrix by entity IDs.

# Arguments
- `rm::RelationalMatrix`: The relational matrix
- `a_id::String`: Entity A identifier
- `b_id::String`: Entity B identifier

# Returns
- Value at the specified position
"""
function Base.getindex(rm::RelationalMatrix, a_id::String, b_id::String)
    i = findfirst(==(a_id), rm.entity_a_ids)
    j = findfirst(==(b_id), rm.entity_b_ids)
    isnothing(i) && error("Entity A ID not found: $a_id")
    isnothing(j) && error("Entity B ID not found: $b_id")
    return rm.matrix[i, j]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Index into the relational matrix by numeric indices.

# Arguments
- `rm::RelationalMatrix`: The relational matrix
- `i::Int`: Row index
- `j::Int`: Column index

# Returns
- Value at the specified position
"""
function Base.getindex(rm::RelationalMatrix, i::Int, j::Int)
    return rm.matrix[i, j]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Pretty print for RelationalMatrix.
"""
function Base.show(io::IO, rm::RelationalMatrix{T}) where {T}
    n_a, n_b = size(rm)
    pct = round(100 * coverage(rm), digits = 1)
    print(io,
        "RelationalMatrix{$T}: $(n_a) $(rm.entity_a_name) × $(n_b) $(rm.entity_b_name), $(pct)% filled")
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert RelationalMatrix to a DataFrame in long format.

# Arguments
- `rm::RelationalMatrix{T}`: Input relational matrix
- `include_missing::Bool`: Include rows for missing values (default: false)

# Returns
- `DataFrames.DataFrame`: Long-format DataFrame with entity_a, entity_b, value columns
"""
function relational_matrix_to_long(rm::RelationalMatrix{T}; include_missing::Bool = false) where {T}
    ## Build column vectors
    a_col = String[]
    b_col = String[]
    v_col = T[]

    for (i, a_id) in enumerate(rm.entity_a_ids)
        for (j, b_id) in enumerate(rm.entity_b_ids)
            val = rm.matrix[i, j]
            is_missing = ismissing(val) || (val isa Number && isnan(val))

            if include_missing || !is_missing
                push!(a_col, a_id)
                push!(b_col, b_id)
                push!(v_col, val)
            end
        end
    end

    return DataFrames.DataFrame(
        Symbol(rm.entity_a_name) => a_col,
        Symbol(rm.entity_b_name) => b_col,
        Symbol(rm.measurement_name) => v_col
    )
end

# =============================================================================
# Distance Functions for RelationalMatrix
# =============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute Gower distance directly from a RelationalMatrix.

Rows are treated as samples, columns as features.

# Arguments
- `rm::RelationalMatrix`: Input relational matrix
- `kwargs...`: Additional arguments passed to gower_distance

# Returns
- `Matrix{Float64}`: Pairwise distance matrix between entity_a samples
"""
function gower_distance(rm::RelationalMatrix; kwargs...)
    return gower_distance(rm.matrix; kwargs...)
end

# =============================================================================
# Relational Clustering Pipeline
# =============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

End-to-end relational clustering pipeline.

Clusters entity_a based on their profiles across entity_b (or vice versa).

# Arguments
- `rm::RelationalMatrix{T}`: Input relational matrix
- `cluster_by::Symbol`: `:entity_a` or `:entity_b` (which dimension to cluster). Default: `:entity_a`
- `distance_metric::Symbol`: Distance metric (`:euclidean`, `:cosine`, `:jaccard`, `:gower`). Default: `:euclidean`
- `imputation::Symbol`: How to handle missing values (`:max`, `:max_observed`, `:median`, `:mean`). Default: `:max_observed`
- `ks::Union{Nothing, AbstractRange}`: K values to test (auto if nothing)
- `min_k::Union{Nothing, Int}`: Minimum k (used if ks is nothing)
- `max_k::Union{Nothing, Int}`: Maximum k (used if ks is nothing)
- `plot_backend::Symbol`: Visualization backend (`:cairomakie` or `:statsplots`). Default: `:cairomakie`

# Returns
- `NamedTuple` with:
  - `relational_matrix`: Original input
  - `distance_matrix`: Computed/imputed distance matrix
  - `hcl`: Hierarchical clustering result
  - `optimal_k`: Inferred optimal cluster count
  - `assignments`: Cluster assignments
  - `rankings`: Vector{ClusterRanking}
  - `rankings_df`: DataFrame of rankings
  - `summary_df`: DataFrame of cluster summaries
  - `figure`: Silhouette plot

# Example
```julia
result = Mycelia.relational_clustering_pipeline(
    rm;
    cluster_by=:entity_a,
    distance_metric=:euclidean,
    ks=2:10
)
```
"""
function relational_clustering_pipeline(
        rm::RelationalMatrix{T};
        cluster_by::Symbol = :entity_a,
        distance_metric::Symbol = :euclidean,
        imputation::Symbol = :max_observed,
        ks::Union{Nothing, AbstractRange} = nothing,
        min_k::Union{Nothing, Int} = nothing,
        max_k::Union{Nothing, Int} = nothing,
        plot_backend::Symbol = :cairomakie
) where {T}
    println("=" ^ 60)
    println("Relational Clustering Pipeline")
    println("=" ^ 60)

    ## Determine which dimension to cluster
    if cluster_by == :entity_a
        profile_matrix = rm.matrix
        entity_ids = rm.entity_a_ids
        entity_name = rm.entity_a_name
        feature_name = rm.entity_b_name
    elseif cluster_by == :entity_b
        profile_matrix = permutedims(rm.matrix)
        entity_ids = rm.entity_b_ids
        entity_name = rm.entity_b_name
        feature_name = rm.entity_a_name
    else
        throw(ArgumentError("cluster_by must be :entity_a or :entity_b, got: $cluster_by"))
    end

    n_samples, n_features = size(profile_matrix)
    println("\n[1/5] Clustering $n_samples $(entity_name) by $n_features $(feature_name) profiles")

    ## Step 2: Compute distance matrix based on metric
    println("\n[2/5] Computing $distance_metric distances...")

    ## Convert to column-major for pairwise_distance_matrix (samples as columns)
    profile_matrix_t = permutedims(profile_matrix)

    raw_distance = if distance_metric == :euclidean
        pairwise_distance_matrix(
            profile_matrix_t;
            dist_func = Distances.euclidean,
            show_progress = true,
            progress_desc = "Computing Euclidean distances"
        )
    elseif distance_metric == :cosine
        pairwise_distance_matrix(
            profile_matrix_t;
            dist_func = Distances.cosine_dist,
            show_progress = true,
            progress_desc = "Computing cosine distances"
        )
    elseif distance_metric == :jaccard
        ## Binarize the matrix
        binary_matrix = BitMatrix(profile_matrix_t .> 0)
        binary_matrix_to_jaccard_distance_matrix(binary_matrix)
    elseif distance_metric == :gower
        ## Gower expects samples as rows
        gower_distance(profile_matrix)
    else
        throw(ArgumentError("Unsupported distance metric: $distance_metric. Use :euclidean, :cosine, :jaccard, or :gower"))
    end

    ## Ensure exact symmetry (hclust requires issymmetric to be true)
    raw_distance = Matrix(LinearAlgebra.Symmetric(raw_distance))
    for i in 1:size(raw_distance, 1)
        raw_distance[i, i] = 0.0
    end

    ## Step 3: Validate and impute missing values
    println("\n[3/5] Validating and imputing distances...")
    validation = validate_distance_matrix(raw_distance)
    if validation.has_missing
        ## Convert symbol to ImputationMethod enum
        imputation_method = if imputation == :max
            IMPUTE_MAX
        elseif imputation == :max_observed
            IMPUTE_MAX_OBSERVED
        elseif imputation == :median
            IMPUTE_MEDIAN
        elseif imputation == :mean
            IMPUTE_MEAN
        else
            throw(ArgumentError("Unknown imputation method: $imputation. Use :max, :max_observed, :median, or :mean"))
        end
        println("  Found $(validation.n_missing_pairs) missing pairs, imputing with $imputation")
        distance_matrix = impute_distances(raw_distance; method = imputation_method)
    else
        distance_matrix = raw_distance
    end

    ## Normalize if needed (skip if max distance is 0 to avoid division by zero)
    max_dist = maximum(filter(x -> !isnan(x) && !isinf(x), vec(distance_matrix)))
    if max_dist > 0
        distance_matrix = normalize_distance_matrix(distance_matrix)
    end

    ## Ensure exact symmetry for clustering algorithms (hclust requires issymmetric to be true)
    ## Use LinearAlgebra.Symmetric to guarantee exact symmetry, then convert back to Matrix
    distance_matrix = Matrix(LinearAlgebra.Symmetric(distance_matrix))
    ## Ensure diagonal is exactly zero
    for i in 1:size(distance_matrix, 1)
        distance_matrix[i, i] = 0.0
    end

    final_validation = validate_distance_matrix(distance_matrix)
    if !final_validation.is_valid
        @warn "Distance matrix validation issues: $(final_validation.violations)"
    end

    ## Step 4: Identify optimal clusters
    println("\n[4/5] Identifying optimal cluster count...")
    cluster_result = identify_optimal_number_of_clusters(
        distance_matrix;
        ks = ks,
        min_k = min_k,
        max_k = max_k,
        plot_backend = plot_backend
    )

    optimal_k = cluster_result.optimal_number_of_clusters
    assignments = cluster_result.assignments

    ## Cluster size statistics
    sizes = values(StatsBase.countmap(assignments))
    println("  Optimal k: $optimal_k")
    println("  Cluster sizes: min=$(minimum(sizes)), max=$(maximum(sizes)), median=$(Statistics.median(collect(sizes)))")

    ## Step 5: Rank cluster members
    println("\n[5/5] Ranking cluster members...")
    rankings = rank_cluster_members(
        distance_matrix,
        assignments,
        entity_ids;
        metric_type = :distance
    )

    rankings_df = rankings_to_dataframe(rankings)
    summary_df = cluster_summary(rankings, distance_matrix)

    println("\n" * "=" ^ 60)
    println("Pipeline complete!")
    println("  - $(n_samples) $(entity_name) clustered into $(optimal_k) groups")
    println("  - Rankings DataFrame: $(DataFrames.nrow(rankings_df)) rows")
    println("  - Summary DataFrame: $(DataFrames.nrow(summary_df)) clusters")
    println("=" ^ 60)

    return (
        relational_matrix = rm,
        distance_matrix = distance_matrix,
        distance_matrix_raw = raw_distance,
        hcl = cluster_result.hcl,
        optimal_k = optimal_k,
        assignments = assignments,
        rankings = rankings,
        rankings_df = rankings_df,
        summary_df = summary_df,
        figure = cluster_result.figure,
        silhouette_scores = cluster_result.silhouette_scores,
        ks = cluster_result.ks
    )
end
