"""
    pcoa_coordinates(
        distance_matrix::AbstractMatrix{<:Real};
        maxoutdim::Int = 2,
        samples_as_rows::Bool = true
    )

Compute Principal Coordinates Analysis (PCoA) coordinates from a precomputed
distance matrix using `pcoa_from_dist`.

# Arguments
- `distance_matrix`: Square distance matrix
- `maxoutdim`: Number of output dimensions to compute
- `samples_as_rows`: When `true` (default), return an `n_samples × maxoutdim`
  matrix suitable for plotting/data-frame workflows. When `false`, return the
  native `maxoutdim × n_samples` layout from `pcoa_from_dist`.

# Returns
A NamedTuple with fields:
- `model`: the fitted `MultivariateStats.MDS` model
- `coordinates`: the embedded coordinates in the requested orientation
"""
function pcoa_coordinates(
    distance_matrix::AbstractMatrix{<:Real};
    maxoutdim::Int = 2,
    samples_as_rows::Bool = true
)
    pcoa_result = pcoa_from_dist(distance_matrix; maxoutdim = maxoutdim)
    coordinates = samples_as_rows ? Matrix(transpose(pcoa_result.coordinates)) : pcoa_result.coordinates
    return (
        model = pcoa_result.model,
        coordinates = coordinates
    )
end

function _distance_matrix_sums_of_squares(
    distance_matrix::AbstractMatrix{<:Real},
    groups::AbstractVector,
    group_levels::AbstractVector
)
    ss_total = sum(abs2, distance_matrix) / (2 * length(groups))
    ss_within = 0.0

    for group in group_levels
        group_idx = findall(isequal(group), groups)
        n_group = length(group_idx)
        n_group <= 1 && continue

        pairwise_ss = 0.0
        for i in 1:(n_group - 1)
            idx_i = group_idx[i]
            for j in (i + 1):n_group
                idx_j = group_idx[j]
                pairwise_ss += abs2(distance_matrix[idx_i, idx_j])
            end
        end
        ss_within += pairwise_ss / n_group
    end

    return (
        ss_total = ss_total,
        ss_within = ss_within,
        ss_between = ss_total - ss_within
    )
end

"""
    permanova(
        distance_matrix::AbstractMatrix{<:Real},
        groups::AbstractVector;
        n_perm::Integer = 9999,
        rng = Random.default_rng()
    )

Run a one-factor PERMANOVA on a precomputed distance matrix.

# Arguments
- `distance_matrix`: Square symmetric distance matrix
- `groups`: Group assignment for each sample in the matrix order. Accepts any
  string-like or categorical vector type.
- `n_perm`: Number of label permutations
- `rng`: Random number generator used for permutations

# Returns
A NamedTuple with fields:
- `F`
- `p`
- `R2`
- `n_perm`
- `df_between`
- `df_within`
- `ss_between`
- `ss_within`
- `ss_total`
- `groups`
- `group_sizes`
"""
function permanova(
    distance_matrix::AbstractMatrix{<:Real},
    groups::AbstractVector;
    n_perm::Integer = 9999,
    rng = Random.default_rng()
)
    n_samples = size(distance_matrix, 1)
    n_samples == size(distance_matrix, 2) ||
        throw(ArgumentError("distance_matrix must be square, got size $(size(distance_matrix))"))
    length(groups) == n_samples ||
        throw(ArgumentError("length(groups) ($(length(groups))) must equal size(distance_matrix, 1) ($n_samples)"))
    n_perm >= 1 || throw(ArgumentError("n_perm must be at least 1, got $n_perm"))
    all(isfinite, distance_matrix) ||
        throw(ArgumentError("distance_matrix must contain only finite values"))
    LinearAlgebra.issymmetric(distance_matrix) ||
        throw(ArgumentError("distance_matrix must be symmetric"))
    any(ismissing, groups) &&
        throw(ArgumentError("groups must not contain missing values"))

    group_vector = collect(groups)
    group_levels = unique(group_vector)
    length(group_levels) >= 2 ||
        throw(ArgumentError("permanova requires at least 2 groups, got $(length(group_levels))"))
    n_samples > length(group_levels) ||
        throw(ArgumentError("permanova requires more samples than groups"))

    observed_ss = _distance_matrix_sums_of_squares(distance_matrix, group_vector, group_levels)
    observed_ss.ss_total > 0 ||
        throw(ArgumentError("distance_matrix has zero total variation; PERMANOVA is undefined"))

    df_between = length(group_levels) - 1
    df_within = n_samples - length(group_levels)
    f_observed = (observed_ss.ss_between / df_between) / (observed_ss.ss_within / df_within)

    permuted_groups = copy(group_vector)
    n_greater = 0
    for _ in 1:n_perm
        Random.shuffle!(rng, permuted_groups)
        permuted_ss = _distance_matrix_sums_of_squares(distance_matrix, permuted_groups, group_levels)
        f_permuted = (permuted_ss.ss_between / df_between) / (permuted_ss.ss_within / df_within)
        if f_permuted >= f_observed
            n_greater += 1
        end
    end

    return (
        F = f_observed,
        p = (n_greater + 1) / (n_perm + 1),
        R2 = observed_ss.ss_between / observed_ss.ss_total,
        n_perm = Int(n_perm),
        df_between = df_between,
        df_within = df_within,
        ss_between = observed_ss.ss_between,
        ss_within = observed_ss.ss_within,
        ss_total = observed_ss.ss_total,
        groups = collect(group_levels),
        group_sizes = [count(isequal(group), group_vector) for group in group_levels]
    )
end

"""
    adjusted_rand_index(labels_a, labels_b)

Compute the Adjusted Rand Index (ARI) between two cluster label assignments.

# Arguments
- `labels_a`: First cluster/label assignment vector
- `labels_b`: Second cluster/label assignment vector

# Returns
Adjusted Rand Index as a `Float64`, where:
- `1.0` indicates perfect agreement
- `0.0` indicates chance-level agreement
- negative values indicate worse-than-chance agreement
"""
function adjusted_rand_index(labels_a, labels_b)
    n = length(labels_a)
    n == length(labels_b) ||
        throw(ArgumentError("labels_a and labels_b must have the same length"))
    n >= 2 || throw(ArgumentError("adjusted_rand_index requires at least 2 observations"))

    contingency_map = Dict{Tuple{Any, Any}, Int}()
    for (a, b) in zip(labels_a, labels_b)
        key = (a, b)
        contingency_map[key] = get(contingency_map, key, 0) + 1
    end

    a_counts = StatsBase.countmap(labels_a)
    b_counts = StatsBase.countmap(labels_b)
    comb2(x) = x * (x - 1) / 2

    sum_nij = sum(comb2(v) for v in values(contingency_map))
    sum_ai = sum(comb2(v) for v in values(a_counts))
    sum_bj = sum(comb2(v) for v in values(b_counts))
    expected = sum_ai * sum_bj / comb2(n)
    max_index = 0.5 * (sum_ai + sum_bj)
    denominator = max_index - expected

    denominator == 0.0 && return 0.0
    return (sum_nij - expected) / denominator
end
