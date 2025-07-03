# src/distance_metrics.jl
"""
Compute the Jaccard distance between columns of a binary matrix.

# Arguments
- `M::AbstractMatrix{<:Integer}`: Binary matrix where rows are features and columns are samples

# Returns
- `Matrix{Float64}`: Symmetric distance matrix with Jaccard distances
"""
function jaccard_distance(M::AbstractMatrix{<:Integer})
    n_samples = size(M, 2)
    D = zeros(Float64, n_samples, n_samples)
    for i in 1:n_samples
        for j in i+1:n_samples
            intersection = sum(M[:, i] .& M[:, j])
            union = sum(M[:, i] .| M[:, j])
            D[i, j] = 1 - intersection / union
            D[j, i] = D[i, j]
        end
    end
    return D
end

"""
Compute the Bray-Curtis distance between columns of a count matrix.

# Arguments
- `M::AbstractMatrix{<:Integer}`: Count matrix where rows are features and columns are samples

# Returns
- `Matrix{Float64}`: Symmetric distance matrix with Bray-Curtis distances
"""
function bray_curtis_distance(M::AbstractMatrix{<:Integer})
    n_samples = size(M, 2)
    D = zeros(Float64, n_samples, n_samples)
    for i in 1:n_samples
        for j in i+1:n_samples
            sum_abs_diff = sum(abs.(M[:, i] - M[:, j]))
            sum_total = sum(M[:, i]) + sum(M[:, j])
            D[i, j] = sum_abs_diff / sum_total
            D[j, i] = D[i, j]
        end
    end
    return D
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the cosine similarity between two k-mer count dictionaries.

# Arguments
- `kmer_counts_1::Dict{String,Int}`: First dictionary mapping k-mer sequences to their counts
- `kmer_counts_2::Dict{String,Int}`: Second dictionary mapping k-mer sequences to their counts

# Returns
- `Float64`: Cosine distance between the two k-mer count vectors, in range [0,1]
  where 0 indicates identical distributions and 1 indicates maximum dissimilarity

# Details
Converts k-mer count dictionaries into vectors using a unified set of keys,
then computes cosine distance. Missing k-mers are treated as count 0.
Result is invariant to input order and total counts (normalized internally).
"""
function kmer_counts_to_cosine_similarity(kmer_counts_1, kmer_counts_2)
    sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
    a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    # Distances.cosine_dist(a, b) == Distances.cosine_dist(b, a) == Distances.cosine_dist(a ./ sum(a), b ./ sum(b)) == Distances.cosine_dist(b ./ sum(b), a ./ sum(a))
    return Distances.cosine_dist(a, b)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Jensen-Shannon divergence between two k-mer frequency distributions.

# Arguments
- `kmer_counts_1`: Dictionary mapping k-mers to their counts in first sequence
- `kmer_counts_2`: Dictionary mapping k-mers to their counts in second sequence

# Returns
- Normalized Jensen-Shannon divergence score between 0 and 1, where:
  - 0 indicates identical distributions
  - 1 indicates maximally different distributions

# Notes
- The measure is symmetric: JS(P||Q) = JS(Q||P)
- Counts are automatically normalized to probability distributions
"""
function kmer_counts_to_js_divergence(kmer_counts_1, kmer_counts_2)
    sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
    a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
    a_norm = a ./ sum(a)
    b_norm = b ./ sum(b)
    # Distances.js_divergence(a ./ sum(a), b ./ sum(b)) == Distances.js_divergence(b ./ sum(b), a ./ sum(a))
    # Distances.js_divergence(a, b) != Distances.js_divergence(a ./ sum(a), b ./ sum(b))
    return Distances.js_divergence(a_norm, b_norm)
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the Jaccard similarity coefficient between two sets.

The Jaccard similarity is defined as the size of the intersection divided by the size
of the union of two sets:

    J(A,B) = |A ∩ B| / |A ∪ B|

# Arguments
- `set1`: First set for comparison
- `set2`: Second set for comparison

# Returns
- `Float64`: A value between 0.0 and 1.0, where:
  * 1.0 indicates identical sets
  * 0.0 indicates completely disjoint sets
"""
function jaccard_similarity(set1, set2)
    return length(intersect(set1, set2)) / length(union(set1, set2))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate the Jaccard distance between two sets, which is the complement of the Jaccard similarity.

The Jaccard distance is defined as:
``J_d(A,B) = 1 - J_s(A,B) = 1 - \\frac{|A ∩ B|}{|A ∪ B|}``

# Arguments
- `set1`: First set to compare
- `set2`: Second set to compare

# Returns
- `Float64`: A value in [0,1] where 0 indicates identical sets and 1 indicates disjoint sets
"""
function jaccard_distance(set1, set2)
    return 1.0 - jaccard_similarity(set1, set2)
end

# function kmer_counts_to_jaccard(kmer_counts_1::AbstractDict{Kmers.DNAKmer{k}, Int64}, kmer_counts_2::AbstractDict{Kmers.DNAKmer{k}, Int64}) where k
#     # sorted_shared_keys = sort(collect(union(keys(kmer_counts_1), keys(kmer_counts_2))))
#     # a = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
#     # b = [get(kmer_counts_1, kmer, 0) for kmer in sorted_shared_keys]
#     # a_indices = findall(a .> 0)
#     # b_indices = findall(b .> 0)
#     # return Distances.jaccard(a_indices, b_indices)
#     return jaccard(collect(keys(kmer_counts_1)), collect(keys(kmer_counts_2)))
# end