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

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

Normalize a distance matrix by dividing all elements by the maximum non-NaN value.

# Arguments
- `distance_matrix`: A matrix of distance values that may contain `NaN`, `nothing`, or `missing` values

# Returns
- Normalized distance matrix with values scaled to [0, 1] range

# Details
- Filters out `NaN`, `nothing`, and `missing` values when finding the maximum
- All elements are divided by the same maximum value to preserve relative distances
- If all values are NaN/nothing/missing, may return NaN values
"""
function normalize_distance_matrix(distance_matrix)
    max_non_nan_value = maximum(filter(x -> !isnan(x) && !isnothing(x) && !ismissing(x), vec(distance_matrix)))
    return distance_matrix ./ max_non_nan_value
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
# where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)
# """
# function count_matrix_to_probability_matrix(
#         counts_matrix,
#         probability_matrix_file = replace(counts_matrix_file, ".bin" => ".probability_matrix.bin")
#     )
#     probability_matrix = Mmap.mmap(probability_matrix_file, Array{Float64, 2}, size(counts_matrix))
#     if !isfile(probability_matrix_file)
#         println("creating new probability matrix $probability_matrix_file")
#         # probability_matrix .= count_matrix_to_probability_matrix(counts_matrix)
#         for (i, col) in enumerate(eachcol(counts_matrix))
#             probability_matrix[:, i] .= col ./ sum(col)
#         end
#     else
#         println("probability matrix found $probability_matrix_file")
#     end
#     return probability_matrix, probability_matrix_file
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert a matrix of counts into a probability matrix by normalizing each column to sum to 1.0.

# Arguments
- `counts_matrix::Matrix{<:Number}`: Input matrix where each column represents counts/frequencies

# Returns
- `Matrix{Float64}`: Probability matrix where each column sums to 1.0
"""
function count_matrix_to_probability_matrix(counts_matrix)
    probability_matrix = zeros(size(counts_matrix))
    for (i, col) in enumerate(eachcol(counts_matrix))
        probability_matrix[:, i] .= col ./ sum(col)
    end
    return probability_matrix
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)

Convert a distance matrix into a Newick tree format using UPGMA hierarchical clustering.

# Arguments
- `distance_matrix`: Square matrix of pairwise distances between entities
- `labels`: Vector of labels corresponding to the entities in the distance matrix
- `outfile`: Path where the Newick tree file will be written

# Returns
- Path to the generated Newick tree file

# Details
Performs hierarchical clustering using the UPGMA (average linkage) method and 
converts the resulting dendrogram into Newick tree format. The branch lengths 
in the tree represent the heights from the clustering.
"""
function distance_matrix_to_newick(;distance_matrix, labels, outfile)
    # phage_names = phage_host_table[indices, :name]
    # this is equivalent to UPGMA
    tree = Clustering.hclust(distance_matrix, linkage=:average, branchorder=:optimal)
    # reference_phage_indices = findall(x -> x in reference_phages, phage_names)
    newick = Dict()
    for row in 1:size(tree.merges, 1)
        left, right = tree.merges[row, :]
        if left < 0
            l = string(labels[abs(left)])
        else
            l = newick[left]
        end
        if right < 0
            r = string(labels[abs(right)])
        else
            r = newick[right]
        end
        height = tree.heights[row]
        newick[row] = "($l:$height, $r:$height)"
    end
    open(outfile, "w") do io
        println(io, newick[size(tree.merges, 1)] * ";")
    end
    return outfile
end

"""
    pairwise_distance_matrix(
        matrix;
        dist_func = Distances.euclidean,
        show_progress = true,
        progress_desc = "Computing distances"
    )

Compute a symmetric pairwise distance matrix between columns of `matrix` using the supplied distance function.

# Arguments
- `matrix`: Column-major matrix (features as rows, entities as columns)
- `dist_func`: Function of the form `f(a, b)` returning the distance between two vectors (default: `Distances.euclidean`)
- `show_progress`: Display progress bar if true (default: true)
- `progress_desc`: Progress bar description (default: "Computing distances")

# Returns
- Symmetric N×N matrix of pairwise distances between columns (entities)
"""
function pairwise_distance_matrix(
    matrix;
    dist_func = Distances.euclidean,
    show_progress = true,
    progress_desc = "Computing distances"
)
    n_entities = size(matrix, 2)
    distance_matrix = zeros(n_entities, n_entities)
    total_pairs = n_entities * (n_entities - 1) ÷ 2

    progress = show_progress ? ProgressMeter.Progress(total_pairs, desc = progress_desc, dt = 0.1) : nothing

    Threads.@threads for entity_1_index in 1:n_entities
        for entity_2_index in entity_1_index+1:n_entities
            a = matrix[:, entity_1_index]
            b = matrix[:, entity_2_index]
            dist = dist_func(a, b)
            distance_matrix[entity_1_index, entity_2_index] = dist
            distance_matrix[entity_2_index, entity_1_index] = dist
            if show_progress && Threads.threadid() == 1
                ProgressMeter.next!(progress)
            end
        end
    end

    if show_progress
        ProgressMeter.finish!(progress)
    end
    return distance_matrix
end

# Wrapper functions

"""
    frequency_matrix_to_euclidean_distance_matrix(counts_table)

Pairwise Euclidean distance between columns of `counts_table`.
"""
frequency_matrix_to_euclidean_distance_matrix(counts_table) =
    pairwise_distance_matrix(counts_table; dist_func = Distances.euclidean, progress_desc = "Euclidean distances")

"""
    frequency_matrix_to_cosine_distance_matrix(probability_matrix)

Pairwise cosine distance between columns of `probability_matrix`.
"""
frequency_matrix_to_cosine_distance_matrix(probability_matrix) =
    pairwise_distance_matrix(probability_matrix; dist_func = Distances.cosine_dist, progress_desc = "Cosine distances")

"""
    binary_matrix_to_jaccard_distance_matrix(binary_matrix::Union{BitMatrix, Matrix{Bool}})

Pairwise Jaccard distance between columns of a binary matrix (BitMatrix or Matrix{Bool}).
Throws an error if the input is not strictly a binary matrix.
"""
function binary_matrix_to_jaccard_distance_matrix(binary_matrix::Union{BitMatrix, Matrix{Bool}})
    return pairwise_distance_matrix(binary_matrix; dist_func = Distances.jaccard, progress_desc = "Jaccard distances")
end

"""
    frequency_matrix_to_jaccard_distance_matrix(matrix)

Thresholds the input matrix at `>0` to obtain a binary matrix, then computes pairwise Jaccard distance between columns.
Accepts any numeric matrix.
"""
function frequency_matrix_to_jaccard_distance_matrix(matrix)
    # Convert to binary (Bool) matrix, then call the binary version
    binary_matrix = matrix .> 0
    return binary_matrix_to_jaccard_distance_matrix(binary_matrix)
end

"""
    frequency_matrix_to_bray_curtis_distance_matrix(counts_table)

Pairwise Bray-Curtis distance between columns of `counts_table`.
"""
frequency_matrix_to_bray_curtis_distance_matrix(counts_table) =
    pairwise_distance_matrix(counts_table; dist_func = Distances.braycurtis, progress_desc = "Bray-Curtis distances")

"""
    frequency_matrix_to_jensen_shannon_distance_matrix(probability_matrix)

Pairwise Jensen-Shannon divergence between columns of `probability_matrix`.

# Arguments
- `probability_matrix`: Matrix where each column is a probability distribution (sums to 1.0).

# Returns
- Symmetric matrix of Jensen-Shannon divergence values between columns.
"""
function frequency_matrix_to_jensen_shannon_distance_matrix(probability_matrix)
    return pairwise_distance_matrix(
        probability_matrix;
        dist_func = Distances.js_divergence,
        progress_desc = "Jensen-Shannon divergence"
    )
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# DEPRECATED: THIS WAS THE MEASURE WITH THE LEAST AGREEMENT TO EXISTING MEASURES LIKE BLAST AND % AVERAGE NUCLEOTIDE IDENTITY
# Create distance matrix from a column-major counts matrix (features as rows and entities as columns)
# where distance is a proportional to total feature count magnitude (size) and cosine similarity (relative frequency)
# """
# function counts_matrix_to_size_normalized_cosine_distance_matrix(counts_table)
#     n_entities = size(counts_table, 2)
#     distance_matrix = zeros(n_entities, n_entities)
#     for entity_1_index in 1:n_entities
#         for entity_2_index in entity_1_index+1:n_entities
#             a = counts_table[:, entity_1_index]
#             b = counts_table[:, entity_2_index]
#             sa = sum(a)
#             sb = sum(b)
#             size_dist = 1-(min(sa, sb)/max(sa, sb))
#             cosine_dist = Distances.cosine_dist(a, b)
#             distances = filter(x -> x > 0, (size_dist, cosine_dist))
#             if isempty(distances)
#                 dist = 0.0
#             else
#                 dist = reduce(*, distances)
#             end
#             distance_matrix[entity_1_index, entity_2_index] = 
#                 distance_matrix[entity_2_index, entity_1_index] = dist
#         end
#     end
#     return distance_matrix
# end

# """
# Create euclidean distance matrix from a column-major counts matrix (features as rows and entities as columns)
# where distance is a proportional to total feature count magnitude (size)
# """
# function frequency_matrix_to_euclidean_distance_matrix(counts_table)
#     n_entities = size(counts_table, 2)
#     distance_matrix = zeros(n_entities, n_entities)
#     ProgressMeter.@showprogress for entity_1_index in 1:n_entities
#         for entity_2_index in entity_1_index+1:n_entities
#             a = counts_table[:, entity_1_index]
#             b = counts_table[:, entity_2_index]
#             distance_matrix[entity_1_index, entity_2_index] = 
#                 distance_matrix[entity_2_index, entity_1_index] = 
#                 Distances.euclidean(a, b)
#         end
#     end
#     return distance_matrix
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Create a Euclidean distance matrix from a column-major counts matrix
# (features as rows and entities as columns), where distance is proportional
# to total feature count magnitude (size).

# Compute pairwise Euclidean distances between entity profiles in a counts matrix.

# # Arguments
# - `counts_table`: A matrix where rows represent features and columns represent entities (column-major format).
#   Each column contains the feature counts/frequencies for one entity.

# # Returns
# - A symmetric N×N matrix of Euclidean distances between each pair of entities, where N is the number of entities.

# # Details
# - Parallelized computation using multi-threading
# - Progress tracking via ProgressMeter
# - Memory efficient: only upper triangle is computed, then mirrored
# - Distance between entities increases with total feature magnitude differences
# """
# function frequency_matrix_to_euclidean_distance_matrix(counts_table)
#     n_entities = size(counts_table, 2)
#     distance_matrix = zeros(n_entities, n_entities)

#     # Initialize a thread-safe progress meter
#     progress = ProgressMeter.Progress(n_entities * (n_entities - 1) ÷ 2, desc = "Computing distances", dt = 0.1)

#     Threads.@threads for entity_1_index in 1:n_entities
#         for entity_2_index in entity_1_index+1:n_entities
#             a = counts_table[:, entity_1_index]
#             b = counts_table[:, entity_2_index]
#             dist = Distances.euclidean(a, b)
#             distance_matrix[entity_1_index, entity_2_index] = dist
#             distance_matrix[entity_2_index, entity_1_index] = dist
#             ProgressMeter.next!(progress)  # Update the progress meter
#         end
#     end

#     ProgressMeter.finish!(progress)  # Ensure the progress meter completes
#     return distance_matrix
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Create cosine distance matrix from a column-major counts matrix (features as rows and entities as columns)
# where distance is a proportional to cosine similarity (relative frequency)

# Compute pairwise cosine distances between entities based on their feature distributions.

# # Arguments
# - `probability_matrix`: Column-major matrix where rows represent features and columns represent entities.
#   Each column should contain frequency/probability values for one entity.

# # Returns
# - Symmetric matrix of size (n_entities × n_entities) containing pairwise cosine distances.
#   Distance values range from 0 (identical distributions) to 1 (orthogonal distributions).

# # Details
# - Computes upper triangle and mirrors to lower triangle for efficiency
# - Uses `Distances.cosine_dist` for the core computation
# - Time complexity is O(n²) where n is the number of entities
# """
# function frequency_matrix_to_cosine_distance_matrix(probability_matrix)
#     n_entities = size(probability_matrix, 2)
#     distance_matrix = zeros(n_entities, n_entities)
#     for entity_1_index in 1:n_entities
#         for entity_2_index in entity_1_index+1:n_entities
#             a = probability_matrix[:, entity_1_index]
#             b = probability_matrix[:, entity_2_index]
#             distance_matrix[entity_1_index, entity_2_index] = 
#                 distance_matrix[entity_2_index, entity_1_index] = 
#                 Distances.cosine_dist(a, b)
#         end
#     end
#     return distance_matrix
# end

# didn't work
# function frequency_matrix_to_euclidean_distance_matrix(counts_table)
#     n_entities = size(counts_table, 2)
#     distance_matrix = zeros(n_entities, n_entities)
#     progress = ProgressMeter.Progress(n_entities)
#     reenrantlock = ReentrantLock()
#     Threads.@threads for entity_1_index in 1:n_entities
#         lock(reenrantlock) do
#             ProgressMeter.next!(progress)
#         end
#         for entity_2_index in entity_1_index+1:n_entities
#             a = counts_table[:, entity_1_index]
#             b = counts_table[:, entity_2_index]
#             distance_matrix[entity_1_index, entity_2_index] = 
#                 distance_matrix[entity_2_index, entity_1_index] = 
#                 Distances.euclidean(a, b)
#         end
#     end
#     return distance_matrix
# end

# function expected_kmer_frequencies(seq::BioSequences.BioSequence{A}, k::Int) where A<:BioSequences.Alphabet
#     k_minus_1_mer_counts = Dict{BioSequences.BioSequence{A}, Int}()
#     for i in 1:(length(seq) - k + 1)
#         k_minus_1_mer = seq[i:(i + k - 2)]
#         k_minus_1_mer_counts[k_minus_1_mer] = get(k_minus_1_mer_counts, k_minus_1_mer, 0) + 1
#     end

#     nucleotide_counts = Dict{BioSequences.BioSequence{A}, Int}()
#     for nucleotide in seq
#         nucleotide_seq = BioSequences.BioSequence{A}([nucleotide])
#         nucleotide_counts[nucleotide_seq] = get(nucleotide_counts, nucleotide_seq, 0) + 1
#     end

#     total_k_minus_1_mers = sum(values(k_minus_1_mer_counts))
#     expected_frequencies = Dict{BioSequences.BioSequence{A}, Float64}()
#     for (k_minus_1_mer, count) in k_minus_1_mer_counts
#         last_nucleotide = k_minus_1_mer[end]
#         last_nucleotide_seq = BioSequences.BioSequence{A}([last_nucleotide])
#         expected_frequencies[k_minus_1_mer] = (count / total_k_minus_1_mers) * (nucleotide_counts[last_nucleotide_seq] / length(seq))
#     end

#     return expected_frequencies
# end

# function precompute_expected_frequencies(sequences::Vector{BioSequences.BioSequence{A}}, k::Int) where A<:BioSequences.Alphabet
#     expected_freqs = Vector{Dict{BioSequences.BioSequence{A}, Float64}}(undef, length(sequences))
#     for (i, seq) in enumerate(sequences)
#         expected_freqs[i] = expected_kmer_frequencies(seq, k)
#     end
#     return expected_freqs
# end

# function D2_star(obs_freq1::Dict{BioSequences.BioSequence{A}, Float64},
#                  obs_freq2::Dict{BioSequences.BioSequence{A}, Float64},
#                  exp_freq::Dict{BioSequences.BioSequence{A}, Float64}) where A<:BioSequences.Alphabet
#     numerator = sum((obs_freq1[kmer] - exp_freq[kmer]) * (obs_freq2[kmer] - exp_freq[kmer]) for kmer in keys(exp_freq))
#     denominator = sqrt(sum((obs_freq1[kmer] - exp_freq[kmer])^2 for kmer in keys(exp_freq)) *
#                        sum((obs_freq2[kmer] - exp_freq[kmer])^2 for kmer in keys(exp_freq)))
#     return numerator / denominator
# end

# function d2_star_normalized(obs_freq1::Dict{BioSequences.BioSequence{A}, Float64},
#                             obs_freq2::Dict{BioSequences.BioSequence{A}, Float64},
#                             exp_freq::Dict{BioSequences.BioSequence{A}, Float64}) where A<:BioSequences.Alphabet
#     numerator = sum((obs_freq1[kmer] - exp_freq[kmer]) * (obs_freq2[kmer] - exp_freq[kmer]) for kmer in keys(exp_freq))
#     var1 = sum((obs_freq1[kmer] - exp_freq[kmer])^2 for kmer in keys(exp_freq))
#     var2 = sum((obs_freq2[kmer] - exp_freq[kmer])^2 for kmer in keys(exp_freq))
#     denominator = sqrt(var1 * var2)
#     return numerator / denominator
# end

# function compute_observed_frequencies(sequences::Vector{BioSequences.BioSequence{A}}, k::Int) where A<:BioSequences.Alphabet
#     observed_freqs = Vector{Dict{BioSequences.BioSequence{A}, Float64}}(undef, length(sequences))
#     for (i, seq) in enumerate(sequences)
#         observed_freqs[i] = observed_kmer_frequencies(seq, k)
#     end
#     return observed_freqs
# end

# function compute_pairwise_statistics(observed_freqs::Vector{Dict{BioSequences.BioSequence{A}, Float64}},
#                                      expected_freqs::Vector{Dict{BioSequences.BioSequence{A}, Float64}}) where A<:BioSequences.Alphabet
#     n = length(observed_freqs)
#     d2_star_matrix = Matrix{Float64}(undef, n, n)
#     d2_star_norm_matrix = Matrix{Float64}(undef, n, n)
#     for i in 1:n
#         for j in i:n
#             d2_star_value = D2_star(observed_freqs[i], observed_freqs[j], expected_freqs[i])
#             d2_star_norm_value = d2_star_normalized(observed_freqs[i], observed_freqs[j], expected_freqs[i])
#             d2_star_matrix[i, j] = d2_star_value
#             d2_star_matrix[j, i] = d2_star_value
#             d2_star_norm_matrix[i, j] = d2_star_norm_value
#             d2_star_norm_matrix[j, i] = d2_star_norm_value
#         end
#     end
#     return d2_star_matrix, d2_star_norm_matrix
# end