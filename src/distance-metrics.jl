# src/distance_metrics.jl

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Generate a Mash sketch from a list of FASTA files.

# Arguments
- `input_fasta_list::Vector{String}`: Paths to input FASTA files
- `output_sketch_path_prefix::String`: Prefix for the generated `.msh` file (default: "mash_sketch")
- `threads::Int`: Number of threads to use (default: `get_default_threads()`)

# Returns
- `String`: Path to the generated Mash sketch file

# Details
- Uses the Mash conda environment
- Skips sketching if the output already exists
"""
function generate_mash_sketch(;input_fasta_list::Vector{String}, output_sketch_path_prefix::String="mash_sketch", threads=get_default_threads())
    sketch_output = output_sketch_path_prefix * ".msh"
    if !isfile(sketch_output)
        mash_input_sketch_list_file = tempname() * ".mash_input.txt"
        open(mash_input_sketch_list_file, "w") do io
            for x in input_fasta_list
                println(io, x)
            end
        end
        # Use the conda environment runner to execute the mash command
        # with live streaming output
        # @info "Running mash sketch on input files: $(input_fasta_list)"
        # @info "Output sketch will be saved to: $sketch_output"
        # Run the mash sketch command with the input list file and output prefix
        Mycelia.add_bioconda_env("mash")
        run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash sketch -l $(mash_input_sketch_list_file) -p $(threads) -o $(output_sketch_path_prefix)`)
        rm(mash_input_sketch_list_file)
    end
    @assert isfile(sketch_output) "Mash sketch file was not created: $sketch_output"
    return sketch_output
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute a pairwise Mash distance matrix for a set of FASTA files.

# Arguments
- `input_fasta_list::Vector{String}`: Paths to input FASTA files
- `cleanup_sketch_output::Bool`: Remove the generated sketch after computing distances (default: true)
- `threads::Int`: Number of threads to use (default: `get_default_threads()`)

# Returns
- `Matrix{Float64}`: Symmetric matrix of Mash distances

# Details
- Reuses or generates a Mash sketch before computing distances
- Utilizes the Mash conda environment for both sketching and distance calculation
"""
function pairwise_mash_distance_matrix(;
    input_fasta_list::Vector{String},
    cleanup_sketch_output::Bool = true,
    threads=get_default_threads()
)
    mash_sketch_file = generate_mash_sketch(input_fasta_list=input_fasta_list)
    pairwise_mash_results = CSV.read(
        open(`$(Mycelia.CONDA_RUNNER) run --live-stream -n mash mash dist -p $(threads) $(mash_sketch_file) $(mash_sketch_file)`),
        DataFrames.DataFrame,
        delim='\t',
        header = ["reference", "query", "mash_distance", "p_value", "matching_hashes"]
    )
    if cleanup_sketch_output && isfile(mash_sketch_file)
        rm(mash_sketch_file)
    end

    fastx_to_index_map = Dict(fasta => i for (i, fasta) in enumerate(input_fasta_list))
    mash_distance_matrix = zeros(length(input_fasta_list), length(input_fasta_list))
    for row in DataFrames.eachrow(pairwise_mash_results)
        matrix_row = fastx_to_index_map[row["reference"]]
        matrix_column = fastx_to_index_map[row["query"]]
        if matrix_row == matrix_column
            @assert row["mash_distance"] == 0.0
        else
            mash_distance_matrix[matrix_row, matrix_column] = row["mash_distance"]
        end
    end
    return mash_distance_matrix
end


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
    b = [get(kmer_counts_2, kmer, 0) for kmer in sorted_shared_keys]
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
    b = [get(kmer_counts_2, kmer, 0) for kmer in sorted_shared_keys]
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
            a = @view matrix[:, entity_1_index]
            b = @view matrix[:, entity_2_index]
            dist = dist_func(a, b)
            distance_matrix[entity_1_index, entity_2_index] = dist
            distance_matrix[entity_2_index, entity_1_index] = dist
            if show_progress
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


# =============================================================================
# Distance Matrix Imputation and Validation Functions
# =============================================================================

"""
Imputation methods for handling missing values in distance calculations.

- `IMPUTE_MAX`: Replace missing with theoretical maximum distance (1.0)
- `IMPUTE_MAX_OBSERVED`: Replace with maximum observed non-missing distance
- `IMPUTE_MEDIAN`: Replace with median of non-missing distances
- `IMPUTE_MEAN`: Replace with mean of non-missing distances
"""
@enum ImputationMethod begin
    IMPUTE_MAX
    IMPUTE_MAX_OBSERVED
    IMPUTE_MEDIAN
    IMPUTE_MEAN
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Impute missing values (NaN, Inf) in a distance matrix.

# Arguments
- `distance_matrix::Matrix{Float64}`: Input distance matrix potentially containing NaN/Inf values
- `method::ImputationMethod`: Imputation strategy (default: IMPUTE_MAX_OBSERVED)

# Returns
- `Matrix{Float64}`: Distance matrix with imputed values

# Details
The imputation is symmetric: the same value is used for D[i,j] and D[j,i].
Diagonal values (self-distances) are always set to 0.0.
"""
function impute_distances(
    distance_matrix::Matrix{Float64};
    method::ImputationMethod = IMPUTE_MAX_OBSERVED
)
    result = copy(distance_matrix)
    n = size(result, 1)

    ## Identify missing values (NaN or Inf)
    missing_mask = isnan.(result) .| isinf.(result)

    if !any(missing_mask)
        return result
    end

    ## Extract valid (non-missing, non-diagonal) values
    valid_values = Float64[]
    for i in 1:n
        for j in (i+1):n
            val = result[i, j]
            if !isnan(val) && !isinf(val)
                push!(valid_values, val)
            end
        end
    end

    if isempty(valid_values)
        error("Cannot impute: all off-diagonal values are missing")
    end

    ## Compute fill value based on method
    fill_value = if method == IMPUTE_MAX
        1.0
    elseif method == IMPUTE_MAX_OBSERVED
        maximum(valid_values)
    elseif method == IMPUTE_MEDIAN
        Statistics.median(valid_values)
    elseif method == IMPUTE_MEAN
        Statistics.mean(valid_values)
    end

    ## Count imputed pairs
    n_imputed = 0

    ## Apply imputation symmetrically
    for i in 1:n
        for j in (i+1):n
            if isnan(result[i, j]) || isinf(result[i, j])
                result[i, j] = fill_value
                result[j, i] = fill_value
                n_imputed += 1
            end
        end
    end

    ## Ensure diagonal is zero
    for i in 1:n
        result[i, i] = 0.0
    end

    @info "Imputed $n_imputed NaN/Inf pairs with $(method) value: $(round(fill_value, digits=4))"

    return result
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute Gower distance between rows of a matrix with missing value handling.

Gower distance is the mean absolute difference normalized by feature ranges,
computed only over features where both observations have non-missing values.

# Arguments
- `matrix::AbstractMatrix`: n_samples × n_features matrix (may contain NaN for missing values)
- `feature_types::Union{Nothing, Vector{Symbol}}`: Type per feature (:numeric, :binary, :categorical).
  Default: all :numeric
- `feature_ranges::Union{Nothing, Vector{Float64}}`: Pre-computed ranges for numeric features.
  Default: computed from data
- `min_shared_features::Int`: Minimum shared non-missing features required (default: 1)

# Returns
- `Matrix{Float64}`: n_samples × n_samples distance matrix. NaN where insufficient shared features.

# Properties
- D[i,i] = 0 (diagonal is zero)
- D[i,j] = D[j,i] (symmetric)
- D[i,j] ∈ [0, 1] for valid pairs
- D[i,j] = NaN if shared features < min_shared_features

# Reference
Gower, J. C. (1971). A general coefficient of similarity and some of its properties.
Biometrics, 27(4), 857-871.
"""
function gower_distance(
    matrix::AbstractMatrix;
    feature_types::Union{Nothing, Vector{Symbol}} = nothing,
    feature_ranges::Union{Nothing, Vector{Float64}} = nothing,
    min_shared_features::Int = 1
)
    n_samples, n_features = size(matrix)

    ## Default: all numeric features
    if isnothing(feature_types)
        feature_types = fill(:numeric, n_features)
    end

    if length(feature_types) != n_features
        throw(ArgumentError("feature_types length ($(length(feature_types))) must match n_features ($n_features)"))
    end

    ## Compute ranges for numeric features if not provided
    if isnothing(feature_ranges)
        ranges = zeros(Float64, n_features)
        for j in 1:n_features
            if feature_types[j] == :numeric
                vals = filter(x -> !ismissing(x) && !(x isa Number && isnan(x)), matrix[:, j])
                if isempty(vals)
                    ranges[j] = 1.0  ## No data for this feature
                else
                    r = maximum(vals) - minimum(vals)
                    ranges[j] = r == 0.0 ? 1.0 : r  ## Avoid division by zero
                end
            else
                ranges[j] = 1.0  ## Binary/categorical: range is 1
            end
        end
    else
        if length(feature_ranges) != n_features
            throw(ArgumentError("feature_ranges length ($(length(feature_ranges))) must match n_features ($n_features)"))
        end
        ranges = feature_ranges
    end

    distance_matrix = zeros(Float64, n_samples, n_samples)

    for i in 1:n_samples
        for j in (i+1):n_samples
            sum_dist = 0.0
            count = 0

            for k in 1:n_features
                xi = matrix[i, k]
                xj = matrix[j, k]

                ## Skip if either is missing
                xi_missing = ismissing(xi) || (xi isa Number && isnan(xi))
                xj_missing = ismissing(xj) || (xj isa Number && isnan(xj))

                if xi_missing || xj_missing
                    continue
                end

                count += 1

                if feature_types[k] == :numeric
                    sum_dist += abs(xi - xj) / ranges[k]
                else  ## :binary or :categorical
                    sum_dist += (xi != xj) ? 1.0 : 0.0
                end
            end

            if count >= min_shared_features
                distance_matrix[i, j] = sum_dist / count
            else
                distance_matrix[i, j] = NaN
            end
            distance_matrix[j, i] = distance_matrix[i, j]
        end
    end

    return distance_matrix
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Convert distances to rank-normalized values.

Ranks are normalized to [0, 1] where 0 is the smallest distance and 1 is the largest.

# Arguments
- `distance_matrix::Matrix{Float64}`: Input symmetric distance matrix
- `handle_ties::Symbol`: How to handle ties (:average, :min, :max). Default: :average

# Returns
- `Matrix{Float64}`: Rank-normalized distance matrix with values in [0, 1]
"""
function normalize_ranks(
    distance_matrix::Matrix{Float64};
    handle_ties::Symbol = :average
)
    n = size(distance_matrix, 1)
    @assert size(distance_matrix, 2) == n "Distance matrix must be square"

    ## Extract upper triangle values and their indices
    upper_indices = [(i, j) for i in 1:n for j in (i+1):n]
    values = [distance_matrix[i, j] for (i, j) in upper_indices]

    ## Handle NaN values: they stay as NaN in the output
    valid_mask = .!isnan.(values)

    if !any(valid_mask)
        ## All NaN: return as-is with zeros on diagonal
        result = copy(distance_matrix)
        for i in 1:n
            result[i, i] = 0.0
        end
        return result
    end

    ## Compute ranks for valid values
    valid_values = values[valid_mask]
    valid_indices = findall(valid_mask)

    ## Use StatsBase for proper rank computation with tie handling
    ranked = if handle_ties == :average
        StatsBase.tiedrank(valid_values)
    elseif handle_ties == :min
        StatsBase.competerank(valid_values)
    elseif handle_ties == :max
        ## For max rank, compute compete rank and adjust
        n_valid = length(valid_values)
        sorted_order = sortperm(valid_values)
        ranks = zeros(n_valid)
        i = 1
        while i <= n_valid
            j = i
            while j < n_valid && valid_values[sorted_order[j]] == valid_values[sorted_order[j+1]]
                j += 1
            end
            ## All tied values get the maximum rank in the group
            for k in i:j
                ranks[sorted_order[k]] = Float64(j)
            end
            i = j + 1
        end
        ranks
    else
        throw(ArgumentError("handle_ties must be :average, :min, or :max"))
    end

    ## Normalize ranks to [0, 1]
    max_rank = maximum(ranked)
    min_rank = minimum(ranked)
    if max_rank == min_rank
        normalized_ranks = fill(0.5, length(ranked))
    else
        normalized_ranks = (ranked .- min_rank) ./ (max_rank - min_rank)
    end

    ## Build result matrix
    result = fill(NaN, n, n)

    for (k, idx) in enumerate(valid_indices)
        i, j = upper_indices[idx]
        result[i, j] = normalized_ranks[k]
        result[j, i] = normalized_ranks[k]
    end

    ## Diagonal is always 0
    for i in 1:n
        result[i, i] = 0.0
    end

    return result
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

Validate properties of a distance matrix.

# Arguments
- `distance_matrix::Matrix{Float64}`: Input distance matrix
- `check_symmetry::Bool`: Check if matrix is symmetric (default: true)
- `check_diagonal::Bool`: Check if diagonal is zero (default: true)
- `check_nonnegative::Bool`: Check if all values are non-negative (default: true)
- `tolerance::Float64`: Tolerance for floating-point comparisons (default: 1e-10)

# Returns
- `NamedTuple`: Validation results with any violations

# Example
```julia
D = [0.0 0.5 0.8; 0.5 0.0 0.3; 0.8 0.3 0.0]
result = Mycelia.validate_distance_matrix(D)
result.is_valid  # true
```
"""
function validate_distance_matrix(
    distance_matrix::Matrix{Float64};
    check_symmetry::Bool = true,
    check_diagonal::Bool = true,
    check_nonnegative::Bool = true,
    tolerance::Float64 = 1e-10
)
    n = size(distance_matrix, 1)

    ## Check square
    is_square = size(distance_matrix, 2) == n
    if !is_square
        return (
            is_valid = false,
            violations = ["Matrix is not square: $(size(distance_matrix))"],
            n_samples = n,
            has_missing = any(isnan, distance_matrix),
            n_missing = count(isnan, distance_matrix)
        )
    end

    violations = String[]

    ## Check symmetry
    is_symmetric = true
    if check_symmetry
        for i in 1:n
            for j in (i+1):n
                d_ij = distance_matrix[i, j]
                d_ji = distance_matrix[j, i]
                ## Both NaN is symmetric, otherwise compare values
                if isnan(d_ij) && isnan(d_ji)
                    continue
                elseif isnan(d_ij) || isnan(d_ji) || abs(d_ij - d_ji) > tolerance
                    is_symmetric = false
                    push!(violations, "Matrix is not symmetric at ($i, $j): $(d_ij) != $(d_ji)")
                    break
                end
            end
            if !is_symmetric
                break
            end
        end
    end

    ## Check diagonal
    is_diag_zero = true
    if check_diagonal
        for i in 1:n
            if abs(distance_matrix[i, i]) > tolerance
                is_diag_zero = false
                push!(violations, "Diagonal contains non-zero values at ($i, $i): $(distance_matrix[i, i])")
                break
            end
        end
    end

    ## Check non-negative
    is_nonneg = true
    if check_nonnegative
        for i in 1:n
            for j in 1:n
                val = distance_matrix[i, j]
                if !isnan(val) && val < -tolerance
                    is_nonneg = false
                    push!(violations, "Matrix contains negative values at ($i, $j): $(val)")
                    break
                end
            end
            if !is_nonneg
                break
            end
        end
    end

    ## Count NaN values
    n_nan = count(isnan, distance_matrix)
    ## Off-diagonal NaN pairs (each appears twice in the matrix)
    n_nan_pairs = 0
    for i in 1:n
        for j in (i+1):n
            if isnan(distance_matrix[i, j])
                n_nan_pairs += 1
            end
        end
    end

    is_valid = is_square && (!check_symmetry || is_symmetric) &&
               (!check_diagonal || is_diag_zero) && (!check_nonnegative || is_nonneg)

    return (
        is_valid = is_valid,
        violations = violations,
        n_samples = n,
        has_missing = n_nan > 0,
        n_missing = n_nan,
        n_missing_pairs = n_nan_pairs
    )
end

# =============================================================================
# Alpha Diversity Metrics
# =============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Shannon diversity index (H') for an abundance vector.

The Shannon diversity index measures the entropy of species proportions:
    H' = -∑(pᵢ × ln(pᵢ))

where pᵢ is the proportion of species i.

# Arguments
- `counts::AbstractVector{<:Real}`: Vector of species abundances/counts

# Returns
- `Float64`: Shannon diversity index (using natural log). Higher values indicate
  greater diversity. Returns 0.0 for empty or all-zero vectors.

# Notes
- Uses natural logarithm (base e) by default
- Zero counts are excluded from the calculation
- For a community with S equally abundant species, H' = ln(S)
"""
function shannon_diversity(counts::AbstractVector{<:Real})
    total = sum(counts)
    total == 0 && return 0.0

    proportions = counts ./ total
    # Filter out zeros to avoid log(0)
    nonzero = filter(x -> x > 0, proportions)
    return -sum(p * log(p) for p in nonzero)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Simpson's diversity index (1 - D) for an abundance vector.

Simpson's diversity index represents the probability that two randomly
selected individuals belong to different species:
    1 - D = 1 - ∑(pᵢ²)

where pᵢ is the proportion of species i.

# Arguments
- `counts::AbstractVector{<:Real}`: Vector of species abundances/counts

# Returns
- `Float64`: Simpson's diversity index in range [0, 1]. Higher values
  indicate greater diversity. Returns 0.0 for empty or all-zero vectors.

# Notes
- Returns 1 - Simpson's dominance (D), which is more intuitive
- A value close to 1 indicates high diversity (many species, evenly distributed)
- A value close to 0 indicates low diversity (few species or dominance)
"""
function simpsons_diversity(counts::AbstractVector{<:Real})
    total = sum(counts)
    total == 0 && return 0.0

    proportions = counts ./ total
    return 1 - sum(proportions .^ 2)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate species richness (number of taxa with non-zero abundance).

# Arguments
- `counts::AbstractVector{<:Real}`: Vector of species abundances/counts

# Returns
- `Int`: Number of species/taxa present (count > 0)
"""
function species_richness(counts::AbstractVector{<:Real})
    return count(x -> x > 0, counts)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Pielou's evenness index (J').

Pielou's evenness measures how evenly individuals are distributed among species:
    J' = H' / ln(S)

where H' is Shannon diversity and S is species richness.

# Arguments
- `counts::AbstractVector{<:Real}`: Vector of species abundances/counts

# Returns
- `Float64`: Pielou's evenness in range [0, 1]. A value of 1 indicates
  perfect evenness (all species equally abundant). Returns 0.0 if richness ≤ 1.

# Notes
- Requires at least 2 species for meaningful calculation
- Based on Shannon diversity using natural logarithm
"""
function pielous_evenness(counts::AbstractVector{<:Real})
    S = species_richness(counts)
    S <= 1 && return 0.0
    H = shannon_diversity(counts)
    return H / log(S)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate multiple alpha diversity metrics for each sample in an abundance matrix.

# Arguments
- `abundance_matrix::AbstractMatrix`: Taxa × samples abundance matrix (rows = taxa, columns = samples)
- `sample_names::AbstractVector`: Names/identifiers for each sample (column)

# Returns
- `DataFrames.DataFrame`: DataFrame with columns:
  - `sample`: Sample identifier
  - `shannon`: Shannon diversity index
  - `simpsons`: Simpson's diversity index
  - `richness`: Species richness (count of non-zero taxa)
  - `evenness`: Pielou's evenness index

# Example
```julia
# Taxa (rows) × Samples (columns)
abundance = [10 5 0; 20 15 30; 5 10 10]
samples = ["Sample1", "Sample2", "Sample3"]
alpha_df = calculate_alpha_diversity(abundance, samples)
```
"""
function calculate_alpha_diversity(abundance_matrix::AbstractMatrix, sample_names::AbstractVector)
    n_samples = size(abundance_matrix, 2)

    alpha_df = DataFrames.DataFrame(
        sample = collect(sample_names),
        shannon = zeros(n_samples),
        simpsons = zeros(n_samples),
        richness = zeros(Int, n_samples),
        evenness = zeros(n_samples)
    )

    for (i, sample) in enumerate(sample_names)
        counts = abundance_matrix[:, i]
        alpha_df[i, :shannon] = shannon_diversity(counts)
        alpha_df[i, :simpsons] = simpsons_diversity(counts)
        alpha_df[i, :richness] = species_richness(counts)
        alpha_df[i, :evenness] = pielous_evenness(counts)
    end

    return alpha_df
end

# =============================================================================
# Beta Diversity Metrics (Pairwise)
# =============================================================================

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Bray-Curtis dissimilarity between two abundance vectors.

Bray-Curtis dissimilarity is defined as:
    BC = ∑|aᵢ - bᵢ| / (∑aᵢ + ∑bᵢ)

This is equivalent to: 1 - 2×∑min(aᵢ, bᵢ) / (∑aᵢ + ∑bᵢ)

# Arguments
- `a::AbstractVector{<:Real}`: First abundance vector
- `b::AbstractVector{<:Real}`: Second abundance vector

# Returns
- `Float64`: Bray-Curtis dissimilarity in range [0, 1].
  0 indicates identical compositions, 1 indicates completely different.

# Notes
- Suitable for abundance/count data (not presence/absence)
- Sensitive to differences in both composition and magnitude
"""
function bray_curtis_dissimilarity(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    sum_total = sum(a) + sum(b)
    sum_total == 0 && return 0.0
    sum_abs_diff = sum(abs.(a .- b))
    return sum_abs_diff / sum_total
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate Jaccard distance between two abundance vectors (presence/absence).

Jaccard distance is defined as:
    J_d = 1 - |A ∩ B| / |A ∪ B|

where A and B are sets of taxa present (count > 0).

# Arguments
- `a::AbstractVector{<:Real}`: First abundance vector
- `b::AbstractVector{<:Real}`: Second abundance vector

# Returns
- `Float64`: Jaccard distance in range [0, 1].
  0 indicates identical presence/absence, 1 indicates no shared taxa.

# Notes
- Based on presence/absence only; ignores abundance magnitudes
- Use `bray_curtis_dissimilarity` when abundances matter
"""
function jaccard_distance_vectors(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    a_present = a .> 0
    b_present = b .> 0

    intersection = sum(a_present .& b_present)
    union_count = sum(a_present .| b_present)

    union_count == 0 && return 0.0
    return 1 - (intersection / union_count)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Calculate pairwise beta diversity matrix from an abundance matrix.

# Arguments
- `abundance_matrix::AbstractMatrix`: Taxa × samples abundance matrix
- `sample_names::AbstractVector`: Names/identifiers for each sample
- `metric::Symbol`: Distance metric to use:
  - `:bray_curtis` (default): Bray-Curtis dissimilarity
  - `:jaccard`: Jaccard distance (presence/absence)

# Returns
- Named tuple with:
  - `distance_matrix::Matrix{Float64}`: Symmetric pairwise distance matrix
  - `sample_names::Vector`: Sample identifiers (for row/column labels)

# Example
```julia
abundance = rand(0:100, 50, 10)  # 50 taxa, 10 samples
samples = ["S\$i" for i in 1:10]
beta = calculate_beta_diversity(abundance, samples, metric=:bray_curtis)
# beta.distance_matrix is a 10×10 symmetric matrix
```
"""
function calculate_beta_diversity(
    abundance_matrix::AbstractMatrix,
    sample_names::AbstractVector;
    metric::Symbol=:bray_curtis
)
    dist_func = if metric == :bray_curtis
        bray_curtis_dissimilarity
    elseif metric == :jaccard
        jaccard_distance_vectors
    else
        error("Unknown metric: $metric. Use :bray_curtis or :jaccard")
    end

    n_samples = size(abundance_matrix, 2)
    dist_matrix = zeros(n_samples, n_samples)

    for i in 1:n_samples
        for j in (i+1):n_samples
            d = dist_func(abundance_matrix[:, i], abundance_matrix[:, j])
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d
        end
    end

    return (distance_matrix=dist_matrix, sample_names=collect(sample_names))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Wrapper around `bray_curtis_distance` for consistency with other matrix functions.

# Arguments
- `M::AbstractMatrix`: Taxa × samples count matrix (rows = features, columns = samples)

# Returns
- `Matrix{Float64}`: Symmetric Bray-Curtis distance matrix

# Notes
- This is an alias for `bray_curtis_distance(M)` for API consistency
- Input matrix is automatically converted to appropriate integer type
"""
function frequency_matrix_to_bray_curtis_distance_matrix(M::AbstractMatrix)
    # Upcast to Int64 to prevent overflow with UInt16 etc.
    return bray_curtis_distance(Int64.(M))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

End-to-end beta diversity analysis: compute distance matrix, run PCoA, return plottable DataFrame.

This is a convenience wrapper that combines:
1. Bray-Curtis (or Jaccard) distance matrix calculation
2. PCoA ordination
3. DataFrame conversion for plotting

# Arguments
- `abundance_matrix::AbstractMatrix`: Taxa × samples abundance matrix
- `sample_names::AbstractVector`: Names/identifiers for each sample
- `metric::Symbol`: Distance metric (:bray_curtis or :jaccard, default: :bray_curtis)
- `maxoutdim::Int`: Number of PCoA dimensions to compute (default: 3)
- `metadata::Union{Nothing, DataFrames.DataFrame}`: Optional sample metadata to join

# Returns
- Named tuple with:
  - `distance_matrix`: Pairwise distance matrix
  - `pcoa`: PCoA result (model + coordinates)
  - `pcoa_df`: DataFrame ready for plotting (sample, PC1, PC2, PC3, ...)
  - `variance_explained`: Proportion of variance explained by each PC

# Example
```julia
result = Mycelia.beta_diversity_pcoa(abundance, samples, metric=:bray_curtis)
StatsPlots.scatter(result.pcoa_df.PC1, result.pcoa_df.PC2,
    title="PCoA of Bray-Curtis Distances",
    xlabel="PC1 (\$(round(result.variance_explained[1]*100, digits=1))%)",
    ylabel="PC2 (\$(round(result.variance_explained[2]*100, digits=1))%)"
)
```
"""
function beta_diversity_pcoa(
    abundance_matrix::AbstractMatrix,
    sample_names::AbstractVector;
    metric::Symbol=:bray_curtis,
    maxoutdim::Int=3,
    metadata::Union{Nothing, DataFrames.DataFrame}=nothing
)
    # Calculate distance matrix
    beta = calculate_beta_diversity(abundance_matrix, sample_names; metric=metric)

    # Run PCoA
    pcoa = pcoa_from_dist(beta.distance_matrix; maxoutdim=maxoutdim)

    # Calculate variance explained
    eigenvalues = MultivariateStats.eigvals(pcoa.model)
    total_var = sum(abs.(eigenvalues))
    variance_explained = abs.(eigenvalues) ./ total_var

    # Convert to DataFrame
    pcoa_df = pcoa_to_dataframe(pcoa, sample_names; metadata=metadata)

    return (
        distance_matrix = beta.distance_matrix,
        pcoa = pcoa,
        pcoa_df = pcoa_df,
        variance_explained = variance_explained
    )
end
