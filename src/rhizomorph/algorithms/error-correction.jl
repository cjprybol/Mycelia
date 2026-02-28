# Graph-based error correction algorithms.
#
# Migrated from panlingual-metagraphs notebook 12.
# Provides greedy sequence error correction using k-mer graph structure.

"""
    hamming_distance(s1::AbstractString, s2::AbstractString)

Compute Hamming distance between two strings of equal length.

# Returns
- `Int`: Number of positions where characters differ
"""
function hamming_distance(s1::AbstractString, s2::AbstractString)
    return sum(c1 != c2 for (c1, c2) in zip(s1, s2))
end

"""
    extract_kmers(seq::AbstractString, k::Int)

Extract all k-mers from a sequence using a sliding window.

# Returns
- `Vector{String}`: K-mers in order of occurrence
"""
function extract_kmers(seq::AbstractString, k::Int)
    s = string(seq)
    return [s[i:(i + k - 1)] for i in 1:(length(s) - k + 1)]
end

"""
    correct_sequence_greedy(query, k, kmer_counts, edge_counts)

Greedy error correction using k-mer graph structure.

For each position, finds the k-mer with highest coverage that extends from
the previous corrected k-mer, allowing single-substitution corrections.

# Arguments
- `query::AbstractString`: Sequence to correct
- `k::Int`: K-mer size
- `kmer_counts::Dict{String, <:Number}`: K-mer to observation count
- `edge_counts::Dict{Tuple{String, String}, <:Number}`: Edge to weight

# Returns
- `Tuple{String, Vector{Int}}`: (corrected sequence, indices of corrected k-mer positions)
"""
function correct_sequence_greedy(
        query::AbstractString,
        k::Int,
        kmer_counts::Dict,
        edge_counts::Dict
)
    query_kmers = extract_kmers(string(query), k)
    corrected_kmers = String[]
    corrections_made = Int[]

    for (i, qkmer) in enumerate(query_kmers)
        if i == 1
            if haskey(kmer_counts, qkmer)
                push!(corrected_kmers, qkmer)
            else
                best_kmer = qkmer
                best_count = 0
                for (kmer, count) in kmer_counts
                    if hamming_distance(kmer, qkmer) <= 1 && count > best_count
                        best_kmer = kmer
                        best_count = count
                    end
                end
                push!(corrected_kmers, best_kmer)
                if best_kmer != qkmer
                    push!(corrections_made, i)
                end
            end
        else
            prev_kmer = corrected_kmers[end]
            best_kmer = qkmer
            best_count = 0

            for (kmer, count) in kmer_counts
                edge = (prev_kmer, kmer)
                if haskey(edge_counts, edge)
                    edge_weight = edge_counts[edge]
                    if edge_weight > best_count
                        best_kmer = kmer
                        best_count = edge_weight
                    end
                end
            end

            push!(corrected_kmers, best_kmer)
            if best_kmer != qkmer
                push!(corrections_made, i)
            end
        end
    end

    if isempty(corrected_kmers)
        return string(query), corrections_made
    end

    corrected_seq = corrected_kmers[1]
    for i in 2:length(corrected_kmers)
        corrected_seq *= corrected_kmers[i][end:end]
    end

    return corrected_seq, corrections_made
end

"""
    correct_sequence_greedy(query, k, graph)

Greedy error correction using a Rhizomorph MetaGraph.
Extracts k-mer counts and edge weights from the graph structure.

# Arguments
- `query::AbstractString`: Sequence to correct
- `k::Int`: K-mer size
- `graph::MetaGraphsNext.MetaGraph`: Rhizomorph graph with evidence data

# Returns
- `Tuple{String, Vector{Int}}`: (corrected sequence, indices of corrected k-mer positions)
"""
function correct_sequence_greedy(
        query::AbstractString,
        k::Int,
        graph::MetaGraphsNext.MetaGraph
)
    kmer_counts = Dict{String, Int}()
    for label in MetaGraphsNext.labels(graph)
        key = string(label)
        kmer_counts[key] = try
            count_evidence(graph[label])
        catch
            1
        end
    end

    edge_counts = Dict{Tuple{String, String}, Int}()
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        key = (string(src), string(dst))
        edge_counts[key] = try
            count_evidence(graph[src, dst])
        catch
            edata = graph[src, dst]
            if edata isa StrandWeightedEdgeData
                round(Int, edata.weight)
            else
                1
            end
        end
    end

    return correct_sequence_greedy(query, k, kmer_counts, edge_counts)
end
