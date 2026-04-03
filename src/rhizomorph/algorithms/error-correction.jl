# Graph-based error correction algorithms.
#
# Migrated from panlingual-metagraphs notebook 12.
# Provides greedy sequence error correction using k-mer graph structure.

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Compute the number of mismatched character positions between two strings.

# Arguments
- `s1::AbstractString`: First string to compare.
- `s2::AbstractString`: Second string to compare.

# Returns
- `Int`: Number of zipped character pairs whose values differ.

# Example
```julia
distance = Mycelia.Rhizomorph.hamming_distance("ABC", "ABD")
```
"""
function hamming_distance(s1::AbstractString, s2::AbstractString)
    return sum(c1 != c2 for (c1, c2) in zip(s1, s2))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Extract all k-mers from a sequence using a sliding window.

# Arguments
- `seq::AbstractString`: Input sequence to scan.
- `k::Int`: Length of each k-mer window.

# Returns
- `Vector{String}`: K-mers in order of occurrence. Returns an empty vector when `k > length(seq)`.

# Example
```julia
kmers = Mycelia.Rhizomorph.extract_kmers("ABCDE", 3)
```
"""
function extract_kmers(seq::AbstractString, k::Int)
    s = string(seq)
    return [s[i:(i + k - 1)] for i in 1:(length(s) - k + 1)]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

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

# Example
```julia
kmer_counts = Dict("ABC" => 3, "BCD" => 3, "CDE" => 3)
edge_counts = Dict(("ABC", "BCD") => 2, ("BCD", "CDE") => 2)
corrected, corrected_positions = Mycelia.Rhizomorph.correct_sequence_greedy(
    "ABXDE",
    3,
    kmer_counts,
    edge_counts,
)
```
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
$(DocStringExtensions.TYPEDSIGNATURES)

Greedy error correction using a Rhizomorph MetaGraph.

Extracts k-mer counts and edge weights from the graph structure before
delegating to the dictionary-backed method.

# Arguments
- `query::AbstractString`: Sequence to correct
- `k::Int`: K-mer size
- `graph::MetaGraphsNext.MetaGraph`: Rhizomorph graph with evidence data

# Returns
- `Tuple{String, Vector{Int}}`: (corrected sequence, indices of corrected k-mer positions)

# Example
```julia
corrected, corrected_positions = Mycelia.Rhizomorph.correct_sequence_greedy(
    "ACGTA",
    3,
    graph,
)
```
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
