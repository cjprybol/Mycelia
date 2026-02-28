# Quality metrics for evaluating generated sequences.
#
# Provides GC content, k-mer profiles, character frequency profiles,
# and composite generation quality evaluation.

"""
    gc_content(seq::AbstractString)

Compute GC content of a nucleotide sequence (case-insensitive).

# Returns
- `Float64`: Fraction of G and C characters (0.0 to 1.0). Returns 0.0 for empty sequences.
"""
function gc_content(seq::AbstractString)
    n = length(seq)
    if n == 0
        return 0.0
    end
    gc_count = count(c -> c in ('G', 'C', 'g', 'c'), seq)
    return gc_count / n
end

"""
    kmer_frequency_profile(seq::AbstractString, k::Int)

Compute normalized k-mer frequency profile using a sliding window.

# Returns
- `Dict{String, Float64}`: K-mer to normalized frequency (sums to 1.0).
  Returns empty Dict for sequences shorter than k.
"""
function kmer_frequency_profile(seq::AbstractString, k::Int)
    s = string(seq)
    n = length(s)
    if n < k
        return Dict{String, Float64}()
    end

    counts = Dict{String, Int}()
    total = 0
    for i in 1:(n - k + 1)
        kmer = s[i:(i + k - 1)]
        counts[kmer] = get(counts, kmer, 0) + 1
        total += 1
    end

    if total == 0
        return Dict{String, Float64}()
    end

    profile = Dict{String, Float64}()
    for (kmer, count) in counts
        profile[kmer] = count / total
    end

    return profile
end

"""
    char_frequency_profile(text::AbstractString)

Compute normalized character frequency profile.

# Returns
- `Dict{Char, Float64}`: Character to normalized frequency (sums to 1.0).
  Returns empty Dict for empty strings.
"""
function char_frequency_profile(text::AbstractString)
    n = length(text)
    if n == 0
        return Dict{Char, Float64}()
    end

    counts = Dict{Char, Int}()
    for c in text
        counts[c] = get(counts, c, 0) + 1
    end

    profile = Dict{Char, Float64}()
    for (c, count) in counts
        profile[c] = count / n
    end

    return profile
end

"""
    evaluate_generation_quality(generated, training_graph; metrics=[:gc_content, :kmer_divergence, :zipf, :diversity])

Compare generated sequences against training corpus distributions derived from the graph.

# Metrics
- `:gc_content` — Mean GC content of generated sequences
- `:kmer_divergence` — Jensen-Shannon divergence between generated and training k-mer profiles
- `:zipf` — Zipf exponent of generated character frequencies
- `:diversity` — Fraction of unique k-mers in generated sequences vs graph vertices

# Returns
- `Dict{Symbol, Float64}`: Metric name to value
"""
function evaluate_generation_quality(
        generated::Vector{<:AbstractString},
        training_graph::MetaGraphsNext.MetaGraph;
        metrics::Vector{Symbol} = [:gc_content, :kmer_divergence, :zipf, :diversity]
)
    results = Dict{Symbol, Float64}()

    if isempty(generated)
        for m in metrics
            results[m] = 0.0
        end
        return results
    end

    # Determine k from first vertex label
    labels = collect(MetaGraphsNext.labels(training_graph))
    k = isempty(labels) ? 3 : length(string(first(labels)))

    if :gc_content in metrics
        gc_values = [gc_content(string(seq)) for seq in generated]
        results[:gc_content] = Statistics.mean(gc_values)
    end

    if :kmer_divergence in metrics
        # Build generated k-mer profile from all sequences
        gen_counts = Dict{String, Int}()
        gen_total = 0
        for seq in generated
            s = string(seq)
            for i in 1:(length(s) - k + 1)
                kmer = s[i:(i + k - 1)]
                gen_counts[kmer] = get(gen_counts, kmer, 0) + 1
                gen_total += 1
            end
        end
        gen_profile = Dict{String, Float64}()
        if gen_total > 0
            for (kmer, count) in gen_counts
                gen_profile[kmer] = count / gen_total
            end
        end

        # Build training profile from graph vertices
        train_counts = Dict{String, Int}()
        train_total = 0
        for label in labels
            key = string(label)
            count = try
                count_evidence(training_graph[label])
            catch
                1
            end
            train_counts[key] = count
            train_total += count
        end
        train_profile = Dict{String, Float64}()
        if train_total > 0
            for (kmer, count) in train_counts
                train_profile[kmer] = count / train_total
            end
        end

        results[:kmer_divergence] = jensen_shannon_divergence(gen_profile, train_profile)
    end

    if :zipf in metrics
        # Character frequencies from all generated sequences
        char_counts = Dict{Char, Int}()
        for seq in generated
            for c in string(seq)
                char_counts[c] = get(char_counts, c, 0) + 1
            end
        end
        freq_vec = collect(values(char_counts))
        results[:zipf] = isempty(freq_vec) ? 0.0 : estimate_zipf_exponent(freq_vec)
    end

    if :diversity in metrics
        # Unique k-mers in generated sequences vs training graph vertices
        gen_kmers = Set{String}()
        for seq in generated
            s = string(seq)
            for i in 1:(length(s) - k + 1)
                push!(gen_kmers, s[i:(i + k - 1)])
            end
        end
        n_training = length(labels)
        results[:diversity] = n_training > 0 ? min(1.0, length(gen_kmers) / n_training) :
                              0.0
    end

    return results
end
