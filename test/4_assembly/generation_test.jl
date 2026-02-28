# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/generation_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/generation_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Smoke tests for the five new Rhizomorph modules:
#   algorithms/generation.jl     - Batch sequence generation via random walks
#   algorithms/metrics.jl        - Graph metrics (centrality, modularity, Betti)
#   algorithms/error-correction.jl - Greedy error correction
#   analysis/information-theory.jl - Entropy, divergence, Zipf
#   analysis/sequence-quality.jl   - GC content, k-mer profiles, evaluation

import Test
import Mycelia
import MetaGraphsNext

# ============================================================================
# Helper: build a small 3-gram string graph for testing
# ============================================================================

function _build_test_ngram_graph()
    # "ABCDEF" with n=3 → vertices: ABC, BCD, CDE, DEF
    # edges: ABC→BCD, BCD→CDE, CDE→DEF (linear chain)
    strings = ["ABCDEF"]
    return Mycelia.Rhizomorph.build_ngram_graph(strings, 3; dataset_id = "gen_test")
end

function _build_test_dna_ngram_graph()
    # "ATCGATCG" with n=3 → 6 vertices, 5 edges (some may overlap)
    strings = ["ATCGATCG"]
    return Mycelia.Rhizomorph.build_ngram_graph(strings, 3; dataset_id = "dna_gen_test")
end

# ============================================================================
# 1. Generation module tests
# ============================================================================

Test.@testset "Generation - generate_sequences" begin
    graph = _build_test_ngram_graph()

    # Convert to weighted graph (required by generation API)
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

    # Basic generation: 5 sequences, short walk
    results = Mycelia.Rhizomorph.generate_sequences(weighted, 5; walk_length = 10, seed = 42)

    Test.@test length(results) > 0
    Test.@test length(results) <= 5

    for r in results
        Test.@test haskey(r, :sequence)
        Test.@test haskey(r, :path)
        Test.@test haskey(r, :walk_probability)
        Test.@test haskey(r, :length)
        Test.@test r.walk_probability > 0.0
        Test.@test r.length > 0
    end
end

Test.@testset "Generation - select_start_vertices" begin
    graph = _build_test_ngram_graph()
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

    for strategy in [:degree_weighted, :uniform, :largest_component]
        vertices = Mycelia.Rhizomorph.select_start_vertices(weighted, 3; strategy = strategy)
        Test.@test length(vertices) == 3

        all_labels = Set(MetaGraphsNext.labels(weighted))
        for v in vertices
            Test.@test v in all_labels
        end
    end
end

Test.@testset "Generation - compute_sequence_likelihood" begin
    graph = _build_test_ngram_graph()
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

    # A valid traversal of the graph should have finite likelihood
    ll = Mycelia.Rhizomorph.compute_sequence_likelihood("ABCDEF", weighted; k = 3)
    Test.@test isfinite(ll)
    Test.@test ll <= 0.0  # log-likelihood is non-positive

    # A sequence too short for any transition → 0.0
    ll_short = Mycelia.Rhizomorph.compute_sequence_likelihood("AB", weighted; k = 3)
    Test.@test ll_short == 0.0

    # A sequence with impossible k-mer transition → -Inf
    ll_impossible = Mycelia.Rhizomorph.compute_sequence_likelihood("XYZXYZ", weighted; k = 3)
    Test.@test ll_impossible == -Inf
end

Test.@testset "Generation - temperature control" begin
    graph = _build_test_ngram_graph()
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

    # Very low temperature should produce near-deterministic walks
    results_cold = Mycelia.Rhizomorph.generate_sequences(
        weighted, 3; walk_length = 10, seed = 42, temperature = 0.01)
    Test.@test length(results_cold) > 0

    # High temperature should still produce valid sequences
    results_hot = Mycelia.Rhizomorph.generate_sequences(
        weighted, 3; walk_length = 10, seed = 42, temperature = 5.0)
    Test.@test length(results_hot) > 0

    # Zero temperature should throw
    Test.@test_throws ArgumentError Mycelia.Rhizomorph.generate_sequences(
        weighted, 1; walk_length = 5, temperature = 0.0)
end

# ============================================================================
# 2. Information theory module tests
# ============================================================================

Test.@testset "Information Theory - shannon_entropy" begin
    # Fair coin: entropy = 1.0 bit
    Test.@test Mycelia.Rhizomorph.shannon_entropy([0.5, 0.5]) ≈ 1.0

    # Certain outcome: entropy = 0.0
    Test.@test Mycelia.Rhizomorph.shannon_entropy([1.0]) ≈ 0.0

    # Uniform over 4 outcomes: entropy = 2.0 bits
    Test.@test Mycelia.Rhizomorph.shannon_entropy([0.25, 0.25, 0.25, 0.25]) ≈ 2.0

    # Zero probabilities are handled safely
    Test.@test Mycelia.Rhizomorph.shannon_entropy([0.5, 0.5, 0.0]) ≈ 1.0
end

Test.@testset "Information Theory - jensen_shannon_divergence" begin
    p = Dict("A" => 0.5, "B" => 0.5)
    q = Dict("A" => 0.5, "B" => 0.5)

    # Identical distributions → JSD = 0
    Test.@test Mycelia.Rhizomorph.jensen_shannon_divergence(p, q) ≈ 0.0 atol = 1e-10

    # Different distributions → JSD > 0
    r = Dict("A" => 1.0, "B" => 0.0)
    jsd = Mycelia.Rhizomorph.jensen_shannon_divergence(p, r)
    Test.@test jsd > 0.0
    Test.@test jsd <= 1.0  # JSD is bounded by log2(2) = 1.0 for two distributions
end

Test.@testset "Information Theory - kl_divergence" begin
    p = Dict("A" => 0.5, "B" => 0.5)
    q = Dict("A" => 0.5, "B" => 0.5)

    # Identical distributions → KL = 0
    Test.@test Mycelia.Rhizomorph.kl_divergence(p, q) ≈ 0.0 atol = 1e-10
end

Test.@testset "Information Theory - estimate_zipf_exponent" begin
    # Perfect Zipf with exponent ~1: frequencies = [100, 50, 33, 25, 20]
    freqs = [100, 50, 33, 25, 20]
    exponent = Mycelia.Rhizomorph.estimate_zipf_exponent(freqs)
    Test.@test isfinite(exponent)
    Test.@test exponent > 0.0  # Zipf exponent should be positive
end

Test.@testset "Information Theory - graph entropy" begin
    graph = _build_test_ngram_graph()
    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

    # Mean transition entropy should be non-negative
    mean_entropy = Mycelia.Rhizomorph.graph_transition_entropy(weighted)
    Test.@test mean_entropy >= 0.0

    # Vertex transition entropy for a specific vertex
    labels = collect(MetaGraphsNext.labels(weighted))
    if !isempty(labels)
        vte = Mycelia.Rhizomorph.vertex_transition_entropy(weighted, first(labels))
        Test.@test vte >= 0.0
    end
end

Test.@testset "Information Theory - vertex/edge distributions" begin
    graph = _build_test_ngram_graph()

    label_freqs = Mycelia.Rhizomorph.vertex_label_frequencies(graph)
    Test.@test !isempty(label_freqs)
    Test.@test all(v -> v > 0, values(label_freqs))

    weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
    edge_weights = Mycelia.Rhizomorph.edge_weight_distribution(weighted)
    Test.@test !isempty(edge_weights)
    Test.@test all(w -> w > 0.0, edge_weights)
end

# ============================================================================
# 3. Sequence quality module tests
# ============================================================================

Test.@testset "Sequence Quality - gc_content" begin
    Test.@test Mycelia.Rhizomorph.gc_content("GCGC") ≈ 1.0
    Test.@test Mycelia.Rhizomorph.gc_content("ATAT") ≈ 0.0
    Test.@test Mycelia.Rhizomorph.gc_content("ATGC") ≈ 0.5
    Test.@test Mycelia.Rhizomorph.gc_content("") ≈ 0.0
    # Case insensitive
    Test.@test Mycelia.Rhizomorph.gc_content("gcgc") ≈ 1.0
end

Test.@testset "Sequence Quality - kmer_frequency_profile" begin
    profile = Mycelia.Rhizomorph.kmer_frequency_profile("AABB", 2)
    Test.@test haskey(profile, "AA")
    Test.@test haskey(profile, "AB")
    Test.@test haskey(profile, "BB")
    Test.@test length(profile) == 3
    # Frequencies should sum to ~1.0
    Test.@test sum(values(profile)) ≈ 1.0

    # Sequence shorter than k → empty
    empty_profile = Mycelia.Rhizomorph.kmer_frequency_profile("A", 3)
    Test.@test isempty(empty_profile)
end

Test.@testset "Sequence Quality - char_frequency_profile" begin
    profile = Mycelia.Rhizomorph.char_frequency_profile("AABC")
    Test.@test profile['A'] ≈ 0.5
    Test.@test profile['B'] ≈ 0.25
    Test.@test profile['C'] ≈ 0.25
    Test.@test sum(values(profile)) ≈ 1.0

    # Empty string → empty dict
    Test.@test isempty(Mycelia.Rhizomorph.char_frequency_profile(""))
end

Test.@testset "Sequence Quality - evaluate_generation_quality" begin
    graph = _build_test_dna_ngram_graph()

    generated = ["ATCGATCG", "ATCGATC"]

    quality = Mycelia.Rhizomorph.evaluate_generation_quality(
        generated, graph; metrics = [:gc_content, :kmer_divergence, :zipf, :diversity])

    Test.@test haskey(quality, :gc_content)
    Test.@test haskey(quality, :kmer_divergence)
    Test.@test haskey(quality, :zipf)
    Test.@test haskey(quality, :diversity)

    Test.@test quality[:gc_content] >= 0.0
    Test.@test quality[:gc_content] <= 1.0
    Test.@test quality[:kmer_divergence] >= 0.0
    Test.@test quality[:diversity] >= 0.0
end

# ============================================================================
# 4. Metrics module tests
# ============================================================================

Test.@testset "Metrics - eigenvector_centrality" begin
    graph = _build_test_ngram_graph()
    centrality = Mycelia.Rhizomorph.compute_eigenvector_centrality(graph)

    Test.@test length(centrality) == MetaGraphsNext.nv(graph)
    Test.@test all(c -> c >= 0.0, centrality)
end

Test.@testset "Metrics - modularity" begin
    graph = _build_test_ngram_graph()
    mod_val, communities = Mycelia.Rhizomorph.compute_modularity(graph)

    # Modularity is between -0.5 and 1.0 for undirected; can vary for directed
    Test.@test isfinite(mod_val)
    Test.@test length(communities) == MetaGraphsNext.nv(graph)
end

Test.@testset "Metrics - closeness_centrality" begin
    graph = _build_test_ngram_graph()
    closeness = Mycelia.Rhizomorph.compute_closeness_centrality(graph)

    Test.@test length(closeness) == MetaGraphsNext.nv(graph)
    Test.@test all(c -> c >= 0.0, closeness)
end

Test.@testset "Metrics - betti_numbers" begin
    graph = _build_test_ngram_graph()
    beta_0, beta_1 = Mycelia.Rhizomorph.compute_betti_numbers(graph)

    # Linear chain: 1 connected component, 0 independent cycles
    Test.@test beta_0 >= 1
    Test.@test beta_1 >= 0
end

# ============================================================================
# 5. Error correction module tests
# ============================================================================

Test.@testset "Error Correction - hamming_distance" begin
    Test.@test Mycelia.Rhizomorph.hamming_distance("ABC", "ABC") == 0
    Test.@test Mycelia.Rhizomorph.hamming_distance("ABC", "ABD") == 1
    Test.@test Mycelia.Rhizomorph.hamming_distance("ABC", "XYZ") == 3
end

Test.@testset "Error Correction - extract_kmers" begin
    kmers = Mycelia.Rhizomorph.extract_kmers("ABCDE", 3)
    Test.@test kmers == ["ABC", "BCD", "CDE"]
    Test.@test isempty(Mycelia.Rhizomorph.extract_kmers("AB", 3))
end

Test.@testset "Error Correction - correct_sequence_greedy (Dict)" begin
    # Build k-mer and edge count dicts from a known-good sequence
    good_seq = "ABCDEF"
    k = 3
    kmers = Mycelia.Rhizomorph.extract_kmers(good_seq, k)
    kmer_counts = Dict{String, Int}()
    for kmer in kmers
        kmer_counts[kmer] = get(kmer_counts, kmer, 0) + 1
    end

    # Edge keys are Tuple{String, String} (not underscore-separated strings)
    edge_counts = Dict{Tuple{String, String}, Int}()
    for i in 1:(length(kmers) - 1)
        edge_key = (kmers[i], kmers[i + 1])
        edge_counts[edge_key] = get(edge_counts, edge_key, 0) + 1
    end

    # Introduce a 1-character error
    bad_seq = "ABXDEF"
    corrected_seq,
    correction_positions = Mycelia.Rhizomorph.correct_sequence_greedy(
        bad_seq, k, kmer_counts, edge_counts)
    Test.@test corrected_seq isa AbstractString
    Test.@test length(corrected_seq) == length(bad_seq)
    Test.@test correction_positions isa Vector{Int}
end
