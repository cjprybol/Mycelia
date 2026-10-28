# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/oracle_strand_profiles_test.jl")'
# ```

# Brute-force oracle tests for the reduced memory profiles (:lightweight and
# :ultralight) of directed double-strand and undirected canonical k-mer graphs.
#
# Every testset builds a graph with
# `Mycelia.Rhizomorph.build_kmer_graph(records, k; mode, memory_profile)` and
# compares its vertex coverage and edge weights against
# `StrandOracle.doublestrand_expected` / `StrandOracle.canonical_expected`.
# For reduced data `count_evidence` returns `total_count`, so the comparison
# is against the aggregated counts the reduced types store.
#
# The same convention as the full-evidence oracle tests applies: an
# observation that is its own reverse complement counts twice, once per
# strand. A palindromic k-mer has coverage 2 * fwd(x); a self-reverse-complement
# edge (b == rc(a)) has weight 2 * fwdE(a, b) in the double-strand graph and
# contributes 2 * fwdE(a, b) to the canonical self-loop (canon(a), canon(a)).

import Test
import Random
import Mycelia
import FASTX
import MetaGraphsNext

if !isdefined(Main, :StrandOracle)
    include(joinpath(@__DIR__, "strand_oracle.jl"))
end

const PROFILE_ORACLE_PROFILES = (:lightweight, :ultralight)
const PROFILE_ORACLE_MODES = (:doublestrand, :canonical)

const PROFILE_ORACLE_HAND_CASES = [
    ("hand: self-reverse-complement edge (ACG→CGT)", ["ACGTTACGT"], 3),
    ("hand: distinct edge and RC-edge pair", ["AACGTT"], 3),
    ("hand: palindromic k-mer at even k", ["AACGTT"], 4),
    ("hand: two records sharing k-mers across strands", ["ACCGTA", "TACGGT"], 3)
]

function profile_oracle_records(seqs::Vector{String})::Vector{FASTX.FASTA.Record}
    return [FASTX.FASTA.Record("oracle_read_$(i)", seq) for (i, seq) in enumerate(seqs)]
end

function profile_observed_and_expected(seqs::Vector{String}, k::Int,
        profile::Symbol, mode::Symbol)::Tuple{Tuple, Tuple}
    graph = Mycelia.Rhizomorph.build_kmer_graph(profile_oracle_records(seqs), k;
        mode = mode, memory_profile = profile)
    if mode == :doublestrand
        observed = StrandOracle.graph_counts(graph; undirected = false)
        expected = StrandOracle.doublestrand_expected(seqs, k)
    elseif mode == :canonical
        observed = StrandOracle.graph_counts(graph; undirected = true)
        expected = StrandOracle.canonical_expected(seqs, k)
    else
        error("unsupported mode $(mode)")
    end
    return observed, expected
end

function profile_matches_oracle(seqs::Vector{String}, k::Int,
        profile::Symbol, mode::Symbol)::Bool
    observed, expected = profile_observed_and_expected(seqs, k, profile, mode)
    return observed == expected
end

function test_profile_oracle(seqs::Vector{String}, k::Int,
        profile::Symbol, mode::Symbol)::Nothing
    (observed_vertices,
        observed_edges),
    (expected_vertices,
        expected_edges) = profile_observed_and_expected(
        seqs, k, profile, mode)
    Test.@test observed_vertices == expected_vertices
    Test.@test observed_edges == expected_edges
    return nothing
end

"""
    reduced_data_snapshot(data)::Vector{Any}

Deep copy of every field of a reduced vertex or edge data value. The reduced
types are mutable structs without a structural `==`, so the snapshot compares
their field values (counts, per-dataset dicts, observation sets) instead.
"""
function reduced_data_snapshot(data)::Vector{Any}
    return Any[deepcopy(getfield(data, name)) for name in fieldnames(typeof(data))]
end

function reduced_graph_snapshot(graph)::Tuple{Dict, Dict}
    vertex_data = Dict(string(label) => reduced_data_snapshot(graph[label])
    for label in MetaGraphsNext.labels(graph))
    edge_data = Dict((string(src), string(dst)) => reduced_data_snapshot(graph[src, dst])
    for (src, dst) in MetaGraphsNext.edge_labels(graph))
    return vertex_data, edge_data
end

Test.@testset "Reduced memory profiles match brute-force strand oracle" begin
    for profile in PROFILE_ORACLE_PROFILES, mode in PROFILE_ORACLE_MODES

        Test.@testset "$(profile) $(mode)" begin
            for (name, seqs, k) in PROFILE_ORACLE_HAND_CASES
                Test.@testset "$(name)" begin
                    test_profile_oracle(seqs, k, profile, mode)
                end
            end

            Test.@testset "randomized, 50 cases" begin
                rng = Random.MersenneTwister(20261028)
                failures = Tuple{Vector{String}, Int}[]
                for i in 1:50
                    k = rand(rng, [3, 4, 5, 7, 11])
                    n = rand(rng, 20:200)
                    n_seqs = rand(rng, 1:2)
                    seqs = [StrandOracle.random_non_rc_closed(rng, n, k)
                            for _ in 1:n_seqs]
                    if !profile_matches_oracle(seqs, k, profile, mode)
                        push!(failures, (seqs, k))
                    end
                end
                if !isempty(failures)
                    @info "randomized reduced-profile oracle failures" profile mode n_failures=length(failures) failing_k=sort(unique(last.(failures))) first_failure=first(failures)
                end
                Test.@test isempty(failures)
            end
        end
    end

    for profile in PROFILE_ORACLE_PROFILES
        Test.@testset "$(profile): input single-strand graph is not mutated" begin
            for seqs in (["ACCGTA", "TACGGT"], ["AACGTT"])
                ss = Mycelia.Rhizomorph.build_kmer_graph(profile_oracle_records(seqs), 3;
                    mode = :singlestrand, memory_profile = profile)
                counts_before = StrandOracle.graph_counts(ss; undirected = false)
                data_before = reduced_graph_snapshot(ss)
                Mycelia.Rhizomorph.convert_to_doublestrand(ss)
                Mycelia.Rhizomorph.convert_to_canonical(ss)
                Test.@test StrandOracle.graph_counts(ss; undirected = false) ==
                           counts_before
                Test.@test reduced_graph_snapshot(ss) == data_before
            end
        end
    end
end
