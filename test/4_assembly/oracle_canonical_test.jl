# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/oracle_canonical_test.jl")'
# ```

# Brute-force oracle tests for undirected canonical k-mer graphs.
#
# Every testset builds a graph with
# `Mycelia.Rhizomorph.build_kmer_graph(records, k; mode = :canonical)` and
# compares its vertex coverage and edge weights (evidence-entry counts) against
# `StrandOracle.canonical_expected`, which is computed with plain String
# operations only. Edge keys are unordered pairs of canonical k-mers.
#
# One convention for both strand modes: an observation that is its own reverse
# complement counts twice, once per strand. A palindromic k-mer has coverage
# 2 * fwd(x), and a self-reverse-complement edge (b == rc(a), e.g. ACG -> CGT)
# contributes 2 * fwdE(a, b) to the self-loop (canon(a), canon(a)).

import Test
import Random
import Mycelia
import FASTX
import Graphs
import MetaGraphsNext

if !isdefined(Main, :StrandOracle)
    include(joinpath(@__DIR__, "strand_oracle.jl"))
end

function canonical_oracle_records(seqs::Vector{String})::Vector{FASTX.FASTA.Record}
    return [FASTX.FASTA.Record("oracle_read_$(i)", seq) for (i, seq) in enumerate(seqs)]
end

function build_canonical_oracle_graph(
        seqs::Vector{String}, k::Int)::MetaGraphsNext.MetaGraph
    return Mycelia.Rhizomorph.build_kmer_graph(canonical_oracle_records(seqs), k;
        mode = :canonical)
end

function canonical_matches_oracle(seqs::Vector{String}, k::Int)::Bool
    graph = build_canonical_oracle_graph(seqs, k)
    observed_vertices,
    observed_edges = StrandOracle.graph_counts(graph; undirected = true)
    expected_vertices, expected_edges = StrandOracle.canonical_expected(seqs, k)
    return observed_vertices == expected_vertices && observed_edges == expected_edges
end

function test_canonical_oracle(seqs::Vector{String}, k::Int)::Nothing
    graph = build_canonical_oracle_graph(seqs, k)
    observed_vertices,
    observed_edges = StrandOracle.graph_counts(graph; undirected = true)
    expected_vertices, expected_edges = StrandOracle.canonical_expected(seqs, k)
    Test.@test observed_vertices == expected_vertices
    Test.@test observed_edges == expected_edges
    return nothing
end

function canonical_evidence_snapshot(graph::MetaGraphsNext.MetaGraph)::Tuple{
        Dict{String}, Dict{Tuple{String, String}}}
    vertex_evidence = Dict(string(label) => deepcopy(graph[label].evidence)
    for label in MetaGraphsNext.labels(graph))
    edge_evidence = Dict((string(src), string(dst)) => deepcopy(graph[src, dst].evidence)
    for (src, dst) in MetaGraphsNext.edge_labels(graph))
    return vertex_evidence, edge_evidence
end

function unordered_label_pairs(
        graph::MetaGraphsNext.MetaGraph)::Set{Tuple{String, String}}
    pairs = Set{Tuple{String, String}}()
    for (src, dst) in MetaGraphsNext.edge_labels(graph)
        s = string(src)
        d = string(dst)
        push!(pairs, d < s ? (d, s) : (s, d))
    end
    return pairs
end

function gfa_segment_depths(path::AbstractString)::Dict{String, Float64}
    depths = Dict{String, Float64}()
    for line in eachline(path)
        startswith(line, "S\t") || continue
        fields = split(line, '\t')
        dp_fields = filter(field -> startswith(field, "DP:f:"), fields)
        depths[String(fields[3])] = parse(Float64, only(dp_fields)[6:end])
    end
    return depths
end

const CANONICAL_HAND_CASES = [
    ("hand: self-reverse-complement edge (ACG→CGT)", ["ACGTTACGT"], 3),
    ("hand: distinct edge and RC-edge pair", ["AACGTT"], 3),
    ("hand: palindromic k-mer at even k", ["AACGTT"], 4),
    ("hand: two records sharing k-mers across strands", ["ACCGTA", "TACGGT"], 3)
]

Test.@testset "Canonical k-mer graph matches brute-force oracle" begin
    for (name, seqs, k) in CANONICAL_HAND_CASES
        Test.@testset "$(name)" begin
            test_canonical_oracle(seqs, k)
        end
    end

    Test.@testset "randomized, 200 cases" begin
        rng = Random.MersenneTwister(20261027)
        failures = Tuple{Vector{String}, Int}[]
        for i in 1:200
            k = rand(rng, [3, 4, 5, 7, 11])
            n = rand(rng, 20:200)
            n_seqs = rand(rng, 1:2)
            seqs = [StrandOracle.random_non_rc_closed(rng, n, k) for _ in 1:n_seqs]
            if !canonical_matches_oracle(seqs, k)
                push!(failures, (seqs, k))
            end
        end
        if !isempty(failures)
            @info "randomized canonical oracle failures" n_failures=length(failures) failing_k=sort(unique(last.(failures))) first_failure=first(failures)
        end
        Test.@test isempty(failures)
    end

    Test.@testset "phiX174 reference" begin
        genome_info = Mycelia.get_test_genome_fasta(use_ncbi = true)
        try
            if genome_info.source != :ncbi
                @warn "phiX174 download unavailable; got $(genome_info.source) genome, marking broken"
                Test.@test_broken genome_info.source == :ncbi
            else
                seqs = String[]
                for record in Mycelia.open_fastx(genome_info.fasta)
                    push!(seqs, uppercase(String(FASTX.sequence(record))))
                end
                Test.@test length(seqs) == 1
                Test.@test all(seq -> all(c -> c in "ACGT", seq), seqs)
                for k in (11, 21)
                    test_canonical_oracle(seqs, k)
                end
            end
        finally
            genome_info.cleanup()
        end
    end

    Test.@testset "input single-strand graph is not mutated" begin
        # (["AACGTT"], 4) has the palindromic k-mer ACGT, exercising the
        # palindrome merge branch.
        for (seqs, k) in ((["ACCGTA", "TACGGT"], 3), (["AACGTT"], 3), (["AACGTT"], 4))
            ss = Mycelia.Rhizomorph.build_kmer_graph(canonical_oracle_records(seqs), k;
                mode = :singlestrand)
            counts_before = StrandOracle.graph_counts(ss; undirected = false)
            evidence_before = canonical_evidence_snapshot(ss)
            Mycelia.Rhizomorph.convert_to_canonical(ss)
            Test.@test StrandOracle.graph_counts(ss; undirected = false) == counts_before
            Test.@test canonical_evidence_snapshot(ss) == evidence_before
        end
    end

    Test.@testset "graph is undirected" begin
        for (_, seqs, k) in CANONICAL_HAND_CASES
            Test.@test Graphs.is_directed(build_canonical_oracle_graph(seqs, k)) == false
        end
    end

    Test.@testset "canonical GFA round-trip keeps the edge set" begin
        for (_, seqs, k) in CANONICAL_HAND_CASES
            graph = build_canonical_oracle_graph(seqs, k)
            expected_vertices, expected_edges = StrandOracle.canonical_expected(seqs, k)
            mktempdir() do dir
                path = joinpath(dir, "canonical.gfa")
                Mycelia.Rhizomorph.write_gfa_canonical(graph, path)
                kmer_type = typeof(first(MetaGraphsNext.labels(graph)))
                reread = Mycelia.Rhizomorph.read_gfa_next(path, kmer_type,
                    Mycelia.Rhizomorph.DoubleStrand)
                Test.@test Set(string.(MetaGraphsNext.labels(reread))) ==
                           Set(keys(expected_vertices))
                Test.@test unordered_label_pairs(reread) == Set(keys(expected_edges))
                Test.@test gfa_segment_depths(path) ==
                           Dict(x => Float64(c) for (x, c) in expected_vertices)
            end
        end
    end
end
