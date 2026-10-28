# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/oracle_doublestrand_test.jl")'
# ```

# Brute-force oracle tests for directed double-strand k-mer graphs.
#
# Every testset builds a graph with
# `Mycelia.Rhizomorph.build_kmer_graph(records, k; mode = :doublestrand)` and
# compares its vertex coverage and edge weights (evidence-entry counts) against
# `StrandOracle.doublestrand_expected`, which is computed with plain String
# operations only.
#
# Self-reverse-complement edges (b == rc(a), e.g. ACG -> CGT) are counted twice
# here: the forward observation and its reverse-complement image land on the
# same directed edge, so its weight is 2 * fwdE(a, b). In the canonical
# (undirected) oracle the same edge keeps weight fwdE(a, b), because the
# canonical graph records each forward observation once. Palindromic k-mers are
# likewise counted on both strands here (coverage 2 * fwd(x)).

import Test
import Random
import Mycelia
import FASTX
import MetaGraphsNext

if !isdefined(Main, :StrandOracle)
    include(joinpath(@__DIR__, "strand_oracle.jl"))
end

function oracle_records(seqs::Vector{String})::Vector{FASTX.FASTA.Record}
    return [FASTX.FASTA.Record("oracle_read_$(i)", seq) for (i, seq) in enumerate(seqs)]
end

function doublestrand_matches_oracle(seqs::Vector{String}, k::Int)::Bool
    graph = Mycelia.Rhizomorph.build_kmer_graph(oracle_records(seqs), k;
        mode = :doublestrand)
    observed_vertices,
    observed_edges = StrandOracle.graph_counts(graph; undirected = false)
    expected_vertices, expected_edges = StrandOracle.doublestrand_expected(seqs, k)
    return observed_vertices == expected_vertices && observed_edges == expected_edges
end

function test_doublestrand_oracle(seqs::Vector{String}, k::Int)
    graph = Mycelia.Rhizomorph.build_kmer_graph(oracle_records(seqs), k;
        mode = :doublestrand)
    observed_vertices,
    observed_edges = StrandOracle.graph_counts(graph; undirected = false)
    expected_vertices, expected_edges = StrandOracle.doublestrand_expected(seqs, k)
    Test.@test observed_vertices == expected_vertices
    Test.@test observed_edges == expected_edges
    return nothing
end

function evidence_snapshot(graph)
    vertex_evidence = Dict(string(label) => deepcopy(graph[label].evidence)
    for label in MetaGraphsNext.labels(graph))
    edge_evidence = Dict((string(src), string(dst)) => deepcopy(graph[src, dst].evidence)
    for (src, dst) in MetaGraphsNext.edge_labels(graph))
    return vertex_evidence, edge_evidence
end

Test.@testset "Double-strand k-mer graph matches brute-force oracle" begin
    Test.@testset "hand: self-reverse-complement edge (ACG→CGT)" begin
        test_doublestrand_oracle(["ACGTTACGT"], 3)
    end

    Test.@testset "hand: distinct edge and RC-edge pair" begin
        test_doublestrand_oracle(["AACGTT"], 3)
    end

    Test.@testset "hand: palindromic k-mer at even k" begin
        test_doublestrand_oracle(["AACGTT"], 4)
    end

    Test.@testset "hand: two records sharing k-mers across strands" begin
        test_doublestrand_oracle(["ACCGTA", "TACGGT"], 3)
    end

    Test.@testset "randomized, 200 cases" begin
        rng = Random.MersenneTwister(20261026)
        failures = Tuple{Vector{String}, Int}[]
        for i in 1:200
            k = rand(rng, [3, 4, 5, 7, 11])
            n = rand(rng, 20:200)
            n_seqs = rand(rng, 1:2)
            seqs = [StrandOracle.random_non_rc_closed(rng, n, k) for _ in 1:n_seqs]
            if !doublestrand_matches_oracle(seqs, k)
                push!(failures, (seqs, k))
            end
        end
        if !isempty(failures)
            @info "randomized double-strand oracle failures" n_failures=length(failures) first_failure=first(failures)
        end
        Test.@test isempty(failures)
    end

    Test.@testset "phiX174 reference" begin
        genome_info = Mycelia.get_test_genome_fasta(use_ncbi = true)
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
                test_doublestrand_oracle(seqs, k)
            end
        end
        genome_info.cleanup()
    end

    Test.@testset "input single-strand graph is not mutated" begin
        for seqs in (["ACCGTA", "TACGGT"], ["AACGTT"])
            ss = Mycelia.Rhizomorph.build_kmer_graph(oracle_records(seqs), 3;
                mode = :singlestrand)
            counts_before = StrandOracle.graph_counts(ss; undirected = false)
            evidence_before = evidence_snapshot(ss)
            Mycelia.Rhizomorph.convert_to_doublestrand(ss)
            Test.@test StrandOracle.graph_counts(ss; undirected = false) == counts_before
            Test.@test evidence_snapshot(ss) == evidence_before
        end
    end
end
