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
# One convention for both strand modes: an observation that is its own reverse
# complement counts twice, once per strand. A self-reverse-complement edge
# (b == rc(a), e.g. ACG -> CGT) has weight 2 * fwdE(a, b), because the forward
# observation and its reverse-complement image land on the same directed edge,
# and a palindromic k-mer has coverage 2 * fwd(x). The canonical oracle applies
# the same rule to the self-loop (canon(a), canon(a)) and to palindromes.

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

function test_doublestrand_oracle(seqs::Vector{String}, k::Int)::Nothing
    graph = Mycelia.Rhizomorph.build_kmer_graph(oracle_records(seqs), k;
        mode = :doublestrand)
    observed_vertices,
    observed_edges = StrandOracle.graph_counts(graph; undirected = false)
    expected_vertices, expected_edges = StrandOracle.doublestrand_expected(seqs, k)
    Test.@test observed_vertices == expected_vertices
    Test.@test observed_edges == expected_edges
    return nothing
end

function evidence_snapshot(graph::MetaGraphsNext.MetaGraph)::Tuple{
        Dict{String}, Dict{Tuple{String, String}}}
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

    Test.@testset "seeded random 5 kb genome" begin
        # Always runs (no network), at assembly-scale k. A single random read
        # never observes both strands, odd k has no palindromic k-mers, and a
        # self-reverse-complement edge is ~1e-3 likely, so two pieces are added
        # to make the case able to fail on the strand-counting defects:
        #   - a palindromic (k + 1)-mer p * rc(p), length(p) = (k + 1) / 2,
        #     spliced mid-sequence: its two k-mers form a guaranteed
        #     self-reverse-complement edge (b == rc(a)), exercising the
        #     "counts twice" rule for edges;
        #   - a second read rc(seq[1001:2000]), so both strands of a 1 kb region
        #     are observed, exercising summation of an observed edge/vertex with
        #     its observed reverse-complement image.
        # The remaining sequence stays non-RC-closed (fixed seeds).
        for k in (21, 31)
            rng = Random.MersenneTwister(20261030)
            base = StrandOracle.random_non_rc_closed(rng, 5000, k)
            p = String(rand(Random.MersenneTwister(20261031 + k), collect("ACGT"),
                (k + 1) ÷ 2))
            palindrome = p * StrandOracle.revcomp(p)
            Test.@test length(palindrome) == k + 1
            Test.@test StrandOracle.revcomp(palindrome) == palindrome
            seq = base[1:2500] * palindrome * base[2501:end]
            reads = [seq, StrandOracle.revcomp(seq[1001:2000])]
            test_doublestrand_oracle(reads, k)
        end
    end

    # When the NCBI download is unavailable (e.g. CI without the
    # ncbi-datasets-cli environment), get_test_genome_fasta returns a simulated
    # genome of the same length; the oracle comparison runs on whichever genome
    # was returned, and the testset name records which one (:ncbi = real phiX).
    genome_info = Mycelia.get_test_genome_fasta(use_ncbi = true)
    Test.@testset "phiX174 genome ($(genome_info.source))" begin
        try
            @info "double-strand oracle genome source" source=genome_info.source fasta=genome_info.fasta
            seqs = String[]
            for record in Mycelia.open_fastx(genome_info.fasta)
                push!(seqs, uppercase(String(FASTX.sequence(record))))
            end
            Test.@test length(seqs) == 1
            Test.@test all(seq -> all(c -> c in "ACGT", seq), seqs)
            for k in (11, 21)
                test_doublestrand_oracle(seqs, k)
            end
        finally
            genome_info.cleanup()
        end
    end

    Test.@testset "input single-strand graph is not mutated" begin
        # (["AACGTT"], 4) has the palindromic k-mer ACGT, so its reverse-
        # complement evidence merges into the copied forward vertex.
        for (seqs, k) in ((["ACCGTA", "TACGGT"], 3), (["AACGTT"], 3), (["AACGTT"], 4))
            ss = Mycelia.Rhizomorph.build_kmer_graph(oracle_records(seqs), k;
                mode = :singlestrand)
            counts_before = StrandOracle.graph_counts(ss; undirected = false)
            evidence_before = evidence_snapshot(ss)
            Mycelia.Rhizomorph.convert_to_doublestrand(ss)
            Test.@test StrandOracle.graph_counts(ss; undirected = false) == counts_before
            Test.@test evidence_snapshot(ss) == evidence_before
        end
    end
end
