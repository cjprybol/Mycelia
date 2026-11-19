# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/oracle_qualmer_strand_test.jl")'
# ```

# Count-level brute-force oracle tests for full-evidence qualmer graphs in the
# directed double-strand and undirected canonical modes.
#
# Every testset builds a graph with
# `Mycelia.Rhizomorph.build_qualmer_graph(records, k; mode)` from FASTQ records
# carrying a constant quality string, and compares its vertex coverage and edge
# weights (evidence-entry counts) against `StrandOracle.doublestrand_expected` /
# `StrandOracle.canonical_expected`. The cases contain a self-reverse-complement
# edge (ACG -> CGT at k = 3) and a palindromic k-mer (ACGT at k = 4), which
# count twice under the shared strand convention.
#
# Only counts are compared. The orientation of the quality vectors on
# reverse-strand evidence is not asserted here.

import Test
import Mycelia
import FASTX
import MetaGraphsNext

if !isdefined(Main, :StrandOracle)
    include(joinpath(@__DIR__, "strand_oracle.jl"))
end

const QUALMER_ORACLE_CASES = [
    (["ACGT"], 3),
    (["ACGT"], 4),
    (["AACGTT"], 3),
    (["AACGTT"], 4)
]

function qualmer_oracle_records(seqs::Vector{String})::Vector{FASTX.FASTQ.Record}
    return [FASTX.FASTQ.Record("qualmer_oracle_read_$(i)", seq, repeat("I", length(seq)))
            for (i, seq) in enumerate(seqs)]
end

Test.@testset "Qualmer strand modes match brute-force oracle (counts)" begin
    for (seqs, k) in QUALMER_ORACLE_CASES
        Test.@testset "$(only(seqs)) k=$(k) :doublestrand" begin
            graph = Mycelia.Rhizomorph.build_qualmer_graph(
                qualmer_oracle_records(seqs), k; mode = :doublestrand)
            observed_vertices,
            observed_edges = StrandOracle.graph_counts(graph; undirected = false)
            expected_vertices,
            expected_edges = StrandOracle.doublestrand_expected(seqs, k)
            Test.@test observed_vertices == expected_vertices
            Test.@test observed_edges == expected_edges
        end

        Test.@testset "$(only(seqs)) k=$(k) :canonical" begin
            graph = Mycelia.Rhizomorph.build_qualmer_graph(
                qualmer_oracle_records(seqs), k; mode = :canonical)
            observed_vertices,
            observed_edges = StrandOracle.graph_counts(graph; undirected = true)
            expected_vertices, expected_edges = StrandOracle.canonical_expected(seqs, k)
            Test.@test observed_vertices == expected_vertices
            Test.@test observed_edges == expected_edges
        end
    end
end

"""
    qualmer_varying_quality_records(seqs::Vector{String})::Vector{FASTX.FASTQ.Record}

FASTQ records whose quality string varies along the read, so an in-place
reversal of a stored quality vector would change the evidence snapshot.
"""
function qualmer_varying_quality_records(
        seqs::Vector{String})::Vector{FASTX.FASTQ.Record}
    return [FASTX.FASTQ.Record("qualmer_mutation_read_$(i)", seq,
                String([Char(40 + (j % 30)) for j in 1:length(seq)]))
            for (i, seq) in enumerate(seqs)]
end

function qualmer_graph_snapshot(graph::MetaGraphsNext.MetaGraph)::Tuple{
        Set{String}, Set{Tuple{String, String}},
        Dict{String}, Dict{Tuple{String, String}}}
    vertex_labels = Set(string(label) for label in MetaGraphsNext.labels(graph))
    edge_set = Set((string(src), string(dst))
    for (src, dst) in MetaGraphsNext.edge_labels(graph))
    vertex_evidence = Dict(string(label) => deepcopy(graph[label].evidence)
    for label in MetaGraphsNext.labels(graph))
    edge_evidence = Dict((string(src), string(dst)) => deepcopy(graph[src, dst].evidence)
    for (src, dst) in MetaGraphsNext.edge_labels(graph))
    return vertex_labels, edge_set, vertex_evidence, edge_evidence
end

Test.@testset "Qualmer strand conversions do not mutate the input graph" begin
    # (["AACGTT"], 3) has the self-reverse-complement edge ACG -> CGT and
    # (["AACGTT"], 4) the palindromic k-mer ACGT, so both merge-into-self
    # branches of the converters are exercised.
    for (seqs, k) in ((["ACCGTA", "TACGGT"], 3), (["AACGTT"], 3), (["AACGTT"], 4))
        ss = Mycelia.Rhizomorph.build_qualmer_graph(
            qualmer_varying_quality_records(seqs), k; mode = :singlestrand)
        counts_before = StrandOracle.graph_counts(ss; undirected = false)
        snapshot_before = qualmer_graph_snapshot(ss)

        Mycelia.Rhizomorph.convert_to_doublestrand(ss)
        Test.@test StrandOracle.graph_counts(ss; undirected = false) == counts_before
        Test.@test qualmer_graph_snapshot(ss) == snapshot_before

        Mycelia.Rhizomorph.convert_to_canonical(ss)
        Test.@test StrandOracle.graph_counts(ss; undirected = false) == counts_before
        Test.@test qualmer_graph_snapshot(ss) == snapshot_before
    end
end
