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
