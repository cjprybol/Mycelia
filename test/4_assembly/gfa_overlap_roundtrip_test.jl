# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/gfa_overlap_roundtrip_test.jl")'
# ```

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX

Test.@testset "GFA link overlaps survive a write/read round-trip" begin
    # Odd exact suffix-prefix overlaps (the OLC builder only considers odd lengths):
    #   a -> b share GTTTG (5), b -> c share CCAAC (5); no other pair overlaps by >= 3.
    records = [
        FASTX.FASTA.Record("ovl_a", BioSequences.LongDNA{4}("ACGTACGTTTG")),
        FASTX.FASTA.Record("ovl_b", BioSequences.LongDNA{4}("GTTTGGCCAAC")),
        FASTX.FASTA.Record("ovl_c", BioSequences.LongDNA{4}("CCAACTAGG"))
    ]

    function overlap_lengths(graph::MetaGraphsNext.MetaGraph)::Vector{Int}
        return sort([graph[src, dst].overlap_length
                     for (src, dst) in MetaGraphsNext.edge_labels(graph)])
    end

    graph = Mycelia.Rhizomorph.build_fasta_graph(
        records; dataset_id = "gfa_overlap_roundtrip", min_overlap = 3)

    # Precondition: the fixture must produce exactly two edges, each with overlap 5.
    Test.@test Graphs.ne(graph.graph) == 2
    Test.@test overlap_lengths(graph) == [5, 5]

    mktempdir() do tmpdir
        gfa_file = joinpath(tmpdir, "overlap_roundtrip.gfa")
        Mycelia.Rhizomorph.write_gfa_next(graph, gfa_file)

        link_lines = filter(line -> startswith(line, "L\t"), readlines(gfa_file))
        Test.@test length(link_lines) == 2
        Test.@test all(split(line, '\t')[6] == "5M" for line in link_lines)

        restored = Mycelia.Rhizomorph.read_gfa_next(gfa_file)
        Test.@test Graphs.ne(restored.graph) == 2
        Test.@test overlap_lengths(restored) == [5, 5]
    end
end
