# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/gfa_io_next.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/gfa_io_next.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences
import FASTX
import Kmers

Test.@testset "GFA I/O Next-Generation Tests (Rhizomorph)" begin
    Test.@testset "GFA Writing" begin
        records = [FASTX.FASTA.Record("sample", "ATCG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records,
            3;
            dataset_id="gfa_write",
            mode=:doublestrand,
        )

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "test.gfa")
            result_file = Mycelia.Rhizomorph.write_gfa_next(graph, gfa_file)

            Test.@test result_file == gfa_file
            Test.@test isfile(gfa_file)

            content = read(gfa_file, String)
            lines = split(strip(content), '\n')

            Test.@test startswith(lines[1], "H\t")
            segment_lines = filter(line -> startswith(line, "S\t"), lines)
            Test.@test !isempty(segment_lines)
            Test.@test all(occursin("DP:f:", line) for line in segment_lines)

            link_lines = filter(line -> startswith(line, "L\t"), lines)
            Test.@test !isempty(link_lines)
            link_fields = split(first(link_lines), '\t')
            Test.@test link_fields[3] == "+"
            Test.@test link_fields[5] == "+"
            Test.@test link_fields[6] == "2M"
        end
    end

    Test.@testset "GFA Reading" begin
        gfa_content = """
        H	VN:Z:1.0
        S	1	ATC	DP:f:2.0
        S	2	TCG	DP:f:1.0
        L	1	+	2	+	2M
        """

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "test_input.gfa")
            write(gfa_file, gfa_content)

            graph = Mycelia.Rhizomorph.read_gfa_next(
                gfa_file,
                Kmers.DNAKmer{3},
                Mycelia.Rhizomorph.DoubleStrand,
            )

            Test.@test graph isa MetaGraphsNext.MetaGraph
            labels = collect(MetaGraphsNext.labels(graph))
            expected_labels = Set([Kmers.DNAKmer{3}("ATC"), Kmers.DNAKmer{3}("TCG")])
            Test.@test Set(labels) == expected_labels

            edge_data = graph[Kmers.DNAKmer{3}("ATC"), Kmers.DNAKmer{3}("TCG")]
            Test.@test edge_data isa Mycelia.Rhizomorph.KmerEdgeData
        end
    end

    Test.@testset "GFA Round-trip" begin
        seq1 = FASTX.FASTA.Record("test1", "ATCG")
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            [seq1],
            3;
            dataset_id="gfa_rt",
            mode=:doublestrand,
        )

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "roundtrip.gfa")
            Mycelia.Rhizomorph.write_gfa_next(graph, gfa_file)
            restored = Mycelia.Rhizomorph.read_gfa_next(
                gfa_file,
                Kmers.DNAKmer{3},
                Mycelia.Rhizomorph.DoubleStrand,
            )

            Test.@test length(MetaGraphsNext.labels(restored)) == length(MetaGraphsNext.labels(graph))
            Test.@test length(collect(MetaGraphsNext.edge_labels(restored))) == length(collect(MetaGraphsNext.edge_labels(graph)))
        end
    end

    Test.@testset "Graph Modes in GFA I/O" begin
        gfa_content = """
        H	VN:Z:1.0
        S	1	ATC
        S	2	GAT
        L	1	+	2	+	2M
        """

        mktempdir() do tmpdir
            gfa_file = joinpath(tmpdir, "modes_test.gfa")
            write(gfa_file, gfa_content)

            single_graph = Mycelia.Rhizomorph.read_gfa_next(
                gfa_file,
                Kmers.DNAKmer{3},
                Mycelia.Rhizomorph.SingleStrand,
            )
            double_graph = Mycelia.Rhizomorph.read_gfa_next(
                gfa_file,
                Kmers.DNAKmer{3},
                Mycelia.Rhizomorph.DoubleStrand,
            )

            for graph in (single_graph, double_graph)
                Test.@test graph isa MetaGraphsNext.MetaGraph
                Test.@test all(label isa Kmers.DNAKmer{3} for label in MetaGraphsNext.labels(graph))
            end
        end
    end

    Test.@testset "Error Handling" begin
        Test.@test_throws SystemError Mycelia.Rhizomorph.read_gfa_next("nonexistent.gfa", Kmers.DNAKmer{3}, Mycelia.Rhizomorph.DoubleStrand)

        empty_graph = MetaGraphsNext.MetaGraph(
            MetaGraphsNext.DiGraph();
            label_type=Kmers.DNAKmer{3},
            vertex_data_type=Mycelia.Rhizomorph.KmerVertexData{Kmers.DNAKmer{3}},
            edge_data_type=Mycelia.Rhizomorph.KmerEdgeData,
        )

        mktempdir() do tmpdir
            empty_gfa = joinpath(tmpdir, "empty.gfa")
            result = Mycelia.Rhizomorph.write_gfa_next(empty_graph, empty_gfa)
            Test.@test result == empty_gfa
            Test.@test isfile(empty_gfa)

            content = read(empty_gfa, String)
            lines = filter(!isempty, split(strip(content), '\n'))
            Test.@test length(lines) == 1
            Test.@test startswith(lines[1], "H\t")
        end
    end
end
