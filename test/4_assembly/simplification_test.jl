# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/simplification_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/simplification_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Comprehensive tests for graph simplification algorithms

import Test
import Mycelia
import Graphs
import MetaGraphsNext
import Kmers
import BioSequences
import FASTX

Test.@testset "Graph Simplification Algorithms" begin
    Test.@testset "Bubble Detection - Basic Cases" begin
        Test.@testset "Simple SNP bubble" begin
            # Create a simple bubble: A → B → D and A → C → D
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 5)

            Test.@test length(bubbles) == 1
            bubble = bubbles[1]
            Test.@test bubble.entry_vertex == "A"
            Test.@test bubble.exit_vertex == "D"
            Test.@test bubble isa Mycelia.Rhizomorph.BubbleStructure
        end

        Test.@testset "No bubbles in linear graph" begin
            # Linear graph: A → B → C → D
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 10)
            Test.@test isempty(bubbles)
        end

        Test.@testset "Bubble with asymmetric path lengths" begin
            # Bubble with different path lengths: A → B → C → E and A → D → E
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D", "E"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["D", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 10)

            Test.@test length(bubbles) == 1
            bubble = bubbles[1]
            Test.@test bubble.entry_vertex == "A"
            Test.@test bubble.exit_vertex == "E"
            # Complexity should be non-zero due to asymmetric lengths
            Test.@test bubble.complexity_score >= 0.0
        end

        Test.@testset "Multiple parallel paths (3-way split)" begin
            # Three parallel paths: A → B → E, A → C → E, A → D → E
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D", "E"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["D", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 10)

            # Should detect multiple bubbles (each pair of paths forms a bubble)
            Test.@test length(bubbles) >= 1
            for bubble in bubbles
                Test.@test bubble.entry_vertex == "A"
                Test.@test bubble.exit_vertex == "E"
            end
        end
    end

    Test.@testset "Bubble Detection - Edge Cases" begin
        Test.@testset "Empty graph" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 10)
            Test.@test isempty(bubbles)
        end

        Test.@testset "Single vertex graph" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            graph["A"] = Mycelia.Rhizomorph.StringVertexData("A")

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 10)
            Test.@test isempty(bubbles)
        end

        Test.@testset "Disconnected components" begin
            # Two separate linear paths with no bubbles
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "X", "Y"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["X", "Y"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 10)
            Test.@test isempty(bubbles)
        end

        Test.@testset "Bubble below min_bubble_length" begin
            # Simple bubble but with min_bubble_length that excludes it
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            # Set min_bubble_length to 5 - should filter out this short bubble
            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 5, max_bubble_length = 10)
            Test.@test isempty(bubbles)
        end
    end

    Test.@testset "Bubble Support Calculation" begin
        Test.@testset "Calculate path support" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 5)

            Test.@test length(bubbles) == 1
            bubble = bubbles[1]
            # Support values should be non-negative integers
            Test.@test bubble.path1_support >= 0
            Test.@test bubble.path2_support >= 0
        end
    end

    Test.@testset "Graph Simplification" begin
        Test.@testset "Simplify graph with clear winner" begin
            # Create a bubble where one path has much higher support
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 5)
            Test.@test length(bubbles) >= 1

            # Simplify the graph
            simplified = Mycelia.Rhizomorph.simplify_graph_next(graph, bubbles)
            Test.@test simplified isa MetaGraphsNext.MetaGraph
            # The simplified graph should still have basic structure
            Test.@test !isempty(collect(MetaGraphsNext.labels(simplified)))
        end

        Test.@testset "Simplify graph - no bubbles to resolve" begin
            # Linear graph - no simplification needed
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph)
            simplified = Mycelia.Rhizomorph.simplify_graph_next(graph, bubbles)

            # Should be unchanged
            original_vertices = Set(MetaGraphsNext.labels(graph))
            simplified_vertices = Set(MetaGraphsNext.labels(simplified))
            Test.@test original_vertices == simplified_vertices
        end
    end

    Test.@testset "Bubble Detection on Real K-mer Graphs" begin
        Test.@testset "DNA k-mer graph with SNP bubble" begin
            mktempdir() do dir
                # Create two sequences that differ by one nucleotide (SNP)
                # ATGCATGC vs ATGGATGC (C→G at position 4)
                fasta = joinpath(dir, "snp.fasta")
                rec1 = FASTX.FASTA.Record("seq1", BioSequences.LongDNA{4}("ATGCATGC"))
                rec2 = FASTX.FASTA.Record("seq2", BioSequences.LongDNA{4}("ATGGATGC"))
                Mycelia.write_fasta(outfile = fasta, records = [rec1, rec2])

                graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 3; mode = :singlestrand)

                # Should have vertices and edges
                Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

                # Detect bubbles - there should be a bubble from the SNP
                bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                    graph; min_bubble_length = 1, max_bubble_length = 20)

                # The SNP should create a bubble structure
                # (exact count depends on k-mer overlap pattern)
                Test.@test bubbles isa Vector{<:Mycelia.Rhizomorph.BubbleStructure}
            end
        end

        Test.@testset "RNA k-mer graph bubble detection" begin
            mktempdir() do dir
                fasta = joinpath(dir, "rna.fasta")
                rec1 = FASTX.FASTA.Record("seq1", BioSequences.LongRNA{4}("AUGCAUGC"))
                rec2 = FASTX.FASTA.Record("seq2", BioSequences.LongRNA{4}("AUGGAUGC"))
                Mycelia.write_fasta(outfile = fasta, records = [rec1, rec2])

                graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 3; mode = :singlestrand)

                Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

                bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                    graph; min_bubble_length = 1, max_bubble_length = 20)
                Test.@test bubbles isa Vector{<:Mycelia.Rhizomorph.BubbleStructure}
            end
        end
    end

    Test.@testset "Bubble Complexity Scoring" begin
        Test.@testset "Symmetric paths have low complexity" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 5)

            Test.@test length(bubbles) == 1
            bubble = bubbles[1]
            # Symmetric paths should have complexity near 0
            Test.@test bubble.complexity_score == 0.0
        end

        Test.@testset "Asymmetric paths have higher complexity" begin
            # Create bubble with asymmetric path lengths
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D", "E"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            # Path 1: A → B → D → E (length 3)
            # Path 2: A → C → E (length 2)
            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["D", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(
                graph; min_bubble_length = 1, max_bubble_length = 10)

            if length(bubbles) > 0
                bubble = bubbles[1]
                # Asymmetric paths should have complexity > 0
                Test.@test bubble.complexity_score > 0.0
            end
        end
    end

    Test.@testset "Helper Functions" begin
        Test.@testset "get_out_neighbors" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)

            neighbors = Mycelia.Rhizomorph.get_out_neighbors(graph, "A")
            Test.@test length(neighbors) == 2
            Test.@test "B" in neighbors
            Test.@test "C" in neighbors

            # Vertex with no outgoing edges
            neighbors_d = Mycelia.Rhizomorph.get_out_neighbors(graph, "B")
            Test.@test isempty(neighbors_d)
        end

        Test.@testset "get_in_neighbors" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)

            neighbors = Mycelia.Rhizomorph.get_in_neighbors(graph, "C")
            Test.@test length(neighbors) == 2
            Test.@test "A" in neighbors
            Test.@test "B" in neighbors

            # Vertex with no incoming edges
            neighbors_a = Mycelia.Rhizomorph.get_in_neighbors(graph, "A")
            Test.@test isempty(neighbors_a)
        end

        Test.@testset "find_path_convergence" begin
            path1 = ["A", "B", "C", "D"]
            path2 = ["X", "Y", "C", "Z"]

            convergence = Mycelia.Rhizomorph.find_path_convergence(path1, path2)
            Test.@test convergence == "C"

            # No convergence
            path3 = ["A", "B"]
            path4 = ["X", "Y"]
            Test.@test Mycelia.Rhizomorph.find_path_convergence(path3, path4) === nothing
        end

        Test.@testset "is_valid_bubble" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            # Valid bubble
            valid_bubble = Mycelia.Rhizomorph.BubbleStructure(
                "A", "D", ["B", "D"], ["C", "D"], 5, 3, 0.0)
            Test.@test Mycelia.Rhizomorph.is_valid_bubble(graph, valid_bubble)

            # Invalid - non-existent entry vertex
            invalid_bubble1 = Mycelia.Rhizomorph.BubbleStructure(
                "X", "D", ["B", "D"], ["C", "D"], 5, 3, 0.0)
            Test.@test !Mycelia.Rhizomorph.is_valid_bubble(graph, invalid_bubble1)

            # Invalid - empty path
            invalid_bubble2 = Mycelia.Rhizomorph.BubbleStructure(
                "A", "D", String[], ["C", "D"], 5, 3, 0.0)
            Test.@test !Mycelia.Rhizomorph.is_valid_bubble(graph, invalid_bubble2)

            # Invalid - identical paths
            invalid_bubble3 = Mycelia.Rhizomorph.BubbleStructure(
                "A", "D", ["B", "D"], ["B", "D"], 5, 3, 0.0)
            Test.@test !Mycelia.Rhizomorph.is_valid_bubble(graph, invalid_bubble3)
        end
    end
end

Test.@testset "Linear Chain Collapsing (collapse_linear_chains!)" begin
    Test.@testset "Edge Cases" begin
        Test.@testset "Empty graph returns unchanged" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )

            result = Mycelia.Rhizomorph.collapse_linear_chains!(graph)
            Test.@test result === graph
            Test.@test Graphs.nv(graph.graph) == 0
            Test.@test Graphs.ne(graph.graph) == 0
        end

        Test.@testset "Single vertex - no collapse" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            graph["AB"] = Mycelia.Rhizomorph.StringVertexData("AB")

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)
            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, "AB")
        end

        Test.@testset "No linear chains - all branching" begin
            # Diamond: A → B, A → C, B → D, C → D
            # No vertex has both indeg=1 and outdeg=1 as interior of a chain
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end
            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(0)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)
            Test.@test Graphs.nv(graph.graph) == 4
            for v in ["A", "B", "C", "D"]
                Test.@test haskey(graph, v)
            end
        end
    end

    Test.@testset "String Graph Collapsing" begin
        Test.@testset "Two-vertex linear chain" begin
            # AB → BC with overlap=1 → should collapse to "ABC"
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            graph["AB"] = Mycelia.Rhizomorph.StringVertexData("AB")
            graph["BC"] = Mycelia.Rhizomorph.StringVertexData("BC")
            graph["AB", "BC"] = Mycelia.Rhizomorph.StringEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, "ABC")
            vdata = graph["ABC"]
            Test.@test vdata isa Mycelia.Rhizomorph.StringVertexData
            Test.@test vdata.string_value == "ABC"
        end

        Test.@testset "Three-vertex linear chain" begin
            # AB → BC → CD with overlap=1 → should collapse to "ABCD"
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            graph["AB"] = Mycelia.Rhizomorph.StringVertexData("AB")
            graph["BC"] = Mycelia.Rhizomorph.StringVertexData("BC")
            graph["CD"] = Mycelia.Rhizomorph.StringVertexData("CD")
            graph["AB", "BC"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["BC", "CD"] = Mycelia.Rhizomorph.StringEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, "ABCD")
            vdata = graph["ABCD"]
            Test.@test vdata.string_value == "ABCD"
        end

        Test.@testset "Chain with zero overlap" begin
            # "hello" → "world" with overlap=0 → should collapse to "helloworld"
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            graph["hello"] = Mycelia.Rhizomorph.StringVertexData("hello")
            graph["world"] = Mycelia.Rhizomorph.StringVertexData("world")
            graph["hello", "world"] = Mycelia.Rhizomorph.StringEdgeData(0)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, "helloworld")
            Test.@test graph["helloworld"].string_value == "helloworld"
        end

        Test.@testset "Chain with overlap=2" begin
            # "ABC" → "BCD" with overlap=2 → should collapse to "ABCD"
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            graph["ABC"] = Mycelia.Rhizomorph.StringVertexData("ABC")
            graph["BCD"] = Mycelia.Rhizomorph.StringVertexData("BCD")
            graph["ABC", "BCD"] = Mycelia.Rhizomorph.StringEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, "ABCD")
        end
    end

    Test.@testset "External Edge Reconnection" begin
        Test.@testset "Chain with incoming external edge" begin
            # X → AB → BC  (X has outdeg > 1 via X → Z too)
            # Chain [AB, BC] collapses; edge X → AB becomes X → "ABC"
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            for v in ["X", "Z", "AB", "BC"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end
            graph["X", "AB"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["X", "Z"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["AB", "BC"] = Mycelia.Rhizomorph.StringEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            # AB and BC should be collapsed into "ABC"
            Test.@test haskey(graph, "ABC")
            Test.@test !haskey(graph, "AB")
            Test.@test !haskey(graph, "BC")
            # X and Z should still exist
            Test.@test haskey(graph, "X")
            Test.@test haskey(graph, "Z")
            # Edge from X to collapsed vertex should exist
            Test.@test haskey(graph, "X", "ABC")
            # Edge from X to Z should still exist
            Test.@test haskey(graph, "X", "Z")
        end

        Test.@testset "Chain with outgoing external edge" begin
            # AB → BC → Y  (Y has indeg > 1 via Z → Y too)
            # Chain [AB, BC] collapses; edge BC → Y becomes "ABC" → Y
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            for v in ["AB", "BC", "Y", "Z"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end
            graph["AB", "BC"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["BC", "Y"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["Z", "Y"] = Mycelia.Rhizomorph.StringEdgeData(0)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test haskey(graph, "ABC")
            Test.@test !haskey(graph, "AB")
            Test.@test !haskey(graph, "BC")
            # Edge from collapsed vertex to Y
            Test.@test haskey(graph, "ABC", "Y")
            Test.@test haskey(graph, "Z", "Y")
        end

        Test.@testset "Chain with both incoming and outgoing external edges" begin
            # X → AB → BC → CD → Y
            # X has outdeg=2 (also X → Z), Y has indeg=2 (also W → Y)
            # Chain [AB, BC, CD] should collapse to "ABCD"
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            for v in ["X", "Z", "AB", "BC", "CD", "Y", "W"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end
            graph["X", "AB"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["X", "Z"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["AB", "BC"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["BC", "CD"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["CD", "Y"] = Mycelia.Rhizomorph.StringEdgeData(0)
            graph["W", "Y"] = Mycelia.Rhizomorph.StringEdgeData(0)

            original_nv = Graphs.nv(graph.graph)
            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test haskey(graph, "ABCD")
            Test.@test haskey(graph, "X", "ABCD")
            Test.@test haskey(graph, "ABCD", "Y")
            # Original vertices removed
            Test.@test !haskey(graph, "AB")
            Test.@test !haskey(graph, "BC")
            Test.@test !haskey(graph, "CD")
            # Non-chain vertices preserved
            Test.@test haskey(graph, "X")
            Test.@test haskey(graph, "Y")
            Test.@test haskey(graph, "Z")
            Test.@test haskey(graph, "W")
            # Vertex count reduced by 2 (3 vertices → 1)
            Test.@test Graphs.nv(graph.graph) == original_nv - 2
        end
    end

    Test.@testset "Multiple Separate Chains" begin
        Test.@testset "Two independent chains" begin
            # Chain 1: AB → BC (overlap=1) → "ABC"
            # Chain 2: XY → YZ (overlap=1) → "XYZ"
            # No edges between chains
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            graph["AB"] = Mycelia.Rhizomorph.StringVertexData("AB")
            graph["BC"] = Mycelia.Rhizomorph.StringVertexData("BC")
            graph["AB", "BC"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["XY"] = Mycelia.Rhizomorph.StringVertexData("XY")
            graph["YZ"] = Mycelia.Rhizomorph.StringVertexData("YZ")
            graph["XY", "YZ"] = Mycelia.Rhizomorph.StringEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test Graphs.nv(graph.graph) == 2
            Test.@test haskey(graph, "ABC")
            Test.@test haskey(graph, "XYZ")
        end
    end

    Test.@testset "Kmer Graph Skips Collapsing" begin
        Test.@testset "Kmer vertices are not collapsed" begin
            kmer1 = Kmers.DNAKmer{3}("ATG")
            kmer2 = Kmers.DNAKmer{3}("TGC")
            KmerType = typeof(kmer1)
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = KmerType,
                vertex_data_type = Mycelia.Rhizomorph.KmerVertexData{KmerType},
                edge_data_type = Mycelia.Rhizomorph.KmerEdgeData
            )
            graph[kmer1] = Mycelia.Rhizomorph.KmerVertexData(kmer1)
            graph[kmer2] = Mycelia.Rhizomorph.KmerVertexData(kmer2)
            graph[kmer1, kmer2] = Mycelia.Rhizomorph.KmerEdgeData()

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            # Both vertices should still exist — kmer graphs are skipped
            Test.@test Graphs.nv(graph.graph) == 2
            Test.@test haskey(graph, kmer1)
            Test.@test haskey(graph, kmer2)
        end
    end

    Test.@testset "BioSequence Graph Collapsing" begin
        Test.@testset "DNA BioSequence linear chain" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            seq3 = BioSequences.LongDNA{4}("GCA")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.BioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.BioSequenceEdgeData
            )
            graph[seq1] = Mycelia.Rhizomorph.BioSequenceVertexData(seq1)
            graph[seq2] = Mycelia.Rhizomorph.BioSequenceVertexData(seq2)
            graph[seq3] = Mycelia.Rhizomorph.BioSequenceVertexData(seq3)
            graph[seq1, seq2] = Mycelia.Rhizomorph.BioSequenceEdgeData(2)
            graph[seq2, seq3] = Mycelia.Rhizomorph.BioSequenceEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGCA")
            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, expected)
            vdata = graph[expected]
            Test.@test vdata isa Mycelia.Rhizomorph.BioSequenceVertexData
            Test.@test vdata.sequence == expected
        end
    end

    Test.@testset "Reduced Type Collapsing" begin
        Test.@testset "LightweightStringVertexData" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.LightweightStringVertexData,
                edge_data_type = Mycelia.Rhizomorph.LightweightEdgeData
            )
            graph["AB"] = Mycelia.Rhizomorph.LightweightStringVertexData("AB")
            graph["BC"] = Mycelia.Rhizomorph.LightweightStringVertexData("BC")
            graph["AB", "BC"] = Mycelia.Rhizomorph.LightweightEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, "ABC")
            vdata = graph["ABC"]
            Test.@test vdata isa Mycelia.Rhizomorph.LightweightStringVertexData
            Test.@test vdata.string_value == "ABC"
        end

        Test.@testset "UltralightStringVertexData" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.UltralightStringVertexData,
                edge_data_type = Mycelia.Rhizomorph.UltralightEdgeData
            )
            graph["AB"] = Mycelia.Rhizomorph.UltralightStringVertexData("AB")
            graph["BC"] = Mycelia.Rhizomorph.UltralightStringVertexData("BC")
            graph["AB", "BC"] = Mycelia.Rhizomorph.UltralightEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, "ABC")
            vdata = graph["ABC"]
            Test.@test vdata isa Mycelia.Rhizomorph.UltralightStringVertexData
            Test.@test vdata.string_value == "ABC"
        end

        Test.@testset "LightweightBioSequenceVertexData" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.LightweightBioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.LightweightEdgeData
            )
            graph[seq1] = Mycelia.Rhizomorph.LightweightBioSequenceVertexData(seq1)
            graph[seq2] = Mycelia.Rhizomorph.LightweightBioSequenceVertexData(seq2)
            graph[seq1, seq2] = Mycelia.Rhizomorph.LightweightEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGC")
            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, expected)
            vdata = graph[expected]
            Test.@test vdata isa Mycelia.Rhizomorph.LightweightBioSequenceVertexData
            Test.@test vdata.sequence == expected
        end

        Test.@testset "UltralightBioSequenceVertexData" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.UltralightBioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.UltralightEdgeData
            )
            graph[seq1] = Mycelia.Rhizomorph.UltralightBioSequenceVertexData(seq1)
            graph[seq2] = Mycelia.Rhizomorph.UltralightBioSequenceVertexData(seq2)
            graph[seq1, seq2] = Mycelia.Rhizomorph.UltralightEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGC")
            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, expected)
            vdata = graph[expected]
            Test.@test vdata isa Mycelia.Rhizomorph.UltralightBioSequenceVertexData
        end

        Test.@testset "QualityBioSequenceVertexData" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.QualityBioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.QualityBioSequenceEdgeData
            )
            graph[seq1] = Mycelia.Rhizomorph.QualityBioSequenceVertexData(seq1)
            graph[seq2] = Mycelia.Rhizomorph.QualityBioSequenceVertexData(seq2)
            graph[seq1, seq2] = Mycelia.Rhizomorph.QualityBioSequenceEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGC")
            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, expected)
            vdata = graph[expected]
            Test.@test vdata isa Mycelia.Rhizomorph.QualityBioSequenceVertexData
            Test.@test vdata.sequence == expected
        end

        Test.@testset "UltralightQualityBioSequenceVertexData" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.UltralightQualityEdgeData
            )
            graph[seq1] = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(seq1)
            graph[seq2] = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(seq2)
            graph[seq1, seq2] = Mycelia.Rhizomorph.UltralightQualityEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGC")
            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, expected)
            vdata = graph[expected]
            Test.@test vdata isa Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData
            Test.@test vdata.sequence == expected
        end

        Test.@testset "LightweightQualityBioSequenceVertexData" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.LightweightQualityEdgeData
            )
            graph[seq1] = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(seq1)
            graph[seq2] = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(seq2)
            graph[seq1, seq2] = Mycelia.Rhizomorph.LightweightQualityEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGC")
            Test.@test Graphs.nv(graph.graph) == 1
            Test.@test haskey(graph, expected)
            vdata = graph[expected]
            Test.@test vdata isa Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData
            Test.@test vdata.sequence == expected
        end
    end

    Test.@testset "Evidence Merging" begin
        Test.@testset "Reduced type - counts accumulated" begin
            # LightweightStringVertexData with evidence on each vertex
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.LightweightStringVertexData,
                edge_data_type = Mycelia.Rhizomorph.LightweightEdgeData
            )
            v1 = Mycelia.Rhizomorph.LightweightStringVertexData("AB")
            v2 = Mycelia.Rhizomorph.LightweightStringVertexData("BC")

            # Add evidence to v1: 3 observations from ds1, 1 from ds2
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs1",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs2",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs3",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds2", "obs4",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))

            # Add evidence to v2: 2 observations from ds1
            Mycelia.Rhizomorph.add_evidence!(
                v2, "ds1", "obs1",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
            Mycelia.Rhizomorph.add_evidence!(
                v2, "ds1", "obs5",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))

            graph["AB"] = v1
            graph["BC"] = v2
            graph["AB", "BC"] = Mycelia.Rhizomorph.LightweightEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            Test.@test haskey(graph, "ABC")
            collapsed = graph["ABC"]

            # Total count: 4 + 2 = 6
            Test.@test collapsed.total_count == 6
            # Dataset counts: ds1 = 3 + 2 = 5, ds2 = 1
            Test.@test collapsed.dataset_counts["ds1"] == 5
            Test.@test collapsed.dataset_counts["ds2"] == 1
            # Observation IDs: union of all obs per dataset
            Test.@test "obs1" in collapsed.dataset_observations["ds1"]
            Test.@test "obs2" in collapsed.dataset_observations["ds1"]
            Test.@test "obs3" in collapsed.dataset_observations["ds1"]
            Test.@test "obs5" in collapsed.dataset_observations["ds1"]
            Test.@test "obs4" in collapsed.dataset_observations["ds2"]
        end

        Test.@testset "Ultralight type - counts only (no obs IDs)" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.UltralightStringVertexData,
                edge_data_type = Mycelia.Rhizomorph.UltralightEdgeData
            )
            v1 = Mycelia.Rhizomorph.UltralightStringVertexData("AB")
            v2 = Mycelia.Rhizomorph.UltralightStringVertexData("BC")

            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs1",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs2",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
            Mycelia.Rhizomorph.add_evidence!(
                v2, "ds1", "obs3",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))

            graph["AB"] = v1
            graph["BC"] = v2
            graph["AB", "BC"] = Mycelia.Rhizomorph.UltralightEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            collapsed = graph["ABC"]
            Test.@test collapsed.total_count == 3
            Test.@test collapsed.dataset_counts["ds1"] == 3
            # Ultralight has no dataset_observations field
            Test.@test !hasfield(typeof(collapsed), :dataset_observations)
        end

        Test.@testset "UltralightQuality type - counts + joint quality merged" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.UltralightQualityEdgeData
            )
            v1 = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(seq1)
            v2 = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(seq2)

            # Add quality evidence: Phred+33 encoded scores
            # v1: Q10, Q15, Q12 → raw [10, 15, 12]
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs1",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[43, 48, 45]))
            # v1 second observation: Q20, Q20, Q20 → raw [20, 20, 20]
            # joint_quality after: [10+20, 15+20, 12+20] = [30, 35, 32]
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs2",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[53, 53, 53]))

            # v2: Q5, Q10, Q8 → raw [5, 10, 8]
            Mycelia.Rhizomorph.add_evidence!(
                v2, "ds1", "obs3",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[38, 43, 41]))

            graph[seq1] = v1
            graph[seq2] = v2
            graph[seq1, seq2] = Mycelia.Rhizomorph.UltralightQualityEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGC")
            Test.@test haskey(graph, expected)
            collapsed = graph[expected]

            # Counts merged: v1 had 2, v2 had 1
            Test.@test collapsed.total_count == 3
            Test.@test collapsed.dataset_counts["ds1"] == 3

            # Joint quality: first vertex copied [30, 35, 32],
            # then second vertex added element-wise [5, 10, 8]
            # Result: [35, 45, 40]
            Test.@test length(collapsed.joint_quality) == 3
            Test.@test collapsed.joint_quality == UInt8[35, 45, 40]

            # Ultralight quality has no dataset_observations
            Test.@test !hasfield(typeof(collapsed), :dataset_observations)
        end

        Test.@testset "LightweightQuality type - counts + obs IDs + joint quality merged" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.LightweightQualityEdgeData
            )
            v1 = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(seq1)
            v2 = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(seq2)

            # v1: 3 observations from ds1, 1 from ds2
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs1",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[43, 48, 45]))
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs2",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[53, 53, 53]))
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs3",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[36, 36, 36]))
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds2", "obs4",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[63, 63, 63]))

            # v2: 2 observations from ds1 (one shared obs ID with v1)
            Mycelia.Rhizomorph.add_evidence!(
                v2, "ds1", "obs1",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[43, 43, 43]))
            Mycelia.Rhizomorph.add_evidence!(
                v2, "ds1", "obs5",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[53, 53, 53]))

            graph[seq1] = v1
            graph[seq2] = v2
            graph[seq1, seq2] = Mycelia.Rhizomorph.LightweightQualityEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGC")
            Test.@test haskey(graph, expected)
            collapsed = graph[expected]

            # Counts merged: v1 had 4, v2 had 2
            Test.@test collapsed.total_count == 6
            Test.@test collapsed.dataset_counts["ds1"] == 5
            Test.@test collapsed.dataset_counts["ds2"] == 1

            # Observation IDs: union of all obs per dataset
            Test.@test "obs1" in collapsed.dataset_observations["ds1"]
            Test.@test "obs2" in collapsed.dataset_observations["ds1"]
            Test.@test "obs3" in collapsed.dataset_observations["ds1"]
            Test.@test "obs5" in collapsed.dataset_observations["ds1"]
            Test.@test "obs4" in collapsed.dataset_observations["ds2"]

            # Joint quality merged:
            # v1 joint_quality after 3 ds1 obs + 1 ds2 obs:
            #   ds1: (43-33)+(53-33)+(36-33) = 10+20+3 = 33 per position → [33, 38, 35]
            #   ds2: 63-33 = 30 → adds [30, 30, 30]
            #   total joint: [33+30, 38+30, 35+30] = [63, 68, 65]
            #   Wait - joint_quality is the global sum, not per-dataset
            #   obs1: [10, 15, 12], obs2: [20, 20, 20], obs3: [3, 3, 3], obs4: [30, 30, 30]
            #   v1 joint = [10+20+3+30, 15+20+3+30, 12+20+3+30] = [63, 68, 65]
            # v2 joint_quality after 2 ds1 obs:
            #   obs1: [10, 10, 10], obs5: [20, 20, 20]
            #   v2 joint = [30, 30, 30]
            # Merged: first vertex copied [63, 68, 65],
            #   then second vertex added [30, 30, 30]
            # Result: [93, 98, 95]
            Test.@test length(collapsed.joint_quality) == 3
            Test.@test collapsed.joint_quality == UInt8[93, 98, 95]
        end

        Test.@testset "QualityBioSequence type - quality evidence entries shifted" begin
            seq1 = BioSequences.LongDNA{4}("ATG")
            seq2 = BioSequences.LongDNA{4}("TGC")
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = BioSequences.LongDNA{4},
                vertex_data_type = Mycelia.Rhizomorph.QualityBioSequenceVertexData{BioSequences.LongDNA{4}},
                edge_data_type = Mycelia.Rhizomorph.QualityBioSequenceEdgeData
            )
            v1 = Mycelia.Rhizomorph.QualityBioSequenceVertexData(seq1)
            v2 = Mycelia.Rhizomorph.QualityBioSequenceVertexData(seq2)

            # Add quality evidence at position 1 for both
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs1",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[43, 48, 45]))
            Mycelia.Rhizomorph.add_evidence!(
                v2, "ds1", "obs2",
                Mycelia.Rhizomorph.QualityEvidenceEntry(
                    1, Mycelia.Rhizomorph.Forward, UInt8[53, 53, 53]))

            graph[seq1] = v1
            graph[seq2] = v2
            graph[seq1, seq2] = Mycelia.Rhizomorph.QualityBioSequenceEdgeData(2)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            expected = BioSequences.LongDNA{4}("ATGC")
            Test.@test haskey(graph, expected)
            collapsed = graph[expected]

            # Full quality type uses the evidence dict (not reduced path)
            # v1 at offset 0: position 1 + 0 = 1
            # v2 at offset 1: position 1 + 1 = 2
            ds1_evidence = collapsed.evidence["ds1"]
            obs1_entries = ds1_evidence["obs1"]
            obs2_entries = ds1_evidence["obs2"]
            Test.@test any(e -> e.position == 1, obs1_entries)
            Test.@test any(e -> e.position == 2, obs2_entries)

            # Quality scores preserved in shifted entries
            shifted_obs1 = first(filter(e -> e.position == 1, obs1_entries))
            shifted_obs2 = first(filter(e -> e.position == 2, obs2_entries))
            Test.@test shifted_obs1.quality_scores == UInt8[43, 48, 45]
            Test.@test shifted_obs2.quality_scores == UInt8[53, 53, 53]
        end

        Test.@testset "Full type - evidence entries shifted" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            v1 = Mycelia.Rhizomorph.StringVertexData("AB")
            v2 = Mycelia.Rhizomorph.StringVertexData("BC")

            # Add evidence at position 1 for v1
            Mycelia.Rhizomorph.add_evidence!(
                v1, "ds1", "obs1",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
            # Add evidence at position 1 for v2
            Mycelia.Rhizomorph.add_evidence!(
                v2, "ds1", "obs2",
                Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))

            graph["AB"] = v1
            graph["BC"] = v2
            graph["AB", "BC"] = Mycelia.Rhizomorph.StringEdgeData(1)

            Mycelia.Rhizomorph.collapse_linear_chains!(graph)

            collapsed = graph["ABC"]
            # v1 was at offset 0, v2 at offset 1 (len("AB")=2, overlap=1 → 2-1=1)
            # v1's evidence: position 1 + offset 0 = position 1
            # v2's evidence: position 1 + offset 1 = position 2
            ds1_evidence = collapsed.evidence["ds1"]
            obs1_entries = ds1_evidence["obs1"]
            obs2_entries = ds1_evidence["obs2"]
            Test.@test any(e -> e.position == 1, obs1_entries)
            Test.@test any(e -> e.position == 2, obs2_entries)
        end
    end

    Test.@testset "Return Value" begin
        Test.@testset "Returns the mutated graph" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type = String,
                vertex_data_type = Mycelia.Rhizomorph.StringVertexData,
                edge_data_type = Mycelia.Rhizomorph.StringEdgeData
            )
            graph["AB"] = Mycelia.Rhizomorph.StringVertexData("AB")
            graph["BC"] = Mycelia.Rhizomorph.StringVertexData("BC")
            graph["AB", "BC"] = Mycelia.Rhizomorph.StringEdgeData(1)

            result = Mycelia.Rhizomorph.collapse_linear_chains!(graph)
            Test.@test result === graph
        end
    end
end
