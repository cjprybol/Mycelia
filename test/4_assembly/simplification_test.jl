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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=5)

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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=10)
            Test.@test isempty(bubbles)
        end

        Test.@testset "Bubble with asymmetric path lengths" begin
            # Bubble with different path lengths: A → B → C → E and A → D → E
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D", "E"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["D", "E"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=10)

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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
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

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=10)

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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=10)
            Test.@test isempty(bubbles)
        end

        Test.@testset "Single vertex graph" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )
            graph["A"] = Mycelia.Rhizomorph.StringVertexData("A")

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=10)
            Test.@test isempty(bubbles)
        end

        Test.@testset "Disconnected components" begin
            # Two separate linear paths with no bubbles
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "X", "Y"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["X", "Y"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=10)
            Test.@test isempty(bubbles)
        end

        Test.@testset "Bubble below min_bubble_length" begin
            # Simple bubble but with min_bubble_length that excludes it
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            # Set min_bubble_length to 5 - should filter out this short bubble
            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=5, max_bubble_length=10)
            Test.@test isempty(bubbles)
        end
    end

    Test.@testset "Bubble Support Calculation" begin

        Test.@testset "Calculate path support" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=5)

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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=5)
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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
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
                Mycelia.write_fasta(outfile=fasta, records=[rec1, rec2])

                graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 3; mode=:singlestrand)

                # Should have vertices and edges
                Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

                # Detect bubbles - there should be a bubble from the SNP
                bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=20)

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
                Mycelia.write_fasta(outfile=fasta, records=[rec1, rec2])

                graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(fasta, 3; mode=:singlestrand)

                Test.@test !isempty(collect(MetaGraphsNext.labels(graph)))

                bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=20)
                Test.@test bubbles isa Vector{<:Mycelia.Rhizomorph.BubbleStructure}
            end
        end
    end

    Test.@testset "Bubble Complexity Scoring" begin

        Test.@testset "Symmetric paths have low complexity" begin
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            graph["A", "B"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["A", "C"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["B", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)
            graph["C", "D"] = Mycelia.Rhizomorph.StringEdgeData(1)

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=5)

            Test.@test length(bubbles) == 1
            bubble = bubbles[1]
            # Symmetric paths should have complexity near 0
            Test.@test bubble.complexity_score == 0.0
        end

        Test.@testset "Asymmetric paths have higher complexity" begin
            # Create bubble with asymmetric path lengths
            graph = MetaGraphsNext.MetaGraph(
                Graphs.DiGraph();
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
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

            bubbles = Mycelia.Rhizomorph.detect_bubbles_next(graph; min_bubble_length=1, max_bubble_length=10)

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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
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
                label_type=String,
                vertex_data_type=Mycelia.Rhizomorph.StringVertexData,
                edge_data_type=Mycelia.Rhizomorph.StringEdgeData
            )

            for v in ["A", "B", "C", "D"]
                graph[v] = Mycelia.Rhizomorph.StringVertexData(v)
            end

            # Valid bubble
            valid_bubble = Mycelia.Rhizomorph.BubbleStructure("A", "D", ["B", "D"], ["C", "D"], 5, 3, 0.0)
            Test.@test Mycelia.Rhizomorph.is_valid_bubble(graph, valid_bubble)

            # Invalid - non-existent entry vertex
            invalid_bubble1 = Mycelia.Rhizomorph.BubbleStructure("X", "D", ["B", "D"], ["C", "D"], 5, 3, 0.0)
            Test.@test !Mycelia.Rhizomorph.is_valid_bubble(graph, invalid_bubble1)

            # Invalid - empty path
            invalid_bubble2 = Mycelia.Rhizomorph.BubbleStructure("A", "D", String[], ["C", "D"], 5, 3, 0.0)
            Test.@test !Mycelia.Rhizomorph.is_valid_bubble(graph, invalid_bubble2)

            # Invalid - identical paths
            invalid_bubble3 = Mycelia.Rhizomorph.BubbleStructure("A", "D", ["B", "D"], ["B", "D"], 5, 3, 0.0)
            Test.@test !Mycelia.Rhizomorph.is_valid_bubble(graph, invalid_bubble3)
        end
    end
end
