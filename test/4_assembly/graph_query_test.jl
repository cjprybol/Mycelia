# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/graph_query_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/graph_query_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Graph Query and Traversal Functions Test
# Tests for graph query and traversal helper functions as specified
# in the rhizomorph graph ecosystem planning document.
# These functions provide efficient querying of graph properties
# and traversal of paths through the graph.

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "Graph Query Functions - Basic Properties" begin

    Test.@testset "vertex_count - Single Vertex" begin
        # Create graph with single record, then check function works
        records = [FASTX.FASTA.Record("read_001", "ATG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        count = Mycelia.Rhizomorph.vertex_count(graph)
        Test.@test count == 1  # One 3-mer: ATG
    end

    Test.@testset "vertex_count - Non-empty Graph" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCGATCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        count = Mycelia.Rhizomorph.vertex_count(graph)
        Test.@test count == 7  # 7 unique 3-mers in ATGCGATCG
    end

    Test.@testset "edge_count - No Edges" begin
        # Single k-mer has no edges
        records = [FASTX.FASTA.Record("read_001", "ATG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        count = Mycelia.Rhizomorph.edge_count(graph)
        Test.@test count == 0
    end

    Test.@testset "edge_count - Linear Path" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        count = Mycelia.Rhizomorph.edge_count(graph)
        Test.@test count == 1  # ATG -> TGC
    end

    Test.@testset "has_vertex - Existing K-mer" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, kmer)
    end

    Test.@testset "has_vertex - Non-existent K-mer" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("GGG")
        Test.@test !Mycelia.Rhizomorph.has_vertex(graph, kmer)
    end

    Test.@testset "has_edge - Existing Edge" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        src = Kmers.DNAKmer{3}("ATG")
        dst = Kmers.DNAKmer{3}("TGC")
        Test.@test Mycelia.Rhizomorph.has_edge(graph, src, dst)
    end

    Test.@testset "has_edge - Non-existent Edge" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        src = Kmers.DNAKmer{3}("ATG")
        dst = Kmers.DNAKmer{3}("GCA")
        Test.@test !Mycelia.Rhizomorph.has_edge(graph, src, dst)
    end
end

Test.@testset "Graph Query Functions - Vertex Properties" begin

    Test.@testset "get_vertex_data - Valid Vertex" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        Test.@test !isnothing(vertex_data)
        Test.@test vertex_data.Kmer == kmer
        Test.@test haskey(vertex_data.evidence, "test")
    end

    Test.@testset "get_vertex_data - Invalid Vertex" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("GGG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        Test.@test isnothing(vertex_data)
    end

    Test.@testset "get_edge_data - Valid Edge" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        src = Kmers.DNAKmer{3}("ATG")
        dst = Kmers.DNAKmer{3}("TGC")
        edge_data = Mycelia.Rhizomorph.get_edge_data(graph, src, dst)

        Test.@test !isnothing(edge_data)
        Test.@test haskey(edge_data.evidence, "test")
    end

    Test.@testset "get_edge_data - Invalid Edge" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        src = Kmers.DNAKmer{3}("ATG")
        dst = Kmers.DNAKmer{3}("GCA")
        edge_data = Mycelia.Rhizomorph.get_edge_data(graph, src, dst)

        Test.@test isnothing(edge_data)
    end

    Test.@testset "get_vertex_observation_count" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCG"),
            FASTX.FASTA.Record("read_002", "ATGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        count = Mycelia.Rhizomorph.get_vertex_observation_count(graph, kmer)

        Test.@test count == 2  # Observed in both reads
    end

    Test.@testset "get_edge_observation_count" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        src = Kmers.DNAKmer{3}("ATG")
        dst = Kmers.DNAKmer{3}("TGC")
        count = Mycelia.Rhizomorph.get_edge_observation_count(graph, src, dst)

        Test.@test count == 2  # Edge observed in both reads
    end
end

Test.@testset "Graph Traversal Functions - Neighbors" begin

    Test.@testset "get_outgoing_neighbors - Single Neighbor" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        neighbors = Mycelia.Rhizomorph.get_outgoing_neighbors(graph, kmer)

        Test.@test length(neighbors) == 1
        Test.@test Kmers.DNAKmer{3}("TGC") in neighbors
    end

    Test.@testset "get_outgoing_neighbors - Multiple Neighbors" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "ATGA")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        neighbors = Mycelia.Rhizomorph.get_outgoing_neighbors(graph, kmer)

        Test.@test length(neighbors) == 2
        Test.@test Kmers.DNAKmer{3}("TGC") in neighbors
        Test.@test Kmers.DNAKmer{3}("TGA") in neighbors
    end

    Test.@testset "get_outgoing_neighbors - No Neighbors" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        neighbors = Mycelia.Rhizomorph.get_outgoing_neighbors(graph, kmer)

        Test.@test length(neighbors) == 0
    end

    Test.@testset "get_incoming_neighbors - Single Neighbor" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("TGC")
        neighbors = Mycelia.Rhizomorph.get_incoming_neighbors(graph, kmer)

        Test.@test length(neighbors) == 1
        Test.@test Kmers.DNAKmer{3}("ATG") in neighbors
    end

    Test.@testset "get_incoming_neighbors - Multiple Neighbors" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "GTGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("TGC")
        neighbors = Mycelia.Rhizomorph.get_incoming_neighbors(graph, kmer)

        Test.@test length(neighbors) == 2
        Test.@test Kmers.DNAKmer{3}("ATG") in neighbors
        Test.@test Kmers.DNAKmer{3}("GTG") in neighbors
    end

    Test.@testset "get_outgoing_degree - Branching" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "ATGA")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        degree = Mycelia.Rhizomorph.get_outgoing_degree(graph, kmer)

        Test.@test degree == 2  # Branches to TGC and TGA
    end

    Test.@testset "get_incoming_degree - Joining" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "GTGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("TGC")
        degree = Mycelia.Rhizomorph.get_incoming_degree(graph, kmer)

        Test.@test degree == 2  # Joined from ATG and GTG
    end
end

Test.@testset "Graph Traversal Functions - Path Finding" begin

    Test.@testset "is_linear_path - True" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("TGC")
        # TGC has one incoming (ATG) and one outgoing (GCG)
        Test.@test Mycelia.Rhizomorph.is_linear_path(graph, kmer)
    end

    Test.@testset "is_linear_path - Branch Point" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "ATGA")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        # ATG has two outgoing neighbors
        Test.@test !Mycelia.Rhizomorph.is_linear_path(graph, kmer)
    end

    Test.@testset "is_source_vertex - True" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("ATG")
        # ATG has no incoming neighbors
        Test.@test Mycelia.Rhizomorph.is_source_vertex(graph, kmer)
    end

    Test.@testset "is_sink_vertex - True" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        kmer = Kmers.DNAKmer{3}("TGC")
        # TGC has no outgoing neighbors
        Test.@test Mycelia.Rhizomorph.is_sink_vertex(graph, kmer)
    end

    Test.@testset "get_all_sources" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "GGAT")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        sources = Mycelia.Rhizomorph.get_all_sources(graph)

        # ATG and GGA are source vertices (no incoming edges)
        Test.@test length(sources) == 2
        Test.@test Kmers.DNAKmer{3}("ATG") in sources
        Test.@test Kmers.DNAKmer{3}("GGA") in sources
    end

    Test.@testset "get_all_sinks" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "GGAT")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        sinks = Mycelia.Rhizomorph.get_all_sinks(graph)

        # TGC and GAT are sink vertices (no outgoing edges)
        Test.@test length(sinks) == 2
        Test.@test Kmers.DNAKmer{3}("TGC") in sinks
        Test.@test Kmers.DNAKmer{3}("GAT") in sinks
    end
end

Test.@testset "Graph Traversal Functions - Simple Paths" begin

    Test.@testset "extend_unitig_forward - Simple Path" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCGATCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        start_kmer = Kmers.DNAKmer{3}("ATG")
        path = Mycelia.Rhizomorph.extend_unitig_forward(graph, start_kmer)

        # Should extend through linear path until branch/sink
        Test.@test length(path) >= 1
        Test.@test path[1] == start_kmer
    end

    Test.@testset "extend_unitig_forward - Stop at Branch" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC"),
            FASTX.FASTA.Record("read_002", "TGCA"),
            FASTX.FASTA.Record("read_003", "TGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        start_kmer = Kmers.DNAKmer{3}("ATG")
        path = Mycelia.Rhizomorph.extend_unitig_forward(graph, start_kmer)

        # Should stop at TGC which has two outgoing neighbors (GCA and GCG)
        Test.@test path[end] == Kmers.DNAKmer{3}("TGC")
    end

    Test.@testset "extend_unitig_backward - Simple Path" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCGATCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        start_kmer = Kmers.DNAKmer{3}("TCG")
        path = Mycelia.Rhizomorph.extend_unitig_backward(graph, start_kmer)

        # Should extend backward through linear path
        Test.@test length(path) >= 1
        Test.@test path[1] == start_kmer
    end

    Test.@testset "get_maximal_unitig - Full Path" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        middle_kmer = Kmers.DNAKmer{3}("TGC")
        unitig = Mycelia.Rhizomorph.get_maximal_unitig(graph, middle_kmer)

        # Should extend both directions to get full unitig
        Test.@test length(unitig) == 3  # ATG, TGC, GCG
        Test.@test unitig[1] == Kmers.DNAKmer{3}("ATG")
        Test.@test unitig[2] == Kmers.DNAKmer{3}("TGC")
        Test.@test unitig[3] == Kmers.DNAKmer{3}("GCG")
    end

    Test.@testset "assemble_path_sequence - Reconstruct Sequence" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        path = [
            Kmers.DNAKmer{3}("ATG"),
            Kmers.DNAKmer{3}("TGC"),
            Kmers.DNAKmer{3}("GCG")
        ]

        sequence = Mycelia.Rhizomorph.assemble_path_sequence(path)

        # Overlapping k-mers should reconstruct original sequence
        Test.@test sequence == BioSequences.LongDNA{4}("ATGCG")
    end
end

Test.@testset "Graph Analysis Functions - Connected Components" begin

    Test.@testset "get_all_vertices - Single Vertex" begin
        records = [FASTX.FASTA.Record("read_001", "ATG")]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        vertices = Mycelia.Rhizomorph.get_all_vertices(graph)
        Test.@test length(vertices) == 1
    end

    Test.@testset "get_all_vertices - Non-empty Graph" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGC")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        vertices = Mycelia.Rhizomorph.get_all_vertices(graph)
        Test.@test length(vertices) == 2  # ATG and TGC
    end

    Test.@testset "filter_vertices_by_observation_count" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCGATCG"),
            FASTX.FASTA.Record("read_002", "ATGCGATCG"),
            FASTX.FASTA.Record("read_003", "GGGG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        # Filter to vertices seen at least twice
        high_coverage = Mycelia.Rhizomorph.filter_vertices_by_observation_count(graph, 2)

        # Should exclude GGG which is only seen once
        Test.@test !(Kmers.DNAKmer{3}("GGG") in high_coverage)
        Test.@test Kmers.DNAKmer{3}("ATG") in high_coverage
    end

    Test.@testset "get_graph_statistics - Basic Stats" begin
        records = [
            FASTX.FASTA.Record("read_001", "ATGCGATCG")
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(records, 3; dataset_id="test")

        stats = Mycelia.Rhizomorph.get_graph_statistics(graph)

        Test.@test haskey(stats, :vertex_count)
        Test.@test haskey(stats, :edge_count)
        Test.@test haskey(stats, :source_count)
        Test.@test haskey(stats, :sink_count)
        Test.@test stats[:vertex_count] > 0
        Test.@test stats[:edge_count] > 0
    end
end

println("✓ Graph query and traversal tests defined (will fail until implementation)")
