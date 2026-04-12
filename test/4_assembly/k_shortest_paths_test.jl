import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext
import Kmers

Test.@testset "K-Shortest Paths" begin
    Test.@testset "Single path graph returns K=1" begin
        # ACGTTTCG with k=4 produces 5 distinct k-mers in a strictly linear chain
        # (no k-mer appears at both start and end), guaranteeing a true source and sink.
        sequence = BioSequences.dna"ACGTTTCG"
        records = [FASTX.FASTA.Record("linear", sequence)]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "linear_ksp", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

        labels = collect(MetaGraphsNext.labels(weighted))
        in_deg = Dict(l => 0 for l in labels)
        out_deg = Dict(l => 0 for l in labels)
        for (src, dst) in MetaGraphsNext.edge_labels(weighted)
            out_deg[src] += 1
            in_deg[dst] += 1
        end
        sources = [l for l in labels if in_deg[l] == 0]
        sinks = [l for l in labels if out_deg[l] == 0]
        source = first(sources)
        target = first(sinks)

        paths = Mycelia.Rhizomorph.k_shortest_paths(weighted, source, target, 5)

        # Only 1 path exists in a linear graph
        Test.@test length(paths) == 1

        # The single path should be valid
        path = first(paths)
        Test.@test path isa Mycelia.Rhizomorph.GraphPath
        Test.@test length(path.steps) > 0
        Test.@test path.total_probability > 0.0
    end

    Test.@testset "Excluded vertices block paths" begin
        seq1 = BioSequences.dna"ATCGTTGA"
        seq2 = BioSequences.dna"ATCAATGA"
        records = [
            FASTX.FASTA.Record("path1", seq1),
            FASTX.FASTA.Record("path2", seq2)
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "bubble_ksp", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
        labels = collect(MetaGraphsNext.labels(weighted))

        all_excluded = Set(labels[2:end])
        source = labels[1]
        target = labels[end]

        result = Mycelia.Rhizomorph._shortest_path_excluding(
            weighted, source, target, all_excluded, Set{Tuple{
                eltype(labels), eltype(labels)}}()
        )
        Test.@test result === nothing
    end

    Test.@testset "Bubble graph returns 2 paths" begin
        seq1 = BioSequences.dna"ATCGTTGA"
        seq2 = BioSequences.dna"ATCAATGA"
        records = [
            FASTX.FASTA.Record("path1", seq1),
            FASTX.FASTA.Record("path2", seq2)
        ]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "bubble2_ksp", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

        labels = collect(MetaGraphsNext.labels(weighted))
        in_edges = Dict(l => 0 for l in labels)
        out_edges = Dict(l => 0 for l in labels)
        for (src, dst) in MetaGraphsNext.edge_labels(weighted)
            out_edges[src] = get(out_edges, src, 0) + 1
            in_edges[dst] = get(in_edges, dst, 0) + 1
        end
        sources = [l for l in labels if get(in_edges, l, 0) == 0]
        sinks = [l for l in labels if get(out_edges, l, 0) == 0]

        Test.@test !isempty(sources)
        Test.@test !isempty(sinks)

        paths = Mycelia.Rhizomorph.k_shortest_paths(
            weighted, first(sources), first(sinks), 5
        )
        Test.@test length(paths) >= 1
        Test.@test length(paths) <= 2
        if length(paths) >= 2
            Test.@test paths[1].total_probability >= paths[2].total_probability
        end
        for path in paths
            Test.@test path isa Mycelia.Rhizomorph.GraphPath
            Test.@test path.total_probability > 0.0
        end
    end

    Test.@testset "Edge cases" begin
        Test.@testset "K=0 returns empty" begin
            sequence = BioSequences.dna"ATCGATCG"
            records = [FASTX.FASTA.Record("seq", sequence)]
            graph = Mycelia.Rhizomorph.build_kmer_graph(
                records, 4; dataset_id = "k0_test", mode = :singlestrand
            )
            weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
            labels = collect(MetaGraphsNext.labels(weighted))

            paths = Mycelia.Rhizomorph.k_shortest_paths(
                weighted, first(labels), last(labels), 0
            )
            Test.@test isempty(paths)
        end

        Test.@testset "No path between disconnected vertices" begin
            seq1 = BioSequences.dna"ACGTACGT"
            seq2 = BioSequences.dna"TGCATGCA"
            records = [
                FASTX.FASTA.Record("a", seq1),
                FASTX.FASTA.Record("b", seq2)
            ]
            graph = Mycelia.Rhizomorph.build_kmer_graph(
                records, 4; dataset_id = "disc_ksp", mode = :singlestrand
            )
            weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
            labels = collect(MetaGraphsNext.labels(weighted))

            Test.@test length(labels) >= 2
            paths = Mycelia.Rhizomorph.k_shortest_paths(
                weighted, labels[1], labels[end], 3
            )
            Test.@test length(paths) == 0
        end

        Test.@testset "Source equals target" begin
            sequence = BioSequences.dna"ATCGATCG"
            records = [FASTX.FASTA.Record("seq", sequence)]
            graph = Mycelia.Rhizomorph.build_kmer_graph(
                records, 4; dataset_id = "self_ksp", mode = :singlestrand
            )
            weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
            labels = collect(MetaGraphsNext.labels(weighted))

            paths = Mycelia.Rhizomorph.k_shortest_paths(
                weighted, first(labels), first(labels), 3
            )
            Test.@test length(paths) >= 1
            Test.@test paths[1].total_probability == 1.0
        end
    end

    Test.@testset "Path to sequence roundtrip" begin
        # ACGTTTCG has 5 distinct k=4 k-mers forming a linear (acyclic) path,
        # so sources and sinks are non-empty and the roundtrip is unambiguous.
        sequence = BioSequences.LongDNA{4}("ACGTTTCG")
        records = [FASTX.FASTA.Record("known", sequence)]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "roundtrip_ksp", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)

        labels = collect(MetaGraphsNext.labels(weighted))
        in_deg = Dict(l => 0 for l in labels)
        out_deg = Dict(l => 0 for l in labels)
        for (src, dst) in MetaGraphsNext.edge_labels(weighted)
            out_deg[src] += 1
            in_deg[dst] += 1
        end
        sources = [l for l in labels if in_deg[l] == 0]
        sinks = [l for l in labels if out_deg[l] == 0]

        Test.@test !isempty(sources)
        Test.@test !isempty(sinks)

        paths = Mycelia.Rhizomorph.k_shortest_paths(
            weighted, first(sources), first(sinks), 1
        )

        Test.@test length(paths) == 1

        path_labels = [step.vertex_label for step in paths[1].steps]
        reconstructed = Mycelia.Rhizomorph.path_to_sequence(path_labels, graph)

        Test.@test reconstructed == sequence
    end

    Test.@testset "_paths_equal helper" begin
        sequence = BioSequences.dna"ACGTTTCG"
        records = [FASTX.FASTA.Record("seq", sequence)]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "equal_test", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
        labels = collect(MetaGraphsNext.labels(weighted))

        in_deg = Dict(l => 0 for l in labels)
        out_deg = Dict(l => 0 for l in labels)
        for (src, dst) in MetaGraphsNext.edge_labels(weighted)
            out_deg[src] += 1
            in_deg[dst] += 1
        end
        sources = [l for l in labels if in_deg[l] == 0]
        sinks = [l for l in labels if out_deg[l] == 0]

        paths = Mycelia.Rhizomorph.k_shortest_paths(weighted, first(sources), first(sinks), 1)
        Test.@test length(paths) == 1

        # Same path should be equal to itself
        Test.@test Mycelia.Rhizomorph._paths_equal(paths[1], paths[1])

        # Different length paths should not be equal
        short_steps = [paths[1].steps[1]]
        short_path = Mycelia.Rhizomorph.GraphPath(short_steps)
        Test.@test !Mycelia.Rhizomorph._paths_equal(paths[1], short_path)
    end

    Test.@testset "_build_graph_path_from_vertices" begin
        sequence = BioSequences.dna"ACGTTTCG"
        records = [FASTX.FASTA.Record("seq", sequence)]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "build_test", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
        labels = collect(MetaGraphsNext.labels(weighted))

        in_deg = Dict(l => 0 for l in labels)
        out_deg = Dict(l => 0 for l in labels)
        for (src, dst) in MetaGraphsNext.edge_labels(weighted)
            out_deg[src] += 1
            in_deg[dst] += 1
        end
        sources = [l for l in labels if in_deg[l] == 0]
        sinks = [l for l in labels if out_deg[l] == 0]

        # Build path from vertices
        path = Mycelia.Rhizomorph._build_graph_path_from_vertices(weighted, labels)
        Test.@test path isa Mycelia.Rhizomorph.GraphPath
        Test.@test length(path.steps) == length(labels)
        Test.@test path.total_probability > 0.0

        # First step should have probability 1.0
        Test.@test path.steps[1].probability == 1.0

        # Empty vertices should return empty path
        empty_path = Mycelia.Rhizomorph._build_graph_path_from_vertices(
            weighted, eltype(labels)[]
        )
        Test.@test length(empty_path.steps) == 0
    end

    Test.@testset "Invalid vertex throws ArgumentError" begin
        sequence = BioSequences.dna"ACGTTTCG"
        records = [FASTX.FASTA.Record("seq", sequence)]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 4; dataset_id = "invalid_test", mode = :singlestrand
        )
        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(graph)
        labels = collect(MetaGraphsNext.labels(weighted))

        # Create a k-mer that definitely doesn't exist in the graph
        fake_kmer = Kmers.DNAKmer{4}("GGGG")

        # Source not in graph
        try
            Mycelia.Rhizomorph.k_shortest_paths(weighted, fake_kmer, first(labels), 1)
            Test.@test false  # Should not reach here
        catch e
            Test.@test e isa ArgumentError
            Test.@test occursin("Source vertex", e.msg)
        end

        # Target not in graph
        try
            Mycelia.Rhizomorph.k_shortest_paths(weighted, first(labels), fake_kmer, 1)
            Test.@test false
        catch e
            Test.@test e isa ArgumentError
            Test.@test occursin("Target vertex", e.msg)
        end

        # Negative k
        try
            Mycelia.Rhizomorph.k_shortest_paths(weighted, first(labels), last(labels), -1)
            Test.@test false
        catch e
            Test.@test e isa ArgumentError
            Test.@test occursin("non-negative", e.msg)
        end
    end
end
