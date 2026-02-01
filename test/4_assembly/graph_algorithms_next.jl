# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/graph_algorithms_next.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/graph_algorithms_next.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import MetaGraphsNext
import Graphs
import FASTX
import BioSequences
import Kmers

function bubble_graph()
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = Mycelia.Rhizomorph.KmerVertexData,
        edge_data_type = Mycelia.Rhizomorph.KmerEdgeData
    )

    vertices = ["AAA", "AAT", "ATC", "TCG", "CGA", "GAT", "ATG", "TGC", "TIP"]
    for v in vertices
        graph[v] = Mycelia.Rhizomorph.KmerVertexData(v)
        Mycelia.Rhizomorph.add_evidence!(graph[v], "ds", "obs",
            Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
    end

    edges = [
        ("AAA", "AAT"),
        ("AAT", "ATC"),
        ("ATC", "TCG"),
        ("ATC", "CGA"),
        ("TCG", "GAT"),
        ("CGA", "GAT"),
        ("GAT", "ATG"),
        ("ATG", "TGC")
    ]
    for (src, dst) in edges
        graph[src, dst] = Mycelia.Rhizomorph.KmerEdgeData()
        Mycelia.Rhizomorph.add_evidence!(graph[src, dst], "ds", "obs",
            Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Forward))
    end

    return graph
end

function cycle_graph()
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type = String,
        vertex_data_type = Mycelia.Rhizomorph.KmerVertexData,
        edge_data_type = Mycelia.Rhizomorph.KmerEdgeData
    )

    for v in ["ATC", "TCG", "CGA"]
        graph[v] = Mycelia.Rhizomorph.KmerVertexData(v)
        Mycelia.Rhizomorph.add_evidence!(graph[v], "cycle", "obs",
            Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
    end

    for (src, dst) in [("ATC", "TCG"), ("TCG", "CGA"), ("CGA", "ATC")]
        graph[src, dst] = Mycelia.Rhizomorph.KmerEdgeData()
        Mycelia.Rhizomorph.add_evidence!(graph[src, dst], "cycle", "obs",
            Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Forward))
    end

    return graph
end

Test.@testset "Rhizomorph Graph Algorithms" begin
    Test.@testset "Eulerian paths" begin
        g = cycle_graph()
        paths = Mycelia.Rhizomorph.find_eulerian_paths_next(g)
        Test.@test length(paths) == 1
        Test.@test length(first(paths)) == 4  # cycle returns to start

        start_idx = MetaGraphsNext.code_for(g, "ATC")
        dfs_indices = Mycelia.Rhizomorph._find_eulerian_path_dfs(g.graph, start_idx)
        Test.@test !isempty(dfs_indices)
        dfs_labels = [MetaGraphsNext.label_for(g, idx) for idx in dfs_indices]
        Test.@test dfs_labels == first(paths)

        manual_path = ["ATC", "TCG", "CGA"]
        seq = Mycelia.Rhizomorph.path_to_sequence(manual_path, g)
        Test.@test seq == "ATCGA"

        empty_graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph();
            label_type = String,
            vertex_data_type = Mycelia.Rhizomorph.KmerVertexData,
            edge_data_type = Mycelia.Rhizomorph.KmerEdgeData
        )
        Test.@test isempty(Mycelia.Rhizomorph.find_eulerian_paths_next(empty_graph))
    end

    Test.@testset "Bubble detection and support" begin
        g = bubble_graph()
        bubbles = Mycelia.Rhizomorph.detect_bubbles_next(g, min_bubble_length = 1, max_bubble_length = 5)
        Test.@test !isempty(bubbles)
        Test.@test any(b.entry_vertex == "ATC" && b.exit_vertex == "GAT" for b in bubbles)

        support = Mycelia.Rhizomorph.calculate_path_support(g, ["AAA", "AAT", "ATC"])
        Test.@test support == 3  # one evidence entry per vertex

        none = Mycelia.Rhizomorph.detect_bubbles_next(g, min_bubble_length = 10, max_bubble_length = 2)
        Test.@test isempty(none)

        empty_graph = MetaGraphsNext.MetaGraph(
            Graphs.DiGraph();
            label_type = String,
            vertex_data_type = Mycelia.Rhizomorph.KmerVertexData,
            edge_data_type = Mycelia.Rhizomorph.KmerEdgeData
        )
        Test.@test isempty(Mycelia.Rhizomorph.detect_bubbles_next(empty_graph))
        Test.@test Mycelia.Rhizomorph.calculate_path_support(empty_graph, String[]) == 0

        boosted = bubble_graph()
        Mycelia.Rhizomorph.add_evidence!(boosted["TCG"], "ds", "extra_obs",
            Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(boosted["GAT"], "ds", "extra_obs",
            Mycelia.Rhizomorph.EvidenceEntry(11, Mycelia.Rhizomorph.Forward))
        high_support = Mycelia.Rhizomorph.calculate_path_support(boosted, [
            "ATC", "TCG", "GAT"])
        low_support = Mycelia.Rhizomorph.calculate_path_support(boosted, [
            "ATC", "CGA", "GAT"])
        Test.@test high_support > low_support
    end

    Test.@testset "Bubble helper utilities" begin
        g = bubble_graph()
        path1 = Mycelia.Rhizomorph.find_limited_path(g, "TCG", 5)
        path2 = Mycelia.Rhizomorph.find_limited_path(g, "CGA", 5)
        Test.@test last(path1) == "TGC"
        Test.@test last(path2) == "TGC"

        convergence = Mycelia.Rhizomorph.find_path_convergence(path1, path2)
        Test.@test convergence == "GAT"

        complexity = Mycelia.Rhizomorph.calculate_bubble_complexity(path1, path2)
        Test.@test 0.0 <= complexity <= 1.0

        bubble = Mycelia.Rhizomorph.BubbleStructure(
            "ATC", "GAT", path1, path2, 4, 3, complexity)
        deduped = Mycelia.Rhizomorph.remove_duplicate_bubbles([bubble, bubble])
        Test.@test length(deduped) == 1
    end

    Test.@testset "Neighbor helpers" begin
        g = bubble_graph()
        outs = Mycelia.Rhizomorph.get_out_neighbors(g, "ATC")
        ins = Mycelia.Rhizomorph.get_in_neighbors(g, "GAT")
        Test.@test length(outs) == 2
        Test.@test length(ins) == 2
        Test.@test Set(outs) == Set(["TCG", "CGA"])
        Test.@test Set(ins) == Set(["TCG", "CGA"])
        Test.@test isempty(Mycelia.Rhizomorph.get_out_neighbors(g, "TGC"))
        Test.@test isempty(Mycelia.Rhizomorph.get_out_neighbors(g, "missing"))
    end

    Test.@testset "Tip removal" begin
        g = bubble_graph()
        Test.@test haskey(g, "TIP")
        Mycelia.Rhizomorph.remove_tips!(g; min_support = 1)
        Test.@test !haskey(g, "TIP")

        supported = bubble_graph()
        Mycelia.Rhizomorph.add_evidence!(supported["TIP"], "ds", "keep_me",
            Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.remove_tips!(supported; min_support = 1)
        Test.@test haskey(supported, "TIP")
    end

    Test.@testset "Path validation and sequences" begin
        g = cycle_graph()
        valid_path = ["ATC", "TCG"]
        Test.@test Mycelia.Rhizomorph.is_valid_path(g, valid_path)
        Test.@test !Mycelia.Rhizomorph.is_valid_path(g, ["ATC", "AAA"])

        # GraphPath-based sequence reconstruction
        steps = [
            Mycelia.Rhizomorph.WalkStep("ATC", Mycelia.Rhizomorph.Forward, 1.0, 1.0),
            Mycelia.Rhizomorph.WalkStep("TCG", Mycelia.Rhizomorph.Forward, 1.0, 1.0)
        ]
        gp = Mycelia.Rhizomorph.GraphPath(steps)
        seq = Mycelia.Rhizomorph.path_to_sequence(gp, g)
        Test.@test seq == "ATCG"
    end

    Test.@testset "GraphPath reverse strand reconstruction" begin
        reads = [FASTX.FASTA.Record("dna", BioSequences.dna"ATCG")]
        kgraph = Mycelia.Rhizomorph.build_kmer_graph(reads, 3; dataset_id = "rev_path", mode = :singlestrand)
        path = [Kmers.DNAKmer{3}("ATC"), Kmers.DNAKmer{3}("TCG")]

        reverse_steps = [
            Mycelia.Rhizomorph.WalkStep(path[1], Mycelia.Rhizomorph.Reverse, 0.5, 0.5),
            Mycelia.Rhizomorph.WalkStep(path[2], Mycelia.Rhizomorph.Reverse, 1.0, 1.0)
        ]
        reverse_path = Mycelia.Rhizomorph.GraphPath(reverse_steps)
        reverse_seq = Mycelia.Rhizomorph.path_to_sequence(reverse_path, kgraph)

        Test.@test reverse_seq isa BioSequences.LongDNA
        Test.@test startswith(string(reverse_seq), "GAT")
    end

    Test.@testset "Simplification smoke test" begin
        g = bubble_graph()
        bubbles = Mycelia.Rhizomorph.detect_bubbles_next(g)
        simplified = Mycelia.Rhizomorph.simplify_graph_next(g, bubbles)
        Test.@test simplified isa typeof(g)
        Test.@test length(MetaGraphsNext.labels(simplified)) <=
                   length(MetaGraphsNext.labels(g))
    end
end
