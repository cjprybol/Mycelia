import Test
import Mycelia
import MetaGraphsNext
import Graphs
import BioSequences

Test.@testset "Graph cleanup" begin
    graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type=Int,
        vertex_data_type=Dict{Symbol, Any},
        edge_data_type=Nothing,
    )

    graph[1] = Dict(:coverage => 1.0, :sequence => BioSequences.LongDNA{4}("AAAA"))
    graph[2] = Dict(:coverage => 10.0, :sequence => BioSequences.LongDNA{4}("ATGC"))
    graph[3] = Dict(:coverage => 1.0, :sequence => BioSequences.LongDNA{4}("ATGC"))

    Graphs.add_edge!(graph, 1, 2)
    Graphs.add_edge!(graph, 2, 3)

    cleaned_graph, stats = Mycelia.statistical_tip_clipping(
        graph;
        min_coverage_threshold=1,
        std_dev_multiplier=3.0,
        preserve_high_quality=true,
    )

    Test.@test !MetaGraphsNext.haskey(cleaned_graph, 1)
    Test.@test MetaGraphsNext.haskey(cleaned_graph, 3)
    Test.@test stats[:tips_removed] == 1
    Test.@test stats[:high_quality_preserved] == 1

    component_graph = MetaGraphsNext.MetaGraph(
        Graphs.DiGraph();
        label_type=Int,
        vertex_data_type=Dict{Symbol, Any},
        edge_data_type=Nothing,
    )

    component_graph[1] = Dict(:coverage => 5.0)
    component_graph[2] = Dict(:coverage => 5.0)
    component_graph[3] = Dict(:coverage => 5.0)

    Graphs.add_edge!(component_graph, 1, 2)

    components = Mycelia.find_connected_components(component_graph)

    Test.@test length(components) == 2
    Test.@test any(component -> Set(component) == Set([1, 2]), components)
    Test.@test any(component -> Set(component) == Set([3]), components)
end
