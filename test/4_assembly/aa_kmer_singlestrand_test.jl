# Amino Acid K-mer Graph - SingleStrand Mode Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/aa_kmer_singlestrand_test.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "AA K-mer SingleStrand Graph (Rhizomorph)" begin
    dataset_id = "aa_single"
    test_aa = BioSequences.aa"MKTV"
    reads = [FASTX.FASTA.Record("test", test_aa)]

    graph = Mycelia.Rhizomorph.build_kmer_graph(
        reads,
        3;
        dataset_id=dataset_id,
        mode=:singlestrand,
    )

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test Set(vertices) == Set([Kmers.AAKmer{3}("MKT"), Kmers.AAKmer{3}("KTV")])
    Test.@test MetaGraphsNext.ne(graph) == 1

    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
        Test.@test Mycelia.Rhizomorph.count_evidence(vertex_data) == 1
        Test.@test haskey(vertex_data.evidence, dataset_id)
    end

    edge_label = only(collect(MetaGraphsNext.edge_labels(graph)))
    edge_data = graph[edge_label...]
    Test.@test edge_data isa Mycelia.Rhizomorph.KmerEdgeData
    Test.@test Mycelia.Rhizomorph.compute_edge_weight(edge_data) == 1

    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
    Test.@test !isempty(paths)

    walk_steps = [
        Mycelia.Rhizomorph.WalkStep(vertex, Mycelia.Rhizomorph.Forward, 1.0, Float64(i))
        for (i, vertex) in enumerate(first(paths))
    ]
    graph_path = Mycelia.Rhizomorph.GraphPath(walk_steps)
    reconstructed = Mycelia.Rhizomorph.path_to_sequence(graph_path, graph)
    Test.@test string(reconstructed) == "MKTV"
end
