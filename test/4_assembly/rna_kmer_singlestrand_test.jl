# RNA K-mer Graph - SingleStrand Mode Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/rna_kmer_singlestrand_test.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import MetaGraphsNext

Test.@testset "RNA K-mer SingleStrand Graph (Rhizomorph)" begin
    dataset_id = "rna_single"
    test_rna = BioSequences.rna"AUCGAUCG"
    reads = [FASTX.FASTA.Record("test", test_rna)]
    observation_id = FASTX.identifier(first(reads))

    graph = Mycelia.Rhizomorph.build_kmer_graph(
        reads,
        3;
        dataset_id=dataset_id,
        mode=:singlestrand,
    )

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test Set(vertices) == Set([
        BioSequences.RNAKmer{3}("AUC"),
        BioSequences.RNAKmer{3}("UCG"),
        BioSequences.RNAKmer{3}("CGA"),
        BioSequences.RNAKmer{3}("GAU"),
    ])
    Test.@test MetaGraphsNext.ne(graph) == 4

    for (vertex_label, expected_positions) in Dict(
        BioSequences.RNAKmer{3}("AUC") => [1, 5],
        BioSequences.RNAKmer{3}("UCG") => [2, 6],
        BioSequences.RNAKmer{3}("CGA") => [3],
        BioSequences.RNAKmer{3}("GAU") => [4],
    )
        vertex_data = graph[vertex_label]
        Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
        Test.@test Mycelia.Rhizomorph.count_evidence(vertex_data) == length(expected_positions)

        dataset_evidence = Mycelia.Rhizomorph.get_dataset_evidence(vertex_data, dataset_id)
        Test.@test dataset_evidence !== nothing
        observation_evidence = Mycelia.Rhizomorph.get_observation_evidence(vertex_data, dataset_id, observation_id)
        Test.@test observation_evidence !== nothing

        positions = sort(entry.position for entry in observation_evidence)
        Test.@test positions == expected_positions
        Test.@test all(entry.strand == Mycelia.Rhizomorph.Forward for entry in observation_evidence)
    end

    expected_edge_positions = Dict(
        (BioSequences.RNAKmer{3}("AUC"), BioSequences.RNAKmer{3}("UCG")) => [(1, 2), (5, 6)],
        (BioSequences.RNAKmer{3}("UCG"), BioSequences.RNAKmer{3}("CGA")) => [(2, 3)],
        (BioSequences.RNAKmer{3}("CGA"), BioSequences.RNAKmer{3}("GAU")) => [(3, 4)],
        (BioSequences.RNAKmer{3}("GAU"), BioSequences.RNAKmer{3}("AUC")) => [(4, 5)],
    )

    for edge_label in MetaGraphsNext.edge_labels(graph)
        src, dst = edge_label
        edge_data = graph[src, dst]
        observation_evidence = Mycelia.Rhizomorph.get_observation_evidence(edge_data, dataset_id, observation_id)
        Test.@test observation_evidence !== nothing
        pairs = sort((entry.from_position, entry.to_position) for entry in observation_evidence)
        Test.@test pairs == sort(expected_edge_positions[(src, dst)])
        Test.@test all(entry.strand == Mycelia.Rhizomorph.Forward for entry in observation_evidence)
    end

    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
    Test.@test !isempty(paths)

    walk_steps = [
        Mycelia.Rhizomorph.WalkStep(label, Mycelia.Rhizomorph.Forward, 1.0, Float64(i))
        for (i, label) in enumerate(first(paths))
    ]
    graph_path = Mycelia.Rhizomorph.GraphPath(walk_steps)

    reconstructed = Mycelia.Rhizomorph.path_to_sequence(graph_path, graph)
    Test.@test reconstructed == "AUCGAUCG"
    Test.@test !occursin('T', reconstructed)
end
