# DNA K-mer Graph - SingleStrand Mode Test
#
# Run with: julia --project=. -e 'include("test/4_assembly/dna_kmer_singlestrand_test.jl")'

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "DNA K-mer SingleStrand Graph (Rhizomorph)" begin
    dataset_id = "dna_single"
    test_dna = BioSequences.dna"ATCGATCG"
    reads = [FASTX.FASTA.Record("test", test_dna)]

    graph = Mycelia.Rhizomorph.build_kmer_graph(
        reads,
        3;
        dataset_id=dataset_id,
        mode=:singlestrand,
    )

    vertices = collect(MetaGraphsNext.labels(graph))
    Test.@test Set(vertices) == Set([
        Kmers.DNAKmer{3}("ATC"),
        Kmers.DNAKmer{3}("TCG"),
        Kmers.DNAKmer{3}("CGA"),
        Kmers.DNAKmer{3}("GAT"),
    ])
    Test.@test MetaGraphsNext.ne(graph) == 4

    for vertex_label in vertices
        vertex_data = graph[vertex_label]
        Test.@test vertex_data isa Mycelia.Rhizomorph.KmerVertexData
        Test.@test Mycelia.Rhizomorph.count_evidence(vertex_data) >= 1
        Test.@test haskey(vertex_data.evidence, dataset_id)
    end

    paths = Mycelia.Rhizomorph.find_eulerian_paths_next(graph)
    Test.@test !isempty(paths)

    reconstructed = Mycelia.Rhizomorph.path_to_sequence(first(paths), graph)
    Test.@test string(reconstructed) == "ATCGATC"
end
