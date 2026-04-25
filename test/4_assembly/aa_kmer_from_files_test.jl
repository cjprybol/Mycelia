import Test
import Mycelia
import FASTX
import Kmers

function write_aa_fasta(path::AbstractString, entries::Vector{Tuple{String, String}})
    open(path, "w") do io
        for (id, sequence) in entries
            println(io, ">$id")
            println(io, sequence)
        end
    end
    return path
end

Test.@testset "AA K-mer Graph From Files" begin
    Test.@testset "multi-file singlestrand merges AA observations" begin
        mktempdir() do dir
            path_a = write_aa_fasta(joinpath(dir, "sample_a.faa"), [("prot_a", "MKVLW")])
            path_b = write_aa_fasta(joinpath(dir, "sample_b.faa"), [("prot_b", "MKVLY")])

            graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [path_a, path_b],
                3;
                mode = :singlestrand
            )

            kmer = Kmers.AAKmer{3}("MKV")
            Test.@test Mycelia.Rhizomorph.has_vertex(graph, kmer)

            vertex_data = graph[kmer]
            Test.@test Set(keys(vertex_data.evidence)) == Set(["sample_a", "sample_b"])
            Test.@test haskey(vertex_data.evidence["sample_a"], "prot_a")
            Test.@test haskey(vertex_data.evidence["sample_b"], "prot_b")
        end
    end

    Test.@testset "incremental add works for amino acid graphs" begin
        graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand(
            [FASTX.FASTA.Record("prot_a", "MKVLW")],
            3;
            dataset_id = "sample_a"
        )

        Mycelia.Rhizomorph.add_observations_to_graph!(
            graph,
            [FASTX.FASTA.Record("prot_b", "MKVLY")],
            3;
            dataset_id = "sample_b"
        )

        kmer = Kmers.AAKmer{3}("MKV")
        Test.@test Mycelia.Rhizomorph.get_vertex_observation_count(graph, kmer) == 2
        Test.@test Set(keys(graph[kmer].evidence)) == Set(["sample_a", "sample_b"])
    end

    Test.@testset "AA multi-file mode guards trigger after merge" begin
        mktempdir() do dir
            path_a = write_aa_fasta(joinpath(dir, "sample_a.faa"), [("prot_a", "MKVLW")])
            path_b = write_aa_fasta(joinpath(dir, "sample_b.faa"), [("prot_b", "MKVLY")])

            Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [path_a, path_b],
                3;
                mode = :doublestrand
            )
            Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [path_a, path_b],
                3;
                mode = :canonical
            )
        end
    end
end
