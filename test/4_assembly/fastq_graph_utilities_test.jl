import Test
import Mycelia
import FASTX
import MetaGraphsNext

Test.@testset "FASTQ Graph Utilities" begin
    mktempdir() do dir
        fastq_path = joinpath(dir, "reads.fastq")
        open(fastq_path, "w") do io
            write(io, "@read1\n")
            write(io, "ATCG\n")
            write(io, "+\n")
            write(io, "IIII\n")
        end

        graph = Mycelia.Rhizomorph.build_fastq_graph_from_file(fastq_path; min_overlap = 3)
        stats = Mycelia.Rhizomorph.get_fastq_graph_statistics(graph)
        Test.@test stats[:num_vertices] == 1
        Test.@test stats[:num_edges] == 0
        Test.@test stats[:sequence_type] == "DNA"
        Test.@test stats[:mean_quality] !== nothing

        records = Mycelia.Rhizomorph.fastq_graph_to_records(graph, "graph")
        Test.@test length(records) == 1
        Test.@test FASTX.FASTQ.identifier(records[1]) == "graph_1"

        fastq_path2 = joinpath(dir, "reads2.fastq")
        open(fastq_path2, "w") do io
            write(io, "@read2\n")
            write(io, "ATGG\n")
            write(io, "+\n")
            write(io, "IIII\n")
        end

        graph_multi = Mycelia.Rhizomorph.build_fastq_graph_from_files(
            [fastq_path, fastq_path2]; min_overlap = 3)
        Test.@test length(MetaGraphsNext.labels(graph_multi)) >= 1

        fasta_path = joinpath(dir, "reads.fasta")
        open(fasta_path, "w") do io
            write(io, ">read1\n")
            write(io, "ATCG\n")
        end
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_fastq_graph_from_file(fasta_path)
    end
end
