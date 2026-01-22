import Test
import Mycelia
import FASTX
import MetaGraphsNext

Test.@testset "FASTA Graph Utilities" begin
    mktempdir() do dir
        fasta_path = joinpath(dir, "reads.fna")
        open(fasta_path, "w") do io
            write(io, ">read1\n")
            write(io, "ATCG\n")
            write(io, ">read2\n")
            write(io, "TCGA\n")
        end

        graph = Mycelia.Rhizomorph.build_fasta_graph_from_file(fasta_path; min_overlap=3)
        stats = Mycelia.Rhizomorph.get_fasta_graph_statistics(graph)
        Test.@test stats[:num_vertices] == 2
        Test.@test stats[:num_edges] == 1
        Test.@test stats[:sequence_type] == "DNA"

        fasta_path2 = joinpath(dir, "reads2.fna")
        open(fasta_path2, "w") do io
            write(io, ">read3\n")
            write(io, "TCGA\n")
            write(io, ">read4\n")
            write(io, "CGAT\n")
        end

        graph_multi = Mycelia.Rhizomorph.build_fasta_graph_from_files([fasta_path, fasta_path2]; min_overlap=3)
        Test.@test length(MetaGraphsNext.labels(graph_multi)) >= 2

        fasta_path3 = joinpath(dir, "reads3.faa")
        open(fasta_path3, "w") do io
            write(io, ">prot1\n")
            write(io, "ACDE\n")
        end
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_fasta_graph_from_files([fasta_path, fasta_path3])
        Test.@test_throws ArgumentError Mycelia.Rhizomorph.build_fasta_graph_from_files([fasta_path]; type_hint=:RNA)

        empty_fasta = joinpath(dir, "empty.fasta")
        open(empty_fasta, "w") do _ end
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_fasta_graph_from_file(empty_fasta)

        fastq_path = joinpath(dir, "reads.fastq")
        open(fastq_path, "w") do io
            write(io, "@read1\n")
            write(io, "ATCG\n")
            write(io, "+\n")
            write(io, "IIII\n")
        end
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_fasta_graph_from_file(fastq_path)
    end
end
