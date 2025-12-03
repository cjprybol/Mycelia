# Rhizomorph doublestrand conversion via file-based builders

import Test
import Mycelia
import BioSequences
import FASTX

Test.@testset "Rhizomorph doublestrand from files" begin
    # K-mer graph from FASTA files with mode=:doublestrand
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        open(fasta, "w") do io
            println(io, ">seq1")
            println(io, "ATGCAT")
        end
        graph = Mycelia.Rhizomorph.build_kmer_graph_from_files([fasta], 3; mode=:doublestrand)

        # Expect both forward and reverse-complement k-mers
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, BioSequences.DNAKmer{3}("ATG"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, BioSequences.DNAKmer{3}("CAT"))
    end

    # Qualmer graph from FASTQ files with mode=:doublestrand
    mktempdir() do dir
        fastq = joinpath(dir, "reads.fastq")
        record = FASTX.FASTQ.Record("read1", "ATGCAT", "IIIIII")
        open(FASTX.FASTQ.Writer, fastq) do writer
            write(writer, record)
        end
        graph = Mycelia.Rhizomorph.build_qualmer_graph_from_files([fastq], 3; mode=:doublestrand)

        Test.@test Mycelia.Rhizomorph.has_vertex(graph, BioSequences.DNAKmer{3}("ATG"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, BioSequences.DNAKmer{3}("CAT"))
    end
end

println("✓ Rhizomorph doublestrand from files tests completed")
