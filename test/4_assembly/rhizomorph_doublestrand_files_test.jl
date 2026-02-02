# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_doublestrand_files_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rhizomorph_doublestrand_files_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Rhizomorph doublestrand conversion via file-based builders

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "Rhizomorph doublestrand from files" begin
    # K-mer graph from FASTA files with mode=:doublestrand
    mktempdir() do dir
        fasta = joinpath(dir, "reads.fasta")
        open(fasta, "w") do io
            println(io, ">seq1")
            println(io, "ATGCAT")
        end
        graph = Mycelia.Rhizomorph.build_kmer_graph_from_files([fasta], 3; mode = :doublestrand)

        # Expect both forward and reverse-complement k-mers
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("ATG"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("CAT"))
    end

    # Qualmer graph from FASTQ files with mode=:doublestrand
    mktempdir() do dir
        fastq = joinpath(dir, "reads.fastq")
        record = FASTX.FASTQ.Record("read1", "ATGCAT", "IIIIII")
        open(FASTX.FASTQ.Writer, fastq) do writer
            write(writer, record)
        end
        graph = Mycelia.Rhizomorph.build_qualmer_graph_from_files([fastq], 3; mode = :doublestrand)

        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("ATG"))
        Test.@test Mycelia.Rhizomorph.has_vertex(graph, Kmers.DNAKmer{3}("CAT"))
    end
end

println("✓ Rhizomorph doublestrand from files tests completed")
