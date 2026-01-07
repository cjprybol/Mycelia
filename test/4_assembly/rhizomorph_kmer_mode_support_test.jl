# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_kmer_mode_support_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rhizomorph_kmer_mode_support_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import Graphs
import MetaGraphsNext

Test.@testset "Rhizomorph k-mer graph modes" begin
    mktempdir() do dir
        file1 = joinpath(dir, "forward.fasta")
        file2 = joinpath(dir, "reverse.fasta")

        # Forward and reverse-complement content to exercise canonical merge
        rec1 = FASTX.FASTA.Record("fwd", BioSequences.LongDNA{4}("ATGCA"))
        rec2 = FASTX.FASTA.Record("rev", BioSequences.LongDNA{4}("TGCAT"))
        Mycelia.write_fasta(outfile=file1, records=[rec1])
        Mycelia.write_fasta(outfile=file2, records=[rec2])

        graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
            [file1, file2], 3; mode=:canonical
        )

        # Canonical graph should collapse RC pairs: ATG↔CAT, TGC↔GCA
        Test.@test MetaGraphsNext.nv(graph) == 2
        Test.@test Set(MetaGraphsNext.labels(graph)) ==
            Set([Kmers.DNAKmer{3}("ATG"), Kmers.DNAKmer{3}("GCA")])

        # Evidence from both files should be present on canonical vertices
        atg_vertex = graph[Kmers.DNAKmer{3}("ATG")]
        Test.@test Set(keys(atg_vertex.evidence)) == Set(["forward", "reverse"])

        # Canonical graph should be undirected
        Test.@test !Graphs.is_directed(graph.graph)
        Test.@test Mycelia.Rhizomorph.has_edge(graph, Kmers.DNAKmer{3}("ATG"), Kmers.DNAKmer{3}("GCA"))
        Test.@test Mycelia.Rhizomorph.has_edge(graph, Kmers.DNAKmer{3}("GCA"), Kmers.DNAKmer{3}("ATG"))  # symmetric in undirected graph
    end

    mktempdir() do dir
        file1 = joinpath(dir, "forward.fasta")
        file2 = joinpath(dir, "reverse.fasta")

        rec1 = FASTX.FASTA.Record("fwd", BioSequences.LongDNA{4}("ATGCA"))
        rec2 = FASTX.FASTA.Record("rev", BioSequences.LongDNA{4}("TGCAT"))
        Mycelia.write_fasta(outfile=file1, records=[rec1])
        Mycelia.write_fasta(outfile=file2, records=[rec2])

        graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
            [file1, file2], 3; mode=:doublestrand
        )

        # Doublestrand graph should carry both orientations and be directed
        Test.@test Graphs.is_directed(graph.graph)

        atg = Kmers.DNAKmer{3}("ATG")
        tgc = Kmers.DNAKmer{3}("TGC")
        rc_atg = BioSequences.reverse_complement(atg)
        rc_tgc = BioSequences.reverse_complement(tgc)

        Test.@test haskey(graph, atg)
        Test.@test haskey(graph, rc_atg)
        # Forward edge present
        Test.@test Mycelia.Rhizomorph.has_edge(graph, atg, tgc)
        # Reverse-complement edge should be reversed in direction
        Test.@test Mycelia.Rhizomorph.has_edge(graph, rc_tgc, rc_atg)
    end

    Test.@testset "Singlestrand graphs are directed" begin
        rec = FASTX.FASTA.Record("fwd", BioSequences.LongDNA{4}("ATGCA"))
        graph = Mycelia.Rhizomorph.build_kmer_graph([rec], 3; mode=:singlestrand)
        Test.@test Graphs.is_directed(graph.graph)
    end

    Test.@testset "Canonical mode not allowed for amino acids" begin
        aa_rec = FASTX.FASTA.Record("aa1", BioSequences.LongAA("ACDEFG"))
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph([aa_rec], 3; mode=:canonical)
    end
end
