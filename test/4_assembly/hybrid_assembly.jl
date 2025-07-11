# Hybrid Assembly tests
import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
Test.@testset "Hybrid Assembly" begin
    Test.@testset "Assembly Core" begin
        kmers = [BioSequences.LongDNA{2}("AAA"), BioSequences.LongDNA{2}("AAT"), BioSequences.LongDNA{2}("ATG")]
        index = Mycelia.get_kmer_index(kmers, kmers[2])
        Test.@test index == 2
        seq = Mycelia.kmer_path_to_sequence(kmers)
        Test.@test String(seq) == "AAATG"
    end
    Test.@testset "Contig Overlap Graph Integrity" begin
        obs = Mycelia.ngrams("ACGT", 2)
        Test.@test obs == ["AC", "CG", "GT"]
        g = Mycelia.string_to_ngram_graph(s="ACGT", n=2)
        Test.@test Graphs.nv(g) == 3
        Test.@test Graphs.ne(g) == 2
    end
end
