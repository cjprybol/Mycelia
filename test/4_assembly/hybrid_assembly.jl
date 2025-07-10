# Hybrid Assembly tests
@testset "Hybrid Assembly" begin
    @testset "Assembly Core" begin
        kmers = [BioSequences.LongDNA{2}("AAA"), BioSequences.LongDNA{2}("AAT"), BioSequences.LongDNA{2}("ATG")]
        index = Mycelia.get_kmer_index(kmers, kmers[2])
        @test index == 2
        seq = Mycelia.kmer_path_to_sequence(kmers)
        @test String(seq) == "AAATG"
    end
    @testset "Contig Overlap Graph Integrity" begin
        obs = Mycelia.ngrams("ACGT", 2)
        @test obs == ["AC", "CG", "GT"]
        g = Mycelia.string_to_ngram_graph(s="ACGT", n=2)
        @test Graphs.nv(g) == 3
        @test Graphs.ne(g) == 2
    end
end
