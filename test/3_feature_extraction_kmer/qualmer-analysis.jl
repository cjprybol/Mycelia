# test/test_qualmer_analysis.jl
using Test
import BioSequences
import FASTX
import Mycelia
import Kmers

@testset "Qualmer Analysis" begin
    
    @testset "Qualmer Types" begin
        # Test DNA qualmer creation and basic operations
        @testset "DNAQualmer" begin
            dna_seq = BioSequences.LongDNA{4}("ACTG")
            kmer = Kmers.DNAKmer{4}(dna_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            
            # Test constructor
            qmer = Mycelia.Qualmer(kmer, quals)
            @test qmer.kmer == kmer
            @test qmer.qualities == quals
            
            # Test methods
            @test length(qmer) == 4
            @test qmer[1] == (kmer[1], quals[1])
            @test qmer[4] == (kmer[4], quals[4])
            
            # Test equality
            different_quals = (0x20, 0x21, 0x22, 0x24)
            qmer3 = Mycelia.Qualmer(kmer, different_quals)
            @test qmer != qmer3
            
            # Test error cases
            @test_throws AssertionError Mycelia.Qualmer(kmer, (0x20, 0x21, 0x22)) # Wrong length
        end
        
        # Test RNA qualmer
        @testset "RNAQualmer" begin
            rna_seq = BioSequences.LongRNA{4}("ACUG")
            kmer = Kmers.RNAKmer{4}(rna_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            
            qmer = Mycelia.Qualmer(kmer, quals)
            @test qmer.kmer == kmer
            @test qmer.qualities == quals
        end
        
        # Test AA qualmer
        @testset "AAQualmer" begin
            aa_seq = BioSequences.LongAA("MWKL")
            kmer = Kmers.AAKmer{4}(aa_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            
            qmer = Mycelia.Qualmer(kmer, quals)
            @test qmer.kmer == kmer
            @test qmer.qualities == quals
        end
    end
    
    @testset "Canonical Representation" begin
        @testset "DNA Canonical" begin
            # Forward case (already canonical)
            dna_seq1 = BioSequences.LongDNA{4}("ACTG")
            kmer1 = Kmers.DNAKmer{4}(dna_seq1)
            quals1 = (0x20, 0x21, 0x22, 0x23)
            qmer1 = Mycelia.Qualmer(kmer1, quals1)
            
            canon_qmer1 = Mycelia.canonical(qmer1)
            @test canon_qmer1.kmer == kmer1  # Already canonical
            @test canon_qmer1.qualities == quals1  # Qualities unchanged
            
            # Reverse complement case
            dna_seq2 = BioSequences.LongDNA{4}("CAGT")  # RC of ACTG
            kmer2 = Kmers.DNAKmer{4}(dna_seq2)
            quals2 = (0x30, 0x31, 0x32, 0x33)
            qmer2 = Mycelia.Qualmer(kmer2, quals2)
            
            canon_qmer2 = Mycelia.canonical(qmer2)
            @test canon_qmer2.kmer == BioSequences.canonical(kmer2)
            @test canon_qmer2.qualities == (0x33, 0x32, 0x31, 0x30)  # Reversed qualities
        end
        
        @testset "RNA Canonical" begin
            # Test both regular and reverse-complement cases
            rna_seq = BioSequences.LongRNA{4}("ACGU")
            kmer = Kmers.RNAKmer{4}(rna_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            qmer = Mycelia.Qualmer(kmer, quals)
            
            canon_qmer = Mycelia.canonical(qmer)
            @test canon_qmer.kmer == kmer  # Already canonical
            @test canon_qmer.qualities !== quals
            @test canon_qmer.qualities == reverse(quals)
        end
        
        @testset "AA Canonical" begin
            # Amino acids don't have canonical form (no RC)
            aa_seq = BioSequences.LongAA("MWKL")
            kmer = Kmers.AAKmer{4}(aa_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            qmer = Mycelia.Qualmer(kmer, quals)
            
            canon_qmer = Mycelia.canonical(qmer)
            @test canon_qmer === qmer  # Should be identical
        end
    end
    
    @testset "Qualmer Iteration" begin
        @testset "DNA Qualmers" begin
            # Create test sequence and quality
            dna_seq = BioSequences.LongDNA{4}("ACTGCATGCAATGC")
            quality = UInt8[0x20 + i for i in 1:length(dna_seq)]
            k = 5
            
            # Test qualmers_fw
            fw_qualmers = collect(Mycelia.qualmers_fw(dna_seq, quality, Val(k)))
            @test length(fw_qualmers) == length(dna_seq) - k + 1
            @test fw_qualmers[1][2] == 1  # First position
            @test fw_qualmers[1][1].kmer == Kmers.Kmer{BioSequences.DNAAlphabet{4}, k}(dna_seq[1:k])
            @test fw_qualmers[1][1].qualities == Tuple(quality[1:k])
            
            # Test qualmers_unambiguous with sequence containing ambiguous bases
            ambig_seq = BioSequences.LongDNA{4}("ACTGNATGC")
            ambig_qual = UInt8[0x20 + i for i in 1:length(ambig_seq)]
            
            unambig_qualmers = collect(Mycelia.qualmers_unambiguous(ambig_seq, ambig_qual, Val(3)))
            @test length(unambig_qualmers) < length(ambig_seq) - 3 + 1  # Should skip ambiguous regions
            
            # Test qualmers_canonical
            canonical_qualmers = collect(Mycelia.qualmers_canonical(dna_seq, quality, Val(k)))
            @test length(canonical_qualmers) == length(dna_seq) - k + 1
        end
        
        @testset "RNA Qualmers" begin
            # Similar tests for RNA
            rna_seq = BioSequences.LongRNA{4}("ACUGCAUGCAAUGC")
            quality = UInt8[0x20 + i for i in 1:length(rna_seq)]
            k = 5
            
            fw_qualmers = collect(Mycelia.qualmers_fw(rna_seq, quality, Val(k)))
            @test length(fw_qualmers) == length(rna_seq) - k + 1
        end
        
        @testset "AA Qualmers" begin
            # Tests for amino acid sequences
            aa_seq = BioSequences.LongAA("MWKLVPGKEC")
            quality = UInt8[0x20 + i for i in 1:length(aa_seq)]
            k = 3
            
            fw_qualmers = collect(Mycelia.qualmers_fw(aa_seq, quality, Val(k)))
            @test length(fw_qualmers) == length(aa_seq) - k + 1
            
            # AA sequences don't have ambiguous characters handled specially
            unambig_qualmers = collect(Mycelia.qualmers_unambiguous(aa_seq, quality, Val(k)))
            @test length(unambig_qualmers) == length(aa_seq) - k + 1
        end
    end
    
    @testset "FASTQ Record Integration" begin
        # Create a test FASTQ record
        header = "test_read"
        seq_str = "ACTGCATGCAATGC"
        quality_str = "IIIIHHHHGGGFFF"
        
        record = FASTX.FASTQ.Record(header, seq_str, quality_str)
        
        @testset "Extract Qualmers from FASTQ" begin
            k = 5
            
            fw_qualmers = collect(Mycelia.qualmers_fw(record, k))
            @test length(fw_qualmers) == length(seq_str) - k + 1
            
            canonical_qualmers = collect(Mycelia.qualmers_canonical(record, k))
            @test length(canonical_qualmers) == length(seq_str) - k + 1
            
            # Test unambiguous filtering
            ambig_record = FASTX.FASTQ.Record(header, "ACTGNCTGA", "IIIIIIIII")
            unambig_qualmers = collect(Mycelia.qualmers_unambiguous(ambig_record, 3))
            @test length(unambig_qualmers) < length(FASTX.sequence(ambig_record)) - 3 + 1
            
            # Test combined canonical and unambiguous
            unambig_canonical = collect(Mycelia.qualmers_unambiguous_canonical(record, k))
            @test length(unambig_canonical) == length(seq_str) - k + 1
        end
    end
    
    @testset "Edge Cases" begin
        # Test empty sequence
        empty_dna = BioSequences.LongDNA{4}("")
        empty_qual = UInt8[]
        @test isempty(collect(Mycelia.qualmers_fw(empty_dna, empty_qual, Val(3))))
        
        # Test sequence shorter than k
        short_dna = BioSequences.LongDNA{4}("AC")
        short_qual = UInt8[0x20, 0x21]
        @test isempty(collect(Mycelia.qualmers_fw(short_dna, short_qual, Val(3))))
        
        # Test sequence exactly length k
        exact_dna = BioSequences.LongDNA{4}("ACTG")
        exact_qual = UInt8[0x20, 0x21, 0x22, 0x23]
        exact_qualmers = collect(Mycelia.qualmers_fw(exact_dna, exact_qual, Val(4)))
        @test length(exact_qualmers) == 1
    end

    @testset "Simple Qualmer Examples" begin
        # Basic Qualmer construction for each alphabet
        dna_kmer = Kmers.DNAKmer{3}(BioSequences.LongDNA{3}("ACG"))
        rna_kmer = Kmers.RNAKmer{3}(BioSequences.LongRNA{3}("ACG"))
        aa_kmer  = Kmers.AAKmer{3}(BioSequences.LongAA("ACD"))
        qual = (0x01, 0x02, 0x03)

        dna_q = Mycelia.Qualmer(dna_kmer, qual)
        rna_q = Mycelia.Qualmer(rna_kmer, qual)
        aa_q  = Mycelia.Qualmer(aa_kmer,  qual)

        @test dna_q.kmer == dna_kmer && dna_q.qualities == qual
        @test rna_q.kmer == rna_kmer && rna_q.qualities == qual
        @test aa_q.kmer  == aa_kmer  && aa_q.qualities  == qual

        # Canonical with reverse complement for DNA and RNA
        rc_dna = Kmers.DNAKmer{3}(BioSequences.LongDNA{3}("GAT"))
        rc_q   = (0x10, 0x11, 0x12)
        canon_dna = Mycelia.canonical(Mycelia.Qualmer(rc_dna, rc_q))
        @test canon_dna.kmer == BioSequences.canonical(rc_dna)
        @test canon_dna.qualities == (0x12, 0x11, 0x10)

        rc_rna = Kmers.RNAKmer{3}(BioSequences.LongRNA{3}("GAU"))
        canon_rna = Mycelia.canonical(Mycelia.Qualmer(rc_rna, rc_q))
        @test canon_rna.kmer == BioSequences.canonical(rc_rna)
        @test canon_rna.qualities == (0x12, 0x11, 0x10)

        # Iterator checks with small sequences
        seq = BioSequences.LongDNA{4}("ACGT")
        q   = UInt8[1,2,3,4]

        fw = collect(Mycelia.qualmers_fw(seq, q, Val(2)))
        @test length(fw) == 3
        @test fw[1] == (Mycelia.Qualmer(Kmers.DNAKmer{2}(BioSequences.LongDNA{2}("AC")), (0x01,0x02)), 1)

        unambig = collect(Mycelia.qualmers_unambiguous(seq, q, Val(2)))
        @test length(unambig) == 3
        @test unambig[3] == (Mycelia.Qualmer(Kmers.DNAKmer{2}(BioSequences.LongDNA{2}("GT")), (0x03,0x04)), 3)

        canon = collect(Mycelia.qualmers_canonical(seq, q, Val(2)))
        @test length(canon) == 3
        @test canon[3] == (Mycelia.Qualmer(BioSequences.canonical(Kmers.DNAKmer{2}(BioSequences.LongDNA{2}("GT"))), (0x04,0x03)), 3)
    end
end
