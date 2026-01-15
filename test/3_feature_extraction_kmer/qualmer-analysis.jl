# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/3_feature_extraction_kmer/qualmer-analysis.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/3_feature_extraction_kmer/qualmer-analysis.jl", "test/3_feature_extraction_kmer", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers
import MetaGraphsNext

Test.@testset "Qualmer Analysis" begin
    
    Test.@testset "Qualmer Types" begin
        # Test DNA qualmer creation and basic operations
        Test.@testset "DNAQualmer" begin
            dna_seq = BioSequences.LongDNA{4}("ACTG")
            kmer = Kmers.DNAKmer{4}(dna_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            
            # Test constructor
            qmer = Mycelia.Qualmer(kmer, quals)
            Test.@test qmer.kmer == kmer
            Test.@test qmer.qualities == quals
            
            # Test methods
            Test.@test length(qmer) == 4
            Test.@test qmer[1] == (kmer[1], quals[1])
            Test.@test qmer[4] == (kmer[4], quals[4])
            
            # Test equality
            different_quals = (0x20, 0x21, 0x22, 0x24)
            qmer3 = Mycelia.Qualmer(kmer, different_quals)
            Test.@test qmer != qmer3
            
            # Test error cases
            Test.@test_throws AssertionError Mycelia.Qualmer(kmer, (0x20, 0x21, 0x22)) # Wrong length
        end
        
        # Test RNA qualmer
        Test.@testset "RNAQualmer" begin
            rna_seq = BioSequences.LongRNA{4}("ACUG")
            kmer = Kmers.RNAKmer{4}(rna_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            
            qmer = Mycelia.Qualmer(kmer, quals)
            Test.@test qmer.kmer == kmer
            Test.@test qmer.qualities == quals
        end
        
        # Test AA qualmer
        Test.@testset "AAQualmer" begin
            aa_seq = BioSequences.LongAA("MWKL")
            kmer = Kmers.AAKmer{4}(aa_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            
            qmer = Mycelia.Qualmer(kmer, quals)
            Test.@test qmer.kmer == kmer
            Test.@test qmer.qualities == quals
        end
    end
    
    Test.@testset "Canonical Representation" begin
        Test.@testset "DNA Canonical" begin
            # Forward case (already canonical)
            dna_seq1 = BioSequences.LongDNA{4}("ACTG")
            kmer1 = Kmers.DNAKmer{4}(dna_seq1)
            quals1 = (0x20, 0x21, 0x22, 0x23)
            qmer1 = Mycelia.Qualmer(kmer1, quals1)
            
            canon_qmer1 = Mycelia.canonical(qmer1)
            Test.@test canon_qmer1.kmer == kmer1  # Already canonical
            Test.@test canon_qmer1.qualities == quals1  # Qualities unchanged
            
            # Reverse complement case
            dna_seq2 = BioSequences.LongDNA{4}("CAGT")  # RC of ACTG
            kmer2 = Kmers.DNAKmer{4}(dna_seq2)
            quals2 = (0x30, 0x31, 0x32, 0x33)
            qmer2 = Mycelia.Qualmer(kmer2, quals2)
            
            canon_qmer2 = Mycelia.canonical(qmer2)
            Test.@test canon_qmer2.kmer == BioSequences.canonical(kmer2)
            Test.@test canon_qmer2.qualities == (0x33, 0x32, 0x31, 0x30)  # Reversed qualities
        end
        
        Test.@testset "RNA Canonical" begin
            # Test both regular and reverse-complement cases
            rna_seq = BioSequences.LongRNA{4}("ACGU")
            kmer = Kmers.RNAKmer{4}(rna_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            qmer = Mycelia.Qualmer(kmer, quals)
            
            canon_qmer = Mycelia.canonical(qmer)
            Test.@test canon_qmer.kmer == kmer  # Already canonical
            Test.@test canon_qmer.qualities == quals  # Qualities should be unchanged for canonical k-mer
            
            # Test with a sequence that needs reverse complement
            rna_seq_rc = BioSequences.LongRNA{4}("UGCC")  # RC would be GGCA, so canonical should be GGCA
            kmer_rc = Kmers.RNAKmer{4}(rna_seq_rc)
            quals_rc = (0x30, 0x31, 0x32, 0x33)
            qmer_rc = Mycelia.Qualmer(kmer_rc, quals_rc)
            
            canon_qmer_rc = Mycelia.canonical(qmer_rc)
            Test.@test canon_qmer_rc.kmer == BioSequences.canonical(kmer_rc)  # Should be canonical
            # If UGCC gets reverse-complemented to GGCA, qualities should be reversed
            if canon_qmer_rc.kmer != kmer_rc
                Test.@test canon_qmer_rc.qualities == (0x33, 0x32, 0x31, 0x30)  # Qualities should be reversed
            else
                Test.@test canon_qmer_rc.qualities == quals_rc  # Qualities unchanged if already canonical
            end
        end
        
        Test.@testset "AA Canonical" begin
            # Amino acids don't have canonical form (no RC)
            aa_seq = BioSequences.LongAA("MWKL")
            kmer = Kmers.AAKmer{4}(aa_seq)
            quals = (0x20, 0x21, 0x22, 0x23)
            qmer = Mycelia.Qualmer(kmer, quals)
            
            canon_qmer = Mycelia.canonical(qmer)
            Test.@test canon_qmer === qmer  # Should be identical
        end
    end
    
    Test.@testset "Qualmer Iteration" begin
        Test.@testset "DNA Qualmers" begin
            # Create test sequence and quality
            dna_seq = BioSequences.LongDNA{4}("ACTGCATGCAATGC")
            quality = UInt8[0x20 + i for i in 1:length(dna_seq)]
            k = 5
            
            # Test qualmers_fw
            fw_qualmers = collect(Mycelia.qualmers_fw(dna_seq, quality, Val(k)))
            Test.@test length(fw_qualmers) == length(dna_seq) - k + 1
            Test.@test fw_qualmers[1][2] == 1  # First position
            Test.@test fw_qualmers[1][1].kmer == Kmers.Kmer{BioSequences.DNAAlphabet{4}, k}(dna_seq[1:k])
            Test.@test fw_qualmers[1][1].qualities == Tuple(quality[1:k])
            
            # Test qualmers_unambiguous with sequence containing ambiguous bases
            ambig_seq = BioSequences.LongDNA{4}("ACTGNATGC")
            ambig_qual = UInt8[0x20 + i for i in 1:length(ambig_seq)]
            
            unambig_qualmers = collect(Mycelia.qualmers_unambiguous(ambig_seq, ambig_qual, Val(3)))
            Test.@test length(unambig_qualmers) < length(ambig_seq) - 3 + 1  # Should skip ambiguous regions
            
            # Test qualmers_canonical
            canonical_qualmers = collect(Mycelia.qualmers_canonical(dna_seq, quality, Val(k)))
            Test.@test length(canonical_qualmers) == length(dna_seq) - k + 1
        end
        
        Test.@testset "RNA Qualmers" begin
            # Similar tests for RNA
            rna_seq = BioSequences.LongRNA{4}("ACUGCAUGCAAUGC")
            quality = UInt8[0x20 + i for i in 1:length(rna_seq)]
            k = 5
            
            fw_qualmers = collect(Mycelia.qualmers_fw(rna_seq, quality, Val(k)))
            Test.@test length(fw_qualmers) == length(rna_seq) - k + 1
        end
        
        Test.@testset "AA Qualmers" begin
            # Tests for amino acid sequences
            aa_seq = BioSequences.LongAA("MWKLVPGKEC")
            quality = UInt8[0x20 + i for i in 1:length(aa_seq)]
            k = 3
            
            fw_qualmers = collect(Mycelia.qualmers_fw(aa_seq, quality, Val(k)))
            Test.@test length(fw_qualmers) == length(aa_seq) - k + 1
            
            # AA sequences don't have ambiguous characters handled specially
            unambig_qualmers = collect(Mycelia.qualmers_unambiguous(aa_seq, quality, Val(k)))
            Test.@test length(unambig_qualmers) == length(aa_seq) - k + 1
        end
    end
    
    Test.@testset "FASTQ Record Integration" begin
        # Create a test FASTQ record
        header = "test_read"
        seq_str = "ACTGCATGCAATGC"
        quality_str = "IIIIHHHHGGGFFF"
        
        record = FASTX.FASTQ.Record(header, seq_str, quality_str)
        
        Test.@testset "Extract Qualmers from FASTQ" begin
            k = 5
            
            fw_qualmers = collect(Mycelia.qualmers_fw(record, k))
            Test.@test length(fw_qualmers) == length(seq_str) - k + 1
            
            canonical_qualmers = collect(Mycelia.qualmers_canonical(record, k))
            Test.@test length(canonical_qualmers) == length(seq_str) - k + 1
            
            # Test unambiguous filtering
            ambig_record = FASTX.FASTQ.Record(header, "ACTGNCTGA", "IIIIIIIII")
            unambig_qualmers = collect(Mycelia.qualmers_unambiguous(ambig_record, 3))
            Test.@test length(unambig_qualmers) < length(FASTX.sequence(ambig_record)) - 3 + 1
            
            # Test combined canonical and unambiguous
            unambig_canonical = collect(Mycelia.qualmers_unambiguous_canonical(record, k))
            Test.@test length(unambig_canonical) == length(seq_str) - k + 1
        end
    end
    
    Test.@testset "Edge Cases" begin
        # Test empty sequence
        empty_dna = BioSequences.LongDNA{4}("")
        empty_qual = UInt8[]
        Test.@test isempty(collect(Mycelia.qualmers_fw(empty_dna, empty_qual, Val(3))))
        
        # Test sequence shorter than k
        short_dna = BioSequences.LongDNA{4}("AC")
        short_qual = UInt8[0x20, 0x21]
        Test.@test isempty(collect(Mycelia.qualmers_fw(short_dna, short_qual, Val(3))))
        
        # Test sequence exactly length k
        exact_dna = BioSequences.LongDNA{4}("ACTG")
        exact_qual = UInt8[0x20, 0x21, 0x22, 0x23]
        exact_qualmers = collect(Mycelia.qualmers_fw(exact_dna, exact_qual, Val(4)))
        Test.@test length(exact_qualmers) == 1
    end

    Test.@testset "Simple Qualmer Examples" begin
        # Basic Qualmer construction for each alphabet
        dna_kmer = Kmers.DNAKmer{3}(BioSequences.LongDNA{4}("ACG"))
        rna_kmer = Kmers.RNAKmer{3}(BioSequences.LongRNA{4}("ACG"))
        aa_kmer  = Kmers.AAKmer{3}(BioSequences.LongAA("ACD"))
        qual = (0x01, 0x02, 0x03)

        dna_q = Mycelia.Qualmer(dna_kmer, qual)
        rna_q = Mycelia.Qualmer(rna_kmer, qual)
        aa_q  = Mycelia.Qualmer(aa_kmer,  qual)

        Test.@test dna_q.kmer == dna_kmer && dna_q.qualities == qual
        Test.@test rna_q.kmer == rna_kmer && rna_q.qualities == qual
        Test.@test aa_q.kmer  == aa_kmer  && aa_q.qualities  == qual

        # Canonical with reverse complement for DNA and RNA
        rc_dna = Kmers.DNAKmer{3}(BioSequences.LongDNA{4}("GAT"))
        rc_q   = (0x10, 0x11, 0x12)
        canon_dna = Mycelia.canonical(Mycelia.Qualmer(rc_dna, rc_q))
        Test.@test canon_dna.kmer == BioSequences.canonical(rc_dna)
        Test.@test canon_dna.qualities == (0x12, 0x11, 0x10)

        rc_rna = Kmers.RNAKmer{3}(BioSequences.LongRNA{4}("GAU"))
        canon_rna = Mycelia.canonical(Mycelia.Qualmer(rc_rna, rc_q))
        Test.@test canon_rna.kmer == BioSequences.canonical(rc_rna)
        Test.@test canon_rna.qualities == (0x12, 0x11, 0x10)

        # Iterator checks with small sequences
        seq = BioSequences.LongDNA{4}("ACGT")
        q   = UInt8[1,2,3,4]

        fw = collect(Mycelia.qualmers_fw(seq, q, Val(2)))
        Test.@test length(fw) == 3
        # Extract the actual k-mer from the sequence for comparison
        expected_kmer = first(Kmers.FwKmers{BioSequences.DNAAlphabet{4}, 2}(seq))
        Test.@test fw[1] == (Mycelia.Qualmer(expected_kmer, (0x01,0x02)), 1)

        unambig = collect(Mycelia.qualmers_unambiguous(seq, q, Val(2)))
        Test.@test length(unambig) == 3
        # UnambiguousDNAMers returns 2-bit k-mers, so create the expected k-mer using that iterator
        third_unambig_kmer = collect(Kmers.UnambiguousDNAMers{2}(seq))[3][1]  # (kmer, pos) tuple
        Test.@test unambig[3] == (Mycelia.Qualmer(third_unambig_kmer, (0x03,0x04)), 3)

        canon = collect(Mycelia.qualmers_canonical(seq, q, Val(2)))
        Test.@test length(canon) == 3
        # The canonical of GT should be AC with reversed qualities
        third_kmer_fw = collect(Kmers.FwKmers{BioSequences.DNAAlphabet{4}, 2}(seq))[3]
        canonical_third = BioSequences.canonical(third_kmer_fw)
        Test.@test canon[3] == (Mycelia.Qualmer(canonical_third, (0x04,0x03)), 3)
    end
    
    Test.@testset "Probability Calculations" begin
        Test.@testset "PHRED Conversion" begin
            # Test PHRED to probability conversion
            Test.@test Mycelia.phred_to_probability(UInt8(10)) ≈ 0.9 atol=0.01  # Q10 = 90% accuracy
            Test.@test Mycelia.phred_to_probability(UInt8(20)) ≈ 0.99 atol=0.001  # Q20 = 99% accuracy
            Test.@test Mycelia.phred_to_probability(UInt8(30)) ≈ 0.999 atol=0.0001  # Q30 = 99.9% accuracy
            
            # Test PHRED to error probability conversion
            Test.@test Mycelia.phred_to_error_probability(UInt8(10)) ≈ 0.1 atol=0.01
            Test.@test Mycelia.phred_to_error_probability(UInt8(20)) ≈ 0.01 atol=0.001
        end
        
        Test.@testset "Qualmer Probability" begin
            # Create test qualmer with known quality scores
            dna_seq = BioSequences.LongDNA{4}("ACTG")
            kmer = Kmers.DNAKmer{4}(dna_seq)
            high_quals = (UInt8(30), UInt8(30), UInt8(30), UInt8(30))  # High quality
            low_quals = (UInt8(10), UInt8(10), UInt8(10), UInt8(10))   # Low quality
            
            high_qmer = Mycelia.Qualmer(kmer, high_quals)
            low_qmer = Mycelia.Qualmer(kmer, low_quals)
            
            # High quality should have higher correctness probability
            high_prob = Mycelia.qualmer_correctness_probability(high_qmer)
            low_prob = Mycelia.qualmer_correctness_probability(low_qmer)
            Test.@test high_prob > low_prob
            Test.@test high_prob > 0.99  # Should be very high for Q30
            Test.@test low_prob < 0.9   # Should be lower for Q10
        end
        
        Test.@testset "Joint Probability" begin
            # Create multiple observations of the same k-mer
            dna_seq = BioSequences.LongDNA{4}("ACTG")
            kmer = Kmers.DNAKmer{4}(dna_seq)
            quals1 = (UInt8(30), UInt8(25), UInt8(20), UInt8(15))
            quals2 = (UInt8(25), UInt8(30), UInt8(25), UInt8(20))
            
            qmer1 = Mycelia.Qualmer(kmer, quals1)
            qmer2 = Mycelia.Qualmer(kmer, quals2)
            
            # Test joint probability calculation
            joint_prob = Mycelia.joint_qualmer_probability([qmer1, qmer2])
            Test.@test joint_prob >= 0.0
            Test.@test joint_prob <= 1.0
            
            # Test position-wise joint probability
            pos_joint_prob = Mycelia.position_wise_joint_probability([qmer1, qmer2])
            Test.@test pos_joint_prob >= 0.0
            Test.@test pos_joint_prob <= 1.0
            
            # Test empty vector
            Test.@test Mycelia.joint_qualmer_probability(Mycelia.Qualmer[]) == 0.0
        end
    end
    
    Test.@testset "Graph Construction and Analysis" begin
        # Create test FASTQ records
        records = [
            FASTX.FASTQ.Record("read1", "ACTGCATGCA", "IIIIIIIIII"),
            FASTX.FASTQ.Record("read2", "TGCATGCAAT", "HHHHHHHHHH"),
            FASTX.FASTQ.Record("read3", "CATGCAATGC", "GGGGGGGGGG")
        ]
        
        Test.@testset "Basic Graph Construction" begin
            # Build a qualmer graph
            graph = Mycelia.build_qualmer_graph(records, k=4)
            Test.@test isa(graph, MetaGraphsNext.MetaGraph)
            
            # Check that graph has vertices
            vertices = collect(MetaGraphsNext.labels(graph))
            Test.@test length(vertices) > 0
            
            # Get basic statistics
            stats = Mycelia.get_qualmer_statistics(graph)
            Test.@test haskey(stats, "num_vertices")
            Test.@test haskey(stats, "num_edges")
            Test.@test haskey(stats, "mean_coverage")
            Test.@test stats["num_vertices"] >= 0
            Test.@test stats["mean_coverage"] >= 0.0
        end
        
        Test.@testset "Quality Metrics" begin
            graph = Mycelia.build_qualmer_graph(records, k=3, min_coverage=1)
            
            # Calculate assembly quality metrics
            quality_metrics = Mycelia.calculate_assembly_quality_metrics(graph)
            Test.@test quality_metrics.mean_coverage >= 0.0
            Test.@test quality_metrics.mean_quality >= 0.0
            Test.@test quality_metrics.mean_confidence >= 0.0
            Test.@test quality_metrics.low_confidence_fraction >= 0.0
            Test.@test quality_metrics.low_confidence_fraction <= 1.0
            Test.@test quality_metrics.total_kmers >= 0
        end
        
        Test.@testset "Error Identification" begin
            graph = Mycelia.build_qualmer_graph(records, k=3, min_coverage=1)
            
            # Identify potential errors with lenient thresholds
            error_labels = Mycelia.identify_potential_errors(graph, 
                min_coverage=1, min_quality=5.0, min_confidence=0.1)
            Test.@test isa(error_labels, Vector)
            # All labels should be k-mers (Kmers.Kmer types)
            if !isempty(error_labels)
                Test.@test all(v -> isa(v, Kmers.Kmer), error_labels)
            end
        end
    end
    
    Test.@testset "Additional Edge Cases and Coverage" begin
        Test.@testset "Quality Score Edge Cases" begin
            # Test with extreme quality scores
            Test.@test Mycelia.phred_to_probability(UInt8(0)) ≈ 0.0 atol=0.01  # Very low quality
            Test.@test Mycelia.phred_to_probability(UInt8(60)) ≈ 1.0 atol=0.000002  # Very high quality
            
            # Test with mixed quality qualmer
            dna_seq = BioSequences.LongDNA{4}("ACTG")
            kmer = Kmers.DNAKmer{4}(dna_seq)
            mixed_quals = (UInt8(0), UInt8(60), UInt8(0), UInt8(60))
            mixed_qmer = Mycelia.Qualmer(kmer, mixed_quals)
            
            prob = Mycelia.qualmer_correctness_probability(mixed_qmer)
            Test.@test prob >= 0.0
            Test.@test prob <= 1.0
        end
        
        Test.@testset "Qualmer Constructor Edge Cases" begin
            # Test constructor with AbstractVector input
            dna_seq = BioSequences.LongDNA{4}("ACG")
            kmer = Kmers.DNAKmer{3}(dna_seq)
            qual_vec = [30, 25, 20]  # Int vector
            
            qmer = Mycelia.Qualmer(kmer, qual_vec)
            Test.@test qmer.qualities == (0x1e, 0x19, 0x14)  # Converted to UInt8
            
            # Test assertion error for mismatched lengths
            Test.@test_throws AssertionError Mycelia.Qualmer(kmer, [30, 25])  # Too short
        end
        
        Test.@testset "Empty Graph Handling" begin
            # Test statistics on empty graph
            empty_graph = MetaGraphsNext.MetaGraph(
                MetaGraphsNext.DiGraph(),
                label_type=Kmers.DNAKmer{3},
                vertex_data_type=Mycelia.QualmerVertexData,
                edge_data_type=Mycelia.QualmerEdgeData,
                weight_function=Mycelia.edge_data_weight,
                default_weight=0.0
            )
            
            stats = Mycelia.get_qualmer_statistics(empty_graph)
            Test.@test stats["num_vertices"] == 0
            Test.@test stats["num_edges"] == 0
            Test.@test stats["mean_coverage"] == 0.0
        end
    end
end
