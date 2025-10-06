# ```bash
# julia --project=. --color=yes -e 'include("test/2_preprocessing_qc/constants.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run from the Mycelia base directory:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/constants.jl", "test/2_preprocessing_qc", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia

# Constants and regex tests
Test.@testset "Constants" begin
    Test.@testset "FASTQ regex" begin
        hypothetical_fastq_files = [
            "sample1.fastq",
            "sample2.fq",
        ]
        for hypothetical_fastq_file in hypothetical_fastq_files
            Test.@test occursin(Mycelia.FASTQ_REGEX, hypothetical_fastq_file)
            Test.@test occursin(Mycelia.FASTQ_REGEX, hypothetical_fastq_file * ".gz")
        end
        invalid_fastq_files = [
            "sample.fast",
            "sample.fastq1",
            "notfastq.txt",
            "sample.fastq.gz.old",
            "reads.fq.gz.bak",
        ]
        for invalid_fastq_file in invalid_fastq_files
            Test.@test !occursin(Mycelia.FASTQ_REGEX, invalid_fastq_file)
        end
    end
    Test.@testset "FASTA regex" begin
        hypothetical_fasta_files = [
            "genome1.fasta",
            "fasta-sequences.fas",
            "fasta.fa",
            "fasta-nucleic-acid.fna",
            "open-reading-frames.ffn",
            "fasta-amino-acid.faa",
            "multi-protein-fasta.mpfa",
            "transcriptome8.frn",
        ]
        for hypothetical_fasta_file in hypothetical_fasta_files
            Test.@test occursin(Mycelia.FASTA_REGEX, hypothetical_fasta_file)
            Test.@test occursin(Mycelia.FASTA_REGEX, hypothetical_fasta_file * ".gz")
        end
        invalid_fasta_files = [
            "genome.fasta1",
            "file.fasta.txt",
            "notfasta.fa.gz.old",
            "genome.fa.bz2",
            "fasta.fna.gz.tmp",
        ]
        for invalid_fasta_file in invalid_fasta_files
            Test.@test !occursin(Mycelia.FASTA_REGEX, invalid_fasta_file)
        end
    end
    Test.@testset "VCF regex" begin
        hypothetical_vcf_files = [
            "variants1.vcf",
            "variants2.vcf.gz",
        ]
        for hypothetical_vcf_file in hypothetical_vcf_files
            Test.@test occursin(Mycelia.VCF_REGEX, hypothetical_vcf_file)
        end
        invalid_vcf_files = [
            "variants.vcf1",
            "vcf_variants.txt",
            "sample.vcf.gz.old",
            "genotypes.vcf.bz2",
            "variants.vcf.gz.bak",
        ]
        for invalid_vcf_file in invalid_vcf_files
            Test.@test !occursin(Mycelia.VCF_REGEX, invalid_vcf_file)
        end
    end
    Test.@testset "XAM regex" begin
        hypothetical_xam_files = [
            "alignment1.sam",
            "alignment2.sam.gz",
            "alignment3.bam",
            "alignment4.cram"
        ]
        for hypothetical_xam_file in hypothetical_xam_files
            Test.@test occursin(Mycelia.XAM_REGEX, hypothetical_xam_file)
        end
        invalid_xam_files = [
            "alignment.sam1",
            "alignment.bam.gz",
            "xam_alignment.txt",
            "alignment.sam.gz.old",
            "alignment.cram.gz",
            "alignment.bam1",
        ]
        for invalid_xam_file in invalid_xam_files
            Test.@test !occursin(Mycelia.XAM_REGEX, invalid_xam_file)
        end
    end

    Test.@testset "Alphabet constants" begin
        ## Test that ACGTN and ACGUN constants are defined
        Test.@test isdefined(Mycelia, :ACGTN_DNA_SYMBOLS)
        Test.@test isdefined(Mycelia, :ACGUN_RNA_SYMBOLS)
        Test.@test isdefined(Mycelia, :ACGTN_DNA_CHARSET)
        Test.@test isdefined(Mycelia, :ACGUN_RNA_CHARSET)
        
        ## Test exact values of alphabet constants
        Test.@testset "Constant values" begin
            ## Test that ACGTN contains exactly A, C, G, T, N (case insensitive)
            expected_acgtn = Set(['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'])
            Test.@test Mycelia.ACGTN_DNA_CHARSET == expected_acgtn
            
            ## Test that ACGUN contains exactly A, C, G, U, N (case insensitive)
            expected_acgun = Set(['A', 'C', 'G', 'U', 'N', 'a', 'c', 'g', 'u', 'n'])
            Test.@test Mycelia.ACGUN_RNA_CHARSET == expected_acgun
            
            ## Test unambiguous DNA contains exactly A, C, G, T (case insensitive)
            expected_unambiguous_dna = Set(['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'])
            Test.@test Mycelia.UNAMBIGUOUS_DNA_CHARSET == expected_unambiguous_dna
            
            ## Test unambiguous RNA contains exactly A, C, G, U (case insensitive)
            expected_unambiguous_rna = Set(['A', 'C', 'G', 'U', 'a', 'c', 'g', 'u'])
            Test.@test Mycelia.UNAMBIGUOUS_RNA_CHARSET == expected_unambiguous_rna
            
            ## Test that ambiguous alphabets include gaps and additional symbols
            Test.@test '-' ∈ Mycelia.AMBIGUOUS_DNA_CHARSET
            Test.@test '-' ∈ Mycelia.AMBIGUOUS_RNA_CHARSET
            Test.@test '-' ∈ Mycelia.AMBIGUOUS_AA_CHARSET
            
            ## Test that ambiguous alphabets include standard IUPAC codes
            Test.@test 'R' ∈ Mycelia.AMBIGUOUS_DNA_CHARSET  # A or G
            Test.@test 'Y' ∈ Mycelia.AMBIGUOUS_DNA_CHARSET  # C or T
            Test.@test 'R' ∈ Mycelia.AMBIGUOUS_RNA_CHARSET  # A or G
            Test.@test 'Y' ∈ Mycelia.AMBIGUOUS_RNA_CHARSET  # C or U
        end

        ## Test subset relationships
        Test.@testset "Subset relationships" begin
            ## DNA hierarchy: unambiguous ⊆ ACGTN ⊆ ambiguous
            Test.@test issubset(Mycelia.UNAMBIGUOUS_DNA_CHARSET, Mycelia.ACGTN_DNA_CHARSET)
            Test.@test issubset(Mycelia.ACGTN_DNA_CHARSET, Mycelia.AMBIGUOUS_DNA_CHARSET)
            Test.@test issubset(Mycelia.UNAMBIGUOUS_DNA_CHARSET, Mycelia.AMBIGUOUS_DNA_CHARSET)
            
            ## RNA hierarchy: unambiguous ⊆ ACGUN ⊆ ambiguous
            Test.@test issubset(Mycelia.UNAMBIGUOUS_RNA_CHARSET, Mycelia.ACGUN_RNA_CHARSET)
            Test.@test issubset(Mycelia.ACGUN_RNA_CHARSET, Mycelia.AMBIGUOUS_RNA_CHARSET)
            Test.@test issubset(Mycelia.UNAMBIGUOUS_RNA_CHARSET, Mycelia.AMBIGUOUS_RNA_CHARSET)
            
            ## AA hierarchy: unambiguous ⊆ ambiguous
            Test.@test issubset(Mycelia.UNAMBIGUOUS_AA_CHARSET, Mycelia.AMBIGUOUS_AA_CHARSET)
        end

        ## Test alphabet size ordering (smallest to largest sets)
        ## This verifies the detect_alphabet search order is from most specific to most general
        Test.@testset "Size ordering" begin
            ## Unambiguous DNA and RNA should have the same size (4 bases each, case insensitive = 8 chars)
            Test.@test length(Mycelia.UNAMBIGUOUS_DNA_CHARSET) == length(Mycelia.UNAMBIGUOUS_RNA_CHARSET) == 8
            
            ## ACGTN should be larger than unambiguous DNA but smaller than full ambiguous
            Test.@test length(Mycelia.UNAMBIGUOUS_DNA_CHARSET) < length(Mycelia.ACGTN_DNA_CHARSET)
            Test.@test length(Mycelia.ACGTN_DNA_CHARSET) < length(Mycelia.AMBIGUOUS_DNA_CHARSET)
            
            ## ACGUN should be larger than unambiguous RNA but smaller than full ambiguous
            Test.@test length(Mycelia.UNAMBIGUOUS_RNA_CHARSET) < length(Mycelia.ACGUN_RNA_CHARSET)
            Test.@test length(Mycelia.ACGUN_RNA_CHARSET) < length(Mycelia.AMBIGUOUS_RNA_CHARSET)
            
            ## Verify that AA alphabet is larger than nucleotide alphabets (more amino acids than nucleotides)
            Test.@test length(Mycelia.UNAMBIGUOUS_DNA_CHARSET) < length(Mycelia.UNAMBIGUOUS_AA_CHARSET)
            Test.@test length(Mycelia.UNAMBIGUOUS_RNA_CHARSET) < length(Mycelia.UNAMBIGUOUS_AA_CHARSET)
            
            ## Full ambiguous alphabets should be the largest
            Test.@test length(Mycelia.UNAMBIGUOUS_AA_CHARSET) < length(Mycelia.AMBIGUOUS_AA_CHARSET)
            
            ## Specific size expectations
            Test.@test length(Mycelia.ACGTN_DNA_CHARSET) == 10  # A,C,G,T,N × 2 cases
            Test.@test length(Mycelia.ACGUN_RNA_CHARSET) == 10  # A,C,G,U,N × 2 cases
        end
    end
end
