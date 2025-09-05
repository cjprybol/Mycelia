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
end
