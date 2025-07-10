# Constants and regex tests
@testset "Constants" begin
    @testset "FASTQ regex" begin
        hypothetical_fastq_files = [
            "sample1.fastq",
            "sample2.fq",
        ]
        for hypothetical_fastq_file in hypothetical_fastq_files
            @test occursin(Mycelia.FASTQ_REGEX, hypothetical_fastq_file)
            @test occursin(Mycelia.FASTQ_REGEX, hypothetical_fastq_file * ".gz")
        end
    end
    @testset "FASTA regex" begin
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
            @test occursin(Mycelia.FASTA_REGEX, hypothetical_fasta_file)
            @test occursin(Mycelia.FASTA_REGEX, hypothetical_fasta_file * ".gz")
        end
    end
    @testset "VCF regex" begin
        hypothetical_vcf_files = [
            "variants1.vcf",
            "variants2.vcf.gz",
        ]
        for hypothetical_vcf_file in hypothetical_vcf_files
            @test occursin(Mycelia.VCF_REGEX, hypothetical_vcf_file)
        end
    end
    @testset "XAM regex" begin
        hypothetical_xam_files = [
            "alignment1.sam",
            "alignment2.sam.gz",
            "alignment3.bam",
            "alignment4.cram"
        ]
        for hypothetical_xam_file in hypothetical_xam_files
            @test occursin(Mycelia.XAM_REGEX, hypothetical_xam_file)
        end
    end
end
