# Constants and regex tests
@testset "Constants" begin
    @testset "FASTQ regex" begin
        hypothethical_fastq_files = [
            "sample1.fastq",
            "sample2.fq",
        ]
        for hypothethical_fastq_file in hypothethical_fastq_files
            @test occursin(Mycelia.FASTQ_REGEX, hypothethical_fastq_file)
            @test occursin(Mycelia.FASTQ_REGEX, hypothethical_fastq_file * ".gz")
        end
    end
    @testset "FASTA regex" begin
        hypothethical_fasta_files = [
            "genome1.fasta",
            "fasta-sequences.fas",
            "fasta.fa",
            "fasta-nucleic-acid.fna",
            "open-reading-frames.ffn",
            "fasta-amino-acid.faa",
            "multi-protein-fasta.mpfa",
            "transcriptome8.frn",
        ]
        for hypothethical_fasta_file in hypothethical_fasta_files
            @test occursin(Mycelia.FASTA_REGEX, hypothethical_fasta_file)
            @test occursin(Mycelia.FASTA_REGEX, hypothethical_fasta_file * ".gz")
        end
    end
    @testset "VCF regex" begin
        hypothethical_vcf_files = [
            "variants1.vcf",
            "variants2.vcf.gz",
        ]
        for hypothethical_vcf_file in hypothethical_vcf_files
            @test occursin(Mycelia.VCF_REGEX, hypothethical_vcf_file)
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
