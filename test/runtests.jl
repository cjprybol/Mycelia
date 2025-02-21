import Pkg
if isinteractive()
    Pkg.activate("..")
end

using Revise
using Test
import Mycelia

import FASTX
import DocStringExtensions
import BioSequences
import Dates
import Random
import SHA
import CSV
import uCSV

const SEED = 42
const phiX174_accession_id = "NC_001422.1"
const phiX174_assembly_id = "GCF_000819615.1"

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

# Simulation: for multi-omics data generation for benchmarking or testing.
@testset "FASTA simulation and acquisition" begin
    @testset "dna record" begin
        dna_record = Mycelia.random_fasta_record(moltype=:DNA, seed=SEED, L = 10)
        # @test FASTX.identifier(dna_record) == "e6614e938383857305089158322604fb9c07d0aac05567844d01f1517ce82ddc"
        @test FASTX.sequence(dna_record) == "CCGCCGCTCA"
    end

    @testset "rna record" begin
        rna_record = Mycelia.random_fasta_record(moltype=:RNA, seed=SEED, L = 10)
        # @test FASTX.identifier(rna_record) == "a4ebfc1dd4051dd9db816c42c5fb8f528ae0b5387d7a8ccad67995e2c31ee48d"
        @test FASTX.sequence(rna_record) == "CCGCCGCUCA"
    end

    @testset "aa record" begin
        aa_record = Mycelia.random_fasta_record(moltype=:AA, seed=SEED, L = 10)
        # @test FASTX.identifier(aa_record) == "cdb40a05b841c44e3757de6a54a1fd011536fce187c7e8c9c1ed36c415db8344"
        @test FASTX.sequence(aa_record) == "VATAGWWITI"
    end
    
    @testset "virus phiX174" begin
        genome_result = Mycelia.download_genome_by_accession(accession=phiX174_accession_id)
        @test basename(genome_result) == phiX174_accession_id * ".fna.gz"
        @test Mycelia.get_base_extension(genome_result) == ".fna.gz"
        rm(genome_result)
        # @test Mycelia.sha256_file(genome_result) == "765354d42319e4a350c93b09260bf911864263653f828ad97b5c35eed950591a"
        
        phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="gff3,rna,cds,protein,genome,seq-report")
        @test basename(phiX174_assembly_dataset.genome) == phiX174_assembly_id * "_ViralProj14015_genomic.fna"
        #TODO add tests for remaining items? do that elsewhere?
        @test Mycelia.get_base_extension(phiX174_assembly_dataset.genome) == ".fna"
        @test Mycelia.get_base_extension(phiX174_assembly_dataset.protein) == ".faa"
        rm(phiX174_assembly_id, recursive=true)
    end

    @testset "bacteria-like" begin
        @test 1 + 1 == 2
    end

    @testset "protist-like" begin
        @test 1 + 1 == 2
    end

    @testset "fungi-like" begin
        @test 1 + 1 == 2
    end

    @testset "plant-like" begin
        @test 1 + 1 == 2
    end

    @testset "animal-like" begin
        @test 1 + 1 == 2
    end

    @testset "microbiome" begin
        @test 1 + 1 == 2
    end
end

@testset "sequence IO" begin
    @testset "detect alphabet type" begin
        @test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA))) == :DNA
        @test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA))) == :RNA
        @test Mycelia.detect_alphabet(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA))) == :AA
    end

    @testset "autoconvert sequences" begin
        @test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA)))) <: BioSequences.LongDNA{4}
        @test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA)))) <: BioSequences.LongRNA{4}
        @test typeof(Mycelia.convert_sequence(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA)))) <: BioSequences.LongAA
    end

    @testset "detect sequence extension" begin
        @test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:DNA))) == ".fna"
        @test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:RNA))) == ".frn"
        @test Mycelia.detect_sequence_extension(FASTX.sequence(Mycelia.random_fasta_record(L=100, moltype=:AA))) == ".faa"
    end
end

@testset "SRA downloading" begin
    outdir = mkpath("test-sra-tools")
    @testset "download SRA short read data" begin
        # https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR31271002&display=metadata
        srr_identifier = "SRR31271002"
        prefetch_results = Mycelia.prefetch(SRR=srr_identifier, outdir=outdir)
        @test prefetch_results.directory == "test-sra-tools/$(srr_identifier)"
        @test prefetch_results.archive == "test-sra-tools/$(srr_identifier)/$(srr_identifier).sra"
        @test filesize(prefetch_results.archive) == 18295489
        
        fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
        @test fasterq_dump_result.forward_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier)_1.fastq.gz"
        @test fasterq_dump_result.reverse_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier)_2.fastq.gz"
        @test ismissing(fasterq_dump_result.unpaired_reads)
        @test filesize(fasterq_dump_result.forward_reads) == 11637534
        @test filesize(fasterq_dump_result.reverse_reads) == 11908008
    end
    
    @testset "download SRA long read data" begin
        # https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR31812976&display=metadata
        srr_identifier = "SRR31812976"
        prefetch_results = Mycelia.prefetch(SRR=srr_identifier, outdir=outdir)
        @test prefetch_results.directory == "$(outdir)/$(srr_identifier)"
        @test prefetch_results.archive == "$(outdir)/$(srr_identifier)/$(srr_identifier).sra"
        @test filesize(prefetch_results.archive) == 33392268

        fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
        @test ismissing(fasterq_dump_result.forward_reads)
        @test ismissing(fasterq_dump_result.reverse_reads)
        @test fasterq_dump_result.unpaired_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier).fastq.gz"
        @test filesize(fasterq_dump_result.unpaired_reads) == 34353450
    end
    rm(outdir, recursive=true)
end

# PRE-PROCESSING & READ QC Tests
@testset "Preprocessing" begin
    @testset "FASTX stats" begin
        # https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR31812976&display=metadata
        srr_identifier = "SRR31812976"
        outdir = mkpath("fastx-stats-test")
        fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
        @test fasterq_dump_result.unpaired_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier).fastq.gz"
        table = Mycelia.fastx_stats(fasterq_dump_result.unpaired_reads)
        io = IOBuffer()
        CSV.write(io, table)
        @test String(take!(io)) == 
        """
        file,format,type,num_seqs,sum_len,min_len,avg_len,max_len,Q1,Q2,Q3,sum_gap,N50,N50_num,Q20(%),Q30(%),AvgQual,GC(%),sum_n,N90
        fastx-stats-test/SRR31812976/SRR31812976.fastq.gz,FASTQ,DNA,58982,34921275,80,592.1,2315,341.0,544.0,811.0,0,743,608,50.26,9.75,11.4,39.62,0,329
        """
        rm(outdir, recursive=true)
    end

    @testset "fastx2normalized_table" begin
        genome_result = Mycelia.download_genome_by_accession(accession=phiX174_accession_id)
        fastx_table = Mycelia.fastx2normalized_table(genome_result)
        io = IOBuffer()
        CSV.write(io, fastx_table, delim='\t')
        @test String(take!(io)) == 
        """
        fastx_path	fastx_sha256	record_identifier	record_description	record_sha256	record_quality	record_alphabet	record_type	mean_record_quality	median_record_quality	record_length	record_sequence
        NC_001422.1.fna.gz	ae3c70dc4d180aab0aa8d6ab81c0048998575c05e402639ffe6b99e82815f953	NC_001422.1	NC_001422.1 Escherichia phage phiX174, complete genome	97038c7e1edea2297667d7f0426ba942b322c74cb30e072ec66ba47f9c0448d0		ACGT	DNA			5386	GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGATGCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTTTCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTTCGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTACGGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTACGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAGTGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAAGGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGCCGCTTTTCTTGGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACAACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTTTCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAACCTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTTGATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACAGACCTATAAACATTCTGTGCCGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGACTAAAGAGATTCAGTACCTTAACGCTAAAGGTGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTGTATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGAAGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCGCCACCATGATTATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGATTGAGTGTGAGGTTATAACGCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGCTTAGGAGTTTAATCATGTTTCAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGTTCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTATATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTGAATGGTCGCCATGATGGTGGTTATTATACCGTCAAGGACTGTGTGACTATTGACGTCCTTCCCCGTACGCCGGGCAATAACGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAATCAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTGCTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAGGCGGTCAAAAAGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCTGGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACCCTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGATAAGTTGCTTGATTTGGTTGGACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTATCTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGGTTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTTCACGCCAGAATACGAAAGACCAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTATGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCAAACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGACTTAGTTCATCAGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTTCTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAGCTGTTGCCGATACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCGTCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTCCAAATCTTGGAGGCTTTTTTATGGTTCGTTCTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTATTGAGGCTTGTGGCATTTCTACTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGCATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATGTTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGAATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTATAGACCACCGCCCCGAAGGGGACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCCCTCTTAAGGATATTCGCGATGAGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATTGCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAGGCTCATGCTGATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTTATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGGTCGCAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGCCGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTCCTACAGGTAGCGTTGACCCTAATTTTGGTCGTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAGCCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATATGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACTTCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGCAGAACGCAAAAAGAGAGATGAGATTGAGGCTGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACCTGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA
        """
        rm(genome_result)
    end
    
    @testset "Read Quality Control" begin
        @test true # https://github.com/OpenGene/fastp
        @test true # https://github.com/FelixKrueger/TrimGalore
        @test true # https://github.com/rrwick/Filtlong
        @test true # https://github.com/OpenGene/fastplong
        @test true # https://github.com/wdecoster/chopper
        # Example: test that adapter trimming and quality filtering work.
        # result = MyceliaAssembly.preprocess_reads("test_data/reads.fastq")
        # @test length(result.filtered_reads) > 0
        # @test result.mean_quality ≥ 30
        @test true  # placeholder
    end
    @testset "Read Statistics" begin
        # Example: test that estimated community composition is within expected bounds.
        # comp = MyceliaAssembly.analyze_community("test_data/reads.fastq")
        # @test 0.8 <= comp["expected_coverage"] <= 1.2
        @test true  # placeholder
    end
end

@testset "tool integration" begin
    @testset "padloc" begin
        ecoli_k12_accession = "GCF_000005845.2"
        result = Mycelia.ncbi_genome_download_accession(accession=ecoli_k12_accession, include_string="genome")
        # @show result.genome
        padloc_result = Mycelia.run_padloc(result.genome)
        @test isfile(padloc_result)
        rm(ecoli_k12_accession, recursive=true)
    end
end

# Data ingestion & normalization: for multi-omics data loading, sanity checks, QC.
# Data Ingestion & Normalization
# - Loading multi-omics data
# - Sanity checks and quality control (QC)
# - Subsampling reads (e.g., subsample_reads_seqtk)
# @testset "Data ingestion & normalization" begin
#     @testset "fasta2normalized_table" begin
#         temp_dir = mktempdir()
#         for alphabet in [:DNA, :AA, :RNA]
#             fasta_record = Mycelia.random_fasta_record(moltype=alphabet, seed=SEED, L = 100)
#             extension = Mycelia.detect_sequence_extension(fasta_record)
#             fasta_file = joinpath(temp_dir, "test" * extension)
#             @test Mycelia.get_base_extension(fasta_file) == extension
#             Mycelia.write_fasta(outfile=fasta_file, records=[fasta_record])
#             result_file = Mycelia.fasta2normalized_table(fasta_file, normalize_name=true)
#             if alphabet == :DNA
#                 @test basename(result_file) == "b58f92d3ad303b329300debba8a8b4e9ebe6fb68fabc1ff7b2aefc1e267b64fb" * extension * ".tsv.gz"
#             elseif alphabet == :RNA
#                 @test basename(result_file) == "dad9c65109b858e38bdc34463807729136ecc102bfbd0247dc00bd4811098dab" * extension * ".tsv.gz"
#             elseif alphabet == :AA
#                 @test basename(result_file) == "5abeeb482b1af9097f71b81d53c72f98c0a69716b6fbfda1a8a855c33b3b0888" * extension * ".tsv.gz"
#             end
#         end
#         rm(temp_dir, recursive=true)
#     end
# end

@testset "FASTQ simulation" begin
    @testset "Illumina" begin
        @test 1 + 1 == 2
    end

    @testset "Ultima" begin
        @test 1 + 1 == 2
    end

    @testset "Nanopore" begin
        @test 1 + 1 == 2
    end

    @testset "PacBio" begin
        @test 1 + 1 == 2
    end

    @testset "multi-entity, even coverage" begin
        @test 1 + 1 == 2
    end

    @testset "multi-entity, log-distributed coverage" begin
        @test 1 + 1 == 2
    end
end

# assembly modules
@testset "assembly modules" begin
    @testset "1. Pre‐processing & Read QC" begin
    end
    @testset "2. k‑mer Analysis" begin
    end
    @testset "3. Hybrid Assembly (Mutual Support Strategy)" begin
    end
    @testset "4. Assembly merging" begin
    end
    @testset "5. Polishing & Error Correction" begin
    end
    @testset "6. Strain resolution" begin
    end
    @testset "7. Validation & Quality Control" begin
    end
end

# PANGENOME/PANPROTEOME AND K-MER ANALYSIS Tests
@testset "Reference Graph and K-mer Analysis" begin
    @testset "Pangenome Construction" begin
        # Example: verify that a reference pangenome graph is built
        # pg_graph = MyceliaAssembly.build_pangenome(["ref1.fasta", "ref2.fasta"])
        # @test typeof(pg_graph) <: AbstractGraph
        @test true  # placeholder
    end
    @testset "Optimal K-mer Selection" begin
        # Example: determine the best k-mer length from given reads
        # best_k = MyceliaAssembly.select_optimal_k("test_data/reads.fastq")
        # @test best_k isa Int
        # @test 21 ≤ best_k ≤ 127
        @test true  # placeholder
    end
end

# HYBRID ASSEMBLY Tests
@testset "Hybrid Assembly" begin
    @testset "Assembly Core" begin
        # Example: run hybrid assembly on a small simulated dataset.
        # assembly = MyceliaAssembly.hybrid_assemble("test_data/reads.fastq"; long_read="test_data/long.fastq")
        # @test length(assembly.contigs) > 0
        # @test assembly.N50 > 5000
        @test true  # placeholder
    end
    @testset "Contig Overlap Graph Integrity" begin
        # Example: check that the overlap graph correctly represents strain variants.
        # graph = MyceliaAssembly.get_overlap_graph(assembly)
        # @test MyceliaAssembly.validate_overlap_graph(graph) == true
        @test true  # placeholder
    end
end

# ASSEMBLY MERGE Tests
@testset "Assembly Merging" begin
    @testset "Contig Merging" begin
        # Example: test that QuickMerge improves contiguity compared to input assemblies.
        # merged = MyceliaAssembly.quick_merge("assembly1.fasta", "assembly2.fasta")
        # @test merged.N50 > max(assembly1.N50, assembly2.N50)
        @test true  # placeholder
    end
end

# POLISHING Tests
@testset "Assembly Polishing" begin
    @testset "Error Correction" begin
        # Example: simulate a polishing step and verify base accuracy improvement.
        # polished = MyceliaAssembly.polish_assembly("raw_assembly.fasta", reads="test_data/reads.fastq")
        # accuracy_before = MyceliaAssembly.evaluate_accuracy("raw_assembly.fasta", "reference.fasta")
        # accuracy_after = MyceliaAssembly.evaluate_accuracy(polished, "reference.fasta")
        # @test accuracy_after > accuracy_before
        @test true  # placeholder
    end
end

# STRAIN RESOLUTION Tests
@testset "Strain Resolution" begin
    @testset "Strain-aware Reassembly" begin
        # Example: test that strain-specific contigs or haplotigs are generated.
        # strains = MyceliaAssembly.resolve_strains("polished_assembly.fasta", reads="test_data/long.fastq")
        # @test length(strains) >= 2  # expect at least two strains in a mixed sample
        @test true  # placeholder
    end
end

# VALIDATION & QUALITY ASSESSMENT Tests
@testset "Assembly Validation" begin
    @testset "Reference-Free Validation" begin
        # merqury
        # ALE
        # CGAL
        # read mapping stats
        # Example: run QUAST analysis and ensure basic quality metric
        @test true  # placeholder
    end
    @testset "Reference-Based Validation" begin
        # Example: run MetaQUAST analysis and ensure basic quality metrics.
        # stats = MyceliaAssembly.validate_with_metaquast("final_assembly.fasta", reference="ref_genome.fasta")
        # @test stats.NGA50 > 10000
        @test true  # placeholder
    end
    @testset "Marker Gene Completeness" begin
        # Example: check completeness with CheckM.
        # quality = MyceliaAssembly.check_assembly_quality("final_assembly.fasta")
        # @test quality.completeness >= 90
        # @test quality.contamination <= 5
        @test true  # placeholder
    end
end

# Assembly & consensus generation: for combining reads under diverse coverage patterns.
# Assembly & Consensus Generation
# - Combining reads under diverse coverage patterns
# - Generating consensus sequences
# - Polishing assemblies (e.g., polish_fastq)
@testset "probabilistic ensemble assembly" begin
    @testset "Illumina" begin
        @test 1 + 1 == 2
    end

    @testset "Ultima" begin
        @test 1 + 1 == 2
    end

    @testset "Nanopore" begin
        @test 1 + 1 == 2
    end

    @testset "PacBio" begin
        @test 1 + 1 == 2
    end

    @testset "multi-entity, even coverage" begin
        @test 1 + 1 == 2
    end

    @testset "multi-entity, log-distributed coverage" begin
        @test 1 + 1 == 2
    end

    @testset "multi-platform" begin
        @test 1 + 1 == 2
    end
end

# Multi-omics alignment & mapping: for consistent coordinate systems or annotation references.
# Multi-omics Alignment & Mapping
# - Consistent coordinate systems or annotation references
# - Mapping reads to reference genomes (e.g., minimap2)
# - Handling different sequencing technologies (e.g., PacBio, ONT)
@testset "Multi-omics alignment & mapping" begin
end

# Graph-based integration: for building and merging genomic/transcriptomic/proteomic graphs.
# Graph-based Integration
# - Building and merging genomic/transcriptomic/proteomic graphs
# - Graph traversal and pathfinding (e.g., take_a_walk)
# - Graph-based clustering and visualization (e.g., draw_radial_tree)
@testset "pangenome construction" begin
    @testset "?" begin
        @test 1 + 1 == 2
    end
end

@testset "pantranscriptome construction" begin
    @testset "?" begin
        @test 1 + 1 == 2
    end
end

@testset "panproteome construction" begin
    @testset "?" begin
        @test 1 + 1 == 2
    end
end

# Annotation & feature extraction: for identifying and comparing variants, functional domains, or transcripts.
# Annotation & Feature Extraction
# - Identifying and comparing variants, functional domains, or transcripts
# - Annotating sequences (e.g., annotate_fasta)
# - Parsing and handling GFF files (e.g., parse_xam)
@testset "Annotation & feature extraction" begin
    @testset "?" begin
        @test 1 + 1 == 2
    end
end

# Comparative analyses: for visualizing similarities/differences across strains or conditions.
# Comparative Analyses
# - Visualizing similarities/differences across strains or conditions
# - Comparative genomics and transcriptomics
# - Clustering and dimensionality reduction (e.g., heirarchically_cluster_distance_matrix)
@testset "Comparative Analyses" begin
    @testset "?" begin
        @test 1 + 1 == 2
    end
end

# Sequence Classification
# - Classifying sequences based on k-mer content
# - Using machine learning models for classification
# - Handling large-scale sequence data efficiently
@testset "sequence classification" begin
    @testset "?" begin
        @test 1 + 1 == 2
    end
end

@testset "phylogenetic analyses" begin
    @testset "?" begin
        @test 1 + 1 == 2
    end
end

# Advanced metadata handling: for storing and retrieving detailed lineage or ontology data.
# Advanced Metadata Handling
# - Storing and retrieving detailed lineage or ontology data
# - Handling taxonomic data (e.g., list_full_taxonomy)
# - Integrating external databases and resources (e.g., NCBI, UniProt)

# Multi-omics Integration
# - Integrating data from genomics, transcriptomics, proteomics, and metabolomics
# - Cross-omics correlation and analysis
# - Building comprehensive multi-omics models

# Visualization & Reporting
# - Generating visual reports for analysis results
# - Interactive visualization of multi-omics data
# - Exporting results in various formats (e.g., PNG, SVG)
