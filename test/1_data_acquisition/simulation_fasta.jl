# FASTA simulation and acquisition tests

import Pkg
Pkg.activate("..")
using Test
import Mycelia
import FASTX
import Random
import CodecZlib

const SEED = 42

@testset "FASTA simulation and acquisition" begin
    @testset "get_base_extension" begin
        @test Mycelia.get_base_extension("foo.fasta.gz") == ".fasta"
        for ext in ["fna", "faa", "frn"]
            @test Mycelia.get_base_extension("foo.$ext.gz") == "." * ext
        end

        @test Mycelia.get_base_extension("foo.fasta") == ".fasta"
        for ext in ["fna", "faa", "frn"]
            @test Mycelia.get_base_extension("foo.$ext") == "." * ext
        end
    end
    
    @testset "random_fasta_record" begin
        for molecule in [:DNA, :RNA, :AA]
            a = Mycelia.random_fasta_record(moltype=molecule, seed=42, L=10)
            b = Mycelia.random_fasta_record(moltype=molecule, seed=42, L=10)
            @test typeof(a) == typeof(b) == FASTX.FASTA.Record
            @test length(FASTX.sequence(a)) == 10
            @test FASTX.sequence(a) == FASTX.sequence(b)
            @test FASTX.identifier(a) != FASTX.identifier(b)
            if molecule == :DNA
                @test FASTX.sequence(a) == "CCGCCGCTCA"
            elseif molecule == :RNA
                @test FASTX.sequence(a) == "CCGCCGCUCA"
            elseif molecule == :AA
                @test FASTX.sequence(a) == "VATAGWWITI"
            end
        end
    end
    @testset "dna record" begin
        dna_record = Mycelia.random_fasta_record(moltype=:DNA, seed=SEED, L = 10)
        @test FASTX.sequence(dna_record) == "CCGCCGCTCA"
    end
    @testset "rna record" begin
        rna_record = Mycelia.random_fasta_record(moltype=:RNA, seed=SEED, L = 10)
        @test FASTX.sequence(rna_record) == "CCGCCGCUCA"
    end
    @testset "aa record" begin
        aa_record = Mycelia.random_fasta_record(moltype=:AA, seed=SEED, L = 10)
        @test FASTX.sequence(aa_record) == "VATAGWWITI"
    end
    @testset "deterministic sequence" begin
        rec1 = Mycelia.random_fasta_record(seed=SEED, L=15)
        rec2 = Mycelia.random_fasta_record(seed=SEED, L=15)
        @test FASTX.sequence(rec1) == FASTX.sequence(rec2)
        @test FASTX.identifier(rec1) != FASTX.identifier(rec2)
    end
    @testset "write_fasta" begin
        record = Mycelia.random_fasta_record(seed=SEED, L=10)
        fasta_file = "temp_write_fasta.fna.gz"
        try
            result = Mycelia.write_fasta(outfile=fasta_file, records=[record])
            @test result == fasta_file
            @test isfile(fasta_file)
            FASTX.FASTA.Reader(CodecZlib.GzipDecompressorStream(open(fasta_file))) do reader
                rec = first(reader)
                @test FASTX.sequence(rec) == FASTX.sequence(record)
                @test FASTX.identifier(rec) == FASTX.identifier(record)
            end
        finally
            isfile(fasta_file) && rm(fasta_file)
        end
    end
    @testset "virus phiX174" begin
        genome_result = Mycelia.download_genome_by_accession(accession=phiX174_accession_id)
        @test basename(genome_result) == phiX174_accession_id * ".fna.gz"
        @test Mycelia.get_base_extension(genome_result) == ".fna.gz"
        rm(genome_result)
        phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="gff3,rna,cds,protein,genome,seq-report")
        @test basename(phiX174_assembly_dataset.genome) == phiX174_assembly_id * "_ViralProj14015_genomic.fna"
        @test Mycelia.get_base_extension(phiX174_assembly_dataset.genome) == ".fna"
        @test Mycelia.get_base_extension(phiX174_assembly_dataset.protein) == ".faa"
        rm(phiX174_assembly_id, recursive=true)
    end
end
