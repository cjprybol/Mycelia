# FASTA simulation and acquisition tests

import Pkg
Pkg.activate("..")
using Test
import Mycelia
import FASTX
import Random

const SEED = 42

@testset "FASTA simulation and acquisition" begin
    @testset "get_base_extension" begin
        @test Mycelia.get_base_extension("foo.fasta.gz") == ".fasta"
        @test Mycelia.get_base_extension("foo.fna.gz") == ".fna"
        @test Mycelia.get_base_extension("foo.faa.gz") == ".faa"
        @test Mycelia.get_base_extension("foo.frn.gz") == ".frn"
    
        @test Mycelia.get_base_extension("foo.fasta") == ".fasta"
        @test Mycelia.get_base_extension("foo.fna") == ".fna"
        @test Mycelia.get_base_extension("foo.faa") == ".faa"
        @test Mycelia.get_base_extension("foo.frn") == ".frn"
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
