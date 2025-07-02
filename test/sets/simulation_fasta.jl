# FASTA simulation and acquisition tests
const SEED = 42
const phiX174_accession_id = "NC_001422.1"
const phiX174_assembly_id = "GCF_000819615.1"

@testset "FASTA simulation and acquisition" begin
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
