# FASTA simulation and acquisition tests
const SEED = 42
const phiX174_accession_id = "NC_001422.1"
const phiX174_assembly_id = "GCF_000819615.1"

@testset "FASTA simulation and acquisition" begin
    @testset "random_fasta_record" begin
        # Test correct sequence length and alphabet for DNA, RNA, AA
        # Example: dna_record = Mycelia.random_fasta_record(moltype=:DNA, seed=42, L=10)
        # @test typeof(dna_record) == FASTX.FASTA.Record
        # @test length(FASTX.sequence(dna_record)) == 10
    end
    @testset "download_genome_by_accession" begin
        # Test file download and format
        # Example: result = Mycelia.download_genome_by_accession(accession="NC_001422.1")
        # @test isfile(result)
        # @test result endswith ".fna.gz"
    end
    @testset "ncbi_genome_download_accession" begin
        # Test all expected files are present
        # Example: result = Mycelia.ncbi_genome_download_accession(accession="GCF_000819615.1", include_string="genome,protein")
        # @test isfile(result.genome)
        # @test isfile(result.protein)
    end
    @testset "get_base_extension" begin
        # @test Mycelia.get_base_extension("foo.fna.gz") == ".fna.gz"
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
