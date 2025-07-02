# NCBI Reference Genome Download tests
@testset "NCBI Reference Genome Download" begin
    # Placeholder: Add actual NCBI genome download tests here
    @test true
end

@testset "download_genome_by_accession" begin
    # Test downloading a known genome by accession
    # Example: result = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    # @test isfile(result)
    # @test result endswith ".fna.gz"
    # Clean up downloaded file(s)
end

@testset "ncbi_genome_download_accession" begin
    # Test downloading a genome and associated files by assembly accession
    # Example: result = Mycelia.ncbi_genome_download_accession(accession="GCF_000819615.1", include_string="genome,protein")
    # @test isfile(result.genome)
    # @test isfile(result.protein)
    # Clean up downloaded files/directories
end

@testset "get_base_extension" begin
    # @test Mycelia.get_base_extension("foo.fna.gz") == ".fna.gz"
end

@testset "random_fasta_record" begin
    # Example: record = Mycelia.random_fasta_record(moltype=:DNA, seed=42, L=10)
    # @test typeof(record) == FASTX.FASTA.Record
    # @test length(FASTX.sequence(record)) == 10
end
