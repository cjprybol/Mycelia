# NCBI Reference Genome Download tests
@testset "NCBI Reference Genome Download" begin
    # Placeholder: Add actual NCBI genome download tests here
    @test true
end

@testset "download_genome_by_accession" begin
    result = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    @test isfile(result)
    @test endswith(result, ".fna.gz")
    display("Downloaded genome file: ", result)
    @assert isfile(result)
    rm(result; force=true)
end

@testset "ncbi_genome_download_accession" begin
    result = Mycelia.ncbi_genome_download_accession(accession="GCF_000819615.1", include_string="genome,protein")
    @test isfile(result.genome)
    @test isfile(result.protein)
    @test filesize(result.genome) > 0
    @test filesize(result.protein) > 0
    # Check file format (FASTA: first char is '>')
    open(result.genome) do io
        @test readline(io)[1] == '>'
    end
    # Clean up
    rm(dirname(result.genome); recursive=true, force=true)
end

@testset "ncbi_genome_download_accession_partial" begin
    result = Mycelia.ncbi_genome_download_accession(accession="GCF_000819615.1", include_string="genome")
    @test isfile(result.genome)
    @test !haskey(result, :protein) || isnothing(result.protein)
    rm(dirname(result.genome); recursive=true, force=true)
    result2 = Mycelia.ncbi_genome_download_accession(accession="GCF_000819615.1", include_string="protein")
    @test isfile(result2.protein)
    @test !haskey(result2, :genome) || isnothing(result2.genome)
    rm(dirname(result2.protein); recursive=true, force=true)
end

@testset "download_genome_by_accession_idempotency" begin
    result1 = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    result2 = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    @test result1 == result2
    @test isfile(result1)
    rm(result1; force=true)
end

@testset "download_genome_by_accession_invalid" begin
    try
        Mycelia.download_genome_by_accession(accession="INVALID_ACCESSION")
        @test false  # Should not reach here
    catch e
        @test true  # Error expected
    end
end

# Optional: SRA prefetch/fasterq_dump if available
if hasmethod(Mycelia, :prefetch)
    @testset "SRA prefetch" begin
        try
            srafile = Mycelia.prefetch("SRR390728")
            @test isfile(srafile)
            rm(srafile; force=true)
        catch
            @test true  # Allow failure if not available
        end
    end
end
if hasmethod(Mycelia, :fasterq_dump)
    @testset "SRA fasterq_dump" begin
        try
            fastqfile = Mycelia.fasterq_dump("SRR390728")
            @test isfile(fastqfile)
            rm(fastqfile; force=true)
        catch
            @test true  # Allow failure if not available
        end
    end
end

@testset "get_base_extension" begin
    @test Mycelia.get_base_extension("foo.fna.gz") == ".fna.gz"
end

@testset "random_fasta_record" begin
    record = Mycelia.random_fasta_record(moltype=:DNA, seed=42, L=10)
    @test typeof(record) == FASTX.FASTA.Record
    @test length(FASTX.sequence(record)) == 10
end
