@testset "download_genome_by_accession" begin
    result = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    @test isfile(result)
    @test endswith(result, ".fna.gz")
    @assert isfile(result)
    rm(result; force=true)
end

@testset "ncbi_genome_download_accession" begin
    accession = "GCF_000819615.1"
    result = Mycelia.ncbi_genome_download_accession(accession=accession)
    @test isdir(result.directory)
    @test isfile(result.genome) && filesize(result.genome) > 0
    @test isfile(result.protein) && filesize(result.protein) > 0
    @test isfile(result.cds) && filesize(result.cds) > 0
    @test isfile(result.gff3) && filesize(result.gff3) > 0
    @test isfile(result.seqreport) && filesize(result.seqreport) > 0
    # Check file format (FASTA: first char is '>')
    open(result.genome) do io
        @test readline(io)[1] == '>'
    end
    # Clean up
    @assert isdir(accession)
    rm(accession, recursive=true, force=true)
end

@testset "ncbi_genome_download_accession_partial" begin
    accession = "GCF_000819615.1"
    if isdir(accession)
        rm(accession, recursive=true, force=true)
    end
    result = Mycelia.ncbi_genome_download_accession(accession=accession, include_string="genome")
    @test isfile(result.genome)
    @test !haskey(result, :protein) || isnothing(result.protein)
    rm(accession, recursive=true, force=true)
    result2 = Mycelia.ncbi_genome_download_accession(accession=accession, include_string="protein")
    @test isfile(result2.protein)
    @test !haskey(result2, :genome) || isnothing(result2.genome)
    rm(accession, recursive=true, force=true)
end

@testset "download_genome_by_accession_idempotency" begin
    result1 = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    result2 = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    @test result1 == result2
    @test isfile(result1)
    rm(result1; force=true)
end

@testset "SRA prefetch and fasterq_dump stepwise" begin
    sra_public_fastq_metadata = CSV.read(CodecZlib.GzipDecompressorStream(open("metadata/20250702.sra-public-fastqs.csv.gz")), DataFrames.DataFrame)
    SRR_identifier = sort!(sra_public_fastq_metadata, "Bytes")[1, "Run"]
    prefetch_result = Mycelia.prefetch(SRR=SRR_identifier)
    @test basename(prefetch_result.directory) == SRR_identifier
    @test basename(prefetch_result.archive) == SRR_identifier * ".sra"
    @test isdir(prefetch_result.directory)
    @test isfile(prefetch_result.archive)
    fasterq_result = Mycelia.fasterq_dump(srr_identifier=SRR_identifier)
    @test ismissing(fasterq_result.forward_reads)
    @test ismissing(fasterq_result.reverse_reads)
    @test basename(fasterq_result.unpaired_reads) == "$(SRR_identifier).fastq.gz"
    @test isfile(fasterq_result.unpaired_reads) && filesize(fasterq_result.unpaired_reads) > 0
    rm(SRR_identifier, recursive=true)
end

@testset "SRA fasterq_dump standalone" begin
    sra_public_fastq_metadata = CSV.read(CodecZlib.GzipDecompressorStream(open("metadata/20250702.sra-public-fastqs.csv.gz")), DataFrames.DataFrame)
    SRR_identifier = sort!(sra_public_fastq_metadata, "Bytes")[1, "Run"]
    fasterq_result = Mycelia.fasterq_dump(srr_identifier=SRR_identifier)
    @test ismissing(fasterq_result.forward_reads)
    @test ismissing(fasterq_result.reverse_reads)
    @test basename(fasterq_result.unpaired_reads) == "$(SRR_identifier).fastq.gz"
    @test isfile(fasterq_result.unpaired_reads) && filesize(fasterq_result.unpaired_reads) > 0
    rm(SRR_identifier, recursive=true)
end

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
