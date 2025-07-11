import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia

Test.@testset "download_genome_by_accession" begin
    # const phiX174_accession_id = "NC_001422.1"
    result = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    Test.@test isfile(result)
    Test.@test endswith(result, ".fna.gz")
    @assert isfile(result)
    rm(result; force=true)
end

Test.@testset "ncbi_genome_download_accession" begin
    # const phiX174_assembly_id = "GCF_000819615.1"
    accession = "GCF_000819615.1"
    result = Mycelia.ncbi_genome_download_accession(accession=accession)
    Test.@test isdir(result.directory)
    Test.@test isfile(result.genome) && filesize(result.genome) > 0
    Test.@test isfile(result.protein) && filesize(result.protein) > 0
    Test.@test isfile(result.cds) && filesize(result.cds) > 0
    Test.@test isfile(result.gff3) && filesize(result.gff3) > 0
    Test.@test isfile(result.seqreport) && filesize(result.seqreport) > 0
    # Check file format (FASTA: first char is '>')
    open(result.genome) do io
        Test.@test readline(io)[1] == '>'
    end
    # Clean up
    @assert isdir(accession)
    rm(accession, recursive=true, force=true)
end

Test.@testset "ncbi_genome_download_accession_partial" begin
    accession = "GCF_000819615.1"
    if isdir(accession)
        rm(accession, recursive=true, force=true)
    end
    result = Mycelia.ncbi_genome_download_accession(accession=accession, include_string="genome")
    Test.@test isfile(result.genome)
    Test.@test !haskey(result, :protein) || isnothing(result.protein)
    rm(accession, recursive=true, force=true)
    result2 = Mycelia.ncbi_genome_download_accession(accession=accession, include_string="protein")
    Test.@test isfile(result2.protein)
    Test.@test !haskey(result2, :genome) || isnothing(result2.genome)
    rm(accession, recursive=true, force=true)
end

Test.@testset "download_genome_by_accession_idempotency" begin
    result1 = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    result2 = Mycelia.download_genome_by_accession(accession="NC_001422.1")
    Test.@test result1 == result2
    Test.@test isfile(result1)
    rm(result1; force=true)
end

Test.@testset "SRA prefetch and fasterq_dump stepwise" begin
    sra_public_fastq_metadata = CSV.read(CodecZlib.GzipDecompressorStream(open("metadata/20250702.sra-public-fastqs.csv.gz")), DataFrames.DataFrame)
    SRR_identifier = sort!(sra_public_fastq_metadata, "Bytes")[1, "Run"]
    prefetch_result = Mycelia.prefetch(SRR=SRR_identifier)
    Test.@test basename(prefetch_result.directory) == SRR_identifier
    Test.@test basename(prefetch_result.archive) == SRR_identifier * ".sra"
    Test.@test isdir(prefetch_result.directory)
    Test.@test isfile(prefetch_result.archive)
    fasterq_result = Mycelia.fasterq_dump(srr_identifier=SRR_identifier)
    Test.@test ismissing(fasterq_result.forward_reads)
    Test.@test ismissing(fasterq_result.reverse_reads)
    Test.@test basename(fasterq_result.unpaired_reads) == "$(SRR_identifier).fastq.gz"
    Test.@test isfile(fasterq_result.unpaired_reads) && filesize(fasterq_result.unpaired_reads) > 0
    rm(SRR_identifier, recursive=true)
end

Test.@testset "SRA fasterq_dump standalone" begin
    sra_public_fastq_metadata = CSV.read(CodecZlib.GzipDecompressorStream(open("metadata/20250702.sra-public-fastqs.csv.gz")), DataFrames.DataFrame)
    SRR_identifier = sort!(sra_public_fastq_metadata, "Bytes")[1, "Run"]
    fasterq_result = Mycelia.fasterq_dump(srr_identifier=SRR_identifier)
    Test.@test ismissing(fasterq_result.forward_reads)
    Test.@test ismissing(fasterq_result.reverse_reads)
    Test.@test basename(fasterq_result.unpaired_reads) == "$(SRR_identifier).fastq.gz"
    Test.@test isfile(fasterq_result.unpaired_reads) && filesize(fasterq_result.unpaired_reads) > 0
    rm(SRR_identifier, recursive=true)
end

# SRA downloading tests
Test.@testset "SRA downloading" begin
    outdir = mkpath("test-sra-tools")
    Test.@testset "download SRA short read data" begin
        srr_identifier = "SRR31271002"
        prefetch_results = Mycelia.prefetch(SRR=srr_identifier, outdir=outdir)
        Test.@test prefetch_results.directory == "test-sra-tools/$(srr_identifier)"
        Test.@test prefetch_results.archive == "test-sra-tools/$(srr_identifier)/$(srr_identifier).sra"
        Test.@test filesize(prefetch_results.archive) == 18295489
        fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
        Test.@test fasterq_dump_result.forward_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier)_1.fastq.gz"
        Test.@test fasterq_dump_result.reverse_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier)_2.fastq.gz"
        Test.@test ismissing(fasterq_dump_result.unpaired_reads)
        Test.@test filesize(fasterq_dump_result.forward_reads) == 11637534
        Test.@test filesize(fasterq_dump_result.reverse_reads) == 11908008
    end
    Test.@testset "download SRA long read data" begin
        srr_identifier = "SRR31812976"
        prefetch_results = Mycelia.prefetch(SRR=srr_identifier, outdir=outdir)
        Test.@test prefetch_results.directory == "$(outdir)/$(srr_identifier)"
        Test.@test prefetch_results.archive == "$(outdir)/$(srr_identifier)/$(srr_identifier).sra"
        Test.@test filesize(prefetch_results.archive) == 33392268
        fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
        Test.@test ismissing(fasterq_dump_result.forward_reads)
        Test.@test ismissing(fasterq_dump_result.reverse_reads)
        Test.@test fasterq_dump_result.unpaired_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier).fastq.gz"
        Test.@test filesize(fasterq_dump_result.unpaired_reads) == 34353450
    end
    rm(outdir, recursive=true)
end
