# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/1_data_acquisition/ncbi_download.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/1_data_acquisition/ncbi_download.jl", "test/1_data_acquisition", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import CSV
import CodecZlib
import DataFrames

Test.@testset "datasets CLI version" begin
    version = Mycelia.datasets_cli_version()
    Test.@test !isempty(strip(version))
end

Test.@testset "datasets genome summary (E. coli)" begin
    summary = Mycelia.datasets_genome_summary(taxon="562", limit=1)
    Test.@test isa(summary, DataFrames.DataFrame)
    Test.@test DataFrames.nrow(summary) > 0
end

Test.@testset "datasets dehydrated download + rehydrate" begin
    accession = "GCF_000819615.1"
    workdir = mktempdir()
    try
        result = Mycelia.datasets_download_genome(accession;
            input_type="accession",
            include=["genome"],
            dehydrated=true,
            extract=true,
            outdir=workdir
        )
        Test.@test isfile(result.zip_path)
        Test.@test isdir(result.directory)
        pre_files = Mycelia.recursively_list_files(result.directory)
        Test.@test !any(x -> occursin(r"\.fna(\.gz)?$", x), pre_files)
        Mycelia.datasets_rehydrate(result.directory; gzip=true, max_workers=2)
        post_files = Mycelia.recursively_list_files(result.directory)
        Test.@test any(x -> occursin(r"\.fna(\.gz)?$", x), post_files)
    finally
        rm(workdir, recursive=true, force=true)
    end
end

Test.@testset "datasets gene download (BRCA1)" begin
    workdir = mktempdir()
    try
        result = Mycelia.datasets_download_gene("BRCA1";
            input_type="symbol",
            taxon="human",
            include=["gene"],
            extract=true,
            outdir=workdir
        )
        Test.@test isfile(result.zip_path)
        Test.@test isdir(result.directory)
        files = Mycelia.recursively_list_files(result.directory)
        Test.@test any(x -> occursin(r"gene\.fna(\.gz)?$", basename(x)), files)
    finally
        rm(workdir, recursive=true, force=true)
    end
end

Test.@testset "download_genome_by_accession" begin
    ## const phiX174_accession_id = "NC_001422.1"
    workdir = mktempdir()
    try
        result = Mycelia.download_genome_by_accession(accession="NC_001422", outdir=workdir)
        Test.@test isfile(result)
        Test.@test endswith(result, ".fna.gz")
        @assert isfile(result)
    finally
        rm(workdir, recursive=true, force=true)
    end
end

Test.@testset "ncbi_genome_download_accession" begin
    ## const phiX174_assembly_id = "GCF_000819615.1"
    ## need to use the versioned identifier - unversioned will add it in breaking the search path
    accession = "GCF_000819615.1"
    workdir = mktempdir()
    try
        result = Mycelia.ncbi_genome_download_accession(accession=accession, outdir=workdir)
        Test.@test isdir(result.directory)
        Test.@test isfile(result.genome) && filesize(result.genome) > 0
        Test.@test isfile(result.protein) && filesize(result.protein) > 0
        Test.@test isfile(result.cds) && filesize(result.cds) > 0
        Test.@test isfile(result.gff3) && filesize(result.gff3) > 0
        Test.@test isfile(result.seqreport) && filesize(result.seqreport) > 0
        ## Check file format (FASTA: first char is '>')
        open(result.genome) do io
            Test.@test readline(io)[1] == '>'
        end
        Test.@test isdir(joinpath(workdir, accession))
    finally
        rm(workdir, recursive=true, force=true)
    end
end

Test.@testset "ncbi_genome_download_accession_partial" begin
    accession = "GCF_000819615.1"
    workdir = mktempdir()
    try
        result = Mycelia.ncbi_genome_download_accession(accession=accession, include_string="genome", outdir=workdir)
        Test.@test isfile(result.genome)
        Test.@test !haskey(result, :protein) || isnothing(result.protein)
        rm(joinpath(workdir, accession), recursive=true, force=true)
        result2 = Mycelia.ncbi_genome_download_accession(accession=accession, include_string="protein", outdir=workdir)
        Test.@test isfile(result2.protein)
        Test.@test !haskey(result2, :genome) || isnothing(result2.genome)
    finally
        rm(workdir, recursive=true, force=true)
    end
end

Test.@testset "download_genome_by_accession_idempotency" begin
    workdir = mktempdir()
    try
        result1 = Mycelia.download_genome_by_accession(accession="NC_001422.1", outdir=workdir)
        result2 = Mycelia.download_genome_by_accession(accession="NC_001422.1", outdir=workdir)
        Test.@test result1 == result2
        Test.@test isfile(result1)
    finally
        rm(workdir, recursive=true, force=true)
    end
end

Test.@testset "SRA prefetch and fasterq_dump stepwise" begin
    sra_public_fastq_metadata_file = joinpath(@__DIR__, "..", "metadata", "20250702.sra-public-fastqs.csv.gz")
    if !isfile(sra_public_fastq_metadata_file)
        @info "Skipping SRA prefetch test; metadata missing at $(sra_public_fastq_metadata_file)"
        return
    end
    sra_public_fastq_metadata = CSV.read(CodecZlib.GzipDecompressorStream(open(sra_public_fastq_metadata_file)), DataFrames.DataFrame)
    SRR_identifier = sort!(sra_public_fastq_metadata, "Bytes")[1, "Run"]
    outdir = mktempdir()
    try
        prefetch_result = Mycelia.prefetch(SRR=SRR_identifier, outdir=outdir)
        Test.@test basename(prefetch_result.directory) == SRR_identifier
        Test.@test basename(prefetch_result.archive) == SRR_identifier * ".sra"
        Test.@test isdir(prefetch_result.directory)
        Test.@test isfile(prefetch_result.archive)
        fasterq_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=SRR_identifier)
        Test.@test ismissing(fasterq_result.forward_reads)
        Test.@test ismissing(fasterq_result.reverse_reads)
        Test.@test basename(fasterq_result.unpaired_reads) == "$(SRR_identifier).fastq.gz"
        Test.@test isfile(fasterq_result.unpaired_reads) && filesize(fasterq_result.unpaired_reads) > 0
    finally
        rm(outdir, recursive=true, force=true)
    end
end

Test.@testset "SRA fasterq_dump standalone" begin
    sra_public_fastq_metadata_file = joinpath(@__DIR__, "..", "metadata", "20250702.sra-public-fastqs.csv.gz")
    if !isfile(sra_public_fastq_metadata_file)
        @info "Skipping SRA fasterq_dump test; metadata missing at $(sra_public_fastq_metadata_file)"
        return
    end
    sra_public_fastq_metadata = CSV.read(CodecZlib.GzipDecompressorStream(open(sra_public_fastq_metadata_file)), DataFrames.DataFrame)
    SRR_identifier = sort!(sra_public_fastq_metadata, "Bytes")[1, "Run"]
    outdir = mktempdir()
    try
        fasterq_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=SRR_identifier)
        Test.@test ismissing(fasterq_result.forward_reads)
        Test.@test ismissing(fasterq_result.reverse_reads)
        Test.@test basename(fasterq_result.unpaired_reads) == "$(SRR_identifier).fastq.gz"
        Test.@test isfile(fasterq_result.unpaired_reads) && filesize(fasterq_result.unpaired_reads) > 0
    finally
        rm(outdir, recursive=true, force=true)
    end
end

# SRA downloading tests
Test.@testset "SRA downloading" begin
    outdir = mktempdir()
    try
        Test.@testset "download SRA short read data" begin
            srr_identifier = "SRR31271002"
            prefetch_results = Mycelia.prefetch(SRR=srr_identifier, outdir=outdir)
            Test.@test prefetch_results.directory == joinpath(outdir, srr_identifier)
            Test.@test prefetch_results.archive == joinpath(outdir, srr_identifier, "$(srr_identifier).sra")
            Test.@test filesize(prefetch_results.archive) == 18295489
            fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
            Test.@test fasterq_dump_result.forward_reads == joinpath(outdir, srr_identifier, "$(srr_identifier)_1.fastq.gz")
            Test.@test fasterq_dump_result.reverse_reads == joinpath(outdir, srr_identifier, "$(srr_identifier)_2.fastq.gz")
            Test.@test ismissing(fasterq_dump_result.unpaired_reads)
            Test.@test filesize(fasterq_dump_result.forward_reads) == 11637534
            Test.@test filesize(fasterq_dump_result.reverse_reads) == 11908008
        end
        Test.@testset "download SRA long read data" begin
            srr_identifier = "SRR31812976"
            prefetch_results = Mycelia.prefetch(SRR=srr_identifier, outdir=outdir)
            Test.@test prefetch_results.directory == joinpath(outdir, srr_identifier)
            Test.@test prefetch_results.archive == joinpath(outdir, srr_identifier, "$(srr_identifier).sra")
            Test.@test filesize(prefetch_results.archive) == 33392268
            fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
            Test.@test ismissing(fasterq_dump_result.forward_reads)
            Test.@test ismissing(fasterq_dump_result.reverse_reads)
            Test.@test fasterq_dump_result.unpaired_reads == joinpath(outdir, srr_identifier, "$(srr_identifier).fastq.gz")
            Test.@test filesize(fasterq_dump_result.unpaired_reads) == 34353450
        end
    finally
        rm(outdir, recursive=true, force=true)
    end
end
