# SRA downloading tests
@testset "SRA downloading" begin
    outdir = mkpath("test-sra-tools")
    @testset "download SRA short read data" begin
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
