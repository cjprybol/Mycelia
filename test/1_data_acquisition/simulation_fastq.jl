# FASTQ simulation tests
const phiX174_assembly_id = "GCF_000819615.1"

@testset "FASTQ simulation" begin
    @testset "Illumina" begin
        if isdir(phiX174_assembly_id)
            rm(phiX174_assembly_id, recursive=true)
        end
        phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
        read_simulation_result = Mycelia.simulate_illumina_paired_reads(in_fasta = phiX174_assembly_dataset.genome, coverage=10)
        @test isfile(read_simulation_result.forward_reads)
        @test isfile(read_simulation_result.reverse_reads)
        @test isfile(read_simulation_result.sam)
        @test isfile(read_simulation_result.error_free_sam)
        rm(phiX174_assembly_id, recursive=true)
    end
    @testset "Ultima" begin
        @test 1 + 1 == 2
    end
    @testset "Nanopore" begin
        @test 1 + 1 == 2
    end
    @testset "PacBio" begin
        @test 1 + 1 == 2
    end
    @testset "multi-entity, even coverage" begin
        @test 1 + 1 == 2
    end
    @testset "multi-entity, log-distributed coverage" begin
        @test 1 + 1 == 2
    end
end
