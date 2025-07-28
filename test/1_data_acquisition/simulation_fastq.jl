# FASTQ simulation tests

# import Pkg
# if isinteractive()
#     Pkg.activate("..")
# end
import Test
import Mycelia
import FASTX
import Random

const phiX174_assembly_id = "GCF_000819615.1"

Test.@testset "FASTQ simulation" begin
    Test.@testset "Illumina" begin
        if isdir(phiX174_assembly_id)
            rm(phiX174_assembly_id, recursive=true)
        end
        phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
        read_simulation_result = Mycelia.simulate_illumina_paired_reads(in_fasta = phiX174_assembly_dataset.genome, coverage=10)
        Test.@test isfile(read_simulation_result.forward_reads)
        Test.@test isfile(read_simulation_result.reverse_reads)
        Test.@test isfile(read_simulation_result.sam)
        Test.@test isfile(read_simulation_result.error_free_sam)
        rm(phiX174_assembly_id, recursive=true)
    end
    Test.@testset "Ultima" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "Nanopore" begin
        if isdir(phiX174_assembly_id)
            rm(phiX174_assembly_id, recursive=true)
        end
        phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
        read_simulation_result = Mycelia.simulate_pacbio_reads(in_fasta = phiX174_assembly_dataset.genome, coverage=10)
    end
    Test.@testset "PacBio" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "multi-entity, even coverage" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "multi-entity, log-distributed coverage" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "simulate_pacbio_reads" begin
        # Test output file exists and is gzipped
        # Example: result = Mycelia.simulate_pacbio_reads(fasta="test.fna", quantity="10x")
        # Test.@test isfile(result)
        # Test.@test endswith(result, ".fq.gz")
    end
    Test.@testset "simulate_nanopore_reads" begin
        # Test output file exists and is gzipped
        # Example: result = Mycelia.simulate_nanopore_reads(fasta="test.fna", quantity="10x")
        # Test.@test isfile(result)
        # Test.@test endswith(result, ".fq.gz")
    end
end
