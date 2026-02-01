import Test
import Mycelia
import FASTX
import StableRNGs

Test.@testset "Minimap Paired-End Mapping" begin
    rng = StableRNGs.StableRNG(1234)
    ref_record = Mycelia.random_fasta_record(L = 1000, seed = rand(rng, 0:typemax(Int)))
    ref_fasta = tempname() * ".fasta"
    writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
    write(writer, ref_record)
    close(writer)

    fastq_files = Mycelia.simulate_illumina_reads(
        fasta = ref_fasta,
        read_count = 2,
        rndSeed = rand(rng, 0:typemax(Int))
    )
    mem_gb = Int(Sys.free_memory()) / 1e9
    threads = 1

    Test.@testset "Header Toggle (Indexless)" begin
        default_cmd = Mycelia.minimap_map_paired_end(
            fasta = ref_fasta,
            forward = fastq_files.forward_reads,
            reverse = fastq_files.reverse_reads,
            as_string = true,
            mapping_type = "sr",
            mem_gb = mem_gb,
            threads = threads
        ).cmd
        Test.@test occursin("--no-header", default_cmd)

        keep_cmd = Mycelia.minimap_map_paired_end(
            fasta = ref_fasta,
            forward = fastq_files.forward_reads,
            reverse = fastq_files.reverse_reads,
            as_string = true,
            mapping_type = "sr",
            mem_gb = mem_gb,
            threads = threads,
            keep_header = true
        ).cmd
        Test.@test !occursin("--no-header", keep_cmd)
    end

    Test.@testset "Command Execution" begin
        map_result = Mycelia.minimap_map_paired_end(
            fasta = ref_fasta,
            forward = fastq_files.forward_reads,
            reverse = fastq_files.reverse_reads,
            mapping_type = "sr",
            mem_gb = mem_gb,
            threads = threads
        )
        run(map_result.cmd)
        Test.@test isfile(map_result.outfile)
        rm(map_result.outfile, force = true)
    end

    index_result = Mycelia.minimap_index(
        fasta = ref_fasta,
        mapping_type = "sr",
        mem_gb = mem_gb,
        threads = threads
    )
    run(index_result.cmd)

    Test.@testset "Header Toggle (With Index)" begin
        default_cmd = Mycelia.minimap_map_paired_end_with_index(
            forward = fastq_files.forward_reads,
            reverse = fastq_files.reverse_reads,
            index_file = index_result.outfile,
            as_string = true,
            mem_gb = mem_gb,
            threads = threads
        ).cmd
        Test.@test occursin("--no-header", default_cmd)

        keep_cmd = Mycelia.minimap_map_paired_end_with_index(
            forward = fastq_files.forward_reads,
            reverse = fastq_files.reverse_reads,
            index_file = index_result.outfile,
            as_string = true,
            mem_gb = mem_gb,
            threads = threads,
            keep_header = true
        ).cmd
        Test.@test !occursin("--no-header", keep_cmd)
    end

    Test.@testset "Minimap Map With Index Header Toggle" begin
        default_cmd = Mycelia.minimap_map_with_index(
            mapping_type = "sr",
            fastq = fastq_files.forward_reads,
            index_file = index_result.outfile,
            as_string = true,
            mem_gb = mem_gb,
            threads = threads
        ).cmd
        Test.@test occursin("--no-header", default_cmd)

        keep_cmd = Mycelia.minimap_map_with_index(
            mapping_type = "sr",
            fastq = fastq_files.forward_reads,
            index_file = index_result.outfile,
            as_string = true,
            mem_gb = mem_gb,
            threads = threads,
            keep_header = true
        ).cmd
        Test.@test !occursin("--no-header", keep_cmd)
    end

    rm(index_result.outfile, force = true)
    rm(ref_fasta, force = true)
    rm(fastq_files.forward_reads, force = true)
    rm(fastq_files.reverse_reads, force = true)
    if !isnothing(fastq_files.sam)
        rm(fastq_files.sam, force = true)
    end
    if !isnothing(fastq_files.error_free_sam)
        rm(fastq_files.error_free_sam, force = true)
    end
end
