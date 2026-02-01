# Legacy third-party assembler tests: metagenome.
import Test
import Mycelia
import Graphs
import Random
import StableRNGs
import BioSequences
import FASTX

threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

Test.@testset "Short Read Metagenomic Assembly" begin
    mktempdir() do dir
        # Create small reference genomes for metagenomic short read assembly
        meta_ref_fasta = joinpath(dir, "metagenomic_ref.fasta")
        rng_meta1 = StableRNGs.StableRNG(456)
        rng_meta2 = StableRNGs.StableRNG(457)
        meta_genome1 = BioSequences.randdnaseq(rng_meta1, 4000)  # 4kb genome
        meta_genome2 = BioSequences.randdnaseq(rng_meta2, 3000)  # 3kb genome

        # Create FASTA records and write using Mycelia.write_fasta
        meta_fasta_record1 = FASTX.FASTA.Record("test_metagenomic_genome_1", meta_genome1)
        meta_fasta_record2 = FASTX.FASTA.Record("test_metagenomic_genome_2", meta_genome2)
        Mycelia.write_fasta(outfile = meta_ref_fasta, records = [
            meta_fasta_record1, meta_fasta_record2])

        # Simulate Illumina paired-end reads with coverage for metagenomic assembly
        meta_simulated_reads = Mycelia.simulate_illumina_reads(
            fasta = meta_ref_fasta, coverage = 10, rndSeed = 456, quiet = true)

        # Extract the paired read files from the result
        meta_fastq1 = meta_simulated_reads.forward_reads
        meta_fastq2 = meta_simulated_reads.reverse_reads

        Test.@testset "Short Read Metagenomic Assembly - MEGAHIT" begin
            # Test MEGAHIT - clean up any existing directory first
            megahit_outdir = joinpath(dir, "megahit_assembly")
            if isdir(megahit_outdir)
                rm(megahit_outdir, recursive = true)
            end
            try
                result = Mycelia.run_megahit(fastq1 = meta_fastq1, fastq2 = meta_fastq2,
                    outdir = megahit_outdir, k_list = "21,29", threads = threads)
                Test.@test result.outdir == megahit_outdir
                Test.@test result.contigs == joinpath(megahit_outdir, "final.contigs.fa")
                Test.@test isfile(result.contigs)
                # Clean up after test
                rm(megahit_outdir, recursive = true, force = true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     MEGAHIT assembly failed due to resource constraints.
                #     Current test: 5kb genome, 15x coverage, ~75kb total sequence data

                #     Required resources for MEGAHIT:
                #     - Memory: ~1-2GB RAM minimum
                #     - CPU: 1-4 cores recommended
                #     - Disk: ~100MB temporary space

                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "MEGAHIT test skipped - insufficient resources"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(megahit_outdir, recursive = true, force = true)
                rethrow(e)
            end
        end

        Test.@testset "Short Read Metagenomic Assembly - metaSPAdes" begin

            # Test metaSPAdes - clean up any existing directory first  
            spades_outdir = joinpath(dir, "metaspades_assembly")
            if isdir(spades_outdir)
                rm(spades_outdir, recursive = true)
            end
            try
                result = Mycelia.run_metaspades(
                    fastq1 = meta_fastq1, fastq2 = meta_fastq2, outdir = spades_outdir,
                    k_list = "21", threads = threads, only_assembler = true)
                Test.@test result.outdir == spades_outdir
                Test.@test result.contigs == joinpath(spades_outdir, "contigs.fasta")
                Test.@test result.scaffolds == joinpath(spades_outdir, "scaffolds.fasta")
                Test.@test isfile(result.contigs)
                Test.@test isfile(result.scaffolds)
                # Clean up after test
                rm(spades_outdir, recursive = true, force = true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     metaSPAdes assembly failed due to resource constraints.
                #     Current test: 5kb genome, 15x coverage, ~75kb total sequence data

                #     Required resources for metaSPAdes:
                #     - Memory: ~2-4GB RAM minimum
                #     - CPU: 1-4 cores recommended
                #     - Disk: ~200MB temporary space

                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "metaSPAdes test skipped - insufficient resources"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(spades_outdir, recursive = true, force = true)
                rethrow(e)
            end
        end

        # Test.@testset "Short Read Metagenomic Assembly - Ray Meta" begin
        #     # Test Ray Meta - clean up any existing directory first  
        #     ray_meta_outdir = joinpath(dir, "ray_meta_assembly")
        #     if isdir(ray_meta_outdir)
        #         rm(ray_meta_outdir, recursive=true)
        #     end
        #     try
        #         result = Mycelia.run_ray_meta([meta_fastq1, meta_fastq2], outdir=ray_meta_outdir, k=25)
        #         Test.@test result.outdir == ray_meta_outdir
        #         Test.@test result.contigs == joinpath(ray_meta_outdir, "Contigs.fasta")
        #         Test.@test result.scaffolds == joinpath(ray_meta_outdir, "Scaffolds.fasta")
        #         Test.@test isfile(result.contigs)
        #         Test.@test isfile(result.scaffolds)
        #         # Clean up after test
        #         rm(ray_meta_outdir, recursive=true, force=true)
        #     catch e
        #         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
        #             @warn """
        #             Ray Meta assembly failed due to resource constraints.
        #             Current test: 4kb genome, 12x coverage, ~48kb total sequence data

        #             Required resources for Ray Meta:
        #             - Memory: ~2-4GB RAM minimum
        #             - CPU: 1-8 cores recommended (distributed metagenomic assembler)
        #             - Disk: ~500MB temporary space
        #             - Note: Ray Meta uses community detection for metagenomic samples

        #             To fix: Increase available memory or reduce test genome size further.
        #             """
        #             Test.@test_skip "Ray Meta test skipped - insufficient resources"
        #         else
        #             rethrow(e)
        #         end
        #         # Clean up on failure
        #         rm(ray_meta_outdir, recursive=true, force=true)
        #     end
        # end

        # Test.@testset "Short Read Metagenomic Assembly - Meta-IDBA" begin
        # end

        # not working
        # Test.@testset "Short Read Metagenomic Assembly - MetaVelvet" begin
        #     # Test metavelvet - clean up any existing directory first  
        #     metavelvet_outdir = joinpath(dir, "metavelvet_assembly")
        #     if isdir(metavelvet_outdir)
        #         rm(metavelvet_outdir, recursive=true)
        #     end
        #     try
        #         result = Mycelia.run_metavelvet(meta_fastq1, fastq2=meta_fastq2, outdir=metavelvet_outdir, k=25)
        #         Test.@test result.outdir == metavelvet_outdir
        #         Test.@test result.contigs == joinpath(metavelvet_outdir, "meta-velvetg.contigs.fa")
        #         Test.@test isfile(result.contigs)
        #         # Clean up after test
        #         rm(metavelvet_outdir, recursive=true, force=true)
        #     catch e
        #         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
        #             @warn """
        #             metavelvet assembly failed due to resource constraints.
        #             Current test: 4kb genome, 12x coverage, ~48kb total sequence data

        #             Required resources for metavelvet:
        #             - Memory: ~1-3GB RAM minimum
        #             - CPU: 1-4 cores recommended
        #             - Disk: ~200MB temporary space
        #             - Note: metavelvet uses velveth + meta-velvetg for metagenomic data

        #             To fix: Increase available memory or reduce test genome size further.
        #             """
        #             Test.@test_skip "metavelvet test skipped - insufficient resources"
        #         else
        #             rethrow(e)
        #         end
        #         # Clean up on failure
        #         rm(metavelvet_outdir, recursive=true, force=true)
        #     end
        # end
    end
end

Test.@testset "Long Read Metagenomic Assembly" begin
    mktempdir() do dir
        # Create reference genomes for long read metagenomic assembly (larger, more complex metagenome)
        meta_long_ref_fasta = joinpath(dir, "metagenomic_long_ref.fasta")

        # Create 3 diverse genomes of moderate sizes to balance complexity with memory usage
        rng_meta_long1 = StableRNGs.StableRNG(789)
        rng_meta_long2 = StableRNGs.StableRNG(790)
        rng_meta_long3 = StableRNGs.StableRNG(791)

        meta_long_genome1 = BioSequences.randdnaseq(rng_meta_long1, 12000)  # 12kb genome
        meta_long_genome2 = BioSequences.randdnaseq(rng_meta_long2, 10000)  # 10kb genome
        meta_long_genome3 = BioSequences.randdnaseq(rng_meta_long3, 8000)  # 8kb genome

        # Create FASTA records and write using Mycelia.write_fasta
        meta_long_fasta_record1 = FASTX.FASTA.Record("test_metagenomic_long_genome_1", meta_long_genome1)
        meta_long_fasta_record2 = FASTX.FASTA.Record("test_metagenomic_long_genome_2", meta_long_genome2)
        meta_long_fasta_record3 = FASTX.FASTA.Record("test_metagenomic_long_genome_3", meta_long_genome3)
        Mycelia.write_fasta(outfile = meta_long_ref_fasta,
            records = [
                meta_long_fasta_record1, meta_long_fasta_record2, meta_long_fasta_record3])

        # Simulate PacBio reads for metagenomic assembly (reuse for all metagenomic tests)
        meta_long_simulated_reads = Mycelia.simulate_pacbio_reads(
            fasta = meta_long_ref_fasta, quantity = "10x", quiet = true)

        # Decompress for assemblers
        meta_long_fastq = joinpath(dir, "meta_long_reads.fq")
        run(pipeline(`gunzip -c $(meta_long_simulated_reads)`, meta_long_fastq))

        Test.@testset "Long Read Metagenome Assembly - hifiasm-meta" begin
            hifiasm_meta_outdir = joinpath(dir, "hifiasm_meta_assembly")
            if isdir(hifiasm_meta_outdir)
                rm(hifiasm_meta_outdir, recursive = true)
            end
            try
                result = Mycelia.run_hifiasm_meta(
                    fastq = meta_long_fastq, outdir = hifiasm_meta_outdir,
                    bloom_filter = 0, read_selection = true, threads = threads)
                Test.@test result.outdir == hifiasm_meta_outdir
                expected_prefix = joinpath(hifiasm_meta_outdir, basename(meta_long_fastq) *
                                                                ".hifiasm_meta")
                Test.@test result.hifiasm_outprefix == expected_prefix
                rm(hifiasm_meta_outdir, recursive = true, force = true)
            catch e
                @error "hifiasm-meta test failed." exception=(e, catch_backtrace())
                Test.@test false
                rm(hifiasm_meta_outdir, recursive = true, force = true)
            end
        end

        Test.@testset "Long Read Metagenomic Assembly - metaFlye" begin
            # Create a smaller, specific dataset for metaFlye to manage memory usage
            metaflye_ref_fasta = joinpath(dir, "metaflye_specific_ref.fasta")
            rng_meta_flye1 = StableRNGs.StableRNG(800)
            rng_meta_flye2 = StableRNGs.StableRNG(801)
            meta_flye_genome1 = BioSequences.randdnaseq(rng_meta_flye1, 4000) # 4kb
            meta_flye_genome2 = BioSequences.randdnaseq(rng_meta_flye2, 3000) # 3kb
            meta_flye_record1 = FASTX.FASTA.Record("meta_flye_genome_1", meta_flye_genome1)
            meta_flye_record2 = FASTX.FASTA.Record("meta_flye_genome_2", meta_flye_genome2)
            Mycelia.write_fasta(outfile = metaflye_ref_fasta,
                records = [meta_flye_record1, meta_flye_record2])

            # Simulate reads with moderate coverage to ensure overlaps in tiny genomes
            meta_flye_simulated_reads = Mycelia.simulate_pacbio_reads(
                fasta = metaflye_ref_fasta, quantity = "30x", quiet = true)
            meta_flye_fastq = joinpath(dir, "meta_flye_reads.fq")
            run(pipeline(`gunzip -c $(meta_flye_simulated_reads)`, meta_flye_fastq))

            # Test metaFlye - clean up any existing directory first
            metaflye_outdir = joinpath(dir, "metaflye_assembly")
            if isdir(metaflye_outdir)
                rm(metaflye_outdir, recursive = true)
            end
            try
                result = Mycelia.run_metaflye(
                    fastq = meta_flye_fastq,
                    outdir = metaflye_outdir,
                    genome_size = "7k",
                    read_type = "pacbio-hifi",
                    min_overlap = 1000,
                    threads = threads
                )
                Test.@test result.outdir == metaflye_outdir
                Test.@test result.assembly == joinpath(metaflye_outdir, "assembly.fasta")
                Test.@test isfile(result.assembly)
                # Clean up after test
                rm(metaflye_outdir, recursive = true, force = true)
            catch e
                # ... (error handling)
                rm(metaflye_outdir, recursive = true, force = true)
                rethrow(e)
            end
        end

        Test.@testset "Long Read Metagenomic Assembly - metamdbg HiFi" begin
            # Test metamdbg with HiFi reads (PacBio)
            metamdbg_hifi_outdir = joinpath(dir, "metamdbg_hifi_assembly")
            if isdir(metamdbg_hifi_outdir)
                rm(metamdbg_hifi_outdir, recursive = true)
            end
            try
                result = Mycelia.run_metamdbg(
                    hifi_reads = meta_long_fastq, outdir = metamdbg_hifi_outdir,
                    abundance_min = 2, threads = threads)
                Test.@test result.outdir == metamdbg_hifi_outdir
                Test.@test !isempty(result.contigs)
                Test.@test !isempty(result.graph)
                Test.@test isfile(result.contigs)
                # Clean up after test
                rm(metamdbg_hifi_outdir, recursive = true, force = true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     metamdbg HiFi assembly failed due to resource constraints.
                #     Current test: 14kb genomes (8kb + 6kb), 12x coverage, ~168kb total sequence data

                #     Required resources for metamdbg:
                #     - Memory: ~2-4GB RAM minimum
                #     - CPU: 1-8 cores recommended
                #     - Disk: ~500MB temporary space
                #     - Note: metamdbg specializes in metagenomic long-read assembly with minimizers

                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "metamdbg HiFi test skipped - insufficient resources"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(metamdbg_hifi_outdir, recursive = true, force = true)
                rethrow(e)
            end
        end

        Test.@testset "Long Read Metagenomic Assembly - metamdbg ONT" begin
            # Test metamdbg with ONT reads - create ONT-specific reads
            # Simulate nanopore reads for ONT test
            ont_simulated_reads = Mycelia.simulate_nanopore_reads(
                fasta = meta_long_ref_fasta, quantity = "10x", quiet = true)
            ont_fastq = joinpath(dir, "meta_ont_reads.fq")
            run(pipeline(`gunzip -c $(ont_simulated_reads)`, ont_fastq))

            metamdbg_ont_outdir = joinpath(dir, "metamdbg_ont_assembly")
            if isdir(metamdbg_ont_outdir)
                rm(metamdbg_ont_outdir, recursive = true)
            end
            try
                result = Mycelia.run_metamdbg(
                    ont_reads = ont_fastq, outdir = metamdbg_ont_outdir,
                    abundance_min = 2, threads = threads)
                Test.@test result.outdir == metamdbg_ont_outdir
                Test.@test !isempty(result.contigs)
                Test.@test !isempty(result.graph)
                Test.@test isfile(result.contigs)
                # Clean up after test
                rm(metamdbg_ont_outdir, recursive = true, force = true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     metamdbg ONT assembly failed due to resource constraints.
                #     Current test: 14kb genomes (8kb + 6kb), 10x coverage, ~140kb total sequence data

                #     Required resources for metamdbg:
                #     - Memory: ~2-4GB RAM minimum
                #     - CPU: 1-8 cores recommended
                #     - Disk: ~500MB temporary space
                #     - Note: metamdbg specializes in metagenomic long-read assembly with minimizers

                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "metamdbg ONT test skipped - insufficient resources"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(metamdbg_ont_outdir, recursive = true, force = true)
                rethrow(e)
            end
        end
    end
end
