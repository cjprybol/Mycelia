# Legacy third-party assembler tests: short read isolate.
import Test
import Mycelia
import Graphs
import Random
import StableRNGs
import BioSequences
import FASTX

threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

Test.@testset "Short Read Isolate Assembly" begin
    mktempdir() do dir
        # Create small reference genome for isolate short read assembly
        ref_fasta = joinpath(dir, "isolate_ref.fasta")
        rng = StableRNGs.StableRNG(123)
        isolate_genome = BioSequences.randdnaseq(rng, 5000)  # 5kb genome for minimal memory usage
        
        # Create FASTA record and write using Mycelia.write_fasta
        fasta_record = FASTX.FASTA.Record("test_isolate_genome", isolate_genome)
        Mycelia.write_fasta(outfile=ref_fasta, records=[fasta_record])
        
        # Simulate Illumina paired-end reads with coverage for isolate assembly
        simulated_reads = Mycelia.simulate_illumina_reads(fasta=ref_fasta, coverage=10, rndSeed=123, quiet=true)
        
        # Extract the paired read files from the result
        fastq1 = simulated_reads.forward_reads
        fastq2 = simulated_reads.reverse_reads

        Test.@testset "Short Read Isolate Assembly - Spades" begin    
            # Test SPAdes - clean up any existing directory first
            spades_outdir = joinpath(dir, "spades_assembly")
            if isdir(spades_outdir)
                rm(spades_outdir, recursive=true)
            end
            try
                result = Mycelia.run_spades(fastq1=fastq1, fastq2=fastq2, outdir=spades_outdir, k_list="21", threads=threads, only_assembler=true)
                Test.@test result.outdir == spades_outdir
                Test.@test result.contigs == joinpath(spades_outdir, "contigs.fasta")
                Test.@test result.scaffolds == joinpath(spades_outdir, "scaffolds.fasta")
                Test.@test isfile(result.contigs)
                Test.@test isfile(result.scaffolds)
                # Clean up after test
                rm(spades_outdir, recursive=true, force=true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     SPAdes assembly failed due to resource constraints.
                #     Current test: 5kb genome, 15x coverage, ~75kb total sequence data
                    
                #     Required resources for SPAdes:
                #     - Memory: ~2-4GB RAM minimum
                #     - CPU: 1-8 cores recommended
                #     - Disk: ~500MB temporary space
                #     - Note: SPAdes is designed for single isolate genome assembly
                    
                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "SPAdes test skipped - insufficient resources"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(spades_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end

        Test.@testset "Short Read Isolate Assembly - SKESA" begin
            # Test SKESA - clean up any existing directory first  
            skesa_outdir = joinpath(dir, "skesa_assembly")
            if isdir(skesa_outdir)
                rm(skesa_outdir, recursive=true)
            end
            try
                result = Mycelia.run_skesa(fastq1=fastq1, fastq2=fastq2, outdir=skesa_outdir, threads=threads)
                Test.@test result.outdir == skesa_outdir
                Test.@test result.contigs == joinpath(skesa_outdir, "contigs.fa")
                Test.@test isfile(result.contigs)
                # Clean up after test
                rm(skesa_outdir, recursive=true, force=true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     SKESA assembly failed due to resource constraints.
                #     Current test: 5kb genome, 15x coverage, ~75kb total sequence data
                    
                #     Required resources for SKESA:
                #     - Memory: ~500MB-1GB RAM minimum
                #     - CPU: 1-4 cores recommended
                #     - Disk: ~50MB temporary space
                    
                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "SKESA test skipped - insufficient resources"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(skesa_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end

        # Test.@testset "Short Read Isolate Assembly - Ray" begin
        #     # Test Ray - clean up any existing directory first  
        #     ray_outdir = joinpath(dir, "ray_assembly")
        #     if isdir(ray_outdir)
        #         rm(ray_outdir, recursive=true)
        #     end
        #     try
        #         result = Mycelia.run_ray([fastq1, fastq2], outdir=ray_outdir, k=25)
        #         Test.@test result.outdir == ray_outdir
        #         Test.@test result.contigs == joinpath(ray_outdir, "Contigs.fasta")
        #         Test.@test result.scaffolds == joinpath(ray_outdir, "Scaffolds.fasta")
        #         Test.@test isfile(result.contigs)
        #         Test.@test isfile(result.scaffolds)
        #         # Clean up after test
        #         rm(ray_outdir, recursive=true, force=true)
        #     catch e
        #         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
        #             @warn """
        #             Ray assembly failed due to resource constraints.
        #             Current test: 5kb genome, 15x coverage, ~75kb total sequence data
                    
        #             Required resources for Ray:
        #             - Memory: ~1-3GB RAM minimum
        #             - CPU: 1-8 cores recommended (distributed assembler)
        #             - Disk: ~200MB temporary space
        #             - Note: Ray uses distributed de Bruijn graph approach
                    
        #             To fix: Increase available memory or reduce test genome size further.
        #             """
        #             Test.@test_skip "Ray test skipped - insufficient resources"
        #         else
        #             rethrow(e)
        #         end
        #         # Clean up on failure
        #         rm(ray_outdir, recursive=true, force=true)
        #     end
        # end

        # Test.@testset "Short Read Isolate Assembly - IDBA-UD" begin
        #     # Test IDBA-UD - clean up any existing directory first  
        #     idba_outdir = joinpath(dir, "idba_ud_assembly")
        #     if isdir(idba_outdir)
        #         rm(idba_outdir, recursive=true)
        #     end
        #     try
        #         result = Mycelia.run_idba_ud(fastq1=meta_fastq1, fastq2=meta_fastq2, outdir=idba_outdir)
        #         Test.@test result.outdir == idba_outdir
        #         Test.@test result.contigs == joinpath(idba_outdir, "contig.fa")
        #         Test.@test result.scaffolds == joinpath(idba_outdir, "scaffold.fa")
        #         Test.@test isfile(result.contigs)
        #         Test.@test isfile(result.scaffolds)
        #         # Clean up after test
        #         rm(idba_outdir, recursive=true, force=true)
        #     catch e
        #         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
        #             @warn """
        #             IDBA-UD assembly failed due to resource constraints.
        #             Current test: 5kb genome, 15x coverage, ~75kb total sequence data
                    
        #             Required resources for IDBA-UD:
        #             - Memory: ~1-2GB RAM minimum
        #             - CPU: 1-4 cores recommended
        #             - Disk: ~100MB temporary space
        #             - Note: IDBA-UD requires paired-end reads for metagenomic assembly
                    
        #             To fix: Increase available memory or reduce test genome size further.
        #             """
        #             Test.@test_skip "IDBA-UD test skipped - insufficient resources"
        #         else
        #             rethrow(e)
        #         end
        #         # Clean up on failure
        #         rm(idba_outdir, recursive=true, force=true)
        #     end
        # end

        Test.@testset "Short Read Isolate Assembly - Velvet" begin
            # Test Velvet - clean up any existing directory first  
            velvet_outdir = joinpath(dir, "velvet_assembly")
            if isdir(velvet_outdir)
                rm(velvet_outdir, recursive=true)
            end
            try
                result = Mycelia.run_velvet(fastq1, fastq2=fastq2, outdir=velvet_outdir, k=25)
                Test.@test result.outdir == velvet_outdir
                Test.@test result.contigs == joinpath(velvet_outdir, "contigs.fa")
                Test.@test isfile(result.contigs)
                # Clean up after test
                rm(velvet_outdir, recursive=true, force=true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     Velvet assembly failed due to resource constraints.
                #     Current test: 5kb genome, 15x coverage, ~75kb total sequence data
                    
                #     Required resources for Velvet:
                #     - Memory: ~500MB-2GB RAM minimum
                #     - CPU: 1-4 cores recommended
                #     - Disk: ~100MB temporary space
                #     - Note: Velvet uses classic de Bruijn graph approach with velveth+velvetg
                    
                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "Velvet test skipped - insufficient resources"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(velvet_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end
    end
end
