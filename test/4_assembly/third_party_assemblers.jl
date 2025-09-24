# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/4_assembly/third_party_assemblers.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/third_party_assemblers.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import Graphs
import Random
import StableRNGs
import BioSequences
import FASTX

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
        simulated_reads = Mycelia.simulate_illumina_reads(fasta=ref_fasta, coverage=15, quiet=true)
        
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
                result = Mycelia.run_spades(fastq1=fastq1, fastq2=fastq2, outdir=spades_outdir)
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
                result = Mycelia.run_skesa(fastq1=fastq1, fastq2=fastq2, outdir=skesa_outdir)
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
        Mycelia.write_fasta(outfile=meta_ref_fasta, records=[meta_fasta_record1, meta_fasta_record2])
        
        # Simulate Illumina paired-end reads with coverage for metagenomic assembly
        meta_simulated_reads = Mycelia.simulate_illumina_reads(fasta=meta_ref_fasta, coverage=12, quiet=true)
        
        # Extract the paired read files from the result
        meta_fastq1 = meta_simulated_reads.forward_reads
        meta_fastq2 = meta_simulated_reads.reverse_reads

        Test.@testset "Short Read Metagenomic Assembly - MEGAHIT" begin
            # Test MEGAHIT - clean up any existing directory first
            megahit_outdir = joinpath(dir, "megahit_assembly")
            if isdir(megahit_outdir)
                rm(megahit_outdir, recursive=true)
            end
            try
                result = Mycelia.run_megahit(fastq1=meta_fastq1, fastq2=meta_fastq2, outdir=megahit_outdir)
                Test.@test result.outdir == megahit_outdir
                Test.@test result.contigs == joinpath(megahit_outdir, "final.contigs.fa")
                Test.@test isfile(result.contigs)
                # Clean up after test
                rm(megahit_outdir, recursive=true, force=true)
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
                rm(megahit_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end

        Test.@testset "Short Read Metagenomic Assembly - metaSPAdes" begin
    
            # Test metaSPAdes - clean up any existing directory first  
            spades_outdir = joinpath(dir, "metaspades_assembly")
            if isdir(spades_outdir)
                rm(spades_outdir, recursive=true)
            end
            try
                result = Mycelia.run_metaspades(fastq1=meta_fastq1, fastq2=meta_fastq2, outdir=spades_outdir)
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
                rm(spades_outdir, recursive=true, force=true)
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
    

Test.@testset "Long Read Isolate Assembly" begin
    mktempdir() do dir
        Test.@testset "Long Read Isolate Assembly - Flye" begin
            # Test Flye - create reference genome and simulate ONT reads
            flye_ref_fasta = joinpath(dir, "flye_isolate_ref.fasta")
            rng_flye = StableRNGs.StableRNG(456)
            flye_genome = BioSequences.randdnaseq(rng_flye, 50000)  # 50kb genome for better assembly success
            
            # Create FASTA record and write using Mycelia.write_fasta
            flye_fasta_record = FASTX.FASTA.Record("flye_isolate_genome", flye_genome)
            Mycelia.write_fasta(outfile=flye_ref_fasta, records=[flye_fasta_record])
            
            # Simulate nanopore reads with 20x coverage (higher coverage for better assembly success)
            flye_simulated_reads = Mycelia.simulate_nanopore_reads(fasta=flye_ref_fasta, quantity="20x", quiet=true)
            
            # Decompress for flye
            flye_fastq = joinpath(dir, "flye_reads.fq")
            run(pipeline(`gunzip -c $(flye_simulated_reads)`, flye_fastq))
        
            # Test Flye - clean up any existing directory first
            flye_outdir = joinpath(dir, "flye_assembly")
            if isdir(flye_outdir)
                rm(flye_outdir, recursive=true)
            end
            try
                result = Mycelia.run_flye(fastq=flye_fastq, outdir=flye_outdir, read_type="nano-raw")
                Test.@test result.outdir == flye_outdir
                Test.@test result.assembly == joinpath(flye_outdir, "assembly.fasta")
                Test.@test isfile(result.assembly)
                # Clean up after test
                rm(flye_outdir, recursive=true, force=true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     Flye assembly failed due to resource constraints.
                #     Current test: 50kb genome, 20x coverage, ~1MB total sequence data

                #     Required resources for Flye:
                #     - Memory: ~1-3GB RAM minimum
                #     - CPU: 1-4 cores recommended
                #     - Disk: ~500MB temporary space

                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "Flye test skipped - insufficient resources"
                # elseif contains(string(e), "No contigs were assembled") || contains(string(e), "Pipeline aborted") || contains(string(e), "assembly failed")
                #     @warn """
                #     Flye assembly failed to generate contigs despite increased parameters.
                #     Current test: 50kb genome, 20x coverage, ~1MB total sequence data

                #     This suggests the simulated reads may not be suitable for assembly.
                #     Consider adjusting read simulation parameters or using different test data.
                #     """
                #     Test.@test_skip "Flye test skipped - assembly failed to generate contigs"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(flye_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end

        Test.@testset "Long Read Isolate Assembly - Canu" begin

            # Test Canu - requires minimum 10x coverage
            # Create reference genome and simulate PacBio reads
            canu_ref_fasta = joinpath(dir, "canu_reference.fasta")
            rng_canu = StableRNGs.StableRNG(42)
            canu_genome = BioSequences.randdnaseq(rng_canu, 25000)  # 25kb genome - smaller for memory
            
            # Create FASTA record and write using Mycelia.write_fasta
            canu_fasta_record = FASTX.FASTA.Record("canu_test_genome", canu_genome)
            Mycelia.write_fasta(outfile=canu_ref_fasta, records=[canu_fasta_record])
            
            # Simulate PacBio reads with 15x coverage (sufficient for Canu's 10x minimum requirement)
            canu_simulated_reads = Mycelia.simulate_pacbio_reads(fasta=canu_ref_fasta, quantity="15x", quiet=true)
            
            # Decompress for canu (canu expects uncompressed fastq)
            canu_fastq = joinpath(dir, "canu_reads.fq")
            run(pipeline(`gunzip -c $(canu_simulated_reads)`, canu_fastq))
            
            # Test Canu - clean up any existing directory first
            canu_outdir = joinpath(dir, "canu_assembly")
            if isdir(canu_outdir)
                rm(canu_outdir, recursive=true)
            end
            try
                result = Mycelia.run_canu(fastq=canu_fastq, outdir=canu_outdir, genome_size="25k", stopOnLowCoverage=8)
                Test.@test result.outdir == canu_outdir
                Test.@test result.assembly == joinpath(canu_outdir, "canu_reads.contigs.fasta")
                Test.@test isfile(result.assembly)
                # Clean up after test
                rm(canu_outdir, recursive=true, force=true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed") || contains(string(e), "coverage")
                #     @warn """
                #     Canu assembly failed due to resource constraints or insufficient coverage.
                #     Current test: 25kb genome, 15x coverage, ~375kb total sequence data
                    
                #     Required resources for Canu:
                #     - Memory: ~4-8GB RAM minimum  
                #     - CPU: 4-8 cores recommended
                #     - Disk: ~1GB temporary space
                #     - Coverage: minimum 10x, prefer 20x+
                    
                #     To fix: Increase available memory/CPU or use smaller genome with higher coverage.
                #     """
                #     Test.@test_skip "Canu test skipped - insufficient resources"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(canu_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end

        Test.@testset "Long Read Isolate Assembly - Hifiasm" begin
            # Check if running in CI environment with resource constraints
            is_ci = haskey(ENV, "CI") || haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "TRAVIS") || haskey(ENV, "CIRCLECI")
            
            if is_ci
                # Skip Hifiasm tests in CI due to memory requirements
                @warn """
                Hifiasm test skipped in CI environment due to resource constraints.
                
                Hifiasm memory requirements:
                - Minimum: ~2-4GB RAM for even small test genomes
                - GitHub CI runners: 16GB total RAM but shared with system processes
                - Risk: Hifiasm's algorithmic overhead can exceed available memory
                
                This test runs locally but is automatically disabled in CI/CD environments
                to ensure reliable builds. To test Hifiasm functionality:
                
                Local testing:
                julia --project=. -e 'include("test/4_assembly/third_party_assemblers.jl")'
                
                Extended testing:
                julia --project=. run_extended_tests.jl tutorials
                """
                Test.@test_skip "hifiasm test skipped in CI environment - resource constraints"
            else
                # Test hifiasm - create reference genome and simulate HiFi reads
                hifiasm_ref_fasta = joinpath(dir, "hifiasm_reference.fasta")
                rng_hifiasm = StableRNGs.StableRNG(333)
                hifiasm_genome = BioSequences.randdnaseq(rng_hifiasm, 15000)  # Reduced from 20kb to 15kb for better CI compatibility
                
                # Create FASTA record and write using Mycelia.write_fasta
                hifiasm_fasta_record = FASTX.FASTA.Record("hifiasm_test_genome", hifiasm_genome)
                Mycelia.write_fasta(outfile=hifiasm_ref_fasta, records=[hifiasm_fasta_record])
                
                # Simulate PacBio reads with 10x coverage (reduced from 12x for memory efficiency)
                hifiasm_simulated_reads = Mycelia.simulate_pacbio_reads(fasta=hifiasm_ref_fasta, quantity="10x", quiet=true)
                
                # Decompress for hifiasm
                hifiasm_fastq = joinpath(dir, "hifiasm_reads.fq")
                run(pipeline(`gunzip -c $(hifiasm_simulated_reads)`, hifiasm_fastq))
                
                # Test hifiasm - clean up any existing directory first
                hifiasm_outdir = joinpath(dir, "hifiasm_assembly")
                if isdir(hifiasm_outdir)
                    rm(hifiasm_outdir, recursive=true)
                end
                try
                    result = Mycelia.run_hifiasm(fastq=hifiasm_fastq, outdir=hifiasm_outdir, bloom_filter=0)
                    Test.@test result.outdir == hifiasm_outdir
                    expected_prefix = joinpath(hifiasm_outdir, basename(hifiasm_fastq) * ".hifiasm")
                    Test.@test result.hifiasm_outprefix == expected_prefix
                    # Clean up after test
                    rm(hifiasm_outdir, recursive=true, force=true)
                catch e
                    # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                    #     @warn """
                    #     hifiasm assembly failed due to resource constraints.
                    #     Current test: 15kb genome, 10x coverage, ~150kb total sequence data
                        
                    #     Required resources for hifiasm:
                    #     - Memory: ~2-4GB RAM minimum
                    #     - CPU: 1-4 cores recommended
                    #     - Disk: ~500MB temporary space
                        
                    #     To fix: Increase available memory or reduce test genome size further.
                    #     """
                    #     Test.@test_skip "hifiasm test skipped - insufficient resources"
                    # else
                    #     rethrow(e)
                    # end
                    # Clean up on failure
                    rm(hifiasm_outdir, recursive=true, force=true)
                    rethrow(e)
                end
            end
        end
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

        meta_long_genome1 = BioSequences.randdnaseq(rng_meta_long1, 35000)  # 35kb genome
        meta_long_genome2 = BioSequences.randdnaseq(rng_meta_long2, 30000)  # 30kb genome
        meta_long_genome3 = BioSequences.randdnaseq(rng_meta_long3, 25000)  # 25kb genome

        # Create FASTA records and write using Mycelia.write_fasta
        meta_long_fasta_record1 = FASTX.FASTA.Record("test_metagenomic_long_genome_1", meta_long_genome1)
        meta_long_fasta_record2 = FASTX.FASTA.Record("test_metagenomic_long_genome_2", meta_long_genome2)
        meta_long_fasta_record3 = FASTX.FASTA.Record("test_metagenomic_long_genome_3", meta_long_genome3)
        Mycelia.write_fasta(outfile=meta_long_ref_fasta, records=[meta_long_fasta_record1, meta_long_fasta_record2, meta_long_fasta_record3])
        
        # Simulate PacBio reads for metagenomic assembly (reuse for all metagenomic tests)
        meta_long_simulated_reads = Mycelia.simulate_pacbio_reads(fasta=meta_long_ref_fasta, quantity="20x", quiet=true)
        
        # Decompress for assemblers
        meta_long_fastq = joinpath(dir, "meta_long_reads.fq")
        run(pipeline(`gunzip -c $(meta_long_simulated_reads)`, meta_long_fastq))
            
        # takes too long to run relative to other tools
        # Test.@testset "Long Read Metagenomic Assembly - hifiasm-meta" begin
        #     # Check if running in CI environment with resource constraints
        #     is_ci = haskey(ENV, "CI") || haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "TRAVIS") || haskey(ENV, "CIRCLECI")
            
        #     # if is_ci
        #     #     # Skip hifiasm-meta tests in CI due to memory requirements
        #     #     @warn """
        #     #     Hifiasm-meta test skipped in CI environment due to resource constraints.
                
        #     #     Hifiasm-meta memory requirements:
        #     #     - Minimum: ~3-6GB RAM for even small test genomes
        #     #     - GitHub CI runners: 16GB total RAM but shared with system processes
        #     #     - Risk: Metagenomic assembly algorithms have higher memory overhead
                
        #     #     This test runs locally but is automatically disabled in CI/CD environments
        #     #     to ensure reliable builds. To test hifiasm-meta functionality:
                
        #     #     Local testing:
        #     #     julia --project=. -e 'include("test/4_assembly/third_party_assemblers.jl")'
                
        #     #     Extended testing:
        #     #     julia --project=. run_extended_tests.jl tutorials
        #     #     """
        #     #     Test.@test_skip "hifiasm-meta test skipped in CI environment - resource constraints"
        #     # else
        #     # Test hifiasm-meta - metagenomic assembly
        #     hifiasm_meta_outdir = joinpath(dir, "hifiasm_meta_assembly")
        #     if isdir(hifiasm_meta_outdir)
        #         rm(hifiasm_meta_outdir, recursive=true)
        #     end
        #     try
        #         result = Mycelia.run_hifiasm_meta(fastq=meta_long_fastq, outdir=hifiasm_meta_outdir, bloom_filter=0, read_selection=true)
        #         Test.@test result.outdir == hifiasm_meta_outdir
        #         expected_prefix = joinpath(hifiasm_meta_outdir, basename(meta_long_fastq) * ".hifiasm_meta")
        #         Test.@test result.hifiasm_outprefix == expected_prefix
        #         # Clean up after test
        #         rm(hifiasm_meta_outdir, recursive=true, force=true)
        #     catch e
        #         # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed") || contains(string(e), "Segmentation fault")
        #         #     @warn """
        #         #     hifiasm-meta assembly failed due to resource constraints.
        #         #     Current test: 90kb total genome (3 genomes), 20x coverage, ~1.8MB total sequence data
                    
        #         #     Required resources for hifiasm-meta:
        #         #     - Memory: ~3-6GB RAM minimum (higher than regular hifiasm)
        #         #     - CPU: 1-4 cores recommended
        #         #     - Disk: ~1GB temporary space
        #         #     - Note: hifiasm-meta is optimized for metagenomic strain resolution
                    
        #         #     To fix: Increase available memory or reduce test genome size further.
        #         #     """
        #         #     Test.@test_skip "hifiasm-meta test skipped - insufficient resources"
        #         # else
        #         #     rethrow(e)
        #         # end
        #         # Clean up on failure
        #         rm(hifiasm_meta_outdir, recursive=true, force=true)
        #         rethrow(e)
        #     end
        # end
        # # end

        Test.@testset "Long Read Metagenomic Assembly - metaFlye" begin
            
            # Test metaFlye - clean up any existing directory first
            metaflye_outdir = joinpath(dir, "metaflye_assembly")
            if isdir(metaflye_outdir)
                rm(metaflye_outdir, recursive=true)
            end
            try
                result = Mycelia.run_metaflye(fastq=meta_long_fastq, outdir=metaflye_outdir, read_type="pacbio-raw")
                Test.@test result.outdir == metaflye_outdir
                Test.@test result.assembly == joinpath(metaflye_outdir, "assembly.fasta")
                Test.@test isfile(result.assembly)
                # Clean up after test
                rm(metaflye_outdir, recursive=true, force=true)
            catch e
                # if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
                #     @warn """
                #     metaFlye assembly failed due to resource constraints.
                #     Current test: 90kb total genome (3 genomes), 20x coverage, ~1.8MB total sequence data

                #     Required resources for metaFlye:
                #     - Memory: ~2-4GB RAM minimum (higher than regular Flye)
                #     - CPU: 1-4 cores recommended
                #     - Disk: ~1GB temporary space
                #     - Note: metaFlye is designed for metagenomic data with uneven coverage

                #     To fix: Increase available memory or reduce test genome size further.
                #     """
                #     Test.@test_skip "metaFlye test skipped - insufficient resources"
                # elseif contains(string(e), "No contigs were assembled") || contains(string(e), "Pipeline aborted") || contains(string(e), "assembly failed")
                #     @warn """
                #     metaFlye assembly failed to generate contigs despite increased parameters.
                #     Current test: 90kb total genome (3 genomes), 20x coverage, ~1.8MB total sequence data

                #     This suggests the simulated reads may not be suitable for metagenomic assembly.
                #     Consider adjusting read simulation parameters or using different test data.
                #     """
                #     Test.@test_skip "metaFlye test skipped - assembly failed to generate contigs"
                # else
                #     rethrow(e)
                # end
                # Clean up on failure
                rm(metaflye_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end

        Test.@testset "Long Read Metagenomic Assembly - metamdbg HiFi" begin
            # Test metamdbg with HiFi reads (PacBio)
            metamdbg_hifi_outdir = joinpath(dir, "metamdbg_hifi_assembly")
            if isdir(metamdbg_hifi_outdir)
                rm(metamdbg_hifi_outdir, recursive=true)
            end
            try
                result = Mycelia.run_metamdbg(hifi_reads=meta_long_fastq, outdir=metamdbg_hifi_outdir, abundance_min=2, threads=4)
                Test.@test result.outdir == metamdbg_hifi_outdir
                Test.@test !isempty(result.contigs)
                Test.@test !isempty(result.graph)
                Test.@test isfile(result.contigs)
                # Clean up after test
                rm(metamdbg_hifi_outdir, recursive=true, force=true)
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
                rm(metamdbg_hifi_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end

        Test.@testset "Long Read Metagenomic Assembly - metamdbg ONT" begin
            # Test metamdbg with ONT reads - create ONT-specific reads
            # Simulate nanopore reads for ONT test
            ont_simulated_reads = Mycelia.simulate_nanopore_reads(fasta=meta_long_ref_fasta, quantity="10x", quiet=true)
            ont_fastq = joinpath(dir, "meta_ont_reads.fq")
            run(pipeline(`gunzip -c $(ont_simulated_reads)`, ont_fastq))
            
            metamdbg_ont_outdir = joinpath(dir, "metamdbg_ont_assembly")
            if isdir(metamdbg_ont_outdir)
                rm(metamdbg_ont_outdir, recursive=true)
            end
            try
                result = Mycelia.run_metamdbg(ont_reads=ont_fastq, outdir=metamdbg_ont_outdir, abundance_min=2, threads=4)
                Test.@test result.outdir == metamdbg_ont_outdir
                Test.@test !isempty(result.contigs)
                Test.@test !isempty(result.graph)
                Test.@test isfile(result.contigs)
                # Clean up after test
                rm(metamdbg_ont_outdir, recursive=true, force=true)
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
                rm(metamdbg_ont_outdir, recursive=true, force=true)
                rethrow(e)
            end
        end
    end
end
    
#     Test.@testset "5. Hybrid Assembly" begin
#         Test.@testset "5a. Hybrid Isolate Assembly" begin
#             mktempdir() do dir
#                 # Create reference genome for hybrid isolate assembly
#                 hybrid_isolate_ref_fasta = joinpath(dir, "hybrid_isolate_ref.fasta")
#                 rng_hybrid_isolate = StableRNGs.StableRNG(789)
#                 hybrid_isolate_genome = BioSequences.randdnaseq(rng_hybrid_isolate, 3000)  # 3kb genome for absolute minimal memory
                
#                 # Create FASTA record and write using Mycelia.write_fasta
#                 hybrid_isolate_fasta_record = FASTX.FASTA.Record("hybrid_isolate_test_genome", hybrid_isolate_genome)
#                 Mycelia.write_fasta(outfile=hybrid_isolate_ref_fasta, records=[hybrid_isolate_fasta_record])
                
#                 # Simulate Illumina short reads with 10x coverage
#                 isolate_short_reads = Mycelia.simulate_illumina_reads(fasta=hybrid_isolate_ref_fasta, coverage=10)
                
#                 # Simulate nanopore long reads with 5x coverage
#                 isolate_long_reads_gz = Mycelia.simulate_nanopore_reads(fasta=hybrid_isolate_ref_fasta, quantity="5x")
                
#                 # Decompress long reads for unicycler
#                 isolate_long_reads = joinpath(dir, "isolate_long_reads.fq")
#                 run(pipeline(`gunzip -c $(isolate_long_reads_gz)`, isolate_long_reads))
            
#             # Test Unicycler - clean up any existing directory first
#             unicycler_outdir = joinpath(dir, "unicycler_assembly")
#             if isdir(unicycler_outdir)
#                 rm(unicycler_outdir, recursive=true)
#             end
#             try
#                 result = Mycelia.run_unicycler(short_1=isolate_short_reads.forward_reads, short_2=isolate_short_reads.reverse_reads, long_reads=isolate_long_reads, outdir=unicycler_outdir)
#                 Test.@test result.outdir == unicycler_outdir
#                 Test.@test result.assembly == joinpath(unicycler_outdir, "assembly.fasta")
#                 Test.@test isfile(result.assembly)
#                 # Clean up after test
#                 rm(unicycler_outdir, recursive=true, force=true)
#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed") || contains(string(e), "failed")
#                     @warn """
#                     Unicycler hybrid assembly failed due to resource constraints.
#                     Current test: 3kb genome, 10x short reads + 5x long reads, ~45kb total sequence data
                    
#                     Required resources for Unicycler:
#                     - Memory: ~4-8GB RAM minimum (SPAdes is memory-intensive)
#                     - CPU: 4-8 cores recommended  
#                     - Disk: ~1-2GB temporary space
#                     - Note: Unicycler internally runs SPAdes which has high memory requirements
                    
#                     To fix: This assembler requires significant computational resources.
#                     Consider running on a machine with more memory, or skip hybrid assembly testing.
                    
#                     Alternative: Use Flye or Canu for long-read-only assembly if resources are limited.
#                     """
#                     Test.@test_skip "Unicycler test skipped - insufficient resources (memory-intensive hybrid assembler)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(unicycler_outdir, recursive=true, force=true)
#             end
#         end
        
#         Test.@testset "5b. Hybrid Metagenomic Assembly" begin
#             mktempdir() do dir
#                 # Create reference genome for hybrid metagenomic assembly
#                 hybrid_meta_ref_fasta = joinpath(dir, "hybrid_meta_ref.fasta")
#                 rng_hybrid_meta = StableRNGs.StableRNG(987)
#                 hybrid_meta_genome = BioSequences.randdnaseq(rng_hybrid_meta, 3000)  # 3kb genome for metagenomic hybrid
                
#                 # Create FASTA record and write using Mycelia.write_fasta
#                 hybrid_meta_fasta_record = FASTX.FASTA.Record("hybrid_meta_test_genome", hybrid_meta_genome)
#                 Mycelia.write_fasta(outfile=hybrid_meta_ref_fasta, records=[hybrid_meta_fasta_record])
                
#                 # Simulate Illumina short reads with 10x coverage
#                 meta_short_reads = Mycelia.simulate_illumina_reads(fasta=hybrid_meta_ref_fasta, coverage=10)
                
#                 # Simulate nanopore long reads with 5x coverage
#                 meta_long_reads_gz = Mycelia.simulate_nanopore_reads(fasta=hybrid_meta_ref_fasta, quantity="5x")
                
#                 # Decompress long reads for HyLight
#                 meta_long_reads = joinpath(dir, "meta_long_reads.fq")
#                 run(pipeline(`gunzip -c $(meta_long_reads_gz)`, meta_long_reads))
            
#             # Test HyLight - hybrid strain-resolved metagenomic assembly
#             hylight_outdir = joinpath(dir, "hylight_assembly")
#             if isdir(hylight_outdir)
#                 rm(hylight_outdir, recursive=true)
#             end
#             try
#                 result = Mycelia.run_hylight(meta_short_reads.forward_reads, meta_short_reads.reverse_reads, meta_long_reads, outdir=hylight_outdir)
#                 Test.@test result.outdir == hylight_outdir
#                 Test.@test result.strain_assemblies == joinpath(hylight_outdir, "strain_assemblies")
#                 Test.@test isdir(result.strain_assemblies)
#                 # Clean up after test
#                 rm(hylight_outdir, recursive=true, force=true)
#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed") || contains(string(e), "failed")
#                     @warn """
#                     HyLight hybrid assembly failed due to resource constraints.
#                     Current test: 3kb genome, 10x short reads + 5x long reads, ~45kb total sequence data
                    
#                     Required resources for HyLight:
#                     - Memory: ~6-12GB RAM minimum (strain-resolved metagenomic hybrid assembly)
#                     - CPU: 4-16 cores recommended
#                     - Disk: ~2-4GB temporary space
#                     - Note: HyLight performs strain resolution on hybrid assemblies
                    
#                     To fix: This is a very resource-intensive assembler requiring significant computational resources.
#                     Consider running on a high-memory machine or skip strain-resolved hybrid assembly testing.
#                     """
#                     Test.@test_skip "HyLight test skipped - insufficient resources (very memory-intensive strain-resolved hybrid assembler)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(hylight_outdir, recursive=true, force=true)
#             end
#         end
#     end
#     Test.@testset "6. Probabilistic Assembly (Mycelia)" begin
#         mktempdir() do dir
#             fastq = joinpath(dir, "reads.fq")
#             open(fastq, "w") do io
#                 println(io, "@r1")
#                 println(io, "ACGTACGTACGTACGTACGT")
#                 println(io, "+")
#                 println(io, "IIIIIIIIIIIIIIIIIIII")
#                 println(io, "@r2")
#                 println(io, "CGTACGTACGTACGTACGTA")
#                 println(io, "+")
#                 println(io, "IIIIIIIIIIIIIIIIIIII")
#             end
            
#             # Test string graph building
#             graph = Mycelia.string_to_ngram_graph(s="ACGTACGTACGTACGTACGT", n=5)
#             Test.@test Graphs.nv(graph) > 0
            
#             # Test Viterbi error correction functions exist
#             Test.@test hasmethod(Mycelia.viterbi_maximum_likelihood_traversals, (Any,))
#             Test.@test isdefined(Mycelia, :polish_fastq)
#         end
#     end
    
#     Test.@testset "7. Assembly merging" begin
#     end
#     Test.@testset "8. Polishing & Error Correction" begin
#         mktempdir() do dir
#             # Create reference genome and assembly for polishing testing
#             polish_ref_fasta = joinpath(dir, "polish_ref.fasta")
#             rng_polish = StableRNGs.StableRNG(222)
#             polish_genome = BioSequences.randdnaseq(rng_polish, 4000)  # 4kb genome for polishing
            
#             # Create FASTA record and write using Mycelia.write_fasta
#             polish_fasta_record = FASTX.FASTA.Record("polish_test_genome", polish_genome)
#             Mycelia.write_fasta(outfile=polish_ref_fasta, records=[polish_fasta_record])
            
#             # Simulate PacBio reads for polishing (Apollo works with PacBio data)
#             polish_simulated_reads = Mycelia.simulate_pacbio_reads(fasta=polish_ref_fasta, quantity="15x")
            
#             # Decompress reads for polishing
#             polish_fastq = joinpath(dir, "polish_reads.fq")
#             run(pipeline(`gunzip -c $(polish_simulated_reads)`, polish_fastq))
            
#             # Generate a draft assembly using Flye (create something to polish)
#             draft_assembly_outdir = joinpath(dir, "draft_assembly")
#             draft_assembly_fasta = joinpath(draft_assembly_outdir, "assembly.fasta")
            
#             try
#                 # Generate draft assembly for polishing
#                 Mycelia.run_flye(fastq=polish_fastq, outdir=draft_assembly_outdir, genome_size="4k", read_type="pacbio-raw")
                
#                 if isfile(draft_assembly_fasta)
#                     # Test Apollo - clean up any existing directory first  
#                     apollo_outdir = joinpath(dir, "apollo_polishing")
#                     if isdir(apollo_outdir)
#                         rm(apollo_outdir, recursive=true)
#                     end
                    
#                     try
#                         result = Mycelia.run_apollo(draft_assembly_fasta, polish_fastq, outdir=apollo_outdir)
#                         Test.@test result.outdir == apollo_outdir
#                         Test.@test result.polished_assembly == joinpath(apollo_outdir, basename(draft_assembly_fasta, ".fasta") * "_polished.fasta")
#                         Test.@test isfile(result.polished_assembly)
#                         # Clean up after test
#                         rm(apollo_outdir, recursive=true, force=true)
#                     catch e
#                         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                             @warn """
#                             Apollo polishing failed due to resource constraints.
#                             Current test: 4kb genome, 15x coverage, ~60kb total sequence data
                            
#                             Required resources for Apollo:
#                             - Memory: ~2-4GB RAM minimum (includes minimap2 + samtools + HMM polishing)
#                             - CPU: 2-8 cores recommended
#                             - Disk: ~500MB temporary space
#                             - Note: Apollo performs HMM-based assembly polishing
                            
#                             To fix: Increase available memory or reduce genome/coverage complexity.
#                             """
#                             Test.@test_skip "Apollo test skipped - insufficient resources"
#                         else
#                             rethrow(e)
#                         end
#                         # Clean up on failure
#                         rm(apollo_outdir, recursive=true, force=true)
#                     end
#                 else
#                     @warn "Draft assembly not generated by Flye - skipping Apollo test"
#                     Test.@test_skip "Apollo test skipped - no draft assembly available"
#                 end
                
#                 # Clean up draft assembly
#                 rm(draft_assembly_outdir, recursive=true, force=true)
                
#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                     @warn """
#                     Flye assembly for Apollo test failed due to resource constraints.
#                     Cannot generate draft assembly needed for polishing testing.
#                     """
#                     Test.@test_skip "Apollo test skipped - cannot generate draft assembly (insufficient resources)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(draft_assembly_outdir, recursive=true, force=true)
#             end
            
#             # Test Homopolish - can reuse the same setup as Apollo
#             try
#                 # Generate draft assembly for Homopolish testing (reuse same setup)
#                 Mycelia.run_flye(fastq=polish_fastq, outdir=draft_assembly_outdir, genome_size="4k", read_type="pacbio-raw")
                
#                 if isfile(draft_assembly_fasta)
#                     # Test Homopolish - clean up any existing directory first  
#                     homopolish_outdir = joinpath(dir, "homopolish_polishing")
#                     if isdir(homopolish_outdir)
#                         rm(homopolish_outdir, recursive=true)
#                     end
                    
#                     try
#                         result = Mycelia.run_homopolish(draft_assembly_fasta, polish_fastq, outdir=homopolish_outdir)
#                         Test.@test result.outdir == homopolish_outdir
#                         Test.@test result.polished_assembly == joinpath(homopolish_outdir, basename(draft_assembly_fasta, ".fasta") * "_homopolished.fasta")
#                         Test.@test isfile(result.polished_assembly)
#                         # Clean up after test
#                         rm(homopolish_outdir, recursive=true, force=true)
#                     catch e
#                         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                             @warn """
#                             Homopolish polishing failed due to resource constraints.
#                             Current test: 4kb genome, 15x coverage, ~60kb total sequence data
                            
#                             Required resources for Homopolish:
#                             - Memory: ~2-3GB RAM minimum (reference-based homopolymer correction)
#                             - CPU: 1-8 cores recommended
#                             - Disk: ~500MB temporary space
#                             - Note: Homopolish performs reference-based homopolymer error correction
                            
#                             To fix: Increase available memory or reduce genome/coverage complexity.
#                             """
#                             Test.@test_skip "Homopolish test skipped - insufficient resources"
#                         else
#                             rethrow(e)
#                         end
#                         # Clean up on failure
#                         rm(homopolish_outdir, recursive=true, force=true)
#                     end
#                 else
#                     @warn "Draft assembly not generated by Flye - skipping Homopolish test"
#                     Test.@test_skip "Homopolish test skipped - no draft assembly available"
#                 end
                
#                 # Clean up draft assembly
#                 rm(draft_assembly_outdir, recursive=true, force=true)
                
#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                     @warn """
#                     Flye assembly for Homopolish test failed due to resource constraints.
#                     Cannot generate draft assembly needed for polishing testing.
#                     """
#                     Test.@test_skip "Homopolish test skipped - cannot generate draft assembly (insufficient resources)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(draft_assembly_outdir, recursive=true, force=true)
#             end
#         end
#     end
#     Test.@testset "6. Strain resolution" begin
#         mktempdir() do dir
#             # Create two related reference strains with realistic variations
#             base_ref_fasta = joinpath(dir, "base_strain.fasta")
#             variant_ref_fasta = joinpath(dir, "variant_strain.fasta")
            
#             rng_strain = StableRNGs.StableRNG(111)
            
#             # Create base strain (5kb genome)
#             base_genome = BioSequences.randdnaseq(rng_strain, 5000)
#             base_fasta_record = FASTX.FASTA.Record("base_strain", base_genome)
#             Mycelia.write_fasta(outfile=base_ref_fasta, records=[base_fasta_record])
            
#             # Create variant strain by introducing realistic variations
#             variant_genome = copy(base_genome)
            
#             # Introduce SNVs every 500bp (1% divergence)
#             for i in 500:500:length(variant_genome)
#                 if i <= length(variant_genome)
#                     original_base = variant_genome[i]
#                     # Change to a different base
#                     new_bases = filter(b -> b != original_base, [BioSequences.DNA_A, BioSequences.DNA_T, BioSequences.DNA_G, BioSequences.DNA_C])
#                     variant_genome[i] = rand(rng_strain, new_bases)
#                 end
#             end
            
#             # Introduce small indels (deletions of 1-3bp every 1000bp)
#             positions_to_delete = collect(1000:1000:length(variant_genome)-10)
#             for pos in reverse(positions_to_delete)  # reverse to maintain positions
#                 if pos + 2 <= length(variant_genome)
#                     deleteat!(variant_genome, pos:pos+1)  # delete 2bp
#                 end
#             end
            
#             variant_fasta_record = FASTX.FASTA.Record("variant_strain", variant_genome)
#             Mycelia.write_fasta(outfile=variant_ref_fasta, records=[variant_fasta_record])
            
#             # Simulate reads from each strain with uneven coverage (3:1 ratio)
#             base_reads = Mycelia.simulate_nanopore_reads(fasta=base_ref_fasta, quantity="12x")  # Higher coverage
#             variant_reads = Mycelia.simulate_nanopore_reads(fasta=variant_ref_fasta, quantity="4x")  # Lower coverage
            
#             # Decompress and combine reads for mixed community
#             base_fastq = joinpath(dir, "base_strain_reads.fq")
#             variant_fastq = joinpath(dir, "variant_strain_reads.fq")
#             mixed_fastq = joinpath(dir, "mixed_strain_reads.fq")
            
#             run(pipeline(`gunzip -c $(base_reads)`, base_fastq))
#             run(pipeline(`gunzip -c $(variant_reads)`, variant_fastq))
            
#             # Combine reads to simulate mixed community
#             run(pipeline(`cat $(base_fastq) $(variant_fastq)`, mixed_fastq))
            
#             # Generate assembly graph using metaFlye (metagenomic long-read assembler)
#             metaflye_for_graph_outdir = joinpath(dir, "metaflye_for_graph")
#             assembly_graph_gfa = joinpath(metaflye_for_graph_outdir, "assembly_graph.gfa")
            
#             try
#                 # Run metaFlye to generate assembly graph from mixed strain reads
#                 Mycelia.run_metaflye(fastq=mixed_fastq, outdir=metaflye_for_graph_outdir, genome_size="5k", read_type="nano-raw")
                
#                 # Test STRONG - clean up any existing directory first  
#                 strong_outdir = joinpath(dir, "strong_assembly")
#                 if isdir(strong_outdir)
#                     rm(strong_outdir, recursive=true)
#                 end
                
#                 if isfile(assembly_graph_gfa)
#                     try
#                         result = Mycelia.run_strong(assembly_graph_gfa, mixed_fastq, outdir=strong_outdir, nb_strains=2)
#                         Test.@test result.outdir == strong_outdir
#                         Test.@test result.strain_unitigs == joinpath(strong_outdir, "strain_unitigs.fasta")
#                         # Clean up after test
#                         rm(strong_outdir, recursive=true, force=true)
#                     catch e
#                         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                             @warn """
#                             STRONG strain resolution failed due to resource constraints.
#                             Current test: 2 strains (5kb each), mixed 12x+4x coverage, ~80kb total sequence data
                            
#                             Required resources for STRONG:
#                             - Memory: ~3-6GB RAM minimum (strain resolution is compute-intensive)
#                             - CPU: 2-8 cores recommended
#                             - Disk: ~1GB temporary space
#                             - Note: STRONG performs strain-aware assembly graph traversal
                            
#                             To fix: Increase available memory or reduce genome/coverage complexity.
#                             """
#                             Test.@test_skip "STRONG test skipped - insufficient resources"
#                         else
#                             rethrow(e)
#                         end
#                         # Clean up on failure
#                         rm(strong_outdir, recursive=true, force=true)
#                     end
#                 else
#                     @warn "Assembly graph not generated by metaFlye - skipping STRONG test"
#                     Test.@test_skip "STRONG test skipped - no assembly graph available"
#                 end
                
#                 # Clean up metaFlye output
#                 rm(metaflye_for_graph_outdir, recursive=true, force=true)
                
#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                     @warn """
#                     metaFlye assembly for STRONG test failed due to resource constraints.
#                     Cannot generate assembly graph needed for strain resolution testing.
                    
#                     This test requires:
#                     1. metaFlye assembly to generate GFA graph (~3-4GB RAM)
#                     2. STRONG strain resolution on the graph (~3-6GB RAM)
                    
#                     Total resources needed: ~5-8GB RAM, 2-8 cores
#                     """
#                     Test.@test_skip "STRONG test skipped - cannot generate assembly graph (insufficient resources)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(metaflye_for_graph_outdir, recursive=true, force=true)
#             end
            
#             # Test Strainy - requires an assembly FASTA file
#             # We can reuse the metaFlye assembly if it was generated successfully above
#             metaflye_assembly = joinpath(dir, "metaflye_for_strainy", "assembly.fasta")
            
#             # Generate assembly using metaFlye for Strainy testing
#             metaflye_for_strainy_outdir = joinpath(dir, "metaflye_for_strainy")
            
#             try
#                 # Run metaFlye to generate assembly FASTA for Strainy
#                 Mycelia.run_metaflye(fastq=mixed_fastq, outdir=metaflye_for_strainy_outdir, genome_size="5k", read_type="nano-raw")
                
#                 if isfile(metaflye_assembly)
#                     # Test Strainy - clean up any existing directory first  
#                     strainy_outdir = joinpath(dir, "strainy_assembly")
#                     if isdir(strainy_outdir)
#                         rm(strainy_outdir, recursive=true)
#                     end
                    
#                     try
#                         result = Mycelia.run_strainy(metaflye_assembly, mixed_fastq, outdir=strainy_outdir, mode="phase")
#                         Test.@test result.outdir == strainy_outdir
#                         Test.@test result.strain_assemblies == joinpath(strainy_outdir, "strain_assemblies.fasta")
#                         # Clean up after test
#                         rm(strainy_outdir, recursive=true, force=true)
#                     catch e
#                         if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                             @warn """
#                             Strainy strain phasing failed due to resource constraints.
#                             Current test: 2 strains (5kb each), mixed 12x+4x coverage, ~80kb total sequence data
                            
#                             Required resources for Strainy:
#                             - Memory: ~2-4GB RAM minimum (includes minimap2 + samtools)
#                             - CPU: 2-8 cores recommended
#                             - Disk: ~1GB temporary space
#                             - Note: Strainy performs strain phasing from long reads mapped to assembly
                            
#                             To fix: Increase available memory or reduce genome/coverage complexity.
#                             """
#                             Test.@test_skip "Strainy test skipped - insufficient resources"
#                         else
#                             rethrow(e)
#                         end
#                         # Clean up on failure
#                         rm(strainy_outdir, recursive=true, force=true)
#                     end
#                 else
#                     @warn "Assembly not generated by metaFlye - skipping Strainy test"
#                     Test.@test_skip "Strainy test skipped - no assembly available"
#                 end
                
#                 # Clean up metaFlye output for Strainy
#                 rm(metaflye_for_strainy_outdir, recursive=true, force=true)
                
#             catch e
#                 if isa(e, ProcessFailedException) || contains(string(e), "memory") || contains(string(e), "Memory") || contains(string(e), "killed")
#                     @warn """
#                     metaFlye assembly for Strainy test failed due to resource constraints.
#                     Cannot generate assembly needed for strain phasing testing.
#                     """
#                     Test.@test_skip "Strainy test skipped - cannot generate assembly (insufficient resources)"
#                 else
#                     rethrow(e)
#                 end
#                 # Clean up on failure
#                 rm(metaflye_for_strainy_outdir, recursive=true, force=true)
#             end
#         end
#     end
#     Test.@testset "9. Validation & Quality Control" begin
#     end
# end
