# Legacy third-party assembler tests: long read isolate.
import Test
import Mycelia
import Graphs
import Random
import StableRNGs
import BioSequences
import FASTX

threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

Test.@testset "Long Read Isolate Assembly" begin
    mktempdir() do dir
        Test.@testset "Long Read Isolate Assembly - Flye" begin
            # Test Flye - create reference genome and simulate ONT reads
            flye_ref_fasta = joinpath(dir, "flye_isolate_ref.fasta")
            rng_flye = StableRNGs.StableRNG(456)
            flye_genome = BioSequences.randdnaseq(rng_flye, 10000)  # 10kb genome for better assembly success
            
            # Create FASTA record and write using Mycelia.write_fasta
            flye_fasta_record = FASTX.FASTA.Record("flye_isolate_genome", flye_genome)
            Mycelia.write_fasta(outfile=flye_ref_fasta, records=[flye_fasta_record])

            # Simulate nanopore reads with 30x coverage (higher coverage for better assembly success)
            # flye_simulated_reads = Mycelia.simulate_nanopore_reads(fasta=flye_ref_fasta, quantity="15x", quiet=true)
            flye_simulated_reads = Mycelia.simulate_nanopore_reads(fasta=flye_ref_fasta, quantity="20x", quiet=true, seed=456)
            
            # Decompress for flye
            flye_fastq = joinpath(dir, "flye_reads.fq")
            run(pipeline(`gunzip -c $(flye_simulated_reads)`, flye_fastq))
        
            # Test Flye - clean up any existing directory first
            flye_outdir = joinpath(dir, "flye_assembly")
            if isdir(flye_outdir)
                rm(flye_outdir, recursive=true)
            end
            try
                result = Mycelia.run_flye(fastq=flye_fastq, outdir=flye_outdir, read_type="nano-hq", threads=threads)
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
            canu_genome = BioSequences.randdnaseq(rng_canu, 15000)  # 15kb genome - smaller for memory
            
            # Create FASTA record and write using Mycelia.write_fasta
            canu_fasta_record = FASTX.FASTA.Record("canu_test_genome", canu_genome)
            Mycelia.write_fasta(outfile=canu_ref_fasta, records=[canu_fasta_record])
            
            # Simulate PacBio reads with 15x coverage (sufficient for Canu's 10x minimum requirement)
            canu_simulated_reads = Mycelia.simulate_pacbio_reads(fasta=canu_ref_fasta, quantity="12x", quiet=true, seed=42)
            
            # Decompress for canu (canu expects uncompressed fastq)
            canu_fastq = joinpath(dir, "canu_reads.fq")
            run(pipeline(`gunzip -c $(canu_simulated_reads)`, canu_fastq))
            
            # Test Canu - clean up any existing directory first
            canu_outdir = joinpath(dir, "canu_assembly")
            if isdir(canu_outdir)
                rm(canu_outdir, recursive=true)
            end
            try
                result = Mycelia.run_canu(fastq=canu_fastq, outdir=canu_outdir, genome_size="15k", stopOnLowCoverage=8, threads=threads)
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
            else
                # Test hifiasm - create reference genome and simulate HiFi reads
                hifiasm_ref_fasta = joinpath(dir, "hifiasm_reference.fasta")
                rng_hifiasm = StableRNGs.StableRNG(333)
                hifiasm_genome = BioSequences.randdnaseq(rng_hifiasm, 10000)  # Reduced from 20kb to 10kb for faster runs
                
                # Create FASTA record and write using Mycelia.write_fasta
                hifiasm_fasta_record = FASTX.FASTA.Record("hifiasm_test_genome", hifiasm_genome)
                Mycelia.write_fasta(outfile=hifiasm_ref_fasta, records=[hifiasm_fasta_record])
                
                # Simulate PacBio reads with 10x coverage (reduced from 12x for memory efficiency)
                hifiasm_simulated_reads = Mycelia.simulate_pacbio_reads(fasta=hifiasm_ref_fasta, quantity="10x", quiet=true, seed=333)
                
                # Decompress for hifiasm
                hifiasm_fastq = joinpath(dir, "hifiasm_reads.fq")
                run(pipeline(`gunzip -c $(hifiasm_simulated_reads)`, hifiasm_fastq))
                
                # Test hifiasm - clean up any existing directory first
                hifiasm_outdir = joinpath(dir, "hifiasm_assembly")
                if isdir(hifiasm_outdir)
                    rm(hifiasm_outdir, recursive=true)
                end
                try
                    result = Mycelia.run_hifiasm(fastq=hifiasm_fastq, outdir=hifiasm_outdir, bloom_filter=0, threads=threads)
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
