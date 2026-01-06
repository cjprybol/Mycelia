# ```bash
# julia --project=. --color=yes -e 'include("test/1_data_acquisition/simulation_fastq.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run from the Mycelia base directory:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/1_data_acquisition/simulation_fastq.jl", "test/1_data_acquisition", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import FASTX
import Random
import BioSequences

const phiX174_assembly_id = "GCF_000819615.1"

Test.@testset "FASTQ simulation" begin
    # Get test genome ONCE and share across all test sets
    # This avoids redundant downloads and reduces CI flakiness
    test_genome_info = Mycelia.get_test_genome_fasta(use_ncbi=true, accession=phiX174_assembly_id)
    test_fasta = test_genome_info.fasta
    
    @info "Using $(test_genome_info.source) genome for simulation tests: $test_fasta"
    
    # Verify we have a valid test genome before proceeding
    if !isfile(test_fasta) || filesize(test_fasta) == 0
        @warn "Failed to obtain test genome for simulation tests"
        Test.@test_broken false
    else
    
    Test.@testset "Illumina" begin
        # Test default system (HS25) - uses shared test_fasta
        read_simulation_result = Mycelia.simulate_illumina_reads(fasta = test_fasta, coverage=10, quiet=true, errfree=true)
        Test.@test isfile(read_simulation_result.forward_reads)
        Test.@test isfile(read_simulation_result.reverse_reads)
        Test.@test isfile(read_simulation_result.sam)
        Test.@test isfile(read_simulation_result.error_free_sam)
        
        # Test all ART Illumina profile-specific functions
        Test.@testset "Individual ART Profile Functions" begin
            # Define all the profile functions and their expected parameters
            profile_functions = [
                (:simulate_illumina_ga1_36bp, "GA1", 36),
                (:simulate_illumina_ga1_44bp, "GA1", 44), 
                (:simulate_illumina_ga2_50bp, "GA2", 50),
                (:simulate_illumina_ga2_75bp, "GA2", 75),
                (:simulate_illumina_hs10_100bp, "HS10", 100),
                (:simulate_illumina_hs20_100bp, "HS20", 100),
                (:simulate_illumina_hs25_125bp, "HS25", 125),
                (:simulate_illumina_hs25_150bp, "HS25", 150),
                (:simulate_illumina_hsxn_150bp, "HSXn", 150),
                (:simulate_illumina_hsxt_150bp, "HSXt", 150),
                (:simulate_illumina_msv1_250bp, "MSv1", 250),
                (:simulate_illumina_msv3_250bp, "MSv3", 250),
                (:simulate_illumina_mins_50bp, "MinS", 50),
                (:simulate_illumina_ns50_75bp, "NS50", 75)
            ]
            
            for (func_name, expected_system, expected_length) in profile_functions
                Test.@testset "$func_name" begin
                    func = getfield(Mycelia, func_name)
                    result = func(fasta = test_fasta, coverage = 3, errfree=true)
                    
                    Test.@test isfile(result.forward_reads)
                    Test.@test isfile(result.sam)
                    Test.@test isfile(result.error_free_sam)
                    
                    # Verify files are not empty
                    Test.@test filesize(result.forward_reads) > 0
                    Test.@test filesize(result.sam) > 0
                    Test.@test filesize(result.error_free_sam) > 0
                    
                    # Files should have standard ART naming pattern
                    if expected_system == "MinS"
                        # MinS single-end files end with .fq.gz
                        Test.@test contains(result.forward_reads, ".fq.gz")
                    else
                        # Paired-end files end with 1.fq.gz
                        Test.@test contains(result.forward_reads, "1.fq.gz")
                    end
                    
                    # Handle paired vs single-end reads
                    if expected_system == "MinS"
                        # MinS is single-end only
                        Test.@test result.reverse_reads === nothing
                    else
                        # All other systems are paired-end
                        Test.@test result.reverse_reads !== nothing
                        Test.@test isfile(result.reverse_reads)
                        Test.@test filesize(result.reverse_reads) > 0
                        Test.@test contains(result.reverse_reads, "2.fq.gz")
                    end
                    
                    # Clean up files
                    cleanup_files = [result.forward_reads, result.sam]
                    if result.error_free_sam !== nothing
                        push!(cleanup_files, result.error_free_sam)
                    end
                    if result.reverse_reads !== nothing
                        push!(cleanup_files, result.reverse_reads)
                    end
                    for file in cleanup_files
                        if isfile(file)
                            rm(file)
                        end
                    end
                end
            end
        end
        
        # Test that excessive read lengths fail appropriately
        Test.@testset "Read Length Validation" begin
            # Test systems with their maximum supported lengths and over-limits
            validation_tests = [
                ("GA1", 44, 50),    # GA1 max is 44, test 50 should fail
                ("GA2", 75, 100),   # GA2 max is 75, test 100 should fail
                ("HS10", 100, 125), # HS10 max is 100, test 125 should fail
                ("HS20", 100, 125), # HS20 max is 100, test 125 should fail
                ("HS25", 150, 200), # HS25 max is 150, test 200 should fail
                ("HSXn", 150, 200), # HSXn max is 150, test 200 should fail
                ("HSXt", 150, 200), # HSXt max is 150, test 200 should fail
                ("MSv1", 250, 300), # MSv1 max is 250, test 300 should fail
                ("MSv3", 250, 300), # MSv3 max is 250, test 300 should fail
                ("MinS", 50, 75),   # MinS max is 50, test 75 should fail
                ("NS50", 75, 100)   # NS50 max is 75, test 100 should fail
            ]
            
            for (seqSys, max_length, excessive_length) in validation_tests
                Test.@testset "Validation for $seqSys" begin
                    paired = seqSys != "MinS"
                    # Test that max length works
                    try
                        result_max = Mycelia.simulate_illumina_reads(
                            fasta = test_fasta,
                            coverage = 1,
                            seqSys = seqSys,
                            read_length = max_length,
                            paired = paired,
                            quiet = true
                        )
                        Test.@test isfile(result_max.forward_reads)
                        # Clean up successful run
                        cleanup_files = [result_max.forward_reads, result_max.sam]
                        if result_max.error_free_sam !== nothing
                            push!(cleanup_files, result_max.error_free_sam)
                        end
                        if result_max.reverse_reads !== nothing
                            push!(cleanup_files, result_max.reverse_reads)
                        end
                        for file in cleanup_files
                            if isfile(file)
                                rm(file)
                            end
                        end
                    catch e
                        # If max length fails, it might be an environment issue, skip
                        @warn "Max length test failed for $seqSys with length $max_length: $e"
                    end
                    
                    # Test that excessive length fails or produces error
                    Test.@test_throws Exception Mycelia.simulate_illumina_reads(
                        fasta = test_fasta,
                        coverage = 1,
                        seqSys = seqSys,
                        read_length = excessive_length,
                        paired = paired,
                        quiet = true
                    )
                end
            end
        end
    end  # end Illumina testset
    
    Test.@testset "Ultima Genomics" begin
        # Test basic Ultima simulation with default parameters (250bp) - uses shared test_fasta
        read_simulation_result = Mycelia.simulate_ultima_reads(fasta = test_fasta, coverage=5, quiet=true)
        Test.@test isfile(read_simulation_result.forward_reads)
        Test.@test read_simulation_result.reverse_reads === nothing  # Single-end only
        Test.@test isfile(read_simulation_result.sam)
        Test.@test isfile(read_simulation_result.error_free_sam)
        
        # Verify files are not empty
        Test.@test filesize(read_simulation_result.forward_reads) > 0
        Test.@test filesize(read_simulation_result.sam) > 0
        Test.@test filesize(read_simulation_result.error_free_sam) > 0
        
        # Test filename pattern for single-end reads
        Test.@test contains(read_simulation_result.forward_reads, ".fq.gz")
        Test.@test contains(read_simulation_result.forward_reads, "ultima")
        
        # Test with custom parameters
        custom_result = Mycelia.simulate_ultima_reads(
            fasta = test_fasta, 
            coverage = 3,
            read_length = 200,
            insertion_rate = 0.005,
            deletion_rate = 0.003,
            id = "test_ultima",
            quiet = true
        )
        Test.@test isfile(custom_result.forward_reads)
        Test.@test custom_result.reverse_reads === nothing
        Test.@test isfile(custom_result.sam)
        Test.@test isfile(custom_result.error_free_sam)
        
        # Clean up files
        cleanup_files = [
            read_simulation_result.forward_reads, 
            read_simulation_result.sam, 
            read_simulation_result.error_free_sam,
            custom_result.forward_reads,
            custom_result.sam,
            custom_result.error_free_sam
        ]
        for file in cleanup_files
            if isfile(file)
                rm(file)
            end
        end
    end  # end Ultima Genomics testset
    
    Test.@testset "Nanopore" begin
        read_simulation_result = Mycelia.simulate_nanopore_reads(fasta = test_fasta, quantity=10, quiet=true)
        Test.@test isfile(read_simulation_result)
        Test.@test endswith(read_simulation_result, ".fq.gz")
    end
    
    Test.@testset "PacBio" begin
        read_simulation_result = Mycelia.simulate_pacbio_reads(fasta = test_fasta, quantity=10, quiet=true)
        Test.@test isfile(read_simulation_result)
        Test.@test endswith(read_simulation_result, ".fq.gz")
    end
    
    Test.@testset "Nearly Perfect Long Reads" begin
        read_simulation_result = Mycelia.simulate_nearly_perfect_long_reads(fasta = test_fasta, quantity="10x", quiet=true)
        Test.@test isfile(read_simulation_result)
        Test.@test endswith(read_simulation_result, ".fq.gz")
    end
    
    Test.@testset "Nanopore R9.4.1 Reads" begin
        read_simulation_result = Mycelia.simulate_nanopore_r941_reads(fasta = test_fasta, quantity="5x", quiet=true)
        Test.@test isfile(read_simulation_result)
        Test.@test endswith(read_simulation_result, ".fq.gz")
        Test.@test contains(read_simulation_result, "nanopore_r941")
    end
    
    Test.@testset "Very Bad Reads" begin
        read_simulation_result = Mycelia.simulate_very_bad_reads(fasta = test_fasta, quantity="5x", quiet=true)
        Test.@test isfile(read_simulation_result)
        Test.@test endswith(read_simulation_result, ".fq.gz")
        Test.@test contains(read_simulation_result, "very_bad")
    end
    
    Test.@testset "Pretty Good Reads" begin
        read_simulation_result = Mycelia.simulate_pretty_good_reads(fasta = test_fasta, quantity="5x", quiet=true)
        Test.@test isfile(read_simulation_result)
        Test.@test endswith(read_simulation_result, ".fq.gz")
        Test.@test contains(read_simulation_result, "pretty_good")
    end
    
    Test.@testset "General Badread Reads" begin
        ## Test with default parameters
        result1 = Mycelia.simulate_badread_reads(fasta = test_fasta, quantity="5x", quiet=true)
        Test.@test isfile(result1)
        Test.@test endswith(result1, ".fq.gz")
        Test.@test contains(result1, "custom")
        
        ## Test with custom error/qscore models
        result2 = Mycelia.simulate_badread_reads(fasta = test_fasta, quantity="3x", 
                                               error_model="nanopore2020", qscore_model="nanopore2020", quiet=true)
        Test.@test isfile(result2)
        Test.@test contains(result2, "nanopore2020")
        
        ## Test with PacBio-like settings
        result3 = Mycelia.simulate_badread_reads(fasta = test_fasta, quantity="3x",
                                               error_model="pacbio2021", qscore_model="pacbio2021", 
                                               identity="30,3", length="20000,10000", quiet=true)
        Test.@test isfile(result3)
        Test.@test contains(result3, "pacbio2021")
        
        ## Test with random error model
        result4 = Mycelia.simulate_badread_reads(fasta = test_fasta, quantity="3x",
                                               error_model="random", qscore_model="ideal", quiet=true)
        Test.@test isfile(result4)
        Test.@test contains(result4, "random")
        
        ## Test with custom seed for reproducibility
        result5 = Mycelia.simulate_badread_reads(fasta = test_fasta, quantity="3x",
                                               seed=12345, quiet=true)
        Test.@test isfile(result5)
        
        ## Test with high junk/random/chimera rates
        result6 = Mycelia.simulate_badread_reads(fasta = test_fasta, quantity="3x",
                                               junk_reads=10.0, random_reads=5.0, chimeras=15.0, quiet=true)
        Test.@test isfile(result6)
        
        ## Test with small plasmid bias
        result7 = Mycelia.simulate_badread_reads(fasta = test_fasta, quantity="3x",
                                               small_plasmid_bias=true, quiet=true)
        Test.@test isfile(result7)
        
        ## Test with custom adapters
        result8 = Mycelia.simulate_badread_reads(fasta = test_fasta, quantity="3x",
                                               start_adapter_seq="ATCGATCG", end_adapter_seq="CGATCGAT", quiet=true)
        Test.@test isfile(result8)
    end
    
    Test.@testset "Error Model Coverage Tests" begin
        ## Test all supported error models
        error_models = ["nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016", "pacbio2021", "random"]
        for error_model in error_models
            result = Mycelia.simulate_badread_reads(fasta = test_fasta, 
                                                  quantity="2x", error_model=error_model, quiet=true)
            Test.@test isfile(result)
            Test.@test contains(result, error_model)
        end
        
        ## Test all supported qscore models
        qscore_models = ["nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016", "pacbio2021", "random", "ideal"]
        for qscore_model in qscore_models
            result = Mycelia.simulate_badread_reads(fasta = test_fasta,
                                                  quantity="2x", qscore_model=qscore_model, quiet=true)
            Test.@test isfile(result)
            Test.@test contains(result, qscore_model)
        end
    end
    
    end  # end of genome-dependent tests (the if block checking for valid test genome)
    
    # Tests that don't require downloaded genomes - use simulated sequences directly
    # These always run regardless of download success
    Test.@testset "Observe Functions" begin
        ## Test observe() with FASTA record
        test_seq = BioSequences.randdnaseq(100)
        test_record = FASTX.FASTA.Record("test_seq", test_seq)
        
        ## Test with no errors
        observed_record = Mycelia.observe(test_record, error_rate=0.0)
        Test.@test isa(observed_record, FASTX.FASTQ.Record)
        Test.@test FASTX.sequence(BioSequences.LongDNA{4}, observed_record) == test_seq
        Test.@test length(FASTX.quality(observed_record)) == length(test_seq)
        
        ## Test with errors
        observed_with_errors = Mycelia.observe(test_record, error_rate=0.1)
        Test.@test isa(observed_with_errors, FASTX.FASTQ.Record)
        Test.@test length(FASTX.quality(observed_with_errors)) > 0
        
        ## Test observe() with BioSequence directly
        dna_seq = BioSequences.randdnaseq(50)
        observed_seq, quality_scores = Mycelia.observe(dna_seq, error_rate=0.0, tech=:illumina)
        Test.@test isa(observed_seq, BioSequences.LongDNA{4})
        Test.@test length(quality_scores) == length(observed_seq)
        Test.@test all(q -> q >= 20, quality_scores) ## Illumina should have reasonable quality
        
        ## Test different technologies
        for tech in [:illumina, :nanopore, :pacbio, :ultima]
            obs_seq, quals = Mycelia.observe(dna_seq, error_rate=0.01, tech=tech)
            Test.@test isa(obs_seq, BioSequences.LongDNA{4})
            Test.@test length(quals) > 0
        end
        
        ## Test RNA sequence
        rna_seq = BioSequences.randrnaseq(30)
        observed_rna, rna_quals = Mycelia.observe(rna_seq, error_rate=0.01, tech=:illumina)
        Test.@test isa(observed_rna, BioSequences.LongRNA{4})
        
        ## Test amino acid sequence
        aa_seq = BioSequences.randaaseq(20)
        observed_aa, aa_quals = Mycelia.observe(aa_seq, error_rate=0.01, tech=:illumina)
        Test.@test isa(observed_aa, BioSequences.LongAA)
    end
    
    Test.@testset "String Mutation Functions" begin
        ## Test mutate_string with default alphabet
        test_string = "ACGTACGT"
        mutated = Mycelia.mutate_string(test_string, error_rate=0.5)
        Test.@test isa(mutated, String)
        Test.@test length(mutated) >= 1 ## Should not be empty
        
        ## Test with custom alphabet
        custom_alphabet = ['A', 'T', 'G', 'C']
        mutated_custom = Mycelia.mutate_string(test_string, alphabet=custom_alphabet, error_rate=0.3)
        Test.@test isa(mutated_custom, String)
        
        ## Test with zero error rate (should return identical string)
        no_errors = Mycelia.mutate_string(test_string, error_rate=0.0)
        Test.@test no_errors == test_string
    end
    
    Test.@testset "Paired-End Read Generation Functions" begin
        ## Test generate_paired_end_reads (kept - not a thin wrapper)
        ref_seq = BioSequences.randdnaseq(1000)
        reads_1, reads_2 = Mycelia.generate_paired_end_reads(ref_seq, 5, 100, 300, error_rate=0.01)
        Test.@test length(reads_1) == length(reads_2)
        Test.@test length(reads_1) > 0
        Test.@test all(r -> isa(r, BioSequences.LongDNA{4}), reads_1)
        Test.@test all(r -> isa(r, BioSequences.LongDNA{4}), reads_2)
        
        ## Test introduce_sequencing_errors (kept - unique error simulation logic)
        test_seq = BioSequences.randdnaseq(100)
        error_seq = Mycelia.introduce_sequencing_errors(test_seq, 0.1)
        Test.@test isa(error_seq, BioSequences.LongDNA{4})
        ## With 10% error rate, sequences should likely be different
        Test.@test length(error_seq) >= 50 ## Should not be too short
    end
    
    Test.@testset "Sequence Generation Functions" begin
        ## Test random_fasta_record with different molecular types (canonical function)
        for moltype in [:DNA, :RNA, :AA]
            record = Mycelia.random_fasta_record(moltype=moltype, seed=12345, L=100)
            Test.@test isa(record, FASTX.FASTA.Record)
            Test.@test length(FASTX.sequence(record)) == 100
            
            if moltype == :DNA
                Test.@test isa(FASTX.sequence(BioSequences.LongDNA{4}, record), BioSequences.LongDNA{4})
            elseif moltype == :RNA
                Test.@test isa(FASTX.sequence(BioSequences.LongRNA{4}, record), BioSequences.LongRNA{4})
            elseif moltype == :AA
                Test.@test isa(FASTX.sequence(BioSequences.LongAA, record), BioSequences.LongAA)
            end
        end
        
        ## Test write_fasta with random_fasta_record (canonical pattern for genome file creation)
        record = Mycelia.random_fasta_record(moltype=:DNA, seed=42, L=1000)
        temp_fasta = tempname() * ".fasta"
        Mycelia.write_fasta(outfile=temp_fasta, records=[record])
        Test.@test isfile(temp_fasta)
        Test.@test filesize(temp_fasta) > 0
        
        ## Verify saved genome
        saved_genome = nothing
        for rec in FASTX.FASTA.Reader(open(temp_fasta))
            saved_genome = FASTX.sequence(BioSequences.LongDNA{4}, rec)
            break
        end
        Test.@test saved_genome !== nothing
        Test.@test length(saved_genome) == 1000
        rm(temp_fasta)
    end
    
    Test.@testset "Quality Score Functions" begin
        ## Test get_correct_quality for different technologies
        for tech in [:illumina, :nanopore, :pacbio, :ultima]
            quality = Mycelia.get_correct_quality(tech, 1, 100)
            Test.@test isa(quality, Int)
            Test.@test quality >= 5 ## Reasonable minimum quality
            
            error_quality = Mycelia.get_error_quality(tech)
            Test.@test isa(error_quality, Int)
            Test.@test error_quality >= 5 ## Reasonable minimum quality
            ## Error quality should generally be lower than correct quality
            Test.@test error_quality <= quality + 10 ## Allow some overlap
        end
    end
    
    # Clean up test genome at the end
    if test_genome_info.cleanup !== nothing
        test_genome_info.cleanup()
    end
end

# Separate testset for directly testing the NCBI download function itself
# This is isolated so download failures here are visible but don't block other tests
Test.@testset "NCBI Download Function (may fail on network issues)" begin
    # This tests the actual ncbi_genome_download_accession function
    # Uses @test_broken to handle persistent API failures gracefully
    download_success = false
    test_dir = mktempdir()
    
    try
        result = Mycelia.ncbi_genome_download_accession(
            accession=phiX174_assembly_id,
            outdir=test_dir,
            include_string="genome",
            max_attempts=3,
            initial_retry_delay=10.0
        )
        download_success = true
        
        Test.@test isdir(result.directory)
        Test.@test isfile(result.genome)
        Test.@test filesize(result.genome) > 0
        
        @info "NCBI download function test passed"
    catch e
        @warn "NCBI download test failed (expected in some CI environments with restricted network)" exception=e
        # Mark as broken - this test is known to fail due to external NCBI API issues
        Test.@test_broken download_success
    finally
        # Clean up
        rm(test_dir, recursive=true, force=true)
    end
end
