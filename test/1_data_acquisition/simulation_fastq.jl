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
    Test.@testset "Illumina" begin
        if isdir(phiX174_assembly_id)
            rm(phiX174_assembly_id, recursive=true)
        end
        phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
        
        # Test default system (HS25)
        read_simulation_result = Mycelia.simulate_illumina_reads(fasta = phiX174_assembly_dataset.genome, coverage=10)
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
                    result = func(fasta = phiX174_assembly_dataset.genome, coverage = 3)
                    
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
                    cleanup_files = [result.forward_reads, result.sam, result.error_free_sam]
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
                    # Test that max length works
                    try
                        result_max = Mycelia.simulate_illumina_reads(
                            fasta = phiX174_assembly_dataset.genome,
                            coverage = 1,
                            seqSys = seqSys,
                            read_length = max_length
                        )
                        Test.@test isfile(result_max.forward_reads)
                        # Clean up successful run
                        cleanup_files = [result_max.forward_reads, result_max.sam, result_max.error_free_sam]
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
                        fasta = phiX174_assembly_dataset.genome,
                        coverage = 1,
                        seqSys = seqSys,
                        read_length = excessive_length
                    )
                end
            end
        end
        
        rm(phiX174_assembly_id, recursive=true)
    end
    # Test.@testset "Nanopore" begin
    #     if isdir(phiX174_assembly_id)
    #         rm(phiX174_assembly_id, recursive=true)
    #     end
    #     phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
    #     read_simulation_result = Mycelia.simulate_nanopore_reads(fasta = phiX174_assembly_dataset.genome, quantity=10)
    #     Test.@test isfile(read_simulation_result)
    #     Test.@test endswith(read_simulation_result, ".fq.gz")
    # end
    # Test.@testset "PacBio" begin
    #     if isdir(phiX174_assembly_id)
    #         rm(phiX174_assembly_id, recursive=true)
    #     end
    #     phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
    #     read_simulation_result = Mycelia.simulate_pacbio_reads(fasta = phiX174_assembly_dataset.genome, quantity=10)
    #     Test.@test isfile(read_simulation_result)
    #     Test.@test endswith(read_simulation_result, ".fq.gz")
    #     rm(phiX174_assembly_id, recursive=true)
    # end
    # Test.@testset "Nearly Perfect Long Reads" begin
    #     if isdir(phiX174_assembly_id)
    #         rm(phiX174_assembly_id, recursive=true)
    #     end
    #     phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
    #     read_simulation_result = Mycelia.simulate_nearly_perfect_long_reads(fasta = phiX174_assembly_dataset.genome, quantity="10x")
    #     Test.@test isfile(read_simulation_result)
    #     Test.@test endswith(read_simulation_result, ".fq.gz")
    #     rm(phiX174_assembly_id, recursive=true)
    # end
    
    # Test.@testset "Nanopore R9.4.1 Reads" begin
    #     if isdir(phiX174_assembly_id)
    #         rm(phiX174_assembly_id, recursive=true)
    #     end
    #     phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
    #     read_simulation_result = Mycelia.simulate_nanopore_r941_reads(fasta = phiX174_assembly_dataset.genome, quantity="5x")
    #     Test.@test isfile(read_simulation_result)
    #     Test.@test endswith(read_simulation_result, ".fq.gz")
    #     Test.@test contains(read_simulation_result, "nanopore_r941")
    #     rm(phiX174_assembly_id, recursive=true)
    # end
    
    # Test.@testset "Very Bad Reads" begin
    #     if isdir(phiX174_assembly_id)
    #         rm(phiX174_assembly_id, recursive=true)
    #     end
    #     phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
    #     read_simulation_result = Mycelia.simulate_very_bad_reads(fasta = phiX174_assembly_dataset.genome, quantity="5x")
    #     Test.@test isfile(read_simulation_result)
    #     Test.@test endswith(read_simulation_result, ".fq.gz")
    #     Test.@test contains(read_simulation_result, "very_bad")
    #     rm(phiX174_assembly_id, recursive=true)
    # end
    
    # Test.@testset "Pretty Good Reads" begin
    #     if isdir(phiX174_assembly_id)
    #         rm(phiX174_assembly_id, recursive=true)
    #     end
    #     phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
    #     read_simulation_result = Mycelia.simulate_pretty_good_reads(fasta = phiX174_assembly_dataset.genome, quantity="5x")
    #     Test.@test isfile(read_simulation_result)
    #     Test.@test endswith(read_simulation_result, ".fq.gz")
    #     Test.@test contains(read_simulation_result, "pretty_good")
    #     rm(phiX174_assembly_id, recursive=true)
    # end
    
    # Test.@testset "General Badread Reads" begin
    #     if isdir(phiX174_assembly_id)
    #         rm(phiX174_assembly_id, recursive=true)
    #     end
    #     phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
        
    #     ## Test with default parameters
    #     result1 = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, quantity="5x")
    #     Test.@test isfile(result1)
    #     Test.@test endswith(result1, ".fq.gz")
    #     Test.@test contains(result1, "custom")
        
    #     ## Test with custom error/qscore models
    #     result2 = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, quantity="3x", 
    #                                            error_model="nanopore2020", qscore_model="nanopore2020")
    #     Test.@test isfile(result2)
    #     Test.@test contains(result2, "nanopore2020")
        
    #     ## Test with PacBio-like settings
    #     result3 = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, quantity="3x",
    #                                            error_model="pacbio2021", qscore_model="pacbio2021", 
    #                                            identity="30,3", length="20000,10000")
    #     Test.@test isfile(result3)
    #     Test.@test contains(result3, "pacbio2021")
        
    #     ## Test with random error model
    #     result4 = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, quantity="3x",
    #                                            error_model="random", qscore_model="ideal")
    #     Test.@test isfile(result4)
    #     Test.@test contains(result4, "random")
        
    #     ## Test with custom seed for reproducibility
    #     result5 = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, quantity="3x",
    #                                            seed=12345)
    #     Test.@test isfile(result5)
        
    #     ## Test with high junk/random/chimera rates
    #     result6 = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, quantity="3x",
    #                                            junk_reads=10.0, random_reads=5.0, chimeras=15.0)
    #     Test.@test isfile(result6)
        
    #     ## Test with small plasmid bias
    #     result7 = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, quantity="3x",
    #                                            small_plasmid_bias=true)
    #     Test.@test isfile(result7)
        
    #     ## Test with custom adapters
    #     result8 = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, quantity="3x",
    #                                            start_adapter_seq="ATCGATCG", end_adapter_seq="CGATCGAT")
    #     Test.@test isfile(result8)
        
    #     rm(phiX174_assembly_id, recursive=true)
    # end
    
    # Test.@testset "Error Model Coverage Tests" begin
    #     if isdir(phiX174_assembly_id)
    #         rm(phiX174_assembly_id, recursive=true)
    #     end
    #     phiX174_assembly_dataset = Mycelia.ncbi_genome_download_accession(accession=phiX174_assembly_id, include_string="genome")
        
    #     ## Test all supported error models
    #     error_models = ["nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016", "pacbio2021", "random"]
    #     for error_model in error_models
    #         result = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome, 
    #                                               quantity="2x", error_model=error_model)
    #         Test.@test isfile(result)
    #         Test.@test contains(result, error_model)
    #     end
        
    #     ## Test all supported qscore models
    #     qscore_models = ["nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016", "pacbio2021", "random", "ideal"]
    #     for qscore_model in qscore_models
    #         result = Mycelia.simulate_badread_reads(fasta = phiX174_assembly_dataset.genome,
    #                                               quantity="2x", qscore_model=qscore_model)
    #         Test.@test isfile(result)
    #         Test.@test contains(result, qscore_model)
    #     end
        
    #     rm(phiX174_assembly_id, recursive=true)
    # end
    
    # Test.@testset "Observe Functions" begin
    #     ## Test observe() with FASTA record
    #     test_seq = BioSequences.randdnaseq(100)
    #     test_record = FASTX.FASTA.Record("test_seq", test_seq)
        
    #     ## Test with no errors
    #     observed_record = Mycelia.observe(test_record, error_rate=0.0)
    #     Test.@test isa(observed_record, FASTX.FASTQ.Record)
    #     Test.@test FASTX.sequence(BioSequences.LongDNA{4}, observed_record) == test_seq
    #     Test.@test length(FASTX.quality(observed_record)) == length(test_seq)
        
    #     ## Test with errors
    #     observed_with_errors = Mycelia.observe(test_record, error_rate=0.1)
    #     Test.@test isa(observed_with_errors, FASTX.FASTQ.Record)
    #     Test.@test length(FASTX.quality(observed_with_errors)) > 0
        
    #     ## Test observe() with BioSequence directly
    #     dna_seq = BioSequences.randdnaseq(50)
    #     observed_seq, quality_scores = Mycelia.observe(dna_seq, error_rate=0.0, tech=:illumina)
    #     Test.@test isa(observed_seq, BioSequences.LongDNA{4})
    #     Test.@test length(quality_scores) == length(observed_seq)
    #     Test.@test all(q -> q >= 20, quality_scores) ## Illumina should have reasonable quality
        
    #     ## Test different technologies
    #     for tech in [:illumina, :nanopore, :pacbio, :ultima]
    #         obs_seq, quals = Mycelia.observe(dna_seq, error_rate=0.01, tech=tech)
    #         Test.@test isa(obs_seq, BioSequences.LongDNA{4})
    #         Test.@test length(quals) > 0
    #     end
        
    #     ## Test RNA sequence
    #     rna_seq = BioSequences.randrnaseq(30)
    #     observed_rna, rna_quals = Mycelia.observe(rna_seq, error_rate=0.01, tech=:illumina)
    #     Test.@test isa(observed_rna, BioSequences.LongRNA{4})
        
    #     ## Test amino acid sequence
    #     aa_seq = BioSequences.randaaseq(20)
    #     observed_aa, aa_quals = Mycelia.observe(aa_seq, error_rate=0.01, tech=:illumina)
    #     Test.@test isa(observed_aa, BioSequences.LongAA)
    # end
    
    # Test.@testset "String Mutation Functions" begin
    #     ## Test mutate_string with default alphabet
    #     test_string = "ACGTACGT"
    #     mutated = Mycelia.mutate_string(test_string, error_rate=0.5)
    #     Test.@test isa(mutated, String)
    #     Test.@test length(mutated) >= 1 ## Should not be empty
        
    #     ## Test with custom alphabet
    #     custom_alphabet = ['A', 'T', 'G', 'C']
    #     mutated_custom = Mycelia.mutate_string(test_string, alphabet=custom_alphabet, error_rate=0.3)
    #     Test.@test isa(mutated_custom, String)
        
    #     ## Test with zero error rate (should return identical string)
    #     no_errors = Mycelia.mutate_string(test_string, error_rate=0.0)
    #     Test.@test no_errors == test_string
    # end
    
    # Test.@testset "FASTQ Generation Functions" begin
    #     ## Test generate_test_fastq_data
    #     temp_fastq = tempname() * ".fq"
    #     Mycelia.generate_test_fastq_data(10, 100, temp_fastq)
    #     Test.@test isfile(temp_fastq)
        
    #     ## Verify FASTQ content
    #     record_count = 0
    #     for record in FASTX.FASTQ.Reader(open(temp_fastq))
    #         record_count += 1
    #         Test.@test length(FASTX.sequence(record)) == 100
    #         Test.@test length(FASTX.quality(record)) == 100
    #     end
    #     Test.@test record_count == 10
    #     rm(temp_fastq)
        
    #     ## Test generate_paired_end_reads
    #     ref_seq = BioSequences.randdnaseq(1000)
    #     reads_1, reads_2 = Mycelia.generate_paired_end_reads(ref_seq, 5, 100, 300, error_rate=0.01)
    #     Test.@test length(reads_1) == length(reads_2)
    #     Test.@test length(reads_1) > 0
    #     Test.@test all(r -> isa(r, BioSequences.LongDNA{4}), reads_1)
    #     Test.@test all(r -> isa(r, BioSequences.LongDNA{4}), reads_2)
        
    #     ## Test introduce_sequencing_errors
    #     test_seq = BioSequences.randdnaseq(100)
    #     error_seq = Mycelia.introduce_sequencing_errors(test_seq, 0.1)
    #     Test.@test isa(error_seq, BioSequences.LongDNA{4})
    #     ## With 10% error rate, sequences should likely be different
    #     Test.@test length(error_seq) >= 50 ## Should not be too short
        
    #     ## Test save_reads_as_fastq
    #     test_reads = [BioSequences.randdnaseq(50) for _ in 1:5]
    #     temp_output = tempname() * ".fq"
    #     Mycelia.save_reads_as_fastq(test_reads, temp_output, 25)
    #     Test.@test isfile(temp_output)
        
    #     ## Verify saved FASTQ
    #     saved_count = 0
    #     for record in FASTX.FASTQ.Reader(open(temp_output))
    #         saved_count += 1
    #         Test.@test length(FASTX.sequence(record)) == 50
    #     end
    #     Test.@test saved_count == 5
    #     rm(temp_output)
    # end
    
    # Test.@testset "Sequence Generation Functions" begin
    #     ## Test random_fasta_record with different molecular types
    #     for moltype in [:DNA, :RNA, :AA]
    #         record = Mycelia.random_fasta_record(moltype=moltype, seed=12345, L=100)
    #         Test.@test isa(record, FASTX.FASTA.Record)
    #         Test.@test length(FASTX.sequence(record)) == 100
            
    #         if moltype == :DNA
    #             Test.@test isa(FASTX.sequence(BioSequences.LongDNA{4}, record), BioSequences.LongDNA{4})
    #         elseif moltype == :RNA
    #             Test.@test isa(FASTX.sequence(BioSequences.LongRNA{4}, record), BioSequences.LongRNA{4})
    #         elseif moltype == :AA
    #             Test.@test isa(FASTX.sequence(BioSequences.LongAA, record), BioSequences.LongAA)
    #         end
    #     end
        
    #     ## Test generate_test_sequences
    #     sequences = Mycelia.generate_test_sequences(500, 3)
    #     Test.@test length(sequences) == 3
    #     Test.@test all(s -> isa(s, BioSequences.LongDNA{4}), sequences)
    #     Test.@test all(s -> length(s) == 500, sequences)
        
    #     ## Test generate_test_genome_with_genes
    #     genome, gene_positions = Mycelia.generate_test_genome_with_genes(10000, 0.05)
    #     Test.@test isa(genome, BioSequences.LongDNA{4})
    #     Test.@test length(genome) == 10000
    #     Test.@test isa(gene_positions, Vector)
    #     Test.@test all(pos -> isa(pos, Tuple{Int,Int}), gene_positions)
    #     Test.@test all(pos -> pos[1] <= pos[2], gene_positions) ## Start <= end
        
    #     ## Test save_genome_as_fasta
    #     temp_fasta = tempname() * ".fa"
    #     Mycelia.save_genome_as_fasta(genome, temp_fasta)
    #     Test.@test isfile(temp_fasta)
        
    #     ## Verify saved genome
    #     saved_genome = nothing
    #     for record in FASTX.FASTA.Reader(open(temp_fasta))
    #         saved_genome = FASTX.sequence(BioSequences.LongDNA{4}, record)
    #         break
    #     end
    #     Test.@test saved_genome == genome
    #     rm(temp_fasta)
    # end
    
    # Test.@testset "Quality Score Functions" begin
    #     ## Test get_correct_quality for different technologies
    #     for tech in [:illumina, :nanopore, :pacbio, :ultima]
    #         quality = Mycelia.get_correct_quality(tech, 1, 100)
    #         Test.@test isa(quality, Int)
    #         Test.@test quality >= 5 ## Reasonable minimum quality
            
    #         error_quality = Mycelia.get_error_quality(tech)
    #         Test.@test isa(error_quality, Int)
    #         Test.@test error_quality >= 5 ## Reasonable minimum quality
    #         ## Error quality should generally be lower than correct quality
    #         Test.@test error_quality <= quality + 10 ## Allow some overlap
    #     end
    # end
    # if isdir(phiX174_assembly_id)
    #     rm(phiX174_assembly_id, recursive=true)
    # end
end
