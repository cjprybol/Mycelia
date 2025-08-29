# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/2_preprocessing_qc/preprocessing.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/2_preprocessing_qc/preprocessing.jl", "test/2_preprocessing_qc", execute=false)'
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
import CSV
import CodecZlib
import DataFrames
import Statistics

const phiX174_accession_id = "NC_001422.1"

Test.@testset "Preprocessing" begin
    Test.@testset "FASTX stats" begin
        srr_identifier = "SRR31812976"
        outdir = mkpath("fastx-stats-test")
        fasterq_dump_result = Mycelia.fasterq_dump(outdir=outdir, srr_identifier=srr_identifier)
        Test.@test fasterq_dump_result.unpaired_reads == "$(outdir)/$(srr_identifier)/$(srr_identifier).fastq.gz"
        table = Mycelia.fastx_stats(fasterq_dump_result.unpaired_reads)
        io = IOBuffer()
        CSV.write(io, table)
        ## The exact values can vary, so test the structure and types
        csv_string = String(take!(io))
        lines = split(strip(csv_string), '\n')
        Test.@test length(lines) == 2  ## header + data row
        
        ## Check header
        header = lines[1]
        Test.@test startswith(header, "file,format,type,num_seqs,sum_len")
        
        ## Check data values are reasonable
        data = lines[2]
        fields = split(data, ',')
        Test.@test fields[1] == "fastx-stats-test/SRR31812976/SRR31812976.fastq.gz"
        Test.@test fields[2] == "FASTQ"
        Test.@test fields[3] == "DNA"
        Test.@test parse(Int, fields[4]) > 50000  ## reasonable number of sequences
        rm(outdir, recursive=true)
    end
    Test.@testset "fastx2normalized_table" begin
        genome_result = Mycelia.download_genome_by_accession(accession=phiX174_accession_id)
        fastx_table = Mycelia.fastx2normalized_table(genome_result)
        
        ## Test that we get expected columns
        expected_cols = [
            "fastx_path", "human_readable_id", "genome_hash", "sequence_hash",
            "genome_identifier", "sequence_identifier", "record_identifier",
            "record_description", "record_length", "record_alphabet", "record_type",
            "mean_record_quality", "median_record_quality", "record_quality",
            "record_sequence"
        ]
        Test.@test names(fastx_table) == expected_cols
        
        ## Test basic properties
        Test.@test DataFrames.nrow(fastx_table) == 1  ## phiX174 is single sequence
        Test.@test fastx_table.record_length[1] == 5386  ## known phiX174 length
        Test.@test fastx_table.record_alphabet[1] == "ACGT"  ## DNA alphabet
        Test.@test fastx_table.record_type[1] == "DNA"
        Test.@test ismissing(fastx_table.record_quality[1])  ## FASTA has no quality
        Test.@test ismissing(fastx_table.mean_record_quality[1])
        Test.@test ismissing(fastx_table.median_record_quality[1])
        
        rm(genome_result)
    end
    
    Test.@testset "normalized_table2fastx" begin
        mktempdir() do dir
            ## Create test data - both FASTA and FASTQ
            fasta_path = joinpath(dir, "test.fasta")  
            fastq_path = joinpath(dir, "test.fastq")
            
            ## Create proper FASTA file
            open(fasta_path, "w") do io
                println(io, ">seq1")
                println(io, "ACGTACGT")
                println(io, ">seq2") 
                println(io, "GGGGCCCC")
            end
            
            ## Create proper FASTQ file  
            open(fastq_path, "w") do io
                println(io, "@seq1")
                println(io, "ACGTACGT")
                println(io, "+")
                println(io, "IIIIIIII")
                println(io, "@seq2")
                println(io, "GGGGCCCC") 
                println(io, "+")
                println(io, "HHHHHHHH")
            end
            
            ## Test FASTA round-trip
            fasta_table = Mycelia.fastx2normalized_table(fasta_path)
            out_fasta = Mycelia.normalized_table2fastx(fasta_table; output_dir=dir, output_basename="test_out")
            Test.@test isfile(out_fasta)
            Test.@test endswith(out_fasta, ".fna")
            
            ## Verify content by reading back
            fasta_roundtrip = Mycelia.fastx2normalized_table(out_fasta)  
            Test.@test DataFrames.nrow(fasta_roundtrip) == 2
            Test.@test Set(fasta_roundtrip.record_sequence) == Set(fasta_table.record_sequence)
            
            ## Test FASTQ round-trip
            fastq_table = Mycelia.fastx2normalized_table(fastq_path)
            out_fastq = Mycelia.normalized_table2fastx(fastq_table; output_dir=dir, output_basename="test_out_fq")
            Test.@test isfile(out_fastq)
            Test.@test endswith(out_fastq, ".fq")
            
            ## Verify content by reading back
            fastq_roundtrip = Mycelia.fastx2normalized_table(out_fastq)
            Test.@test DataFrames.nrow(fastq_roundtrip) == 2  
            Test.@test Set(fastq_roundtrip.record_sequence) == Set(fastq_table.record_sequence)
            Test.@test !any(ismissing, fastq_roundtrip.mean_record_quality)  ## FASTQ should have quality
            
            ## Test gzip compression
            out_fasta_gz = Mycelia.normalized_table2fastx(fasta_table; output_dir=dir, output_basename="test_gz", gzip=true)
            Test.@test isfile(out_fasta_gz)
            Test.@test endswith(out_fasta_gz, ".fna.gz")
        end
    end
    Test.@testset "Read Quality Control" begin
        mktempdir() do dir
            ## Create test FASTQ files for testing preprocessing functions
            fastq_single = joinpath(dir, "test_single.fastq.gz")
            fastq_R1 = joinpath(dir, "test_R1.fastq.gz")
            fastq_R2 = joinpath(dir, "test_R2.fastq.gz")
            
            ## Create a simple test FASTQ content (single-end)
            fastq_content = """@seq1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@seq2  
GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
"""
            
            ## Write gzipped test files
            open(fastq_single, "w") do f
                stream = CodecZlib.GzipCompressorStream(f)
                write(stream, fastq_content)
                close(stream)
            end
            open(fastq_R1, "w") do f
                stream = CodecZlib.GzipCompressorStream(f)
                write(stream, fastq_content)  
                close(stream)
            end
            open(fastq_R2, "w") do f
                stream = CodecZlib.GzipCompressorStream(f)
                write(stream, fastq_content)
                close(stream)
            end
            
            ## Test fastp (short reads)
            Test.@testset "fastp" begin
                result = Mycelia.qc_filter_short_reads_fastp(
                    forward_reads = fastq_R1,
                    reverse_reads = fastq_R2
                )
                Test.@test isfile(result.out_forward)
                Test.@test isfile(result.out_reverse) 
                Test.@test isfile(result.json)
                Test.@test isfile(result.html)
            end
            
            ## Test trim_galore (paired-end)
            Test.@testset "trim_galore" begin
                result = Mycelia.trim_galore_paired(
                    forward_reads = fastq_R1,
                    reverse_reads = fastq_R2,
                    outdir = dir
                )
                Test.@test isfile(result.trimmed_forward)
                Test.@test isfile(result.trimmed_reverse)
            end
            
            ## Test filtlong (long reads)
            Test.@testset "filtlong" begin
                out_fastq = Mycelia.qc_filter_long_reads_filtlong(
                    in_fastq = fastq_single
                )
                Test.@test isfile(out_fastq)
                Test.@test endswith(out_fastq, ".filtlong.fq.gz")
            end
            
            ## Test fastplong (long reads)
            Test.@testset "fastplong" begin
                out_fastq = Mycelia.qc_filter_long_reads_fastplong(
                    in_fastq = fastq_single
                )
                Test.@test isfile(out_fastq)
                Test.@test endswith(out_fastq, ".fastplong.fq.gz")
            end
            
            ## Test chopper (long reads)
            Test.@testset "chopper" begin
                out_fastq = Mycelia.qc_filter_long_reads_chopper(
                    in_fastq = fastq_single
                )
                Test.@test isfile(out_fastq)
                Test.@test endswith(out_fastq, ".chopper.fq.gz")
            end
        end
    end
    Test.@testset "Read Statistics" begin
        mktempdir() do dir
            ## Create test FASTQ file for statistics
            fastq_path = joinpath(dir, "stats_test.fastq")
            
            ## Create FASTQ with known properties  
            open(fastq_path, "w") do io
                ## First read: ACGT (length=4, high quality)
                println(io, "@read1")
                println(io, "ACGT")
                println(io, "+")
                println(io, "IIII")  ## Quality 40 for each base
                
                ## Second read: GGCC (length=4, lower quality)
                println(io, "@read2")
                println(io, "GGCC")
                println(io, "+")
                println(io, "HHHH")  ## Quality 39 for each base
                
                ## Third read: longer sequence
                println(io, "@read3")
                println(io, "ACGTACGTACGT")  ## length=12
                println(io, "+")
                println(io, "IIIIIIIIIIII")
            end
            
            ## Test fastx_stats function
            Test.@testset "fastx_stats" begin
                stats = Mycelia.fastx_stats(fastq_path)
                Test.@test DataFrames.nrow(stats) == 1  ## One row for the file
                Test.@test stats.num_seqs[1] == 3  ## 3 sequences
                Test.@test stats.sum_len[1] == 20  ## 4+4+12 = 20 total bases
                Test.@test stats.min_len[1] == 4  ## shortest read
                Test.@test stats.max_len[1] == 12  ## longest read
                Test.@test stats.format[1] == "FASTQ"
                Test.@test stats.type[1] == "DNA"
            end
            
            ## Test assess_duplication_rates function
            Test.@testset "assess_duplication_rates" begin
                results_file = Mycelia.assess_duplication_rates(fastq_path)
                Test.@test isfile(results_file)
                Test.@test endswith(results_file, ".duplication_rates.tsv")
                
                ## Read and check the results
                results_df = CSV.read(results_file, DataFrames.DataFrame, delim='\t')
                Test.@test results_df.total_records[1] == 3
                Test.@test results_df.total_unique_observations[1] == 3  ## All sequences unique
                Test.@test results_df.percent_duplication_rate[1] ≈ 0.0  ## No duplicates
            end
            
            ## Test count_records function
            Test.@testset "count_records" begin
                count = Mycelia.count_records(fastq_path)
                Test.@test count == 3
            end
            
            ## Test determine_read_lengths function
            Test.@testset "determine_read_lengths" begin
                lengths = Mycelia.determine_read_lengths(fastq_path)
                Test.@test length(lengths) == 3
                Test.@test Set(lengths) == Set([4, 4, 12])  ## Expected lengths
            end
            
            ## Test quality conversion functions
            Test.@testset "quality_conversions" begin
                ## Test q-value to error rate conversion
                Test.@test Mycelia.q_value_to_error_rate(10) ≈ 0.1
                Test.@test Mycelia.q_value_to_error_rate(20) ≈ 0.01
                Test.@test Mycelia.q_value_to_error_rate(30) ≈ 0.001
                
                ## Test error rate to q-value conversion  
                Test.@test Mycelia.error_rate_to_q_value(0.1) ≈ 10.0
                Test.@test Mycelia.error_rate_to_q_value(0.01) ≈ 20.0
                Test.@test Mycelia.error_rate_to_q_value(0.001) ≈ 30.0
            end
        end
    end
    Test.@testset "fastx2normalized_table - compressed and raw input" begin
        mktempdir() do dir
            # Prepare test records
            fasta_file = joinpath(dir, "test-normalized.fasta")
            fastq_file = joinpath(dir, "test-normalized.fastq")
            gz_fasta = fasta_file * ".gz"
            gz_fastq = fastq_file * ".gz"
            
            # Write FASTA with simple sequences
            open(fasta_file, "w") do io
                println(io, ">seq1")
                println(io, "ACGTACGTAC")
                println(io, ">seq2")
                println(io, "GGGGCCCCAA")
            end
            
            # Write FASTQ with same sequences + quality
            open(fastq_file, "w") do io
                println(io, "@seq1")
                println(io, "ACGTACGTAC")
                println(io, "+")
                println(io, "IIIIIIIIII")
                println(io, "@seq2")
                println(io, "GGGGCCCCAA")
                println(io, "+")
                println(io, "HHHHHHHHHH")
            end
            
            # Compress files
            run(pipeline(`gzip -c $fasta_file`, stdout=gz_fasta))
            run(pipeline(`gzip -c $fastq_file`, stdout=gz_fastq))

            # Test raw and gzipped FASTA
            for file in (fasta_file, gz_fasta)
                table = Mycelia.fastx2normalized_table(file)
                Test.@test DataFrames.nrow(table) == 2
                # Case-insensitive: sequences should match after uppercasing
                seqs = [uppercase(row.record_sequence) for row in DataFrames.eachrow(table)]
                Test.@test length(unique(seqs)) == 2
            end

            # Test raw and gzipped FASTQ
            for file in (fastq_file, gz_fastq)
                table = Mycelia.fastx2normalized_table(file)
                Test.@test DataFrames.nrow(table) == 2
                # Case-sensitive: sequences should match exactly
                seqs = [row.record_sequence for row in DataFrames.eachrow(table)]
                Test.@test length(unique(seqs)) == 2
            end
        end
    end

    Test.@testset "fastx2normalized_table - identifier normalization" begin
        mktempdir() do dir
            # Two records with different identifiers but same sequence (case-insensitive)
            fasta_file = joinpath(dir, "test-idnorm.fasta")
            open(fasta_file, "w") do io
                println(io, ">id1")
                println(io, "acgtACGT")
                println(io, ">id2")
                println(io, "ACGTacgt")
            end
            table = Mycelia.fastx2normalized_table(fasta_file)
            seqs = [uppercase(row.record_sequence) for row in DataFrames.eachrow(table)]
            Test.@test length(unique(seqs)) == 1
        end
    end

    Test.@testset "human_readable_id automatic assignment" begin
        mktempdir() do dir
            ## Test cases for automatic human readable ID extraction
            test_cases = [
                # (filename, expected_id, description)
                ("simple.fasta", "simple", "Simple filename"),
                ("NC_001422.1.fna", "NC_001422.1", "GenBank accession with version"),
                ("GCF_000005825.2_ASM582v2_genomic.fna", "GCF_000005825.2", "RefSeq assembly (avoid meaningless prefix)"),
                ("sample_data_v2.fastq", "sample_data_v2", "Descriptive name with underscores"),
                ("long-file-name.fq", "long-file-name", "Hyphenated name"),
                ("my_genome_assembly_v1_final.fasta", "my_genome", "Long name with meaningful truncation"),
                ("E_coli_K12.fasta", "E_coli_K12", "Species name with underscore"),
            ]
            
            for (filename, expected_id, description) in test_cases
                Test.@testset "$(description)" begin
                    filepath = joinpath(dir, filename)
                    
                    ## Create appropriate file format based on extension
                    open(filepath, "w") do io
                        if endswith(filename, ".fastq") || endswith(filename, ".fq")
                            ## Create FASTQ format
                            println(io, "@test_seq")
                            println(io, "ACGTACGT")
                            println(io, "+")
                            println(io, "IIIIIIII")
                        else
                            ## Create FASTA format  
                            println(io, ">test_seq")
                            println(io, "ACGTACGT")
                        end
                    end
                    
                    ## Test automatic ID extraction
                    table = Mycelia.fastx2normalized_table(filepath)
                    Test.@test table.human_readable_id[1] == expected_id
                    Test.@test startswith(table.genome_identifier[1], expected_id * "_")
                    Test.@test startswith(table.sequence_identifier[1], expected_id * "_")
                end
            end
            
            ## Test explicit human_readable_id override
            Test.@testset "explicit human_readable_id override" begin
                filepath = joinpath(dir, "auto_name.fasta")
                open(filepath, "w") do io
                    println(io, ">test_seq")
                    println(io, "ACGTACGT")
                end
                
                table = Mycelia.fastx2normalized_table(filepath; human_readable_id="custom_id")
                Test.@test table.human_readable_id[1] == "custom_id"
                Test.@test startswith(table.genome_identifier[1], "custom_id_")
            end
            
            ## Test force_truncate functionality
            Test.@testset "force_truncate functionality" begin
                ## Create a filename that has no meaningful delimiters and exceeds 16 chars
                filepath = joinpath(dir, "verylongfilenamewithoutdelimitersthatexceedssixteencharacters.fasta")
                open(filepath, "w") do io
                    println(io, ">test_seq")
                    println(io, "ACGTACGT")
                end
                
                ## Should work with force_truncate=true
                table = Mycelia.fastx2normalized_table(filepath; force_truncate=true)
                Test.@test length(table.human_readable_id[1]) <= 16
                
                ## Should error without force_truncate (test this by catching the error)
                Test.@test_throws ErrorException Mycelia.fastx2normalized_table(filepath; force_truncate=false)
            end
            
            ## Test length validation
            Test.@testset "human_readable_id length validation" begin
                filepath = joinpath(dir, "test.fasta")
                open(filepath, "w") do io
                    println(io, ">test_seq")
                    println(io, "ACGTACGT")
                end
                
                ## Should error for ID > 16 chars without force_truncate
                Test.@test_throws ErrorException Mycelia.fastx2normalized_table(
                    filepath; 
                    human_readable_id="this_id_is_way_too_long_and_exceeds_sixteen_chars"
                )
                
                ## Should work with force_truncate
                table = Mycelia.fastx2normalized_table(
                    filepath; 
                    human_readable_id="this_id_is_way_too_long", 
                    force_truncate=true
                )
                Test.@test length(table.human_readable_id[1]) == 16
                Test.@test table.human_readable_id[1] == "this_id_is_way_t"
            end
            
            ## Test identifier propagation through the hierarchy
            Test.@testset "identifier hierarchy propagation" begin
                filepath = joinpath(dir, "hierarchy_test.fasta")
                open(filepath, "w") do io
                    println(io, ">seq1")
                    println(io, "ACGTACGT")
                    println(io, ">seq2")
                    println(io, "GGGGCCCC")
                end
                
                table = Mycelia.fastx2normalized_table(filepath; human_readable_id="test_genome")
                
                ## Check that all identifiers have proper hierarchy
                Test.@test all(table.human_readable_id .== "test_genome")
                Test.@test all(length.(table.genome_hash) .== 16)  ## Base58 encoded, 16 chars
                Test.@test all(length.(table.sequence_hash) .== 16)  ## Base58 encoded, 16 chars
                
                ## genome_identifier = human_readable_id + "_" + genome_hash
                for row in DataFrames.eachrow(table)
                    Test.@test row.genome_identifier == "test_genome_" * row.genome_hash
                    Test.@test row.sequence_identifier == row.genome_identifier * "_" * row.sequence_hash
                    Test.@test length(row.sequence_identifier) <= 50  ## NCBI limit
                end
                
                ## All sequences should have the same genome_identifier but different sequence_identifier
                Test.@test length(unique(table.genome_identifier)) == 1
                Test.@test length(unique(table.sequence_identifier)) == 2
            end
        end
    end

    Test.@testset "fastx_stats and fastx2normalized_table - fixtures" begin
        mktempdir() do dir
            fasta_path = joinpath(dir, "small.fasta")
            fastq_path = joinpath(dir, "small.fastq")
            
            ## Create proper FASTA file
            open(fasta_path, "w") do io
                println(io, ">seq1")
                println(io, "ACGT")
                println(io, ">seq2")
                println(io, "GGGG")
            end
            
            ## Create proper FASTQ file
            open(fastq_path, "w") do io
                println(io, "@seq1")
                println(io, "ACGT")
                println(io, "+")
                println(io, "IIII")
                println(io, "@seq2")
                println(io, "GGGG")
                println(io, "+")
                println(io, "IIII")
            end

            expected_cols = [
                "file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len",
                "max_len", "Q1", "Q2", "Q3", "sum_gap", "N50", "N50_num", "Q20(%)",
                "Q30(%)", "AvgQual", "GC(%)", "sum_n", "N90",
            ]

            stats_fasta = Mycelia.fastx_stats(fasta_path)
            Test.@test names(stats_fasta) == expected_cols
            Test.@test stats_fasta.num_seqs[1] == 2
            Test.@test stats_fasta.sum_len[1] == 8
            Test.@test isapprox(stats_fasta[1, "GC(%)"], 75.0; atol=0.01)

            stats_fastq = Mycelia.fastx_stats(fastq_path)
            Test.@test stats_fastq.num_seqs[1] == 2
            Test.@test isapprox(stats_fastq.AvgQual[1], 40.0; atol=0.01)

            ## Updated column names for fastx2normalized_table
            norm_cols = [
                "fastx_path", "human_readable_id", "genome_hash", "sequence_hash",
                "genome_identifier", "sequence_identifier", "record_identifier",
                "record_description", "record_length", "record_alphabet", "record_type",
                "mean_record_quality", "median_record_quality", "record_quality",
                "record_sequence"
            ]

            nfasta = Mycelia.fastx2normalized_table(fasta_path)
            nfastq = Mycelia.fastx2normalized_table(fastq_path)
            Test.@test names(nfasta) == norm_cols
            Test.@test names(nfastq) == norm_cols
            Test.@test DataFrames.nrow(nfasta) == 2
            gc_val = sum(count(c -> c in ['G','C'], row.record_sequence) for row in DataFrames.eachrow(nfasta)) / sum(nfasta.record_length) * 100
            Test.@test isapprox(gc_val, 75.0; atol=0.01)
            Test.@test isapprox(Statistics.mean(skipmissing(nfastq.mean_record_quality)), 40.0; atol=0.01)
        end
    end
end
