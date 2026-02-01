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
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import FASTX
import CSV
import CodecZlib
import DataFrames
import Statistics
import JSON

const phiX174_accession_id = "NC_001422.1"

Test.@testset "Preprocessing" begin
    Test.@testset "FASTX stats" begin
        srr_identifier = "SRR31812976"
        outdir = mkpath("fastx-stats-test")
        fasterq_dump_result = Mycelia.fasterq_dump(outdir = outdir, srr_identifier = srr_identifier)
        Test.@test fasterq_dump_result.unpaired_reads ==
                   "$(outdir)/$(srr_identifier)/$(srr_identifier).fastq.gz"
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
        rm(outdir, recursive = true)
    end
    Test.@testset "fastx2normalized_table" begin
        genome_result = Mycelia.download_genome_by_accession(accession = phiX174_accession_id)
        fastx_table = Mycelia.fastx2normalized_table(genome_result)

        ## Test that we get expected columns
        expected_cols = [
            "fastx_path", "human_readable_id", "fastx_hash", "sequence_hash",
            "fastx_identifier", "sequence_identifier", "record_identifier",
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
            out_fasta = Mycelia.normalized_table2fastx(
                fasta_table; output_dir = dir, output_basename = "test_out")
            Test.@test isfile(out_fasta)
            Test.@test endswith(out_fasta, ".fna")

            ## Verify content by reading back
            fasta_roundtrip = Mycelia.fastx2normalized_table(out_fasta)
            Test.@test DataFrames.nrow(fasta_roundtrip) == 2
            Test.@test Set(fasta_roundtrip.record_sequence) ==
                       Set(fasta_table.record_sequence)

            ## Test FASTQ round-trip
            fastq_table = Mycelia.fastx2normalized_table(fastq_path)
            out_fastq = Mycelia.normalized_table2fastx(
                fastq_table; output_dir = dir, output_basename = "test_out_fq")
            Test.@test isfile(out_fastq)
            Test.@test endswith(out_fastq, ".fq")

            ## Verify content by reading back
            fastq_roundtrip = Mycelia.fastx2normalized_table(out_fastq)
            Test.@test DataFrames.nrow(fastq_roundtrip) == 2
            Test.@test Set(fastq_roundtrip.record_sequence) ==
                       Set(fastq_table.record_sequence)
            Test.@test !any(ismissing, fastq_roundtrip.mean_record_quality)  ## FASTQ should have quality

            ## Test gzip compression
            out_fasta_gz = Mycelia.normalized_table2fastx(
                fasta_table; output_dir = dir, output_basename = "test_gz", gzip = true)
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
                ## Test default behavior (automatic dedup logic)
                Test.@testset "fastp default (auto dedup)" begin
                    result = Mycelia.qc_filter_short_reads_fastp(
                        forward_reads = fastq_R1,
                        reverse_reads = fastq_R2
                    )
                    Test.@test isfile(result.out_forward)
                    Test.@test isfile(result.out_reverse)
                    Test.@test isfile(result.json)
                    Test.@test isfile(result.html)
                end

                ## Test with explicit deduplication disabled
                Test.@testset "fastp dedup disabled" begin
                    result = Mycelia.qc_filter_short_reads_fastp(
                        forward_reads = fastq_R1,
                        reverse_reads = fastq_R2,
                        enable_dedup = false
                    )
                    Test.@test isfile(result.out_forward)
                    Test.@test isfile(result.out_reverse)
                    Test.@test isfile(result.json)
                    Test.@test isfile(result.html)
                end

                ## Test with explicit deduplication enabled
                Test.@testset "fastp dedup enabled" begin
                    result = Mycelia.qc_filter_short_reads_fastp(
                        forward_reads = fastq_R1,
                        reverse_reads = fastq_R2,
                        enable_dedup = true
                    )
                    Test.@test isfile(result.out_forward)
                    Test.@test isfile(result.out_reverse)
                    Test.@test isfile(result.json)
                    Test.@test isfile(result.html)
                end

                ## Test deduplication functionality with actual duplicated sequences
                Test.@testset "fastp deduplication validation" begin
                    ## Create unique temporary directories for each test to avoid conflicts
                    dup_test_dir = mktempdir(dir)

                    ## Create test files with duplicated sequences using Mycelia.write_fastq
                    dup_fastq_R1 = joinpath(dup_test_dir, "test_dup_R1.fastq.gz")
                    dup_fastq_R2 = joinpath(dup_test_dir, "test_dup_R2.fastq.gz")

                    try
                        ## Create FASTQ records with duplicates:
                        ## - Same sequence, different identifiers (should be deduplicated)
                        ## - Same sequence, same identifier (should be deduplicated)  
                        ## - Unique sequences (should be retained)
                        duplicate_records = [
                            ## First occurrence of duplicate sequence
                            FASTX.FASTQ.Record("seq1",
                                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
                            ## Same sequence, different identifier (duplicate)
                            FASTX.FASTQ.Record("seq2_different_id",
                                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
                            ## Same sequence, same identifier (duplicate)
                            FASTX.FASTQ.Record("seq1",
                                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
                            ## Unique sequence 1
                            FASTX.FASTQ.Record("seq3_unique",
                                "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
                                "JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"),
                            ## Unique sequence 2
                            FASTX.FASTQ.Record("seq4_another_unique",
                                "GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC",
                                "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK")
                        ]

                        ## Write duplicate test files using Mycelia.write_fastq
                        Mycelia.write_fastq(records = duplicate_records, filename = dup_fastq_R1)
                        Mycelia.write_fastq(records = duplicate_records, filename = dup_fastq_R2)

                        ## Test with deduplication disabled - should keep all sequences
                        Test.@testset "dedup disabled - keeps duplicates" begin
                            ## Use unique output paths to avoid caching
                            no_dedup_out_forward = joinpath(dup_test_dir, "test_no_dedup_R1.fastp.1.fq.gz")
                            no_dedup_out_reverse = joinpath(dup_test_dir, "test_no_dedup_R2.fastp.2.fq.gz")
                            no_dedup_json = joinpath(dup_test_dir, "test_no_dedup_R.fastp_report.json")
                            no_dedup_html = joinpath(dup_test_dir, "test_no_dedup_R.fastp_report.html")

                            result_no_dedup = nothing
                            reader = nothing
                            try
                                result_no_dedup = Mycelia.qc_filter_short_reads_fastp(
                                    forward_reads = dup_fastq_R1,
                                    reverse_reads = dup_fastq_R2,
                                    out_forward = no_dedup_out_forward,
                                    out_reverse = no_dedup_out_reverse,
                                    json = no_dedup_json,
                                    html = no_dedup_html,
                                    enable_dedup = false
                                )

                                ## Count sequences in output (should have all 5 sequences per file)
                                seq_count_r1 = 0
                                reader = Mycelia.open_fastx(result_no_dedup.out_forward)
                                for record in reader
                                    seq_count_r1 += 1
                                end
                                Test.@test seq_count_r1 == 5  ## All sequences should be retained
                            finally
                                ## Cleanup reader
                                reader !== nothing && close(reader)
                            end
                        end

                        ## Test with deduplication enabled - should remove duplicates
                        Test.@testset "dedup enabled - removes duplicates" begin
                            ## Use different output paths to avoid file caching
                            dedup_out_forward = joinpath(dup_test_dir, "test_dedup_R1.fastp.1.fq.gz")
                            dedup_out_reverse = joinpath(dup_test_dir, "test_dedup_R2.fastp.2.fq.gz")
                            dedup_json = joinpath(dup_test_dir, "test_dedup_R.fastp_report.json")
                            dedup_html = joinpath(dup_test_dir, "test_dedup_R.fastp_report.html")

                            result_with_dedup = nothing
                            reader = nothing
                            try
                                result_with_dedup = Mycelia.qc_filter_short_reads_fastp(
                                    forward_reads = dup_fastq_R1,
                                    reverse_reads = dup_fastq_R2,
                                    out_forward = dedup_out_forward,
                                    out_reverse = dedup_out_reverse,
                                    json = dedup_json,
                                    html = dedup_html,
                                    enable_dedup = true
                                )

                                ## Count sequences in output (should have fewer sequences after dedup)
                                seq_count_r1 = 0
                                reader = Mycelia.open_fastx(result_with_dedup.out_forward)
                                for record in reader
                                    seq_count_r1 += 1
                                end

                                ## Should have removed duplicates - expecting 3 unique sequences:
                                ## seq1 (first occurrence), seq3_unique, seq4_another_unique
                                Test.@test seq_count_r1 < 5  ## Should have fewer than original
                                Test.@test seq_count_r1 >= 3  ## Should retain at least the 3 unique sequences

                                ## Verify the JSON report contains deduplication stats
                                Test.@test isfile(result_with_dedup.json)
                                json_content = JSON.parse(read(result_with_dedup.json, String))
                                Test.@test haskey(json_content, "duplication")  ## fastp should report duplication stats
                            finally
                                ## Cleanup reader
                                reader !== nothing && close(reader)
                            end
                        end
                    finally
                        ## Cleanup test directory and all its contents
                        isdir(dup_test_dir) &&
                            rm(dup_test_dir, recursive = true, force = true)
                    end
                end
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
                result = Mycelia.qc_filter_long_reads_fastplong(
                    in_fastq = fastq_single
                )
                Test.@test isfile(result.out_fastq)
                Test.@test isfile(result.html_report)
                Test.@test isfile(result.json_report)
                Test.@test endswith(result.out_fastq, ".fastplong.fq.gz")
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
                results_df = CSV.read(results_file, DataFrames.DataFrame, delim = '\t')
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
            run(pipeline(`gzip -c $fasta_file`, stdout = gz_fasta))
            run(pipeline(`gzip -c $fastq_file`, stdout = gz_fastq))

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
                ("GCF_000005825.2_ASM582v2_genomic.fna", "GCF_000005825.2",
                    "RefSeq assembly (avoid meaningless prefix)"),
                ("sample_data_v2.fastq", "sample_data_v2",
                    "Descriptive name with underscores"),
                ("long-file-name.fq", "long-file-name", "Hyphenated name"),
                ("my_genome_assembly_v1_final.fasta", "my_genome",
                    "Long name with meaningful truncation"),
                ("E_coli_K12.fasta", "E_coli_K12", "Species name with underscore")
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
                    Test.@test startswith(table.fastx_identifier[1], expected_id * "_")
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

                table = Mycelia.fastx2normalized_table(filepath; human_readable_id = "custom_id")
                Test.@test table.human_readable_id[1] == "custom_id"
                Test.@test startswith(table.fastx_identifier[1], "custom_id_")
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
                table = Mycelia.fastx2normalized_table(filepath; force_truncate = true)
                Test.@test length(table.human_readable_id[1]) <= 16

                ## Should error without force_truncate (test this by catching the error)
                Test.@test_throws ErrorException Mycelia.fastx2normalized_table(
                    filepath; force_truncate = false)
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
                    human_readable_id = "this_id_is_way_too_long_and_exceeds_sixteen_chars"
                )

                ## Should work with force_truncate
                table = Mycelia.fastx2normalized_table(
                    filepath;
                    human_readable_id = "this_id_is_way_too_long",
                    force_truncate = true
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

                table = Mycelia.fastx2normalized_table(filepath; human_readable_id = "test_genome")

                ## Check that all identifiers have proper hierarchy
                Test.@test all(table.human_readable_id .== "test_genome")
                Test.@test all(length.(table.fastx_hash) .== 16)  ## Base58 encoded, 16 chars
                Test.@test all(length.(table.sequence_hash) .== 16)  ## Base58 encoded, 16 chars

                ## fastx_identifier = human_readable_id + "_" + fastx_hash
                for row in DataFrames.eachrow(table)
                    Test.@test row.fastx_identifier == "test_genome_" * row.fastx_hash
                    Test.@test row.sequence_identifier ==
                               row.fastx_identifier * "_" * row.sequence_hash
                    Test.@test length(row.sequence_identifier) <= 50  ## NCBI limit
                end

                ## All sequences should have the same fastx_identifier but different sequence_identifier
                Test.@test length(unique(table.fastx_identifier)) == 1
                Test.@test length(unique(table.sequence_identifier)) == 2
            end
        end
    end

    Test.@testset "Serialization - Arrow and JLD2" begin
        mktempdir() do dir
            ## Create a test FASTA file to get a normalized table
            fasta_path = joinpath(dir, "test.fasta")
            open(fasta_path, "w") do io
                println(io, ">seq1")
                println(io, "ACGTACGTACGTACGT")
                println(io, ">seq2")
                println(io, "GGGGCCCCAAAATTTT")
            end

            ## Get normalized table
            normalized_table = Mycelia.fastx2normalized_table(fasta_path; human_readable_id = "test_genome")

            Test.@testset "Arrow serialization" begin
                ## Test write_arrow and read_arrow round trip
                arrow_file = joinpath(dir, "test.arrow")

                ## Write the table
                written_file = Mycelia.write_arrow(normalized_table; filename = arrow_file)
                Test.@test written_file == arrow_file
                Test.@test isfile(arrow_file)

                ## Read back the table
                read_table = Mycelia.read_arrow(arrow_file)

                ## Check that sanitized version matches (Arrow converts mixed types to strings)
                sanitized_table = Mycelia.sanitize_for_arrow(normalized_table)
                Test.@test isequal(read_table, sanitized_table)

                ## Verify key properties are preserved
                Test.@test DataFrames.nrow(read_table) == DataFrames.nrow(normalized_table)
                Test.@test names(read_table) == names(normalized_table)
                Test.@test read_table.human_readable_id ==
                           normalized_table.human_readable_id
                Test.@test read_table.record_sequence == normalized_table.record_sequence

                ## Test with auto-sanitization (default behavior)
                arrow_file2 = joinpath(dir, "test2.arrow")
                written_file2 = Mycelia.write_arrow(normalized_table; filename = arrow_file2, sanitize = true)
                read_table2 = Mycelia.read_arrow(arrow_file2)
                Test.@test isequal(read_table2, sanitized_table)
            end

            Test.@testset "JLD2 serialization" begin
                ## Test JLD2_write_table and JLD2_read_table
                jld2_file = joinpath(dir, "test.jld2")

                ## Write the table using JLD2
                Mycelia.JLD2_write_table(df = normalized_table, filename = jld2_file)
                Test.@test isfile(jld2_file)

                ## Read back the table
                read_table = Mycelia.JLD2_read_table(jld2_file)

                ## JLD2 should preserve exact data types and values
                Test.@test isequal(read_table, normalized_table)
                Test.@test DataFrames.nrow(read_table) == DataFrames.nrow(normalized_table)
                Test.@test names(read_table) == names(normalized_table)

                ## Test save_df_jld2 and load_df_jld2 functions
                jld2_file2 = joinpath(dir, "test2")  ## Extension will be added automatically by save
                Mycelia.save_df_jld2(df = normalized_table, filename = jld2_file2)
                Test.@test isfile(jld2_file2 * ".jld2")

                ## Need to include .jld2 extension when loading since load doesn't auto-add it
                read_table2 = Mycelia.load_df_jld2(jld2_file2 * ".jld2")
                Test.@test isequal(read_table2, normalized_table)

                ## Test with custom key
                jld2_file3 = joinpath(dir, "test3.jld2")
                Mycelia.save_df_jld2(df = normalized_table, filename = jld2_file3, key = "normalized_genome")
                read_table3 = Mycelia.load_df_jld2(jld2_file3, key = "normalized_genome")
                Test.@test isequal(read_table3, normalized_table)
            end

            Test.@testset "Sanitization behavior" begin
                ## Test sanitize_for_arrow function behavior
                sanitized = Mycelia.sanitize_for_arrow(normalized_table)

                ## Should have same structure but potentially different types
                Test.@test DataFrames.nrow(sanitized) == DataFrames.nrow(normalized_table)
                Test.@test names(sanitized) == names(normalized_table)

                ## String columns should be preserved as strings
                for col in names(sanitized)
                    if eltype(normalized_table[!, col]) <: AbstractString
                        Test.@test eltype(sanitized[!, col]) <: AbstractString
                    end
                end

                ## Check that missing values are handled correctly
                Test.@test sum(ismissing.(sanitized.record_quality)) ==
                           sum(ismissing.(normalized_table.record_quality))
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
                "Q30(%)", "AvgQual", "GC(%)", "sum_n", "N90"
            ]

            stats_fasta = Mycelia.fastx_stats(fasta_path)
            Test.@test names(stats_fasta) == expected_cols
            Test.@test stats_fasta.num_seqs[1] == 2
            Test.@test stats_fasta.sum_len[1] == 8
            Test.@test isapprox(stats_fasta[1, "GC(%)"], 75.0; atol = 0.01)

            stats_fastq = Mycelia.fastx_stats(fastq_path)
            Test.@test stats_fastq.num_seqs[1] == 2
            Test.@test isapprox(stats_fastq.AvgQual[1], 40.0; atol = 0.01)

            ## Updated column names for fastx2normalized_table
            norm_cols = [
                "fastx_path", "human_readable_id", "fastx_hash", "sequence_hash",
                "fastx_identifier", "sequence_identifier", "record_identifier",
                "record_description", "record_length", "record_alphabet", "record_type",
                "mean_record_quality", "median_record_quality", "record_quality",
                "record_sequence"
            ]

            nfasta = Mycelia.fastx2normalized_table(fasta_path)
            nfastq = Mycelia.fastx2normalized_table(fastq_path)
            Test.@test names(nfasta) == norm_cols
            Test.@test names(nfastq) == norm_cols
            Test.@test DataFrames.nrow(nfasta) == 2
            gc_val = sum(count(c -> c in ['G', 'C'], row.record_sequence)
            for row in DataFrames.eachrow(nfasta)) / sum(nfasta.record_length) * 100
            Test.@test isapprox(gc_val, 75.0; atol = 0.01)
            Test.@test isapprox(Statistics.mean(skipmissing(nfastq.mean_record_quality)), 40.0; atol = 0.01)
        end
    end

    Test.@testset "Test Suite Cleanup Verification" begin
        ## Verify no test files are left behind in the project root
        project_root = dirname(dirname(dirname(@__FILE__)))  # Go up from test/2_preprocessing_qc to project root
        workspace_files = readdir(project_root)

        ## Check for common test file patterns that shouldn't be in workspace root
        test_patterns = [
            r"\.fasta$", r"\.fastq$", r"\.fna$", r"\.fq$", r"\.arrow$", r"\.jld2$"]
        test_files = []

        for file in workspace_files
            for pattern in test_patterns
                if occursin(pattern, file)
                    push!(test_files, file)
                end
            end
        end

        if !isempty(test_files)
            @warn "Found lingering test files in workspace root: $(test_files)"
        end
        Test.@test isempty(test_files)

        ## Check that fastx-stats-test directory was cleaned up
        if isdir("/workspaces/Mycelia/fastx-stats-test")
            @warn "fastx-stats-test directory was not cleaned up"
        end
        Test.@test !isdir("/workspaces/Mycelia/fastx-stats-test")

        println("✅ Test suite cleanup verification passed - no lingering files found")
    end

    Test.@testset "biological sequence hashing functions" begin
        Test.@testset "create_sequence_hash unified interface" begin
            test_seq = "ATCGATCGATCG"

            # Test default behavior (Blake3 + Base58)
            default_hash = Mycelia.create_sequence_hash(test_seq)
            Test.@test isa(default_hash, String)
            Test.@test length(default_hash) == 64  # Default Blake3 length for tree-of-life scale

            # Test different hash algorithms
            blake3_hash = Mycelia.create_sequence_hash(test_seq, hash_function = :blake3, encoded_length = 32)
            sha256_hash = Mycelia.create_sequence_hash(test_seq, hash_function = :sha256,
                encoded_length = 32, allow_truncation = true)
            md5_hash = Mycelia.create_sequence_hash(test_seq, hash_function = :md5,
                encoded_length = 20, allow_truncation = true)  # MD5 is 16 bytes, produces ~22 Base58 chars

            Test.@test length(blake3_hash) == 32
            Test.@test length(sha256_hash) == 32
            Test.@test length(md5_hash) == 20

            # All should be different
            Test.@test length(Set([blake3_hash, sha256_hash, md5_hash])) == 3

            # Test different encodings
            hex_hash = Mycelia.create_sequence_hash(test_seq, encoding = :hex, encoded_length = 32)
            base58_hash = Mycelia.create_sequence_hash(test_seq, encoding = :base58, encoded_length = 32)
            base64_hash = Mycelia.create_sequence_hash(test_seq, encoding = :base64, encoded_length = 32)

            Test.@test length(hex_hash) == 32
            Test.@test length(base58_hash) == 32
            Test.@test length(base64_hash) == 32

            # All encodings should be different
            Test.@test length(Set([hex_hash, base58_hash, base64_hash])) == 3
        end

        Test.@testset "create_base58_hash backwards compatibility" begin
            test_seq = "ATCGATCGATCG"

            # Test backwards compatibility with existing function
            old_style_hash = Mycelia.create_base58_hash(test_seq, encoded_length = 32)
            new_style_hash = Mycelia.create_sequence_hash(
                test_seq, hash_function = :blake3, encoding = :base58, encoded_length = 32)

            Test.@test length(old_style_hash) == 32
            Test.@test length(new_style_hash) == 32

            # Should produce the same result (both use Blake3 + Base58)
            Test.@test old_style_hash == new_style_hash

            # Test case normalization in backwards compatibility
            upper_hash = Mycelia.create_base58_hash("ATCG", encoded_length = 16)
            lower_hash = Mycelia.create_base58_hash("atcg", encoded_length = 16)
            Test.@test upper_hash == lower_hash  # Case normalization should be enabled by default
        end

        Test.@testset "generate_joint_sequence_hash with full parameter control" begin
            sequences = ["ATCG", "GCTA", "TTAA", "CCGG"]

            # Test order independence (core feature)
            joint_hash1 = Mycelia.generate_joint_sequence_hash(sequences, encoded_length = 32)
            joint_hash2 = Mycelia.generate_joint_sequence_hash(reverse(sequences), encoded_length = 32)
            Test.@test joint_hash1 == joint_hash2

            # Test with different hash functions
            blake3_joint = Mycelia.generate_joint_sequence_hash(
                sequences, hash_function = :blake3, encoded_length = 32)
            sha256_joint = Mycelia.generate_joint_sequence_hash(
                sequences, hash_function = :sha256,
                encoded_length = 32, allow_truncation = true)

            Test.@test length(blake3_joint) == 32
            Test.@test length(sha256_joint) == 32
            Test.@test blake3_joint != sha256_joint  # Different algorithms should produce different hashes

            # Test with different encodings
            base58_joint = Mycelia.generate_joint_sequence_hash(
                sequences, encoding = :base58, encoded_length = 32)
            hex_joint = Mycelia.generate_joint_sequence_hash(
                sequences, encoding = :hex, encoded_length = 32, allow_truncation = true)

            Test.@test length(base58_joint) == 32
            Test.@test length(hex_joint) == 32
            Test.@test base58_joint != hex_joint  # Different encodings should produce different representations

            # Test case normalization in joint hashing
            mixed_case_seqs = ["atcg", "GCTA", "TtAa"]
            upper_case_seqs = ["ATCG", "GCTA", "TTAA"]

            mixed_joint = Mycelia.generate_joint_sequence_hash(mixed_case_seqs, encoded_length = 24)
            upper_joint = Mycelia.generate_joint_sequence_hash(upper_case_seqs, encoded_length = 24)
            Test.@test mixed_joint == upper_joint  # Case normalization should work

            # Test with case normalization disabled
            mixed_joint_no_norm = Mycelia.generate_joint_sequence_hash(
                mixed_case_seqs, encoded_length = 24, normalize_case = false)
            upper_joint_no_norm = Mycelia.generate_joint_sequence_hash(
                upper_case_seqs, encoded_length = 24, normalize_case = false)
            Test.@test mixed_joint_no_norm != upper_joint_no_norm  # Should be different without case normalization
        end

        Test.@testset "sequence hash integration with fastx2normalized_table" begin
            mktempdir() do dir
                # Create test FASTA file
                fasta_path = joinpath(dir, "hash_test.fasta")
                open(fasta_path, "w") do io
                    println(io, ">seq1")
                    println(io, "ATCGATCGATCG")
                    println(io, ">seq2")
                    println(io, "GGGGCCCCAAAA")
                end

                # Test that fastx2normalized_table produces consistent hashes
                table1 = Mycelia.fastx2normalized_table(fasta_path; human_readable_id = "test_sample1")
                table2 = Mycelia.fastx2normalized_table(fasta_path; human_readable_id = "test_sample2")

                # Sequence hashes should be identical (same sequences)
                Test.@test table1.sequence_hash == table2.sequence_hash

                # fastx hashes should be identical (same sequence content enables deduplication)
                Test.@test table1.fastx_hash == table2.fastx_hash

                # Verify hash lengths are correct (16 characters for preprocessing functions)
                Test.@test all(length.(table1.sequence_hash) .== 16)
                Test.@test all(length.(table1.fastx_hash) .== 16)

                # Test that identical sequences produce identical hashes
                Test.@test length(unique(table1.sequence_hash)) == 2  # Two unique sequences

                # Test case normalization integration
                # Create mixed case version
                fasta_mixed_path = joinpath(dir, "hash_test_mixed.fasta")
                open(fasta_mixed_path, "w") do io
                    println(io, ">seq1")
                    println(io, "atcgatcgatcg")  # lowercase
                    println(io, ">seq2")
                    println(io, "GGGGCCCCAAAA")  # uppercase
                end

                table_mixed = Mycelia.fastx2normalized_table(
                    fasta_mixed_path; human_readable_id = "test_sample1")

                # Sequence hashes should be identical due to case normalization
                Test.@test table1.sequence_hash == table_mixed.sequence_hash
            end
        end

        Test.@testset "hash function edge cases and validation" begin
            # Test empty sequence handling
            Test.@test_throws ErrorException Mycelia.create_sequence_hash("")  # Should handle empty sequences gracefully

            # Test very short sequences
            short_hash = Mycelia.create_sequence_hash("A", encoded_length = 16)
            Test.@test length(short_hash) == 16

            # Test sequences with non-standard characters (should still hash)
            nonstandard_hash = Mycelia.create_sequence_hash("ATCGN", encoded_length = 16)
            Test.@test length(nonstandard_hash) == 16

            # Test that different sequences produce different hashes
            hash_a = Mycelia.create_sequence_hash("AAAA", encoded_length = 16)
            hash_t = Mycelia.create_sequence_hash("TTTT", encoded_length = 16)
            hash_g = Mycelia.create_sequence_hash("GGGG", encoded_length = 16)
            hash_c = Mycelia.create_sequence_hash("CCCC", encoded_length = 16)

            Test.@test length(Set([hash_a, hash_t, hash_g, hash_c])) == 4  # All unique

            # Test joint hashing with edge cases
            Test.@test_throws ErrorException Mycelia.generate_joint_sequence_hash(String[])  # Empty vector

            single_seq_joint = Mycelia.generate_joint_sequence_hash(["ATCG"], encoded_length = 16)
            Test.@test length(single_seq_joint) == 16

            # Test that joint hashing with identical sequences is deterministic
            identical_seqs = ["ATCG", "ATCG", "ATCG"]
            joint1 = Mycelia.generate_joint_sequence_hash(identical_seqs, encoded_length = 16)
            joint2 = Mycelia.generate_joint_sequence_hash(identical_seqs, encoded_length = 16)
            Test.@test joint1 == joint2
        end

        Test.@testset "tree-of-life scale hash defaults" begin
            test_seq = "ATCGATCGATCG"

            # Test that default lengths are suitable for tree-of-life applications
            default_sequence_hash = Mycelia.create_sequence_hash(test_seq)
            Test.@test length(default_sequence_hash) == 64  # 64 characters provides ~380 bits entropy

            default_base58_hash = Mycelia.create_base58_hash(test_seq)
            Test.@test length(default_base58_hash) == 64  # Should match sequence hash default

            # Test that joint sequence hashing also defaults appropriately
            sequences = ["ATCG", "GCTA"]
            default_joint_hash = Mycelia.generate_joint_sequence_hash(sequences)
            Test.@test length(default_joint_hash) == 64  # Should use same default

            # Verify collision resistance by testing many sequences
            many_sequences = ["A" * string(i) * "T" * string(i) for i in 1:100]
            many_hashes = [Mycelia.create_sequence_hash(seq, encoded_length = 32)
                           for seq in many_sequences]
            Test.@test length(Set(many_hashes)) == 100  # All should be unique
        end
    end
end
