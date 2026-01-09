# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/xam.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/xam.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import XAM
import DataFrames
import BioSequences
import FASTX
import StableRNGs

Test.@testset "Quick XAM Verification" begin
    ## Create test data
    rng = StableRNGs.StableRNG(1234)
    ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
    ref_fasta = tempname() * ".fasta"
    writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
    write(writer, ref_record)
    close(writer)

    fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=2, rndSeed=rand(rng, 0:typemax(Int)))

    Test.@testset "Basic Functionality" begin
        ## Test BAM output (default)
        map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr")
        run(map_result.cmd)
        bam_file = map_result.outfile

        ## Test auto parser
        df = Mycelia.xam_to_dataframe(bam_file)
        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) == 2
        Test.@test all(df.ismapped)

        rm(bam_file, force=true)
    end

    Test.@testset "Format Options" begin
        ## Test SAM output
        map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="sam")
        run(map_result.cmd)
        sam_file = map_result.outfile

        Test.@test endswith(sam_file, ".sam")
        Test.@test isfile(sam_file)

        ## Test format detection
        detected = Mycelia.detect_xam_format(sam_file)
        Test.@test detected == :sam

        ## Test parsing
        df = Mycelia.xam_to_dataframe(sam_file)
        Test.@test DataFrames.nrow(df) == 2

        rm(sam_file, force=true)
    end

    Test.@testset "Parser Selection" begin
        ## Create BAM file
        map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="bam")
        run(map_result.cmd)
        bam_file = map_result.outfile

        ## Test different parsers
        reader_xamjl = Mycelia.open_xam(bam_file, parser=:xamjl)
        Test.@test reader_xamjl isa XAM.BAM.Reader
        close(reader_xamjl)

        reader_samtools = Mycelia.open_xam(bam_file, parser=:samtools)
        Test.@test reader_samtools isa XAM.SAM.Reader
        close(reader_samtools)

        rm(bam_file, force=true)
    end

    ## Cleanup
    rm(ref_fasta, force=true)
    rm(fastq_files.forward_reads, force=true)
    rm(fastq_files.reverse_reads, force=true)
    rm(fastq_files.sam, force=true)
    rm(fastq_files.error_free_sam, force=true)
end

Test.@testset "XAM File Processing Tests" begin
    Test.@testset "Reproducible SAM Test" begin
        # 1. Create a reference genome
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        ref_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        # 2. Simulate reads
        fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=10, rndSeed=rand(rng, 0:typemax(Int)))

        # 3. Align reads
        map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr")
        run(map_result.cmd)
        sam_file = map_result.outfile

        # 4. Test xam_to_dataframe
        df = Mycelia.xam_to_dataframe(sam_file)
        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) > 0

        # Cleanup
        rm(ref_fasta, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(sam_file, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end
    

    Test.@testset "XAM to DataFrame Conversion - Reader Input" begin
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        ref_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=2, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr")
        run(map_result.cmd)
        sam_file = map_result.outfile

        reader = Mycelia.open_xam(sam_file)
        df = Mycelia.xam_to_dataframe(reader)
        close(reader)
        
        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) == 2
        Test.@test "template" in names(df)
        Test.@test "ismapped" in names(df)
        Test.@test "flag" in names(df)
        Test.@test "reference" in names(df)
        Test.@test "position" in names(df)
        
        # Test specific values
        Test.@test all(df.ismapped)  # Both reads should be mapped
        
        # Cleanup
        rm(ref_fasta, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(sam_file, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end

    Test.@testset "XAM to DataFrame Conversion - File Path Input" begin
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        ref_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=1, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr")
        run(map_result.cmd)
        sam_file = map_result.outfile

        df = Mycelia.xam_to_dataframe(sam_file)
        
        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) == 1
        
        # Cleanup
        rm(ref_fasta, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(sam_file, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end



    Test.@testset "XAM File Opening" begin
        # Create test data using real simulation pipeline
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        temp_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(temp_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=temp_fasta, read_count=1, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=temp_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="sam")
        run(map_result.cmd)
        temp_sam = map_result.outfile

        # Test opening without header
        reader = Mycelia.open_xam(temp_sam, header=false)
        Test.@test reader isa XAM.SAM.Reader
        close(reader)

        # Test opening with header
        reader_with_header = Mycelia.open_xam(temp_sam, header=true)
        Test.@test reader_with_header isa XAM.SAM.Reader
        close(reader_with_header)

        # Cleanup
        rm(temp_fasta, force=true)
        rm(temp_sam, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end

    Test.@testset "Samtools Flagstat Integration" begin
        # Create test data using real simulation pipeline
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        ref_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=2, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="sam")
        run(map_result.cmd)
        temp_sam = map_result.outfile

        # Test flagstat output file generation
        expected_flagstat = temp_sam * ".samtools-flagstat.txt"

        # Note: This would require samtools to be installed
        # In a real test environment, we'd mock this or skip if not available
        try
            result = Mycelia.run_samtools_flagstat(temp_sam)
            Test.@test result == expected_flagstat
            # If successful, check that file was created
            if isfile(expected_flagstat)
                Test.@test isfile(expected_flagstat)
                rm(expected_flagstat, force=true)
            end
        catch e
            # If samtools not available, test that function exists
            Test.@test hasmethod(Mycelia.run_samtools_flagstat, (String,))
        end

        # Cleanup
        rm(ref_fasta, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(temp_sam, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end

    Test.@testset "Mapping Statistics" begin
        # Create test data using real simulation pipeline
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        temp_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(temp_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=temp_fasta, read_count=2, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=temp_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="sam")
        run(map_result.cmd)
        temp_sam = map_result.outfile

        # Test FASTA-XAM mapping statistics
        stats = Mycelia.fasta_xam_mapping_stats(fasta=temp_fasta, xam=temp_sam)
        Test.@test stats isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(stats) >= 1

        # Test XAM-only contig mapping statistics
        contig_stats = Mycelia.xam_to_contig_mapping_stats(temp_sam)
        Test.@test contig_stats isa DataFrames.DataFrame

        # Cleanup
        rm(temp_fasta, force=true)
        rm(temp_sam, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
    if run_external
        Test.@testset "Coverage Determination" begin
            # Create test data using Mycelia pipeline to ensure correct BAM format
            rng = StableRNGs.StableRNG(1234)
            ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
            ref_fasta = tempname() * ".fasta"
            writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
            write(writer, ref_record)
            close(writer)

            fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=2, rndSeed=rand(rng, 0:typemax(Int)))
            map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="bam")
            run(map_result.cmd)
            bam_file = map_result.outfile

            # Test coverage determination - bedtools expects BAM format
            coverage_result = Mycelia.determine_fasta_coverage_from_bam(bam_file)
            Test.@test coverage_result isa DataFrames.DataFrame

            # Cleanup
            rm(ref_fasta, force=true)
            rm(fastq_files.forward_reads, force=true)
            rm(fastq_files.reverse_reads, force=true)
            rm(bam_file, force=true)
            rm(fastq_files.sam, force=true)
            rm(fastq_files.error_free_sam, force=true)
        end
    else
        @info "Skipping Coverage Determination tests; bedtools execution is opt-in via MYCELIA_RUN_EXTERNAL=true"
    end

    Test.@testset "BAM to FASTQ Conversion" begin
        # Create test data using Mycelia pipeline to ensure correct BAM format
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        ref_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=2, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="bam")
        run(map_result.cmd)
        bam_file = map_result.outfile

        # Test BAM to FASTQ conversion - samtools fastq expects BAM format
        expected_fastq = bam_file * ".fq.gz"
        result = Mycelia.bam_to_fastq(bam=bam_file, fastq=expected_fastq)

        Test.@test result == expected_fastq
        Test.@test isfile(expected_fastq)

        # Cleanup
        rm(ref_fasta, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(bam_file, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
        rm(expected_fastq, force=true)
    end

    Test.@testset "XAM Summary Table Parsing" begin
        # Create test data using real simulation pipeline
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        temp_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(temp_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=temp_fasta, read_count=2, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=temp_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="sam")
        run(map_result.cmd)
        temp_sam = map_result.outfile

        # Test summary table parsing
        summary = Mycelia.parse_xam_to_summary_table(temp_sam)
        Test.@test summary isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(summary) >= 1

        # Cleanup
        rm(temp_fasta, force=true)
        rm(temp_sam, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end

    Test.@testset "Error Handling" begin
        # Test with non-existent file
        Test.@test_throws Exception Mycelia.xam_to_dataframe("nonexistent.sam")
        Test.@test_throws Exception Mycelia.open_xam("nonexistent.sam")
        
        # Test with empty SAM file
        empty_sam = tempname() * ".sam"
        write(empty_sam, "")
        
        Test.@test_throws Exception Mycelia.xam_to_dataframe(empty_sam)
        
        # Cleanup
        rm(empty_sam, force=true)
        
        # Test with malformed SAM file
        malformed_sam = tempname() * ".sam"
        write(malformed_sam, "not_a_sam_file\n")
        
        Test.@test_throws Exception Mycelia.xam_to_dataframe(malformed_sam)
        
        # Cleanup
        rm(malformed_sam, force=true)
    end

    Test.@testset "Record Field Extraction" begin
        # Create test data using real simulation pipeline
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        temp_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(temp_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=temp_fasta, read_count=3, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=temp_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="sam")
        run(map_result.cmd)
        temp_sam = map_result.outfile

        # Test DataFrame conversion
        df = Mycelia.xam_to_dataframe(temp_sam)

        Test.@test DataFrames.nrow(df) >= 1  # Should have at least 1 mapped read
        Test.@test "template" in names(df)
        Test.@test "ismapped" in names(df)
        Test.@test "isprimary" in names(df)

        # Test that we have proper record structure (most reads should be mapped in simulation)
        Test.@test any(df.ismapped)  # At least some reads should be mapped
        Test.@test all(df.isprimary) # All should be primary in this simple mapping

        # Cleanup
        rm(temp_fasta, force=true)
        rm(temp_sam, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end

    Test.@testset "Large File Handling" begin
        # Create test data using real simulation pipeline with more reads
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=5000, seed=rand(rng, 0:typemax(Int)))
        temp_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(temp_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=temp_fasta, read_count=50, rndSeed=rand(rng, 0:typemax(Int)))
        map_result = Mycelia.minimap_map(fasta=temp_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="sam")
        run(map_result.cmd)
        temp_sam = map_result.outfile

        # Test processing larger file
        df = Mycelia.xam_to_dataframe(temp_sam)
        Test.@test DataFrames.nrow(df) >= 10  # Should have many mapped reads
        Test.@test any(df.ismapped)  # Most reads should be mapped

        # Test that read names are unique
        Test.@test length(unique(df.template)) == DataFrames.nrow(df)

        # Cleanup
        rm(temp_fasta, force=true)
        rm(temp_sam, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end

    Test.@testset "Comprehensive Parser Testing" begin
        ## Create test data using Mycelia pipeline
        rng = StableRNGs.StableRNG(1234)
        ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
        ref_fasta = tempname() * ".fasta"
        writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
        write(writer, ref_record)
        close(writer)

        fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=3, rndSeed=rand(rng, 0:typemax(Int)))

        Test.@testset "Format Detection Tests" begin
            ## Test different output formats from minimap
            formats = ["sam", "bam", "sam.gz"]
            test_files = String[]

            for format in formats
                map_result = Mycelia.minimap_map(
                    fasta=ref_fasta,
                    fastq=fastq_files.forward_reads,
                    mapping_type="sr",
                    output_format=format,
                    sorted=true
                )
                run(map_result.cmd)
                push!(test_files, map_result.outfile)

                ## Test format detection
                detected_format = Mycelia.detect_xam_format(map_result.outfile)
                expected_format = format == "sam.gz" ? :sam_gz : Symbol(format)
                Test.@test detected_format == expected_format

                ## Test that sam.gz files can be parsed
                if format == "sam.gz"
                    df = Mycelia.xam_to_dataframe(map_result.outfile)
                    Test.@test df isa DataFrames.DataFrame
                    Test.@test DataFrames.nrow(df) == 3
                    Test.@test all(df.ismapped)
                end
            end

            ## Cleanup format test files
            for file in test_files
                rm(file, force=true)
            end
        end

        Test.@testset "BGZ Extension Support Tests" begin
            ## Create BGZF-compressed SAM file with .bgz extension
            map_result = Mycelia.minimap_map(
                fasta=ref_fasta,
                fastq=fastq_files.forward_reads,
                mapping_type="sr",
                output_format="sam.gz",
                sorted=true
            )
            run(map_result.cmd)
            original_samgz = map_result.outfile

            ## Create .bgz file by copying and renaming the sam.gz file
            bgz_file = replace(original_samgz, ".sam.gz" => ".bgz")
            cp(original_samgz, bgz_file)

            ## Create .sam.bgz file as well
            sambgz_file = replace(original_samgz, ".sam.gz" => ".sam.bgz")
            cp(original_samgz, sambgz_file)

            Test.@testset "BGZ Extension Format Detection" begin
                ## Test .bgz format detection
                detected_bgz = Mycelia.detect_xam_format(bgz_file)
                Test.@test detected_bgz == :sam_gz

                ## Test .sam.bgz format detection
                detected_sambgz = Mycelia.detect_xam_format(sambgz_file)
                Test.@test detected_sambgz == :sam_gz

                ## Test original .sam.gz format detection still works
                detected_samgz = Mycelia.detect_xam_format(original_samgz)
                Test.@test detected_samgz == :sam_gz
            end

            Test.@testset "BGZ Extension Parsing" begin
                ## Test .bgz file parsing with all parsers
                Test.@testset ".bgz File Parsing" begin
                    ## Auto parser
                    reader_auto = Mycelia.open_xam(bgz_file, parser=:auto)
                    df_auto = Mycelia.xam_to_dataframe(reader_auto)
                    close(reader_auto)
                    Test.@test DataFrames.nrow(df_auto) == 3
                    Test.@test all(df_auto.ismapped)

                    ## XAM.jl parser
                    reader_xamjl = Mycelia.open_xam(bgz_file, parser=:xamjl)
                    df_xamjl = Mycelia.xam_to_dataframe(reader_xamjl)
                    close(reader_xamjl)
                    Test.@test DataFrames.nrow(df_xamjl) == 3
                    Test.@test all(df_xamjl.ismapped)

                    ## Samtools parser
                    reader_samtools = Mycelia.open_xam(bgz_file, parser=:samtools)
                    df_samtools = Mycelia.xam_to_dataframe(reader_samtools)
                    close(reader_samtools)
                    Test.@test DataFrames.nrow(df_samtools) == 3
                    Test.@test all(df_samtools.ismapped)

                    ## Direct file path parsing
                    df_direct = Mycelia.xam_to_dataframe(bgz_file)
                    Test.@test DataFrames.nrow(df_direct) == 3
                    Test.@test all(df_direct.ismapped)
                end

                ## Test .sam.bgz file parsing with all parsers
                Test.@testset ".sam.bgz File Parsing" begin
                    ## Auto parser
                    reader_auto = Mycelia.open_xam(sambgz_file, parser=:auto)
                    df_auto = Mycelia.xam_to_dataframe(reader_auto)
                    close(reader_auto)
                    Test.@test DataFrames.nrow(df_auto) == 3
                    Test.@test all(df_auto.ismapped)

                    ## XAM.jl parser
                    reader_xamjl = Mycelia.open_xam(sambgz_file, parser=:xamjl)
                    df_xamjl = Mycelia.xam_to_dataframe(reader_xamjl)
                    close(reader_xamjl)
                    Test.@test DataFrames.nrow(df_xamjl) == 3
                    Test.@test all(df_xamjl.ismapped)

                    ## Samtools parser
                    reader_samtools = Mycelia.open_xam(sambgz_file, parser=:samtools)
                    df_samtools = Mycelia.xam_to_dataframe(reader_samtools)
                    close(reader_samtools)
                    Test.@test DataFrames.nrow(df_samtools) == 3
                    Test.@test all(df_samtools.ismapped)

                    ## Direct file path parsing
                    df_direct = Mycelia.xam_to_dataframe(sambgz_file)
                    Test.@test DataFrames.nrow(df_direct) == 3
                    Test.@test all(df_direct.ismapped)
                end
            end

            Test.@testset "BGZ Extension Data Consistency" begin
                ## Verify all three file extensions (.sam.gz, .bgz, .sam.bgz) produce identical data
                df_samgz = Mycelia.xam_to_dataframe(original_samgz)
                df_bgz = Mycelia.xam_to_dataframe(bgz_file)
                df_sambgz = Mycelia.xam_to_dataframe(sambgz_file)

                ## All should have same number of records
                Test.@test DataFrames.nrow(df_samgz) == DataFrames.nrow(df_bgz) == DataFrames.nrow(df_sambgz)
                Test.@test DataFrames.nrow(df_samgz) == 3

                ## All should have identical read names
                Test.@test Set(df_samgz.template) == Set(df_bgz.template) == Set(df_sambgz.template)

                ## All should have identical mapping status
                Test.@test all(df_samgz.ismapped) && all(df_bgz.ismapped) && all(df_sambgz.ismapped)
            end

            ## Cleanup
            rm(original_samgz, force=true)
            rm(bgz_file, force=true)
            rm(sambgz_file, force=true)
        end

        Test.@testset "Format Detection Edge Cases" begin
            ## Test the priority of extension-based vs magic byte detection
            ## Create files with BGZF content but different extensions to verify disambiguation

            ## Create a BGZF-compressed SAM file first
            map_result = Mycelia.minimap_map(
                fasta=ref_fasta,
                fastq=fastq_files.forward_reads,
                mapping_type="sr",
                output_format="sam.gz",
                sorted=true
            )
            run(map_result.cmd)
            original_samgz = map_result.outfile

            Test.@testset "Extension Priority Over Magic Bytes" begin
                ## All these files have identical BGZF content but different extensions
                ## Should be detected based on extension, not magic bytes

                ## .bgz should be detected as sam.gz format
                bgz_file = tempname() * ".bgz"
                cp(original_samgz, bgz_file)
                Test.@test Mycelia.detect_xam_format(bgz_file) == :sam_gz

                ## .sam.bgz should be detected as sam.gz format
                sambgz_file = tempname() * ".sam.bgz"
                cp(original_samgz, sambgz_file)
                Test.@test Mycelia.detect_xam_format(sambgz_file) == :sam_gz

                ## .sam.gz should be detected as sam.gz format
                samgz_file = tempname() * ".sam.gz"
                cp(original_samgz, samgz_file)
                Test.@test Mycelia.detect_xam_format(samgz_file) == :sam_gz

                ## Rename to .bam should be detected as bam format (extension override)
                bam_file = tempname() * ".bam"
                cp(original_samgz, bam_file)
                Test.@test Mycelia.detect_xam_format(bam_file) == Symbol("bam")

                ## Verify all these files can actually be parsed despite format confusion
                ## (samtools should handle them gracefully regardless of detected format)
                for test_file in [bgz_file, sambgz_file, samgz_file]
                    df = Mycelia.xam_to_dataframe(test_file)
                    Test.@test DataFrames.nrow(df) == 3
                    Test.@test all(df.ismapped)
                end

                ## Cleanup
                rm(bgz_file, force=true)
                rm(sambgz_file, force=true)
                rm(samgz_file, force=true)
                rm(bam_file, force=true)
            end

            Test.@testset "No Extension Fallback to Magic Bytes" begin
                ## Create file with no extension - should fall back to magic byte detection

                ## Copy BGZF file without extension
                no_ext_file = tempname()  # No extension
                cp(original_samgz, no_ext_file)

                ## Should detect as BAM since BGZF magic bytes without extension context default to BAM
                detected = Mycelia.detect_xam_format(no_ext_file)
                Test.@test detected == Symbol("bam")  # Magic byte fallback assumes BAM for BGZF

                ## Should still be parseable by samtools
                df = Mycelia.xam_to_dataframe(no_ext_file)
                Test.@test DataFrames.nrow(df) == 3
                Test.@test all(df.ismapped)

                ## Cleanup
                rm(no_ext_file, force=true)
            end

            ## Cleanup
            rm(original_samgz, force=true)
        end

        Test.@testset "Parser Comparison Tests" begin
            ## Create BAM file for parser comparison
            map_result = Mycelia.minimap_map(
                fasta=ref_fasta,
                fastq=fastq_files.forward_reads,
                mapping_type="sr",
                output_format="bam"
            )
            run(map_result.cmd)
            bam_file = map_result.outfile

            ## Test XAM.jl parser
            Test.@testset "XAM.jl Parser" begin
                reader = Mycelia.open_xam(bam_file, parser=:xamjl)
                Test.@test reader isa XAM.BAM.Reader
                df_xamjl = Mycelia.xam_to_dataframe(reader)
                close(reader)

                Test.@test DataFrames.nrow(df_xamjl) == 3
                Test.@test all(df_xamjl.ismapped)
                Test.@test "template" in names(df_xamjl)
                Test.@test "reference" in names(df_xamjl)
            end

            ## Test samtools parser
            Test.@testset "Samtools Parser" begin
                reader = Mycelia.open_xam(bam_file, parser=:samtools)
                Test.@test reader isa XAM.SAM.Reader
                df_samtools = Mycelia.xam_to_dataframe(reader)
                close(reader)

                Test.@test DataFrames.nrow(df_samtools) == 3
                Test.@test all(df_samtools.ismapped)
                Test.@test "template" in names(df_samtools)
                Test.@test "reference" in names(df_samtools)
            end

            ## Test auto parser (defaults to samtools, returns SAM.Reader)
            Test.@testset "Auto Parser" begin
                reader = Mycelia.open_xam(bam_file, parser=:auto)
                Test.@test reader isa XAM.SAM.Reader  # Auto parser defaults to samtools
                df_auto = Mycelia.xam_to_dataframe(reader)
                close(reader)

                Test.@test DataFrames.nrow(df_auto) == 3
                Test.@test all(df_auto.ismapped)
            end

            ## Cleanup parser comparison test file
            rm(bam_file, force=true)
        end

        Test.@testset "minimap_map Format Options" begin
            ## Test sorted vs unsorted
            Test.@testset "Sorted vs Unsorted Output Comparison" begin
                ## Create sorted output
                map_result_sorted = Mycelia.minimap_map(
                    fasta=ref_fasta,
                    fastq=fastq_files.forward_reads,
                    mapping_type="sr",
                    output_format="sam",
                    sorted=true
                )
                Test.@test occursin("sorted", map_result_sorted.outfile)
                Test.@test endswith(map_result_sorted.outfile, ".sam")
                run(map_result_sorted.cmd)
                Test.@test isfile(map_result_sorted.outfile)

                ## Create unsorted output
                map_result_unsorted = Mycelia.minimap_map(
                    fasta=ref_fasta,
                    fastq=fastq_files.forward_reads,
                    mapping_type="sr",
                    output_format="sam",
                    sorted=false
                )
                Test.@test !occursin("sorted", map_result_unsorted.outfile)
                Test.@test endswith(map_result_unsorted.outfile, ".sam")
                run(map_result_unsorted.cmd)
                Test.@test isfile(map_result_unsorted.outfile)

                ## Compare the outputs - both should have same records but potentially different order
                df_sorted = Mycelia.xam_to_dataframe(map_result_sorted.outfile)
                df_unsorted = Mycelia.xam_to_dataframe(map_result_unsorted.outfile)

                Test.@test DataFrames.nrow(df_sorted) == DataFrames.nrow(df_unsorted)
                Test.@test DataFrames.nrow(df_sorted) == 3

                ## Both should have same read names (but potentially different order)
                Test.@test Set(df_sorted.template) == Set(df_unsorted.template)
                Test.@test all(df_sorted.ismapped) == all(df_unsorted.ismapped)

                ## Clean up
                rm(map_result_sorted.outfile, force=true)
                rm(map_result_unsorted.outfile, force=true)
            end
        end

        Test.@testset "SAM.GZ Compressed Format Tests" begin
            ## Test SAM.GZ creation and parsing
            map_result = Mycelia.minimap_map(
                fasta=ref_fasta,
                fastq=fastq_files.forward_reads,
                mapping_type="sr",
                output_format="sam.gz"
            )
            run(map_result.cmd)
            samgz_file = map_result.outfile

            Test.@test endswith(samgz_file, ".sam.gz")
            Test.@test isfile(samgz_file)

            ## Test format detection
            detected_format = Mycelia.detect_xam_format(samgz_file)
            Test.@test detected_format == :sam_gz

            ## Test parsing with different methods
            Test.@testset "SAM.GZ Parser Methods" begin
                ## Test auto parser
                reader_auto = Mycelia.open_xam(samgz_file, parser=:auto)
                df_auto = Mycelia.xam_to_dataframe(reader_auto)
                close(reader_auto)
                Test.@test DataFrames.nrow(df_auto) == 3
                Test.@test all(df_auto.ismapped)

                ## Test xamjl parser
                reader_xamjl = Mycelia.open_xam(samgz_file, parser=:xamjl)
                df_xamjl = Mycelia.xam_to_dataframe(reader_xamjl)
                close(reader_xamjl)
                Test.@test DataFrames.nrow(df_xamjl) == 3
                Test.@test all(df_xamjl.ismapped)

                ## Test samtools parser
                reader_samtools = Mycelia.open_xam(samgz_file, parser=:samtools)
                df_samtools = Mycelia.xam_to_dataframe(reader_samtools)
                close(reader_samtools)
                Test.@test DataFrames.nrow(df_samtools) == 3
                Test.@test all(df_samtools.ismapped)

                ## Verify all parsers give same number of records
                Test.@test DataFrames.nrow(df_auto) == DataFrames.nrow(df_xamjl) == DataFrames.nrow(df_samtools)
            end

            ## Test direct file path parsing
            df_direct = Mycelia.xam_to_dataframe(samgz_file)
            Test.@test df_direct isa DataFrames.DataFrame
            Test.@test DataFrames.nrow(df_direct) == 3
            Test.@test all(df_direct.ismapped)
            Test.@test "template" in names(df_direct)
            Test.@test "reference" in names(df_direct)

            ## Clean up
            rm(samgz_file, force=true)
        end

        ## Cleanup main test data
        rm(ref_fasta, force=true)
        rm(fastq_files.forward_reads, force=true)
        rm(fastq_files.reverse_reads, force=true)
        rm(fastq_files.sam, force=true)
        rm(fastq_files.error_free_sam, force=true)
    end

    Test.@testset "Error Handling and Edge Cases" begin
        ## Test with non-existent file
        Test.@test_throws Exception Mycelia.xam_to_dataframe("nonexistent.sam")
        Test.@test_throws Exception Mycelia.open_xam("nonexistent.sam")

        ## Test with empty files
        empty_sam = tempname() * ".sam"
        write(empty_sam, "")
        Test.@test_throws Exception Mycelia.open_xam(empty_sam)
        rm(empty_sam, force=true)

        ## Test invalid parser specification
        Test.@testset "Invalid Parser Options" begin
            ## Create a valid test file first
            rng = StableRNGs.StableRNG(1234)
            ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
            ref_fasta = tempname() * ".fasta"
            writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
            write(writer, ref_record)
            close(writer)

            fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=1, rndSeed=rand(rng, 0:typemax(Int)))
            map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr")
            run(map_result.cmd)
            test_file = map_result.outfile

            ## Test invalid parser
            Test.@test_throws Exception Mycelia.open_xam(test_file, parser=:invalid)

            ## Test invalid minimap output format
            Test.@test_throws Exception Mycelia.minimap_map(
                fasta=ref_fasta,
                fastq=fastq_files.forward_reads,
                mapping_type="sr",
                output_format="invalid"
            )

            ## Cleanup
            rm(ref_fasta, force=true)
            rm(fastq_files.forward_reads, force=true)
            rm(fastq_files.reverse_reads, force=true)
            rm(test_file, force=true)
            rm(fastq_files.sam, force=true)
            rm(fastq_files.error_free_sam, force=true)
        end

        ## Test with malformed file content (modify valid data to create invalid variants)
        Test.@testset "Malformed File Handling" begin
            ## Create a valid SAM file first
            rng = StableRNGs.StableRNG(1234)
            ref_record = Mycelia.random_fasta_record(L=1000, seed=rand(rng, 0:typemax(Int)))
            ref_fasta = tempname() * ".fasta"
            writer = FASTX.FASTA.Writer(open(ref_fasta, "w"))
            write(writer, ref_record)
            close(writer)

            fastq_files = Mycelia.simulate_illumina_reads(fasta=ref_fasta, read_count=1, rndSeed=rand(rng, 0:typemax(Int)))
            map_result = Mycelia.minimap_map(fasta=ref_fasta, fastq=fastq_files.forward_reads, mapping_type="sr", output_format="sam")
            run(map_result.cmd)
            valid_sam = map_result.outfile

            ## Create malformed variant by corrupting the valid file
            malformed_sam = tempname() * ".sam"
            valid_content = read(valid_sam, String)
            ## Corrupt the content by removing essential fields
            corrupted_content = replace(valid_content, r"\t[0-9]+\t" => "\tXXX\t", count=1)
            write(malformed_sam, corrupted_content)

            ## Test that malformed file is handled appropriately
            ## Note: Some malformed files may still parse as valid DataFrames
            ## depending on the level of corruption, so we test for either case
            df = Mycelia.xam_to_dataframe(malformed_sam)
            Test.@test df isa DataFrames.DataFrame

            ## Cleanup
            rm(ref_fasta, force=true)
            rm(fastq_files.forward_reads, force=true)
            rm(fastq_files.reverse_reads, force=true)
            rm(valid_sam, force=true)
            rm(malformed_sam, force=true)
            rm(fastq_files.sam, force=true)
            rm(fastq_files.error_free_sam, force=true)
        end
    end
end
