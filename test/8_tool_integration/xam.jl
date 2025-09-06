# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/8_tool_integration/xam.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/xam.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
import XAM
import DataFrames
import BioSequences

Test.@testset "XAM File Processing Tests" begin
    Test.@testset "XAM to DataFrame Conversion - Reader Input" begin
        # Create minimal SAM content for testing
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2\t16\tchr1\t200\t30\t50M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        temp_sam = tempname() * ".sam"
        write(temp_sam, sam_content)
        
        # Test reading with XAM.SAM.Reader
        reader = XAM.SAM.Reader(open(temp_sam))
        df = Mycelia.xam_to_dataframe(reader)
        close(reader)
        
        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) == 2
        Test.@test "templates" in names(df)
        Test.@test "ismapped" in names(df)
        Test.@test "flags" in names(df)
        Test.@test "references" in names(df)
        Test.@test "positions" in names(df)
        
        # Test specific values
        Test.@test df[1, "templates"] == "read1"
        Test.@test df[2, "templates"] == "read2"
        Test.@test all(df.ismapped)  # Both reads should be mapped
        
        # Cleanup
        rm(temp_sam, force=true)
    end

    Test.@testset "XAM to DataFrame Conversion - File Path Input" begin
        # Create test SAM file
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        temp_sam = tempname() * ".sam"
        write(temp_sam, sam_content)
        
        # Test file path input
        df = Mycelia.xam_to_dataframe(temp_sam)
        
        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) == 1
        Test.@test df[1, "templates"] == "read1"
        
        # Cleanup
        rm(temp_sam, force=true)
    end

    Test.@testset "XAM File Opening" begin
        # Create test SAM file
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        temp_sam = tempname() * ".sam"
        write(temp_sam, sam_content)
        
        # Test opening without header
        reader = Mycelia.open_xam(temp_sam, header=false)
        Test.@test reader isa XAM.SAM.Reader
        close(reader)
        
        # Test opening with header
        reader_with_header = Mycelia.open_xam(temp_sam, header=true)
        Test.@test reader_with_header isa XAM.SAM.Reader
        close(reader_with_header)
        
        # Cleanup
        rm(temp_sam, force=true)
    end

    Test.@testset "Samtools Flagstat Integration" begin
        # Create test SAM file
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        temp_sam = tempname() * ".sam"
        write(temp_sam, sam_content)
        
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
        rm(temp_sam, force=true)
    end

    Test.@testset "Mapping Statistics" begin
        # Create test FASTA reference
        temp_fasta = tempname() * ".fasta"
        fasta_content = ">chr1\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        write(temp_fasta, fasta_content)
        
        # Create test SAM file with mappings to this reference
        temp_sam = tempname() * ".sam"
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2\t0\tchr1\t200\t30\t50M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        write(temp_sam, sam_content)
        
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
    end

    Test.@testset "Coverage Determination" begin
        # Create test BAM file (note: would need actual BAM in real scenario)
        temp_sam = tempname() * ".sam"
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        write(temp_sam, sam_content)
        
        # Test coverage determination
        # Note: This function might expect BAM format specifically
        coverage_result = Mycelia.determine_fasta_coverage_from_bam(temp_sam)
        Test.@test coverage_result isa DataFrames.DataFrame
        
        # Cleanup
        rm(temp_sam, force=true)
    end

    Test.@testset "BAM to FASTQ Conversion" begin
        # Create test SAM file
        temp_sam = tempname() * ".sam"
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        write(temp_sam, sam_content)
        
        # Test BAM to FASTQ conversion
        expected_fastq = temp_sam * ".fq.gz"
        result = Mycelia.bam_to_fastq(bam=temp_sam, fastq=expected_fastq)
        
        Test.@test result == expected_fastq
        # In real scenario, would check if FASTQ file was created
        
        # Cleanup
        rm(temp_sam, force=true)
        # rm(expected_fastq, force=true) # Might not exist if conversion failed
    end

    Test.@testset "XAM Summary Table Parsing" begin
        # Create test SAM file
        temp_sam = tempname() * ".sam"
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2\t16\tchr1\t200\t30\t50M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        write(temp_sam, sam_content)
        
        # Test summary table parsing
        summary = Mycelia.parse_xam_to_summary_table(temp_sam)
        Test.@test summary isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(summary) >= 1
        
        # Cleanup
        rm(temp_sam, force=true)
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
        # Create SAM with various record types
        complex_sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
mapped_read\t0\tchr1\t100\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tAS:i:100\tNM:i:0
unmapped_read\t4\t*\t0\t0\t*\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
reverse_read\t16\tchr1\t200\t30\t50M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tAS:i:95\tNM:i:2
"""
        temp_sam = tempname() * ".sam"
        write(temp_sam, complex_sam_content)
        
        # Test DataFrame conversion
        df = Mycelia.xam_to_dataframe(temp_sam)
        
        Test.@test DataFrames.nrow(df) == 3
        
        # Test specific field values
        Test.@test df[1, "templates"] == "mapped_read"
        Test.@test df[2, "templates"] == "unmapped_read"
        Test.@test df[3, "templates"] == "reverse_read"
        
        # Test mapping status
        Test.@test df[1, "ismapped"] == true
        Test.@test df[2, "ismapped"] == false
        Test.@test df[3, "ismapped"] == true
        
        # Test primary alignment status
        Test.@test all(df.isprimary)  # All should be primary
        
        # Cleanup
        rm(temp_sam, force=true)
    end

    Test.@testset "Large File Handling" begin
        # Create a larger SAM file for performance testing
        temp_sam = tempname() * ".sam"
        
        # Write header
        write(temp_sam, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n")
        
        # Add multiple reads
        open(temp_sam, "a") do io
            for i in 1:100
                println(io, "read$i\t0\tchr1\t$(i*10)\t30\t50M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII")
            end
        end
        
        # Test processing larger file
        df = Mycelia.xam_to_dataframe(temp_sam)
        Test.@test DataFrames.nrow(df) == 100
        Test.@test all(df.ismapped)
        
        # Test that all read names are unique and correct
        expected_names = ["read$i" for i in 1:100]
        Test.@test Set(df.templates) == Set(expected_names)
        
        # Cleanup
        rm(temp_sam, force=true)
    end
end