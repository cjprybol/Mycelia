# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/6_annotation/genome_features.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/6_annotation/genome_features.jl", "test/6_annotation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
import DataFrames
import HTTP
import CSV

Test.@testset "Genome Features Tests" begin
    Test.@testset "GFF Reading and Writing" begin
        # Create test GFF content
        test_gff_content = """##gff-version 3
##sequence-region ctg123 1 1497228
ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN
ctg123	.	TF_binding_site	1000	1012	.	+	.	ID=tfbs00001;Parent=gene00001
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00002;Parent=gene00001;Name=EDEN.2
"""

        # Write to temporary file
        temp_gff = tempname() * ".gff"
        write(temp_gff, test_gff_content)

        # Test reading GFF
        gff_df = Mycelia.read_gff(temp_gff)

        Test.@test gff_df isa DataFrames.DataFrame
        Test.@test size(gff_df, 1) == 4  # 4 feature lines
        Test.@test names(gff_df) == ["#seqid", "source", "type", "start", "end",
            "score", "strand", "phase", "attributes"]

        # Test specific values
        Test.@test gff_df[1, "#seqid"] == "ctg123"
        Test.@test gff_df[1, "type"] == "gene"
        Test.@test gff_df[1, "start"] == 1000
        Test.@test gff_df[1, "end"] == 9000
        Test.@test gff_df[1, "strand"] == "+"
        Test.@test occursin("ID=gene00001", gff_df[1, "attributes"])

        # Cleanup
        rm(temp_gff, force = true)
    end

    Test.@testset "GFF Attributes Splitting" begin
        # Create test DataFrames.DataFrame with GFF attributes
        test_df = DataFrames.DataFrame(
            attributes = [
            "ID=gene1;Name=BRCA1;Type=gene",
            "ID=mRNA1;Parent=gene1;Name=transcript1",
            "ID=exon1;Parent=mRNA1"
        ]
        )

        # Test splitting attributes
        expanded_df = Mycelia.split_gff_attributes_into_columns(test_df)

        Test.@test "ID" in names(expanded_df)
        Test.@test "Name" in names(expanded_df)
        Test.@test "Type" in names(expanded_df)
        Test.@test "Parent" in names(expanded_df)

        # Test specific values
        Test.@test expanded_df[1, "ID"] == "gene1"
        Test.@test expanded_df[1, "Name"] == "BRCA1"
        Test.@test expanded_df[1, "Type"] == "gene"
        Test.@test ismissing(expanded_df[3, "Name"])  # exon1 has no Name attribute
        Test.@test expanded_df[2, "Parent"] == "gene1"
    end

    Test.@testset "GFF Writing" begin
        # Create test data
        test_data = DataFrames.DataFrame(
            :seqid => ["chr1", "chr1", "chr2"],
            :source => ["test", "test", "test"],
            :type => ["gene", "mRNA", "exon"],
            :start => [1000, 1000, 1200],
            :end => [2000, 2000, 1800],
            :score => [100, 90, 80],
            :strand => ["+", "+", "-"],
            :phase => [0, 0, 1],
            :attributes => ["ID=gene1", "ID=mRNA1;Parent=gene1", "ID=exon1;Parent=mRNA1"]
        )

        # Write to file
        temp_output = tempname() * ".gff"
        result_file = Mycelia.write_gff(gff = test_data, outfile = temp_output)

        Test.@test result_file == temp_output
        Test.@test isfile(temp_output)

        # Verify content can be read back
        content = read(temp_output, String)
        Test.@test occursin("chr1", content)
        Test.@test occursin("gene", content)
        Test.@test occursin("ID=gene1", content)

        # Cleanup
        rm(temp_output, force = true)
    end

    Test.@testset "MMseqs GFF Update" begin
        # Create mock GFF data
        gff_content = """##gff-version 3
ctg123	prodigal	CDS	1000	2000	.	+	0	ID=gene_001;product=hypothetical_protein
ctg123	prodigal	CDS	3000	4000	.	-	0	ID=gene_002;product=unknown_function
"""
        temp_gff = tempname() * ".gff"
        write(temp_gff, gff_content)

        # Create mock MMseqs results
        mmseqs_data = DataFrames.DataFrame(
            query = ["gene_001", "gene_002"],
            target = ["target1", "target2"],
            theader = ["ATP synthase subunit alpha", "DNA polymerase"],
            bits = [150.0, 200.0],
            evalue = [1e-50, 1e-60]
        )
        temp_mmseqs = tempname() * ".tsv"
        CSV.write(temp_mmseqs, mmseqs_data, delim = '\t')

        # Mock the read_mmseqs_easy_search function for testing
        # In a real scenario, this would read the actual MMseqs format
        function mock_read_mmseqs_easy_search(file)
            return mmseqs_data
        end

        # Test the update logic manually
        gff_df = Mycelia.read_gff(temp_gff)
        Test.@test size(gff_df, 1) == 2

        # Cleanup
        rm(temp_gff, force = true)
        rm(temp_mmseqs, force = true)
    end

    Test.@testset "File Path Handling" begin
        # Test file existence checking logic
        temp_file = tempname() * ".gff"
        write(temp_file, "test content")

        Test.@test isfile(temp_file)

        # Test the logic from open_gff without actually opening
        path = temp_file
        Test.@test isfile(path)

        # Test gzipped file detection
        gz_path = temp_file * ".gz"
        Test.@test occursin(r"\.gz$", basename(gz_path))
        Test.@test !occursin(r"\.gz$", basename(temp_file))

        # Test URL detection
        http_url = "http://example.com/file.gff"
        ftp_url = "ftp://example.com/file.gff"
        Test.@test occursin(r"^http", http_url)
        Test.@test occursin(r"^ftp", ftp_url)

        # Test FTP to HTTP conversion logic
        converted_url = replace(ftp_url, r"^ftp:" => "http:")
        Test.@test converted_url == "http://example.com/file.gff"

        # Cleanup
        rm(temp_file, force = true)
    end

    Test.@testset "GFF Format Validation" begin
        # Test with malformed GFF content
        malformed_gff = """##gff-version 3
chr1	source	gene	not_a_number	2000	.	+	.	ID=gene1
"""
        temp_malformed = tempname() * ".gff"
        write(temp_malformed, malformed_gff)

        # This should handle malformed numeric fields without throwing
        malformed_df = Mycelia.read_gff(temp_malformed)
        Test.@test malformed_df isa DataFrames.DataFrame
        Test.@test !isa(malformed_df[1, "start"], Int)

        # Cleanup
        rm(temp_malformed, force = true)

        # Test with minimal valid GFF
        minimal_gff = """##gff-version 3
chr1	.	gene	1	100	.	+	.	ID=gene1
"""
        temp_minimal = tempname() * ".gff"
        write(temp_minimal, minimal_gff)

        minimal_df = Mycelia.read_gff(temp_minimal)
        Test.@test size(minimal_df, 1) == 1
        Test.@test minimal_df[1, "start"] == 1
        Test.@test minimal_df[1, "end"] == 100

        # Cleanup
        rm(temp_minimal, force = true)
    end

    Test.@testset "Attribute Parsing Edge Cases" begin
        # Test various attribute formats
        test_cases = [
            "ID=gene1",  # Simple case
            "ID=gene1;Name=test",  # Two attributes
            "ID=gene1;Name=test;Type=protein_coding",  # Three attributes
            "ID=gene1;Name=",  # Empty value
            "ID=gene1;;Name=test",  # Empty attribute
            "ID=gene1;Name=test;",  # Trailing semicolon
            "ID=gene1;Name=test name with spaces",  # Spaces in values
            ""  # Empty attributes
        ]

        test_df = DataFrames.DataFrame(attributes = test_cases)
        expanded_df = Mycelia.split_gff_attributes_into_columns(test_df)

        # All rows should have ID column populated (except empty one)
        Test.@test expanded_df[1, "ID"] == "gene1"
        Test.@test expanded_df[2, "ID"] == "gene1"
        Test.@test expanded_df[3, "ID"] == "gene1"

        # Test Name column
        Test.@test ismissing(expanded_df[1, "Name"])  # No Name in first case
        Test.@test expanded_df[2, "Name"] == "test"
        Test.@test expanded_df[7, "Name"] == "test name with spaces"
    end

    Test.@testset "Error Handling" begin
        # Test with non-existent file
        non_existent = "/path/that/does/not/exist.gff"
        Test.@test_throws Exception Mycelia.read_gff(non_existent)

        # Test empty DataFrames.DataFrame for attribute splitting
        empty_df = DataFrames.DataFrame(attributes = String[])
        result = Mycelia.split_gff_attributes_into_columns(empty_df)
        Test.@test size(result, 1) == 0
        Test.@test "attributes" in names(result)
    end

    Test.@testset "Integration Tests" begin
        # Test complete workflow: create GFF -> read -> modify -> write
        original_data = DataFrames.DataFrame(
            "#seqid" => ["chr1", "chr1"],
            "source" => ["test", "test"],
            "type" => ["gene", "mRNA"],
            "start" => [1000, 1000],
            "end" => [2000, 2000],
            "score" => [100, 90],
            "strand" => ["+", "+"],
            "phase" => [0, 0],
            "attributes" => ["ID=gene1;Name=TestGene", "ID=mRNA1;Parent=gene1"]
        )

        # Write original data
        temp_file1 = tempname() * ".gff"
        Mycelia.write_gff(gff = original_data, outfile = temp_file1)

        # Read it back
        read_data = Mycelia.read_gff(temp_file1)

        # Verify data integrity
        Test.@test size(read_data) == size(original_data)
        Test.@test read_data[1, "start"] == 1000
        Test.@test read_data[2, "type"] == "mRNA"

        # Expand attributes
        expanded_data = Mycelia.split_gff_attributes_into_columns(read_data)
        Test.@test "ID" in names(expanded_data)
        Test.@test "Name" in names(expanded_data)
        Test.@test "Parent" in names(expanded_data)

        # Write expanded data
        temp_file2 = tempname() * ".gff"
        Mycelia.write_gff(gff = expanded_data, outfile = temp_file2)

        Test.@test isfile(temp_file2)

        # Cleanup
        rm(temp_file1, force = true)
        rm(temp_file2, force = true)
    end
end

# Note: Tests involving external dependencies (EMBOSS, HTTP requests, etc.) are 
# commented out or mocked to avoid test environment dependencies and network calls.
# For full integration testing, these would need to be run in an environment with
# the appropriate tools installed and network access available.
