# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/8_tool_integration/sourmash_parsers_test.jl")'
# ```

import Test
import Mycelia
import CSV
import DataFrames

Test.@testset "Sourmash Output Parsers" begin
    Test.@testset "parse_sourmash_gather_output" begin
        # Create a temporary test file with sourmash gather CSV format
        test_dir = mktempdir()
        test_file = joinpath(test_dir, "test_gather.csv")

        # Sourmash gather output columns (standard format)
        gather_content = """intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,name,filename,md5,f_match_orig,unique_intersect_bp,gather_result_rank,remaining_bp,query_filename,query_name,query_md5,query_bp
1000000,0.5,0.8,0.45,0.42,2.1,2.0,0.5,genome1,ref.sig,abc123,0.9,950000,0,500000,sample.sig,sample1,def456,2000000
500000,0.25,0.6,0.22,0.20,1.5,1.0,0.3,genome2,ref.sig,ghi789,0.7,480000,1,0,sample.sig,sample1,def456,2000000"""

        open(test_file, "w") do io
            write(io, gather_content)
        end

        # Test parsing
        df = Mycelia.parse_sourmash_gather_output(test_file)

        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) == 2
        Test.@test "intersect_bp" in DataFrames.names(df)
        Test.@test "f_orig_query" in DataFrames.names(df)
        Test.@test "name" in DataFrames.names(df)
        Test.@test df.name[1] == "genome1"
        Test.@test df.intersect_bp[1] == 1000000

        # Test error on missing file
        Test.@test_throws ErrorException Mycelia.parse_sourmash_gather_output("nonexistent.csv")

        rm(test_dir, recursive = true)
    end

    Test.@testset "parse_sourmash_search_output" begin
        # Create a temporary test file with sourmash search CSV format
        test_dir = mktempdir()
        test_file = joinpath(test_dir, "test_search.csv")

        # Sourmash search output columns (standard format)
        search_content = """similarity,name,filename,md5
0.95,genome1,ref1.sig,abc123
0.85,genome2,ref2.sig,def456
0.75,genome3,ref3.sig,ghi789"""

        open(test_file, "w") do io
            write(io, search_content)
        end

        # Test parsing
        df = Mycelia.parse_sourmash_search_output(test_file)

        Test.@test df isa DataFrames.DataFrame
        Test.@test DataFrames.nrow(df) == 3
        Test.@test "similarity" in DataFrames.names(df)
        Test.@test "name" in DataFrames.names(df)
        Test.@test df.similarity[1] == 0.95
        Test.@test df.name[2] == "genome2"

        # Test error on missing file
        Test.@test_throws ErrorException Mycelia.parse_sourmash_search_output("nonexistent.csv")

        rm(test_dir, recursive = true)
    end
end

println("All sourmash parser tests passed!")
