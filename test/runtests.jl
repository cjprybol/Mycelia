# Mycelia Test Suite - Tiered Testing System
#
# Usage:
#   Core tests: julia --project=. -e "import Pkg; Pkg.test()"
#   Benchmarks: julia --project=. benchmarking/run_all_benchmarks.jl
#   Tutorials:  julia --project=. tutorials/run_all_tutorials.jl

println("Running Mycelia tests")

# Helper function to include all test files in a directory
function include_all_tests(dir)
    test_count = 0
    for (root, dirs, files) in walkdir(dir)
        for file in sort(files)
            if endswith(file, ".jl")
                include(joinpath(root, file))
                test_count += 1
            end
        end
    end
    return test_count
end

println("\n=== Running Tests ===")

# Run all tests including bioconda-dependent ones
include_all_tests(joinpath(@__DIR__, "1_data_acquisition"))
include_all_tests(joinpath(@__DIR__, "2_preprocessing_qc"))
include_all_tests(joinpath(@__DIR__, "3_feature_extraction_kmer"))
# include_all_tests(joinpath(@__DIR__, "4_assembly"))
# include_all_tests(joinpath(@__DIR__, "5_validation"))
# include_all_tests(joinpath(@__DIR__, "6_annotation"))
# include_all_tests(joinpath(@__DIR__, "7_comparative_pangenomics"))
# include_all_tests(joinpath(@__DIR__, "8_tool_integration"))

println("\n✅ Tests completed!")