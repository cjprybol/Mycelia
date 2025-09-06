# Mycelia Test Suite - Tiered Testing System
#
# Usage:
#   Quick core tests (Julia-only):     julia --project=. -e "import Pkg; Pkg.test()"
#   Integration tests (with bioconda): julia --project=. test/run_integration_tests.jl
#   Benchmarks:                        julia --project=. benchmarking/run_all_benchmarks.jl
#   Tutorials:                         julia --project=. tutorials/run_all_tutorials.jl

# Determine which test tier to run based on environment variable
const TEST_TIER = get(ENV, "MYCELIA_TEST_TIER", "core")

println("Running Mycelia tests - Tier: $TEST_TIER")

# Helper function to include all test files in a directory
function include_all_tests(dir; exclude_patterns=String[])
    test_count = 0
    for (root, dirs, files) in walkdir(dir)
        for file in sort(files)
            if endswith(file, ".jl")
                # Check if file matches any exclude pattern
                should_exclude = any(pattern -> occursin(pattern, file), exclude_patterns)
                if !should_exclude
                    include(joinpath(root, file))
                    test_count += 1
                end
            end
        end
    end
    return test_count
end

# Define which tests require external dependencies
const BIOCONDA_DEPENDENT_TESTS = [
    "bioconda.jl",
    "tool_integration.jl",
    # Add more bioconda-dependent test files here
]

const LONG_RUNNING_TESTS = [
    "end_to_end_assembly_tests.jl",
    "ensemble_assembly.jl",
    # Add more long-running tests here
]

# Run tests based on tier
if TEST_TIER == "core"
    println("\n=== Running Core Tests (Julia-only) ===")
    println("Excluding bioconda-dependent and long-running tests for speed\n")
    
    # Exclude patterns for core testing
    exclude_patterns = vcat(BIOCONDA_DEPENDENT_TESTS, LONG_RUNNING_TESTS)
    
    # Run core tests in logical order
    println("1. Data Acquisition Tests (core)")
    include_all_tests(joinpath(@__DIR__, "1_data_acquisition"), exclude_patterns=exclude_patterns)
    
    println("\n2. Preprocessing & QC Tests (core)")
    include_all_tests(joinpath(@__DIR__, "2_preprocessing_qc"), exclude_patterns=exclude_patterns)
    
    println("\n3. K-mer Analysis Tests (core)")
    include_all_tests(joinpath(@__DIR__, "3_feature_extraction_kmer"), exclude_patterns=exclude_patterns)
    
    # println("\n4. Assembly Tests (core)")
    # include_all_tests(joinpath(@__DIR__, "4_assembly"), exclude_patterns=exclude_patterns)
    
    # println("\n5. Validation Tests (core)")
    # include_all_tests(joinpath(@__DIR__, "5_validation"), exclude_patterns=exclude_patterns)
    
    # Skip these directories entirely for core tests as they're mostly tool integrations
    # include_all_tests(joinpath(@__DIR__, "6_annotation"), exclude_patterns=exclude_patterns)
    # include_all_tests(joinpath(@__DIR__, "7_comparative_pangenomics"), exclude_patterns=exclude_patterns)
    # include_all_tests(joinpath(@__DIR__, "8_tool_integration"), exclude_patterns=exclude_patterns)
    
elseif TEST_TIER == "integration"
    println("\n=== Running Integration Tests (includes external tools) ===")
    println("This may take longer and requires bioconda environment\n")
    
    # Run all tests including bioconda-dependent ones
    include_all_tests(joinpath(@__DIR__, "1_data_acquisition"))
    include_all_tests(joinpath(@__DIR__, "2_preprocessing_qc"))
    include_all_tests(joinpath(@__DIR__, "3_feature_extraction_kmer"))
    include_all_tests(joinpath(@__DIR__, "4_assembly"))
    include_all_tests(joinpath(@__DIR__, "5_validation"))
    include_all_tests(joinpath(@__DIR__, "6_annotation"))
    include_all_tests(joinpath(@__DIR__, "7_comparative_pangenomics"))
    include_all_tests(joinpath(@__DIR__, "8_tool_integration"))
    
elseif TEST_TIER == "quick"
    println("\n=== Running Quick Smoke Tests ===")
    println("Minimal tests to verify basic functionality\n")
    
    # Just run a few critical tests
    include(joinpath(@__DIR__, "1_data_acquisition", "simulation_fasta.jl"))
    include(joinpath(@__DIR__, "3_feature_extraction_kmer", "kmer_analysis.jl"))
    include(joinpath(@__DIR__, "4_assembly", "intelligent_assembly_tests.jl"))
    
else
    error("Unknown test tier: $TEST_TIER. Valid options: core, integration, quick")
end

println("\n✅ Test tier '$TEST_TIER' completed!")