#!/usr/bin/env julia
# Integration Test Runner for Mycelia
# This runs the full test suite including bioconda-dependent tests
# Suitable for HPC environments with proper conda setup

println("""
================================================================================
Mycelia Integration Test Suite
================================================================================
This will run ALL tests including those requiring external tools via bioconda.
Ensure you have:
- Sufficient memory (>16GB recommended)
- Conda/mamba available
- Internet access for package downloads
- Time (this may take 30+ minutes)
================================================================================
""")

# Set environment to run integration tests
ENV["MYCELIA_TEST_TIER"] = "integration"

# Activate the project
import Pkg
Pkg.activate(dirname(@__DIR__))

# Run the full test suite
println("\nStarting integration tests at $(Dates.now())")
start_time = time()

try
    Pkg.test()
    elapsed = round(time() - start_time, digits=1)
    println("\n✅ All integration tests passed in $elapsed seconds!")
catch e
    elapsed = round(time() - start_time, digits=1)
    println("\n❌ Integration tests failed after $elapsed seconds")
    rethrow(e)
end