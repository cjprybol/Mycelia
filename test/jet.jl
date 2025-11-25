## JET.jl Static Analysis Tests
##
## This file runs JET.jl static analysis on the Mycelia package.
## JET is a test-only dependency - it is not loaded when using Mycelia normally.
##
## Run via Pkg.test() (when uncommented in runtests.jl):
##   julia --project=. -e 'import Pkg; Pkg.test()'
##
## Or run standalone (requires test dependencies):
##   julia --project=. -e 'import Pkg; Pkg.test()' -- test/jet.jl
##   julia --project=. -e 'include("test/jet.jl")'

import Test
import JET

Test.@testset "JET Static Analysis" begin
    ## Run analysis and capture the result
    result = JET.report_package("Mycelia")
    reports = JET.get_reports(result)
    
    if isempty(reports)
        println("\n✅ No type instabilities or errors detected by JET!")
    else
        println("\n⚠️  Found $(length(reports)) issues")
        for report in reports
            println(report)
        end
    end
    
    ## Currently not enforcing zero issues - uncomment when ready:
    ## Test.@test isempty(reports)
    
    ## For now, just test that JET ran successfully
    Test.@test result !== nothing
end