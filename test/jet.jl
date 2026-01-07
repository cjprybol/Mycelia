# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/jet.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/jet.jl", "test", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("..")
## end
## using Revise

## JET.jl Static Analysis Tests
## This file runs JET.jl static analysis on the Mycelia package.
## JET is a test-only dependency - it is not loaded when using Mycelia normally.
## Run via Pkg.test() (when uncommented in runtests.jl):
## Or run standalone (requires test dependencies):

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
