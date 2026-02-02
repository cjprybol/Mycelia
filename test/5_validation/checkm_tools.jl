# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/5_validation/checkm_tools.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/5_validation/checkm_tools.jl", "test/5_validation", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# CheckM, CheckM2, CheckV tool tests

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia

Test.@testset "CheckM Tools" begin
    Test.@testset "Setup Functions" begin
        # Test that setup functions return database paths
        Test.@test isa(Mycelia.setup_checkv(), String)
        Test.@test isa(Mycelia.setup_checkm(), String)
        Test.@test isa(Mycelia.setup_checkm2(), String)
    end

    Test.@testset "FASTA File Detection" begin
        # Test FASTA file detection with existing constants
        Test.@test occursin(Mycelia.FASTA_REGEX, "test.fasta")
        Test.@test occursin(Mycelia.FASTA_REGEX, "test.fa")
        Test.@test !occursin(Mycelia.FASTA_REGEX, "test.txt")
    end
end
