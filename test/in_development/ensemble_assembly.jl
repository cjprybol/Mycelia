# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/in_development/ensemble_assembly.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/in_development/ensemble_assembly.jl", "test/in_development", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Probabilistic ensemble assembly tests

import Pkg
if isinteractive()
    Pkg.activate("..")
end
import Test
import Mycelia
Test.@testset "probabilistic ensemble assembly" begin
    Test.@testset "Illumina" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "Ultima" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "Nanopore" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "PacBio" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "multi-entity, even coverage" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "multi-entity, log-distributed coverage" begin
        Test.@test 1 + 1 == 2
    end
    Test.@testset "multi-platform" begin
        Test.@test 1 + 1 == 2
    end
end
