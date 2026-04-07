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

const RUN_ALL = lowercase(get(ENV, "MYCELIA_RUN_ALL", "false")) == "true"
const RUN_EXTERNAL = RUN_ALL || lowercase(get(ENV, "MYCELIA_RUN_EXTERNAL", "false")) == "true"

function missing_checkm_envs()::Vector{String}
    required_envs = ["checkv", "checkm-genome", "checkm2"]
    return [
        env_name for env_name in required_envs if !Mycelia._check_conda_env_exists(env_name)
    ]
end

Test.@testset "CheckM Tools" begin
    Test.@testset "Setup Functions" begin
        if !RUN_EXTERNAL
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true to run CheckM setup tests"
        else
            missing_envs = missing_checkm_envs()
            if !isempty(missing_envs)
                Test.@test_skip "Requires pre-provisioned conda envs: $(join(missing_envs, ", "))"
            else
                # Test that setup functions return database paths
                Test.@test isa(Mycelia.setup_checkv(), String)
                Test.@test isa(Mycelia.setup_checkm(), String)
                Test.@test isa(Mycelia.setup_checkm2(), String)
            end
        end
    end

    Test.@testset "FASTA File Detection" begin
        # Test FASTA file detection with existing constants
        Test.@test occursin(Mycelia.FASTA_REGEX, "test.fasta")
        Test.@test occursin(Mycelia.FASTA_REGEX, "test.fa")
        Test.@test !occursin(Mycelia.FASTA_REGEX, "test.txt")
    end
end
