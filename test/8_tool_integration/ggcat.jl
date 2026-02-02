# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/ggcat.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/ggcat.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia

Test.@testset "GGCAT Tool Integration" begin
    Test.@testset "Function Availability" begin
        Test.@test isdefined(Mycelia, :install_ggcat)
        Test.@test isdefined(Mycelia, :ggcat_build)
        Test.@test isdefined(Mycelia, :ggcat_query)
        Test.@test Mycelia.GGCAT_ENV_NAME isa String
    end

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        Test.@testset "GGCAT Integration (Requires Installation)" begin
            if !isfile(Mycelia.CONDA_RUNNER)
                Test.@test_skip "Conda runner not available; skipping GGCAT integration tests"
            else
                workdir = mktempdir()
                try
                    input_fasta = joinpath(workdir, "input.fasta")
                    open(input_fasta, "w") do io
                        println(io, ">seq1")
                        println(io, "ATGCATGCATGCATGCATGC")
                        println(io, ">seq2")
                        println(io, "GCATGCATGCATGCATGCAT")
                    end

                    Mycelia.install_ggcat()
                    Test.@test Mycelia.check_bioconda_env_is_installed(Mycelia.GGCAT_ENV_NAME)

                    k_len = 5
                    graph_output = Mycelia.ggcat_build(
                        input_fasta,
                        joinpath(workdir, "graph.fasta.lz4"),
                        k_len;
                        threads = 1,
                        min_multiplicity = 1
                    )
                    Test.@test isfile(graph_output)
                    Test.@test filesize(graph_output) > 0

                    query_fasta = joinpath(workdir, "query.fasta")
                    open(query_fasta, "w") do io
                        println(io, ">q1")
                        println(io, "ATGCATGC")
                    end

                    query_output = Mycelia.ggcat_query(
                        graph_output,
                        query_fasta,
                        joinpath(workdir, "query_results.txt"),
                        k_len;
                        threads = 1
                    )
                    Test.@test isfile(query_output)
                    Test.@test filesize(query_output) > 0
                finally
                    rm(workdir; recursive = true, force = true)
                end
            end
        end
    else
        Test.@testset "GGCAT Integration (Skipped)" begin
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true to run GGCAT integration tests"
        end
    end
end
