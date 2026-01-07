# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/metagraph.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/metagraph.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# MetaGraph wrapper tests
# To run the integration test (requires conda + network):

import Test
import Mycelia

Test.@testset "MetaGraph Tool Integration" begin
    Test.@testset "Function Availability" begin
        Test.@test isdefined(Mycelia, :install_metagraph)
        Test.@test isdefined(Mycelia, :metagraph_cmd)
        Test.@test isdefined(Mycelia, :run_metagraph)
        Test.@test isdefined(Mycelia, :metagraph_build)
        Test.@test isdefined(Mycelia, :metagraph_annotate)
        Test.@test isdefined(Mycelia, :metagraph_query)
        Test.@test isdefined(Mycelia, :metagraph_align)
        Test.@test isdefined(Mycelia, :metagraph_transform)
        Test.@test isdefined(Mycelia, :metagraph_transform_anno)
        Test.@test isdefined(Mycelia, :metagraph_relax_brwt)
        Test.@test isdefined(Mycelia, :metagraph_server_query)
        Test.@test isdefined(Mycelia, :metagraph_assemble)
        Test.@test isdefined(Mycelia, :metagraph_stats)
        Test.@test isdefined(Mycelia, :metagraph_clean)
        Test.@test Mycelia.METAGRAPH_CONDA_ENV isa String
    end

    Test.@testset "Command Construction" begin
        cmd = Mycelia.metagraph_cmd(["--version"])
        Test.@test cmd isa Cmd
        Test.@test string(cmd.exec[1]) == Mycelia.CONDA_RUNNER
        Test.@test "metagraph" in cmd.exec

        cmd_prot = Mycelia.metagraph_cmd(["build"]; executable="metagraph_Protein")
        Test.@test "metagraph_Protein" in cmd_prot.exec
    end

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        Test.@testset "MetaGraph Integration (Requires Conda)" begin
            if !isfile(Mycelia.CONDA_RUNNER) && !haskey(ENV, "CONDA_PREFIX")
                Test.@test_skip "Conda runner not available; skipping MetaGraph integration tests"
            else
                try
                    Mycelia.install_metagraph()
                    Test.@test Mycelia.check_bioconda_env_is_installed(Mycelia.METAGRAPH_CONDA_ENV)
                    Mycelia.run_metagraph(["--version"]; live_stream=false)
                    Test.@test true
                catch
                    Test.@test false
                end
            end
        end

        Test.@testset "MetaGraph Basic Workflow" begin
            if Mycelia.check_bioconda_env_is_installed(Mycelia.METAGRAPH_CONDA_ENV)
                mktempdir() do tmp_dir
                    fasta_path = joinpath(tmp_dir, "test.fasta")
                    open(fasta_path, "w") do io
                        write(io, ">seq1\nACGTACGTACGTACGT\n")
                    end

                    Mycelia.metagraph_build([fasta_path];
                        k=15,
                        outfile_base="graph",
                        out_dir=tmp_dir,
                        verbose=false
                    )

                    dbg_files = filter(path -> endswith(path, ".dbg"), readdir(tmp_dir; join=true))
                    Test.@test !isempty(dbg_files)

                    Mycelia.metagraph_stats(first(dbg_files))
                end
            else
                Test.@test_skip "MetaGraph environment not available; skipping workflow test"
            end
        end
    else
        Test.@testset "MetaGraph Integration (Skipped)" begin
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true to run MetaGraph integration tests"
        end
    end
end
