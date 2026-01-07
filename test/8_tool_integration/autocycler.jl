# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/autocycler.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/autocycler.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# Autocycler wrapper tests
# To run the install smoke test (requires conda + network):
# To run the pipeline (requires installed env + data paths):
#   MYCELIA_RUN_EXTERNAL=true \
#   MYCELIA_AUTOCYCLER_LONG_READS=/path/to/reads.fastq \

import Test
import Mycelia

Test.@testset "Autocycler Wrapper Tests" begin
    Test.@testset "Constants and paths" begin
        Test.@test Mycelia.AUTOCYCLER_ENV_NAME isa String
        Test.@test !isempty(Mycelia.AUTOCYCLER_ENV_NAME)
        Test.@test Mycelia.AUTOCYCLER_SCRIPT_URL isa String
        Test.@test Mycelia.AUTOCYCLER_ENV_URL isa String

        install_dir, script_path, env_file_path = Mycelia._autocycler_paths()
        Test.@test endswith(script_path, "autocycler_full.sh")
        Test.@test endswith(env_file_path, "environment.yml")
        Test.@test endswith(install_dir, joinpath("deps", "autocycler"))
    end

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        Test.@testset "Install Smoke Test" begin
            if isfile(Mycelia.CONDA_RUNNER) || haskey(ENV, "CONDA_PREFIX")
                script_path = Mycelia.install_autocycler()
                Test.@test isfile(script_path)
                Test.@test Mycelia.check_bioconda_env_is_installed(Mycelia.AUTOCYCLER_ENV_NAME)
            else
                Test.@test_skip "Conda not available, skipping Autocycler install smoke test"
            end
        end
    else
        @info "Skipping Autocycler install smoke test; set MYCELIA_RUN_EXTERNAL=true to enable"
    end

    if run_external
        Test.@testset "Pipeline Smoke Test" begin
            long_reads = get(ENV, "MYCELIA_AUTOCYCLER_LONG_READS", "")
            short_reads_1 = get(ENV, "MYCELIA_AUTOCYCLER_SHORT_READS_1", "")
            short_reads_2 = get(ENV, "MYCELIA_AUTOCYCLER_SHORT_READS_2", "")

            if isempty(long_reads)
                Test.@test_skip "Set MYCELIA_AUTOCYCLER_LONG_READS to run Autocycler pipeline"
            elseif (isempty(short_reads_1) && !isempty(short_reads_2)) || (!isempty(short_reads_1) && isempty(short_reads_2))
                Test.@test_skip "Set both MYCELIA_AUTOCYCLER_SHORT_READS_1 and MYCELIA_AUTOCYCLER_SHORT_READS_2 for short reads"
            else
                mktempdir() do tmp_dir
                    output_dir = joinpath(tmp_dir, "autocycler_out")
                    if isempty(short_reads_1)
                        Mycelia.run_autocycler(long_reads=long_reads, out_dir=output_dir)
                    else
                        Mycelia.run_autocycler(
                            long_reads=long_reads,
                            out_dir=output_dir,
                            short_reads_1=short_reads_1,
                            short_reads_2=short_reads_2
                        )
                    end
                    Test.@test isdir(output_dir)
                end
            end
        end
    else
        @info "Skipping Autocycler pipeline test; set MYCELIA_RUN_EXTERNAL=true to enable"
    end
end
