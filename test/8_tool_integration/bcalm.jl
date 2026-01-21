# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/bcalm.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/bcalm.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Bcalm wrapper tests
# To run the integration test (requires conda + network):

import Test
import Mycelia
import FASTX
import MetaGraphsNext

Test.@testset "Bcalm Tool Integration" begin
    Test.@testset "Function Availability" begin
        Test.@test isdefined(Mycelia, :install_bcalm)
        Test.@test isdefined(Mycelia, :run_bcalm)
        Test.@test isdefined(Mycelia, :_bcalm_paths)
        Test.@test Mycelia.BCALM_ENV_NAME isa String
        Test.@test Mycelia.BCALM_GFA_SCRIPT_URL isa String
    end

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        Test.@testset "Bcalm Integration (Requires Installation)" begin
            if !isfile(Mycelia.CONDA_RUNNER)
                Test.@test_skip "Conda runner not available; skipping bcalm integration tests"
            else
                workdir = mktempdir()
                try
                    input1 = joinpath(workdir, "bcalm_reads_1.fasta")
                    input2 = joinpath(workdir, "bcalm_reads_2.fasta")

                    record1 = Mycelia.random_fasta_record(moltype=:DNA, seed=11, L=400)
                    seq1 = FASTX.sequence(record1)
                    records1 = [
                        FASTX.FASTA.Record("read1", seq1),
                        FASTX.FASTA.Record("read2", seq1),
                    ]
                    Mycelia.write_fasta(outfile=input1, records=records1)

                    record2 = Mycelia.random_fasta_record(moltype=:DNA, seed=22, L=400)
                    seq2 = FASTX.sequence(record2)
                    records2 = [
                        FASTX.FASTA.Record("read3", seq2),
                        FASTX.FASTA.Record("read4", seq2),
                    ]
                    Mycelia.write_fasta(outfile=input2, records=records2)

                    Mycelia.install_bcalm()
                    Test.@test Mycelia.check_bioconda_env_is_installed(Mycelia.BCALM_ENV_NAME)
                    _, script_path = Mycelia._bcalm_paths()
                    Test.@test isfile(script_path)

                    out_single = joinpath(workdir, "bcalm_single")
                    result_single = Mycelia.run_bcalm(
                        input1,
                        out_single;
                        kmer_size=21,
                        abundance_min=2,
                        threads=1,
                    )
                    Test.@test isfile(result_single.unitigs)
                    Test.@test filesize(result_single.unitigs) > 0
                    Test.@test isfile(result_single.gfa)

                    gfa_content = read(result_single.gfa, String)
                    gfa_lines = filter(!isempty, split(strip(gfa_content), '\n'))
                    Test.@test !isempty(gfa_lines)
                    Test.@test startswith(gfa_lines[1], "H\t") || startswith(gfa_lines[1], "S\t")

                    graph = Mycelia.Rhizomorph.read_gfa_next(result_single.gfa)
                    Test.@test graph isa MetaGraphsNext.MetaGraph

                    out_multi = joinpath(workdir, "bcalm_multi")
                    result_multi = Mycelia.run_bcalm(
                        [input1, input2],
                        out_multi;
                        kmer_size=21,
                        abundance_min=2,
                        threads=1,
                    )
                    Test.@test isfile(result_multi.unitigs)
                    Test.@test isfile(result_multi.gfa)
                finally
                    rm(workdir; recursive=true, force=true)
                end
            end
        end
    else
        Test.@testset "Bcalm Integration (Skipped)" begin
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true to run bcalm integration tests"
        end
    end
end
