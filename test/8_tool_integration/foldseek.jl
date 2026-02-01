# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_ALL=true MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/foldseek.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/foldseek.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import Printf

Test.@testset "Foldseek Tool Integration" begin
    Test.@testset "Function Availability" begin
        Test.@test isdefined(Mycelia, :install_foldseek)
        Test.@test isdefined(Mycelia, :foldseek_easy_search)
        Test.@test isdefined(Mycelia, :foldseek_easy_cluster)
        Test.@test isdefined(Mycelia, :foldseek_createdb)
        Test.@test isdefined(Mycelia, :foldseek_databases)
    end

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        Test.@testset "Foldseek Integration (Requires Installation)" begin
            if !isfile(Mycelia.CONDA_RUNNER)
                Test.@test_skip "Conda runner not available; skipping Foldseek integration tests"
            else
                workdir = mktempdir()
                try
                    query_dir = joinpath(workdir, "query")
                    target_dir = joinpath(workdir, "target")
                    mkpath(query_dir)
                    mkpath(target_dir)

                    lines = ["HEADER    TEST STRUCTURE"]
                    atom_id = 1
                    for residue_id in 1:20
                        x = 10.0 + residue_id
                        y = 10.0
                        z = 10.0
                        push!(lines,
                            Printf.@sprintf("ATOM  %5d  N   GLY A%4d    %8.3f%8.3f%8.3f  1.00 20.00           N",
                                atom_id, residue_id, x, y, z,))
                        atom_id += 1
                        push!(lines,
                            Printf.@sprintf("ATOM  %5d  CA  GLY A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C",
                                atom_id, residue_id, x + 0.6, y + 0.2, z + 0.1,))
                        atom_id += 1
                        push!(lines,
                            Printf.@sprintf("ATOM  %5d  C   GLY A%4d    %8.3f%8.3f%8.3f  1.00 20.00           C",
                                atom_id, residue_id, x + 1.2, y + 0.4, z + 0.2,))
                        atom_id += 1
                        push!(lines,
                            Printf.@sprintf("ATOM  %5d  O   GLY A%4d    %8.3f%8.3f%8.3f  1.00 20.00           O",
                                atom_id, residue_id, x + 1.8, y + 0.6, z + 0.3,))
                        atom_id += 1
                    end
                    push!(lines, "TER")
                    push!(lines, "END")
                    pdb_content = join(lines, "\n")

                    query_pdb = joinpath(query_dir, "query.pdb")
                    target_pdb = joinpath(target_dir, "target.pdb")
                    write(query_pdb, pdb_content)
                    write(target_pdb, pdb_content)

                    Mycelia.install_foldseek()
                    Test.@test Mycelia.check_bioconda_env_is_installed("foldseek")

                    db_path = joinpath(workdir, "foldseek_db")
                    Mycelia.foldseek_createdb(query_dir, db_path; threads = 1)
                    db_entries = readdir(dirname(db_path))
                    Test.@test any(startswith(entry, basename(db_path))
                    for entry in db_entries)

                    search_tmp = joinpath(workdir, "tmp_search")
                    mkpath(search_tmp)
                    search_output = joinpath(workdir, "foldseek_search.tsv")
                    Mycelia.foldseek_easy_search(
                        query_dir,
                        target_dir,
                        search_output;
                        tmp_dir = search_tmp,
                        threads = 1
                    )
                    Test.@test isfile(search_output)

                    cluster_tmp = joinpath(workdir, "tmp_cluster")
                    mkpath(cluster_tmp)
                    cluster_prefix = joinpath(workdir, "foldseek_cluster")
                    cluster_output = Mycelia.foldseek_easy_cluster(
                        query_dir,
                        cluster_prefix;
                        tmp_dir = cluster_tmp,
                        threads = 1
                    )
                    Test.@test isfile(cluster_output)
                finally
                    rm(workdir; recursive = true, force = true)
                end
            end
        end
    else
        Test.@testset "Foldseek Integration (Skipped)" begin
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true to run Foldseek integration tests"
        end
    end
end
