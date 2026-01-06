import Test
import Mycelia

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

                    pdb_content = join([
                        "HEADER    TEST STRUCTURE",
                        "ATOM      1  N   MET A   1      11.104  13.207   8.276  1.00 20.00           N",
                        "ATOM      2  CA  MET A   1      12.560  13.105   8.556  1.00 20.00           C",
                        "ATOM      3  C   MET A   1      13.047  11.651   8.835  1.00 20.00           C",
                        "ATOM      4  O   MET A   1      12.396  10.679   8.505  1.00 20.00           O",
                        "TER",
                        "END",
                    ], "\n")

                    query_pdb = joinpath(query_dir, "query.pdb")
                    target_pdb = joinpath(target_dir, "target.pdb")
                    write(query_pdb, pdb_content)
                    write(target_pdb, pdb_content)

                    Mycelia.install_foldseek()
                    Test.@test Mycelia.check_bioconda_env_is_installed("foldseek")

                    db_path = joinpath(workdir, "foldseek_db")
                    Mycelia.foldseek_createdb(query_dir, db_path; threads=1)
                    db_entries = readdir(dirname(db_path))
                    Test.@test any(startswith(entry, basename(db_path)) for entry in db_entries)

                    search_tmp = joinpath(workdir, "tmp_search")
                    mkpath(search_tmp)
                    search_output = joinpath(workdir, "foldseek_search.tsv")
                    Mycelia.foldseek_easy_search(
                        query_dir,
                        target_dir,
                        search_output;
                        tmp_dir=search_tmp,
                        threads=1,
                    )
                    Test.@test isfile(search_output)

                    cluster_tmp = joinpath(workdir, "tmp_cluster")
                    mkpath(cluster_tmp)
                    cluster_prefix = joinpath(workdir, "foldseek_cluster")
                    cluster_output = Mycelia.foldseek_easy_cluster(
                        query_dir,
                        cluster_prefix;
                        tmp_dir=cluster_tmp,
                        threads=1,
                    )
                    Test.@test isfile(cluster_output)
                finally
                    rm(workdir; recursive=true, force=true)
                end
            end
        end
    else
        Test.@testset "Foldseek Integration (Skipped)" begin
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true to run Foldseek integration tests"
        end
    end
end
