# From the Mycelia base directory, run the tests with:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/pantools.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/8_tool_integration/pantools.jl", "test/8_tool_integration", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# PanTools wrapper tests
# To run optional CLI smoke tests (requires conda + network):

import Test
import Mycelia

Test.@testset "PanTools Wrapper Tests" begin
    Test.@testset "Types" begin
        db = Mycelia.PantoolsDB("some_db")
        Test.@test db isa Mycelia.PantoolsDB
        Test.@test db.path == "some_db"
        Test.@test Mycelia.pantools_db_path(db) == "some_db"
        Test.@test Mycelia.pantools_db_path("other_db") == "other_db"
    end

    Test.@testset "Java opts" begin
        Test.@test Mycelia.pantools_java_opts() == String[]
        Test.@test Mycelia.pantools_java_opts(xmx="5g") == ["-Xmx5g"]
        Test.@test Mycelia.pantools_java_opts(xms="2g", xmx="5g") == ["-Xms2g", "-Xmx5g"]
    end

    Test.@testset "Command construction" begin
        cmd_ls = Mycelia.pantools_cmd(["--help"]; env="pantools", live_stream=true)
        Test.@test cmd_ls isa Cmd
        Test.@test string(cmd_ls.exec[1]) == Mycelia.CONDA_RUNNER
        Test.@test "run" in cmd_ls.exec
        Test.@test "--live-stream" in cmd_ls.exec
        Test.@test "-n" in cmd_ls.exec
        Test.@test "pantools" in cmd_ls.exec

        cmd_no_ls = Mycelia.pantools_cmd(["--help"]; env="pantools", live_stream=false)
        Test.@test cmd_no_ls isa Cmd
        Test.@test !("--live-stream" in cmd_no_ls.exec)
    end

    Test.@testset "Input file writers" begin
        workdir = mktempdir()
        created = String[]
        try
            fasta1 = joinpath(workdir, "a.fasta")
            fasta2 = joinpath(workdir, "b.fasta")
            gff1 = joinpath(workdir, "a.gff3")
            gff2 = joinpath(workdir, "b.gff3")

            open(fasta1, "w") do io
                write(io, ">a\nACGT\n")
            end
            open(fasta2, "w") do io
                write(io, ">b\nACGT\n")
            end
            open(gff1, "w") do io
                write(io, "##gff-version 3\n")
            end
            open(gff2, "w") do io
                write(io, "##gff-version 3\n")
            end

            genomes_file = Mycelia.write_pantools_genome_locations_file([fasta1, fasta2])
            push!(created, genomes_file)
            Test.@test isfile(genomes_file)
            lines = readlines(genomes_file)
            Test.@test length(lines) == 2
            Test.@test all(startswith.(lines, workdir))

            ann_file = Mycelia.write_pantools_annotation_locations_file(Dict(2 => gff2, 1 => gff1))
            push!(created, ann_file)
            Test.@test isfile(ann_file)
            ann_lines = readlines(ann_file)
            Test.@test ann_lines[1] == "1 $(abspath(gff1))"
            Test.@test ann_lines[2] == "2 $(abspath(gff2))"

            genome_numbers_file = Mycelia.write_pantools_genome_numbers_file([1, 3, 2])
            push!(created, genome_numbers_file)
            Test.@test isfile(genome_numbers_file)
            Test.@test readlines(genome_numbers_file) == ["1", "3", "2"]

            regions_file = Mycelia.write_pantools_regions_file([(1, "contig1", 1, 10), (2, "contig2", 5, 20, '+')])
            push!(created, regions_file)
            Test.@test isfile(regions_file)
            Test.@test readlines(regions_file) == ["1 contig1 1 10", "2 contig2 5 20 +"]

            Test.@test_throws AssertionError Mycelia.write_pantools_regions_file([(1, "c", 1)])
        finally
            for file in created
                isfile(file) && rm(file)
            end
            rm(workdir; recursive=true, force=true)
        end
    end

    # Optional smoke test: install env + check CLI responds
    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        Test.@testset "PanTools CLI Smoke (Requires Installation)" begin
            if isfile(Mycelia.CONDA_RUNNER) || haskey(ENV, "CONDA_PREFIX")
                Mycelia.add_bioconda_env("pantools")
                Test.@test Mycelia.check_bioconda_env_is_installed("pantools")
                Test.@test success(`$(Mycelia.CONDA_RUNNER) run -n pantools pantools --help`)
            else
                Test.@test_skip "Conda not available, skipping PanTools CLI smoke test"
            end
        end
    else
        @info "Skipping PanTools CLI smoke test; set MYCELIA_RUN_EXTERNAL=true to enable"
    end
end
