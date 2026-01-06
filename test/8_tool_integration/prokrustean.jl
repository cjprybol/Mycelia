# Prokrustean wrapper tests
#
# From the Mycelia base directory, run:
#   julia --project=. -e 'include("test/8_tool_integration/prokrustean.jl")'
#
# To run the integration test (requires build tools + conda + network):
#   MYCELIA_RUN_EXTERNAL=true julia --project=. -e 'include("test/8_tool_integration/prokrustean.jl")'

import Test
import Mycelia

Test.@testset "Prokrustean Tool Integration" begin
    Test.@testset "Function Availability" begin
        Test.@test isdefined(Mycelia, :install_prokrustean)
        Test.@test isdefined(Mycelia, :prokrustean_build_graph)
        Test.@test isdefined(Mycelia, :prokrustean_kmer_count)
        Test.@test isdefined(Mycelia, :prokrustean_unitig_count)
        Test.@test isdefined(Mycelia, :prokrustean_braycurtis)
        Test.@test isdefined(Mycelia, :prokrustean_overlap)
    end

    run_all = get(ENV, "MYCELIA_RUN_ALL", "false") == "true"
    run_external = run_all || get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"

    if run_external
        Test.@testset "Prokrustean Integration (Requires Build Tools)" begin
            missing_tools = String[]
            for tool in ("git", "cmake", "make")
                if Sys.which(tool) === nothing
                    push!(missing_tools, tool)
                end
            end

            if !isempty(missing_tools)
                Test.@test_skip "Missing build tools: $(join(missing_tools, ", "))"
            elseif !isfile(Mycelia.CONDA_RUNNER)
                Test.@test_skip "Conda runner not available; skipping ropebwt2 integration"
            else
                workdir = mktempdir()
                try
                    Mycelia.add_bioconda_env("ropebwt2")
                    if !Mycelia.check_bioconda_env_is_installed("ropebwt2")
                        Test.@test_skip "ropebwt2 environment not available"
                    else
                        fasta_path = joinpath(workdir, "test.fasta")
                        open(fasta_path, "w") do io
                            write(io, ">seq1\nACGTACGTACGTACGT\n")
                            write(io, ">seq2\nTACGTACGTACGTACG\n")
                            write(io, ">seq3\nACGTACGT\n")
                        end

                        bwt_path = joinpath(workdir, "test.bwt")
                        run(pipeline(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ropebwt2 ropebwt2 $fasta_path`, stdout=bwt_path))
                        Test.@test isfile(bwt_path)

                        install_dir = joinpath(workdir, "prokrustean_install")
                        build_ok = true
                        try
                            Mycelia.install_prokrustean(dest_dir=install_dir, force=true)
                            Test.@test isfile(joinpath(install_dir, "build", "prokrustean"))
                        catch e
                            Test.@test_skip "Prokrustean build failed: $(e)"
                            build_ok = false
                        end

                        if build_ok
                            graph_path = joinpath(workdir, "test.prokrustean")
                            Mycelia.prokrustean_build_graph(bwt_path, graph_path; kmin=5, install_dir=install_dir)
                            Test.@test isfile(graph_path)

                            kmer_counts_file = joinpath(workdir, "kmer_counts.txt")
                            Mycelia.prokrustean_kmer_count(graph_path, kmer_counts_file; install_dir=install_dir)
                            Test.@test isfile(kmer_counts_file)
                            Test.@test filesize(kmer_counts_file) > 0

                            unitig_counts_file = joinpath(workdir, "unitig_counts.txt")
                            Mycelia.prokrustean_unitig_count(graph_path, unitig_counts_file; install_dir=install_dir)
                            Test.@test isfile(unitig_counts_file)
                            lines = readlines(unitig_counts_file)
                            Test.@test !isempty(lines)
                            Test.@test all(line -> tryparse(Float64, split(strip(line))[1]) !== nothing, lines)

                            overlap_file = joinpath(workdir, "overlap.txt")
                            Mycelia.prokrustean_overlap(graph_path, overlap_file; install_dir=install_dir)
                            Test.@test isfile(overlap_file)
                        end
                    end
                finally
                    rm(workdir; recursive=true, force=true)
                end
            end
        end
    else
        Test.@testset "Prokrustean Integration (Skipped)" begin
            Test.@test_skip "Set MYCELIA_RUN_EXTERNAL=true to run prokrustean integration tests"
        end
    end
end
