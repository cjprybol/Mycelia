# From the Mycelia base directory, run these tests with:
#
# ```bash
# julia --project=. --compiled-modules=no -e \
#   'include("test/8_tool_integration/autocycler.jl")'
# ```
#
# A real smoke test is opt-in because Autocycler is compute intensive:
#
# ```bash
# MYCELIA_RUN_EXTERNAL=true \
# MYCELIA_AUTOCYCLER_LONG_READS=/path/to/long.fastq \
# MYCELIA_AUTOCYCLER_SHORT_READS_1=/path/to/R1.fastq \
# MYCELIA_AUTOCYCLER_SHORT_READS_2=/path/to/R2.fastq \
# julia --project=. --compiled-modules=no -e \
#   'include("test/8_tool_integration/autocycler.jl")'
# ```

import Mycelia
import Test

function _autocycler_test_runner!(step::NamedTuple)::Nothing
    for output_path in step.expected_outputs
        mkpath(dirname(output_path))
        if endswith(output_path, ".gfa")
            write(
                output_path,
                "H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\tdp:f:20.0\n",
            )
        elseif endswith(output_path, ".fasta") ||
               endswith(output_path, ".fa")
            write(output_path, ">authoritative_consensus\nTTAACCGG\n")
        elseif endswith(output_path, ".sam")
            write(output_path, "@HD\tVN:1.6\n")
        else
            write(output_path, "test artifact\n")
        end
    end
    return nothing
end

function _autocycler_test_graph_only_runner!(step::NamedTuple)::Nothing
    for output_path in step.expected_outputs
        if endswith(output_path, ".gfa")
            mkpath(dirname(output_path))
            write(output_path, "H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\n")
        end
    end
    return nothing
end

function _autocycler_test_missing_output_runner(step::NamedTuple)::Nothing
    return nothing
end

function _autocycler_test_error(function_to_run::Function)::Exception
    try
        function_to_run()
    catch error
        return error
    end
    Base.error("Expected function to throw")
end

Test.@testset "Autocycler wrapper" begin
    Test.@testset "Constants, paths, and bundled environment" begin
        Test.@test Mycelia.AUTOCYCLER_ENV_NAME == "autocycler"
        Test.@test Mycelia.AUTOCYCLER_READ_TYPES == (
            "ont_r9",
            "ont_r10",
            "pacbio_clr",
            "pacbio_hifi",
        )
        Test.@test !isempty(Mycelia.AUTOCYCLER_SCRIPT_URL)
        Test.@test !isempty(Mycelia.AUTOCYCLER_ENV_URL)

        install_dir, script_path, env_file_path = Mycelia._autocycler_paths()
        Test.@test endswith(script_path, "autocycler_full.sh")
        Test.@test endswith(env_file_path, "environment.yml")
        Test.@test endswith(install_dir, joinpath("deps", "autocycler"))

        environment_spec = read(env_file_path, String)
        for dependency in ("bwa", "parallel", "polypolish", "pypolca", "sed")
            Test.@test occursin("  - $(dependency)", environment_spec)
        end
    end

    Test.@testset "Exact long-read-only upstream command" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            out_dir = joinpath(temp_dir, "results")
            script_path = joinpath(temp_dir, "autocycler_full.sh")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(script_path, "#!/usr/bin/env bash\n")

            plan = Mycelia._autocycler_command_plan(
                long_reads,
                out_dir;
                threads = 8,
                jobs = 2,
                read_type = "pacbio_hifi",
                script_path = script_path,
                conda_runner = "/test/conda",
            )
            command = only(plan.steps).command
            Test.@test command.exec == String[
                "/test/conda",
                "run",
                "--live-stream",
                "-n",
                "autocycler",
                "bash",
                abspath(script_path),
                abspath(long_reads),
                "8",
                "2",
                "pacbio_hifi",
            ]
            Test.@test command.dir == abspath(out_dir)
            Test.@test plan.graph == joinpath(
                abspath(out_dir),
                "autocycler_out",
                "consensus_assembly.gfa",
            )
            Test.@test plan.assembly == joinpath(
                abspath(out_dir),
                "autocycler_out",
                "consensus_assembly.fasta",
            )
            Test.@test only(plan.steps).expected_outputs == String[
                plan.graph,
                plan.assembly,
            ]
            Test.@test !ispath(out_dir)
        end
    end

    Test.@testset "Parameter and output validation fails loudly" begin
        invalid_type = _autocycler_test_error() do
            Mycelia._autocycler_command_plan(
                "reads.fastq",
                "out";
                read_type = "illumina",
            )
        end
        Test.@test invalid_type isa ArgumentError
        Test.@test occursin("read_type must be one of", sprint(showerror, invalid_type))
        Test.@test occursin("illumina", sprint(showerror, invalid_type))

        invalid_threads = _autocycler_test_error() do
            Mycelia._autocycler_command_plan(
                "reads.fastq",
                "out";
                threads = 0,
            )
        end
        Test.@test invalid_threads isa ArgumentError
        Test.@test occursin("threads must be positive", sprint(showerror, invalid_threads))

        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            nonempty_out_dir = joinpath(temp_dir, "nonempty")
            mkpath(nonempty_out_dir)
            write(joinpath(nonempty_out_dir, "owned.txt"), "keep\n")

            nonempty_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    nonempty_out_dir;
                    dependency_checker = () -> nothing,
                    runner = _autocycler_test_runner!,
                )
            end
            Test.@test nonempty_error isa ArgumentError
            Test.@test occursin(
                "output directory must be empty",
                sprint(showerror, nonempty_error),
            )
            Test.@test isfile(joinpath(nonempty_out_dir, "owned.txt"))

            missing_output_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "missing-output");
                    dependency_checker = () -> nothing,
                    runner = _autocycler_test_missing_output_runner,
                )
            end
            Test.@test missing_output_error isa ErrorException
            Test.@test occursin(
                "did not create a nonempty artifact",
                sprint(showerror, missing_output_error),
            )
            Test.@test occursin(
                "consensus_assembly.gfa",
                sprint(showerror, missing_output_error),
            )

            missing_fasta_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    long_reads,
                    joinpath(temp_dir, "missing-fasta");
                    dependency_checker = () -> nothing,
                    runner = _autocycler_test_graph_only_runner!,
                )
            end
            Test.@test missing_fasta_error isa ErrorException
            Test.@test occursin(
                "did not create a nonempty artifact",
                sprint(showerror, missing_fasta_error),
            )
            Test.@test occursin(
                "consensus_assembly.fasta",
                sprint(showerror, missing_fasta_error),
            )
        end
    end

    Test.@testset "Autocycler returns authoritative FASTA and GFA" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            dependency_checks = Ref(0)
            result = Mycelia._run_autocycler(
                long_reads,
                joinpath(temp_dir, "autocycler-run");
                threads = 4,
                jobs = 3,
                read_type = "ont_r9",
                dependency_checker = () -> begin
                    dependency_checks[] += 1
                    return nothing
                end,
                runner = _autocycler_test_runner!,
            )
            Test.@test dependency_checks[] == 1
            Test.@test isfile(result.graph)
            Test.@test isfile(result.assembly)
            Test.@test read(result.assembly, String) ==
                       ">authoritative_consensus\nTTAACCGG\n"
            Test.@test read(result.graph, String) ==
                       "H\tVN:Z:1.0\nS\tcontig_1\tACGTACGT\tdp:f:20.0\n"
            Test.@test result.outdir == abspath(joinpath(temp_dir, "autocycler-run"))
        end
    end

    Test.@testset "Paired-short polishing command plan" begin
        mktempdir() do temp_dir
            assembly = joinpath(temp_dir, "draft.fasta")
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            out_dir = joinpath(temp_dir, "results")
            write(assembly, ">contig\nACGT\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")

            plan = Mycelia._autocycler_polishing_command_plan(
                assembly,
                short_reads_1,
                short_reads_2,
                out_dir;
                threads = 6,
                polypolish_careful = true,
                conda_runner = "/test/conda",
            )
            Test.@test Tuple(step.name for step in plan.steps) == (
                :bwa_index,
                :bwa_mem_1,
                :bwa_mem_2,
                :polypolish_filter,
                :polypolish,
                :pypolca,
            )

            bwa_mem_1 = plan.steps[2].command.exec
            bwa_mem_2 = plan.steps[3].command.exec
            Test.@test "-a" in bwa_mem_1
            Test.@test "-a" in bwa_mem_2
            Test.@test bwa_mem_1[end] == abspath(short_reads_1)
            Test.@test bwa_mem_2[end] == abspath(short_reads_2)
            Test.@test !(abspath(short_reads_2) in bwa_mem_1)
            Test.@test !(abspath(short_reads_1) in bwa_mem_2)

            filter_arguments = plan.steps[4].command.exec
            Test.@test "--in1" in filter_arguments
            Test.@test "--in2" in filter_arguments
            Test.@test "--out1" in filter_arguments
            Test.@test "--out2" in filter_arguments
            Test.@test "--careful" in plan.steps[5].command.exec
            Test.@test "--careful" in plan.steps[6].command.exec
            Test.@test endswith(
                plan.assembly,
                joinpath("pypolca", "autocycler_polished_corrected.fasta"),
            )

            noncareful_plan = Mycelia._autocycler_polishing_command_plan(
                assembly,
                short_reads_1,
                short_reads_2,
                out_dir;
                polypolish_careful = false,
            )
            Test.@test !("--careful" in noncareful_plan.steps[5].command.exec)
            Test.@test "--careful" in noncareful_plan.steps[6].command.exec
            Test.@test !ispath(out_dir)
        end
    end

    Test.@testset "Full fake polished workflow preserves every stage" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")
            dependency_checks = Ref(0)

            result = Mycelia._run_autocycler_polished(
                long_reads,
                short_reads_1,
                short_reads_2,
                joinpath(temp_dir, "polished-run");
                threads = 4,
                jobs = 2,
                read_type = "ont_r10",
                dependency_checker = () -> begin
                    dependency_checks[] += 1
                    return nothing
                end,
                runner = _autocycler_test_runner!,
            )
            Test.@test dependency_checks[] == 1
            Test.@test isfile(result.graph)
            Test.@test isfile(result.autocycler_assembly)
            Test.@test isfile(result.polypolish_assembly)
            Test.@test isfile(result.assembly)
            Test.@test isfile(result.pypolca_report)
            Test.@test result.assembly != result.autocycler_assembly
            Test.@test result.assembly != result.polypolish_assembly
        end
    end

    Test.@testset "Paired-short inputs fail before long assembly starts" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            runner_calls = Ref(0)
            runner = (step::NamedTuple) -> begin
                runner_calls[] += 1
                return nothing
            end

            missing_pair_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    joinpath(temp_dir, "missing-R2.fastq"),
                    joinpath(temp_dir, "should-not-run");
                    dependency_checker = () -> nothing,
                    runner = runner,
                )
            end
            Test.@test missing_pair_error isa ArgumentError
            Test.@test occursin(
                "Paired short-read R2 FASTQ not found",
                sprint(showerror, missing_pair_error),
            )
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "should-not-run"))
        end
    end
end

function _autocycler_test_env_enabled(name::AbstractString)::Bool
    value = Base.get(ENV, name, "false")
    return Base.lowercase(Base.strip(value)) == "true"
end

run_all = _autocycler_test_env_enabled("MYCELIA_RUN_ALL")
run_external = run_all || _autocycler_test_env_enabled("MYCELIA_RUN_EXTERNAL")

if run_external
    Test.@testset "Autocycler gated real smoke" begin
        long_reads = get(ENV, "MYCELIA_AUTOCYCLER_LONG_READS", "")
        short_reads_1 = get(ENV, "MYCELIA_AUTOCYCLER_SHORT_READS_1", "")
        short_reads_2 = get(ENV, "MYCELIA_AUTOCYCLER_SHORT_READS_2", "")
        read_type = get(ENV, "MYCELIA_AUTOCYCLER_READ_TYPE", "ont_r10")

        if isempty(long_reads)
            Test.@test_skip "Set MYCELIA_AUTOCYCLER_LONG_READS for a real smoke test"
        elseif xor(isempty(short_reads_1), isempty(short_reads_2))
            Test.@test_skip "Set both Autocycler paired-short read variables or neither"
        else
            Mycelia.install_autocycler()
            mktempdir() do temp_dir
                if isempty(short_reads_1)
                    result = Mycelia.run_autocycler(
                        long_reads = long_reads,
                        out_dir = joinpath(temp_dir, "autocycler"),
                        read_type = read_type,
                    )
                    Test.@test isfile(result.graph)
                    Test.@test filesize(result.graph) > 0
                    Test.@test isfile(result.assembly)
                    Test.@test filesize(result.assembly) > 0
                else
                    result = Mycelia.run_autocycler_polished(
                        long_reads = long_reads,
                        short_reads_1 = short_reads_1,
                        short_reads_2 = short_reads_2,
                        out_dir = joinpath(temp_dir, "autocycler-polished"),
                        read_type = read_type,
                    )
                    Test.@test isfile(result.graph)
                    Test.@test filesize(result.graph) > 0
                    Test.@test isfile(result.assembly)
                    Test.@test filesize(result.assembly) > 0
                end
            end
        end
    end
else
    @info "Skipping Autocycler real smoke; set MYCELIA_RUN_EXTERNAL=true to enable"
end
