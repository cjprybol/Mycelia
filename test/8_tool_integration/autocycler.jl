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

function _autocycler_test_compatible_packages()::Dict{String, String}
    return Dict(
        specification.name =>
            (specification.constraint == :present ? "present" : specification.version)
        for specification in Mycelia.AUTOCYCLER_REQUIRED_PACKAGE_SPECS
    )
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
        Test.@test occursin(
            Mycelia.AUTOCYCLER_SCRIPT_REVISION,
            Mycelia.AUTOCYCLER_SCRIPT_URL,
        )
        Test.@test length(Mycelia.AUTOCYCLER_SCRIPT_SHA256) == 64
        install_dir, script_path, env_file_path = Mycelia._autocycler_paths()
        Test.@test endswith(script_path, "autocycler_full.sh")
        Test.@test endswith(env_file_path, "environment.yml")
        Test.@test endswith(install_dir, joinpath("deps", "autocycler"))

        environment_spec = read(env_file_path, String)
        Test.@test occursin("  - autocycler=0.5.2", environment_spec)
        dependency_section = split(environment_spec, "dependencies:"; limit = 2)[2]
        dependency_names = Set{String}()
        for line in split(dependency_section, '\n')
            dependency_match = match(r"^  - ([A-Za-z0-9_-]+)", line)
            dependency_match === nothing || push!(
                dependency_names,
                something(only(dependency_match.captures)),
            )
        end
        Test.@test dependency_names == Set(Mycelia.AUTOCYCLER_REQUIRED_PACKAGES)
        for specification in Mycelia.AUTOCYCLER_REQUIRED_PACKAGE_SPECS
            constraint = if specification.constraint == :exact
                "="
            elseif specification.constraint == :minimum
                ">="
            else
                ""
            end
            Test.@test occursin(
                "  - $(specification.name)$(constraint)$(specification.version)",
                environment_spec,
            )
        end
        Test.@test Mycelia._autocycler_script_is_verified(script_path)
        compatible_packages = _autocycler_test_compatible_packages()
        toolchain = Mycelia._autocycler_toolchain_metadata(compatible_packages)
        Test.@test Set(keys(toolchain)) == Set([
            "autocycler_script_revision",
            "autocycler_script_sha256",
            "environment_spec_sha256",
            "packages",
        ])
        Test.@test toolchain["autocycler_script_revision"] ==
                   Mycelia.AUTOCYCLER_SCRIPT_REVISION
        Test.@test toolchain["autocycler_script_sha256"] ==
                   Mycelia.AUTOCYCLER_SCRIPT_SHA256
        Test.@test toolchain["environment_spec_sha256"] ==
                   Mycelia._autocycler_sha256(env_file_path)
        Test.@test toolchain["packages"] == compatible_packages
    end

    Test.@testset "Pinned script and stale-environment preflight" begin
        mktempdir() do temp_dir
            script_path = joinpath(temp_dir, "autocycler_full.sh")
            downloader = function (
                    url::AbstractString,
                    destination::AbstractString,
            )
                Test.@test url == Mycelia.AUTOCYCLER_SCRIPT_URL
                write(destination, "tampered script\n")
                return String(destination)
            end
            checksum_error = _autocycler_test_error() do
                Mycelia._install_verified_autocycler_script!(
                    script_path;
                    downloader,
                )
            end
            Test.@test checksum_error isa ErrorException
            Test.@test occursin(
                "checksum mismatch",
                sprint(showerror, checksum_error),
            )
            Test.@test !isfile(script_path)
        end

        mktempdir() do temp_dir
            installed_script = joinpath(temp_dir, "autocycler_full.sh")
            bundled_script = Mycelia._autocycler_paths()[2]
            downloader = function (
                    url::AbstractString,
                    destination::AbstractString,
            )
                Test.@test url == Mycelia.AUTOCYCLER_SCRIPT_URL
                cp(bundled_script, destination; force = true)
                return String(destination)
            end
            Test.@test Mycelia._install_verified_autocycler_script!(
                installed_script;
                downloader,
            ) == abspath(installed_script)
            Test.@test Mycelia._autocycler_script_is_verified(installed_script)
        end

        parsed_versions = Mycelia._autocycler_environment_packages(;
            conda_runner = "/test/conda",
            command_reader = command -> begin
                Test.@test command.exec == String[
                    "/test/conda",
                    "list",
                    "-n",
                    "autocycler",
                    "--json",
                ]
                return "[{\"name\":\"autocycler\",\"version\":\"0.5.2\"}]"
            end,
        )
        Test.@test parsed_versions == Dict("autocycler" => "0.5.2")

        mktempdir() do temp_dir
            paths = (
                temp_dir,
                joinpath(temp_dir, "autocycler_full.sh"),
                joinpath(temp_dir, "environment.yml"),
            )
            missing_environment_error = _autocycler_test_error() do
                Mycelia.install_autocycler(; paths)
            end
            Test.@test missing_environment_error isa ErrorException
            Test.@test occursin(
                "environment file is missing or empty",
                sprint(showerror, missing_environment_error),
            )
            write(paths[3], "")
            empty_environment_error = _autocycler_test_error() do
                Mycelia.install_autocycler(; paths)
            end
            Test.@test empty_environment_error isa ErrorException
        end

        inspection_calls = Ref(0)
        package_inspector = function ()
            inspection_calls[] += 1
            packages = _autocycler_test_compatible_packages()
            if inspection_calls[] == 1
                packages["autocycler"] = "0.5.1"
                delete!(packages, "canu")
            end
            return packages
        end
        installer_forces = Bool[]
        installer = function (; force::Bool = false)
            push!(installer_forces, force)
            return nothing
        end
        versions = Mycelia._ensure_autocycler_packages!(
            package_inspector,
            installer,
        )
        Test.@test inspection_calls[] == 2
        Test.@test installer_forces == [true]
        Test.@test versions["autocycler"] == "0.5.2"
        Test.@test haskey(versions, "canu")

        too_old = _autocycler_test_compatible_packages()
        too_old["flye"] = "2.9.5"
        Test.@test any(
            issue -> occursin("flye must be at least 2.9.6", issue),
            Mycelia._autocycler_package_issues(too_old),
        )
        necat_without_update = _autocycler_test_compatible_packages()
        necat_without_update["necat"] = "0.0.1"
        Test.@test any(
            issue -> occursin("necat must be at least 0.0.1_update20200803", issue),
            Mycelia._autocycler_package_issues(necat_without_update),
        )
        Test.@test Mycelia._autocycler_version_at_least(
            "0.0.1_update20210101",
            "0.0.1_update20200803",
        )

        always_stale = () -> Dict(
            "autocycler" => "0.5.2",
            "bwa" => "0.7.17",
        )
        stale_error = _autocycler_test_error() do
            Mycelia._ensure_autocycler_packages!(always_stale, installer)
        end
        Test.@test stale_error isa ErrorException
        Test.@test occursin(
            "missing or incompatible required packages",
            sprint(showerror, stale_error),
        )
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

    Test.@testset "Step execution preserves exception semantics" begin
        step = (
            name = :synthetic,
            command = `true`,
            stdout = nothing,
            expected_outputs = String[],
        )

        interrupted = _autocycler_test_error() do
            Mycelia._execute_autocycler_steps(
                (step,);
                runner = _ -> throw(InterruptException()),
            )
        end
        Test.@test interrupted isa InterruptException

        runner_bug = ErrorException("synthetic runner bug")
        preserved_bug = _autocycler_test_error() do
            Mycelia._execute_autocycler_steps(
                (step,);
                runner = _ -> throw(runner_bug),
            )
        end
        Test.@test preserved_bug === runner_bug

        process_error = _autocycler_test_error() do
            Mycelia._execute_autocycler_steps(
                (step,);
                runner = _ -> Base.run(`false`),
            )
        end
        Test.@test process_error isa ErrorException
        Test.@test occursin(
            "workflow command synthetic failed",
            sprint(showerror, process_error),
        )
    end

    Test.@testset "Autocycler returns authoritative FASTA and GFA" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            dependency_checks = Ref(0)
            expected_toolchain = Dict(
                "autocycler_script_revision" => "test-revision",
            )
            result = Mycelia._run_autocycler(
                long_reads,
                joinpath(temp_dir, "autocycler-run");
                threads = 4,
                jobs = 3,
                read_type = "ont_r9",
                dependency_checker = () -> begin
                    dependency_checks[] += 1
                    return expected_toolchain
                end,
                runner = _autocycler_test_runner!,
            )
            Test.@test dependency_checks[] == 1
            Test.@test result.toolchain == expected_toolchain
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
            Test.@test plan.steps[1].expected_outputs == String[
                "$(abspath(assembly)).$(extension)" for
                extension in ("amb", "ann", "bwt", "pac", "sa")
            ]
            Test.@test plan.bwa_index_files == plan.steps[1].expected_outputs

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
            Test.@test result.toolchain === nothing
            Test.@test isempty(result.intermediates)
            Test.@test !any(isfile, String[
                "$(result.autocycler_assembly).$(extension)" for
                extension in ("amb", "ann", "bwt", "pac", "sa")
            ])

            retained = Mycelia._run_autocycler_polished(
                long_reads,
                short_reads_1,
                short_reads_2,
                joinpath(temp_dir, "retained-polishing-run");
                dependency_checker = () -> nothing,
                runner = _autocycler_test_runner!,
                keep_intermediates = true,
            )
            Test.@test length(retained.intermediates) == 9
            Test.@test all(isfile, retained.intermediates)
        end
    end

    Test.@testset "Failed polishing cleans only route-owned intermediates" begin
        mktempdir() do temp_dir
            long_reads = joinpath(temp_dir, "long.fastq")
            short_reads_1 = joinpath(temp_dir, "R1.fastq")
            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            write(long_reads, "@long\nACGT\n+\nIIII\n")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")

            for keep_intermediates in (false, true)
                out_dir = joinpath(temp_dir, "failure-$(keep_intermediates)")
                marker_error = ErrorException("synthetic pypolca failure")
                runner = function (step::NamedTuple)
                    _autocycler_test_runner!(step)
                    step.name == :pypolca && throw(marker_error)
                    return nothing
                end
                observed_error = _autocycler_test_error() do
                    Mycelia._run_autocycler_polished(
                        long_reads,
                        short_reads_1,
                        short_reads_2,
                        out_dir;
                        keep_intermediates,
                        dependency_checker = () -> nothing,
                        runner,
                    )
                end
                Test.@test observed_error === marker_error

                autocycler_assembly = joinpath(
                    out_dir,
                    "autocycler_out",
                    "consensus_assembly.fasta",
                )
                polishing_plan = Mycelia._autocycler_polishing_command_plan(
                    autocycler_assembly,
                    short_reads_1,
                    short_reads_2,
                    out_dir,
                )
                Test.@test all(
                    path -> isfile(path) == keep_intermediates,
                    polishing_plan.intermediate_files,
                )
                Test.@test isfile(autocycler_assembly)
                Test.@test isfile(joinpath(
                    out_dir,
                    "autocycler_out",
                    "consensus_assembly.gfa",
                ))
                Test.@test isfile(polishing_plan.polypolish_assembly)
                Test.@test isfile(polishing_plan.assembly)
                Test.@test isfile(polishing_plan.pypolca_report)
            end
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

            short_reads_2 = joinpath(temp_dir, "R2.fastq")
            write(short_reads_2, "@different/2\nACGT\n+\nIIII\n")
            dependency_checks = Ref(0)
            mismatch_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "mismatch-should-not-run");
                    dependency_checker = () -> begin
                        dependency_checks[] += 1
                        return nothing
                    end,
                    runner = runner,
                )
            end
            Test.@test mismatch_error isa ArgumentError
            Test.@test occursin("out of sync", sprint(showerror, mismatch_error))
            Test.@test dependency_checks[] == 0
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "mismatch-should-not-run"))

            same_file_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_1,
                    joinpath(temp_dir, "same-file-should-not-run");
                    dependency_checker = () -> nothing,
                    runner = runner,
                )
            end
            Test.@test same_file_error isa ArgumentError
            Test.@test occursin("must be distinct files", sprint(showerror, same_file_error))
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "same-file-should-not-run"))

            write(short_reads_1, "@pair/2\nACGT\n+\nIIII\n")
            write(short_reads_2, "@pair/1\nACGT\n+\nIIII\n")
            reversed_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "reversed-should-not-run");
                    dependency_checker = () -> nothing,
                    runner = runner,
                )
            end
            Test.@test reversed_error isa ArgumentError
            Test.@test occursin(
                "invalid explicit mate roles",
                sprint(showerror, reversed_error),
            )
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "reversed-should-not-run"))

            fasta_pair = joinpath(temp_dir, "not-fastq.fasta")
            write(short_reads_1, "@pair/1\nACGT\n+\nIIII\n")
            write(fasta_pair, ">pair/2\nACGT\n")
            fasta_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    fasta_pair,
                    joinpath(temp_dir, "fasta-should-not-run");
                    dependency_checker = () -> nothing,
                    runner = runner,
                )
            end
            Test.@test fasta_error isa ArgumentError
            Test.@test occursin("must be FASTQ", sprint(showerror, fasta_error))
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "fasta-should-not-run"))

            fasta_long = joinpath(temp_dir, "not-long-fastq.fasta")
            write(fasta_long, ">long\nACGT\n")
            long_dependency_checks = Ref(0)
            long_fasta_error = _autocycler_test_error() do
                Mycelia._run_autocycler(
                    fasta_long,
                    joinpath(temp_dir, "long-fasta-should-not-run");
                    dependency_checker = () -> begin
                        long_dependency_checks[] += 1
                        return nothing
                    end,
                    runner = runner,
                )
            end
            Test.@test long_fasta_error isa ArgumentError
            Test.@test occursin("must be a FASTQ", sprint(showerror, long_fasta_error))
            Test.@test long_dependency_checks[] == 0
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "long-fasta-should-not-run"))

            write(short_reads_2, "@pair/2\nACGT\n+\nIIII\n")
            stale_environment_error = _autocycler_test_error() do
                Mycelia._run_autocycler_polished(
                    long_reads,
                    short_reads_1,
                    short_reads_2,
                    joinpath(temp_dir, "stale-env-should-not-run");
                    dependency_checker = () -> throw(
                        ErrorException("missing required packages: pypolca"),
                    ),
                    runner = runner,
                )
            end
            Test.@test stale_environment_error isa ErrorException
            Test.@test occursin(
                "missing required packages",
                sprint(showerror, stale_environment_error),
            )
            Test.@test runner_calls[] == 0
            Test.@test !ispath(joinpath(temp_dir, "stale-env-should-not-run"))
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
