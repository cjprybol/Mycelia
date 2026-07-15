# Default-CI contract tests for the single-technology metaMDBG wrapper.

import Test
import Mycelia

function _test_metamdbg_error(
        thunk::Function,
        expected_type::Type{<:Exception},
        expected_message::Regex,
)::Nothing
    exception = try
        thunk()
        nothing
    catch caught
        caught
    end
    Test.@test exception isa expected_type
    if exception isa Exception
        Test.@test occursin(expected_message, sprint(showerror, exception))
    end
    return nothing
end

function _write_test_metamdbg_contract!(
        outdir::AbstractString,
        hifi_reads::AbstractString;
        abundance_min::Int = 3,
        graph_k::Int = 21,
)::String
    selected_input = Mycelia._metamdbg_selected_input(
        String(hifi_reads),
        nothing,
    )
    input_contract = Mycelia._metamdbg_input_contract(
        selected_input,
        abundance_min,
    )
    outputs = Mycelia._metamdbg_output_paths(outdir, graph_k)
    return Mycelia._write_metamdbg_contract!(outputs, input_contract)
end

Test.@testset "metaMDBG input and artifact contracts" begin
    mktempdir() do temporary_root
        valid_reads = joinpath(temporary_root, "reads.fastq")
        empty_reads = joinpath(temporary_root, "empty.fastq")
        write(valid_reads, "@read-1\nACGT\n+\nIIII\n")
        touch(empty_reads)

        provisioning_calls = Ref(0)
        dependency_checker = () -> begin
            provisioning_calls[] += 1
            return nothing
        end
        forbidden_runner = function (_command::Cmd)
            error("metaMDBG command must not run")
        end

        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                outdir = joinpath(temporary_root, "neither"),
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ArgumentError,
            r"exactly one input technology",
        )
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                ont_reads = valid_reads,
                outdir = joinpath(temporary_root, "both"),
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ArgumentError,
            r"exactly one input technology",
        )
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = String[],
                outdir = joinpath(temporary_root, "empty-vector"),
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ArgumentError,
            r"must contain at least one non-empty path",
        )
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = joinpath(temporary_root, "missing.fastq"),
                outdir = joinpath(temporary_root, "missing-file"),
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ArgumentError,
            r"input file does not exist",
        )
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = empty_reads,
                outdir = joinpath(temporary_root, "empty-file"),
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ArgumentError,
            r"input file is empty",
        )
        Test.@test provisioning_calls[] == 0
        Test.@test !ispath(joinpath(temporary_root, "missing-file"))

        unprovenanced_outdir = joinpath(temporary_root, "unprovenanced")
        mkpath(unprovenanced_outdir)
        write(
            joinpath(unprovenanced_outdir, "contigs.fasta.gz"),
            "untrusted-output",
        )
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir = unprovenanced_outdir,
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ErrorException,
            r"without its provenance contract marker",
        )
        Test.@test provisioning_calls[] == 0

        Test.@testset "local plain-output normalization and idempotent reuse" begin
            outdir = joinpath(temporary_root, "local")
            command_count = Ref(0)
            local_runner = function (_command::Cmd)
                command_count[] += 1
                if command_count[] == 1
                    write(
                        joinpath(outdir, "contigs.fasta"),
                        ">contig-1\nACGT\n",
                    )
                elseif command_count[] == 2
                    write(
                        joinpath(outdir, "assemblyGraph_k21_4bps.gfa"),
                        "H\tVN:Z:1.0\n",
                    )
                else
                    error("unexpected metaMDBG command")
                end
                return nothing
            end
            local_provisioning_calls = Ref(0)
            result = Mycelia._run_metamdbg(;
                ont_reads = valid_reads,
                outdir,
                graph_k = 21,
                dependency_checker = () -> begin
                    local_provisioning_calls[] += 1
                    return nothing
                end,
                local_runner,
            )

            Test.@test command_count[] == 2
            Test.@test local_provisioning_calls[] == 1
            Test.@test result.contigs == joinpath(outdir, "contigs.fasta.gz")
            Test.@test isfile(result.contigs)
            Test.@test filesize(result.contigs) > 0
            Test.@test isfile(joinpath(outdir, "contigs.fasta"))
            Test.@test result.graph == joinpath(outdir, "assemblyGraph_k21.gfa")
            Test.@test islink(result.graph)
            Test.@test read(result.graph, String) == "H\tVN:Z:1.0\n"
            outputs = Mycelia._metamdbg_output_paths(outdir, 21)
            expected_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(nothing, valid_reads),
                3,
            )
            Test.@test isfile(outputs.contract_marker)
            Test.@test read(outputs.contract_marker, String) ==
                       expected_contract.contents

            reader = Mycelia.open_fastx(result.contigs)
            records = try
                collect(reader)
            finally
                close(reader)
            end
            Test.@test length(records) == 1

            reused = Mycelia._run_metamdbg(;
                ont_reads = valid_reads,
                outdir,
                graph_k = 21,
                dependency_checker = () -> error(
                    "complete metaMDBG reuse must not provision",
                ),
                local_runner = forbidden_runner,
            )
            Test.@test reused == result
            Test.@test command_count[] == 2

            changed_reads = joinpath(temporary_root, "changed-reads.fastq")
            write(changed_reads, "@read-1\nACGT\n+\nIIII\n")
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = changed_reads,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> begin
                        local_provisioning_calls[] += 1
                        return nothing
                    end,
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"existing output contract does not match",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> begin
                        local_provisioning_calls[] += 1
                        return nothing
                    end,
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"existing output contract does not match",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    outdir,
                    abundance_min = 4,
                    graph_k = 21,
                    dependency_checker = () -> begin
                        local_provisioning_calls[] += 1
                        return nothing
                    end,
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"existing output contract does not match",
            )
            open(valid_reads, "a") do io
                write(io, "@read-2\nTGCA\n+\nIIII\n")
            end
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> begin
                        local_provisioning_calls[] += 1
                        return nothing
                    end,
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"existing output contract does not match",
            )
            Test.@test local_provisioning_calls[] == 1
            Test.@test command_count[] == 2
        end

        Test.@testset "exact graph-k discovery and ambiguity failure" begin
            outdir = joinpath(temporary_root, "exact-k")
            mkpath(outdir)
            write(joinpath(outdir, "contigs.fasta.gz"), "compressed-placeholder")
            write(joinpath(outdir, "assemblyGraph_k20_4bps.gfa"), "wrong-k\n")
            _write_test_metamdbg_contract!(outdir, valid_reads)
            graph_calls = Ref(0)
            result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir,
                graph_k = 21,
                dependency_checker = () -> nothing,
                local_runner = function (command::Cmd)
                    graph_calls[] += 1
                    Test.@test occursin(" gfa ", string(command))
                    write(
                        joinpath(outdir, "assemblyGraph_k21_4bps.gfa"),
                        "right-k\n",
                    )
                    return nothing
                end,
            )
            Test.@test graph_calls[] == 1
            Test.@test read(result.graph, String) == "right-k\n"

            ambiguous_outdir = joinpath(temporary_root, "ambiguous")
            mkpath(ambiguous_outdir)
            write(
                joinpath(ambiguous_outdir, "contigs.fasta.gz"),
                "compressed-placeholder",
            )
            write(
                joinpath(ambiguous_outdir, "assemblyGraph_k21_4bps.gfa"),
                "first\n",
            )
            write(
                joinpath(ambiguous_outdir, "assemblyGraph_k21_5bps.gfa"),
                "second\n",
            )
            _write_test_metamdbg_contract!(ambiguous_outdir, valid_reads)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = ambiguous_outdir,
                    graph_k = 21,
                    dependency_checker = () -> error(
                        "ambiguous reuse must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"multiple nonempty graph artifacts for k=21",
            )
        end

        Test.@testset "executor script validates dynamic artifacts" begin
            outdir = joinpath(temporary_root, "executor output")
            outputs = Mycelia._metamdbg_output_paths(outdir, 21)
            input_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            script = Mycelia._metamdbg_executor_script(
                "metaMDBG asm --fixture",
                "metaMDBG gfa --fixture",
                outputs,
                21,
                input_contract,
            )
            Test.@test occursin("assemblyGraph_k21_*bps.gfa", script)
            Test.@test occursin("metaMDBG asm --fixture", script)
            Test.@test occursin("metaMDBG gfa --fixture", script)
            Test.@test occursin("gzip -c", script)
            Test.@test occursin("find_metamdbg_graph", script)
            Test.@test occursin("ln -s", script)
            Test.@test occursin("test -s \"\$contigs_gz\"", script)
            Test.@test occursin("test -s \"\$graph_alias\"", script)
            Test.@test occursin("mycelia_metamdbg_contract.json", script)
            Test.@test occursin("cmp -s", script)
            Test.@test occursin("mv -f", script)
            Test.@test occursin(input_contract.signature, script)
        end

        Test.@testset "collected executor reports normalized planned paths" begin
            outdir = joinpath(temporary_root, "collected")
            executor = Mycelia.CollectExecutor()
            result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir,
                graph_k = 31,
                executor,
                dependency_checker = () -> nothing,
                local_runner = forbidden_runner,
            )
            Test.@test length(executor.jobs) == 1
            Test.@test result.contigs == joinpath(outdir, "contigs.fasta.gz")
            Test.@test result.graph == joinpath(outdir, "assemblyGraph_k31.gfa")
            Test.@test occursin(
                "assemblyGraph_k31_*bps.gfa",
                only(executor.jobs).cmd,
            )
            Test.@test occursin("test -s \"\$contigs_gz\"", only(executor.jobs).cmd)
            Test.@test occursin("test -s \"\$graph_alias\"", only(executor.jobs).cmd)
            Test.@test occursin(
                "mycelia_metamdbg_contract.json",
                only(executor.jobs).cmd,
            )
            Test.@test occursin("cmp -s", only(executor.jobs).cmd)
        end
    end
end
