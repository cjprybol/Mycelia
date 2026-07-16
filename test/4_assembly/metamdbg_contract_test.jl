# Default-CI contract tests for the single-technology metaMDBG wrapper.

import CodecZlib
import JSON
import Test
import Mycelia

struct _TestUnverifiableMetamdbgExecutor <: Mycelia.AbstractExecutor end

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

function _test_metamdbg_toolchain()::NamedTuple
    return Mycelia._require_metamdbg_package_version(NamedTuple[
        (;
            name = "libzlib",
            version = "1.3.1",
            build = "h8359307_2",
            channel = "conda-forge",
        ),
        (;
            name = "metamdbg",
            version = "1.4",
            build = "h43eeafb_2",
            channel = "bioconda",
        ),
    ])
end

function _write_test_metamdbg_gzip_fasta!(
        path::AbstractString;
        contents::AbstractString = ">contig-1\nACGT\n",
)::String
    open(path, "w") do raw_output
        gzip_output = CodecZlib.GzipCompressorStream(raw_output)
        try
            write(gzip_output, contents)
        finally
            close(gzip_output)
        end
    end
    return String(path)
end

function _write_test_metamdbg_gfa!(
        path::AbstractString;
        sequence::AbstractString = "ACGT",
)::String
    write(path, "H\tVN:Z:1.0\nS\tcontig-1\t$(sequence)\n")
    return String(path)
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
        temporary_root = realpath(temporary_root)
        valid_reads = joinpath(temporary_root, "reads.fastq")
        empty_reads = joinpath(temporary_root, "empty.fastq")
        write(valid_reads, "@read-1\nACGT\n+\nIIII\n")
        touch(empty_reads)
        symlink_reads = joinpath(temporary_root, "reads-symlink.fastq")
        hardlink_reads = joinpath(temporary_root, "reads-hardlink.fastq")
        symlink(valid_reads, symlink_reads)
        Base.Filesystem.hardlink(valid_reads, hardlink_reads)

        provisioning_calls = Ref(0)
        dependency_checker = () -> begin
            provisioning_calls[] += 1
            return _test_metamdbg_toolchain()
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
        aliased_inputs = (
            exact = String[valid_reads, valid_reads],
            symlink = String[valid_reads, symlink_reads],
            hardlink = String[valid_reads, hardlink_reads],
        )
        for (alias_type, input_paths) in pairs(aliased_inputs)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = input_paths,
                    outdir = joinpath(temporary_root, "$(alias_type)-alias"),
                    dependency_checker,
                    local_runner = forbidden_runner,
                ),
                ArgumentError,
                r"must refer to physically distinct files",
            )
            Test.@test !ispath(joinpath(
                temporary_root,
                "$(alias_type)-alias",
            ))
        end
        distinct_reads = joinpath(temporary_root, "reads-distinct.fastq")
        write(distinct_reads, "@read-2\nTGCA\n+\nJJJJ\n")
        ordered_paths = String[distinct_reads, valid_reads]
        ordered_input =
            Mycelia._metamdbg_selected_input(ordered_paths, nothing)
        ordered_contract = Mycelia._metamdbg_input_contract(ordered_input, 3)
        Test.@test getproperty.(
            ordered_contract.contract.inputs,
            :path,
        ) == ordered_paths
        Test.@test getproperty.(
            ordered_contract.contract.inputs,
            :sha256,
        ) == Mycelia._metamdbg_sha256.(ordered_paths)
        reversed_contract = Mycelia._metamdbg_input_contract(
            Mycelia._metamdbg_selected_input(reverse(ordered_paths), nothing),
            3,
        )
        Test.@test ordered_contract.signature != reversed_contract.signature
        distinct_executor = Mycelia.CollectExecutor()
        distinct_result = Mycelia._run_metamdbg(;
            hifi_reads = ordered_paths,
            outdir = joinpath(temporary_root, "distinct-input-plan"),
            executor = distinct_executor,
            dependency_checker = () -> error(
                "planned distinct inputs must not provision metaMDBG",
            ),
            local_runner = forbidden_runner,
        )
        Test.@test distinct_result.status == :planned
        distinct_script = only(distinct_executor.jobs).cmd
        distinct_positions = map(ordered_paths) do path
            position = findfirst(path, distinct_script)
            Test.@test position !== nothing
            return first(position)
        end
        Test.@test first(distinct_positions) < last(distinct_positions)
        for fingerprint in getproperty.(
                ordered_contract.contract.inputs,
                :sha256,
        )
            Test.@test occursin(fingerprint, distinct_script)
        end
        distinct_execution_outdir =
            joinpath(temporary_root, "distinct-input-execution")
        distinct_command_count = Ref(0)
        distinct_runner = function (command::Cmd)
            distinct_command_count[] += 1
            if distinct_command_count[] == 1
                input_flag_index = only(findall(
                    ==("--in-hifi"),
                    command.exec,
                ))
                abundance_index = only(findall(
                    ==("--min-abundance"),
                    command.exec,
                ))
                staged_paths = command.exec[
                    (input_flag_index + 1):(abundance_index - 1)
                ]
                Test.@test length(staged_paths) == 2
                Test.@test read.(staged_paths, String) ==
                           read.(ordered_paths, String)
                write(
                    joinpath(distinct_execution_outdir, "contigs.fasta"),
                    ">distinct-contig\nACGT\n",
                )
            elseif distinct_command_count[] == 2
                _write_test_metamdbg_gfa!(joinpath(
                    distinct_execution_outdir,
                    "assemblyGraph_k21_4bps.gfa",
                ))
            else
                error("unexpected distinct-input metaMDBG command")
            end
            return nothing
        end
        distinct_execution = Mycelia._run_metamdbg(;
            hifi_reads = ordered_paths,
            outdir = distinct_execution_outdir,
            dependency_checker = _test_metamdbg_toolchain,
            local_runner = distinct_runner,
        )
        Test.@test distinct_execution.status == :complete
        Test.@test distinct_command_count[] == 2
        executed_contract = JSON.parse(read(
            distinct_execution.contract_marker,
            String,
        ))
        Test.@test getindex.(
            executed_contract["contract"]["inputs"],
            "path",
        ) == ordered_paths
        Test.@test getindex.(
            executed_contract["contract"]["inputs"],
            "sha256",
        ) == Mycelia._metamdbg_sha256.(ordered_paths)
        Test.@test provisioning_calls[] == 0
        Test.@test !ispath(joinpath(temporary_root, "missing-file"))

        unsupported_outdir = joinpath(temporary_root, "unsupported-executor")
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir = unsupported_outdir,
                executor = _TestUnverifiableMetamdbgExecutor(),
                dependency_checker,
                local_runner = forbidden_runner,
                submission_runner = (_job, _executor) -> error(
                    "unsupported executor must fail before submission",
                ),
            ),
            ArgumentError,
            r"refusing an unverifiable custom nonlocal executor",
        )
        Test.@test !ispath(unsupported_outdir)
        Test.@test isempty(
            Mycelia._metamdbg_submission_reservation_paths(unsupported_outdir),
        )

        unprovenanced_outdir = joinpath(temporary_root, "unprovenanced")
        mkpath(unprovenanced_outdir)
        write(
            joinpath(unprovenanced_outdir, "arbitrary.partial"),
            "legacy-or-partial-output",
        )
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir = unprovenanced_outdir,
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ErrorException,
            r"nonempty output root without its regular provenance contract marker",
        )
        Test.@test provisioning_calls[] == 0

        Test.@testset "local plain-output normalization and idempotent reuse" begin
            outdir = joinpath(temporary_root, "local")
            output_lock_path = Mycelia._metamdbg_output_lock_path(outdir)
            command_count = Ref(0)
            local_runner = function (command::Cmd)
                command_count[] += 1
                Test.@test isdir(output_lock_path)
                Test.@test !ispath(joinpath(
                    outdir,
                    "mycelia_metamdbg_contract.json",
                ))
                if command_count[] == 1
                    input_flag_index = only(findall(
                        ==("--in-ont"),
                        command.exec,
                    ))
                    staged_input = command.exec[input_flag_index + 1]
                    Test.@test staged_input != valid_reads
                    Test.@test occursin(
                        ".mycelia-metamdbg-inputs.",
                        staged_input,
                    )
                    Test.@test read(staged_input, String) ==
                               read(valid_reads, String)
                    Test.@test (stat(staged_input).mode & 0o777) == 0o400
                    write(
                        joinpath(outdir, "contigs.fasta"),
                        ">contig-1\nACGT\n",
                    )
                elseif command_count[] == 2
                    _write_test_metamdbg_gfa!(
                        joinpath(outdir, "assemblyGraph_k21_4bps.gfa"),
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
                    Test.@test isdir(output_lock_path)
                    return _test_metamdbg_toolchain()
                end,
                local_runner,
            )

            Test.@test command_count[] == 2
            Test.@test local_provisioning_calls[] == 2
            Test.@test result.status == :complete
            Test.@test result.contigs == joinpath(outdir, "contigs.fasta.gz")
            Test.@test isfile(result.contigs)
            Test.@test filesize(result.contigs) > 0
            Test.@test isfile(joinpath(outdir, "contigs.fasta"))
            Test.@test result.graph == joinpath(outdir, "assemblyGraph_k21.gfa")
            Test.@test islink(result.graph)
            Test.@test read(result.graph, String) ==
                       "H\tVN:Z:1.0\nS\tcontig-1\tACGT\n"
            outputs = Mycelia._metamdbg_output_paths(outdir, 21)
            expected_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(nothing, valid_reads),
                3,
            )
            Test.@test isfile(outputs.contract_marker)
            Test.@test read(outputs.contract_marker, String) ==
                       expected_contract.contents
            expected_input = only(expected_contract.contract.inputs)
            Test.@test expected_contract.contract.schema_version == 4
            Test.@test expected_input.sha256 ==
                       Mycelia._metamdbg_sha256(valid_reads)
            Test.@test !hasproperty(expected_input, :modification_time_ns)
            Test.@test !occursin(
                "modification_time_ns",
                expected_contract.contents,
            )
            Test.@test result.provenance.contract_signature ==
                       expected_contract.signature
            Test.@test result.provenance.metamdbg_version == "1.4"
            Test.@test result.provenance.environment_name ==
                       Mycelia.METAMDBG_ENV_NAME
            Test.@test result.provenance.graph_k == 21
            Test.@test occursin(
                r"^[0-9a-f]{64}$",
                result.provenance.workflow_signature,
            )
            Test.@test result.provenance.package_inventory_sha256 ==
                       _test_metamdbg_toolchain().package_inventory_sha256
            Test.@test result.provenance.package_count == 2
            Test.@test result.provenance.contigs_size_bytes ==
                       filesize(result.contigs)
            Test.@test result.provenance.graph_size_bytes ==
                       filesize(realpath(result.graph))
            Test.@test result.completion_marker == outputs.completion_marker
            Test.@test isfile(outputs.completion_marker)
            completion_record = JSON.parse(read(
                outputs.completion_marker,
                String,
            ))
            Test.@test completion_record["manifest"]["workflow"]["graph_k"] ==
                       21
            Test.@test completion_record["manifest"]["artifacts"]["contigs"]["sha256"] ==
                       Mycelia._metamdbg_sha256(result.contigs)
            Test.@test completion_record["manifest"]["artifacts"]["graph"]["sha256"] ==
                       Mycelia._metamdbg_sha256(realpath(result.graph))
            Test.@test !ispath(output_lock_path)

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

            contigs_backup = joinpath(temporary_root, "contigs-backup.fasta.gz")
            cp(result.contigs, contigs_backup)
            _write_test_metamdbg_gzip_fasta!(
                result.contigs;
                contents = ">contig-1\nTGCA\n",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> error(
                        "artifact-replaced reuse must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"completion manifest does not match",
            )
            cp(contigs_backup, result.contigs; force = true)

            graph_source = realpath(result.graph)
            graph_backup = read(graph_source, String)
            _write_test_metamdbg_gfa!(graph_source; sequence = "TGCA")
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> error(
                        "graph-replaced reuse must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"completion manifest does not match",
            )
            write(graph_source, graph_backup)

            changed_reads = joinpath(temporary_root, "changed-reads.fastq")
            write(changed_reads, "@read-1\nACGT\n+\nIIII\n")
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = changed_reads,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> begin
                        local_provisioning_calls[] += 1
                        return _test_metamdbg_toolchain()
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
                        return _test_metamdbg_toolchain()
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
                        return _test_metamdbg_toolchain()
                    end,
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"existing output contract does not match",
            )

            metadata_reference = joinpath(
                temporary_root,
                "reads-original-metadata.fastq",
            )
            write(metadata_reference, read(valid_reads, String))
            run(`touch -r $(valid_reads) $(metadata_reference)`)
            write(valid_reads, "@read-1\nTGCA\n+\nIIII\n")
            run(`touch -r $(metadata_reference) $(valid_reads)`)
            changed_content_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(nothing, valid_reads),
                3,
            )
            changed_content_input = only(changed_content_contract.contract.inputs)
            expected_snapshot = only(expected_contract.input_snapshot)
            changed_content_snapshot =
                only(changed_content_contract.input_snapshot)
            Test.@test changed_content_input.size_bytes == expected_input.size_bytes
            Test.@test changed_content_snapshot.modification_time_ns ==
                       expected_snapshot.modification_time_ns
            Test.@test changed_content_input.sha256 != expected_input.sha256
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> error(
                        "content-changed reuse must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"existing output contract does not match",
            )
            write(valid_reads, "@read-1\nACGT\n+\nIIII\n")
            run(`touch -r $(metadata_reference) $(valid_reads)`)
            restored_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(nothing, valid_reads),
                3,
            )
            Test.@test restored_contract.contents == expected_contract.contents

            run(`touch -t 203001010101 $(valid_reads)`)
            touched_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(nothing, valid_reads),
                3,
            )
            Test.@test touched_contract.contents == expected_contract.contents
            Test.@test touched_contract.signature == expected_contract.signature
            Test.@test touched_contract.input_snapshot !=
                       expected_contract.input_snapshot
            run(`touch -r $(metadata_reference) $(valid_reads)`)

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
                        return _test_metamdbg_toolchain()
                    end,
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"existing output contract does not match",
            )
            Test.@test local_provisioning_calls[] == 2
            Test.@test command_count[] == 2
        end

        Test.@testset "local output root identity rejects symlink swaps" begin
            outdir = joinpath(temporary_root, "local-output-root-swap")
            outputs = Mycelia._metamdbg_output_paths(outdir, 21)
            output_lock_path = Mycelia._metamdbg_output_lock_path(outdir)
            attack_target =
                joinpath(temporary_root, "local-output-root-swap-target")
            mkpath(attack_target)
            command_count = Ref(0)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = function (_command::Cmd)
                        command_count[] += 1
                        rm(outdir; recursive = true)
                        symlink(attack_target, outdir)
                        write(
                            joinpath(attack_target, "contigs.fasta"),
                            ">escaped-contig\nACGT\n",
                        )
                        return nothing
                    end,
                ),
                ErrorException,
                r"output root.*changed",
            )
            Test.@test command_count[] == 1
            Test.@test islink(outdir)
            Test.@test isfile(joinpath(attack_target, "contigs.fasta"))
            Test.@test !ispath(joinpath(
                attack_target,
                basename(outputs.contract_marker),
            ))
            Test.@test !ispath(joinpath(
                attack_target,
                basename(outputs.completion_marker),
            ))
            Test.@test !ispath(output_lock_path)
        end

        Test.@testset "local input mutation never stamps stale contract" begin
            mutation_reads = joinpath(temporary_root, "local-mutation.fastq")
            write(mutation_reads, "@local-mutation\nACGT\n+\nIIII\n")
            outdir = joinpath(temporary_root, "local-mutation-output")
            command_count = Ref(0)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = mutation_reads,
                    outdir,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = function (_command::Cmd)
                        command_count[] += 1
                        if command_count[] == 1
                            write(
                                joinpath(outdir, "contigs.fasta"),
                                ">contig-1\nACGT\n",
                            )
                        elseif command_count[] == 2
                            _write_test_metamdbg_gfa!(joinpath(
                                outdir,
                                "assemblyGraph_k21_4bps.gfa",
                            ))
                            write(
                                mutation_reads,
                                "@local-mutation\nTGCA\n+\nIIII\n",
                            )
                        else
                            error("unexpected metaMDBG command")
                        end
                        return nothing
                    end,
                ),
                ErrorException,
                r"input content contract changed",
            )
            outputs = Mycelia._metamdbg_output_paths(outdir, 21)
            Test.@test command_count[] == 2
            Test.@test !ispath(outputs.contract_marker)
            Test.@test !ispath(outputs.completion_marker)
            Test.@test !ispath(Mycelia._metamdbg_output_lock_path(outdir))

            restored_reads =
                joinpath(temporary_root, "local-mutate-restore.fastq")
            original_reads = "@local-mutate-restore\nACGT\n+\nIIII\n"
            write(restored_reads, original_reads)
            restored_metadata =
                joinpath(temporary_root, "local-mutate-restore-metadata.fastq")
            cp(restored_reads, restored_metadata)
            run(`touch -r $(restored_reads) $(restored_metadata)`)
            restored_outdir =
                joinpath(temporary_root, "local-mutate-restore-output")
            restored_commands = Ref(0)
            restored_result = Mycelia._run_metamdbg(;
                hifi_reads = restored_reads,
                outdir = restored_outdir,
                dependency_checker = _test_metamdbg_toolchain,
                local_runner = function (command::Cmd)
                    restored_commands[] += 1
                    if restored_commands[] == 1
                        input_flag_index = only(findall(
                            ==("--in-hifi"),
                            command.exec,
                        ))
                        staged_input = command.exec[input_flag_index + 1]
                        Test.@test staged_input != restored_reads
                        Test.@test read(staged_input, String) == original_reads
                        write(
                            restored_reads,
                            "@local-mutate-restore\nTGCA\n+\nIIII\n",
                        )
                        Test.@test read(staged_input, String) == original_reads
                        write(restored_reads, original_reads)
                        run(`touch -r $(restored_metadata) $(restored_reads)`)
                        write(
                            joinpath(restored_outdir, "contigs.fasta"),
                            ">contig-1\nACGT\n",
                        )
                    else
                        _write_test_metamdbg_gfa!(joinpath(
                            restored_outdir,
                            "assemblyGraph_k21_4bps.gfa",
                        ))
                    end
                    return nothing
                end,
            )
            Test.@test restored_result.status == :complete
            Test.@test restored_commands[] == 2
            Test.@test read(restored_reads, String) == original_reads

            drift_reads = joinpath(temporary_root, "toolchain-drift.fastq")
            write(drift_reads, "@toolchain-drift\nACGT\n+\nIIII\n")
            drift_outdir = joinpath(temporary_root, "toolchain-drift-output")
            dependency_calls = Ref(0)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = drift_reads,
                    outdir = drift_outdir,
                    dependency_checker = () -> begin
                        dependency_calls[] += 1
                        toolchain = _test_metamdbg_toolchain()
                        dependency_calls[] == 1 && return toolchain
                        changed_inventory = copy(toolchain.package_inventory)
                        changed_inventory[1] = merge(
                            changed_inventory[1],
                            (; build = "changed-after-tools"),
                        )
                        return Mycelia._require_metamdbg_package_version(
                            changed_inventory,
                        )
                    end,
                    local_runner = function (_command::Cmd)
                        if !ispath(joinpath(drift_outdir, "contigs.fasta"))
                            write(
                                joinpath(drift_outdir, "contigs.fasta"),
                                ">contig-1\nACGT\n",
                            )
                        else
                            _write_test_metamdbg_gfa!(joinpath(
                                drift_outdir,
                                "assemblyGraph_k21_4bps.gfa",
                            ))
                        end
                        return nothing
                    end,
                ),
                ErrorException,
                r"resolved package inventory changed",
            )
            Test.@test dependency_calls[] == 2
            drift_outputs =
                Mycelia._metamdbg_output_paths(drift_outdir, 21)
            Test.@test !ispath(drift_outputs.contract_marker)
            Test.@test !ispath(drift_outputs.completion_marker)
        end

        Test.@testset "input hashing is bounded at lifecycle boundaries" begin
            local_reads = joinpath(temporary_root, "bounded-local.fastq")
            write(local_reads, "@bounded-local\nACGT\n+\nIIII\n")
            local_outdir = joinpath(temporary_root, "bounded-local-output")
            local_digest_calls = String[]
            local_commands = Ref(0)
            local_result = Mycelia._run_metamdbg(;
                hifi_reads = local_reads,
                outdir = local_outdir,
                dependency_checker = _test_metamdbg_toolchain,
                input_digest_function = function (path::AbstractString)
                    push!(local_digest_calls, String(path))
                    return Mycelia._metamdbg_sha256(path)
                end,
                local_runner = function (_command::Cmd)
                    local_commands[] += 1
                    if local_commands[] == 1
                        write(
                            joinpath(local_outdir, "contigs.fasta"),
                            ">contig-1\nACGT\n",
                        )
                    else
                        _write_test_metamdbg_gfa!(joinpath(
                            local_outdir,
                            "assemblyGraph_k21_4bps.gfa",
                        ))
                    end
                    return nothing
                end,
            )
            Test.@test local_result.status == :complete
            Test.@test local_commands[] == 2
            Test.@test local_digest_calls == String[local_reads, local_reads]

            planned_reads = joinpath(temporary_root, "bounded-planned.fastq")
            write(planned_reads, "@bounded-planned\nACGT\n+\nIIII\n")
            planned_digest_calls = String[]
            planned_result = Mycelia._run_metamdbg(;
                hifi_reads = planned_reads,
                outdir = joinpath(temporary_root, "bounded-planned-output"),
                executor = Mycelia.CollectExecutor(),
                dependency_checker = () -> error(
                    "planned execution must not provision metaMDBG",
                ),
                local_runner = forbidden_runner,
                input_digest_function = function (path::AbstractString)
                    push!(planned_digest_calls, String(path))
                    return Mycelia._metamdbg_sha256(path)
                end,
            )
            Test.@test planned_result.status == :planned
            Test.@test planned_digest_calls == String[planned_reads]

            restored_reads = joinpath(temporary_root, "restored-mtime.fastq")
            restored_reference =
                joinpath(temporary_root, "restored-mtime-reference.fastq")
            write(restored_reads, "@restored-mtime\nACGT\n+\nIIII\n")
            cp(restored_reads, restored_reference)
            run(`touch -r $(restored_reads) $(restored_reference)`)
            restored_outdir =
                joinpath(temporary_root, "restored-mtime-output")
            restored_digest_calls = String[]
            restored_commands = Ref(0)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = restored_reads,
                    outdir = restored_outdir,
                    dependency_checker = _test_metamdbg_toolchain,
                    input_digest_function = function (path::AbstractString)
                        push!(restored_digest_calls, String(path))
                        return Mycelia._metamdbg_sha256(path)
                    end,
                    local_runner = function (_command::Cmd)
                        restored_commands[] += 1
                        if restored_commands[] == 1
                            write(
                                joinpath(restored_outdir, "contigs.fasta"),
                                ">contig-1\nACGT\n",
                            )
                        else
                            _write_test_metamdbg_gfa!(joinpath(
                                restored_outdir,
                                "assemblyGraph_k21_4bps.gfa",
                            ))
                            write(
                                restored_reads,
                                "@restored-mtime\nTGCA\n+\nIIII\n",
                            )
                            run(`touch -r $(restored_reference) $(restored_reads)`)
                        end
                        return nothing
                    end,
                ),
                ErrorException,
                r"input content contract changed",
            )
            Test.@test restored_commands[] == 2
            Test.@test restored_digest_calls ==
                       String[restored_reads, restored_reads]
            Test.@test !ispath(joinpath(
                restored_outdir,
                "mycelia_metamdbg_contract.json",
            ))
        end

        Test.@testset "exact graph-k discovery and ambiguity failure" begin
            partial_outdir =
                joinpath(temporary_root, "contracted-partial-contigs")
            mkpath(partial_outdir)
            _write_test_metamdbg_gzip_fasta!(joinpath(
                partial_outdir,
                "contigs.fasta.gz",
            ))
            _write_test_metamdbg_contract!(partial_outdir, valid_reads)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = partial_outdir,
                    dependency_checker = () -> error(
                        "partial output must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"partial contracted contigs without a realized-stage completion manifest",
            )

            legacy_complete_outdir =
                joinpath(temporary_root, "legacy-unbound-complete")
            mkpath(legacy_complete_outdir)
            _write_test_metamdbg_gzip_fasta!(joinpath(
                legacy_complete_outdir,
                "contigs.fasta.gz",
            ))
            _write_test_metamdbg_gfa!(joinpath(
                legacy_complete_outdir,
                "assemblyGraph_k21_4bps.gfa",
            ))
            _write_test_metamdbg_contract!(
                legacy_complete_outdir,
                valid_reads,
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = legacy_complete_outdir,
                    graph_k = 21,
                    dependency_checker = () -> error(
                        "unbound complete reuse must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"missing its regular, nonempty completion manifest",
            )

            outdir = joinpath(temporary_root, "exact-k")
            mkpath(outdir)
            _write_test_metamdbg_gzip_fasta!(
                joinpath(outdir, "contigs.fasta.gz"),
            )
            _write_test_metamdbg_gfa!(
                joinpath(outdir, "assemblyGraph_k20_4bps.gfa"),
            )
            _write_test_metamdbg_contract!(outdir, valid_reads)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> error(
                        "different graph_k must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"exactly one graph_k lifecycle per output root",
            )

            ambiguous_outdir = joinpath(temporary_root, "ambiguous")
            mkpath(ambiguous_outdir)
            _write_test_metamdbg_gzip_fasta!(
                joinpath(ambiguous_outdir, "contigs.fasta.gz"),
            )
            _write_test_metamdbg_gfa!(
                joinpath(ambiguous_outdir, "assemblyGraph_k21_4bps.gfa"),
            )
            _write_test_metamdbg_gfa!(
                joinpath(ambiguous_outdir, "assemblyGraph_k21_5bps.gfa"),
                sequence = "TGCA",
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
                r"multiple graph artifacts for k=21",
            )

            mixed_k_outdir = joinpath(temporary_root, "mixed-graph-k")
            mkpath(mixed_k_outdir)
            _write_test_metamdbg_gzip_fasta!(joinpath(
                mixed_k_outdir,
                "contigs.fasta.gz",
            ))
            _write_test_metamdbg_gfa!(joinpath(
                mixed_k_outdir,
                "assemblyGraph_k21_4bps.gfa",
            ))
            _write_test_metamdbg_gfa!(
                joinpath(mixed_k_outdir, "assemblyGraph_k31_4bps.gfa");
                sequence = "TGCA",
            )
            _write_test_metamdbg_contract!(mixed_k_outdir, valid_reads)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = mixed_k_outdir,
                    graph_k = 21,
                    dependency_checker = () -> error(
                        "mixed graph_k output must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"exactly one total graph artifact per output root",
            )
        end

        Test.@testset "spec-addressed exact environment and install lock" begin
            _, bundled_environment = Mycelia._metamdbg_paths()
            Test.@test Mycelia._metamdbg_sha256(bundled_environment) ==
                       Mycelia.METAMDBG_ENVIRONMENT_SPEC_SHA256
            Test.@test occursin(
                "metamdbg==1.4",
                read(bundled_environment, String),
            )
            Test.@test occursin(
                "python=3.11",
                read(bundled_environment, String),
            )
            Test.@test Mycelia._METAMDBG_INSTALL_LOCK_STALE_SECONDS > 0
            Test.@test Mycelia._METAMDBG_INSTALL_LOCK_REFRESH_SECONDS > 0
            Test.@test Mycelia._METAMDBG_INSTALL_LOCK_REFRESH_SECONDS <
                       Mycelia._METAMDBG_INSTALL_LOCK_STALE_SECONDS

            inventory = Mycelia._metamdbg_environment_packages(;
                conda_runner = "fixture-conda",
                command_reader = function (command::Cmd)
                    rendered = string(command)
                    Test.@test occursin("fixture-conda list -n", rendered)
                    Test.@test occursin(Mycelia.METAMDBG_ENV_NAME, rendered)
                    return "[" *
                           "{\"name\":\"metamdbg\",\"version\":\"1.4\"," *
                           "\"build_string\":\"h43eeafb_2\"," *
                           "\"channel\":\"bioconda\"}," *
                           "{\"name\":\"libzlib\",\"version\":\"1.3.1\"," *
                           "\"build_string\":\"h8359307_2\"," *
                           "\"channel\":\"conda-forge\"}]"
                end,
            )
            Test.@test inventory ==
                       _test_metamdbg_toolchain().package_inventory
            realized_toolchain =
                Mycelia._require_metamdbg_package_version(inventory)
            Test.@test realized_toolchain == _test_metamdbg_toolchain()
            Test.@test realized_toolchain.package_count == 2
            Test.@test occursin(
                r"^[0-9a-f]{64}$",
                realized_toolchain.package_inventory_sha256,
            )
            changed_build_toolchain =
                Mycelia._require_metamdbg_package_version(NamedTuple[
                    (;
                        name = "libzlib",
                        version = "1.3.1",
                        build = "h8359307_3",
                        channel = "conda-forge",
                    ),
                    (;
                        name = "metamdbg",
                        version = "1.4",
                        build = "h43eeafb_2",
                        channel = "bioconda",
                    ),
                ])
            Test.@test changed_build_toolchain.metamdbg_version ==
                       realized_toolchain.metamdbg_version
            Test.@test changed_build_toolchain.package_inventory_sha256 !=
                       realized_toolchain.package_inventory_sha256

            escaped_inventory = Any[
                Dict(
                    "channel" => "bio\"conda",
                    "build_string" => "h43eeafb_2",
                    "version" => "1.4",
                    "name" => "metamdbg",
                ),
                Dict(
                    "version" => "1.3.1",
                    "name" => "libzlib",
                    "channel" => "conda-forge",
                    "build" => "h8359307_2",
                ),
            ]
            inventory_json = joinpath(
                temporary_root,
                "escaped-reordered-inventory.json",
            )
            inventory_tsv = joinpath(
                temporary_root,
                "escaped-reordered-inventory.tsv",
            )
            serialized_inventory = JSON.json(reverse(escaped_inventory))
            write(
                inventory_json,
                "  \n" * replace(serialized_inventory, "},{" => "},\n {") *
                "\n",
            )
            run(`python3 -c $(Mycelia._metamdbg_runtime_inventory_canonicalizer_python()) $(inventory_json) $(inventory_tsv)`)
            Test.@test read(inventory_tsv, String) ==
                       Mycelia._metamdbg_package_inventory_contents(
                escaped_inventory,
            )
            escaped_changed = deepcopy(escaped_inventory)
            escaped_changed[1]["build_string"] = "h43eeafb_3"
            changed_inventory_json = joinpath(
                temporary_root,
                "changed-build-inventory.json",
            )
            changed_inventory_tsv = joinpath(
                temporary_root,
                "changed-build-inventory.tsv",
            )
            write(changed_inventory_json, JSON.json(escaped_changed))
            run(`python3 -c $(Mycelia._metamdbg_runtime_inventory_canonicalizer_python()) $(changed_inventory_json) $(changed_inventory_tsv)`)
            Test.@test Mycelia._metamdbg_sha256(changed_inventory_tsv) !=
                       Mycelia._metamdbg_sha256(inventory_tsv)
            _test_metamdbg_error(
                () -> Mycelia._require_metamdbg_package_version(Any[
                    Dict("name" => "metamdbg", "version" => "1.4"),
                ]),
                ErrorException,
                r"missing or empty build field",
            )
            _test_metamdbg_error(
                () -> Mycelia._require_metamdbg_package_version(
                    NamedTuple[(;
                        name = "metamdbg",
                        version = "1.3",
                        build = "h43eeafb_2",
                        channel = "bioconda",
                    )],
                ),
                ErrorException,
                r"exactly 1.4",
            )

            install_dir = joinpath(temporary_root, "metamdbg-install")
            environment_path = joinpath(install_dir, "environment.yml")
            mkpath(install_dir)
            cp(bundled_environment, environment_path)
            created = Ref(false)
            lock_calls = Ref(0)
            installed_toolchain = Mycelia._ensure_metamdbg_installed(;
                paths = (install_dir, environment_path),
                environment_checker = function (
                        _runner::AbstractString,
                )
                    return false
                end,
                environment_creator = function (
                        specification::AbstractString,
                        environment_name::AbstractString,
                        _runner::AbstractString;
                        force::Bool,
                )
                    Test.@test specification == environment_path
                    Test.@test environment_name == Mycelia.METAMDBG_ENV_NAME
                    Test.@test !force
                    created[] = true
                    return String(environment_name)
                end,
                package_inspector = function (
                        _runner::AbstractString,
                )
                    return _test_metamdbg_toolchain().package_inventory
                end,
                lock_path = joinpath(temporary_root, "injected-install.pid"),
                lock_runner = function (
                        action::Function,
                        lock_path::AbstractString,
                )
                    lock_calls[] += 1
                    Test.@test endswith(lock_path, "injected-install.pid")
                    return action()
                end,
            )
            Test.@test created[]
            Test.@test lock_calls[] == 1
            Test.@test installed_toolchain == _test_metamdbg_toolchain()

            runner_root_one = joinpath(temporary_root, "conda-root-one")
            runner_root_two = joinpath(temporary_root, "conda-root-two")
            runner_one = joinpath(runner_root_one, "bin", "conda")
            runner_two = joinpath(runner_root_two, "bin", "conda")
            for runner in (runner_one, runner_two)
                mkpath(dirname(runner))
                write(runner, "#!/bin/sh\nexit 0\n")
                chmod(runner, 0o700)
            end
            canonical_runner_one = realpath(runner_one)
            canonical_runner_two = realpath(runner_two)
            runner_one_lock = Mycelia._metamdbg_install_lock_path(runner_one)
            runner_two_lock = Mycelia._metamdbg_install_lock_path(runner_two)
            runner_one_alias = joinpath(temporary_root, "conda-one-alias")
            symlink(runner_one, runner_one_alias)
            Test.@test runner_one_lock == joinpath(
                realpath(runner_root_one),
                ".mycelia-locks",
                "$(Mycelia.METAMDBG_ENV_NAME).pid",
            )
            Test.@test Mycelia._metamdbg_install_lock_path(runner_one_alias) ==
                       runner_one_lock
            Test.@test runner_two_lock != runner_one_lock
            mkpath(joinpath(
                runner_root_one,
                "envs",
                Mycelia.METAMDBG_ENV_NAME,
                "conda-meta",
            ))
            Test.@test Mycelia._metamdbg_environment_is_installed(runner_one)
            Test.@test !Mycelia._metamdbg_environment_is_installed(runner_two)
            original_depot_path = copy(Base.DEPOT_PATH)
            try
                Base.DEPOT_PATH[1] = joinpath(temporary_root, "other-depot")
                Test.@test Mycelia._metamdbg_install_lock_path(runner_one) ==
                           runner_one_lock
            finally
                empty!(Base.DEPOT_PATH)
                append!(Base.DEPOT_PATH, original_depot_path)
            end

            observed_inventory_runner = Ref("")
            runner_inventory = Mycelia._metamdbg_environment_packages(;
                conda_runner = runner_one,
                command_reader = function (command::Cmd)
                    observed_inventory_runner[] = first(command.exec)
                    return "[" *
                           "{\"name\":\"metamdbg\",\"version\":\"1.4\"," *
                           "\"build_string\":\"h43eeafb_2\"," *
                           "\"channel\":\"bioconda\"}," *
                           "{\"name\":\"libzlib\",\"version\":\"1.3.1\"," *
                           "\"build_string\":\"h8359307_2\"," *
                           "\"channel\":\"conda-forge\"}]"
                end,
            )
            Test.@test observed_inventory_runner[] == canonical_runner_one
            Test.@test runner_inventory ==
                       _test_metamdbg_toolchain().package_inventory

            ensured_runners = String[]
            created_with_runner = Ref(false)
            observed_default_lock = Ref("")
            withenv("MYCELIA_CONDA_RUNNER" => runner_one) do
                default_runner_toolchain = Mycelia._ensure_metamdbg_installed(;
                    paths = (install_dir, environment_path),
                    environment_checker = function (
                            runner::AbstractString,
                    )
                        push!(ensured_runners, String(runner))
                        return false
                    end,
                    environment_creator = function (
                            specification::AbstractString,
                            environment_name::AbstractString,
                            runner::AbstractString;
                            force::Bool,
                    )
                        Test.@test specification == environment_path
                        Test.@test environment_name ==
                                   Mycelia.METAMDBG_ENV_NAME
                        Test.@test runner == canonical_runner_one
                        Test.@test !force
                        created_with_runner[] = true
                        return String(environment_name)
                    end,
                    package_inspector = function (
                            runner::AbstractString,
                    )
                        push!(ensured_runners, String(runner))
                        return _test_metamdbg_toolchain().package_inventory
                    end,
                    lock_runner = function (
                            action::Function,
                            lock_path::AbstractString,
                    )
                        observed_default_lock[] = String(lock_path)
                        return action()
                    end,
                )
                Test.@test default_runner_toolchain ==
                           _test_metamdbg_toolchain()
            end
            Test.@test created_with_runner[]
            Test.@test ensured_runners ==
                       String[canonical_runner_one, canonical_runner_one]
            Test.@test observed_default_lock[] == runner_one_lock

            explicit_runner_outdir =
                joinpath(temporary_root, "explicit-runner-local")
            explicit_runner_commands = Cmd[]
            explicit_dependency_runners = String[]
            explicit_runner_result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir = explicit_runner_outdir,
                conda_runner = runner_two,
                dependency_checker = function (
                        runner::AbstractString,
                )
                    push!(explicit_dependency_runners, String(runner))
                    return _test_metamdbg_toolchain()
                end,
                local_runner = function (command::Cmd)
                    push!(explicit_runner_commands, command)
                    Test.@test first(command.exec) == canonical_runner_two
                    if length(explicit_runner_commands) == 1
                        write(
                            joinpath(explicit_runner_outdir, "contigs.fasta"),
                            ">contig-1\nACGT\n",
                        )
                    else
                        _write_test_metamdbg_gfa!(joinpath(
                            explicit_runner_outdir,
                            "assemblyGraph_k21_4bps.gfa",
                        ))
                    end
                    return nothing
                end,
            )
            Test.@test explicit_runner_result.status == :complete
            Test.@test length(explicit_runner_commands) == 2
            Test.@test explicit_dependency_runners ==
                       String[canonical_runner_two, canonical_runner_two]

            default_runner_executor = Mycelia.CollectExecutor()
            withenv("MYCELIA_CONDA_RUNNER" => runner_one) do
                default_runner_result = Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    outdir = joinpath(
                        temporary_root,
                        "default-runner-nonlocal",
                    ),
                    executor = default_runner_executor,
                    dependency_checker = () -> error(
                        "planned runner test must not provision metaMDBG",
                    ),
                    local_runner = forbidden_runner,
                )
                Test.@test default_runner_result.status == :planned
            end
            default_runner_script = only(default_runner_executor.jobs).cmd
            Test.@test occursin(
                "conda_runner='$(canonical_runner_one)'",
                default_runner_script,
            )

            real_lock_path = joinpath(temporary_root, "real-install.pid")
            lock_result = Mycelia._with_metamdbg_install_lock(
                () -> begin
                    Test.@test ispath(real_lock_path)
                    return :locked
                end,
                real_lock_path,
            )
            Test.@test lock_result == :locked
            Test.@test !ispath(real_lock_path)
        end

        Test.@testset "semantic artifact failures never stamp contracts" begin
            empty_identifier_fasta =
                joinpath(temporary_root, "empty-identifier.fasta")
            write(empty_identifier_fasta, ">\nACGT\n")
            _test_metamdbg_error(
                () -> Mycelia._require_valid_metamdbg_fasta(
                    empty_identifier_fasta,
                    "empty-identifier fixture FASTA",
                ),
                ErrorException,
                r"empty FASTA identifier|not valid FASTA",
            )
            duplicate_identifier_fasta =
                joinpath(temporary_root, "duplicate-identifier.fasta")
            write(
                duplicate_identifier_fasta,
                ">duplicate first\nACGT\n>duplicate second\nTGCA\n",
            )
            _test_metamdbg_error(
                () -> Mycelia._require_valid_metamdbg_fasta(
                    duplicate_identifier_fasta,
                    "duplicate-identifier fixture FASTA",
                ),
                ErrorException,
                r"duplicate FASTA identifier",
            )
            whitespace_identifier_fasta =
                joinpath(temporary_root, "whitespace-identifier.fasta")
            write(whitespace_identifier_fasta, "> leading-space\nACGT\n")
            _test_metamdbg_error(
                () -> Mycelia._require_valid_metamdbg_fasta(
                    whitespace_identifier_fasta,
                    "whitespace-identifier fixture FASTA",
                ),
                ErrorException,
                r"empty FASTA identifier|not valid FASTA",
            )

            invalid_dna_fasta = joinpath(temporary_root, "invalid-dna.fasta")
            write(invalid_dna_fasta, ">invalid-dna\nACGTZ\n")
            _test_metamdbg_error(
                () -> Mycelia._require_valid_metamdbg_fasta(
                    invalid_dna_fasta,
                    "invalid fixture FASTA",
                ),
                ErrorException,
                r"invalid DNA",
            )
            invalid_dna_gfa = joinpath(temporary_root, "invalid-dna.gfa")
            _write_test_metamdbg_gfa!(
                invalid_dna_gfa;
                sequence = "ACGTZ",
            )
            _test_metamdbg_error(
                () -> Mycelia._require_valid_metamdbg_gfa(
                    invalid_dna_gfa,
                    "invalid fixture GFA",
                ),
                ErrorException,
                r"invalid DNA",
            )
            gapped_dna_fasta = joinpath(temporary_root, "gapped-dna.fasta")
            write(gapped_dna_fasta, ">gapped-dna\nAC-GT\n")
            _test_metamdbg_error(
                () -> Mycelia._require_valid_metamdbg_fasta(
                    gapped_dna_fasta,
                    "gapped fixture FASTA",
                ),
                ErrorException,
                r"invalid DNA",
            )
            gapped_dna_gfa = joinpath(temporary_root, "gapped-dna.gfa")
            _write_test_metamdbg_gfa!(
                gapped_dna_gfa;
                sequence = "AC-GT",
            )
            _test_metamdbg_error(
                () -> Mycelia._require_valid_metamdbg_gfa(
                    gapped_dna_gfa,
                    "gapped fixture GFA",
                ),
                ErrorException,
                r"invalid DNA",
            )

            structured_gfa = joinpath(temporary_root, "structured.gfa")
            write(
                structured_gfa,
                "H\tVN:Z:1.0\n" *
                "S\tcontig-1\tACGT\n" *
                "S\tcontig-2\tTGCA\n" *
                "S\tcontig,with,comma\tGATTACA\tLN:i:7\n" *
                "L\tcontig-1\t+\tcontig-2\t-\t4M\n" *
                "P\tpath-1\tcontig-1+,contig-2-\t4M\n" *
                "P\tcomma-path\tcontig,with,comma+\t*\n",
            )
            Test.@test Mycelia._require_valid_metamdbg_gfa(
                structured_gfa,
                "structured fixture GFA",
            ) == structured_gfa

            typed_tags_gfa = joinpath(temporary_root, "typed-tags.gfa")
            write(
                typed_tags_gfa,
                "H\tVN:Z:1.2\n" *
                "S\tcontig-1\tACGT\tZA:Z:text\tJA:J:{\"ok\":[1,true]}" *
                "\tHA:H:0AFF\tBc:B:c,-128,127\tBC:B:C,0,255" *
                "\tBs:B:s,-32768,32767\tBS:B:S,0,65535" *
                "\tBi:B:i,-2147483648,2147483647" *
                "\tBI:B:I,0,4294967295\tBf:B:f,-1.5,3e2\n",
            )
            Test.@test Mycelia._require_valid_metamdbg_gfa(
                typed_tags_gfa,
                "typed tags fixture GFA",
            ) == typed_tags_gfa

            many_link_gfa = joinpath(temporary_root, "many-forward-links.gfa")
            open(many_link_gfa, "w") do output
                write(output, "H\tVN:Z:1.0\n")
                for _ in 1:1_000
                    write(output, "L\tcontig-1\t+\tcontig-2\t-\t4M\n")
                end
                write(output, "S\tcontig-1\tACGT\nS\tcontig-2\tTGCA\n")
            end
            Test.@test Mycelia._require_valid_metamdbg_gfa(
                many_link_gfa,
                "many-link fixture GFA",
            ) == many_link_gfa

            small_path_field = join(
                ("segment-$(index)+" for index in 1:2_000),
                ',',
            )
            large_path_field = join(
                ("segment-$(index)+" for index in 1:4_000),
                ',',
            )
            Mycelia._metamdbg_gfa_path_step_identifiers(
                small_path_field,
                1,
                "allocation warmup",
                "fixture.gfa",
            )
            GC.gc()
            small_path_allocations = @allocated begin
                Mycelia._metamdbg_gfa_path_step_identifiers(
                    small_path_field,
                    1,
                    "small allocation fixture",
                    "fixture.gfa",
                )
            end
            GC.gc()
            large_path_allocations = @allocated begin
                Mycelia._metamdbg_gfa_path_step_identifiers(
                    large_path_field,
                    1,
                    "large allocation fixture",
                    "fixture.gfa",
                )
            end
            Test.@test large_path_allocations <=
                       3 * small_path_allocations + 1_000_000
            long_path_gfa = joinpath(temporary_root, "long-path.gfa")
            long_steps = join(fill("contig-1+", 10_000), ',')
            write(
                long_path_gfa,
                "H\tVN:Z:1.0\nS\tcontig-1\tACGT\n" *
                "P\tlong-path\t$(long_steps)\t*\n",
            )
            Test.@test Mycelia._require_valid_metamdbg_gfa(
                long_path_gfa,
                "long path fixture GFA",
            ) == long_path_gfa

            invalid_gfa_semantics = (
                leading_star = (
                    "S\t*segment\tACGT\n",
                    r"invalid GFA segment identifier",
                ),
                leading_equals = (
                    "S\t=segment\tACGT\n",
                    r"invalid GFA segment identifier",
                ),
                shared_namespace = (
                    "S\tshared-name\tACGT\n" *
                    "P\tshared-name\tshared-name+\t*\n",
                    r"duplicate GFA segment/path name",
                ),
                duplicate_tag = (
                    "S\tcontig-1\tACGT\tLN:i:4\tLN:i:5\n",
                    r"duplicate GFA segment optional tag name",
                ),
                invalid_integer_tag = (
                    "S\tcontig-1\tACGT\tLN:i:not-an-integer\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                empty_z_tag = (
                    "S\tcontig-1\tACGT\tZZ:Z:\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                invalid_json_tag = (
                    "S\tcontig-1\tACGT\tJJ:J:{not-json}\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                nonstandard_json_constant = (
                    "S\tcontig-1\tACGT\tJJ:J:NaN\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                odd_hex_tag = (
                    "S\tcontig-1\tACGT\tHH:H:ABC\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                b_signed_range = (
                    "S\tcontig-1\tACGT\tBB:B:c,-129\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                b_unsigned_byte_range = (
                    "S\tcontig-1\tACGT\tBB:B:C,-1\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                b_signed_short_range = (
                    "S\tcontig-1\tACGT\tBB:B:s,32768\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                b_unsigned_short_range = (
                    "S\tcontig-1\tACGT\tBB:B:S,65536\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                b_signed_int_range = (
                    "S\tcontig-1\tACGT\tBB:B:i,2147483648\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                b_unsigned_range = (
                    "S\tcontig-1\tACGT\tBB:B:I,4294967296\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                b_float_range = (
                    "S\tcontig-1\tACGT\tBB:B:f,1e400\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                b_unknown_subtype = (
                    "S\tcontig-1\tACGT\tBB:B:d,1\n",
                    r"invalid or unsupported GFA segment optional tag value",
                ),
                non_gfa1_version = (
                    "H\tVN:Z:2.0\nS\tcontig-1\tACGT\n",
                    r"non-GFA1 VN header",
                ),
                leading_whitespace_pseudo_comment = (
                    " #not-a-comment\nS\tcontig-1\tACGT\n",
                    r"unknown GFA record type",
                ),
            )
            for (failure_type, fixture) in pairs(invalid_gfa_semantics)
                invalid_gfa = joinpath(
                    temporary_root,
                    "gfa-semantics-$(failure_type).gfa",
                )
                write(invalid_gfa, fixture[1])
                _test_metamdbg_error(
                    () -> Mycelia._require_valid_metamdbg_gfa(
                        invalid_gfa,
                        "$(failure_type) fixture GFA",
                    ),
                    ErrorException,
                    fixture[2],
                )
            end

            invalid_structured_gfas = (
                leading_garbage = (
                    "garbage\nS\tcontig-1\tACGT\n",
                    r"unknown GFA record type",
                ),
                malformed_link = (
                    "S\tcontig-1\tACGT\n" *
                    "S\tcontig-2\tTGCA\n" *
                    "L\tcontig-1\t?\tcontig-2\t+\t4M\n",
                    r"invalid GFA link source orientation",
                ),
                dangling_link = (
                    "S\tcontig-1\tACGT\n" *
                    "L\tcontig-1\t+\tmissing\t+\t4M\n",
                    r"dangling GFA link segment reference",
                ),
            )
            for (failure_type, (contents, expected_message)) in
                pairs(invalid_structured_gfas)
                invalid_gfa = joinpath(
                    temporary_root,
                    "$(failure_type).gfa",
                )
                write(invalid_gfa, contents)
                _test_metamdbg_error(
                    () -> Mycelia._require_valid_metamdbg_gfa(
                        invalid_gfa,
                        "$(failure_type) fixture GFA",
                    ),
                    ErrorException,
                    expected_message,
                )
            end

            malformed_fasta_outdir = joinpath(temporary_root, "malformed-fasta")
            fasta_calls = Ref(0)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = malformed_fasta_outdir,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = function (_command::Cmd)
                        fasta_calls[] += 1
                        write(
                            joinpath(malformed_fasta_outdir, "contigs.fasta"),
                            ">empty-contig\n",
                        )
                        return nothing
                    end,
                ),
                ErrorException,
                r"empty FASTA sequence|no FASTA records|not valid FASTA",
            )
            Test.@test fasta_calls[] == 1
            Test.@test !ispath(joinpath(
                malformed_fasta_outdir,
                "mycelia_metamdbg_contract.json",
            ))
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = malformed_fasta_outdir,
                    dependency_checker = () -> error(
                        "uncontracted partial output must fail first",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"nonempty output root without its regular provenance contract marker",
            )

            malformed_graph_outdir = joinpath(temporary_root, "malformed-graph")
            graph_calls = Ref(0)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    outdir = malformed_graph_outdir,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = function (_command::Cmd)
                        graph_calls[] += 1
                        if graph_calls[] == 1
                            write(
                                joinpath(
                                    malformed_graph_outdir,
                                    "contigs.fasta",
                                ),
                                ">contig-1\nACGT\n",
                            )
                        else
                            write(
                                joinpath(
                                    malformed_graph_outdir,
                                    "assemblyGraph_k21_4bps.gfa",
                                ),
                                "H\tVN:Z:1.0\n",
                            )
                        end
                        return nothing
                    end,
                ),
                ErrorException,
                r"no sequence-bearing GFA segments",
            )
            Test.@test graph_calls[] == 2
            Test.@test !ispath(joinpath(
                malformed_graph_outdir,
                "mycelia_metamdbg_contract.json",
            ))

            interrupted_outdir = joinpath(temporary_root, "interrupted")
            interruption = try
                Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = interrupted_outdir,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = function (_command::Cmd)
                        write(
                            joinpath(interrupted_outdir, "contigs.fasta"),
                            ">partial\nACGT\n",
                        )
                        throw(InterruptException())
                    end,
                )
                nothing
            catch caught
                caught
            end
            Test.@test interruption isa InterruptException
            Test.@test !ispath(joinpath(
                interrupted_outdir,
                "mycelia_metamdbg_contract.json",
            ))
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(interrupted_outdir),
            )

            contracted_fasta_outdir =
                joinpath(temporary_root, "contracted-malformed-fasta")
            mkpath(contracted_fasta_outdir)
            _write_test_metamdbg_gzip_fasta!(
                joinpath(contracted_fasta_outdir, "contigs.fasta.gz");
                contents = ">empty-contig\n",
            )
            _write_test_metamdbg_gfa!(
                joinpath(
                    contracted_fasta_outdir,
                    "assemblyGraph_k21_4bps.gfa",
                ),
            )
            _write_test_metamdbg_contract!(
                contracted_fasta_outdir,
                valid_reads,
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = contracted_fasta_outdir,
                    dependency_checker = () -> error(
                        "malformed reuse must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"empty FASTA sequence|no FASTA records|not valid FASTA",
            )

            contracted_graph_outdir =
                joinpath(temporary_root, "contracted-malformed-graph")
            mkpath(contracted_graph_outdir)
            _write_test_metamdbg_gzip_fasta!(
                joinpath(contracted_graph_outdir, "contigs.fasta.gz"),
            )
            _write_test_metamdbg_gfa!(
                joinpath(
                    contracted_graph_outdir,
                    "assemblyGraph_k21_4bps.gfa",
                );
                sequence = "*",
            )
            _write_test_metamdbg_contract!(
                contracted_graph_outdir,
                valid_reads,
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = contracted_graph_outdir,
                    dependency_checker = () -> error(
                        "malformed reuse must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"no sequence for GFA segment",
            )
        end

        Test.@testset "canonical output ownership and lifecycle locks" begin
            physical_parent = joinpath(temporary_root, "physical-parent")
            alias_parent = joinpath(temporary_root, "alias-parent")
            mkpath(physical_parent)
            symlink(physical_parent, alias_parent)
            aliased_outdir = joinpath(alias_parent, "assembly")
            physical_outdir = joinpath(physical_parent, "assembly")
            Test.@test Mycelia._metamdbg_output_paths(
                aliased_outdir,
                21,
            ).outdir == physical_outdir
            Test.@test Mycelia._metamdbg_output_lock_path(aliased_outdir) ==
                       Mycelia._metamdbg_output_lock_path(physical_outdir)

            nested_outdir = joinpath(
                temporary_root,
                "new-parent",
                "nested",
                "assembly",
            )
            Test.@test isempty(
                Mycelia._metamdbg_submission_reservation_paths(nested_outdir),
            )
            nested_result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir = nested_outdir,
                executor = Mycelia.CollectExecutor(),
                dependency_checker = () -> error(
                    "nested planned execution must not provision metaMDBG",
                ),
                local_runner = forbidden_runner,
            )
            Test.@test nested_result.status == :planned
            Test.@test isdir(dirname(nested_outdir))

            remnant_outdir = joinpath(temporary_root, "reservation-remnants")
            remnant_prefix =
                Mycelia._metamdbg_submission_reservation_prefix(remnant_outdir)
            temporary_remnant =
                joinpath(temporary_root, remnant_prefix * "tmp.fixture")
            consumed_remnant =
                joinpath(temporary_root, remnant_prefix * "consumed.fixture")
            mkdir(temporary_remnant)
            mkdir(consumed_remnant)
            Test.@test isempty(
                Mycelia._metamdbg_submission_reservation_paths(remnant_outdir),
            )
            malformed_remnant =
                joinpath(temporary_root, remnant_prefix * "malformed")
            mkdir(malformed_remnant)
            _test_metamdbg_error(
                () -> Mycelia._metamdbg_submission_reservation_paths(
                    remnant_outdir,
                ),
                ErrorException,
                r"malformed submission reservation entry",
            )
            rm(malformed_remnant)
            rm(temporary_remnant)
            rm(consumed_remnant)

            missing_parent_root =
                joinpath(temporary_root, "reclaim-missing-parent")
            missing_parent_outdir = joinpath(missing_parent_root, "assembly")
            missing_parent_outputs =
                Mycelia._metamdbg_output_paths(missing_parent_outdir, 21)
            missing_parent_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            missing_parent_reservation =
                Mycelia._metamdbg_submission_reservation(
                    missing_parent_outputs,
                    missing_parent_contract,
                    21;
                    owner_token = "missing-parent-reclaim-fixture",
                )
            Mycelia._with_metamdbg_output_lock(missing_parent_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    missing_parent_reservation,
                    missing_parent_outdir,
                )
            end
            missing_parent_metadata = (;
                canonical_outdir = missing_parent_reservation.canonical_outdir,
                path = missing_parent_reservation.path,
                workflow_signature =
                    missing_parent_reservation.workflow_signature,
                input_contract_signature =
                    missing_parent_reservation.input_contract_signature,
                graph_k = missing_parent_reservation.graph_k,
                owner_token = missing_parent_reservation.owner_token,
                job_id = "899",
            )
            rm(missing_parent_root; recursive = true)
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    missing_parent_metadata;
                    owner_token = missing_parent_metadata.owner_token,
                    job_id = missing_parent_metadata.job_id,
                    confirm_cancelled = true,
                ),
                ErrorException,
                r"missing or was already consumed",
            )
            Test.@test !ispath(missing_parent_root)

            symlink_target = joinpath(temporary_root, "symlink-target")
            symlink_outdir = joinpath(temporary_root, "symlink-outdir")
            mkpath(symlink_target)
            symlink(symlink_target, symlink_outdir)
            _test_metamdbg_error(
                () -> Mycelia._metamdbg_output_paths(symlink_outdir, 21),
                ArgumentError,
                r"must not be a symbolic link",
            )

            dangling_ancestor = joinpath(temporary_root, "dangling-ancestor")
            symlink(joinpath(temporary_root, "missing-target"), dangling_ancestor)
            _test_metamdbg_error(
                () -> Mycelia._metamdbg_output_paths(
                    joinpath(dangling_ancestor, "assembly"),
                    21,
                ),
                ArgumentError,
                r"dangling symbolic-link ancestor",
            )

            locked_outdir = joinpath(temporary_root, "locked-output")
            locked_path = Mycelia._metamdbg_output_lock_path(locked_outdir)
            mkdir(locked_path)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = locked_outdir,
                    dependency_checker = () -> error(
                        "locked output must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"locked by another lifecycle",
            )
            Test.@test !ispath(locked_outdir)
            rm(locked_path)

            lifecycle_outdir = joinpath(temporary_root, "lifecycle-output")
            lifecycle_lock =
                Mycelia._metamdbg_output_lock_path(lifecycle_outdir)
            lock_result = Mycelia._with_metamdbg_output_lock(
                () -> begin
                    Test.@test isdir(lifecycle_lock)
                    return :locked
                end,
                lifecycle_outdir,
            )
            Test.@test lock_result == :locked
            Test.@test !ispath(lifecycle_lock)

            primary_outdir = joinpath(temporary_root, "primary-failure-output")
            primary_lock = Mycelia._metamdbg_output_lock_path(primary_outdir)
            primary_failure = ErrorException("primary workflow failure")
            cleanup_failure = ErrorException("synthetic lock cleanup failure")
            preserved_failure = try
                Mycelia._with_metamdbg_output_lock(
                    primary_outdir;
                    lock_remover = _path -> throw(cleanup_failure),
                ) do
                    Test.@test isdir(primary_lock)
                    throw(primary_failure)
                end
                nothing
            catch caught
                caught
            end
            Test.@test preserved_failure === primary_failure
            Test.@test isdir(primary_lock)
            rm(primary_lock)

            interrupted_outdir =
                joinpath(temporary_root, "cleanup-interrupted-output")
            interrupted_lock =
                Mycelia._metamdbg_output_lock_path(interrupted_outdir)
            primary_interrupt = InterruptException()
            preserved_interrupt = try
                Mycelia._with_metamdbg_output_lock(
                    interrupted_outdir;
                    lock_remover = _path -> throw(cleanup_failure),
                ) do
                    Test.@test isdir(interrupted_lock)
                    throw(primary_interrupt)
                end
                nothing
            catch caught
                caught
            end
            Test.@test preserved_interrupt === primary_interrupt
            Test.@test isdir(interrupted_lock)
            rm(interrupted_lock)
        end

        Test.@testset "queued nonlocal work blocks local execution" begin
            outdir = joinpath(temporary_root, "queued-blocks-local")
            outputs = Mycelia._metamdbg_output_paths(outdir, 21)
            input_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            reservation = Mycelia._metamdbg_submission_reservation(
                outputs,
                input_contract,
                21;
                owner_token = "queued-blocks-local-fixture",
            )
            Mycelia._with_metamdbg_output_lock(outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    reservation,
                    outdir,
                )
            end
            blocked_digest_calls = Ref(0)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    dependency_checker = () -> error(
                        "queued reservation must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                    input_digest_function = function (_path::AbstractString)
                        blocked_digest_calls[] += 1
                        error("reservation probe must precede input hashing")
                    end,
                ),
                ErrorException,
                r"active nonlocal submission reservation",
            )
            Test.@test blocked_digest_calls[] == 0
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    executor = Mycelia.CollectExecutor(),
                    dependency_checker = () -> error(
                        "queued reservation must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                    input_digest_function = function (_path::AbstractString)
                        blocked_digest_calls[] += 1
                        error("reservation probe must precede input hashing")
                    end,
                ),
                ErrorException,
                r"active nonlocal submission reservation",
            )
            Test.@test blocked_digest_calls[] == 0
            Test.@test isdir(reservation.path)
            Test.@test !ispath(outputs.contract_marker)
            Mycelia._with_metamdbg_output_lock(outdir) do
                Mycelia._remove_metamdbg_submission_reservation!(reservation)
            end
            Test.@test !ispath(reservation.path)
        end

        Test.@testset "shared hierarchical output-root reservations" begin
            same_root = joinpath(temporary_root, "shared-domain-same")
            Mycelia._with_metamdbg_output_domain_lock(same_root) do
                conflict = fetch(@async try
                    Mycelia._with_unicycler_output_lock(same_root) do _reserved
                        :unexpected
                    end
                catch caught
                    caught
                end)
                Test.@test conflict isa ArgumentError
                Test.@test occursin(
                    "already reserved",
                    sprint(showerror, conflict),
                )
            end

            parent_root = joinpath(temporary_root, "shared-domain-parent")
            child_root = joinpath(parent_root, "child")
            Mycelia._with_metamdbg_output_domain_lock(parent_root) do
                conflict = fetch(@async try
                    Mycelia._with_autocycler_output_lock(child_root) do _reserved
                        :unexpected
                    end
                catch caught
                    caught
                end)
                Test.@test conflict isa ArgumentError
                Test.@test occursin(
                    "active ancestor output-root reservation",
                    sprint(showerror, conflict),
                )
            end
            Mycelia._with_metamdbg_output_domain_lock(child_root) do
                conflict = fetch(@async try
                    Mycelia._with_unicycler_output_lock(parent_root) do _reserved
                        :unexpected
                    end
                catch caught
                    caught
                end)
                Test.@test conflict isa ArgumentError
                Test.@test occursin(
                    "active descendant output-root reservation",
                    sprint(showerror, conflict),
                )
            end

            sibling_meta = joinpath(temporary_root, "shared-domain-sibling-a")
            sibling_generic =
                joinpath(temporary_root, "shared-domain-sibling-b")
            Mycelia._with_metamdbg_output_domain_lock(sibling_meta) do
                result = fetch(@async Mycelia._with_unicycler_output_lock(
                    sibling_generic,
                ) do _reserved
                    :sibling_reserved
                end)
                Test.@test result == :sibling_reserved
            end

            input_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            function reserve_output_root!(
                    outdir::String,
                    owner_token::String,
            )::NamedTuple
                outputs = Mycelia._metamdbg_output_paths(outdir, 21)
                reservation = Mycelia._metamdbg_submission_reservation(
                    outputs,
                    input_contract,
                    21;
                    owner_token,
                )
                Mycelia._with_metamdbg_output_domain_lock(outdir) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        reservation,
                        outdir,
                    )
                end
                return reservation
            end

            queued_parent =
                joinpath(temporary_root, "shared-domain-queued-parent")
            queued_parent_reservation = reserve_output_root!(
                queued_parent,
                "shared-domain-queued-parent-owner",
            )
            Test.@test isfile(
                queued_parent_reservation.output_root_reservation_marker,
            )
            Test.@test Mycelia._output_root_reservation_is_active(
                queued_parent_reservation.output_root_reservation_marker,
            )
            _test_metamdbg_error(
                () -> Mycelia._with_unicycler_output_lock(queued_parent) do _
                    :unexpected
                end,
                ArgumentError,
                r"active same-root output-root reservation",
            )
            _test_metamdbg_error(
                () -> Mycelia._with_metamdbg_output_domain_lock(
                    joinpath(queued_parent, "meta-child"),
                ) do
                    :unexpected
                end,
                ArgumentError,
                r"active ancestor output-root reservation",
            )
            queued_sibling = joinpath(
                temporary_root,
                "shared-domain-queued-sibling",
            )
            Test.@test Mycelia._with_autocycler_output_lock(
                queued_sibling,
            ) do _reserved
                :queued_sibling_reserved
            end == :queued_sibling_reserved
            Mycelia._with_metamdbg_output_lock(queued_parent) do
                Mycelia._remove_metamdbg_submission_reservation!(
                    queued_parent_reservation,
                )
            end

            queued_child =
                joinpath(parent_root, "queued-meta-child")
            queued_child_reservation = reserve_output_root!(
                queued_child,
                "shared-domain-queued-child-owner",
            )
            _test_metamdbg_error(
                () -> Mycelia._with_metamdbg_output_domain_lock(
                    parent_root,
                ) do
                    :unexpected
                end,
                ArgumentError,
                r"active descendant output-root reservation",
            )
            Mycelia._with_metamdbg_output_lock(queued_child) do
                Mycelia._remove_metamdbg_submission_reservation!(
                    queued_child_reservation,
                )
            end
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
            Test.@test occursin("umask 077", script)
            Test.@test occursin("mkdir -m 700", script)
            Test.@test occursin("mktemp -d", script)
            Test.@test occursin("validate_fasta_stream", script)
            Test.@test occursin("validate_contigs", script)
            Test.@test occursin("validate_gfa", script)
            Test.@test !occursin("link_from", script)
            Test.@test !occursin("path_reference", script)
            Test.@test !occursin("object_count = split", script)
            Test.@test occursin("json.load", script)
            Test.@test !occursin("lines = list", script)
            Test.@test occursin(
                "for line_number, fields in records(path)",
                script,
            )
            Test.@test occursin("capture_package_inventory", script)
            Test.@test occursin(
                "resolved package inventory changed",
                script,
            )
            Test.@test !occursin("grep -Eq", script)
            Test.@test occursin("reservation_tombstone=", script)
            Test.@test occursin(
                "mv -- \"\$submission_reservation_dir\" " *
                "\"\$reservation_tombstone\"",
                script,
            )
            Test.@test !occursin(
                "rm -- \"\$submission_reservation_contract\"",
                script,
            )
            Test.@test occursin(Mycelia.METAMDBG_ENV_NAME, script)
            Test.@test occursin(Mycelia.METAMDBG_VERSION, script)
            Test.@test occursin(
                Mycelia.METAMDBG_ENVIRONMENT_SPEC_SHA256,
                script,
            )
            Test.@test occursin("mycelia_metamdbg_contract.json", script)
            Test.@test occursin("cmp -s", script)
            Test.@test occursin("mv -n -- \"\$contract_new\"", script)
            Test.@test occursin(input_contract.signature, script)
            Test.@test occursin("sha256_file", script)
            Test.@test occursin(
                only(input_contract.contract.inputs).sha256,
                script,
            )
            script_lines = split(script, '\n')
            lock_position = findfirst(==("lock_acquired=1"), script_lines)
            validation_positions = findall(
                ==("validate_metamdbg_inputs"),
                script_lines,
            )
            reuse_position = findfirst(==("contract_exists=0"), script_lines)
            finalization_position = findfirst(
                ==("if [ \"\$contract_exists\" -eq 1 ]; then"),
                script_lines,
            )
            Test.@test length(validation_positions) == 3
            Test.@test lock_position < first(validation_positions) <
                       reuse_position < validation_positions[2] <
                       last(validation_positions) <
                       finalization_position
            Test.@test !occursin("\$\$", script)
            Test.@test !occursin("\${contigs_gz}.tmp", script)
            Test.@test !occursin("rmdir -- \"\$lock_dir\" || true", script)

            script_path = joinpath(temporary_root, "metamdbg-executor.sh")
            write(script_path, script)
            Test.@test success(`bash -n $(script_path)`)

            fake_conda = joinpath(temporary_root, "fake-conda")
            write(
                fake_conda,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "if [ \"\$1\" = \"run\" ]; then\n" *
                "  shift\n" *
                "  while [ \"\$#\" -gt 0 ] && [ \"\$1\" != \"python\" ]; do shift; done\n" *
                "  [ \"\$#\" -gt 0 ] || exit 90\n" *
                "  shift\n" *
                "  exec python3 \"\$@\"\n" *
                "fi\n" *
                "[ \"\$1\" = \"list\" ] || exit 90\n" *
                "printf '%s\\n' " *
                "'[{\"name\":\"metamdbg\",\"version\":\"1.4\"," *
                "\"build_string\":\"h43eeafb_2\"," *
                "\"channel\":\"bioconda\"}," *
                "{\"name\":\"libzlib\",\"version\":\"1.3.1\"," *
                "\"build_string\":\"h8359307_2\"," *
                "\"channel\":\"conda-forge\"}]'\n",
            )
            chmod(fake_conda, 0o700)
            fake_asm =
                "printf '>contig-1\\nACGT\\n' > \"\$outdir/contigs.fasta\""
            fake_gfa =
                "printf 'H\\tVN:Z:1.0\\n" *
                "S\\tcontig-1\\tACGT\\tJA:J:{\"ok\":true}" *
                "\\tBc:B:c,-128,127\\tBI:B:I,0,4294967295" *
                "\\tBf:B:f,-1.5,3e2\\n" *
                "S\\tcontig-2\\tTGCA\\n" *
                "S\\tcontig,with,comma\\tGATTACA\\tLN:i:7\\n" *
                "L\\tcontig-1\\t+\\tcontig-2\\t-\\t4M\\n" *
                "P\\tpath-1\\tcontig-1+,contig-2-\\t4M\\n" *
                "P\\tcomma-path\\tcontig,with,comma+\\t*\\n' " *
                "> \"\$outdir/assemblyGraph_k21_4bps.gfa\""

            owner_runtime_outdir =
                joinpath(temporary_root, "executor-owner-capability")
            owner_runtime_outputs =
                Mycelia._metamdbg_output_paths(owner_runtime_outdir, 21)
            owner_runtime_reservation =
                Mycelia._metamdbg_submission_reservation(
                    owner_runtime_outputs,
                    input_contract,
                    21;
                    owner_token = "executor-owner-capability-fixture",
                )
            Mycelia._with_metamdbg_output_domain_lock(
                owner_runtime_outdir,
            ) do
                Mycelia._create_metamdbg_submission_reservation!(
                    owner_runtime_reservation,
                    owner_runtime_outdir,
                )
            end
            owner_runtime_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    owner_runtime_reservation,
                    "899",
                )
            foreign_reservation =
                Mycelia._output_root_durable_reservation_path_from_canonical(
                    owner_runtime_outdir,
                    "foreign-owner-proof",
                )
            write(foreign_reservation, "foreign owner\n")
            chmod(foreign_reservation, 0o600)
            owner_assembly_marker =
                joinpath(temporary_root, "owner-capability-assembly-ran")
            owner_runtime_script = Mycelia._metamdbg_executor_script(
                "touch $(Base.shell_escape(owner_assembly_marker)); $(fake_asm)",
                fake_gfa,
                owner_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = owner_runtime_reservation,
            )
            owner_runtime_script_path =
                joinpath(temporary_root, "executor-owner-capability.sh")
            write(owner_runtime_script_path, owner_runtime_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => owner_runtime_reservation.job_id,
            ) do
                success(`bash $(owner_runtime_script_path)`)
            end
            Test.@test isdir(owner_runtime_reservation.path)
            Test.@test isfile(
                owner_runtime_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                owner_runtime_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(owner_assembly_marker)
            rm(foreign_reservation)
            Test.@test Base.withenv(
                "SLURM_JOB_ID" => owner_runtime_reservation.job_id,
            ) do
                success(`bash $(owner_runtime_script_path)`)
            end
            Test.@test !ispath(owner_runtime_reservation.path)
            Test.@test !ispath(
                owner_runtime_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                owner_runtime_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test ispath(owner_assembly_marker)

            swapped_runtime_outdir =
                joinpath(temporary_root, "executor-output-root-swap")
            swapped_runtime_outputs =
                Mycelia._metamdbg_output_paths(swapped_runtime_outdir, 21)
            swapped_runtime_lock =
                Mycelia._metamdbg_output_lock_path(swapped_runtime_outdir)
            swapped_runtime_target =
                joinpath(temporary_root, "executor-output-root-swap-target")
            mkpath(swapped_runtime_target)
            swapped_runtime_reservation =
                Mycelia._metamdbg_submission_reservation(
                    swapped_runtime_outputs,
                    input_contract,
                    21;
                    owner_token = "executor-output-root-swap-fixture",
                )
            Mycelia._with_metamdbg_output_lock(swapped_runtime_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    swapped_runtime_reservation,
                    swapped_runtime_outdir,
                )
            end
            swapped_runtime_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    swapped_runtime_reservation,
                    "900",
                )
            swapped_runtime_asm =
                "rm -rf -- \"\$outdir\"; " *
                "ln -s -- $(Base.shell_escape(swapped_runtime_target)) " *
                "\"\$outdir\"; " *
                "printf '>escaped-contig\\nACGT\\n' > " *
                "\"\$outdir/contigs.fasta\""
            swapped_runtime_script = Mycelia._metamdbg_executor_script(
                swapped_runtime_asm,
                fake_gfa,
                swapped_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = swapped_runtime_reservation,
            )
            swapped_runtime_script_path =
                joinpath(temporary_root, "executor-output-root-swap.sh")
            write(swapped_runtime_script_path, swapped_runtime_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => swapped_runtime_reservation.job_id,
            ) do
                success(`bash $(swapped_runtime_script_path)`)
            end
            Test.@test islink(swapped_runtime_outdir)
            Test.@test isfile(joinpath(
                swapped_runtime_target,
                "contigs.fasta",
            ))
            Test.@test !ispath(swapped_runtime_outputs.contract_marker)
            Test.@test !ispath(swapped_runtime_outputs.completion_marker)
            Test.@test !ispath(swapped_runtime_reservation.path)
            Test.@test !ispath(swapped_runtime_lock)

            mixed_k_runtime_outdir =
                joinpath(temporary_root, "executor-mixed-graph-k")
            mixed_k_runtime_outputs =
                Mycelia._metamdbg_output_paths(mixed_k_runtime_outdir, 21)
            mixed_k_runtime_reservation =
                Mycelia._metamdbg_submission_reservation(
                    mixed_k_runtime_outputs,
                    input_contract,
                    21;
                    owner_token = "executor-mixed-graph-k-fixture",
                )
            Mycelia._with_metamdbg_output_lock(mixed_k_runtime_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    mixed_k_runtime_reservation,
                    mixed_k_runtime_outdir,
                )
            end
            mixed_k_runtime_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    mixed_k_runtime_reservation,
                    "901",
                )
            mixed_k_runtime_gfa = fake_gfa * "\n" *
                                  "printf 'H\\tVN:Z:1.0\\n" *
                                  "S\\tcontig-1\\tACGT\\n' " *
                                  "> \"\$outdir/" *
                                  "assemblyGraph_k31_4bps.gfa\""
            mixed_k_runtime_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                mixed_k_runtime_gfa,
                mixed_k_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = mixed_k_runtime_reservation,
            )
            mixed_k_runtime_script_path =
                joinpath(temporary_root, "executor-mixed-graph-k.sh")
            write(mixed_k_runtime_script_path, mixed_k_runtime_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => mixed_k_runtime_reservation.job_id,
            ) do
                success(`bash $(mixed_k_runtime_script_path)`)
            end
            Test.@test !ispath(mixed_k_runtime_outputs.graph_alias)
            Test.@test !ispath(mixed_k_runtime_outputs.contract_marker)
            Test.@test !ispath(mixed_k_runtime_outputs.completion_marker)
            Test.@test !ispath(mixed_k_runtime_reservation.path)

            unbound_runtime_outdir =
                joinpath(temporary_root, "executor-unbound-reservation")
            unbound_runtime_outputs =
                Mycelia._metamdbg_output_paths(unbound_runtime_outdir, 21)
            unbound_runtime_reservation =
                Mycelia._metamdbg_submission_reservation(
                    unbound_runtime_outputs,
                    input_contract,
                    21;
                    owner_token = "executor-unbound-reservation-fixture",
                )
            Mycelia._with_metamdbg_output_lock(unbound_runtime_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    unbound_runtime_reservation,
                    unbound_runtime_outdir,
                )
            end
            unbound_runtime_marker = joinpath(
                temporary_root,
                "unbound-runtime-assembly-ran",
            )
            unbound_runtime_script = Mycelia._metamdbg_executor_script(
                "touch $(Base.shell_escape(unbound_runtime_marker)); $(fake_asm)",
                fake_gfa,
                unbound_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = unbound_runtime_reservation,
            )
            unbound_runtime_script_path = joinpath(
                temporary_root,
                "executor-unbound-reservation.sh",
            )
            write(unbound_runtime_script_path, unbound_runtime_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => "998",
            ) do
                success(`bash $(unbound_runtime_script_path)`)
            end
            Test.@test isdir(unbound_runtime_reservation.path)
            Test.@test !ispath(unbound_runtime_marker)
            Test.@test !ispath(unbound_runtime_outputs.contract_marker)
            Mycelia._with_metamdbg_output_lock(unbound_runtime_outdir) do
                Mycelia._remove_metamdbg_submission_reservation!(
                    unbound_runtime_reservation,
                )
            end

            bound_runtime_cases = (
                matching = (job_id = "900", succeeds = true),
                missing = (job_id = nothing, succeeds = false),
                mismatched = (job_id = "999", succeeds = false),
            )
            for (case_name, runtime_case) in pairs(bound_runtime_cases)
                bound_runtime_outdir = joinpath(
                    temporary_root,
                    "executor-bound-job-$(case_name)",
                )
                bound_runtime_outputs =
                    Mycelia._metamdbg_output_paths(bound_runtime_outdir, 21)
                bound_runtime_reservation =
                    Mycelia._metamdbg_submission_reservation(
                        bound_runtime_outputs,
                        input_contract,
                        21;
                        owner_token = "executor-bound-job-$(case_name)",
                    )
                Mycelia._with_metamdbg_output_lock(bound_runtime_outdir) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        bound_runtime_reservation,
                        bound_runtime_outdir,
                    )
                end
                bound_runtime_reservation =
                    Mycelia._bind_metamdbg_submission_job!(
                        bound_runtime_reservation,
                        "900",
                    )
                bound_runtime_script = Mycelia._metamdbg_executor_script(
                    fake_asm,
                    fake_gfa,
                    bound_runtime_outputs,
                    21,
                    input_contract;
                    conda_runner = fake_conda,
                    submission_reservation = bound_runtime_reservation,
                )
                bound_runtime_script_path = joinpath(
                    temporary_root,
                    "executor-bound-job-$(case_name).sh",
                )
                write(bound_runtime_script_path, bound_runtime_script)
                runtime_succeeded = withenv(
                    "SLURM_JOB_ID" => runtime_case.job_id,
                ) do
                    success(`bash $(bound_runtime_script_path)`)
                end
                Test.@test runtime_succeeded == runtime_case.succeeds
                if runtime_case.succeeds
                    Test.@test !ispath(bound_runtime_reservation.path)
                    Test.@test isfile(bound_runtime_outputs.completion_marker)
                else
                    Test.@test isdir(bound_runtime_reservation.path)
                    Test.@test !ispath(
                        bound_runtime_outputs.completion_marker,
                    )
                    Mycelia._with_metamdbg_output_lock(
                        bound_runtime_outdir,
                    ) do
                        Mycelia._remove_metamdbg_submission_reservation!(
                            bound_runtime_reservation,
                        )
                    end
                end
            end

            partial_runtime_outdir =
                joinpath(temporary_root, "executor-partial-contracted")
            mkpath(partial_runtime_outdir)
            _write_test_metamdbg_gzip_fasta!(joinpath(
                partial_runtime_outdir,
                "contigs.fasta.gz",
            ))
            _write_test_metamdbg_contract!(
                partial_runtime_outdir,
                valid_reads,
            )
            partial_runtime_outputs = Mycelia._metamdbg_output_paths(
                partial_runtime_outdir,
                21,
            )
            partial_runtime_reservation =
                Mycelia._metamdbg_submission_reservation(
                    partial_runtime_outputs,
                    input_contract,
                    21;
                    owner_token = "executor-partial-contracted-fixture",
                )
            Mycelia._with_metamdbg_output_lock(partial_runtime_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    partial_runtime_reservation,
                    partial_runtime_outdir,
                )
            end
            partial_runtime_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    partial_runtime_reservation,
                    "902",
                )
            partial_runtime_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                partial_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = partial_runtime_reservation,
            )
            partial_runtime_script_path =
                joinpath(temporary_root, "executor-partial-contracted.sh")
            write(partial_runtime_script_path, partial_runtime_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => partial_runtime_reservation.job_id,
            ) do
                success(`bash $(partial_runtime_script_path)`)
            end
            Test.@test !ispath(partial_runtime_outputs.completion_marker)
            Test.@test !ispath(partial_runtime_reservation.path)

            retry_outdir = joinpath(temporary_root, "executor-lock-retry")
            retry_outputs = Mycelia._metamdbg_output_paths(retry_outdir, 21)
            retry_reservation = Mycelia._metamdbg_submission_reservation(
                retry_outputs,
                input_contract,
                21;
                owner_token = "executor-lock-retry-fixture",
            )
            Mycelia._with_metamdbg_output_lock(retry_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    retry_reservation,
                    retry_outdir,
                )
            end
            retry_reservation = Mycelia._bind_metamdbg_submission_job!(
                retry_reservation,
                "903",
            )
            retry_lock = Mycelia._metamdbg_output_lock_path(retry_outdir)
            mkdir(retry_lock)
            retry_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                retry_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = retry_reservation,
                lock_retry_attempts = 100,
                lock_retry_delay_seconds = 0.02,
            )
            Test.@test occursin("while ! mkdir", retry_script)
            retry_script_path =
                joinpath(temporary_root, "executor-lock-retry.sh")
            write(retry_script_path, retry_script)
            lock_releaser = run(
                `bash -c $("sleep 0.2; rmdir -- " * Base.shell_escape(retry_lock))`;
                wait = false,
            )
            Test.@test Base.withenv(
                "SLURM_JOB_ID" => retry_reservation.job_id,
            ) do
                success(`bash $(retry_script_path)`)
            end
            wait(lock_releaser)
            Test.@test success(lock_releaser)
            Test.@test !ispath(retry_reservation.path)
            Test.@test !ispath(retry_lock)
            Test.@test isfile(retry_outputs.contract_marker)

            timeout_outdir = joinpath(temporary_root, "executor-lock-timeout")
            timeout_outputs = Mycelia._metamdbg_output_paths(timeout_outdir, 21)
            timeout_reservation = Mycelia._metamdbg_submission_reservation(
                timeout_outputs,
                input_contract,
                21;
                owner_token = "executor-lock-timeout-fixture",
            )
            Mycelia._with_metamdbg_output_lock(timeout_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    timeout_reservation,
                    timeout_outdir,
                )
            end
            timeout_reservation = Mycelia._bind_metamdbg_submission_job!(
                timeout_reservation,
                "904",
            )
            timeout_lock = Mycelia._metamdbg_output_lock_path(timeout_outdir)
            mkdir(timeout_lock)
            timeout_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                timeout_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = timeout_reservation,
                lock_retry_attempts = 2,
                lock_retry_delay_seconds = 0,
            )
            timeout_script_path =
                joinpath(temporary_root, "executor-lock-timeout.sh")
            write(timeout_script_path, timeout_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => timeout_reservation.job_id,
            ) do
                success(`bash $(timeout_script_path)`)
            end
            Test.@test isdir(timeout_reservation.path)
            Test.@test isdir(timeout_lock)
            Test.@test !ispath(timeout_outputs.contract_marker)
            rm(timeout_lock)
            timeout_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    timeout_outdir,
                ),
            )
            timeout_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    timeout_metadata;
                    owner_token = timeout_metadata.owner_token,
                    job_id = timeout_metadata.job_id,
                    confirm_cancelled = true,
                )
            Test.@test timeout_reclaimed.status == :reclaimed
            Test.@test timeout_reclaimed.recovery_reason == :cancelled
            Test.@test !ispath(timeout_reservation.path)

            invalid_executor_fastas = (
                empty_identifier =
                    "printf '>\\nACGT\\n' > \"\$outdir/contigs.fasta\"",
                duplicate_identifier =
                    "printf '>duplicate first\\nACGT\\n" *
                    ">duplicate second\\nTGCA\\n' " *
                    "> \"\$outdir/contigs.fasta\"",
                whitespace_identifier =
                    "printf '> leading-space\\nACGT\\n' " *
                    "> \"\$outdir/contigs.fasta\"",
            )
            for (failure_index, (failure_type, invalid_asm)) in
                enumerate(pairs(invalid_executor_fastas))
                invalid_outdir = joinpath(
                    temporary_root,
                    "executor-fasta-$(failure_type)-output",
                )
                invalid_outputs =
                    Mycelia._metamdbg_output_paths(invalid_outdir, 21)
                invalid_reservation =
                    Mycelia._metamdbg_submission_reservation(
                        invalid_outputs,
                        input_contract,
                        21;
                        owner_token = "executor-fasta-$(failure_type)-fixture",
                    )
                Mycelia._with_metamdbg_output_lock(invalid_outdir) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        invalid_reservation,
                        invalid_outdir,
                    )
                end
                invalid_reservation =
                    Mycelia._bind_metamdbg_submission_job!(
                        invalid_reservation,
                        string(910 + failure_index),
                    )
                invalid_script = Mycelia._metamdbg_executor_script(
                    invalid_asm,
                    fake_gfa,
                    invalid_outputs,
                    21,
                    input_contract;
                    conda_runner = fake_conda,
                    submission_reservation = invalid_reservation,
                )
                invalid_script_path = joinpath(
                    temporary_root,
                    "executor-fasta-$(failure_type).sh",
                )
                write(invalid_script_path, invalid_script)
                Test.@test !Base.withenv(
                    "SLURM_JOB_ID" => invalid_reservation.job_id,
                ) do
                    success(`bash $(invalid_script_path)`)
                end
                Test.@test !ispath(invalid_outputs.contract_marker)
                Test.@test !ispath(invalid_reservation.path)
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(invalid_outdir),
                )
            end

            invalid_executor_gfas = (
                leading_garbage =
                    "printf 'garbage\\nS\\tcontig-1\\tACGT\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                malformed_link =
                    "printf 'S\\tcontig-1\\tACGT\\n" *
                    "S\\tcontig-2\\tTGCA\\n" *
                    "L\\tcontig-1\\t?\\tcontig-2\\t+\\t4M\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                dangling_link =
                    "printf 'S\\tcontig-1\\tACGT\\n" *
                    "L\\tcontig-1\\t+\\tmissing\\t+\\t4M\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                leading_star =
                    "printf 'S\\t*segment\\tACGT\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                leading_equals =
                    "printf 'S\\t=segment\\tACGT\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                shared_namespace =
                    "printf 'S\\tshared-name\\tACGT\\n" *
                    "P\\tshared-name\\tshared-name+\\t*\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                duplicate_tag =
                    "printf 'S\\tcontig-1\\tACGT\\tLN:i:4\\tLN:i:5\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                invalid_integer_tag =
                    "printf 'S\\tcontig-1\\tACGT\\t" *
                    "LN:i:not-an-integer\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                empty_z_tag =
                    "printf 'S\\tcontig-1\\tACGT\\tZZ:Z:\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                invalid_json_tag =
                    "printf 'S\\tcontig-1\\tACGT\\tJJ:J:{bad}\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                nonstandard_json_constant =
                    "printf 'S\\tcontig-1\\tACGT\\tJJ:J:NaN\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                b_signed_range =
                    "printf 'S\\tcontig-1\\tACGT\\tBB:B:c,-129\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                b_unsigned_range =
                    "printf 'S\\tcontig-1\\tACGT\\tBB:B:I,4294967296\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                b_float_range =
                    "printf 'S\\tcontig-1\\tACGT\\tBB:B:f,1e400\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                non_gfa1_version =
                    "printf 'H\\tVN:Z:2.0\\nS\\tcontig-1\\tACGT\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
                leading_whitespace_pseudo_comment =
                    "printf ' #not-a-comment\\n" *
                    "S\\tcontig-1\\tACGT\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
            )
            for (failure_index, (failure_type, invalid_gfa)) in
                enumerate(pairs(invalid_executor_gfas))
                invalid_outdir = joinpath(
                    temporary_root,
                    "executor-$(failure_type)-output",
                )
                invalid_outputs =
                    Mycelia._metamdbg_output_paths(invalid_outdir, 21)
                invalid_reservation =
                    Mycelia._metamdbg_submission_reservation(
                        invalid_outputs,
                        input_contract,
                        21;
                        owner_token = "executor-$(failure_type)-fixture",
                    )
                Mycelia._with_metamdbg_output_lock(invalid_outdir) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        invalid_reservation,
                        invalid_outdir,
                    )
                end
                invalid_reservation =
                    Mycelia._bind_metamdbg_submission_job!(
                        invalid_reservation,
                        string(930 + failure_index),
                    )
                invalid_script = Mycelia._metamdbg_executor_script(
                    fake_asm,
                    invalid_gfa,
                    invalid_outputs,
                    21,
                    input_contract;
                    conda_runner = fake_conda,
                    submission_reservation = invalid_reservation,
                )
                invalid_script_path = joinpath(
                    temporary_root,
                    "executor-$(failure_type).sh",
                )
                write(invalid_script_path, invalid_script)
                Test.@test !Base.withenv(
                    "SLURM_JOB_ID" => invalid_reservation.job_id,
                ) do
                    success(`bash $(invalid_script_path)`)
                end
                Test.@test !ispath(invalid_outputs.contract_marker)
                Test.@test !ispath(invalid_reservation.path)
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(invalid_outdir),
                )
            end

            for (queued_index, queued_change) in
                enumerate((:mutation, :replacement))
                queued_reads = joinpath(
                    temporary_root,
                    "queued-$(queued_change).fastq",
                )
                write(queued_reads, "@queued\nACGT\n+\nIIII\n")
                queued_mutation_target = queued_reads
                queued_input = if queued_change == :mutation
                    queued_secondary = joinpath(
                        temporary_root,
                        "queued-mutation-secondary.fastq",
                    )
                    write(queued_secondary, "@secondary\nACGT\n+\nIIII\n")
                    queued_mutation_target = queued_secondary
                    String[queued_reads, queued_secondary]
                else
                    queued_reads
                end
                queued_outdir = joinpath(
                    temporary_root,
                    "queued-$(queued_change)-output",
                )
                queued_outputs =
                    Mycelia._metamdbg_output_paths(queued_outdir, 21)
                queued_contract = Mycelia._metamdbg_input_contract(
                    Mycelia._metamdbg_selected_input(queued_input, nothing),
                    3,
                )
                queued_reservation = Mycelia._metamdbg_submission_reservation(
                    queued_outputs,
                    queued_contract,
                    21;
                    owner_token = "queued-$(queued_change)-fixture",
                )
                Mycelia._with_metamdbg_output_lock(queued_outdir) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        queued_reservation,
                        queued_outdir,
                    )
                end
                queued_reservation = Mycelia._bind_metamdbg_submission_job!(
                    queued_reservation,
                    string(950 + queued_index),
                )
                assembly_marker = joinpath(
                    temporary_root,
                    "queued-$(queued_change)-assembly-ran",
                )
                queued_asm =
                    "touch $(Base.shell_escape(assembly_marker)); $(fake_asm)"
                queued_script = Mycelia._metamdbg_executor_script(
                    queued_asm,
                    fake_gfa,
                    queued_outputs,
                    21,
                    queued_contract;
                    conda_runner = fake_conda,
                    submission_reservation = queued_reservation,
                )
                if queued_change == :mutation
                    write(
                        queued_mutation_target,
                        "@secondary\nTGCA\n+\nIIII\n",
                    )
                else
                    replacement_reads = joinpath(
                        temporary_root,
                        "queued-replacement-new.fastq",
                    )
                    write(replacement_reads, "@queued\nTGCA\n+\nIIII\n")
                    mv(replacement_reads, queued_reads; force = true)
                end
                queued_script_path = joinpath(
                    temporary_root,
                    "queued-$(queued_change)-executor.sh",
                )
                write(queued_script_path, queued_script)
                Test.@test !Base.withenv(
                    "SLURM_JOB_ID" => queued_reservation.job_id,
                ) do
                    success(`bash $(queued_script_path)`)
                end
                Test.@test !ispath(assembly_marker)
                Test.@test !ispath(queued_outputs.contract_marker)
                Test.@test !ispath(queued_reservation.path)
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(queued_outdir),
                )
            end

            runtime_mutation_reads =
                joinpath(temporary_root, "runtime-mutation.fastq")
            write(
                runtime_mutation_reads,
                "@runtime-mutation\nACGT\n+\nIIII\n",
            )
            runtime_mutation_outdir =
                joinpath(temporary_root, "runtime-mutation-output")
            runtime_mutation_outputs =
                Mycelia._metamdbg_output_paths(runtime_mutation_outdir, 21)
            runtime_mutation_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(
                    runtime_mutation_reads,
                    nothing,
                ),
                3,
            )
            runtime_mutation_reservation =
                Mycelia._metamdbg_submission_reservation(
                    runtime_mutation_outputs,
                    runtime_mutation_contract,
                    21;
                    owner_token = "runtime-mutation-fixture",
                )
            Mycelia._with_metamdbg_output_lock(
                runtime_mutation_outdir,
            ) do
                Mycelia._create_metamdbg_submission_reservation!(
                    runtime_mutation_reservation,
                    runtime_mutation_outdir,
                )
            end
            runtime_mutation_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    runtime_mutation_reservation,
                    "905",
                )
            runtime_mutation_asm =
                fake_asm * "; printf '%s' " *
                Base.shell_escape("@runtime-mutation\nTGCA\n+\nIIII\n") *
                " > " * Base.shell_escape(runtime_mutation_reads)
            runtime_mutation_script = Mycelia._metamdbg_executor_script(
                runtime_mutation_asm,
                fake_gfa,
                runtime_mutation_outputs,
                21,
                runtime_mutation_contract;
                conda_runner = fake_conda,
                submission_reservation = runtime_mutation_reservation,
            )
            runtime_mutation_script_path = joinpath(
                temporary_root,
                "runtime-mutation-executor.sh",
            )
            write(runtime_mutation_script_path, runtime_mutation_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => runtime_mutation_reservation.job_id,
            ) do
                success(`bash $(runtime_mutation_script_path)`)
            end
            Test.@test !ispath(runtime_mutation_reservation.path)
            Test.@test !ispath(runtime_mutation_outputs.contract_marker)
            Test.@test !ispath(runtime_mutation_outputs.completion_marker)
            Test.@test !ispath(Mycelia._metamdbg_output_lock_path(
                runtime_mutation_outdir,
            ))

            missing_reservation_outdir =
                joinpath(temporary_root, "missing-reservation-output")
            missing_reservation_outputs =
                Mycelia._metamdbg_output_paths(missing_reservation_outdir, 21)
            missing_reservation = Mycelia._metamdbg_submission_reservation(
                missing_reservation_outputs,
                input_contract,
                21;
                owner_token = "missing-reservation-fixture",
            )
            missing_reservation_assembly_marker = joinpath(
                temporary_root,
                "missing-reservation-assembly-ran",
            )
            missing_reservation_script = Mycelia._metamdbg_executor_script(
                "touch " *
                Base.shell_escape(missing_reservation_assembly_marker) *
                "; " * fake_asm,
                fake_gfa,
                missing_reservation_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = missing_reservation,
            )
            missing_reservation_script_path = joinpath(
                temporary_root,
                "missing-reservation-executor.sh",
            )
            write(
                missing_reservation_script_path,
                missing_reservation_script,
            )
            Test.@test !success(`bash $(missing_reservation_script_path)`)
            Test.@test !ispath(missing_reservation_assembly_marker)
            Test.@test !ispath(missing_reservation_outputs.contract_marker)
            Test.@test !ispath(Mycelia._metamdbg_output_lock_path(
                missing_reservation_outdir,
            ))

            executable_outdir = joinpath(temporary_root, "executor-fixture")
            executable_outputs =
                Mycelia._metamdbg_output_paths(executable_outdir, 21)
            executable_reservation = Mycelia._metamdbg_submission_reservation(
                executable_outputs,
                input_contract,
                21;
                owner_token = "successful-executor-fixture",
            )
            Mycelia._with_metamdbg_output_lock(executable_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    executable_reservation,
                    executable_outdir,
                )
            end
            executable_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    executable_reservation,
                    "906",
                )
            executable_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                executable_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = executable_reservation,
            )
            executable_script_path =
                joinpath(temporary_root, "metamdbg-executable-fixture.sh")
            write(executable_script_path, executable_script)
            Test.@test Base.withenv(
                "SLURM_JOB_ID" => executable_reservation.job_id,
            ) do
                success(`bash $(executable_script_path)`)
            end
            Test.@test Mycelia._require_valid_metamdbg_fasta(
                executable_outputs.contigs_gz,
                "executor fixture contigs",
            ) == executable_outputs.contigs_gz
            graph_source = joinpath(
                executable_outdir,
                "assemblyGraph_k21_4bps.gfa",
            )
            Test.@test Mycelia._require_valid_metamdbg_gfa(
                graph_source,
                "executor fixture graph",
            ) == graph_source
            Test.@test read(executable_outputs.contract_marker, String) ==
                       input_contract.contents
            executable_artifacts =
                Mycelia._require_metamdbg_artifacts!(executable_outputs, 21)
            executable_completion =
                Mycelia._require_metamdbg_completion_manifest!(
                    executable_outputs,
                    executable_artifacts,
                    input_contract,
                    21,
                )
            Test.@test executable_completion.manifest.toolchain.package_inventory_sha256 ==
                       _test_metamdbg_toolchain().package_inventory_sha256
            Test.@test !ispath(executable_reservation.path)
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(executable_outdir),
            )
            Test.@test isempty(filter(
                name -> startswith(
                    name,
                    ".$(basename(executable_outdir)).mycelia.",
                ),
                readdir(dirname(executable_outdir)),
            ))

            drift_runtime_outdir =
                joinpath(temporary_root, "executor-inventory-drift")
            drift_runtime_outputs = Mycelia._metamdbg_output_paths(
                drift_runtime_outdir,
                21,
            )
            drift_runtime_reservation =
                Mycelia._metamdbg_submission_reservation(
                    drift_runtime_outputs,
                    input_contract,
                    21;
                    owner_token = "executor-inventory-drift-fixture",
                )
            Mycelia._with_metamdbg_output_lock(drift_runtime_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    drift_runtime_reservation,
                    drift_runtime_outdir,
                )
            end
            drift_runtime_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    drift_runtime_reservation,
                    "907",
                )
            drift_conda = joinpath(temporary_root, "drifting-fake-conda")
            drift_state = joinpath(temporary_root, "drifting-conda-state")
            write(
                drift_conda,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "if [ \"\$1\" = \"run\" ]; then exec $(Base.shell_escape(fake_conda)) \"\$@\"; fi\n" *
                "count=0; [ ! -f $(Base.shell_escape(drift_state)) ] || count=\$(cat $(Base.shell_escape(drift_state)))\n" *
                "count=\$((count + 1)); printf '%s' \"\$count\" > $(Base.shell_escape(drift_state))\n" *
                "build=h8359307_2; [ \"\$count\" -eq 1 ] || build=h8359307_3\n" *
                "printf '[{\"name\":\"metamdbg\",\"version\":\"1.4\",\"build_string\":\"h43eeafb_2\",\"channel\":\"bioconda\"},{\"name\":\"libzlib\",\"version\":\"1.3.1\",\"build_string\":\"%s\",\"channel\":\"conda-forge\"}]\\n' \"\$build\"\n",
            )
            chmod(drift_conda, 0o700)
            drift_runtime_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                drift_runtime_outputs,
                21,
                input_contract;
                conda_runner = drift_conda,
                submission_reservation = drift_runtime_reservation,
            )
            drift_runtime_script_path =
                joinpath(temporary_root, "drifting-runtime.sh")
            write(drift_runtime_script_path, drift_runtime_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => drift_runtime_reservation.job_id,
            ) do
                success(`bash $(drift_runtime_script_path)`)
            end
            Test.@test !ispath(drift_runtime_outputs.contract_marker)
            Test.@test !ispath(drift_runtime_outputs.completion_marker)
            Test.@test !ispath(drift_runtime_reservation.path)

            publication_race_outdir =
                joinpath(temporary_root, "executor-publication-race")
            publication_race_outputs = Mycelia._metamdbg_output_paths(
                publication_race_outdir,
                21,
            )
            publication_race_reservation =
                Mycelia._metamdbg_submission_reservation(
                    publication_race_outputs,
                    input_contract,
                    21;
                    owner_token = "executor-publication-race-fixture",
                )
            Mycelia._with_metamdbg_output_lock(publication_race_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    publication_race_reservation,
                    publication_race_outdir,
                )
            end
            publication_race_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    publication_race_reservation,
                    "908",
                )
            publication_race_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                publication_race_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = publication_race_reservation,
                post_completion_publication_hook =
                    "printf '>contig-1\\nTGCA\\n' | gzip -c > \"\$contigs_gz\"",
            )
            publication_race_script_path =
                joinpath(temporary_root, "publication-race.sh")
            write(publication_race_script_path, publication_race_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => publication_race_reservation.job_id,
            ) do
                success(`bash $(publication_race_script_path)`)
            end
            Test.@test !ispath(publication_race_outputs.completion_marker)
            Test.@test isfile(publication_race_outputs.contract_marker)
            Test.@test !ispath(publication_race_reservation.path)

            sentinel_victim = joinpath(temporary_root, "sentinel-victim")
            write(sentinel_victim, "must-remain-unchanged\n")
            adversarial_reads = joinpath(
                temporary_root,
                "reads;printf PWNED>sentinel-victim;#'\"\$(x).fastq",
            )
            write(adversarial_reads, "@adversarial\nACGT\n+\nIIII\n")
            sentinel_outdir = joinpath(
                temporary_root,
                "executor;printf PWNED>sentinel-victim;#'\"\$(x)",
            )
            sentinel_outputs =
                Mycelia._metamdbg_output_paths(sentinel_outdir, 21)
            sentinel_input_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(adversarial_reads, nothing),
                3,
            )
            sentinel_reservation = Mycelia._metamdbg_submission_reservation(
                sentinel_outputs,
                sentinel_input_contract,
                21;
                owner_token = "sentinel-executor-fixture",
            )
            Mycelia._with_metamdbg_output_lock(sentinel_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    sentinel_reservation,
                    sentinel_outdir,
                )
            end
            sentinel_reservation = Mycelia._bind_metamdbg_submission_job!(
                sentinel_reservation,
                "909",
            )
            sentinel_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                sentinel_outputs,
                21,
                sentinel_input_contract;
                conda_runner = fake_conda,
                submission_reservation = sentinel_reservation,
            )
            sentinel_script_path =
                joinpath(temporary_root, "metamdbg-sentinel-fixture.sh")
            write(sentinel_script_path, sentinel_script)
            sentinel_log = joinpath(temporary_root, "sentinel-runtime.log")
            sentinel_success = Base.withenv(
                "SLURM_JOB_ID" => sentinel_reservation.job_id,
            ) do
                open(sentinel_log, "w") do log_io
                    success(pipeline(
                        Cmd(
                            `bash $(sentinel_script_path)`;
                            dir = temporary_root,
                        );
                        stderr = log_io,
                    ))
                end
            end
            sentinel_success || @info "sentinel runtime failure" log = read(
                sentinel_log,
                String,
            )
            Test.@test sentinel_success
            Test.@test read(sentinel_victim, String) ==
                       "must-remain-unchanged\n"
            Test.@test Mycelia._require_valid_metamdbg_fasta(
                sentinel_outputs.contigs_gz,
                "sentinel fixture contigs",
            ) == sentinel_outputs.contigs_gz
            Test.@test read(sentinel_outputs.contract_marker, String) ==
                       sentinel_input_contract.contents
            sentinel_artifacts =
                Mycelia._require_metamdbg_artifacts!(sentinel_outputs, 21)
            Test.@test Mycelia._require_metamdbg_completion_manifest!(
                sentinel_outputs,
                sentinel_artifacts,
                sentinel_input_contract,
                21,
            ).manifest.workflow.graph_k == 21
            Test.@test !ispath(sentinel_reservation.path)
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(sentinel_outdir),
            )
        end

        Test.@testset "collected executor reports normalized planned paths" begin
            outdir = joinpath(temporary_root, "collected")
            executor = Mycelia.CollectExecutor()
            result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir,
                graph_k = 31,
                executor,
                dependency_checker = () -> error(
                    "planned collection must not provision metaMDBG",
                ),
                local_runner = forbidden_runner,
            )
            Test.@test length(executor.jobs) == 1
            Test.@test result.status == :planned
            Test.@test result.submission == 1
            Test.@test !hasproperty(result, :contigs)
            Test.@test !hasproperty(result, :graph)
            Test.@test result.expected_artifacts.contigs ==
                       joinpath(outdir, "contigs.fasta.gz")
            Test.@test result.expected_artifacts.graph ==
                       joinpath(outdir, "assemblyGraph_k31.gfa")
            Test.@test result.expected_artifacts.completion_marker ==
                       joinpath(outdir, "mycelia_metamdbg_completion.json")
            Test.@test result.provenance.graph_k == 31
            alternate_graph_k = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir,
                graph_k = 32,
                executor = Mycelia.CollectExecutor(),
                dependency_checker = () -> error(
                    "planned collection must not provision metaMDBG",
                ),
                local_runner = forbidden_runner,
            )
            Test.@test alternate_graph_k.provenance.contract_signature ==
                       result.provenance.contract_signature
            Test.@test alternate_graph_k.provenance.workflow_signature !=
                       result.provenance.workflow_signature
            Test.@test alternate_graph_k.provenance.graph_k == 32
            Test.@test occursin(
                "assemblyGraph_k31_*bps.gfa",
                only(executor.jobs).cmd,
            )
            Test.@test occursin("validate_contigs", only(executor.jobs).cmd)
            Test.@test occursin("validate_gfa", only(executor.jobs).cmd)
            Test.@test occursin(
                "mycelia_metamdbg_contract.json",
                only(executor.jobs).cmd,
            )
            Test.@test occursin(
                "mycelia_metamdbg_completion.json",
                only(executor.jobs).cmd,
            )
            Test.@test occursin("cmp -s", only(executor.jobs).cmd)
            Test.@test isempty(
                Mycelia._metamdbg_submission_reservation_paths(outdir),
            )
            Test.@test !ispath(joinpath(
                outdir,
                "mycelia_metamdbg_contract.json",
            ))
            Test.@test occursin(
                "submission_reservation_dir=",
                only(executor.jobs).cmd,
            )
            collected_script_path = joinpath(
                temporary_root,
                "collected-unreserved-executor.sh",
            )
            write(collected_script_path, only(executor.jobs).cmd)
            Test.@test !success(`bash $(collected_script_path)`)
            Test.@test !ispath(Mycelia._metamdbg_output_lock_path(outdir))
            Test.@test !ispath(joinpath(
                outdir,
                "mycelia_metamdbg_contract.json",
            ))
        end

        Test.@testset "dry-run executors never persist reservations" begin
            planned_executors = (
                dry_run = Mycelia.DryRunExecutor(),
                slurm_dry_run = Mycelia.SlurmExecutor(dry_run = true),
            )
            for (label, executor) in pairs(planned_executors)
                outdir = joinpath(temporary_root, "$(label)-output")
                captured_job = Ref{Any}(nothing)
                result = Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    executor,
                    dependency_checker = () -> error(
                        "dry-run execution must not provision metaMDBG",
                    ),
                    local_runner = forbidden_runner,
                    submission_runner = function (
                            job::Mycelia.JobSpec,
                            _executor::Mycelia.AbstractExecutor,
                    )
                        captured_job[] = job
                        return :planned_fixture
                    end,
                )
                Test.@test result.status == :planned
                Test.@test result.submission == :planned_fixture
                Test.@test captured_job[] !== nothing
                Test.@test occursin(
                    "submission_reservation_dir=",
                    captured_job[].cmd,
                )
                Test.@test isempty(
                    Mycelia._metamdbg_submission_reservation_paths(outdir),
                )
                Test.@test !ispath(joinpath(
                    outdir,
                    "mycelia_metamdbg_contract.json",
                ))
            end
        end

        Test.@testset "real submissions reserve exactly once" begin
            outdir = joinpath(temporary_root, "submitted-output")
            submission_calls = Ref(0)
            reserved_path = Ref("")
            lifecycle_events = Symbol[]
            submission_runner = function (
                    job::Mycelia.JobSpec,
                    executor::Mycelia.AbstractExecutor,
            )
                push!(lifecycle_events, :submit_held)
                submission_calls[] += 1
                Test.@test executor isa Mycelia.SlurmExecutor
                Test.@test executor.hold
                Test.@test !executor.dry_run
                reservations =
                    Mycelia._metamdbg_submission_reservation_paths(outdir)
                Test.@test length(reservations) == 1
                reserved_path[] = only(reservations)
                Test.@test isdir(reserved_path[])
                Test.@test occursin(
                    Base.shell_escape(reserved_path[]),
                    job.cmd,
                )
                reservation_record = JSON.parse(read(
                    joinpath(reserved_path[], "contract.json"),
                    String,
                ))
                Test.@test job.job_name ==
                           reservation_record["scheduler_job_name"]
                Test.@test startswith(job.job_name, "custom-job-name-")
                Test.@test endswith(
                    job.job_name,
                    "-mycelia-metamdbg-" *
                    reservation_record["workflow_signature"],
                )
                return Mycelia.SubmitResult(
                    ok = true,
                    dry_run = false,
                    held = true,
                    scheduler_acceptance = :accepted,
                    site = :scg,
                    backend = :sbatch,
                    job_id = "123",
                )
            end
            submission_job_binder = function (
                    reservation::NamedTuple,
                    job_id::AbstractString,
            )
                push!(lifecycle_events, :bind_job_record)
                Test.@test lifecycle_events ==
                           Symbol[:submit_held, :bind_job_record]
                return Mycelia._bind_metamdbg_submission_job_after_submit!(
                    reservation,
                    job_id,
                )
            end
            submission_release_runner = function (job_id::AbstractString)
                push!(lifecycle_events, :release_exact_job)
                Test.@test lifecycle_events == Symbol[
                    :submit_held,
                    :bind_job_record,
                    :release_exact_job,
                ]
                Test.@test String(job_id) == "123"
                Test.@test isfile(joinpath(reserved_path[], "job.json"))
                return String(job_id)
            end
            result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir,
                job_name = " custom job/name ",
                executor = Mycelia.SlurmExecutor(dry_run = false),
                site = :scg,
                dependency_checker = _test_metamdbg_toolchain,
                local_runner = forbidden_runner,
                submission_runner,
                submission_job_binder,
                submission_release_runner,
            )
            Test.@test result.status == :submitted
            Test.@test result.submission.job_id == "123"
            Test.@test result.submission.held
            reservation_metadata = result.submission_reservation
            Test.@test reservation_metadata.path == reserved_path[]
            Test.@test reservation_metadata.job_id == "123"
            Test.@test reservation_metadata.submission_state == :submitted
            Test.@test !isempty(reservation_metadata.owner_token)
            Test.@test occursin(
                reservation_metadata.workflow_signature,
                reservation_metadata.path,
            )
            Test.@test submission_calls[] == 1
            Test.@test lifecycle_events == Symbol[
                :submit_held,
                :bind_job_record,
                :release_exact_job,
            ]
            Test.@test isdir(reserved_path[])
            Test.@test readdir(reserved_path[]) ==
                       String["contract.json", "job.json"]
            submitted_job_marker = joinpath(reserved_path[], "job.json")
            Test.@test (stat(submitted_job_marker).mode & 0o777) == 0o600
            submitted_job_record = JSON.parse(read(
                submitted_job_marker,
                String,
            ))
            Test.@test submitted_job_record["schema_version"] ==
                       Mycelia._METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION
            Test.@test submitted_job_record["workflow_signature"] ==
                       reservation_metadata.workflow_signature
            Test.@test submitted_job_record["owner_token"] ==
                       reservation_metadata.owner_token
            Test.@test submitted_job_record["job_id"] == "123"
            inspected_submission = only(
                Mycelia.inspect_metamdbg_submission_reservations(outdir),
            )
            Test.@test inspected_submission.job_id == "123"
            Test.@test inspected_submission.submission_state == :submitted

            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner,
                    submission_job_binder,
                    submission_release_runner,
                ),
                ErrorException,
                r"active nonlocal submission reservation",
            )
            Test.@test submission_calls[] == 1
            Test.@test isdir(reserved_path[])
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    graph_k = 31,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner,
                    submission_job_binder,
                    submission_release_runner,
                ),
                ErrorException,
                r"active nonlocal submission reservation",
            )
            Test.@test submission_calls[] == 1
            Test.@test isdir(reserved_path[])

            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    reservation_metadata;
                    owner_token = reservation_metadata.owner_token,
                    job_id = reservation_metadata.job_id,
                ),
                ArgumentError,
                r"exactly one of confirm_cancelled=true",
            )
            Test.@test isdir(reserved_path[])
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    reservation_metadata;
                    owner_token = "wrong-owner-token",
                    job_id = reservation_metadata.job_id,
                    confirm_cancelled = true,
                ),
                ErrorException,
                r"owner token does not match",
            )
            Test.@test isdir(reserved_path[])
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    reservation_metadata;
                    owner_token = reservation_metadata.owner_token,
                    job_id = "124",
                    confirm_cancelled = true,
                ),
                ErrorException,
                r"job id does not match",
            )
            Test.@test isdir(reserved_path[])
            reclaimed = Mycelia.reclaim_metamdbg_submission_reservation!(
                reservation_metadata;
                owner_token = reservation_metadata.owner_token,
                job_id = reservation_metadata.job_id,
                confirm_cancelled = true,
            )
            Test.@test reclaimed.status == :reclaimed
            Test.@test reclaimed.job_id == "123"
            Test.@test reclaimed.recovery_reason == :cancelled
            Test.@test !ispath(reserved_path[])

            terminal_outdir =
                joinpath(temporary_root, "terminal-failed-submission")
            terminal_outputs =
                Mycelia._metamdbg_output_paths(terminal_outdir, 21)
            terminal_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            terminal_reservation = Mycelia._metamdbg_submission_reservation(
                terminal_outputs,
                terminal_contract,
                21;
                owner_token = "terminal-failed-owner",
            )
            Mycelia._with_metamdbg_output_lock(terminal_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    terminal_reservation,
                    terminal_outdir,
                )
            end
            terminal_bound_reservation =
                Mycelia._bind_metamdbg_submission_job!(
                    terminal_reservation,
                    "789",
                )
            terminal_metadata = (;
                canonical_outdir = terminal_bound_reservation.canonical_outdir,
                path = terminal_bound_reservation.path,
                workflow_signature =
                    terminal_bound_reservation.workflow_signature,
                input_contract_signature =
                    terminal_bound_reservation.input_contract_signature,
                graph_k = terminal_bound_reservation.graph_k,
                owner_token = terminal_bound_reservation.owner_token,
                job_id = terminal_bound_reservation.job_id,
            )
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    terminal_metadata;
                    owner_token = terminal_metadata.owner_token,
                    job_id = terminal_metadata.job_id,
                    confirm_terminal = :completed,
                ),
                ArgumentError,
                r"accepts only :failed",
            )
            Test.@test isdir(terminal_reservation.path)
            terminal_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    terminal_metadata;
                    owner_token = terminal_metadata.owner_token,
                    job_id = terminal_metadata.job_id,
                    confirm_terminal = :failed,
                )
            Test.@test terminal_reclaimed.status == :reclaimed
            Test.@test terminal_reclaimed.job_id == "789"
            Test.@test terminal_reclaimed.recovery_reason == :failed
            Test.@test !ispath(terminal_reservation.path)

            crashed_outdir =
                joinpath(temporary_root, "pre-submit-process-death")
            crashed_outputs =
                Mycelia._metamdbg_output_paths(crashed_outdir, 21)
            crashed_input_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            crashed_reservation = Mycelia._metamdbg_submission_reservation(
                crashed_outputs,
                crashed_input_contract,
                21;
                owner_token = "pre-submit-owner-capability",
            )
            Mycelia._with_metamdbg_output_lock(crashed_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    crashed_reservation,
                    crashed_outdir,
                )
            end
            Test.@test (stat(crashed_reservation.path).mode & 0o777) == 0o700
            Test.@test (
                stat(crashed_reservation.contract_marker).mode & 0o777
            ) == 0o600
            chmod(crashed_reservation.contract_marker, 0o640)
            _test_metamdbg_error(
                () -> Mycelia.inspect_metamdbg_submission_reservations(
                    crashed_outdir,
                ),
                ErrorException,
                r"must have mode 0600",
            )
            chmod(crashed_reservation.contract_marker, 0o600)
            run(`touch -t 200001010101 $(crashed_reservation.contract_marker)`)
            inspected =
                Mycelia.inspect_metamdbg_submission_reservations(crashed_outdir)
            Test.@test length(inspected) == 1
            inspected_reservation = only(inspected)
            Test.@test inspected_reservation.path == crashed_reservation.path
            Test.@test inspected_reservation.owner_token ==
                       crashed_reservation.owner_token
            Test.@test inspected_reservation.job_id === nothing
            Test.@test inspected_reservation.submission_state == :reserved
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    inspected_reservation;
                    owner_token = "wrong-pre-submit-token",
                    confirm_not_submitted = true,
                ),
                ErrorException,
                r"owner token does not match",
            )
            Test.@test isdir(crashed_reservation.path)
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    inspected_reservation;
                    owner_token = inspected_reservation.owner_token,
                ),
                ArgumentError,
                r"exactly one of confirm_cancelled=true",
            )
            Test.@test isdir(crashed_reservation.path)
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    inspected_reservation;
                    owner_token = inspected_reservation.owner_token,
                    job_id = "777",
                    confirm_not_submitted = true,
                ),
                ArgumentError,
                r"Do not provide a scheduler job_id",
            )
            Test.@test isdir(crashed_reservation.path)
            recovered = Mycelia.reclaim_metamdbg_submission_reservation!(
                inspected_reservation;
                owner_token = inspected_reservation.owner_token,
                confirm_not_submitted = true,
            )
            Test.@test recovered.status == :reclaimed
            Test.@test recovered.job_id === nothing
            Test.@test recovered.recovery_reason == :not_submitted
            Test.@test !ispath(crashed_reservation.path)
            Test.@test isempty(
                Mycelia.inspect_metamdbg_submission_reservations(crashed_outdir),
            )

            replacement_bind_outdir =
                joinpath(temporary_root, "replacement-owner-before-bind")
            replacement_bind_outputs = Mycelia._metamdbg_output_paths(
                replacement_bind_outdir,
                21,
            )
            replacement_bind_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            original_bind_reservation =
                Mycelia._metamdbg_submission_reservation(
                    replacement_bind_outputs,
                    replacement_bind_contract,
                    21;
                    owner_token = "original-bind-owner",
                    job_name = "custom-bind-job",
                )
            replacement_bind_reservation =
                Mycelia._metamdbg_submission_reservation(
                    replacement_bind_outputs,
                    replacement_bind_contract,
                    21;
                    owner_token = "replacement-bind-owner",
                    job_name = "custom-bind-job",
                )
            Mycelia._with_metamdbg_output_lock(replacement_bind_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    original_bind_reservation,
                    replacement_bind_outdir,
                )
            end
            stale_bind_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    replacement_bind_outdir,
                ),
            )
            Mycelia._with_metamdbg_output_lock(replacement_bind_outdir) do
                Mycelia._remove_metamdbg_submission_reservation!(
                    original_bind_reservation,
                )
                Mycelia._create_metamdbg_submission_reservation!(
                    replacement_bind_reservation,
                    replacement_bind_outdir,
                )
            end
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    stale_bind_metadata;
                    owner_token = stale_bind_metadata.owner_token,
                    job_id = "777",
                    confirm_submitted = true,
                ),
                ErrorException,
                r"owner_token changed before its scheduler job could be bound",
            )
            current_bind_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    replacement_bind_outdir,
                ),
            )
            Test.@test current_bind_metadata.owner_token ==
                       replacement_bind_reservation.owner_token
            Test.@test current_bind_metadata.job_id === nothing
            current_bound_metadata =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    current_bind_metadata;
                    owner_token = current_bind_metadata.owner_token,
                    job_id = "778",
                    confirm_submitted = true,
                )
            Test.@test current_bound_metadata.job_id == "778"
            Test.@test current_bound_metadata.scheduler_job_name ==
                       replacement_bind_reservation.scheduler_job_name
            Mycelia.reclaim_metamdbg_submission_reservation!(
                current_bound_metadata;
                owner_token = current_bound_metadata.owner_token,
                job_id = current_bound_metadata.job_id,
                confirm_cancelled = true,
            )

            post_submit_crash_outdir =
                joinpath(temporary_root, "post-submit-pre-bind-death")
            post_submit_crash_outputs =
                Mycelia._metamdbg_output_paths(post_submit_crash_outdir, 21)
            post_submit_crash_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            post_submit_crash_reservation =
                Mycelia._metamdbg_submission_reservation(
                    post_submit_crash_outputs,
                    post_submit_crash_contract,
                    21;
                    owner_token = "post-submit-pre-bind-owner",
                )
            Mycelia._with_metamdbg_output_lock(post_submit_crash_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    post_submit_crash_reservation,
                    post_submit_crash_outdir,
                )
            end
            post_submit_crash_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    post_submit_crash_outdir,
                ),
            )
            Test.@test post_submit_crash_metadata.submission_state == :reserved
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    post_submit_crash_metadata;
                    owner_token = post_submit_crash_metadata.owner_token,
                    job_id = "321",
                ),
                ArgumentError,
                r"confirm_submitted=true",
            )
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    post_submit_crash_metadata;
                    owner_token = "wrong-post-submit-owner",
                    job_id = "321",
                    confirm_submitted = true,
                ),
                ErrorException,
                r"owner token does not match",
            )
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    post_submit_crash_metadata;
                    owner_token = post_submit_crash_metadata.owner_token,
                    job_id = "invalid/job/id",
                    confirm_submitted = true,
                ),
                ArgumentError,
                r"must contain only decimal digits",
            )
            post_submit_bound =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    post_submit_crash_metadata;
                    owner_token = post_submit_crash_metadata.owner_token,
                    job_id = "321",
                    confirm_submitted = true,
                )
            Test.@test post_submit_bound.job_id == "321"
            Test.@test post_submit_bound.submission_state == :submitted
            post_submit_bound_reservation =
                Mycelia._metamdbg_bound_submission_reservation(
                    post_submit_crash_reservation,
                    post_submit_bound.job_id,
                )
            Test.@test readdir(post_submit_bound.path) ==
                       String["contract.json", "job.json"]
            Test.@test (
                stat(post_submit_bound_reservation.job_marker).mode & 0o777
            ) == 0o600
            _test_metamdbg_error(
                () -> Mycelia._bind_metamdbg_submission_job!(
                    post_submit_bound_reservation,
                    "322",
                ),
                ErrorException,
                r"refuses to overwrite an existing durable submission job record",
            )
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    post_submit_crash_metadata;
                    owner_token = post_submit_crash_metadata.owner_token,
                    job_id = "322",
                    confirm_submitted = true,
                ),
                ErrorException,
                r"already has a scheduler job id",
            )
            chmod(post_submit_bound_reservation.job_marker, 0o640)
            _test_metamdbg_error(
                () -> Mycelia.inspect_metamdbg_submission_reservations(
                    post_submit_crash_outdir,
                ),
                ErrorException,
                r"job record must have mode 0600",
            )
            chmod(post_submit_bound_reservation.job_marker, 0o600)
            canonical_job_contents = read(
                post_submit_bound_reservation.job_marker,
                String,
            )
            tampered_job_record = JSON.parse(canonical_job_contents)
            tampered_job_record["owner_token"] = "tampered-owner"
            write(
                post_submit_bound_reservation.job_marker,
                JSON.json(tampered_job_record) * "\n",
            )
            _test_metamdbg_error(
                () -> Mycelia.inspect_metamdbg_submission_reservations(
                    post_submit_crash_outdir,
                ),
                ErrorException,
                r"not canonical or does not match its workflow owner",
            )
            write(
                post_submit_bound_reservation.job_marker,
                canonical_job_contents,
            )
            inspected_post_submit_bound = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    post_submit_crash_outdir,
                ),
            )
            Test.@test inspected_post_submit_bound.job_id ==
                       "321"
            Test.@test inspected_post_submit_bound.submission_state ==
                       :submitted
            recovered_bound_submission =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    inspected_post_submit_bound;
                    owner_token = inspected_post_submit_bound.owner_token,
                    job_id = inspected_post_submit_bound.job_id,
                    confirm_cancelled = true,
                )
            Test.@test recovered_bound_submission.recovery_reason == :cancelled
            Test.@test !ispath(post_submit_bound.path)

            immediate_outdir = joinpath(
                temporary_root,
                "immediate-start-submitted-output",
            )
            immediate_outputs =
                Mycelia._metamdbg_output_paths(immediate_outdir, 21)
            immediate_fake_conda = joinpath(
                temporary_root,
                "immediate-start-fake-conda",
            )
            immediate_staged_observed = joinpath(
                temporary_root,
                "immediate-staged-input-observed.fastq",
            )
            immediate_original_reads = read(valid_reads, String)
            immediate_mutated_reads = replace(
                immediate_original_reads,
                "ACGT" => "TGCA";
                count = 1,
            )
            write(
                immediate_fake_conda,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "if [ \"\$1\" = \"run\" ]; then\n" *
                "  shift\n" *
                "  while [ \"\$#\" -gt 0 ] && [ \"\$1\" != \"python\" ] && [ \"\$1\" != \"metaMDBG\" ]; do shift; done\n" *
                "  [ \"\$#\" -gt 0 ] || exit 90\n" *
                "  tool=\$1; shift\n" *
                "  if [ \"\$tool\" = \"python\" ]; then exec python3 \"\$@\"; fi\n" *
                "  operation=\$1; shift\n" *
                "  output=; staged=; graph_k=\n" *
                "  while [ \"\$#\" -gt 0 ]; do\n" *
                "    case \"\$1\" in\n" *
                "      --out-dir|--assembly-dir) output=\$2; shift 2;;\n" *
                "      --in-hifi|--in-ont) staged=\$2; shift 2;;\n" *
                "      --k) graph_k=\$2; shift 2;;\n" *
                "      *) shift;;\n" *
                "    esac\n" *
                "  done\n" *
                "  if [ \"\$operation\" = \"asm\" ]; then\n" *
                "    [ -f \"\$staged\" ] && [ \"\$staged\" != $(Base.shell_escape(valid_reads)) ]\n" *
                "    printf '%s' $(Base.shell_escape(immediate_mutated_reads)) > $(Base.shell_escape(valid_reads))\n" *
                "    cp -- \"\$staged\" $(Base.shell_escape(immediate_staged_observed))\n" *
                "    printf '%s' $(Base.shell_escape(immediate_original_reads)) > $(Base.shell_escape(valid_reads))\n" *
                "    printf '>contig-1\\nACGT\\n' > \"\$output/contigs.fasta\"\n" *
                "  elif [ \"\$operation\" = \"gfa\" ]; then\n" *
                "    printf 'H\\tVN:Z:1.0\\nS\\tcontig-1\\tACGT\\n' > \"\$output/assemblyGraph_k\${graph_k}_4bps.gfa\"\n" *
                "  else exit 91; fi\n" *
                "  exit 0\n" *
                "fi\n" *
                "[ \"\$1\" = \"list\" ] || exit 90\n" *
                "printf '%s\\n' " *
                "'[{\"name\":\"metamdbg\",\"version\":\"1.4\"," *
                "\"build_string\":\"h43eeafb_2\"," *
                "\"channel\":\"bioconda\"}," *
                "{\"name\":\"libzlib\",\"version\":\"1.3.1\"," *
                "\"build_string\":\"h8359307_2\"," *
                "\"channel\":\"conda-forge\"}]'\n",
            )
            chmod(immediate_fake_conda, 0o700)
            immediate_script_path = joinpath(
                temporary_root,
                "immediate-start-executor.sh",
            )
            immediate_runner = function (
                    job::Mycelia.JobSpec,
                    executor::Mycelia.AbstractExecutor,
            )
                Test.@test executor isa Mycelia.SlurmExecutor
                Test.@test executor.hold
                Test.@test length(
                    Mycelia._metamdbg_submission_reservation_paths(
                        immediate_outdir,
                    ),
                ) == 1
                Test.@test occursin("staged_input_path_1", job.cmd)
                Test.@test occursin(
                    "metaMDBG asm --out-dir \"\$outdir\" --in-hifi " *
                    "\"\$staged_input_path_1\"",
                    job.cmd,
                )
                runtime_script = replace(
                    job.cmd,
                    r"(?m)^conda_runner=.*$" =>
                        "conda_runner=$(Base.shell_escape(immediate_fake_conda))",
                )
                write(immediate_script_path, runtime_script)
                Test.@test !Base.withenv(
                    "SLURM_JOB_ID" => "456",
                ) do
                    success(`bash $(immediate_script_path)`)
                end
                Test.@test length(
                    Mycelia._metamdbg_submission_reservation_paths(
                        immediate_outdir,
                    ),
                ) == 1
                Test.@test !ispath(immediate_outputs.contract_marker)
                return Mycelia.SubmitResult(
                    ok = true,
                    dry_run = false,
                    held = true,
                    scheduler_acceptance = :accepted,
                    site = :scg,
                    backend = :sbatch,
                    job_id = "456",
                )
            end
            immediate_release_runner = function (job_id::AbstractString)
                Test.@test String(job_id) == "456"
                inspected = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        immediate_outdir,
                    ),
                )
                Test.@test inspected.job_id == String(job_id)
                Test.@test inspected.submission_state == :submitted
                Test.@test Base.withenv(
                    "SLURM_JOB_ID" => String(job_id),
                ) do
                    success(`bash $(immediate_script_path)`)
                end
                return String(job_id)
            end
            immediate_result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir = immediate_outdir,
                executor = Mycelia.SlurmExecutor(dry_run = false),
                site = :scg,
                dependency_checker = _test_metamdbg_toolchain,
                local_runner = forbidden_runner,
                submission_runner = immediate_runner,
                submission_release_runner = immediate_release_runner,
            )
            Test.@test immediate_result.status == :submitted
            Test.@test immediate_result.submission.job_id ==
                       "456"
            Test.@test immediate_result.submission_reservation.job_id ==
                       "456"
            Test.@test immediate_result.submission_reservation.submission_state ==
                       :consumed
            Test.@test isfile(immediate_outputs.contract_marker)
            Test.@test isfile(immediate_outputs.completion_marker)
            Test.@test read(immediate_staged_observed, String) ==
                       immediate_original_reads
            Test.@test read(valid_reads, String) == immediate_original_reads
            Test.@test Mycelia._require_valid_metamdbg_fasta(
                immediate_outputs.contigs_gz,
                "immediate-start contigs",
            ) == immediate_outputs.contigs_gz
            Test.@test Mycelia._require_valid_metamdbg_gfa(
                joinpath(
                    immediate_outdir,
                    "assemblyGraph_k21_4bps.gfa",
                ),
                "immediate-start graph",
            ) == joinpath(
                immediate_outdir,
                "assemblyGraph_k21_4bps.gfa",
            )
            immediate_complete = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir = immediate_outdir,
                dependency_checker = () -> error(
                    "runtime-complete reuse must not provision metaMDBG",
                ),
                local_runner = forbidden_runner,
            )
            Test.@test immediate_complete.status == :complete
            Test.@test immediate_complete.provenance.package_inventory_sha256 ==
                       _test_metamdbg_toolchain().package_inventory_sha256
        end

        Test.@testset "submission ambiguity preserves recovery state" begin
            thrown_outdir = joinpath(
                temporary_root,
                "thrown-submission-failure-output",
            )
            submission_failure =
                ErrorException("synthetic submission failure")
            caught_failure = try
                Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = thrown_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = function (
                            _job::Mycelia.JobSpec,
                            _executor::Mycelia.AbstractExecutor,
                    )
                        Test.@test length(
                            Mycelia._metamdbg_submission_reservation_paths(
                                thrown_outdir,
                            ),
                        ) == 1
                        throw(submission_failure)
                    end,
                )
                nothing
            catch caught
                caught
            end
            Test.@test caught_failure isa ErrorException
            Test.@test occursin(
                "ambiguous response after attempting a held SLURM submission",
                sprint(showerror, caught_failure),
            )
            Test.@test occursin(
                "synthetic submission failure",
                sprint(showerror, caught_failure),
            )
            Test.@test length(
                Mycelia._metamdbg_submission_reservation_paths(thrown_outdir),
            ) == 1
            thrown_reservation = only(
                Mycelia.inspect_metamdbg_submission_reservations(thrown_outdir),
            )
            thrown_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    thrown_reservation;
                    owner_token = thrown_reservation.owner_token,
                    confirm_not_submitted = true,
                )
            Test.@test thrown_reclaimed.recovery_reason == :not_submitted

            legacy_failed_result_outdir = joinpath(
                temporary_root,
                "legacy-failed-result-submission-output",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = legacy_failed_result_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor,
                    ) ->
                        Mycelia.SubmitResult(
                            false,
                            false,
                            :scg,
                            :sbatch,
                            nothing,
                            nothing,
                            "sbatch --hold --parsable fixture.sbatch",
                            nothing,
                            "",
                            String[],
                            ["synthetic legacy sbatch failure"],
                        ),
                ),
                ErrorException,
                r"ambiguous response.*synthetic legacy sbatch failure",
            )
            legacy_failed_reservation = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    legacy_failed_result_outdir,
                ),
            )
            Mycelia.reclaim_metamdbg_submission_reservation!(
                legacy_failed_reservation;
                owner_token = legacy_failed_reservation.owner_token,
                confirm_not_submitted = true,
            )

            no_stdout_outdir = joinpath(
                temporary_root,
                "no-stdout-unknown-submission-output",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = no_stdout_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor,
                    ) -> Mycelia.SubmitResult(
                        ok = false,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :unknown,
                        site = :scg,
                        backend = :sbatch,
                        stdout = nothing,
                        errors = ["transport ended without stdout"],
                    ),
                ),
                ErrorException,
                r"ambiguous response.*transport ended without stdout",
            )
            no_stdout_reservation = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    no_stdout_outdir,
                ),
            )
            Mycelia.reclaim_metamdbg_submission_reservation!(
                no_stdout_reservation;
                owner_token = no_stdout_reservation.owner_token,
                confirm_not_submitted = true,
            )

            wrong_backend_outdir = joinpath(
                temporary_root,
                "wrong-backend-submission-output",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = wrong_backend_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor,
                    ) -> Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :not_attempted,
                        site = :scg,
                        backend = :salloc,
                    ),
                ),
                ErrorException,
                r"did not use the sbatch backend",
            )
            Test.@test isempty(
                Mycelia._metamdbg_submission_reservation_paths(
                    wrong_backend_outdir,
                ),
            )

            ambiguous_results = (
                legacy_unheld_success = (
                    result = Mycelia.SubmitResult(
                        true,
                        false,
                        :scg,
                        :sbatch,
                        nothing,
                        nothing,
                        "sbatch --hold --parsable fixture.sbatch",
                        "1001",
                        "1001\n",
                        String[],
                        String[],
                    ),
                    message = r"was not verified as held",
                ),
                missing_job_id = (
                    result = Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                    ),
                    message = r"returned no job id",
                ),
                unheld = (
                    result = Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = false,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = "999",
                    ),
                    message = r"was not verified as held",
                ),
                mismatched_stdout = (
                    result = Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = "12345",
                        stdout = "67890\n",
                    ),
                    message = r"stdout did not contain the exact returned job id",
                ),
            )
            for (label, fixture) in pairs(ambiguous_results)
                outdir = joinpath(temporary_root, "$(label)-submission-output")
                _test_metamdbg_error(
                    () -> Mycelia._run_metamdbg(;
                        hifi_reads = valid_reads,
                        outdir,
                        executor = Mycelia.SlurmExecutor(dry_run = false),
                        site = :scg,
                        dependency_checker = _test_metamdbg_toolchain,
                        local_runner = forbidden_runner,
                        submission_runner = (
                            _job::Mycelia.JobSpec,
                            _executor::Mycelia.AbstractExecutor,
                        ) -> fixture.result,
                    ),
                    ErrorException,
                    fixture.message,
                )
                inspected = only(
                    Mycelia.inspect_metamdbg_submission_reservations(outdir),
                )
                Test.@test inspected.job_id === nothing
                if fixture.result.job_id isa AbstractString
                    bound = Mycelia.bind_metamdbg_submission_reservation_job!(
                        inspected;
                        owner_token = inspected.owner_token,
                        job_id = fixture.result.job_id,
                        confirm_submitted = true,
                    )
                    Mycelia.reclaim_metamdbg_submission_reservation!(
                        bound;
                        owner_token = bound.owner_token,
                        job_id = bound.job_id,
                        confirm_cancelled = true,
                    )
                else
                    Mycelia.reclaim_metamdbg_submission_reservation!(
                        inspected;
                        owner_token = inspected.owner_token,
                        confirm_not_submitted = true,
                    )
                end
                Test.@test isempty(
                    Mycelia._metamdbg_submission_reservation_paths(outdir),
                )
            end

            binder_failure_outdir = joinpath(
                temporary_root,
                "binder-failure-submission-output",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = binder_failure_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor,
                    ) -> Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = "111",
                    ),
                    submission_job_binder = (
                        _reservation::NamedTuple,
                        _job_id::AbstractString,
                    ) -> error("synthetic durable bind failure"),
                    submission_release_runner = _job_id -> error(
                        "release must not run before durable binding",
                    ),
                ),
                ErrorException,
                r"could not durably bind job.json.*synthetic durable bind failure",
            )
            unbound_after_bind_failure = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    binder_failure_outdir,
                ),
            )
            Test.@test unbound_after_bind_failure.job_id === nothing
            recovered_bound =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    unbound_after_bind_failure;
                    owner_token = unbound_after_bind_failure.owner_token,
                    job_id = "111",
                    confirm_submitted = true,
                )
            Mycelia.reclaim_metamdbg_submission_reservation!(
                recovered_bound;
                owner_token = recovered_bound.owner_token,
                job_id = recovered_bound.job_id,
                confirm_cancelled = true,
            )

            release_failure_outdir = joinpath(
                temporary_root,
                "release-failure-submission-output",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = release_failure_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor,
                    ) -> Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = "222",
                    ),
                    submission_release_runner = _job_id -> error(
                        "synthetic scheduler release failure",
                    ),
                ),
                ErrorException,
                r"failed to release.*synthetic scheduler release failure",
            )
            bound_after_release_failure = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    release_failure_outdir,
                ),
            )
            Test.@test bound_after_release_failure.job_id ==
                       "222"
            Mycelia.reclaim_metamdbg_submission_reservation!(
                bound_after_release_failure;
                owner_token = bound_after_release_failure.owner_token,
                job_id = bound_after_release_failure.job_id,
                confirm_cancelled = true,
            )

            replacement_outdir = joinpath(
                temporary_root,
                "replacement-owner-submission-output",
            )
            replacement_outputs =
                Mycelia._metamdbg_output_paths(replacement_outdir, 21)
            replacement_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            replacement_reservation =
                Mycelia._metamdbg_submission_reservation(
                    replacement_outputs,
                    replacement_contract,
                    21;
                    owner_token = "replacement-owner-fixture",
                )
            replacement_failure =
                ErrorException("submission failed after owner replacement")
            caught_replacement_failure = try
                Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = replacement_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = function (
                            _job::Mycelia.JobSpec,
                            _executor::Mycelia.AbstractExecutor,
                    )
                        active_reservation = only(
                            Mycelia._metamdbg_submission_reservation_paths(
                                replacement_outdir,
                            ),
                        )
                        Test.@test active_reservation ==
                                   replacement_reservation.path
                        Mycelia._with_metamdbg_output_lock(
                            replacement_outdir,
                        ) do
                            active_record =
                                Mycelia._metamdbg_submission_reservation_from_path(
                                    active_reservation,
                                    replacement_outdir,
                                )
                            rm(active_record.output_root_reservation_marker)
                            rm(joinpath(
                                active_reservation,
                                Mycelia._METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
                            ))
                            rm(active_reservation)
                            Mycelia._create_metamdbg_submission_reservation!(
                                replacement_reservation,
                                replacement_outdir,
                            )
                        end
                        throw(replacement_failure)
                    end,
                )
                nothing
            catch caught
                caught
            end
            Test.@test caught_replacement_failure isa ErrorException
            Test.@test occursin(
                "ambiguous response after attempting a held SLURM submission",
                sprint(showerror, caught_replacement_failure),
            )
            Test.@test occursin(
                sprint(showerror, replacement_failure),
                sprint(showerror, caught_replacement_failure),
            )
            Test.@test isdir(replacement_reservation.path)
            Test.@test read(
                replacement_reservation.contract_marker,
                String,
            ) == replacement_reservation.contents
            Mycelia._with_metamdbg_output_lock(replacement_outdir) do
                Mycelia._remove_metamdbg_submission_reservation!(
                    replacement_reservation,
                )
            end
            Test.@test !ispath(replacement_reservation.path)
        end
    end
end
