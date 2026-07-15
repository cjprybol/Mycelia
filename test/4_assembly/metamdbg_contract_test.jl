# Default-CI contract tests for the single-technology metaMDBG wrapper.

import CodecZlib
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

function _test_metamdbg_toolchain()::NamedTuple
    return Mycelia._metamdbg_expected_toolchain()
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
        Test.@test provisioning_calls[] == 0
        Test.@test !ispath(joinpath(temporary_root, "missing-file"))

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
            local_runner = function (_command::Cmd)
                command_count[] += 1
                Test.@test isdir(output_lock_path)
                Test.@test !ispath(joinpath(
                    outdir,
                    "mycelia_metamdbg_contract.json",
                ))
                if command_count[] == 1
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
            Test.@test local_provisioning_calls[] == 1
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
            Test.@test local_provisioning_calls[] == 1
            Test.@test command_count[] == 2
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
            Test.@test !ispath(Mycelia._metamdbg_output_lock_path(outdir))
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
            outdir = joinpath(temporary_root, "exact-k")
            mkpath(outdir)
            _write_test_metamdbg_gzip_fasta!(
                joinpath(outdir, "contigs.fasta.gz"),
            )
            _write_test_metamdbg_gfa!(
                joinpath(outdir, "assemblyGraph_k20_4bps.gfa"),
            )
            _write_test_metamdbg_contract!(outdir, valid_reads)
            graph_calls = Ref(0)
            result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir,
                graph_k = 21,
                dependency_checker = _test_metamdbg_toolchain,
                local_runner = function (command::Cmd)
                    graph_calls[] += 1
                    Test.@test occursin(" gfa ", string(command))
                    _write_test_metamdbg_gfa!(
                        joinpath(outdir, "assemblyGraph_k21_4bps.gfa"),
                    )
                    return nothing
                end,
            )
            Test.@test graph_calls[] == 1
            Test.@test read(result.graph, String) ==
                       "H\tVN:Z:1.0\nS\tcontig-1\tACGT\n"

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
        end

        Test.@testset "spec-addressed exact environment and install lock" begin
            _, bundled_environment = Mycelia._metamdbg_paths()
            Test.@test Mycelia._metamdbg_sha256(bundled_environment) ==
                       Mycelia.METAMDBG_ENVIRONMENT_SPEC_SHA256
            Test.@test occursin(
                "metamdbg==1.4",
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
                    return "[{\"name\":\"metamdbg\",\"version\":\"1.4\"}]"
                end,
            )
            Test.@test inventory == Dict("metamdbg" => "1.4")
            Test.@test Mycelia._require_metamdbg_package_version(inventory) ==
                       _test_metamdbg_toolchain()
            _test_metamdbg_error(
                () -> Mycelia._require_metamdbg_package_version(
                    Dict("metamdbg" => "1.3"),
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
                environment_checker = () -> false,
                environment_creator = function (
                        specification::AbstractString,
                        environment_name::AbstractString;
                        force::Bool,
                )
                    Test.@test specification == environment_path
                    Test.@test environment_name == Mycelia.METAMDBG_ENV_NAME
                    Test.@test !force
                    created[] = true
                    return String(environment_name)
                end,
                package_inspector = () -> Dict("metamdbg" => "1.4"),
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
                job_id = "missing-parent-job",
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
            Test.@test occursin("object_count = split", script)
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
            Test.@test length(validation_positions) == 2
            Test.@test lock_position < first(validation_positions) <
                       reuse_position < last(validation_positions) <
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
                "if [ \"\$1\" != \"list\" ]; then exit 90; fi\n" *
                "printf '%s\\n' " *
                "'[{\"name\":\"metamdbg\",\"version\":\"1.4\"}]'\n",
            )
            chmod(fake_conda, 0o700)
            fake_asm =
                "printf '>contig-1\\nACGT\\n' > \"\$outdir/contigs.fasta\""
            fake_gfa =
                "printf 'H\\tVN:Z:1.0\\n" *
                "S\\tcontig-1\\tACGT\\n" *
                "S\\tcontig-2\\tTGCA\\n" *
                "S\\tcontig,with,comma\\tGATTACA\\tLN:i:7\\n" *
                "L\\tcontig-1\\t+\\tcontig-2\\t-\\t4M\\n" *
                "P\\tpath-1\\tcontig-1+,contig-2-\\t4M\\n" *
                "P\\tcomma-path\\tcontig,with,comma+\\t*\\n' " *
                "> \"\$outdir/assemblyGraph_k21_4bps.gfa\""

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
            Test.@test success(`bash $(retry_script_path)`)
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
            Test.@test !success(`bash $(timeout_script_path)`)
            Test.@test isdir(timeout_reservation.path)
            Test.@test isdir(timeout_lock)
            Test.@test !ispath(timeout_outputs.contract_marker)
            rm(timeout_lock)
            Mycelia._with_metamdbg_output_lock(timeout_outdir) do
                Mycelia._remove_metamdbg_submission_reservation!(
                    timeout_reservation,
                )
            end

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
            for (failure_type, invalid_asm) in pairs(invalid_executor_fastas)
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
                Test.@test !success(`bash $(invalid_script_path)`)
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
                leading_whitespace_pseudo_comment =
                    "printf ' #not-a-comment\\n" *
                    "S\\tcontig-1\\tACGT\\n' " *
                    "> \"\$outdir/assemblyGraph_k21_4bps.gfa\"",
            )
            for (failure_type, invalid_gfa) in pairs(invalid_executor_gfas)
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
                Test.@test !success(`bash $(invalid_script_path)`)
                Test.@test !ispath(invalid_outputs.contract_marker)
                Test.@test !ispath(invalid_reservation.path)
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(invalid_outdir),
                )
            end

            for queued_change in (:mutation, :replacement)
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
                Test.@test !success(`bash $(queued_script_path)`)
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
            Test.@test !success(`bash $(runtime_mutation_script_path)`)
            Test.@test !ispath(runtime_mutation_reservation.path)
            Test.@test !ispath(runtime_mutation_outputs.contract_marker)
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
            Test.@test success(`bash $(executable_script_path)`)
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

            sentinel_outdir = joinpath(temporary_root, "executor-sentinel")
            mkpath(sentinel_outdir)
            write(
                joinpath(sentinel_outdir, "contigs.fasta"),
                ">contig-1\nACGT\n",
            )
            _write_test_metamdbg_contract!(sentinel_outdir, valid_reads)
            sentinel_victim = joinpath(temporary_root, "sentinel-victim")
            write(sentinel_victim, "must-remain-unchanged\n")
            predictable_sentinel =
                joinpath(sentinel_outdir, "contigs.fasta.gz.tmp")
            symlink(sentinel_victim, predictable_sentinel)
            sentinel_outputs =
                Mycelia._metamdbg_output_paths(sentinel_outdir, 21)
            sentinel_reservation = Mycelia._metamdbg_submission_reservation(
                sentinel_outputs,
                input_contract,
                21;
                owner_token = "sentinel-executor-fixture",
            )
            Mycelia._with_metamdbg_output_lock(sentinel_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    sentinel_reservation,
                    sentinel_outdir,
                )
            end
            sentinel_script = Mycelia._metamdbg_executor_script(
                "echo 'unexpected assembly invocation' >&2; exit 91",
                fake_gfa,
                sentinel_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = sentinel_reservation,
            )
            sentinel_script_path =
                joinpath(temporary_root, "metamdbg-sentinel-fixture.sh")
            write(sentinel_script_path, sentinel_script)
            Test.@test success(`bash $(sentinel_script_path)`)
            Test.@test read(sentinel_victim, String) ==
                       "must-remain-unchanged\n"
            Test.@test islink(predictable_sentinel)
            Test.@test read(predictable_sentinel, String) ==
                       "must-remain-unchanged\n"
            Test.@test Mycelia._require_valid_metamdbg_fasta(
                sentinel_outputs.contigs_gz,
                "sentinel fixture contigs",
            ) == sentinel_outputs.contigs_gz
            Test.@test read(sentinel_outputs.contract_marker, String) ==
                       input_contract.contents
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
            submission_runner = function (
                    job::Mycelia.JobSpec,
                    _executor::Mycelia.AbstractExecutor,
            )
                submission_calls[] += 1
                reservations =
                    Mycelia._metamdbg_submission_reservation_paths(outdir)
                Test.@test length(reservations) == 1
                reserved_path[] = only(reservations)
                Test.@test isdir(reserved_path[])
                Test.@test occursin(
                    Base.shell_escape(reserved_path[]),
                    job.cmd,
                )
                return Mycelia.SubmitResult(
                    ok = true,
                    dry_run = false,
                    site = :scg,
                    backend = :sbatch,
                    job_id = "fixture-123",
                )
            end
            result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir,
                executor = Mycelia.SlurmExecutor(dry_run = false),
                site = :scg,
                dependency_checker = _test_metamdbg_toolchain,
                local_runner = forbidden_runner,
                submission_runner,
            )
            Test.@test result.status == :submitted
            Test.@test result.submission.job_id == "fixture-123"
            reservation_metadata = result.submission_reservation
            Test.@test reservation_metadata.path == reserved_path[]
            Test.@test reservation_metadata.job_id == "fixture-123"
            Test.@test !isempty(reservation_metadata.owner_token)
            Test.@test occursin(
                reservation_metadata.workflow_signature,
                reservation_metadata.path,
            )
            Test.@test submission_calls[] == 1
            Test.@test isdir(reserved_path[])

            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner,
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
                r"confirm_cancelled=true",
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
                    job_id = "wrong-job-id",
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
            Test.@test reclaimed.job_id == "fixture-123"
            Test.@test !ispath(reserved_path[])

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
            write(
                immediate_fake_conda,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "if [ \"\$1\" != \"list\" ]; then exit 90; fi\n" *
                "printf '%s\\n' " *
                "'[{\"name\":\"metamdbg\",\"version\":\"1.4\"}]'\n",
            )
            chmod(immediate_fake_conda, 0o700)
            immediate_script_path = joinpath(
                temporary_root,
                "immediate-start-executor.sh",
            )
            immediate_runner = function (
                    job::Mycelia.JobSpec,
                    _executor::Mycelia.AbstractExecutor,
            )
                Test.@test length(
                    Mycelia._metamdbg_submission_reservation_paths(
                        immediate_outdir,
                    ),
                ) == 1
                _write_test_metamdbg_gzip_fasta!(
                    immediate_outputs.contigs_gz,
                )
                _write_test_metamdbg_gfa!(joinpath(
                    immediate_outdir,
                    "assemblyGraph_k21_4bps.gfa",
                ))
                _write_test_metamdbg_contract!(
                    immediate_outdir,
                    valid_reads,
                )
                runtime_script = replace(
                    job.cmd,
                    r"(?m)^conda_runner=.*$" =>
                        "conda_runner=$(Base.shell_escape(immediate_fake_conda))",
                )
                write(immediate_script_path, runtime_script)
                Test.@test success(`bash $(immediate_script_path)`)
                Test.@test isempty(
                    Mycelia._metamdbg_submission_reservation_paths(
                        immediate_outdir,
                    ),
                )
                return Mycelia.SubmitResult(
                    ok = true,
                    dry_run = false,
                    site = :scg,
                    backend = :sbatch,
                    job_id = "fixture-immediate-456",
                )
            end
            immediate_result = Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                outdir = immediate_outdir,
                executor = Mycelia.SlurmExecutor(dry_run = false),
                site = :scg,
                dependency_checker = _test_metamdbg_toolchain,
                local_runner = forbidden_runner,
                submission_runner = immediate_runner,
            )
            Test.@test immediate_result.status == :submitted
            Test.@test immediate_result.submission.job_id ==
                       "fixture-immediate-456"
            Test.@test isfile(immediate_outputs.contract_marker)
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
        end

        Test.@testset "submission failures clean owned reservations" begin
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
            Test.@test caught_failure === submission_failure
            Test.@test isempty(
                Mycelia._metamdbg_submission_reservation_paths(thrown_outdir),
            )

            failed_result_outdir = joinpath(
                temporary_root,
                "failed-result-submission-output",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = failed_result_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor,
                    ) ->
                        Mycelia.SubmitResult(
                            ok = false,
                            dry_run = false,
                            site = :scg,
                            backend = :sbatch,
                            errors = ["synthetic sbatch rejection"],
                        ),
                ),
                ErrorException,
                r"submission failed.*synthetic sbatch rejection",
            )
            Test.@test isempty(
                Mycelia._metamdbg_submission_reservation_paths(
                    failed_result_outdir,
                ),
            )

            falsified_results = (
                wrong_backend = (
                    result = Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        site = :scg,
                        backend = :salloc,
                    ),
                    message = r"did not use the sbatch backend",
                ),
                missing_job_id = (
                    result = Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        site = :scg,
                        backend = :sbatch,
                    ),
                    message = r"returned no job id",
                ),
            )
            for (label, fixture) in pairs(falsified_results)
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
                Test.@test isempty(
                    Mycelia._metamdbg_submission_reservation_paths(outdir),
                )
            end

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
            Test.@test caught_replacement_failure === replacement_failure
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
