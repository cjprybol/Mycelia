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
            Test.@test expected_contract.contract.schema_version == 3
            Test.@test expected_input.sha256 ==
                       Mycelia._metamdbg_sha256(valid_reads)
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
            Test.@test changed_content_input.size_bytes == expected_input.size_bytes
            Test.@test changed_content_input.modification_time_ns ==
                       expected_input.modification_time_ns
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
            lock_position = first(something(findfirst(
                "lock_acquired=1",
                script,
            )))
            digest_position = first(something(findfirst(
                "actual_input_sha256_1=",
                script,
            )))
            reuse_position = first(something(findfirst(
                "contract_exists=0",
                script,
            )))
            Test.@test lock_position < digest_position < reuse_position
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
                "printf 'H\\tVN:Z:1.0\\nS\\tcontig-1\\tACGT\\n' " *
                "> \"\$outdir/assemblyGraph_k21_4bps.gfa\""

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
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(queued_outdir),
                )
            end

            executable_outdir = joinpath(temporary_root, "executor-fixture")
            executable_outputs =
                Mycelia._metamdbg_output_paths(executable_outdir, 21)
            executable_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                executable_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
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
            sentinel_script = Mycelia._metamdbg_executor_script(
                "echo 'unexpected assembly invocation' >&2; exit 91",
                fake_gfa,
                sentinel_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
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
        end
    end
end
