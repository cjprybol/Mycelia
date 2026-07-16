# Default-CI lifecycle and provenance tests for the common Unicycler wrapper.
# All external effects are injected; no Conda environment or assembler runs.

import Mycelia
import Test

function _test_unicycler_inventory(;
        unicycler_build::AbstractString = "pyhdfd78af_0",
)::Vector{NamedTuple}
    return NamedTuple[
        (
            name = "python",
            version = "3.11.9",
            build = "h955ad1f_0_cpython",
            channel = "conda-forge",
        ),
        (
            name = "spades",
            version = "4.2.0",
            build = "h5ca1c30_1",
            channel = "bioconda",
        ),
        (
            name = "unicycler",
            version = "0.5.1",
            build = String(unicycler_build),
            channel = "bioconda",
        ),
    ]
end

function _test_unicycler_toolchain(;
        unicycler_build::AbstractString = "pyhdfd78af_0",
)::Dict{String, Any}
    return Mycelia._unicycler_toolchain_metadata(
        _test_unicycler_inventory(; unicycler_build),
    )
end

function _test_unicycler_error(
        thunk::Function,
        expected_message::AbstractString,
        expected_type::Type{<:Exception} = ErrorException,
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

Test.@testset "common Unicycler wrapper contract" begin
    Test.@testset "inventory normalization and detached provenance" begin
        inventory = _test_unicycler_inventory()
        digest = Mycelia._unicycler_package_inventory_sha256(inventory)
        Test.@test occursin(r"^[0-9a-f]{64}$", digest)
        Test.@test Mycelia._unicycler_package_inventory_sha256(
            reverse(inventory),
        ) == digest

        changed_inventory = _test_unicycler_inventory(
            unicycler_build = "pyhdfd78af_1",
        )
        Test.@test Mycelia._unicycler_package_inventory_sha256(
            changed_inventory,
        ) != digest
        Test.@test Mycelia._unicycler_environment_lock_path(
            "/opt/conda-a/bin/conda",
        ) != Mycelia._unicycler_environment_lock_path(
            "/opt/conda-b/bin/conda",
        )

        reported = _test_unicycler_toolchain()
        Test.@test reported["package_inventory_sha256"] ==
                   Mycelia._unicycler_package_inventory_sha256(inventory)
        python_packages = filter(
            package -> package["name"] == "python",
            reported["packages"],
        )
        Test.@test length(python_packages) == 1
        Test.@test only(python_packages) == Dict{String, Any}(
            "name" => "python",
            "version" => "3.11.9",
            "build" => "h955ad1f_0_cpython",
            "channel" => "conda-forge",
        )
        validated = Mycelia._require_unicycler_toolchain_provenance(reported)
        reported_packages = reported["packages"]
        reported_packages[1]["version"] = "mutated-after-validation"
        Test.@test validated["packages"][1]["version"] !=
                   "mutated-after-validation"
        Test.@test validated == _test_unicycler_toolchain()

        tampered = _test_unicycler_toolchain()
        tampered["packages"][1]["version"] = "0"
        _test_unicycler_error(
            () -> Mycelia._require_unicycler_toolchain_provenance(tampered),
            "inventory digest does not match",
        )
        for drift_field in ("build", "channel")
            transitive_drift = _test_unicycler_toolchain()
            python_package = only(filter(
                package -> package["name"] == "python",
                transitive_drift["packages"],
            ))
            python_package[drift_field] *= "-transitive-drift"
            _test_unicycler_error(
                () -> Mycelia._require_unicycler_toolchain_provenance(
                    transitive_drift,
                ),
                "inventory digest does not match",
            )
        end
    end

    Test.@testset "locked fresh local execution reports exact provenance" begin
        mktempdir() do temp_dir
            output_dir = joinpath(temp_dir, "output with spaces")
            lock_path = joinpath(temp_dir, "unicycler-environment.pid")
            output_lock_path = Mycelia._unicycler_output_lock_path(output_dir)
            conda_runner = joinpath(temp_dir, "custom conda", "bin", "conda")
            events = NamedTuple[]
            observe = stage -> push!(events, (;
                stage,
                environment_held = isfile(lock_path),
                output_held = isfile(output_lock_path),
            ))
            result = withenv("MYCELIA_CONDA_RUNNER" => conda_runner) do
                Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/short R1.fastq",
                    short_2 = "/reads/short R2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = output_dir,
                    threads = 3,
                    executor = Mycelia.LocalExecutor(),
                    environment_preparer = runner -> begin
                        Test.@test runner == abspath(conda_runner)
                        observe(:prepare)
                    end,
                    toolchain_inspector = runner -> begin
                        Test.@test runner == abspath(conda_runner)
                        observe(:snapshot)
                        _test_unicycler_toolchain()
                    end,
                    environment_lock_path = lock_path,
                    command_runner = command -> begin
                        observe(:runner)
                        Test.@test first(command.exec) == abspath(conda_runner)
                        Test.@test "/reads/short R1.fastq" in command.exec
                        Test.@test "/reads/short R2.fastq" in command.exec
                        mkpath(output_dir)
                        write(
                            joinpath(output_dir, "assembly.fasta"),
                            ">contig\nACGT\n",
                        )
                        write(
                            joinpath(output_dir, "assembly.gfa"),
                            "H\tVN:Z:1.0\nS\tcontig\tACGT\n",
                        )
                    end,
                )
            end
            Test.@test [event.stage for event in events] ==
                       [:prepare, :snapshot, :runner, :snapshot]
            Test.@test all(event -> event.environment_held, events)
            Test.@test all(event -> event.output_held, events)
            Test.@test !ispath(lock_path)
            Test.@test !ispath(output_lock_path)
            Test.@test result.provenance_status == "realized-local-exact"
            Test.@test result.toolchain == _test_unicycler_toolchain()
            Test.@test result.outdir == abspath(output_dir)
            Test.@test result.requested_outdir == output_dir
            Test.@test result.assembly == joinpath(
                abspath(output_dir),
                "assembly.fasta",
            )
            Test.@test result.graph == joinpath(
                abspath(output_dir),
                "assembly.gfa",
            )
        end
    end

    Test.@testset "failure releases lock and inventory mutation fails loud" begin
        mktempdir() do temp_dir
            lock_path = joinpath(temp_dir, "unicycler-environment.pid")
            output_dir = joinpath(temp_dir, "failure-output")
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = output_dir,
                    environment_preparer = runner -> nothing,
                    toolchain_inspector = runner -> _test_unicycler_toolchain(),
                    environment_lock_path = lock_path,
                    command_runner = command -> error("synthetic runner failure"),
                ),
                "synthetic runner failure",
            )
            Test.@test !ispath(lock_path)
            reacquired = Mycelia._with_unicycler_environment_lock(lock_path) do
                Test.@test isfile(lock_path)
                :reacquired
            end
            Test.@test reacquired == :reacquired
            Test.@test !ispath(lock_path)

            snapshots = Any[
                _test_unicycler_toolchain(),
                _test_unicycler_toolchain(
                    unicycler_build = "pyhdfd78af_1",
                ),
            ]
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = output_dir,
                    environment_preparer = runner -> nothing,
                    toolchain_inspector = runner -> popfirst!(snapshots),
                    environment_lock_path = lock_path,
                    command_runner = command -> begin
                        mkpath(output_dir)
                        write(
                            joinpath(output_dir, "assembly.fasta"),
                            ">contig\nACGT\n",
                        )
                        write(
                            joinpath(output_dir, "assembly.gfa"),
                            "H\tVN:Z:1.0\nS\tcontig\tACGT\n",
                        )
                    end,
                ),
                "package inventory changed while the assembler ran",
            )
            Test.@test isempty(snapshots)
            Test.@test !ispath(lock_path)
        end
    end

    Test.@testset "reuse is unclaimed and nonlocal executors are rejected" begin
        mktempdir() do temp_dir
            lock_path = joinpath(temp_dir, "unicycler-environment.pid")
            reused_dir = joinpath(temp_dir, "reused")
            mkpath(reused_dir)
            write(joinpath(reused_dir, "assembly.fasta"), ">old\nACGT\n")
            write(
                joinpath(reused_dir, "assembly.gfa"),
                "H\tVN:Z:1.0\nS\told\tACGT\n",
            )
            reused = Mycelia._run_unicycler_with_contract(
                short_1 = "/reads/r1.fastq",
                short_2 = "/reads/r2.fastq",
                long_reads = "/reads/long.fastq",
                outdir = reused_dir,
                environment_preparer = runner -> nothing,
                toolchain_inspector = runner -> _test_unicycler_toolchain(),
                environment_lock_path = lock_path,
                command_runner = command -> error("reuse unexpectedly executed"),
            )
            Test.@test reused.toolchain === nothing
            Test.@test reused.provenance_status == "unavailable-reused-output"
            Test.@test reused.outdir == realpath(reused_dir)
            Test.@test reused.requested_outdir == reused_dir
            Test.@test reused.assembly == joinpath(
                realpath(reused_dir),
                "assembly.fasta",
            )
            Test.@test reused.graph == joinpath(
                realpath(reused_dir),
                "assembly.gfa",
            )

            incomplete_dir = joinpath(temp_dir, "incomplete-reuse")
            mkpath(incomplete_dir)
            incomplete_assembly = joinpath(incomplete_dir, "assembly.fasta")
            write(incomplete_assembly, ">old\nACGT\n")
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = incomplete_dir,
                    environment_preparer = runner -> nothing,
                    toolchain_inspector = runner -> _test_unicycler_toolchain(),
                    environment_lock_path = lock_path,
                ),
                "non-empty assembly GFA",
            )
            Test.@test read(incomplete_assembly, String) == ">old\nACGT\n"

            empty_dir = joinpath(temp_dir, "empty-reuse")
            mkpath(empty_dir)
            touch(joinpath(empty_dir, "assembly.fasta"))
            write(
                joinpath(empty_dir, "assembly.gfa"),
                "H\tVN:Z:1.0\nS\told\tACGT\n",
            )
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = empty_dir,
                    environment_preparer = runner -> nothing,
                    toolchain_inspector = runner -> _test_unicycler_toolchain(),
                    environment_lock_path = lock_path,
                ),
                "non-empty assembly FASTA",
            )

            symlink_dir = joinpath(temp_dir, "symlink-reuse")
            mkpath(symlink_dir)
            external_assembly = joinpath(temp_dir, "external.fasta")
            write(external_assembly, ">old\nACGT\n")
            symlink(external_assembly, joinpath(symlink_dir, "assembly.fasta"))
            write(
                joinpath(symlink_dir, "assembly.gfa"),
                "H\tVN:Z:1.0\nS\told\tACGT\n",
            )
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = symlink_dir,
                    environment_preparer = runner -> nothing,
                    toolchain_inspector = runner -> _test_unicycler_toolchain(),
                    environment_lock_path = lock_path,
                ),
                "regular non-symlink file",
            )

            collector = Mycelia.CollectExecutor()
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = joinpath(temp_dir, "submitted"),
                    executor = collector,
                    environment_preparer = runner ->
                        error("unexpected preparation"),
                    toolchain_inspector = runner ->
                        error("unexpected inspection"),
                    environment_lock_path = lock_path,
                ),
                "supports only synchronous local execution",
                ArgumentError,
            )
            Test.@test isempty(collector.jobs)
            Test.@test !ispath(lock_path)
        end
    end

    Test.@testset "public relative output returns canonical safety paths" begin
        mktempdir() do temp_dir
            conda_root = joinpath(temp_dir, "fake-conda")
            conda_runner = joinpath(conda_root, "bin", "conda")
            environment_prefix = joinpath(
                conda_root,
                "envs",
                "unicycler",
            )
            mkpath(dirname(conda_runner))
            mkpath(joinpath(environment_prefix, "conda-meta"))
            inventory = Mycelia.JSON.json(_test_unicycler_inventory())
            script = replace(
                raw"""#!/bin/sh
if [ "$1" = "env" ] && [ "$2" = "list" ]; then
    printf '%s\n' '__ENVIRONMENT_PREFIX__'
    exit 0
fi
if [ "$1" = "list" ] && [ "$2" = "-n" ] && [ "$3" = "unicycler" ]; then
    cat <<'MYCELIA_UNICYCLER_JSON'
__PACKAGE_INVENTORY__
MYCELIA_UNICYCLER_JSON
    exit 0
fi
exit 66
""",
                "__ENVIRONMENT_PREFIX__" => environment_prefix,
                "__PACKAGE_INVENTORY__" => inventory,
            )
            write(conda_runner, script)
            chmod(conda_runner, 0o755)

            cd(temp_dir) do
                requested_outdir = "relative-unicycler-output"
                mkpath(requested_outdir)
                write(
                    joinpath(requested_outdir, "assembly.fasta"),
                    ">contig\nACGT\n",
                )
                write(
                    joinpath(requested_outdir, "assembly.gfa"),
                    "H\tVN:Z:1.0\nS\tcontig\tACGT\n",
                )
                result = withenv(
                    "MYCELIA_CONDA_RUNNER" => conda_runner,
                ) do
                    Mycelia.run_unicycler(
                        short_1 = "reads/r1.fastq",
                        short_2 = "reads/r2.fastq",
                        long_reads = "reads/long.fastq",
                        outdir = requested_outdir,
                    )
                end
                canonical_outdir = realpath(requested_outdir)
                Test.@test result.outdir == canonical_outdir
                Test.@test result.assembly == joinpath(
                    canonical_outdir,
                    "assembly.fasta",
                )
                Test.@test result.graph == joinpath(
                    canonical_outdir,
                    "assembly.gfa",
                )
                Test.@test result.requested_outdir == requested_outdir
                Test.@test result.toolchain === nothing
                Test.@test result.provenance_status ==
                           "unavailable-reused-output"
            end
        end
    end

    Test.@testset "output ownership fails closed" begin
        _test_unicycler_error(
            () -> Mycelia._run_unicycler_with_contract(
                short_1 = "/reads/r1.fastq",
                short_2 = "/reads/r2.fastq",
                long_reads = "/reads/long.fastq",
                outdir = "   ",
                environment_preparer = runner ->
                    error("unexpected preparation"),
            ),
            "outdir must be a non-empty path",
            ArgumentError,
        )

        mktempdir() do temp_dir
            output_dir = joinpath(temp_dir, "incomplete")
            mkpath(output_dir)
            marker = joinpath(output_dir, "caller-owned.txt")
            write(marker, "preserve me")
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = output_dir,
                    environment_preparer = runner -> nothing,
                    toolchain_inspector = runner -> _test_unicycler_toolchain(),
                    environment_lock_path = joinpath(temp_dir, "environment.pid"),
                    command_runner = command -> error("unexpected execution"),
                ),
                "Refusing to remove non-empty Unicycler outdir",
                ArgumentError,
            )
            Test.@test read(marker, String) == "preserve me"
        end

        mktempdir() do temp_dir
            physical_parent = joinpath(temp_dir, "physical")
            alias_parent = joinpath(temp_dir, "alias")
            mkpath(physical_parent)
            symlink(physical_parent, alias_parent)
            physical_outdir = joinpath(physical_parent, "shared-output")
            alias_outdir = joinpath(alias_parent, "shared-output")
            expected_lock_path = Mycelia._unicycler_output_lock_path(
                physical_outdir,
            )
            Test.@test Mycelia._unicycler_output_lock_path(alias_outdir) ==
                       expected_lock_path

            observed_lock_path = Ref("")
            pidlock_runner = function (
                    action::Function,
                    lock_path::AbstractString;
                    stale_age::Real,
                    refresh::Real,
            )
                observed_lock_path[] = String(lock_path)
                Test.@test stale_age > 0
                Test.@test refresh > 0
                return action()
            end
            Test.@test Mycelia._with_unicycler_output_lock(
                alias_outdir;
                pidlock_runner,
            ) do reserved_outdir
                Test.@test reserved_outdir == physical_outdir
                :reserved
            end == :reserved
            Test.@test observed_lock_path[] == expected_lock_path

            alias_result = Mycelia._run_unicycler_with_contract(
                short_1 = "/reads/r1.fastq",
                short_2 = "/reads/r2.fastq",
                long_reads = "/reads/long.fastq",
                outdir = alias_outdir,
                environment_preparer = runner -> nothing,
                toolchain_inspector = runner -> _test_unicycler_toolchain(),
                environment_lock_path = joinpath(
                    temp_dir,
                    "alias-environment.pid",
                ),
                command_runner = command -> begin
                    output_flag = findfirst(==("-o"), command.exec)
                    Test.@test output_flag !== nothing
                    command_outdir = command.exec[something(output_flag) + 1]
                    Test.@test command_outdir == physical_outdir
                    mkpath(command_outdir)
                    write(
                        joinpath(command_outdir, "assembly.fasta"),
                        ">contig\nACGT\n",
                    )
                    write(
                        joinpath(command_outdir, "assembly.gfa"),
                        "H\tVN:Z:1.0\nS\tcontig\tACGT\n",
                    )
                end,
            )
            Test.@test alias_result.outdir == physical_outdir
            Test.@test alias_result.requested_outdir == alias_outdir
            Test.@test alias_result.assembly == joinpath(
                physical_outdir,
                "assembly.fasta",
            )
            Test.@test alias_result.graph == joinpath(
                physical_outdir,
                "assembly.gfa",
            )

            retarget_physical_parent = joinpath(temp_dir, "retarget-physical")
            retarget_alias_parent = joinpath(temp_dir, "retarget-alias")
            retarget_destination = joinpath(temp_dir, "retarget-destination")
            mkpath(retarget_physical_parent)
            mkpath(retarget_destination)
            symlink(retarget_physical_parent, retarget_alias_parent)
            retarget_alias_outdir = joinpath(
                retarget_alias_parent,
                "output",
            )
            retarget_physical_outdir = joinpath(
                retarget_physical_parent,
                "output",
            )
            inspection_count = Ref(0)
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = retarget_alias_outdir,
                    environment_preparer = runner -> nothing,
                    toolchain_inspector = runner -> begin
                        inspection_count[] += 1
                        if inspection_count[] == 2
                            rm(retarget_alias_parent)
                            symlink(
                                retarget_destination,
                                retarget_alias_parent,
                            )
                        end
                        return _test_unicycler_toolchain()
                    end,
                    environment_lock_path = joinpath(
                        temp_dir,
                        "retarget-environment.pid",
                    ),
                    command_runner = command -> begin
                        mkpath(retarget_physical_outdir)
                        write(
                            joinpath(
                                retarget_physical_outdir,
                                "assembly.fasta",
                            ),
                            ">contig\nACGT\n",
                        )
                        write(
                            joinpath(
                                retarget_physical_outdir,
                                "assembly.gfa",
                            ),
                            "H\tVN:Z:1.0\nS\tcontig\tACGT\n",
                        )
                    end,
                ),
                "outdir changed physical identity after reservation",
                ArgumentError,
            )
            Test.@test inspection_count[] == 2
            Test.@test isfile(joinpath(
                retarget_physical_outdir,
                "assembly.fasta",
            ))

            external_outdir = joinpath(temp_dir, "external")
            mkpath(external_outdir)
            write(
                joinpath(external_outdir, "assembly.fasta"),
                ">external\nACGT\n",
            )
            write(
                joinpath(external_outdir, "assembly.gfa"),
                "H\tVN:Z:1.0\nS\texternal\tACGT\n",
            )
            linked_outdir = joinpath(temp_dir, "linked-output")
            symlink(external_outdir, linked_outdir)
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = linked_outdir,
                    environment_preparer = runner ->
                        error("environment must not be prepared"),
                ),
                "outdir must not be a symbolic link",
                ArgumentError,
            )

            retargeted_outdir = joinpath(temp_dir, "retargeted-output")
            retargeting_lock_runner = function (
                    action::Function,
                    reserved_outdir::AbstractString,
            )
                symlink(external_outdir, retargeted_outdir)
                return action(reserved_outdir)
            end
            _test_unicycler_error(
                () -> Mycelia._run_unicycler_with_contract(
                    short_1 = "/reads/r1.fastq",
                    short_2 = "/reads/r2.fastq",
                    long_reads = "/reads/long.fastq",
                    outdir = retargeted_outdir,
                    environment_preparer = runner ->
                        error("environment must not be prepared"),
                    output_lock_runner = retargeting_lock_runner,
                ),
                "outdir must not be a symbolic link",
                ArgumentError,
            )
            Test.@test read(
                joinpath(external_outdir, "assembly.fasta"),
                String,
            ) == ">external\nACGT\n"
        end
    end
end
