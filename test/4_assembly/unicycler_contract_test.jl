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
            conda_runner = joinpath(temp_dir, "custom conda", "bin", "conda")
            events = NamedTuple[]
            observe = stage -> push!(events, (; stage, held = isfile(lock_path)))
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
            Test.@test all(event -> event.held, events)
            Test.@test !ispath(lock_path)
            Test.@test result.provenance_status == "realized-local-exact"
            Test.@test result.toolchain == _test_unicycler_toolchain()
            Test.@test result.outdir == abspath(output_dir)
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
    end
end
