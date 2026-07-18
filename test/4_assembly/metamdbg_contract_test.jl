# Default-CI contract tests for the single-technology metaMDBG wrapper.

import CodecZlib
import FileWatching
import JSON
import Logging
import Test
import Mycelia

struct _TestUnverifiableMetamdbgExecutor <: Mycelia.AbstractExecutor end

struct _TestInterruptingMetamdbgString <: AbstractString
    interrupt::InterruptException
end

Base.ncodeunits(::_TestInterruptingMetamdbgString)::Int = 1

function Base.codeunit(
        value::_TestInterruptingMetamdbgString,
        index::Integer,
)::UInt8
    index == 1 || throw(BoundsError(value, index))
    return Base.codeunit("x", 1)
end

function Base.isvalid(
        ::_TestInterruptingMetamdbgString,
        index::Integer,
)::Bool
    return index == 1
end

function Base.String(value::_TestInterruptingMetamdbgString)::String
    throw(value.interrupt)
end

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

function _test_metamdbg_output_root_pid_lock_identity(
        lock_path::AbstractString,
)::NamedTuple
    return Mycelia._metamdbg_output_root_pid_lock_snapshot(lock_path).identity
end

function _test_metamdbg_file_publisher_hardening(
        publisher::Function,
        marker::AbstractString,
        held_path::AbstractString,
        cleanup_pattern::Regex,
        label::AbstractString,
)::Nothing
    captured_temporary = Ref("")
    replacement_contents = "$(label) replacement must survive\n"
    observed_aba = Test.@test_logs (
        :warn,
        cleanup_pattern,
    ) min_level=Logging.Warn match_mode=:any begin
        try
            publisher(
                (temporary_path, _identity) -> begin
                    captured_temporary[] = temporary_path
                    mv(temporary_path, held_path)
                    write(temporary_path, replacement_contents)
                end,
                (_path, _identity) -> nothing,
            )
            nothing
        catch caught
            caught
        end
    end
    Test.@test observed_aba isa ErrorException
    Test.@test occursin(
        "temporary artifact was replaced",
        sprint(showerror, observed_aba),
    )
    Test.@test read(captured_temporary[], String) == replacement_contents
    Test.@test isfile(held_path)
    Test.@test !ispath(marker)
    rm(captured_temporary[]; force = true)
    rm(held_path; force = true)

    post_hardlink_primary = ErrorException(
        "synthetic $(label) post-hardlink failure",
    )
    observed_post_hardlink = try
        publisher(
            (temporary_path, _identity) ->
                (captured_temporary[] = temporary_path),
            (_path, _identity) -> throw(post_hardlink_primary),
        )
        nothing
    catch caught
        caught
    end
    Test.@test observed_post_hardlink === post_hardlink_primary
    Test.@test !ispath(captured_temporary[])
    Test.@test !ispath(marker)
    return nothing
end

Test.@testset "metaMDBG interrupts retain control-flow identity" begin
    mktempdir() do temporary_root
        child_cleanup_root = joinpath(temporary_root, "child-cleanup")
        child_directory = joinpath(child_cleanup_root, "nested")
        mkpath(child_directory)
        write(joinpath(child_directory, "sentinel"), "retain\n")
        cleanup_directory = Base.Filesystem.open(
            child_cleanup_root,
            Base.JL_O_RDONLY |
            Base.JL_O_DIRECTORY |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
        )
        child_open_interrupt = InterruptException()
        observed_child_open_interrupt = try
            Mycelia._run_cleanup_after_primary_error!(
                () -> Mycelia._remove_metamdbg_quarantined_directory_contents!(
                    cleanup_directory,
                    child_cleanup_root;
                    child_opener = (
                        _directory::Base.Filesystem.File,
                        _entry::AbstractString,
                        _flags::Integer,
                    ) -> throw(child_open_interrupt),
                ),
                ErrorException("synthetic primary failure"),
                "synthetic child cleanup failure",
            )
            nothing
        catch caught
            caught
        end
        close(cleanup_directory)
        Test.@test observed_child_open_interrupt === child_open_interrupt
        Test.@test isfile(joinpath(child_directory, "sentinel"))

        reads = joinpath(temporary_root, "reads.fastq")
        write(reads, "@read-1\nACGT\n+\nIIII\n")
        interrupt_toolchain = () ->
            Mycelia._require_metamdbg_package_version(NamedTuple[
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
        reservation_outdir = joinpath(temporary_root, "reservation-parser")
        reservation_outputs =
            Mycelia._metamdbg_output_paths(reservation_outdir, 21)
        reservation_contract = Mycelia._metamdbg_input_contract(
            Mycelia._metamdbg_selected_input(reads, nothing),
            3,
        )
        reservation = Mycelia._metamdbg_submission_reservation(
            reservation_outputs,
            reservation_contract,
            21;
            owner_token = "interrupt-parser-fixture",
        )
        Mycelia._with_metamdbg_output_lock(reservation_outdir) do
            Mycelia._create_metamdbg_submission_reservation!(
                reservation,
                reservation_outdir,
            )
        end

        contract_parse_interrupt = InterruptException()
        observed_contract_parse_interrupt = try
            Mycelia._metamdbg_submission_reservation_from_path(
                reservation.path,
                reservation_outdir;
                json_parser = _contents -> throw(contract_parse_interrupt),
            )
            nothing
        catch caught
            caught
        end
        Test.@test observed_contract_parse_interrupt === contract_parse_interrupt

        write(reservation.job_marker, "{}\n")
        job_parse_interrupt = InterruptException()
        parse_calls = Ref(0)
        observed_job_parse_interrupt = try
            Mycelia._metamdbg_submission_reservation_from_path(
                reservation.path,
                reservation_outdir;
                json_parser = contents -> begin
                    parse_calls[] += 1
                    parse_calls[] == 1 && return JSON.parse(contents)
                    throw(job_parse_interrupt)
                end,
            )
            nothing
        catch caught
            caught
        end
        Test.@test parse_calls[] == 2
        Test.@test observed_job_parse_interrupt === job_parse_interrupt

        completion_marker = joinpath(temporary_root, "completion.json")
        write(completion_marker, "{}\n")
        completion_parse_interrupt = InterruptException()
        observed_completion_parse_interrupt = try
            Mycelia._require_metamdbg_completion_manifest!(
                (; completion_marker),
                (;),
                (;),
                21;
                json_parser = _contents ->
                    throw(completion_parse_interrupt),
            )
            nothing
        catch caught
            caught
        end
        Test.@test observed_completion_parse_interrupt ===
                   completion_parse_interrupt

        release_interrupt = InterruptException()
        observed_release_interrupt = try
            Mycelia._release_metamdbg_submission_job!(
                (;
                    path = joinpath(temporary_root, "release-reservation"),
                    scheduler_job_name = "interrupt-release",
                ),
                "919",
                nothing,
                _job_id -> throw(release_interrupt),
            )
            nothing
        catch caught
            caught
        end
        Test.@test observed_release_interrupt === release_interrupt

        binder_interrupt = InterruptException()
        binder_outdir = joinpath(temporary_root, "binder-interrupt")
        observed_binder_interrupt = try
            Mycelia._run_metamdbg(;
                hifi_reads = reads,
                outdir = binder_outdir,
                executor = Mycelia.SlurmExecutor(dry_run = false),
                site = :scg,
                dependency_checker = interrupt_toolchain,
                local_runner = _command -> error("local runner is forbidden"),
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
                    job_id = "920",
                    stdout = "920\n",
                ),
                submission_job_binder = (
                    _reservation::NamedTuple,
                    _job_id::AbstractString,
                ) -> throw(binder_interrupt),
                submission_release_runner = _job_id -> error(
                    "release must not run after an interrupted bind",
                ),
            )
            nothing
        catch caught
            caught
        end
        Test.@test observed_binder_interrupt === binder_interrupt

        submission_interrupt = InterruptException()
        submission_outdir = joinpath(temporary_root, "submission-interrupt")
        observed_submission_interrupt = Test.@test_logs (
            :warn,
            r"ambiguous response after attempting",
        ) min_level=Logging.Warn match_mode=:any begin
            try
                Mycelia._run_metamdbg(;
                    hifi_reads = reads,
                    outdir = submission_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = interrupt_toolchain,
                    local_runner = _command -> error(
                        "local runner is forbidden",
                    ),
                    submission_runner = (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor,
                    ) -> throw(submission_interrupt),
                )
                nothing
            catch caught
                caught
            end
        end
        Test.@test observed_submission_interrupt === submission_interrupt
    end
end

Test.@testset "metaMDBG pathname publication is identity bound" begin
    mktempdir() do temporary_root
        reads = joinpath(temporary_root, "reads.fastq")
        write(reads, "@read-1\nACGT\n+\nIIII\n")
        selected_input = Mycelia._metamdbg_selected_input(reads, nothing)
        input_contract =
            Mycelia._metamdbg_input_contract(selected_input, 3)

        marker_outdir = joinpath(temporary_root, "output-marker")
        marker_outputs = Mycelia._metamdbg_output_paths(marker_outdir, 21)
        marker_reservation = Mycelia._metamdbg_submission_reservation(
            marker_outputs,
            input_contract,
            21;
            owner_token = "output-marker-hardening",
        )
        function marker_publisher(
                after_hook::Function,
                post_hook::Function,
        )::String
            return Mycelia._publish_metamdbg_output_root_reservation_marker!(
                marker_reservation;
                after_temporary_binding_hook = after_hook,
                post_hardlink_hook = post_hook,
            )
        end
        _test_metamdbg_file_publisher_hardening(
            marker_publisher,
            marker_reservation.output_root_reservation_marker,
            joinpath(temporary_root, "held-output-marker"),
            r"shared output-root marker cleanup failed",
            "output marker",
        )

        contract_outdir = joinpath(temporary_root, "contract-publisher")
        mkpath(contract_outdir)
        contract_outputs =
            Mycelia._metamdbg_output_paths(contract_outdir, 21)
        function contract_publisher(
                after_hook::Function,
                post_hook::Function,
        )::String
            return Mycelia._write_metamdbg_contract!(
                contract_outputs,
                input_contract;
                after_temporary_binding_hook = after_hook,
                post_hardlink_hook = post_hook,
            )
        end
        _test_metamdbg_file_publisher_hardening(
            contract_publisher,
            contract_outputs.contract_marker,
            joinpath(contract_outdir, "held-contract"),
            r"provenance-contract cleanup failed",
            "provenance contract",
        )

        completion_outdir = joinpath(temporary_root, "completion-publisher")
        mkpath(completion_outdir)
        completion_outputs =
            Mycelia._metamdbg_output_paths(completion_outdir, 21)
        completion = (; contents = "{\"fixture\":true}\n")
        function completion_publisher(
                after_hook::Function,
                post_hook::Function,
        )::String
            return Mycelia._write_metamdbg_completion_manifest!(
                completion_outputs,
                completion;
                after_temporary_binding_hook = after_hook,
                post_hardlink_hook = post_hook,
            )
        end
        _test_metamdbg_file_publisher_hardening(
            completion_publisher,
            completion_outputs.completion_marker,
            joinpath(completion_outdir, "held-completion"),
            r"completion-manifest cleanup failed",
            "completion manifest",
        )

        pending_outdir = joinpath(temporary_root, "pending-publisher")
        pending_outputs = Mycelia._metamdbg_output_paths(pending_outdir, 21)
        pending_reservation = Mycelia._metamdbg_submission_reservation(
            pending_outputs,
            input_contract,
            21;
            owner_token = "pending-hardening",
        )
        Mycelia._with_metamdbg_output_lock(pending_outdir) do
            Mycelia._create_metamdbg_submission_reservation!(
                pending_reservation,
                pending_outdir,
            )
        end
        pending_bound = Mycelia._metamdbg_bound_submission_reservation(
            pending_reservation,
            "931",
        )
        write(pending_bound.pending_job_marker, pending_bound.job_contents)
        chmod(pending_bound.pending_job_marker, 0o600)
        pending = Mycelia._metamdbg_pending_submission_job_record(
            pending_reservation,
            pending_bound.pending_job_marker,
        )
        reservation_identity =
            Mycelia._metamdbg_submission_reservation_identity(
                pending_reservation,
            )
        shared_identity = Mycelia._metamdbg_shared_reservation_identity(
            pending_reservation,
        )
        held_pending = joinpath(temporary_root, "held-pending")
        observed_pending_aba = try
            Mycelia._publish_metamdbg_pending_submission_job_record!(
                pending_reservation,
                pending,
                reservation_identity,
                shared_identity;
                pre_hardlink_hook = (path, _identity) -> begin
                    mv(path, held_pending)
                    write(path, "replacement pending record\n")
                    chmod(path, 0o600)
                end,
            )
            nothing
        catch caught
            caught
        end
        Test.@test observed_pending_aba isa ErrorException
        Test.@test occursin(
            "replaced before exact publication",
            sprint(showerror, observed_pending_aba),
        )
        Test.@test read(pending.path, String) ==
                   "replacement pending record\n"
        Test.@test isfile(held_pending)
        Test.@test !ispath(pending_reservation.job_marker)
        rm(pending.path)
        mv(held_pending, pending.path)

        pending_post_primary = ErrorException(
            "synthetic pending post-hardlink failure",
        )
        observed_pending_post = try
            Mycelia._publish_metamdbg_pending_submission_job_record!(
                pending_reservation,
                pending,
                reservation_identity,
                shared_identity;
                post_hardlink_hook = (_path, _identity) ->
                    throw(pending_post_primary),
            )
            nothing
        catch caught
            caught
        end
        Test.@test observed_pending_post === pending_post_primary
        Test.@test isfile(pending.path)
        Test.@test !ispath(pending_reservation.job_marker)

        for phase in (:pre_rename, :post_rename)
            create_outdir = joinpath(
                temporary_root,
                "create-reservation-$(phase)",
            )
            create_outputs =
                Mycelia._metamdbg_output_paths(create_outdir, 21)
            create_reservation = Mycelia._metamdbg_submission_reservation(
                create_outputs,
                input_contract,
                21;
                owner_token = "create-$(phase)-hardening",
            )
            held_reservation = joinpath(
                temporary_root,
                "held-create-$(phase)",
            )
            create_replacement_path = Ref("")
            observed_create_aba = Test.@test_logs (
                :warn,
                r"failed to clean an incomplete submission reservation",
            ) min_level=Logging.Warn match_mode=:any begin
                try
                    Mycelia._with_metamdbg_output_lock(create_outdir) do
                        Mycelia._create_metamdbg_submission_reservation!(
                            create_reservation,
                            create_outdir;
                            pre_rename_hook = temporary -> begin
                                phase == :pre_rename || return nothing
                                create_replacement_path[] = temporary.path
                                mv(temporary.path, held_reservation)
                                mkpath(temporary.path)
                                write(
                                    joinpath(temporary.path, "replacement"),
                                    "retain\n",
                                )
                                return nothing
                            end,
                            post_rename_hook = published -> begin
                                phase == :post_rename || return nothing
                                create_replacement_path[] = published.path
                                mv(published.path, held_reservation)
                                mkpath(published.path)
                                write(
                                    joinpath(published.path, "replacement"),
                                    "retain\n",
                                )
                                return nothing
                            end,
                        )
                    end
                    nothing
                catch caught
                    caught
                end
            end
            Test.@test observed_create_aba isa ErrorException
            Test.@test occursin(
                phase == :pre_rename ?
                "replaced before publication" :
                "replaced after publication",
                sprint(showerror, observed_create_aba),
            )
            Test.@test read(
                joinpath(create_replacement_path[], "replacement"),
                String,
            ) == "retain\n"
            Test.@test isdir(held_reservation)
            Test.@test !ispath(
                create_reservation.output_root_reservation_marker,
            )
            rm(create_replacement_path[]; recursive = true)
            rm(held_reservation; recursive = true)
        end

        for phase in (:pre_private, :post_private, :post_shared)
            remove_outdir = joinpath(
                temporary_root,
                "remove-reservation-$(phase)",
            )
            remove_outputs =
                Mycelia._metamdbg_output_paths(remove_outdir, 21)
            remove_reservation = Mycelia._metamdbg_submission_reservation(
                remove_outputs,
                input_contract,
                21;
                owner_token = "remove-$(phase)-hardening",
            )
            Mycelia._with_metamdbg_output_lock(remove_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    remove_reservation,
                    remove_outdir,
                )
            end
            held_remove = joinpath(temporary_root, "held-remove-$(phase)")
            replacement_path = phase == :pre_private ?
                               remove_reservation.path :
                               remove_reservation.reclaiming_path
            function swap_reservation(path::AbstractString)::Nothing
                mv(path, held_remove)
                mkpath(path)
                write(joinpath(path, "replacement"), "retain\n")
                return nothing
            end
            observed_remove_aba = try
                Mycelia._remove_metamdbg_submission_reservation!(
                    remove_reservation;
                    pre_private_rename_hook = reservation -> begin
                        phase == :pre_private || return nothing
                        swap_reservation(reservation.path)
                    end,
                    post_private_rename_hook = transition -> begin
                        phase == :post_private || return nothing
                        swap_reservation(transition.path)
                    end,
                    post_shared_release_hook = transition -> begin
                        phase == :post_shared || return nothing
                        swap_reservation(transition.path)
                    end,
                )
                nothing
            catch caught
                caught
            end
            Test.@test observed_remove_aba isa ErrorException
            Test.@test occursin(
                phase == :pre_private ?
                "replaced before its identity-bound reclaim" :
                "reclaim transition was replaced",
                sprint(showerror, observed_remove_aba),
            )
            Test.@test read(
                joinpath(replacement_path, "replacement"),
                String,
            ) == "retain\n"
            Test.@test isdir(held_remove)
            marker_exists = Mycelia._output_root_path_entry_exists(
                remove_reservation.output_root_reservation_marker,
            )
            Test.@test marker_exists == (phase != :post_shared)
            rm(replacement_path; recursive = true)
            mv(held_remove, replacement_path)
            current = Mycelia._metamdbg_submission_reservation_from_path(
                replacement_path,
                remove_outdir,
            )
            Mycelia._remove_metamdbg_submission_reservation!(current)
            Test.@test !ispath(replacement_path)
        end

        graph_outdir = joinpath(temporary_root, "graph-alias")
        mkpath(graph_outdir)
        graph_outputs = Mycelia._metamdbg_output_paths(graph_outdir, 21)
        graph_source = joinpath(
            graph_outdir,
            "assemblyGraph_k21_4bps.gfa",
        )
        write(graph_source, "H\tVN:Z:1.0\nS\tcontig-1\tACGT\n")
        write(graph_outputs.graph_alias, "replacement must survive\n")
        observed_graph_path = try
            Mycelia._normalize_metamdbg_graph!(graph_outputs, 21)
            nothing
        catch caught
            caught
        end
        Test.@test observed_graph_path isa ErrorException
        Test.@test occursin(
            "refuses to overwrite an existing graph-alias path",
            sprint(showerror, observed_graph_path),
        )
        Test.@test read(graph_outputs.graph_alias, String) ==
                   "replacement must survive\n"
        rm(graph_outputs.graph_alias)
        held_graph_alias = graph_outputs.graph_alias * ".original"
        replacement_graph_source =
            joinpath(graph_outdir, "replacement-graph.gfa")
        write(
            replacement_graph_source,
            "H\tVN:Z:1.0\nS\treplacement\tTGCA\n",
        )
        observed_graph_replacement = try
            Mycelia._normalize_metamdbg_graph!(
                graph_outputs,
                21;
                post_symlink_hook = (path, _identity) -> begin
                    mv(path, held_graph_alias)
                    symlink(basename(replacement_graph_source), path)
                    return nothing
                end,
            )
            nothing
        catch caught
            caught
        end
        Test.@test observed_graph_replacement isa ErrorException
        Test.@test occursin(
            "graph alias changed after exact creation",
            sprint(showerror, observed_graph_replacement),
        )
        Test.@test readlink(graph_outputs.graph_alias) ==
                   basename(replacement_graph_source)
        Test.@test readlink(held_graph_alias) == basename(graph_source)
        rm(graph_outputs.graph_alias)
        rm(held_graph_alias)
        symlink(basename(graph_source), graph_outputs.graph_alias)
        Test.@test Mycelia._normalize_metamdbg_graph!(graph_outputs, 21) ==
                   graph_outputs.graph_alias
        Test.@test islink(graph_outputs.graph_alias)
    end
end

Test.@testset "metaMDBG public preflight fails before dependencies" begin
    missing_hifi = joinpath(tempdir(), "mycelia-missing-public-hifi.fastq")
    missing_ont = joinpath(tempdir(), "mycelia-missing-public-ont.fastq")
    _test_metamdbg_error(
        () -> Mycelia.run_metamdbg(),
        ArgumentError,
        r"exactly one input technology",
    )
    _test_metamdbg_error(
        () -> Mycelia.run_metamdbg(;
            hifi_reads = missing_hifi,
            ont_reads = missing_ont,
        ),
        ArgumentError,
        r"exactly one input technology",
    )
    _test_metamdbg_error(
        () -> Mycelia.run_metamdbg(; ont_reads = missing_ont),
        ArgumentError,
        r"ont_r10_4_plus=true",
    )
    _test_metamdbg_error(
        () -> Mycelia.run_metamdbg(;
            hifi_reads = missing_hifi,
            ont_r10_4_plus = true,
        ),
        ArgumentError,
        r"applies only when ont_reads",
    )
    for graph_k in (0, -1)
        _test_metamdbg_error(
            () -> Mycelia.run_metamdbg(;
                hifi_reads = missing_hifi,
                graph_k,
            ),
            ArgumentError,
            r"graph_k must be positive",
        )
    end
end

Test.@testset "metaMDBG durable directory and unlink fences" begin
    mktempdir() do temporary_root
        durable_directory = joinpath(temporary_root, "durable-directory")
        mkdir(durable_directory)
        Test.@test Mycelia._fsync_metamdbg_directory(durable_directory) ===
                   nothing

        symlink_path = joinpath(temporary_root, "directory-symlink")
        symlink(durable_directory, symlink_path)
        _test_metamdbg_error(
            () -> Mycelia._fsync_metamdbg_directory(symlink_path),
            ErrorException,
            r"regular, non-symlink directory",
        )

        pre_fsync_path = joinpath(temporary_root, "pre-fsync-replacement")
        pre_fsync_original = pre_fsync_path * ".original"
        mkdir(pre_fsync_path)
        _test_metamdbg_error(
            () -> Mycelia._fsync_metamdbg_directory(
                pre_fsync_path;
                post_descriptor_open_hook = function (
                        path::AbstractString,
                        _identity::NamedTuple,
                )
                    mv(path, pre_fsync_original)
                    mkdir(path)
                    return nothing
                end,
            ),
            ErrorException,
            r"directory path changed before fsync",
        )
        Test.@test isdir(pre_fsync_path)
        Test.@test isdir(pre_fsync_original)

        post_fsync_path = joinpath(temporary_root, "post-fsync-replacement")
        post_fsync_original = post_fsync_path * ".original"
        mkdir(post_fsync_path)
        _test_metamdbg_error(
            () -> Mycelia._fsync_metamdbg_directory(
                post_fsync_path;
                post_fsync_hook = function (
                        path::AbstractString,
                        _identity::NamedTuple,
                )
                    mv(path, post_fsync_original)
                    mkdir(path)
                    return nothing
                end,
            ),
            ErrorException,
            r"directory path changed during fsync",
        )
        Test.@test isdir(post_fsync_path)
        Test.@test isdir(post_fsync_original)

        pid_parent = joinpath(temporary_root, "pid-parent")
        pid_parent_original = pid_parent * ".original"
        mkdir(pid_parent)
        pid_lock = joinpath(pid_parent, "owner.pid")
        write(pid_lock, "111 fixture-host")
        _test_metamdbg_error(
            () -> Mycelia._metamdbg_output_root_pid_lock_snapshot(
                pid_lock;
                post_parent_descriptor_open_hook = function (
                        parent_path::AbstractString,
                        _identity::NamedTuple,
                )
                    mv(parent_path, pid_parent_original)
                    mkdir(parent_path)
                    write(
                        joinpath(parent_path, basename(pid_lock)),
                        "222 replacement-host",
                    )
                    return nothing
                end,
            ),
            ErrorException,
            r"PID lock parent changed after its descriptor was opened",
        )
        Test.@test read(
            joinpath(pid_parent_original, basename(pid_lock)),
            String,
        ) == "111 fixture-host"
        Test.@test read(pid_lock, String) == "222 replacement-host"

        Test.@test Mycelia._metamdbg_linux_renameat2_syscall_number(
            :x86_64,
        ) == 316
        Test.@test Mycelia._metamdbg_linux_renameat2_syscall_number(
            :aarch64,
        ) == 276
        _test_metamdbg_error(
            () -> Mycelia._metamdbg_linux_renameat2_syscall_number(
                :unverified_architecture,
            ),
            ErrorException,
            r"no verified Linux renameat2 syscall number",
        )
        libc_calls = NamedTuple[]
        syscall_calls = NamedTuple[]
        function missing_renameat2_fixture(
                source_descriptor::Cint,
                source::AbstractString,
                destination_descriptor::Cint,
                destination::AbstractString,
                flags::Cuint,
        )::Cint
            push!(libc_calls, (;
                source_descriptor,
                source = String(source),
                destination_descriptor,
                destination = String(destination),
                flags,
            ))
            error("could not load symbol renameat2")
        end
        function renameat2_syscall_fixture(
                syscall_number::Clong,
                source_descriptor::Cint,
                source::AbstractString,
                destination_descriptor::Cint,
                destination::AbstractString,
                flags::Cuint,
        )::Clong
            push!(syscall_calls, (;
                syscall_number,
                source_descriptor,
                source = String(source),
                destination_descriptor,
                destination = String(destination),
                flags,
            ))
            return Clong(-1)
        end
        fallback_result = Mycelia._metamdbg_linux_renameat2(
            Cint(11),
            "source-entry",
            Cint(12),
            "destination-entry",
            Cuint(1);
            libc_runner = missing_renameat2_fixture,
            syscall_runner = renameat2_syscall_fixture,
        )
        Test.@test fallback_result == Cint(-1)
        Test.@test libc_calls == [(;
            source_descriptor = Cint(11),
            source = "source-entry",
            destination_descriptor = Cint(12),
            destination = "destination-entry",
            flags = Cuint(1),
        )]
        Test.@test syscall_calls == [(;
            syscall_number = Mycelia._metamdbg_linux_renameat2_syscall_number(),
            source_descriptor = Cint(11),
            source = "source-entry",
            destination_descriptor = Cint(12),
            destination = "destination-entry",
            flags = Cuint(1),
        )]

        durable_file = joinpath(temporary_root, "durable-file.json")
        write(durable_file, "fixture\n")
        durable_identity = Mycelia._metamdbg_regular_file_identity(durable_file)
        Mycelia._remove_exact_metamdbg_durable_file!(
            durable_file,
            durable_identity,
        )
        Test.@test !ispath(durable_file)

        ordered_file = joinpath(temporary_root, "ordered-fsync.json")
        write(ordered_file, "ordered\n")
        ordered_identity =
            Mycelia._metamdbg_regular_file_identity(ordered_file)
        fsync_order = Symbol[]
        Mycelia._remove_exact_metamdbg_durable_file!(
            ordered_file,
            ordered_identity;
            post_quarantine_destination_fsync_hook =
                (_path::AbstractString, _quarantine::AbstractString) ->
                    push!(fsync_order, :destination),
            pre_source_parent_fsync_hook =
                (_path::AbstractString, _quarantine::AbstractString) ->
                    push!(fsync_order, :source),
        )
        Test.@test fsync_order == [:destination, :source]
        Test.@test !ispath(ordered_file)

        replacement_file = joinpath(temporary_root, "replacement-file.json")
        replacement_original = replacement_file * ".original"
        write(replacement_file, "original\n")
        original_identity =
            Mycelia._metamdbg_regular_file_identity(replacement_file)
        # Keep the original inode linked so Linux cannot recycle it for the
        # replacement fixture before the identity fence runs.
        mv(replacement_file, replacement_original)
        write(replacement_file, "replacement\n")
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                replacement_file,
                original_identity,
            ),
            ErrorException,
            r"refuses to remove a replacement durable file",
        )
        Test.@test read(replacement_original, String) == "original\n"
        Test.@test read(replacement_file, String) == "replacement\n"

        retained_file = joinpath(temporary_root, "retained-file.json")
        write(retained_file, "retained\n")
        retained_identity = Mycelia._metamdbg_regular_file_identity(retained_file)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                retained_file,
                retained_identity;
                remover = _path -> nothing,
            ),
            ErrorException,
            r"retained fail-closed durable removal evidence",
        )
        Test.@test !ispath(retained_file)
        retained_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(retained_file) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        Test.@test read(
            joinpath(retained_quarantine, "payload"),
            String,
        ) == "retained\n"
        Test.@test isfile(joinpath(retained_quarantine, "manifest.json"))

        close_failure_file = joinpath(
            temporary_root,
            "descriptor-close-failure.json",
        )
        write(close_failure_file, "descriptor cleanup fixture\n")
        close_failure_identity =
            Mycelia._metamdbg_regular_file_identity(close_failure_file)
        close_failure_primary = ErrorException(
            "synthetic durable-entry primary failure",
        )
        close_failure_interrupt = InterruptException()
        descriptor_close_calls = Ref(0)
        observed_descriptor_cleanup = Test.@test_logs (
            :warn,
            r"additional cleanup failure",
        ) min_level=Logging.Warn match_mode=:any begin
            try
                Mycelia._remove_exact_metamdbg_durable_file!(
                    close_failure_file,
                    close_failure_identity;
                    post_quarantine_validation_hook = (_path, _payload) ->
                        throw(close_failure_primary),
                    descriptor_closer = descriptor -> begin
                        descriptor_close_calls[] += 1
                        close(descriptor)
                        descriptor_close_calls[] == 1 && return nothing
                        descriptor_close_calls[] == 4 &&
                            throw(close_failure_interrupt)
                        error("synthetic ordinary descriptor close failure")
                    end,
                )
                nothing
            catch caught
                caught
            end
        end
        Test.@test observed_descriptor_cleanup === close_failure_interrupt
        Test.@test descriptor_close_calls[] == 4
        Test.@test !ispath(close_failure_file)
        close_failure_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(close_failure_file) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        Test.@test isfile(joinpath(close_failure_quarantine, "payload"))
        rm(close_failure_quarantine; recursive = true, force = true)

        collision_file = joinpath(temporary_root, "collision-file.json")
        write(collision_file, "expected\n")
        collision_identity =
            Mycelia._metamdbg_regular_file_identity(collision_file)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                collision_file,
                collision_identity;
                pre_quarantine_rename_hook = function (
                        _path::AbstractString,
                        quarantine::AbstractString,
                )
                    write(joinpath(quarantine, "payload"), "foreign\n")
                    return nothing
                end,
            ),
            ErrorException,
            r"retained fail-closed durable removal evidence.*File exists",
        )
        Test.@test read(collision_file, String) == "expected\n"
        collision_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(collision_file) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        Test.@test read(
            joinpath(collision_quarantine, "payload"),
            String,
        ) == "foreign\n"

        fifo_path = joinpath(temporary_root, "durable-entry.fifo")
        Test.@test ccall(
            :mkfifo,
            Cint,
            (Cstring, Base.Cmode_t),
            fifo_path,
            Base.Cmode_t(0o600),
        ) == 0
        fifo_status = lstat(fifo_path)
        parent = Base.Filesystem.open(
            temporary_root,
            Base.JL_O_RDONLY |
            Base.JL_O_DIRECTORY |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
        )
        try
            Test.@test Mycelia._metamdbg_openat_entry_exists(
                parent,
                basename(fifo_path),
            )
            _test_metamdbg_error(
                () -> Mycelia._metamdbg_openat_entry_identity(
                    parent,
                    basename(fifo_path),
                    fifo_path,
                    false,
                ),
                ErrorException,
                r"changed type",
            )
            _test_metamdbg_error(
                () -> Mycelia._remove_exact_metamdbg_durable_file!(
                    fifo_path,
                    (;
                        path = fifo_path,
                        device = fifo_status.device,
                        inode = fifo_status.inode,
                    ),
                ),
                ErrorException,
                r"changed type",
            )
        finally
            close(parent)
        end

        scan_directory = joinpath(temporary_root, "repeatable-scan")
        mkdir(scan_directory)
        write(joinpath(scan_directory, "first"), "first\n")
        scan_handle = Base.Filesystem.open(
            scan_directory,
            Base.JL_O_RDONLY |
            Base.JL_O_DIRECTORY |
            Base.JL_O_NOFOLLOW |
            Base.JL_O_CLOEXEC,
        )
        try
            Test.@test Mycelia._metamdbg_descriptor_directory_entries(
                scan_handle,
            ) == ["first"]
            write(joinpath(scan_directory, "second"), "second\n")
            Test.@test Mycelia._metamdbg_descriptor_directory_entries(
                scan_handle,
            ) == ["first", "second"]
        finally
            close(scan_handle)
        end

        rename_race_file = joinpath(temporary_root, "rename-race.json")
        rename_race_original = rename_race_file * ".original"
        write(rename_race_file, "expected\n")
        rename_race_identity =
            Mycelia._metamdbg_regular_file_identity(rename_race_file)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                rename_race_file,
                rename_race_identity;
                pre_quarantine_rename_hook = function (
                        path::AbstractString,
                        _quarantine::AbstractString,
                )
                    mv(path, rename_race_original)
                    write(path, "replacement\n")
                    return nothing
                end,
            ),
            ErrorException,
            r"moved a replacement into its fail-closed removal quarantine",
        )
        Test.@test !ispath(rename_race_file)
        Test.@test read(rename_race_original, String) == "expected\n"
        rename_race_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(rename_race_file) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        Test.@test read(
            joinpath(rename_race_quarantine, "payload"),
            String,
        ) == "replacement\n"

        unlink_race_file = joinpath(temporary_root, "unlink-race.json")
        unlink_race_original = unlink_race_file * ".original"
        write(unlink_race_file, "expected\n")
        unlink_race_identity =
            Mycelia._metamdbg_regular_file_identity(unlink_race_file)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                unlink_race_file,
                unlink_race_identity;
                pre_quarantine_unlink_hook = function (
                        _path::AbstractString,
                        payload::AbstractString,
                )
                    mv(payload, unlink_race_original)
                    write(payload, "replacement\n")
                    return nothing
                end,
            ),
            ErrorException,
            r"was replaced before unlink",
        )
        Test.@test !ispath(unlink_race_file)
        Test.@test read(unlink_race_original, String) == "expected\n"
        unlink_race_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(unlink_race_file) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        Test.@test read(
            joinpath(unlink_race_quarantine, "payload"),
            String,
        ) == "replacement\n"

        post_validation_file =
            joinpath(temporary_root, "post-validation-race.json")
        post_validation_original = post_validation_file * ".original"
        write(post_validation_file, "expected\n")
        post_validation_identity =
            Mycelia._metamdbg_regular_file_identity(post_validation_file)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                post_validation_file,
                post_validation_identity;
                post_quarantine_validation_hook = function (
                        _path::AbstractString,
                        payload::AbstractString,
                )
                    mv(payload, post_validation_original)
                    write(payload, "replacement\n")
                    return nothing
                end,
            ),
            ErrorException,
            r"was replaced after descriptor validation",
        )
        Test.@test read(post_validation_original, String) == "expected\n"
        post_validation_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(post_validation_file) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        Test.@test read(
            joinpath(post_validation_quarantine, "payload"),
            String,
        ) == "replacement\n"

        final_unlink_directory =
            joinpath(temporary_root, "final-unlink-race")
        final_unlink_original = final_unlink_directory * ".original"
        mkdir(final_unlink_directory)
        write(joinpath(final_unlink_directory, "entry.txt"), "entry\n")
        final_unlink_identity =
            Mycelia._metamdbg_directory_path_identity(final_unlink_directory)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_directory!(
                final_unlink_directory,
                final_unlink_identity;
                recursive = true,
                pre_quarantine_final_unlink_hook = function (
                        _path::AbstractString,
                        payload::AbstractString,
                )
                    mv(payload, final_unlink_original)
                    mkdir(payload)
                    write(joinpath(payload, "replacement.txt"), "foreign\n")
                    return nothing
                end,
            ),
            ErrorException,
            r"directory was replaced before final unlink",
        )
        Test.@test isdir(final_unlink_original)
        Test.@test isempty(readdir(final_unlink_original))
        final_unlink_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(final_unlink_directory) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        Test.@test read(
            joinpath(
                final_unlink_quarantine,
                "payload",
                "replacement.txt",
            ),
            String,
        ) == "foreign\n"

        nested_preopen_directory =
            joinpath(temporary_root, "nested-child-preopen-race")
        nested_preopen_child =
            joinpath(nested_preopen_directory, "nested")
        mkpath(nested_preopen_child)
        write(joinpath(nested_preopen_child, "original.txt"), "original\n")
        nested_preopen_identity =
            Mycelia._metamdbg_directory_path_identity(
                nested_preopen_directory,
            )
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_directory!(
                nested_preopen_directory,
                nested_preopen_identity;
                recursive = true,
                nested_child_pre_open_hook = function (
                        path::AbstractString,
                        _identity::NamedTuple,
                        entry_kind::Symbol,
                )
                    entry_kind == :directory || return nothing
                    held_path = String(path) * ".original"
                    mv(path, held_path)
                    mkdir(path)
                    write(joinpath(path, "replacement.txt"), "replacement\n")
                    return nothing
                end,
            ),
            ErrorException,
            r"descriptor does not bind the expected payload",
        )
        nested_preopen_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(nested_preopen_directory) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        nested_preopen_payload =
            joinpath(nested_preopen_quarantine, "payload")
        Test.@test read(
            joinpath(nested_preopen_payload, "nested", "replacement.txt"),
            String,
        ) == "replacement\n"
        Test.@test read(
            joinpath(
                nested_preopen_payload,
                "nested.original",
                "original.txt",
            ),
            String,
        ) == "original\n"

        nested_preunlink_directory =
            joinpath(temporary_root, "nested-file-preunlink-race")
        mkdir(nested_preunlink_directory)
        write(
            joinpath(nested_preunlink_directory, "nested.txt"),
            "original\n",
        )
        nested_preunlink_identity =
            Mycelia._metamdbg_directory_path_identity(
                nested_preunlink_directory,
            )
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_directory!(
                nested_preunlink_directory,
                nested_preunlink_identity;
                recursive = true,
                nested_child_pre_unlink_hook = function (
                        _path::AbstractString,
                        quarantine_path::AbstractString,
                        _identity::NamedTuple,
                        entry_kind::Symbol,
                )
                    entry_kind == :file || return nothing
                    held_path = String(quarantine_path) * ".original"
                    mv(quarantine_path, held_path)
                    write(quarantine_path, "replacement\n")
                    return nothing
                end,
            ),
            ErrorException,
            r"quarantined child was replaced before exact unlink",
        )
        nested_preunlink_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(nested_preunlink_directory) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        nested_preunlink_payload =
            joinpath(nested_preunlink_quarantine, "payload")
        nested_preunlink_replacement = only(filter(
            path -> startswith(
                basename(path),
                ".mycelia-child.removing.",
            ) && !endswith(path, ".original"),
            readdir(nested_preunlink_payload; join = true),
        ))
        Test.@test read(nested_preunlink_replacement, String) ==
                   "replacement\n"
        Test.@test read(
            nested_preunlink_replacement * ".original",
            String,
        ) == "original\n"

        quarantine_swap_file =
            joinpath(temporary_root, "quarantine-swap.json")
        write(quarantine_swap_file, "expected\n")
        quarantine_swap_identity =
            Mycelia._metamdbg_regular_file_identity(quarantine_swap_file)
        quarantine_replacement = Ref("")
        quarantine_original = Ref("")
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                quarantine_swap_file,
                quarantine_swap_identity;
                pre_quarantine_directory_unlink_hook = function (
                        _path::AbstractString,
                        quarantine::AbstractString,
                )
                    original = quarantine * ".original"
                    mv(quarantine, original)
                    mkdir(quarantine)
                    chmod(quarantine, 0o700)
                    quarantine_replacement[] = String(quarantine)
                    quarantine_original[] = original
                    return nothing
                end,
            ),
            ErrorException,
            r"removal quarantine was replaced before final unlink",
        )
        Test.@test !ispath(quarantine_swap_file)
        Test.@test isdir(quarantine_replacement[])
        Test.@test isdir(quarantine_original[])
        Test.@test isempty(readdir(quarantine_replacement[]))
        Test.@test isempty(readdir(quarantine_original[]))

        dead_child = run(`sleep 0.1`; wait = false)
        dead_pid = getpid(dead_child)
        wait(dead_child)
        hostname = Base.gethostname()
        Test.@test !FileWatching.Pidfile.isvalidpid(
            hostname,
            Cuint(dead_pid),
        )
        pid_race_file = joinpath(temporary_root, "pid-content-race.pid")
        write(pid_race_file, "$(dead_pid) $(hostname)")
        pid_race_identity =
            Mycelia._metamdbg_regular_file_identity(pid_race_file)
        live_contents = "$(getpid()) $(hostname)"
        _test_metamdbg_error(
            () -> Mycelia._remove_dead_metamdbg_output_root_pid_lock!(
                pid_race_file;
                pre_quarantine_unlink_hook = function (
                        _path::AbstractString,
                        payload::AbstractString,
                )
                    open(payload, "w") do output
                        write(output, live_contents)
                    end
                    return nothing
                end,
            ),
            ErrorException,
            r"quarantined output-root PID lock contents changed",
        )
        Test.@test !ispath(pid_race_file)
        pid_race_quarantine = only(filter(
            path -> startswith(
                basename(path),
                basename(pid_race_file) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))
        pid_payload = joinpath(pid_race_quarantine, "payload")
        Test.@test read(pid_payload, String) == live_contents
        Test.@test stat(pid_payload).inode == pid_race_identity.inode

        nonempty_directory = joinpath(temporary_root, "nonempty-directory")
        mkdir(nonempty_directory)
        write(joinpath(nonempty_directory, "blocker"), "block\n")
        nonempty_identity =
            Mycelia._metamdbg_directory_path_identity(nonempty_directory)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_directory!(
                nonempty_directory,
                nonempty_identity,
            ),
            ErrorException,
            r"refuses nonrecursive removal of a nonempty durable directory",
        )
        Test.@test isdir(nonempty_directory)
        Test.@test read(
            joinpath(nonempty_directory, "blocker"),
            String,
        ) == "block\n"
        Test.@test isempty(filter(
            path -> startswith(
                basename(path),
                basename(nonempty_directory) * ".removing.",
            ),
            readdir(temporary_root; join = true),
        ))

        durable_subdirectory = joinpath(temporary_root, "durable-subdirectory")
        mkdir(durable_subdirectory)
        write(joinpath(durable_subdirectory, "entry.txt"), "entry\n")
        durable_subdirectory_identity =
            Mycelia._metamdbg_directory_path_identity(durable_subdirectory)
        Mycelia._remove_exact_metamdbg_durable_directory!(
            durable_subdirectory,
            durable_subdirectory_identity;
            recursive = true,
        )
        Test.@test !ispath(durable_subdirectory)
    end
end

Test.@testset "metaMDBG reports interrupted removal evidence on restart" begin
    mktempdir() do temporary_root
        fresh_nested_outdir = joinpath(
            temporary_root,
            "missing",
            "child",
            "restart-output",
        )
        Test.@test Mycelia._metamdbg_removal_quarantine_evidence(
            fresh_nested_outdir,
        ) == String[]
        Test.@test Mycelia._require_no_metamdbg_removal_quarantine_evidence!(
            fresh_nested_outdir,
        ) === nothing
    end

    mktempdir() do temporary_root
        outdir = joinpath(temporary_root, "restart-output")
        lock_path = Mycelia._metamdbg_output_lock_path(outdir)
        Base.Filesystem.mkdir(lock_path; mode = 0o700)
        lock_identity =
            Mycelia._metamdbg_directory_path_identity(lock_path)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_directory!(
                lock_path,
                lock_identity;
                pre_quarantine_rename_hook =
                    (_path::AbstractString, _quarantine::AbstractString) ->
                        error("synthetic pre-rename interruption"),
            ),
            ErrorException,
            r"retained fail-closed durable removal evidence.*pre-rename",
        )
        Test.@test isdir(lock_path)
        evidence =
            Mycelia._metamdbg_removal_quarantine_evidence(outdir)
        Test.@test length(evidence) == 1
        Test.@test startswith(
            basename(only(evidence)),
            basename(lock_path) * ".removing.",
        )
        _test_metamdbg_error(
            () -> Mycelia.inspect_metamdbg_submission_reservations(outdir),
            ErrorException,
            r"fail-closed \.removing evidence.*refusing automatic restart",
        )
    end

    mktempdir() do temporary_root
        outdir = joinpath(temporary_root, "nested-restart-output")
        reservation_path = joinpath(
            temporary_root,
            Mycelia._metamdbg_submission_reservation_prefix(outdir) *
            repeat("a", 64),
        )
        Base.Filesystem.mkdir(reservation_path; mode = 0o700)
        orphan = joinpath(
            reservation_path,
            "job.json.removing." * repeat("b", 32),
        )
        Base.Filesystem.mkdir(orphan; mode = 0o700)
        Test.@test Mycelia._metamdbg_removal_quarantine_evidence(outdir) ==
                   [orphan]
        _test_metamdbg_error(
            () -> Mycelia._require_no_active_metamdbg_submission_reservation!(
                outdir,
            ),
            ErrorException,
            r"fail-closed \.removing evidence",
        )
    end

    mktempdir() do temporary_root
        outdir = joinpath(temporary_root, "unrelated-target-output")
        unrelated_file = joinpath(temporary_root, "unrelated.json")
        write(unrelated_file, "unrelated\n")
        unrelated_identity =
            Mycelia._metamdbg_regular_file_identity(unrelated_file)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                unrelated_file,
                unrelated_identity;
                pre_quarantine_unlink_hook =
                    (_path::AbstractString, _payload::AbstractString) ->
                        error("synthetic unrelated interruption"),
            ),
            ErrorException,
            r"retained fail-closed durable removal evidence.*unrelated",
        )
        Test.@test Mycelia._metamdbg_removal_quarantine_evidence(outdir) ==
                   String[]
        Test.@test Mycelia._require_no_metamdbg_removal_quarantine_evidence!(
            outdir,
        ) === nothing
    end

    mktempdir() do temporary_root
        outdir = joinpath(temporary_root, "reservation-scan-race-output")
        reservation_path = joinpath(
            temporary_root,
            Mycelia._metamdbg_submission_reservation_prefix(outdir) *
            repeat("c", 64),
        )
        mkdir(reservation_path)
        _test_metamdbg_error(
            () -> Mycelia._metamdbg_removal_quarantine_evidence(
                outdir;
                post_reservation_scan_hook = function (
                        path::AbstractString,
                        _identity::NamedTuple,
                )
                    mv(path, path * ".original")
                    mkdir(path)
                    mkdir(joinpath(
                        path,
                        "job.json.removing." * repeat("d", 32),
                    ))
                    return nothing
                end,
            ),
            ErrorException,
            r"submission reservation changed during removal-evidence scan",
        )
        Test.@test isdir(reservation_path * ".original")
        Test.@test isdir(reservation_path)
    end

    for mutation in (:add, :remove)
        mktempdir() do temporary_root
            outdir = joinpath(
                temporary_root,
                "reservation-entry-set-$(mutation)-output",
            )
            reservation_path = joinpath(
                temporary_root,
                Mycelia._metamdbg_submission_reservation_prefix(outdir) *
                repeat("e", 64),
            )
            mkdir(reservation_path)
            orphan = joinpath(
                reservation_path,
                "job.json.removing." * repeat("f", 32),
            )
            mutation == :remove && mkdir(orphan)
            reservation_identity = stat(reservation_path)
            _test_metamdbg_error(
                () -> Mycelia._metamdbg_removal_quarantine_evidence(
                    outdir;
                    post_reservation_scan_hook = function (
                            _path::AbstractString,
                            _identity::NamedTuple,
                    )
                        mutation == :add ? mkdir(orphan) : rm(orphan)
                        return nothing
                    end,
                ),
                ErrorException,
                r"submission reservation entry set changed",
            )
            final_reservation_identity = stat(reservation_path)
            Test.@test final_reservation_identity.device ==
                       reservation_identity.device
            Test.@test final_reservation_identity.inode ==
                       reservation_identity.inode
            Test.@test isdir(orphan) == (mutation == :add)
        end
    end

    mktempdir() do temporary_root
        outdir = joinpath(temporary_root, "outdir-entry-set-race-output")
        mkdir(outdir)
        reservation_path = joinpath(
            temporary_root,
            Mycelia._metamdbg_submission_reservation_prefix(outdir) *
            repeat("1", 64),
        )
        mkdir(reservation_path)
        orphan = joinpath(
            outdir,
            "assembly.fasta.removing." * repeat("2", 32),
        )
        _test_metamdbg_error(
            () -> Mycelia._metamdbg_removal_quarantine_evidence(
                outdir;
                post_reservation_scan_hook = function (
                        _path::AbstractString,
                        _identity::NamedTuple,
                )
                    mkdir(orphan)
                    return nothing
                end,
            ),
            ErrorException,
            r"output directory removal-evidence entry set changed",
        )
        Test.@test isdir(orphan)
    end

    mktempdir() do temporary_root
        outdir = joinpath(temporary_root, "parent-entry-set-race-output")
        first_reservation = joinpath(
            temporary_root,
            Mycelia._metamdbg_submission_reservation_prefix(outdir) *
            repeat("3", 64),
        )
        second_reservation = joinpath(
            temporary_root,
            Mycelia._metamdbg_submission_reservation_prefix(outdir) *
            repeat("4", 64),
        )
        mkdir(first_reservation)
        _test_metamdbg_error(
            () -> Mycelia._metamdbg_removal_quarantine_evidence(
                outdir;
                post_reservation_scan_hook = function (
                        _path::AbstractString,
                        _identity::NamedTuple,
                )
                    mkdir(second_reservation)
                    return nothing
                end,
            ),
            ErrorException,
            r"output-root removal-evidence entry set changed",
        )
        Test.@test isdir(first_reservation)
        Test.@test isdir(second_reservation)
    end

    mktempdir() do temporary_root
        reads = joinpath(temporary_root, "contracted-reads.fastq")
        write(reads, "@read-1\nACGT\n+\nIIII\n")
        outdir = joinpath(temporary_root, "contracted-restart-output")
        mkdir(outdir)
        selected_input = Mycelia._metamdbg_selected_input(reads, nothing)
        input_contract = Mycelia._metamdbg_input_contract(selected_input, 3)
        outputs = Mycelia._metamdbg_output_paths(outdir, 21)
        Mycelia._write_metamdbg_contract!(outputs, input_contract)
        Test.@test Mycelia._require_metamdbg_contract!(
            outputs,
            input_contract,
        ) == outputs.contract_marker
        write(outputs.contigs_plain, ">contig-1\nACGT\n")
        contigs_identity =
            Mycelia._metamdbg_regular_file_identity(outputs.contigs_plain)
        _test_metamdbg_error(
            () -> Mycelia._remove_exact_metamdbg_durable_file!(
                outputs.contigs_plain,
                contigs_identity;
                pre_quarantine_unlink_hook =
                    (_path::AbstractString, _payload::AbstractString) ->
                        error("synthetic contracted-output interruption"),
            ),
            ErrorException,
            r"retained fail-closed durable removal evidence.*contracted-output",
        )
        evidence = Mycelia._metamdbg_removal_quarantine_evidence(outdir)
        Test.@test length(evidence) == 1
        Test.@test dirname(only(evidence)) == outdir
        Test.@test startswith(
            basename(only(evidence)),
            basename(outputs.contigs_plain) * ".removing.",
        )
        _test_metamdbg_error(
            () -> Mycelia._require_no_active_metamdbg_submission_reservation!(
                outdir,
            ),
            ErrorException,
            r"fail-closed \.removing evidence",
        )
    end
end

Test.@testset "metaMDBG dead-process recovery rechecks runtime evidence" begin
    mktempdir() do temporary_root
        reads = joinpath(temporary_root, "reads.fastq")
        write(reads, "@read-1\nACGT\n+\nIIII\n")
        outdir = joinpath(temporary_root, "runtime-race-output")
        selected_input = Mycelia._metamdbg_selected_input(reads, nothing)
        input_contract = Mycelia._metamdbg_input_contract(selected_input, 3)
        outputs = Mycelia._metamdbg_output_paths(outdir, 21)
        reservation = Mycelia._metamdbg_submission_reservation(
            outputs,
            input_contract,
            21;
            owner_token = "runtime-recheck-owner",
        )
        Mycelia._with_metamdbg_output_domain_lock(outdir) do
            Mycelia._create_metamdbg_submission_reservation!(
                reservation,
                outdir,
            )
        end
        private_lock = Mycelia._metamdbg_output_lock_path(outdir)
        Base.Filesystem.mkdir(private_lock; mode = 0o700)
        cleanup_reservation =
            Mycelia._publish_metamdbg_lifecycle_cleanup_reservation!(outdir)
        private_identity = Mycelia._metamdbg_output_lock_identity(private_lock)
        cleanup_identity =
            Mycelia._metamdbg_output_lock_identity(cleanup_reservation)
        pid_lock = Mycelia._output_root_reservation_lock_path_from_canonical(
            Mycelia._metamdbg_canonical_output_path(outdir),
        )
        dead_pid = Base.typemax(Cint)
        Test.@test !FileWatching.Pidfile.isvalidpid(
            Base.gethostname(),
            Cuint(dead_pid),
        )
        write(pid_lock, "$(dead_pid) $(Base.gethostname())")
        pending_recovery_called = Ref(false)
        hook_called = Ref(false)
        _test_metamdbg_error(
            () -> Mycelia._inspect_metamdbg_submission_reservations(
                outdir;
                confirm_process_dead = true,
                pending_recovery_function = _outdir -> begin
                    pending_recovery_called[] = true
                    return nothing
                end,
                post_recovery_pid_acquisition_hook = _outdir -> begin
                    hook_called[] = true
                    Base.Filesystem.mkdir(
                        reservation.runtime_output_root_reservation_marker;
                        mode = 0o700,
                    )
                    Mycelia._fsync_metamdbg_directory(dirname(
                        reservation.runtime_output_root_reservation_marker,
                    ))
                    return nothing
                end,
            ),
            ErrorException,
            r"runtime scheduler ownership appeared after the replacement recovery PID",
        )
        Test.@test hook_called[]
        Test.@test !pending_recovery_called[]
        Test.@test Mycelia._metamdbg_output_lock_identity(private_lock) ==
                   private_identity
        Test.@test Mycelia._metamdbg_output_lock_identity(cleanup_reservation) ==
                   cleanup_identity
        Test.@test isdir(reservation.runtime_output_root_reservation_marker)
    end
end

Test.@testset "metaMDBG preserves federated scheduler identity" begin
    mktempdir() do temporary_root
        reads = joinpath(temporary_root, "reads.fastq")
        write(reads, "@read-1\nACGT\n+\nIIII\n")
        results = NamedTuple[]
        for (cluster, pending_state) in (
                ("cluster-a", :partial),
                ("cluster-b", :complete),
            )
            outdir = joinpath(temporary_root, cluster)
            release_calls = NamedTuple[]
            result = Mycelia._run_metamdbg(;
                hifi_reads = reads,
                outdir,
                executor = Mycelia.SlurmExecutor(dry_run = false),
                site = :scg,
                dependency_checker = () ->
                    Mycelia._require_metamdbg_package_version(NamedTuple[
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
                    ]),
                local_runner = _command -> error("local runner is forbidden"),
                submission_runner = function (
                        job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor,
                )
                    Test.@test occursin("SLURM_CLUSTER_NAME", job.cmd)
                    return Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = "42001",
                        job_cluster = cluster,
                        stdout = "42001;$(cluster)\n",
                    )
                end,
                submission_job_binder = function (
                        reservation::NamedTuple,
                        job_id::AbstractString,
                        job_cluster::AbstractString,
                )
                    Test.@test job_cluster == cluster
                    return Mycelia._bind_metamdbg_submission_job_after_submit!(
                        reservation,
                        job_id,
                        job_cluster,
                    )
                end,
                submission_release_runner = function (
                        job_id::AbstractString;
                        job_cluster::Union{Nothing, AbstractString} = nothing,
                )
                    push!(release_calls, (;
                        job_id = String(job_id),
                        job_cluster = job_cluster === nothing ?
                                      nothing : String(job_cluster),
                    ))
                    return String(job_id)
                end,
            )
            Test.@test release_calls == [(;
                job_id = "42001",
                job_cluster = cluster,
            )]
            Test.@test result.submission.job_cluster == cluster
            Test.@test result.submission_reservation.job_cluster == cluster
            metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(outdir),
            )
            Test.@test metadata.job_id == "42001"
            Test.@test metadata.job_cluster == cluster
            job_record = JSON.parse(read(
                joinpath(metadata.path, "job.json"),
                String,
            ))
            Test.@test job_record["job_id"] == "42001"
            Test.@test job_record["job_cluster"] == cluster
            capability =
                Mycelia._metamdbg_output_root_reservation_capability(
                    result.submission_reservation.workflow_signature,
                    result.submission_reservation.owner_token,
                )
            pending_a = Mycelia._metamdbg_pending_submission_job_path(
                outdir,
                capability,
                "42001",
                "cluster-a",
            )
            pending_b = Mycelia._metamdbg_pending_submission_job_path(
                outdir,
                capability,
                "42001",
                "cluster-b",
            )
            Test.@test pending_a != pending_b
            push!(results, result)
            wrong_cluster = cluster == "cluster-a" ? "cluster-b" : "cluster-a"
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    metadata;
                    owner_token = metadata.owner_token,
                    job_id = metadata.job_id,
                    job_cluster = wrong_cluster,
                    confirm_cancelled = true,
                ),
                ErrorException,
                r"scheduler cluster does not match",
            )
            reclaimed = Mycelia.reclaim_metamdbg_submission_reservation!(
                metadata;
                owner_token = metadata.owner_token,
                job_id = metadata.job_id,
                job_cluster = metadata.job_cluster,
                confirm_cancelled = true,
            )
            Test.@test reclaimed.job_cluster == cluster
        end
        Test.@test results[1].submission.job_id == results[2].submission.job_id
        Test.@test results[1].submission.job_cluster !=
                   results[2].submission.job_cluster
    end
end

Test.@testset "metaMDBG recovers schema-1 durable job evidence" begin
    mktempdir() do temporary_root
        reads = joinpath(temporary_root, "reads.fastq")
        write(reads, "@read-1\nACGT\n+\nIIII\n")
        selected_input = Mycelia._metamdbg_selected_input(reads, nothing)
        input_contract = Mycelia._metamdbg_input_contract(selected_input, 3)
        legacy_schema =
            Mycelia._METAMDBG_LEGACY_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION

        committed_outdir = joinpath(temporary_root, "legacy-committed")
        committed_outputs =
            Mycelia._metamdbg_output_paths(committed_outdir, 21)
        committed_reservation = Mycelia._metamdbg_submission_reservation(
            committed_outputs,
            input_contract,
            21;
            owner_token = "legacy-committed-owner",
        )
        Mycelia._with_metamdbg_output_domain_lock(committed_outdir) do
            Mycelia._create_metamdbg_submission_reservation!(
                committed_reservation,
                committed_outdir,
            )
        end
        legacy_committed = Mycelia._metamdbg_bound_submission_reservation(
            committed_reservation,
            "51001",
            nothing,
            legacy_schema,
        )
        open(committed_reservation.job_marker, "w") do output
            write(output, legacy_committed.job_contents)
            chmod(committed_reservation.job_marker, 0o600)
            Mycelia._fsync_metamdbg_file(
                output,
                committed_reservation.job_marker,
            )
        end
        Mycelia._fsync_metamdbg_directory(committed_reservation.path)
        committed_metadata = only(
            Mycelia.inspect_metamdbg_submission_reservations(
                committed_outdir,
            ),
        )
        Test.@test committed_metadata.job_id == "51001"
        Test.@test committed_metadata.job_cluster === nothing
        Test.@test committed_metadata.job_schema_version == legacy_schema
        committed_reclaimed =
            Mycelia.reclaim_metamdbg_submission_reservation!(
                committed_metadata;
                owner_token = committed_metadata.owner_token,
                job_id = committed_metadata.job_id,
                confirm_cancelled = true,
            )
        Test.@test committed_reclaimed.job_cluster === nothing

        for pending_state in (:complete, :partial, :empty)
            outdir = joinpath(temporary_root, "legacy-pending-$(pending_state)")
            outputs = Mycelia._metamdbg_output_paths(outdir, 21)
            reservation = Mycelia._metamdbg_submission_reservation(
                outputs,
                input_contract,
                21;
                owner_token = "legacy-pending-$(pending_state)-owner",
            )
            Mycelia._with_metamdbg_output_domain_lock(outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    reservation,
                    outdir,
                )
            end
            pending_path = Mycelia._metamdbg_pending_submission_job_path(
                outdir,
                reservation.output_root_reservation_capability,
                "51002",
            )
            legacy_contents = Mycelia._metamdbg_submission_job_contents(
                reservation,
                "51002",
                nothing,
                legacy_schema,
            )
            pending_contents = if pending_state == :complete
                legacy_contents
            elseif pending_state == :partial
                schema_marker = findfirst(
                    "\"schema_version\":1",
                    legacy_contents,
                )
                schema_marker isa UnitRange || error(
                    "Legacy schema marker is missing from the test fixture.",
                )
                legacy_contents[1:last(schema_marker)]
            else
                ""
            end
            write(pending_path, pending_contents)
            chmod(pending_path, 0o600)
            Mycelia._fsync_metamdbg_directory(dirname(pending_path))
            pending_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(outdir),
            )
            expected_schema = pending_state == :empty ?
                              Mycelia._METAMDBG_SUBMISSION_RESERVATION_JOB_SCHEMA_VERSION :
                              legacy_schema
            Test.@test pending_metadata.job_id === nothing
            Test.@test pending_metadata.pending_job_id == "51002"
            Test.@test pending_metadata.pending_job_cluster === nothing
            Test.@test pending_metadata.pending_job_schema_version ==
                       expected_schema
            bound = Mycelia.bind_metamdbg_submission_reservation_job!(
                pending_metadata;
                owner_token = pending_metadata.owner_token,
                job_id = pending_metadata.pending_job_id,
                confirm_submitted = true,
            )
            Test.@test bound.job_schema_version == expected_schema
            Test.@test !ispath(pending_path)
            bound_record = JSON.parse(read(
                joinpath(bound.path, "job.json"),
                String,
            ))
            Test.@test bound_record["schema_version"] == expected_schema
            Test.@test !haskey(bound_record, "job_cluster") ==
                       (expected_schema == legacy_schema)
            rebound_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(outdir),
            )
            Test.@test rebound_metadata.job_schema_version == expected_schema
            Mycelia.reclaim_metamdbg_submission_reservation!(
                rebound_metadata;
                owner_token = rebound_metadata.owner_token,
                job_id = rebound_metadata.job_id,
                confirm_cancelled = true,
            )
        end
    end
end

Test.@testset "metaMDBG recovers cluster-qualified pending jobs" begin
    mktempdir() do temporary_root
        reads = joinpath(temporary_root, "reads.fastq")
        write(reads, "@read-1\nACGT\n+\nIIII\n")
        selected_input = Mycelia._metamdbg_selected_input(reads, nothing)
        input_contract = Mycelia._metamdbg_input_contract(selected_input, 3)
        recovered = NamedTuple[]
        for (cluster, pending_state) in (
                ("cluster-a", :complete),
                ("cluster-b", :partial),
        )
            outdir = joinpath(temporary_root, cluster)
            outputs = Mycelia._metamdbg_output_paths(outdir, 21)
            reservation = Mycelia._metamdbg_submission_reservation(
                outputs,
                input_contract,
                21;
                owner_token = "pending-$(cluster)-owner",
            )
            Mycelia._with_metamdbg_output_domain_lock(outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    reservation,
                    outdir,
                )
            end
            pending_path = Mycelia._metamdbg_pending_submission_job_path(
                outdir,
                reservation.output_root_reservation_capability,
                "52001",
                cluster,
            )
            pending_contents = Mycelia._metamdbg_submission_job_contents(
                reservation,
                "52001",
                cluster,
            )
            persisted_contents = if pending_state == :partial
                partial_length = max(1, length(pending_contents) - 12)
                pending_contents[1:partial_length]
            else
                pending_contents
            end
            write(pending_path, persisted_contents)
            chmod(pending_path, 0o600)
            Mycelia._fsync_metamdbg_directory(dirname(pending_path))
            metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(outdir),
            )
            Test.@test metadata.pending_job_id == "52001"
            Test.@test metadata.pending_job_cluster == cluster
            Test.@test metadata.pending_job_complete ==
                       (pending_state == :complete)
            wrong_cluster = cluster == "cluster-a" ? "cluster-b" : "cluster-a"
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    metadata;
                    owner_token = metadata.owner_token,
                    job_id = metadata.pending_job_id,
                    job_cluster = wrong_cluster,
                    confirm_submitted = true,
                ),
                ErrorException,
                r"cluster does not match the exact pending job record",
            )
            bound = Mycelia.bind_metamdbg_submission_reservation_job!(
                metadata;
                owner_token = metadata.owner_token,
                job_id = metadata.pending_job_id,
                job_cluster = cluster,
                confirm_submitted = true,
            )
            Test.@test bound.job_id == "52001"
            Test.@test bound.job_cluster == cluster
            Test.@test !ispath(pending_path)
            push!(recovered, bound)
            rebound_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(outdir),
            )
            Mycelia.reclaim_metamdbg_submission_reservation!(
                rebound_metadata;
                owner_token = rebound_metadata.owner_token,
                job_id = rebound_metadata.job_id,
                job_cluster = rebound_metadata.job_cluster,
                confirm_cancelled = true,
            )
        end
        Test.@test recovered[1].job_id == recovered[2].job_id
        Test.@test recovered[1].job_cluster != recovered[2].job_cluster
    end
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

function _kill_metamdbg_reservation_publication!(
        outdir::AbstractString,
        ready_path::AbstractString,
        phase::Symbol,
        ;
        live_check::Function = (_process::Base.Process) -> nothing,
)::Nothing
    phase in (:pre_rename, :post_rename) || throw(ArgumentError(
        "phase must be :pre_rename or :post_rename.",
    ))
    active_project = Base.active_project()
    active_project isa AbstractString || error(
        "A project is required for the metaMDBG SIGKILL fixture.",
    )
    script_path = joinpath(
        dirname(String(outdir)),
        "metamdbg-$(phase)-sigkill-fixture.jl",
    )
    hook_keyword = phase == :pre_rename ?
                   "pre_rename_hook" : "post_rename_hook"
    script = """
    ENV["LD_LIBRARY_PATH"] = ""
    const Mycelia = Base.require(Main, :Mycelia)

    outdir = $(repr(String(outdir)))
    ready_path = $(repr(String(ready_path)))
    outputs = Mycelia._metamdbg_output_paths(outdir, 21)
    reservation = Mycelia._metamdbg_submission_reservation(
        outputs,
        (; signature = repeat("a", 64)),
        21;
        owner_token = "$(phase)-sigkill-owner",
    )
    publication_hook = function (_reservation::NamedTuple)
        write(ready_path, "ready\\n")
        while true
            sleep(60)
        end
    end
    Mycelia._with_metamdbg_output_domain_lock(outdir) do
        Mycelia._create_metamdbg_submission_reservation!(
            reservation,
            outdir;
            $(hook_keyword) = publication_hook,
        )
    end
    """
    write(script_path, script)
    project_directory = dirname(String(active_project))
    project_command =
        `$(Base.julia_cmd()) --project=$(project_directory)`
    command =
        `$project_command --startup-file=no --compiled-modules=yes $(script_path)`
    process = run(ignorestatus(command); wait = false)
    wait_status = Base.timedwait(
        () -> isfile(ready_path) || !Base.process_running(process),
        600.0;
        pollint = 0.05,
    )
    if wait_status != :ok || !isfile(ready_path)
        Base.process_running(process) && Base.kill(process, Base.SIGKILL)
        wait(process)
        error(
            "metaMDBG $(phase) SIGKILL fixture did not reach its publication " *
            "hook (exit=$(process.exitcode), signal=$(process.termsignal)).",
        )
    end
    try
        live_check(process)
    finally
        Base.process_running(process) && Base.kill(process, Base.SIGKILL)
        wait(process)
    end
    Test.@test !success(process)
    Test.@test process.termsignal == 9
    return nothing
end

function _kill_metamdbg_job_record_publication!(
        outdir::AbstractString,
        reads::AbstractString,
        phase::Symbol,
        ;
        live_check::Function = (
            _process::Base.Process,
            _reservation::NamedTuple,
        ) -> nothing,
)::NamedTuple
    phase in (:pending_name, :pre_publication, :post_publication) || throw(
        ArgumentError(
            "phase must be :pending_name, :pre_publication, or " *
            ":post_publication.",
        ),
    )
    active_project = Base.active_project()
    active_project isa AbstractString || error(
        "A project is required for the metaMDBG job-bind SIGKILL fixture.",
    )
    normalized_outdir = String(outdir)
    outputs = Mycelia._metamdbg_output_paths(normalized_outdir, 21)
    input_contract = Mycelia._metamdbg_input_contract(
        Mycelia._metamdbg_selected_input(String(reads), nothing),
        3,
    )
    reservation = Mycelia._metamdbg_submission_reservation(
        outputs,
        input_contract,
        21;
        owner_token = "$(phase)-job-bind-owner",
    )
    Mycelia._with_metamdbg_output_lock(normalized_outdir) do
        Mycelia._create_metamdbg_submission_reservation!(
            reservation,
            normalized_outdir,
        )
    end
    ready_path = joinpath(
        dirname(normalized_outdir),
        "metamdbg-$(phase)-job-bind.ready",
    )
    script_path = joinpath(
        dirname(normalized_outdir),
        "metamdbg-$(phase)-job-bind.jl",
    )
    hook_keyword = if phase == :pending_name
        "post_pending_job_record_creation_hook"
    elseif phase == :pre_publication
        "pre_job_record_publication_hook"
    else
        "post_job_record_publication_hook"
    end
    script = """
    ENV["LD_LIBRARY_PATH"] = ""
    const Mycelia = Base.require(Main, :Mycelia)

    outdir = $(repr(normalized_outdir))
    ready_path = $(repr(ready_path))
    reservation_path = $(repr(reservation.path))
    reservation = Mycelia._metamdbg_submission_reservation_from_path(
        reservation_path,
        outdir,
    )
    pause_hook = function (_args...)
        write(ready_path, "ready\\n")
        while true
            sleep(60)
        end
    end
    Mycelia._with_metamdbg_output_domain_lock(
        outdir;
        allowed_same_root_locks = (
            reservation.output_root_reservation_marker,
        ),
    ) do
        Mycelia._bind_metamdbg_submission_job!(
            reservation,
            "771";
            $(hook_keyword) = pause_hook,
        )
    end
    """
    write(script_path, script)
    project_directory = dirname(String(active_project))
    process = run(
        `$(Base.julia_cmd()) --project=$(project_directory) --startup-file=no --compiled-modules=yes $(script_path)`;
        wait = false,
    )
    wait_status = Base.timedwait(
        () -> isfile(ready_path) || !Base.process_running(process),
        600.0;
        pollint = 0.05,
    )
    if wait_status != :ok || !isfile(ready_path)
        Base.process_running(process) && Base.kill(process, Base.SIGKILL)
        wait(process)
        error(
            "metaMDBG $(phase) job-bind fixture did not reach its " *
            "publication hook (exit=$(process.exitcode), " *
            "signal=$(process.termsignal)).",
        )
    end
    try
        live_check(process, reservation)
    finally
        Base.process_running(process) && Base.kill(process, Base.SIGKILL)
        wait(process)
    end
    Test.@test !success(process)
    Test.@test process.termsignal == 9
    return reservation
end

function _kill_metamdbg_runtime_transition!(
        outdir::AbstractString,
        reads::AbstractString,
        phase::Symbol,
        ;
        live_check::Function = (
            _process::Base.Process,
            _reservation::NamedTuple,
        ) -> nothing,
)::NamedTuple
    phase in (
        :post_runtime_marker,
        :post_private_rename,
        :post_queued_release,
        :post_runtime_release,
        :post_consumed_rename,
    ) || throw(
        ArgumentError(
            "phase must be a supported runtime transition hook.",
        ),
    )
    normalized_outdir = String(outdir)
    outputs = Mycelia._metamdbg_output_paths(normalized_outdir, 21)
    input_contract = Mycelia._metamdbg_input_contract(
        Mycelia._metamdbg_selected_input(String(reads), nothing),
        3,
    )
    reservation = Mycelia._metamdbg_submission_reservation(
        outputs,
        input_contract,
        21;
        owner_token = "$(phase)-runtime-owner",
    )
    Mycelia._with_metamdbg_output_lock(normalized_outdir) do
        Mycelia._create_metamdbg_submission_reservation!(
            reservation,
            normalized_outdir,
        )
    end
    bound = Mycelia._bind_metamdbg_submission_job!(reservation, "8317")
    fixture_root = dirname(normalized_outdir)
    ready_path = joinpath(fixture_root, "metamdbg-$(phase)-runtime.ready")
    fake_conda = joinpath(fixture_root, "metamdbg-$(phase)-fake-conda")
    write(
        fake_conda,
        "#!/usr/bin/env bash\n" *
        "set -euo pipefail\n" *
        "[ \"\$1\" = run ] || exit 90\n" *
        "shift\n" *
        "while [ \"\$#\" -gt 0 ] && [ \"\$1\" != python ]; do shift; done\n" *
        "[ \"\$#\" -gt 0 ] || exit 91\n" *
        "shift\n" *
        "exec python3 \"\$@\"\n",
    )
    chmod(fake_conda, 0o700)
    pause_hook = "touch $(Base.shell_escape(ready_path)); kill -STOP \$\$"
    script = Mycelia._metamdbg_executor_script(
        "false",
        "false",
        outputs,
        21,
        input_contract;
        conda_runner = fake_conda,
        submission_reservation = bound,
        post_runtime_marker_publication_hook =
            phase == :post_runtime_marker ? pause_hook : nothing,
        post_runtime_private_rename_hook =
            phase == :post_private_rename ? pause_hook : nothing,
        post_queued_marker_release_hook =
            phase == :post_queued_release ? pause_hook : nothing,
        post_runtime_marker_release_hook =
            phase == :post_runtime_release ? pause_hook : nothing,
        post_consumed_owner_rename_hook =
            phase == :post_consumed_rename ? pause_hook : nothing,
        lock_retry_attempts = 1,
        lock_retry_delay_seconds = 0.0,
    )
    script_path = joinpath(fixture_root, "metamdbg-$(phase)-runtime.sh")
    write(script_path, script)
    process = run(
        addenv(`bash $(script_path)`, "SLURM_JOB_ID" => bound.job_id);
        wait = false,
    )
    wait_status = Base.timedwait(
        () -> isfile(ready_path) || !Base.process_running(process),
        120.0;
        pollint = 0.05,
    )
    if wait_status != :ok || !isfile(ready_path)
        Base.process_running(process) && Base.kill(process, Base.SIGKILL)
        wait(process)
        error(
            "metaMDBG $(phase) runtime fixture did not reach its transition " *
            "hook (exit=$(process.exitcode), signal=$(process.termsignal)).",
        )
    end
    try
        live_check(process, bound)
    finally
        Base.process_running(process) && Base.kill(process, Base.SIGKILL)
        wait(process)
    end
    Test.@test !success(process)
    Test.@test process.termsignal == 9
    return bound
end

function _kill_metamdbg_reclaim_transition!(
        outdir::AbstractString,
        reads::AbstractString,
        phase::Symbol,
        ;
        resume_reclaiming::Bool = false,
        live_check::Function = (
            _process::Base.Process,
            _reservation::NamedTuple,
        ) -> nothing,
)::NamedTuple
    phase in (:post_private_rename, :post_shared_release) || throw(
        ArgumentError(
            "phase must be :post_private_rename or :post_shared_release.",
        ),
    )
    active_project = Base.active_project()
    active_project isa AbstractString || error(
        "A project is required for the metaMDBG reclaim SIGKILL fixture.",
    )
    normalized_outdir = String(outdir)
    outputs = Mycelia._metamdbg_output_paths(normalized_outdir, 21)
    input_contract = Mycelia._metamdbg_input_contract(
        Mycelia._metamdbg_selected_input(String(reads), nothing),
        3,
    )
    reservation = Mycelia._metamdbg_submission_reservation(
        outputs,
        input_contract,
        21;
        owner_token = "$(phase)-reclaim-owner",
    )
    if !resume_reclaiming
        Mycelia._with_metamdbg_output_lock(normalized_outdir) do
            Mycelia._create_metamdbg_submission_reservation!(
                reservation,
                normalized_outdir,
            )
        end
    end
    ready_path = joinpath(
        dirname(normalized_outdir),
        "metamdbg-$(phase)-reclaim.ready",
    )
    script_path = joinpath(
        dirname(normalized_outdir),
        "metamdbg-$(phase)-reclaim.jl",
    )
    rm(ready_path; force = true)
    hook_keyword = phase == :post_private_rename ?
                   "post_private_rename_hook" : "post_shared_release_hook"
    script = """
    ENV["LD_LIBRARY_PATH"] = ""
    const Mycelia = Base.require(Main, :Mycelia)

    outdir = $(repr(normalized_outdir))
    ready_path = $(repr(ready_path))
    reservation_path = $(repr(
        resume_reclaiming ? reservation.reclaiming_path : reservation.path,
    ))
    reservation = Mycelia._metamdbg_submission_reservation_from_path(
        reservation_path,
        outdir,
    )
    pause_hook = function (_reservation::NamedTuple)
        write(ready_path, "ready\\n")
        while true
            sleep(60)
        end
    end
    Mycelia._with_metamdbg_output_domain_lock(
        outdir;
        allowed_same_root_locks = (
            reservation.output_root_reservation_marker,
        ),
    ) do
        Mycelia._remove_metamdbg_submission_reservation!(
            reservation;
            $(hook_keyword) = pause_hook,
        )
    end
    """
    write(script_path, script)
    project_directory = dirname(String(active_project))
    process = run(
        `$(Base.julia_cmd()) --project=$(project_directory) --startup-file=no --compiled-modules=yes $(script_path)`;
        wait = false,
    )
    wait_status = Base.timedwait(
        () -> isfile(ready_path) || !Base.process_running(process),
        600.0;
        pollint = 0.05,
    )
    if wait_status != :ok || !isfile(ready_path)
        Base.process_running(process) && Base.kill(process, Base.SIGKILL)
        wait(process)
        error(
            "metaMDBG $(phase) reclaim fixture did not reach its transition " *
            "hook (exit=$(process.exitcode), signal=$(process.termsignal)).",
        )
    end
    try
        live_check(process, reservation)
    finally
        Base.process_running(process) && Base.kill(process, Base.SIGKILL)
        wait(process)
    end
    Test.@test !success(process)
    Test.@test process.termsignal == 9
    return reservation
end

Test.@testset "metaMDBG input and artifact contracts" begin
    mktempdir() do temporary_root
        temporary_root = realpath(temporary_root)
        valid_reads = joinpath(temporary_root, "reads.fastq")
        empty_reads = joinpath(temporary_root, "empty.fastq")
        near_name_max_reads = joinpath(
            temporary_root,
            repeat("n", 244) * ".fastq",
        )
        backslash_reads = joinpath(
            temporary_root,
            "reads\\backslash.fastq",
        )
        write(valid_reads, "@read-1\nACGT\n+\nIIII\n")
        touch(empty_reads)
        write(near_name_max_reads, "@long-name\nACGT\n+\nIIII\n")
        write(backslash_reads, "@backslash\nACGT\n+\nIIII\n")
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
        unattested_ont_outdir =
            joinpath(temporary_root, "unattested-ont-input")
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                ont_reads = valid_reads,
                outdir = unattested_ont_outdir,
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ArgumentError,
            r"ont_r10_4_plus=true.*generic, R9, and unknown ONT",
        )
        Test.@test !ispath(unattested_ont_outdir)
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                ont_r10_4_plus = true,
                outdir = joinpath(temporary_root, "hifi-ont-attestation"),
                dependency_checker,
                local_runner = forbidden_runner,
            ),
            ArgumentError,
            r"applies only when ont_reads is selected",
        )
        for (executor_name, executor) in (
                (:local, Mycelia.LocalExecutor()),
                (:collect, Mycelia.CollectExecutor()),
                (:slurm, Mycelia.SlurmExecutor(dry_run = false)),
        )
            for graph_k in (0, -1)
                invalid_graph_outdir = joinpath(
                    temporary_root,
                    "invalid-graph-$(executor_name)-$(graph_k)",
                )
                side_effect_calls = Ref(0)
                _test_metamdbg_error(
                    () -> Mycelia._run_metamdbg(;
                        hifi_reads = valid_reads,
                        outdir = invalid_graph_outdir,
                        graph_k,
                        executor,
                        dependency_checker = () -> begin
                            side_effect_calls[] += 1
                            return _test_metamdbg_toolchain()
                        end,
                        local_runner = forbidden_runner,
                        submission_runner = (
                            _job::Mycelia.JobSpec,
                            _executor::Mycelia.AbstractExecutor,
                        ) -> begin
                            side_effect_calls[] += 1
                            return nothing
                        end,
                    ),
                    ArgumentError,
                    r"graph_k must be positive",
                )
                Test.@test side_effect_calls[] == 0
                Test.@test !ispath(invalid_graph_outdir)
                if executor isa Mycelia.CollectExecutor
                    Test.@test isempty(executor.jobs)
                end
            end
        end
        _test_metamdbg_error(
            () -> Mycelia._run_metamdbg(;
                hifi_reads = valid_reads,
                ont_reads = valid_reads,
                ont_r10_4_plus = true,
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

        Test.@testset "local staged-input capabilities are identity-bound" begin
            staged = Mycelia._stage_metamdbg_inputs!(
                ordered_input,
                ordered_contract,
                temporary_root,
            )
            Test.@test staged.root_identity ==
                       Mycelia._metamdbg_directory_path_identity(staged.root)
            Test.@test (stat(staged.root).mode & 0o777) == 0o700
            Test.@test isopen(staged.root_descriptor)
            Test.@test read.(staged.selected_input.paths, String) ==
                       read.(ordered_paths, String)
            Test.@test all(
                staged_path ->
                    (stat(staged_path).mode & 0o777) == 0o400,
                staged.selected_input.paths,
            )
            staged_root = staged.root
            Mycelia._cleanup_staged_metamdbg_inputs!(staged)
            Test.@test !ispath(staged_root)
            Test.@test !isopen(staged.root_descriptor)

            injected_copy_root = mktempdir(
                temporary_root;
                prefix = "injected-copy-root-",
                cleanup = false,
            )
            injected_copy_identity =
                Mycelia._metamdbg_directory_path_identity(injected_copy_root)
            injected_copy_directory = Base.Filesystem.open(
                injected_copy_root,
                Base.JL_O_RDONLY |
                Base.JL_O_DIRECTORY |
                Base.JL_O_NOFOLLOW |
                Base.JL_O_CLOEXEC,
            )
            injected_copy_primary = ErrorException(
                "synthetic staged-copy primary failure",
            )
            injected_copy_interrupt = InterruptException()
            injected_copy_close_calls = Ref(0)
            observed_injected_copy_cleanup = Test.@test_logs (
                :warn,
                r"additional cleanup failure",
            ) min_level=Logging.Warn match_mode=:any begin
                try
                    Mycelia._copy_metamdbg_input!(
                        first(ordered_paths),
                        joinpath(injected_copy_root, "copied.fastq"),
                        injected_copy_directory,
                        injected_copy_identity;
                        post_copy_hook = (_source, _destination) ->
                            throw(injected_copy_primary),
                        stream_closer = stream -> begin
                            injected_copy_close_calls[] += 1
                            close(stream)
                            injected_copy_close_calls[] == 1 && error(
                                "synthetic staged-copy ordinary close failure",
                            )
                            throw(injected_copy_interrupt)
                        end,
                    )
                    nothing
                catch caught
                    caught
                end
            end
            Test.@test observed_injected_copy_cleanup ===
                       injected_copy_interrupt
            Test.@test injected_copy_close_calls[] == 2
            close(injected_copy_directory)
            rm(injected_copy_root; recursive = true, force = true)

            injected_root_stage = Mycelia._stage_metamdbg_inputs!(
                ordered_input,
                ordered_contract,
                temporary_root,
            )
            injected_root_primary = ErrorException(
                "synthetic staged-root primary failure",
            )
            observed_injected_root_primary = Test.@test_logs (
                :warn,
                r"root descriptor cleanup failed while preserving",
            ) min_level=Logging.Warn match_mode=:any begin
                try
                    Mycelia._cleanup_staged_metamdbg_inputs!(
                        injected_root_stage;
                        pre_remove_hook = (_root, _identity) ->
                            throw(injected_root_primary),
                        descriptor_closer = descriptor -> begin
                            close(descriptor)
                            error("synthetic staged-root close failure")
                        end,
                    )
                    nothing
                catch caught
                    caught
                end
            end
            Test.@test observed_injected_root_primary === injected_root_primary
            Test.@test !isopen(injected_root_stage.root_descriptor)
            Test.@test isdir(injected_root_stage.root)
            rm(injected_root_stage.root; recursive = true, force = true)

            symlink_target = joinpath(
                temporary_root,
                "staged-input-symlink-target.fastq",
            )
            write(symlink_target, "preserve symlink target\n")
            symlink_root = Ref("")
            symlink_error = try
                Mycelia._stage_metamdbg_inputs!(
                    ordered_input,
                    ordered_contract,
                    temporary_root;
                    pre_child_open_hook = function (
                            root::NamedTuple,
                            destination::AbstractString,
                            input_index::Int,
                    )
                        if input_index == 1
                            symlink_root[] = root.path
                            symlink(symlink_target, destination)
                        end
                        return nothing
                    end,
                )
                nothing
            catch caught
                caught
            end
            Test.@test symlink_error isa SystemError
            Test.@test occursin(
                "openat",
                sprint(showerror, symlink_error),
            )
            Test.@test read(symlink_target, String) ==
                       "preserve symlink target\n"
            Test.@test !ispath(symlink_root[])

            swapped_root = Ref("")
            held_copy_root = joinpath(
                temporary_root,
                "held-staged-input-copy-root",
            )
            replacement_copy_marker = Ref("")
            root_swap_error = Test.@test_logs (
                :warn,
                r"retained staged-input cleanup evidence",
            ) min_level=Logging.Warn begin
                try
                    Mycelia._stage_metamdbg_inputs!(
                        ordered_input,
                        ordered_contract,
                        temporary_root;
                        pre_child_open_hook = function (
                                root::NamedTuple,
                                _destination::AbstractString,
                                input_index::Int,
                        )
                            input_index == 1 || return nothing
                            swapped_root[] = root.path
                            mv(root.path, held_copy_root)
                            mkpath(root.path)
                            replacement_copy_marker[] = joinpath(
                                root.path,
                                "replacement-marker.txt",
                            )
                            write(
                                replacement_copy_marker[],
                                "preserve copy replacement\n",
                            )
                            return nothing
                        end,
                    )
                    nothing
                catch caught
                    caught
                end
            end
            Test.@test root_swap_error isa ErrorException
            Test.@test occursin(
                "staged-input root changed before child creation",
                sprint(showerror, root_swap_error),
            )
            Test.@test read(replacement_copy_marker[], String) ==
                       "preserve copy replacement\n"
            Test.@test isdir(held_copy_root)
            rm(swapped_root[]; recursive = true, force = true)
            rm(held_copy_root; recursive = true, force = true)

            cleanup_swap_stage = Mycelia._stage_metamdbg_inputs!(
                ordered_input,
                ordered_contract,
                temporary_root,
            )
            held_cleanup_root = joinpath(
                temporary_root,
                "held-staged-input-cleanup-root",
            )
            replacement_cleanup_marker = joinpath(
                cleanup_swap_stage.root,
                "replacement-marker.txt",
            )
            cleanup_swap_error = try
                Mycelia._cleanup_staged_metamdbg_inputs!(
                    cleanup_swap_stage;
                    pre_remove_hook = function (
                            root::AbstractString,
                            _identity::NamedTuple,
                    )
                        mv(root, held_cleanup_root)
                        mkpath(root)
                        write(
                            replacement_cleanup_marker,
                            "preserve cleanup replacement\n",
                        )
                        return nothing
                    end,
                )
                nothing
            catch caught
                caught
            end
            Test.@test cleanup_swap_error isa ErrorException
            Test.@test occursin(
                "refuses to remove a replacement durable directory",
                sprint(showerror, cleanup_swap_error),
            )
            Test.@test read(replacement_cleanup_marker, String) ==
                       "preserve cleanup replacement\n"
            Test.@test !isopen(cleanup_swap_stage.root_descriptor)
            Test.@test read.(
                joinpath.(
                    held_cleanup_root,
                    basename.(cleanup_swap_stage.selected_input.paths),
                ),
                String,
            ) == read.(ordered_paths, String)
            rm(cleanup_swap_stage.root; recursive = true, force = true)
            rm(held_cleanup_root; recursive = true, force = true)

            staging_interrupt = InterruptException()
            interrupt_staging_root = Ref("")
            held_interrupt_staging_root = joinpath(
                temporary_root,
                "held-interrupted-staging-root",
            )
            interrupt_replacement_marker = Ref("")
            observed_staging_interrupt = Test.@test_logs (
                :warn,
                r"preserving the primary staging failure",
            ) min_level=Logging.Warn match_mode=:any begin
                try
                    Mycelia._stage_metamdbg_inputs!(
                        ordered_input,
                        ordered_contract,
                        temporary_root;
                        pre_child_open_hook = function (
                                root::NamedTuple,
                                _destination::AbstractString,
                                input_index::Int,
                        )
                            input_index == 1 || return nothing
                            interrupt_staging_root[] = root.path
                            mv(root.path, held_interrupt_staging_root)
                            mkpath(root.path)
                            interrupt_replacement_marker[] = joinpath(
                                root.path,
                                "replacement-marker.txt",
                            )
                            write(
                                interrupt_replacement_marker[],
                                "preserve interrupted replacement\n",
                            )
                            throw(staging_interrupt)
                        end,
                    )
                    nothing
                catch caught
                    caught
                end
            end
            Test.@test observed_staging_interrupt === staging_interrupt
            Test.@test read(interrupt_replacement_marker[], String) ==
                       "preserve interrupted replacement\n"
            Test.@test isdir(held_interrupt_staging_root)
            rm(interrupt_staging_root[]; recursive = true, force = true)
            rm(held_interrupt_staging_root; recursive = true, force = true)

            for (case_label, primary_error) in (
                    (
                        :ordinary,
                        ErrorException("synthetic metaMDBG local primary failure"),
                    ),
                    (:interrupt, InterruptException()),
            )
                failure_outdir = joinpath(
                    temporary_root,
                    "staged-cleanup-preserves-$(case_label)-primary",
                )
                staged_root = Ref("")
                held_staged_root = joinpath(
                    temporary_root,
                    "held-local-$(case_label)-staged-root",
                )
                replacement_marker = Ref("")
                local_runner = function (command::Cmd)
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
                    staged_root[] = dirname(first(staged_paths))
                    mv(staged_root[], held_staged_root)
                    mkpath(staged_root[])
                    replacement_marker[] = joinpath(
                        staged_root[],
                        "replacement-marker.txt",
                    )
                    write(
                        replacement_marker[],
                        "preserve local replacement\n",
                    )
                    throw(primary_error)
                end
                observed_error = Test.@test_logs (
                    :warn,
                    r"preserving the primary local workflow failure",
                ) min_level=Logging.Warn match_mode=:any begin
                    try
                        Mycelia._run_metamdbg(;
                            hifi_reads = ordered_paths,
                            outdir = failure_outdir,
                            dependency_checker = _test_metamdbg_toolchain,
                            local_runner,
                        )
                        nothing
                    catch caught
                        caught
                    end
                end
                Test.@test observed_error === primary_error
                Test.@test read(replacement_marker[], String) ==
                           "preserve local replacement\n"
                Test.@test isdir(held_staged_root)
                rm(staged_root[]; recursive = true, force = true)
                rm(held_staged_root; recursive = true, force = true)
            end
        end

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

        Test.@test ncodeunits(basename(near_name_max_reads)) == 250
        expected_long_staged_name =
            Mycelia._metamdbg_staged_input_name(1, near_name_max_reads)
        Test.@test ncodeunits(expected_long_staged_name) <= 255
        Test.@test !occursin(repeat("n", 244), expected_long_staged_name)
        Test.@test endswith(
            Mycelia._metamdbg_staged_input_name(
                2,
                near_name_max_reads * ".gz",
            ),
            ".gz",
        )
        long_input_outdir = joinpath(temporary_root, "long-input-local")
        long_input_command_count = Ref(0)
        long_input_runner = function (command::Cmd)
            long_input_command_count[] += 1
            if long_input_command_count[] == 1
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
                Test.@test length(staged_paths) == 1
                staged_path = only(staged_paths)
                Test.@test basename(staged_path) == expected_long_staged_name
                Test.@test read(staged_path, String) ==
                           read(near_name_max_reads, String)
                write(
                    joinpath(long_input_outdir, "contigs.fasta"),
                    ">long-input-contig\nACGT\n",
                )
            elseif long_input_command_count[] == 2
                _write_test_metamdbg_gfa!(joinpath(
                    long_input_outdir,
                    "assemblyGraph_k21_4bps.gfa",
                ))
            else
                error("unexpected long-input metaMDBG command")
            end
            return nothing
        end
        long_input_result = Mycelia._run_metamdbg(;
            hifi_reads = near_name_max_reads,
            outdir = long_input_outdir,
            dependency_checker = _test_metamdbg_toolchain,
            local_runner = long_input_runner,
        )
        Test.@test long_input_result.status == :complete
        Test.@test long_input_command_count[] == 2
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
            compression_root = joinpath(temporary_root, "compression-cleanup")
            mkpath(compression_root)
            compression_input = joinpath(compression_root, "plain.fasta")
            compression_output = joinpath(compression_root, "contigs.fasta.gz")
            write(compression_input, ">contig-1\nACGT\n")
            compression_primary =
                ErrorException("synthetic compression primary failure")
            captured_temporary = Ref("")
            held_temporary = joinpath(compression_root, "held-temporary")
            replacement_marker = "replacement survives\n"
            observed_compression_primary = Test.@test_logs (
                :warn,
                r"contig-compression cleanup failed",
            ) min_level=Logging.Warn match_mode=:any begin
                try
                    Mycelia._gzip_metamdbg_contigs!(
                        compression_input,
                        compression_output;
                        after_temporary_binding_hook =
                            (temporary_path, _identity) -> begin
                                captured_temporary[] = temporary_path
                                mv(temporary_path, held_temporary)
                                write(temporary_path, replacement_marker)
                                throw(compression_primary)
                            end,
                    )
                    nothing
                catch caught
                    caught
                end
            end
            Test.@test observed_compression_primary === compression_primary
            Test.@test read(captured_temporary[], String) == replacement_marker
            Test.@test isfile(held_temporary)
            Test.@test !ispath(compression_output)
            rm(captured_temporary[]; force = true)
            rm(held_temporary; force = true)

            captured_success_aba = Ref("")
            held_success_aba = joinpath(
                compression_root,
                "held-success-aba-temporary",
            )
            success_aba_replacement = "successful-action replacement\n"
            observed_success_aba = Test.@test_logs (
                :warn,
                r"contig-compression cleanup failed",
            ) min_level=Logging.Warn match_mode=:any begin
                try
                    Mycelia._gzip_metamdbg_contigs!(
                        compression_input,
                        compression_output;
                        after_temporary_binding_hook =
                            (temporary_path, _identity) -> begin
                                captured_success_aba[] = temporary_path
                                mv(temporary_path, held_success_aba)
                                write(
                                    temporary_path,
                                    success_aba_replacement,
                                )
                            end,
                    )
                    nothing
                catch caught
                    caught
                end
            end
            Test.@test observed_success_aba isa ErrorException
            Test.@test occursin(
                "temporary artifact was replaced",
                sprint(showerror, observed_success_aba),
            )
            Test.@test read(captured_success_aba[], String) ==
                       success_aba_replacement
            Test.@test isfile(held_success_aba)
            Test.@test !ispath(compression_output)
            rm(captured_success_aba[]; force = true)
            rm(held_success_aba; force = true)

            post_hardlink_temporary = Ref("")
            post_hardlink_primary = ErrorException(
                "synthetic gzip post-hardlink failure",
            )
            observed_post_hardlink_primary = try
                Mycelia._gzip_metamdbg_contigs!(
                    compression_input,
                    compression_output;
                    after_temporary_binding_hook =
                        (temporary_path, _identity) ->
                            (post_hardlink_temporary[] = temporary_path),
                    post_hardlink_hook = (_path, _identity) ->
                        throw(post_hardlink_primary),
                )
                nothing
            catch caught
                caught
            end
            Test.@test observed_post_hardlink_primary === post_hardlink_primary
            Test.@test !ispath(post_hardlink_temporary[])
            Test.@test !ispath(compression_output)

            in_place_mutation_temporary = Ref("")
            observed_in_place_mutation = try
                Mycelia._gzip_metamdbg_contigs!(
                    compression_input,
                    compression_output;
                    after_temporary_binding_hook =
                        (temporary_path, _identity) ->
                            (in_place_mutation_temporary[] = temporary_path),
                    post_hardlink_hook = (path, _identity) -> begin
                        raw_output = open(path, "w")
                        replacement_output =
                            CodecZlib.GzipCompressorStream(raw_output)
                        write(
                            replacement_output,
                            ">replacement-contig\nTGCA\n",
                        )
                        close(replacement_output)
                        return nothing
                    end,
                )
                nothing
            catch caught
                caught
            end
            Test.@test observed_in_place_mutation isa ErrorException
            Test.@test occursin(
                "content digest changed after publication",
                sprint(showerror, observed_in_place_mutation),
            )
            Test.@test !ispath(in_place_mutation_temporary[])
            Test.@test !ispath(compression_output)

            fasta_cleanup_failure = ErrorException(
                "synthetic metaMDBG FASTA reader close failure",
            )
            observed_fasta_cleanup = try
                Mycelia._require_valid_metamdbg_fasta(
                    compression_input,
                    "close-injected metaMDBG contigs";
                    reader_closer = reader -> begin
                        close(reader)
                        throw(fasta_cleanup_failure)
                    end,
                )
                nothing
            catch caught
                caught
            end
            Test.@test observed_fasta_cleanup === fasta_cleanup_failure

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
                ont_r10_4_plus = true,
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
                Mycelia._metamdbg_selected_input(
                    nothing,
                    valid_reads,
                    true,
                ),
                3,
            )
            Test.@test isfile(outputs.contract_marker)
            Test.@test read(outputs.contract_marker, String) ==
                       expected_contract.contents
            expected_input = only(expected_contract.contract.inputs)
            Test.@test expected_contract.contract.schema_version == 5
            Test.@test expected_contract.contract.platform_attestation ==
                       "nanopore-r10.4-or-later"
            Test.@test expected_input.sha256 ==
                       Mycelia._metamdbg_sha256(valid_reads)
            Test.@test !hasproperty(expected_input, :modification_time_ns)
            Test.@test !occursin(
                "modification_time_ns",
                expected_contract.contents,
            )
            Test.@test result.provenance.contract_signature ==
                       expected_contract.signature
            Test.@test result.provenance.input_platform_attestation ==
                       "nanopore-r10.4-or-later"
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
            Test.@test result.provenance.package_inventory ==
                       _test_metamdbg_toolchain().package_inventory
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
            Test.@test completion_record["schema_version"] == 2
            Test.@test completion_record["manifest"]["schema_version"] == 2
            Test.@test completion_record["manifest"]["workflow"]["graph_k"] ==
                       21
            recorded_package_inventory = completion_record["manifest"][
                "toolchain"
            ]["package_inventory"]
            Test.@test recorded_package_inventory == JSON.parse(JSON.json(
                _test_metamdbg_toolchain().package_inventory,
            ))
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
                ont_r10_4_plus = true,
                outdir,
                graph_k = 21,
                dependency_checker = () -> error(
                    "complete metaMDBG reuse must not provision",
                ),
                local_runner = forbidden_runner,
            )
            Test.@test reused == result
            Test.@test command_count[] == 2

            completion_backup = read(outputs.completion_marker, String)
            mismatched_inventory_record = deepcopy(completion_record)
            mismatched_inventory_record["manifest"]["toolchain"][
                "package_inventory"
            ][1]["build"] = "synthetic-mismatched-build"
            write(
                outputs.completion_marker,
                JSON.json(mismatched_inventory_record) * "\n",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    ont_r10_4_plus = true,
                    outdir,
                    graph_k = 21,
                    dependency_checker = () -> error(
                        "inventory-mismatched reuse must fail before provisioning",
                    ),
                    local_runner = forbidden_runner,
                ),
                ErrorException,
                r"package inventory does not match its normalized digest",
            )
            write(outputs.completion_marker, completion_backup)

            contigs_backup = joinpath(temporary_root, "contigs-backup.fasta.gz")
            cp(result.contigs, contigs_backup)
            _write_test_metamdbg_gzip_fasta!(
                result.contigs;
                contents = ">contig-1\nTGCA\n",
            )
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    ont_reads = valid_reads,
                    ont_r10_4_plus = true,
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
                    ont_r10_4_plus = true,
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
                    ont_r10_4_plus = true,
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
                    ont_r10_4_plus = true,
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
                Mycelia._metamdbg_selected_input(nothing, valid_reads, true),
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
                    ont_r10_4_plus = true,
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
                Mycelia._metamdbg_selected_input(nothing, valid_reads, true),
                3,
            )
            Test.@test restored_contract.contents == expected_contract.contents

            run(`touch -t 203001010101 $(valid_reads)`)
            touched_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(nothing, valid_reads, true),
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
                    ont_r10_4_plus = true,
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

            for mutation_kind in (:contigs, :graph)
                artifact_reads = joinpath(
                    temporary_root,
                    "local-artifact-$(mutation_kind).fastq",
                )
                write(
                    artifact_reads,
                    "@local-artifact-$(mutation_kind)\nACGT\n+\nIIII\n",
                )
                artifact_outdir = joinpath(
                    temporary_root,
                    "local-artifact-$(mutation_kind)-output",
                )
                artifact_outputs =
                    Mycelia._metamdbg_output_paths(artifact_outdir, 21)
                artifact_original_contigs =
                    artifact_outputs.contigs_gz * ".original"
                artifact_dependency_calls = Ref(0)
                artifact_commands = Ref(0)
                _test_metamdbg_error(
                    () -> Mycelia._run_metamdbg(;
                        hifi_reads = artifact_reads,
                        outdir = artifact_outdir,
                        dependency_checker = () -> begin
                            artifact_dependency_calls[] += 1
                            if artifact_dependency_calls[] == 2
                                if mutation_kind == :contigs
                                    # Keep the original inode linked so Linux
                                    # cannot recycle it for the replacement
                                    # fixture before cleanup's identity fence.
                                    mv(
                                        artifact_outputs.contigs_gz,
                                        artifact_original_contigs,
                                    )
                                    _write_test_metamdbg_gzip_fasta!(
                                        artifact_outputs.contigs_gz;
                                        contents = ">replacement\nTGCA\n",
                                    )
                                else
                                    _write_test_metamdbg_gfa!(joinpath(
                                        artifact_outdir,
                                        "assemblyGraph_k21_4bps.gfa",
                                    ); sequence = "TGCA")
                                end
                            end
                            return _test_metamdbg_toolchain()
                        end,
                        local_runner = function (_command::Cmd)
                            artifact_commands[] += 1
                            if artifact_commands[] == 1
                                write(
                                    artifact_outputs.contigs_plain,
                                    ">contig-1\nACGT\n",
                                )
                            else
                                _write_test_metamdbg_gfa!(joinpath(
                                    artifact_outdir,
                                    "assemblyGraph_k21_4bps.gfa",
                                ))
                            end
                            return nothing
                        end,
                    ),
                    ErrorException,
                    r"artifacts changed after post-execution toolchain validation",
                )
                Test.@test artifact_dependency_calls[] == 2
                Test.@test artifact_commands[] == 2
                Test.@test !ispath(artifact_outputs.graph_alias)
                Test.@test !ispath(artifact_outputs.contract_marker)
                Test.@test !ispath(artifact_outputs.completion_marker)
                Test.@test isfile(artifact_outputs.contigs_plain)
                Test.@test isfile(joinpath(
                    artifact_outdir,
                    "assemblyGraph_k21_4bps.gfa",
                ))
                Test.@test ispath(artifact_outputs.contigs_gz) ==
                           (mutation_kind == :contigs)

                if mutation_kind == :contigs
                    Test.@test isfile(artifact_original_contigs)
                    Test.@test isfile(artifact_outputs.contigs_gz)
                    Test.@test !Base.Filesystem.samefile(
                        artifact_original_contigs,
                        artifact_outputs.contigs_gz,
                    )
                    Test.@test String(Base.transcode(
                        CodecZlib.GzipDecompressor,
                        read(artifact_original_contigs),
                    )) == ">contig-1\nACGT\n"
                    Test.@test String(Base.transcode(
                        CodecZlib.GzipDecompressor,
                        read(artifact_outputs.contigs_gz),
                    )) == ">replacement\nTGCA\n"
                end

                if mutation_kind == :graph
                    _test_metamdbg_error(
                        () -> Mycelia._run_metamdbg(;
                            hifi_reads = artifact_reads,
                            outdir = artifact_outdir,
                            dependency_checker =
                                _test_metamdbg_toolchain,
                            local_runner = forbidden_runner,
                        ),
                        ErrorException,
                        r"nonempty output root without its regular provenance contract",
                    )
                    rm(artifact_outputs.contigs_plain)
                    rm(joinpath(
                        artifact_outdir,
                        "assemblyGraph_k21_4bps.gfa",
                    ))
                    recovery_commands = Ref(0)
                    recovered = Mycelia._run_metamdbg(;
                        hifi_reads = artifact_reads,
                        outdir = artifact_outdir,
                        dependency_checker = _test_metamdbg_toolchain,
                        local_runner = function (_command::Cmd)
                            recovery_commands[] += 1
                            if recovery_commands[] == 1
                                write(
                                    artifact_outputs.contigs_plain,
                                    ">recovered\nACGT\n",
                                )
                            else
                                _write_test_metamdbg_gfa!(joinpath(
                                    artifact_outdir,
                                    "assemblyGraph_k21_4bps.gfa",
                                ))
                            end
                            return nothing
                        end,
                    )
                    Test.@test recovered.status == :complete
                    Test.@test recovery_commands[] == 2
                end
            end

            final_reads =
                joinpath(temporary_root, "local-final-input.fastq")
            write(final_reads, "@local-final-input\nACGT\n+\nIIII\n")
            final_outdir =
                joinpath(temporary_root, "local-final-input-output")
            final_outputs =
                Mycelia._metamdbg_output_paths(final_outdir, 21)
            final_digest_calls = Ref(0)
            final_commands = Ref(0)
            _test_metamdbg_error(
                () -> Mycelia._run_metamdbg(;
                    hifi_reads = final_reads,
                    outdir = final_outdir,
                    dependency_checker = _test_metamdbg_toolchain,
                    input_digest_function = function (
                            path::AbstractString,
                    )
                        final_digest_calls[] += 1
                        if final_digest_calls[] == 3
                            write(
                                path,
                                "@local-final-input\nTGCA\n+\nIIII\n",
                            )
                        end
                        return Mycelia._metamdbg_sha256(path)
                    end,
                    local_runner = function (_command::Cmd)
                        final_commands[] += 1
                        if final_commands[] == 1
                            write(
                                final_outputs.contigs_plain,
                                ">contig-1\nACGT\n",
                            )
                        else
                            _write_test_metamdbg_gfa!(joinpath(
                                final_outdir,
                                "assemblyGraph_k21_4bps.gfa",
                            ))
                        end
                        return nothing
                    end,
                ),
                ErrorException,
                r"input content contract changed",
            )
            Test.@test final_digest_calls[] == 3
            Test.@test final_commands[] == 2
            Test.@test !ispath(final_outputs.contigs_gz)
            Test.@test !ispath(final_outputs.graph_alias)
            Test.@test !ispath(final_outputs.contract_marker)
            Test.@test !ispath(final_outputs.completion_marker)
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
            Test.@test local_digest_calls ==
                       String[local_reads, local_reads, local_reads]

            local_outputs =
                Mycelia._metamdbg_output_paths(local_outdir, 21)
            local_artifacts = (;
                outdir = local_outdir,
                contigs = local_result.contigs,
                graph = local_result.graph,
            )
            artifact_snapshot_digest_calls = Ref(0)
            artifact_snapshot =
                Mycelia._validated_metamdbg_artifact_snapshot(
                    local_artifacts;
                    digest_function = function (input::IO)
                        artifact_snapshot_digest_calls[] += 1
                        return Mycelia._metamdbg_descriptor_sha256(input)
                    end,
                )
            Test.@test artifact_snapshot_digest_calls[] == 4
            unchanged_artifact_digest_calls = Ref(0)
            Mycelia._require_unchanged_metamdbg_artifacts!(
                local_artifacts,
                artifact_snapshot,
                "during bounded test";
                digest_function = function (input::IO)
                    unchanged_artifact_digest_calls[] += 1
                    return Mycelia._metamdbg_descriptor_sha256(input)
                end,
            )
            Test.@test unchanged_artifact_digest_calls[] <= 2

            local_selected_input =
                Mycelia._metamdbg_selected_input(local_reads, nothing)
            local_input_contract = Mycelia._metamdbg_input_contract(
                local_selected_input,
                3,
            )
            fresh_completion_digest_calls = Ref(0)
            fresh_completion_digest_maximum = 0
            fresh_completion = Mycelia._metamdbg_completion_manifest(
                local_outputs,
                local_artifacts,
                local_input_contract,
                21,
                _test_metamdbg_toolchain();
                artifact_snapshot,
                digest_function = function (path::AbstractString)
                    fresh_completion_digest_calls[] += 1
                    return Mycelia._metamdbg_sha256(path)
                end,
            )
            Test.@test fresh_completion_digest_calls[] <=
                       fresh_completion_digest_maximum
            Test.@test fresh_completion.contents ==
                       read(local_outputs.completion_marker, String)

            reuse_completion_digest_calls = Ref(0)
            reuse_completion_digest_maximum = 2
            reused_completion =
                Mycelia._require_metamdbg_completion_manifest!(
                    local_outputs,
                    local_artifacts,
                    local_input_contract,
                    21;
                    digest_function = function (path::AbstractString)
                        reuse_completion_digest_calls[] += 1
                        return Mycelia._metamdbg_sha256(path)
                    end,
                )
            Test.@test reuse_completion_digest_calls[] <=
                       reuse_completion_digest_maximum
            Test.@test reused_completion.contents == fresh_completion.contents

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
            missing_metamdbg_inventory = NamedTuple[(;
                name = "libzlib",
                version = "1.3.1",
                build = "h8359307_2",
                channel = "conda-forge",
            )]
            missing_metamdbg_completion_toolchain = (;
                Mycelia._metamdbg_expected_toolchain()...,
                package_inventory = missing_metamdbg_inventory,
                package_inventory_sha256 =
                    Mycelia._metamdbg_package_inventory_sha256(
                        missing_metamdbg_inventory,
                    ),
                package_count = 1,
            )
            _test_metamdbg_error(
                () -> Mycelia._metamdbg_completion_toolchain_summary(
                    missing_metamdbg_completion_toolchain,
                ),
                ErrorException,
                r"must contain metamdbg exactly 1.4",
            )

            escaped_inventory = Any[
                Dict(
                    "channel" => "bio\"conda\\edge",
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
            inventory_manifest = joinpath(
                temporary_root,
                "escaped-reordered-inventory.manifest.json",
            )
            serialized_inventory = JSON.json(reverse(escaped_inventory))
            write(
                inventory_json,
                "  \n" * replace(serialized_inventory, "},{" => "},\n {") *
                "\n",
            )
            inventory_canonicalizer =
                Mycelia._metamdbg_runtime_inventory_canonicalizer_python()
            run(Cmd(String[
                "python3",
                "-c",
                inventory_canonicalizer,
                inventory_json,
                inventory_tsv,
                inventory_manifest,
            ]))
            Test.@test read(inventory_tsv, String) ==
                       Mycelia._metamdbg_package_inventory_contents(
                escaped_inventory,
            )
            normalized_escaped_inventory =
                Mycelia._metamdbg_normalized_package_inventory(
                    escaped_inventory,
                )
            Test.@test read(inventory_manifest, String) ==
                       JSON.json(normalized_escaped_inventory)
            Test.@test JSON.parse(read(inventory_manifest, String)) ==
                       JSON.parse(JSON.json(normalized_escaped_inventory))
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
            changed_inventory_manifest = joinpath(
                temporary_root,
                "changed-build-inventory.manifest.json",
            )
            write(changed_inventory_json, JSON.json(escaped_changed))
            run(Cmd(String[
                "python3",
                "-c",
                inventory_canonicalizer,
                changed_inventory_json,
                changed_inventory_tsv,
                changed_inventory_manifest,
            ]))
            Test.@test Mycelia._metamdbg_sha256(changed_inventory_tsv) !=
                       Mycelia._metamdbg_sha256(inventory_tsv)
            Test.@test read(changed_inventory_manifest, String) !=
                       read(inventory_manifest, String)
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
                    ont_r10_4_plus = true,
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
                    ont_r10_4_plus = true,
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

            no_pid_recovery_outdir =
                joinpath(temporary_root, "no-pid-recovery-output")
            _test_metamdbg_error(
                () -> Mycelia.inspect_metamdbg_submission_reservations(
                    no_pid_recovery_outdir;
                    confirm_process_dead = true,
                ),
                ErrorException,
                r"requires a pre-existing canonical local PID.*no-op recovery",
            )

            pre_kill_outdir =
                joinpath(temporary_root, "pre-rename-sigkill-output")
            pre_kill_ready =
                joinpath(temporary_root, "pre-rename-sigkill.ready")
            pre_kill_live_check = function (process::Base.Process)
                Test.@test Base.process_running(process)
                private_lock =
                    Mycelia._metamdbg_output_lock_path(pre_kill_outdir)
                cleanup_reservation =
                    Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                        pre_kill_outdir,
                    )
                pid_lock =
                    Mycelia._output_root_reservation_lock_path_from_canonical(
                        Mycelia._metamdbg_canonical_output_path(
                            pre_kill_outdir,
                        ),
                    )
                provisional_path = only(
                    Mycelia._metamdbg_submission_reservation_paths(
                        pre_kill_outdir,
                    ),
                )
                provisional_record =
                    Mycelia._metamdbg_submission_reservation_from_path(
                        provisional_path,
                        pre_kill_outdir;
                        allow_provisional = true,
                    )
                tracked_paths = String[
                    private_lock,
                    cleanup_reservation,
                    pid_lock,
                    provisional_path,
                    provisional_record.output_root_reservation_marker,
                ]
                tracked_identities = map(tracked_paths) do path
                    path_status = stat(path)
                    return (;
                        device = path_status.device,
                        inode = path_status.inode,
                    )
                end
                _test_metamdbg_error(
                    () -> Mycelia.inspect_metamdbg_submission_reservations(
                        pre_kill_outdir;
                        confirm_process_dead = true,
                    ),
                    ErrorException,
                    r"PID lock still names a live or remotely unverifiable process",
                )
                Test.@test all(ispath, tracked_paths)
                observed_identities = map(tracked_paths) do path
                    path_status = stat(path)
                    return (;
                        device = path_status.device,
                        inode = path_status.inode,
                    )
                end
                Test.@test observed_identities == tracked_identities
                return nothing
            end
            _kill_metamdbg_reservation_publication!(
                pre_kill_outdir,
                pre_kill_ready,
                :pre_rename,
                ;
                live_check = pre_kill_live_check,
            )
            pre_kill_private_lock =
                Mycelia._metamdbg_output_lock_path(pre_kill_outdir)
            pre_kill_cleanup_reservation =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    pre_kill_outdir,
                )
            pre_kill_pid_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(pre_kill_outdir),
                )
            Test.@test isdir(pre_kill_private_lock)
            Test.@test isdir(pre_kill_cleanup_reservation)
            Test.@test isfile(pre_kill_pid_lock)
            pre_kill_locked_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    pre_kill_outdir,
                ),
            )
            Test.@test pre_kill_locked_metadata.publication_state ==
                       :provisional
            Test.@test pre_kill_locked_metadata.private_lock_identity !==
                       nothing
            Test.@test pre_kill_locked_metadata.cleanup_reservation_identity !==
                       nothing
            canonical_pre_kill_pid_contents =
                read(pre_kill_pid_lock, String)
            pre_kill_pid_mode = stat(pre_kill_pid_lock).mode & 0o777
            chmod(pre_kill_pid_lock, 0o600)
            write(pre_kill_pid_lock, "")
            chmod(pre_kill_pid_lock, pre_kill_pid_mode)
            _test_metamdbg_error(
                () -> Mycelia.inspect_metamdbg_submission_reservations(
                    pre_kill_outdir;
                    confirm_process_dead = true,
                ),
                ErrorException,
                r"PID lock is empty or malformed",
            )
            Test.@test isdir(pre_kill_private_lock)
            Test.@test isdir(pre_kill_cleanup_reservation)
            Test.@test isfile(pre_kill_pid_lock)
            chmod(pre_kill_pid_lock, 0o600)
            write(pre_kill_pid_lock, "not-a-canonical-pidfile")
            chmod(pre_kill_pid_lock, pre_kill_pid_mode)
            _test_metamdbg_error(
                () -> Mycelia.inspect_metamdbg_submission_reservations(
                    pre_kill_outdir;
                    confirm_process_dead = true,
                ),
                ErrorException,
                r"PID lock is empty or malformed",
            )
            Test.@test isdir(pre_kill_private_lock)
            Test.@test isdir(pre_kill_cleanup_reservation)
            Test.@test isfile(pre_kill_pid_lock)
            chmod(pre_kill_pid_lock, 0o600)
            write(pre_kill_pid_lock, canonical_pre_kill_pid_contents)
            chmod(pre_kill_pid_lock, pre_kill_pid_mode)
            pre_kill_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    pre_kill_outdir,
                    ;
                    confirm_process_dead = true,
                ),
            )
            Test.@test !ispath(pre_kill_private_lock)
            Test.@test !ispath(pre_kill_cleanup_reservation)
            Test.@test !ispath(pre_kill_pid_lock)
            Test.@test pre_kill_metadata.publication_state == :provisional
            Test.@test pre_kill_metadata.submission_state == :reserved
            Test.@test pre_kill_metadata.job_id === nothing
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    pre_kill_metadata;
                    owner_token = pre_kill_metadata.owner_token,
                    job_id = "101",
                    confirm_submitted = true,
                ),
                ArgumentError,
                r"fully published submission reservation",
            )
            Test.@test Mycelia._metamdbg_submission_reservation_path_state(
                pre_kill_metadata.path,
                pre_kill_outdir,
            ) == :provisional
            pre_kill_record =
                Mycelia._metamdbg_submission_reservation_from_path(
                    pre_kill_metadata.path,
                    pre_kill_outdir;
                    allow_provisional = true,
                )
            Test.@test isfile(
                pre_kill_record.output_root_reservation_marker,
            )
            chmod(pre_kill_record.contract_marker, 0o640)
            _test_metamdbg_error(
                () -> Mycelia.inspect_metamdbg_submission_reservations(
                    pre_kill_outdir,
                ),
                ErrorException,
                r"could not validate a provisional.*potentially paired shared",
            )
            chmod(pre_kill_record.contract_marker, 0o600)
            held_directory_marker =
                joinpath(temporary_root, "queued-marker-directory-hold")
            mv(
                pre_kill_record.output_root_reservation_marker,
                held_directory_marker,
            )
            mkdir(pre_kill_record.output_root_reservation_marker)
            chmod(pre_kill_record.contract_marker, 0o640)
            _test_metamdbg_error(
                () -> Mycelia.inspect_metamdbg_submission_reservations(
                    pre_kill_outdir,
                ),
                ErrorException,
                r"could not validate a provisional.*potentially paired shared",
            )
            chmod(pre_kill_record.contract_marker, 0o600)
            rm(pre_kill_record.output_root_reservation_marker)
            mv(
                held_directory_marker,
                pre_kill_record.output_root_reservation_marker,
            )
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    pre_kill_metadata;
                    owner_token = pre_kill_metadata.owner_token,
                    confirm_cancelled = true,
                ),
                ArgumentError,
                r"provisional.*confirm_not_submitted",
            )

            replacement_hold =
                joinpath(temporary_root, "provisional-replacement-hold")
            mkpath(replacement_hold)
            held_private = joinpath(replacement_hold, "private")
            mv(pre_kill_metadata.path, held_private)
            mkpath(pre_kill_metadata.path)
            chmod(pre_kill_metadata.path, 0o700)
            replacement_contract = joinpath(
                pre_kill_metadata.path,
                Mycelia._METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
            )
            write(replacement_contract, pre_kill_record.contents)
            chmod(replacement_contract, 0o600)
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    pre_kill_metadata;
                    owner_token = pre_kill_metadata.owner_token,
                    confirm_not_submitted = true,
                ),
                ErrorException,
                r"reservation was replaced after recovery inspection",
            )
            Test.@test isdir(pre_kill_metadata.path)
            Test.@test isfile(
                pre_kill_record.output_root_reservation_marker,
            )
            rm(pre_kill_metadata.path; recursive = true)
            mv(held_private, pre_kill_metadata.path)

            held_shared = joinpath(replacement_hold, "shared")
            mv(pre_kill_record.output_root_reservation_marker, held_shared)
            write(
                pre_kill_record.output_root_reservation_marker,
                pre_kill_record.output_root_reservation_contents,
            )
            chmod(pre_kill_record.output_root_reservation_marker, 0o600)
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    pre_kill_metadata;
                    owner_token = pre_kill_metadata.owner_token,
                    confirm_not_submitted = true,
                ),
                ErrorException,
                r"shared output-root reservation was replaced after recovery",
            )
            Test.@test isdir(pre_kill_metadata.path)
            Test.@test isfile(
                pre_kill_record.output_root_reservation_marker,
            )
            rm(pre_kill_record.output_root_reservation_marker)
            mv(held_shared, pre_kill_record.output_root_reservation_marker)
            provisional_reclaim =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    pre_kill_metadata;
                    owner_token = pre_kill_metadata.owner_token,
                    confirm_not_submitted = true,
                )
            Test.@test provisional_reclaim.status == :reclaimed
            Test.@test provisional_reclaim.publication_state == :provisional
            Test.@test !ispath(pre_kill_metadata.path)
            Test.@test !ispath(
                pre_kill_record.output_root_reservation_marker,
            )
            Test.@test Mycelia._with_metamdbg_output_domain_lock(
                pre_kill_outdir,
            ) do
                return isempty(
                    Mycelia._metamdbg_submission_reservation_paths(
                        pre_kill_outdir,
                    ),
                )
            end

            post_kill_outdir =
                joinpath(temporary_root, "post-rename-sigkill-output")
            post_kill_ready =
                joinpath(temporary_root, "post-rename-sigkill.ready")
            _kill_metamdbg_reservation_publication!(
                post_kill_outdir,
                post_kill_ready,
                :post_rename,
            )
            post_kill_private_lock =
                Mycelia._metamdbg_output_lock_path(post_kill_outdir)
            post_kill_cleanup_reservation =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    post_kill_outdir,
                )
            post_kill_pid_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(post_kill_outdir),
                )
            Test.@test isdir(post_kill_private_lock)
            Test.@test isdir(post_kill_cleanup_reservation)
            Test.@test isfile(post_kill_pid_lock)
            post_kill_locked_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    post_kill_outdir,
                ),
            )
            Test.@test post_kill_locked_metadata.publication_state ==
                       :published
            Test.@test post_kill_locked_metadata.private_lock_identity !==
                       nothing
            Test.@test post_kill_locked_metadata.cleanup_reservation_identity !==
                       nothing
            post_kill_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    post_kill_outdir,
                    ;
                    confirm_process_dead = true,
                ),
            )
            Test.@test !ispath(post_kill_private_lock)
            Test.@test !ispath(post_kill_cleanup_reservation)
            Test.@test !ispath(post_kill_pid_lock)
            Test.@test post_kill_metadata.publication_state == :published
            Test.@test post_kill_metadata.submission_state == :reserved
            identity_free_post_metadata = (;
                canonical_outdir = post_kill_metadata.canonical_outdir,
                path = post_kill_metadata.path,
                workflow_signature = post_kill_metadata.workflow_signature,
                scheduler_job_name = post_kill_metadata.scheduler_job_name,
                input_contract_signature =
                    post_kill_metadata.input_contract_signature,
                graph_k = post_kill_metadata.graph_k,
                owner_token = post_kill_metadata.owner_token,
                job_id = post_kill_metadata.job_id,
                submission_state = post_kill_metadata.submission_state,
                publication_state = post_kill_metadata.publication_state,
            )
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    identity_free_post_metadata;
                    owner_token = identity_free_post_metadata.owner_token,
                    confirm_not_submitted = true,
                ),
                ArgumentError,
                r"exact private and shared filesystem identities",
            )
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    identity_free_post_metadata;
                    owner_token = identity_free_post_metadata.owner_token,
                    job_id = "102",
                    confirm_submitted = true,
                ),
                ArgumentError,
                r"exact private and shared filesystem identities",
            )
            post_kill_record =
                Mycelia._metamdbg_submission_reservation_from_path(
                    post_kill_metadata.path,
                    post_kill_outdir,
                )
            post_replacement_hold =
                joinpath(temporary_root, "published-replacement-hold")
            mkpath(post_replacement_hold)
            held_post_private = joinpath(post_replacement_hold, "private")
            mv(post_kill_metadata.path, held_post_private)
            mkpath(post_kill_metadata.path)
            chmod(post_kill_metadata.path, 0o700)
            post_replacement_contract = joinpath(
                post_kill_metadata.path,
                Mycelia._METAMDBG_SUBMISSION_RESERVATION_CONTRACT_FILENAME,
            )
            write(post_replacement_contract, post_kill_record.contents)
            chmod(post_replacement_contract, 0o600)
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    post_kill_metadata;
                    owner_token = post_kill_metadata.owner_token,
                    job_id = "103",
                    confirm_submitted = true,
                ),
                ErrorException,
                r"reservation was replaced after recovery inspection",
            )
            Test.@test !ispath(joinpath(
                post_kill_metadata.path,
                Mycelia._METAMDBG_SUBMISSION_RESERVATION_JOB_FILENAME,
            ))
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    post_kill_metadata;
                    owner_token = post_kill_metadata.owner_token,
                    confirm_not_submitted = true,
                ),
                ErrorException,
                r"reservation was replaced after recovery inspection",
            )
            Test.@test isdir(post_kill_metadata.path)
            Test.@test isfile(
                post_kill_record.output_root_reservation_marker,
            )
            rm(post_kill_metadata.path; recursive = true)
            mv(held_post_private, post_kill_metadata.path)

            held_post_shared = joinpath(post_replacement_hold, "shared")
            mv(
                post_kill_record.output_root_reservation_marker,
                held_post_shared,
            )
            write(
                post_kill_record.output_root_reservation_marker,
                post_kill_record.output_root_reservation_contents,
            )
            chmod(post_kill_record.output_root_reservation_marker, 0o600)
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    post_kill_metadata;
                    owner_token = post_kill_metadata.owner_token,
                    confirm_not_submitted = true,
                ),
                ErrorException,
                r"shared output-root reservation was replaced after recovery",
            )
            Test.@test isdir(post_kill_metadata.path)
            Test.@test isfile(
                post_kill_record.output_root_reservation_marker,
            )
            rm(post_kill_record.output_root_reservation_marker)
            mv(
                held_post_shared,
                post_kill_record.output_root_reservation_marker,
            )
            published_reclaim =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    post_kill_metadata;
                    owner_token = post_kill_metadata.owner_token,
                    confirm_not_submitted = true,
                )
            Test.@test published_reclaim.status == :reclaimed
            Test.@test published_reclaim.publication_state == :published
            Test.@test !ispath(post_kill_metadata.path)
            Test.@test !ispath(
                post_kill_record.output_root_reservation_marker,
            )
            Test.@test Mycelia._with_metamdbg_output_domain_lock(
                post_kill_outdir,
            ) do
                return isempty(
                    Mycelia._metamdbg_submission_reservation_paths(
                        post_kill_outdir,
                    ),
                )
            end

            for phase in (:pre_publication, :post_publication)
                bind_cleanup_outdir = joinpath(
                    temporary_root,
                    "$(phase)-job-bind-cleanup-output",
                )
                bind_cleanup_outputs = Mycelia._metamdbg_output_paths(
                    bind_cleanup_outdir,
                    21,
                )
                bind_cleanup_contract = Mycelia._metamdbg_input_contract(
                    Mycelia._metamdbg_selected_input(valid_reads, nothing),
                    3,
                )
                bind_cleanup_reservation =
                    Mycelia._metamdbg_submission_reservation(
                        bind_cleanup_outputs,
                        bind_cleanup_contract,
                        21;
                        owner_token = "$(phase)-bind-cleanup-owner",
                    )
                Mycelia._with_metamdbg_output_lock(
                    bind_cleanup_outdir,
                ) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        bind_cleanup_reservation,
                        bind_cleanup_outdir,
                    )
                end
                parent_entries_before = Set(readdir(temporary_root))
                pre_bind_hook = phase == :pre_publication ?
                                (
                    _reservation::NamedTuple,
                    _path::AbstractString,
                ) -> error("synthetic pre-publication bind failure") :
                                (
                    _reservation::NamedTuple,
                    _path::AbstractString,
                ) -> nothing
                post_bind_hook = phase == :post_publication ?
                                 (_reservation::NamedTuple) -> error(
                    "synthetic post-publication bind failure",
                ) : (_reservation::NamedTuple) -> nothing
                _test_metamdbg_error(
                    () -> Mycelia._with_metamdbg_output_domain_lock(
                        bind_cleanup_outdir;
                        allowed_same_root_locks = (
                            bind_cleanup_reservation.output_root_reservation_marker,
                        ),
                    ) do
                        Mycelia._bind_metamdbg_submission_job!(
                            bind_cleanup_reservation,
                            "770";
                            pre_job_record_publication_hook = pre_bind_hook,
                            post_job_record_publication_hook = post_bind_hook,
                        )
                    end,
                    ErrorException,
                    r"synthetic (pre|post)-publication bind failure",
                )
                bind_cleanup_pending =
                    Mycelia._metamdbg_pending_submission_job_path(
                        bind_cleanup_outdir,
                        bind_cleanup_reservation.output_root_reservation_capability,
                        "770",
                    )
                Test.@test isfile(bind_cleanup_pending) ==
                           (phase == :pre_publication)
                Test.@test isfile(bind_cleanup_reservation.job_marker) ==
                           (phase == :post_publication)
                bind_cleanup_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        bind_cleanup_outdir,
                    ),
                )
                if phase == :pre_publication
                    Test.@test bind_cleanup_metadata.job_id === nothing
                    Test.@test bind_cleanup_metadata.pending_job_id == "770"
                    Test.@test bind_cleanup_metadata.submission_state ==
                               :submission_pending
                    _test_metamdbg_error(
                        () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                            bind_cleanup_metadata;
                            owner_token = bind_cleanup_metadata.owner_token,
                            confirm_not_submitted = true,
                        ),
                        ArgumentError,
                        r"pending scheduler job evidence.*cannot be reclaimed",
                    )
                    bind_cleanup_metadata =
                        Mycelia.bind_metamdbg_submission_reservation_job!(
                            bind_cleanup_metadata;
                            owner_token = bind_cleanup_metadata.owner_token,
                            job_id = "770",
                            confirm_submitted = true,
                        )
                end
                Test.@test bind_cleanup_metadata.job_id == "770"
                Test.@test !ispath(bind_cleanup_pending)
                Test.@test Set(readdir(temporary_root)) ==
                           parent_entries_before
                cleanup_reclaimed =
                    Mycelia.reclaim_metamdbg_submission_reservation!(
                        bind_cleanup_metadata;
                        owner_token = bind_cleanup_metadata.owner_token,
                        job_id = bind_cleanup_metadata.job_id,
                        confirm_cancelled = true,
                    )
                Test.@test cleanup_reclaimed.recovery_reason == :cancelled

                bind_kill_outdir = joinpath(
                    temporary_root,
                    "$(phase)-job-bind-sigkill-output",
                )
                bind_reservation =
                    _kill_metamdbg_job_record_publication!(
                        bind_kill_outdir,
                        valid_reads,
                        phase;
                        live_check = function (
                                process::Base.Process,
                                reservation::NamedTuple,
                        )
                            Test.@test Base.process_running(process)
                            pid_lock =
                                Mycelia._output_root_reservation_lock_path_from_canonical(
                                    bind_kill_outdir,
                                )
                            private_lock =
                                Mycelia._metamdbg_output_lock_path(
                                    bind_kill_outdir,
                                )
                            cleanup_reservation =
                                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                                    bind_kill_outdir,
                                )
                            pending_job_record =
                                Mycelia._metamdbg_pending_submission_job_path(
                                    bind_kill_outdir,
                                    reservation.output_root_reservation_capability,
                                    "771",
                                )
                            tracked_paths = String[
                                reservation.path,
                                reservation.output_root_reservation_marker,
                                pending_job_record,
                                pid_lock,
                                private_lock,
                                cleanup_reservation,
                            ]
                            tracked_identities = map(tracked_paths) do path
                                path_status = stat(path)
                                return (;
                                    device = path_status.device,
                                    inode = path_status.inode,
                                )
                            end
                            Test.@test isfile(reservation.job_marker) ==
                                       (phase == :post_publication)
                            live_metadata = only(
                                Mycelia.inspect_metamdbg_submission_reservations(
                                    bind_kill_outdir,
                                ),
                            )
                            Test.@test live_metadata.pending_job_id == "771"
                            Test.@test live_metadata.pending_job_complete
                            Test.@test live_metadata.submission_state ==
                                       (phase == :post_publication ?
                                        :submission_commit_cleanup_pending :
                                        :submission_pending)
                            _test_metamdbg_error(
                                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                                    live_metadata;
                                    owner_token = live_metadata.owner_token,
                                    confirm_not_submitted = true,
                                ),
                                ArgumentError,
                                r"pending scheduler job evidence.*cannot be reclaimed",
                            )
                            _test_metamdbg_error(
                                () -> Mycelia.inspect_metamdbg_submission_reservations(
                                    bind_kill_outdir;
                                    confirm_process_dead = true,
                                ),
                                ErrorException,
                                r"still names a live or remotely unverifiable process",
                            )
                            Test.@test all(ispath, tracked_paths)
                            observed_identities = map(tracked_paths) do path
                                path_status = stat(path)
                                return (;
                                    device = path_status.device,
                                    inode = path_status.inode,
                                )
                            end
                            Test.@test observed_identities ==
                                       tracked_identities
                            return nothing
                        end,
                    )
                bind_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        bind_kill_outdir;
                        confirm_process_dead = true,
                    ),
                )
                Test.@test bind_metadata.job_id == "771"
                Test.@test bind_metadata.pending_job_id === nothing
                Test.@test bind_metadata.private_lock_identity === nothing
                Test.@test bind_metadata.cleanup_reservation_identity ===
                           nothing
                reclaimed_bind =
                    Mycelia.reclaim_metamdbg_submission_reservation!(
                        bind_metadata;
                        owner_token = bind_metadata.owner_token,
                        job_id = bind_metadata.job_id,
                        confirm_cancelled = true,
                    )
                Test.@test reclaimed_bind.recovery_reason == :cancelled
                Test.@test !ispath(bind_reservation.path)
                Test.@test !ispath(
                    bind_reservation.output_root_reservation_marker,
                )
                Test.@test !ispath(
                    Mycelia._metamdbg_pending_submission_job_path(
                        bind_kill_outdir,
                        bind_reservation.output_root_reservation_capability,
                        "771",
                    ),
                )
                Test.@test isempty(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        bind_kill_outdir,
                    ),
                )
            end

            stale_clean_outdir = joinpath(
                temporary_root,
                "stale-clean-bind-metadata-output",
            )
            stale_clean_outputs =
                Mycelia._metamdbg_output_paths(stale_clean_outdir, 21)
            stale_clean_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            stale_clean_reservation = Mycelia._metamdbg_submission_reservation(
                stale_clean_outputs,
                stale_clean_contract,
                21;
                owner_token = "stale-clean-bind-owner",
            )
            Mycelia._with_metamdbg_output_lock(stale_clean_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    stale_clean_reservation,
                    stale_clean_outdir,
                )
            end
            stale_clean_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    stale_clean_outdir,
                ),
            )
            _test_metamdbg_error(
                () -> Mycelia._with_metamdbg_output_domain_lock(
                    stale_clean_outdir;
                    allowed_same_root_locks = (
                        stale_clean_reservation.output_root_reservation_marker,
                    ),
                ) do
                    Mycelia._bind_metamdbg_submission_job!(
                        stale_clean_reservation,
                        "776";
                        pre_job_record_publication_hook = (
                            _reservation::NamedTuple,
                            _path::AbstractString,
                        ) -> error("synthetic accepted-job bind failure"),
                    )
                end,
                ErrorException,
                r"synthetic accepted-job bind failure",
            )
            stale_pending_path =
                Mycelia._metamdbg_pending_submission_job_path(
                    stale_clean_outdir,
                    stale_clean_reservation.output_root_reservation_capability,
                    "776",
                )
            conflicting_pending_path =
                Mycelia._metamdbg_pending_submission_job_path(
                    stale_clean_outdir,
                    stale_clean_reservation.output_root_reservation_capability,
                    "777",
                )
            _test_metamdbg_error(
                () -> Mycelia.bind_metamdbg_submission_reservation_job!(
                    stale_clean_metadata;
                    owner_token = stale_clean_metadata.owner_token,
                    job_id = "777",
                    confirm_submitted = true,
                ),
                ErrorException,
                r"pending scheduler job evidence appeared after inspection",
            )
            Test.@test isfile(stale_pending_path)
            Test.@test !ispath(conflicting_pending_path)
            Test.@test !ispath(stale_clean_reservation.job_marker)
            exact_pending_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    stale_clean_outdir,
                ),
            )
            exact_pending_bound =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    exact_pending_metadata;
                    owner_token = exact_pending_metadata.owner_token,
                    job_id = "776",
                    confirm_submitted = true,
                )
            stale_clean_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    exact_pending_bound;
                    owner_token = exact_pending_bound.owner_token,
                    job_id = exact_pending_bound.job_id,
                    confirm_cancelled = true,
                )
            Test.@test stale_clean_reclaimed.recovery_reason == :cancelled

            pending_cleanup_outdir = joinpath(
                temporary_root,
                "committed-pending-cleanup-failure-output",
            )
            pending_cleanup_outputs = Mycelia._metamdbg_output_paths(
                pending_cleanup_outdir,
                21,
            )
            pending_cleanup_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            pending_cleanup_reservation =
                Mycelia._metamdbg_submission_reservation(
                    pending_cleanup_outputs,
                    pending_cleanup_contract,
                    21;
                    owner_token = "committed-pending-cleanup-owner",
                )
            Mycelia._with_metamdbg_output_lock(pending_cleanup_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    pending_cleanup_reservation,
                    pending_cleanup_outdir,
                )
            end
            _test_metamdbg_error(
                () -> Mycelia._with_metamdbg_output_domain_lock(
                    pending_cleanup_outdir;
                    allowed_same_root_locks = (
                        pending_cleanup_reservation.output_root_reservation_marker,
                    ),
                ) do
                    Mycelia._bind_metamdbg_submission_job!(
                        pending_cleanup_reservation,
                        "773";
                        pending_job_record_remover = _path -> error(
                            "synthetic pending removal failure",
                        ),
                    )
                end,
                ErrorException,
                r"synthetic pending removal failure",
            )
            pending_cleanup_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    pending_cleanup_outdir,
                ),
            )
            Test.@test pending_cleanup_metadata.job_id == "773"
            Test.@test pending_cleanup_metadata.pending_job_id == "773"
            Test.@test pending_cleanup_metadata.submission_state ==
                       :submission_commit_cleanup_pending
            _test_metamdbg_error(
                () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                    pending_cleanup_metadata;
                    owner_token = pending_cleanup_metadata.owner_token,
                    job_id = pending_cleanup_metadata.job_id,
                    confirm_cancelled = true,
                ),
                ArgumentError,
                r"pending scheduler job evidence.*cannot be reclaimed",
            )
            resumed_pending_cleanup =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    pending_cleanup_metadata;
                    owner_token = pending_cleanup_metadata.owner_token,
                    job_id = pending_cleanup_metadata.job_id,
                    confirm_submitted = true,
                )
            Test.@test resumed_pending_cleanup.job_id == "773"
            Test.@test resumed_pending_cleanup.pending_job_id === nothing
            resumed_pending_cleanup = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    pending_cleanup_outdir,
                ),
            )
            pending_cleanup_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    resumed_pending_cleanup;
                    owner_token = resumed_pending_cleanup.owner_token,
                    job_id = resumed_pending_cleanup.job_id,
                    confirm_cancelled = true,
                )
            Test.@test pending_cleanup_reclaimed.recovery_reason == :cancelled

            replacement_pending_outdir = joinpath(
                temporary_root,
                "pending-same-content-replacement-output",
            )
            replacement_pending_outputs = Mycelia._metamdbg_output_paths(
                replacement_pending_outdir,
                21,
            )
            replacement_pending_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            replacement_pending_reservation =
                Mycelia._metamdbg_submission_reservation(
                    replacement_pending_outputs,
                    replacement_pending_contract,
                    21;
                    owner_token = "pending-replacement-owner",
                )
            Mycelia._with_metamdbg_output_lock(replacement_pending_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    replacement_pending_reservation,
                    replacement_pending_outdir,
                )
            end
            held_original_pending = joinpath(
                temporary_root,
                "held-original-pending-job-record",
            )
            replace_pending_hook = function (
                    bound::NamedTuple,
                    pending_path::AbstractString,
            )
                mv(pending_path, held_original_pending)
                write(pending_path, bound.job_contents)
                chmod(pending_path, 0o600)
                return nothing
            end
            _test_metamdbg_error(
                () -> Mycelia._with_metamdbg_output_domain_lock(
                    replacement_pending_outdir;
                    allowed_same_root_locks = (
                        replacement_pending_reservation.output_root_reservation_marker,
                    ),
                ) do
                    Mycelia._bind_metamdbg_submission_job!(
                        replacement_pending_reservation,
                        "774";
                        pre_job_record_publication_hook = replace_pending_hook,
                    )
                end,
                ErrorException,
                r"pending submission job record identity changed",
            )
            replacement_pending_path =
                Mycelia._metamdbg_pending_submission_job_path(
                    replacement_pending_outdir,
                    replacement_pending_reservation.output_root_reservation_capability,
                    "774",
                )
            Test.@test isfile(replacement_pending_path)
            Test.@test isfile(held_original_pending)
            Test.@test !ispath(replacement_pending_reservation.job_marker)
            rm(replacement_pending_path)
            mv(held_original_pending, replacement_pending_path)
            replacement_pending_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    replacement_pending_outdir,
                ),
            )
            replacement_pending_bound =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    replacement_pending_metadata;
                    owner_token = replacement_pending_metadata.owner_token,
                    job_id = "774",
                    confirm_submitted = true,
                )
            replacement_pending_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    replacement_pending_bound;
                    owner_token = replacement_pending_bound.owner_token,
                    job_id = replacement_pending_bound.job_id,
                    confirm_cancelled = true,
                )
            Test.@test replacement_pending_reclaimed.recovery_reason ==
                       :cancelled

            for replacement_kind in (:owner, :shared)
                snapshot_outdir = joinpath(
                    temporary_root,
                    "pending-snapshot-$(replacement_kind)-replacement-output",
                )
                snapshot_outputs =
                    Mycelia._metamdbg_output_paths(snapshot_outdir, 21)
                snapshot_contract = Mycelia._metamdbg_input_contract(
                    Mycelia._metamdbg_selected_input(valid_reads, nothing),
                    3,
                )
                snapshot_reservation =
                    Mycelia._metamdbg_submission_reservation(
                        snapshot_outputs,
                        snapshot_contract,
                        21;
                        owner_token =
                            "pending-snapshot-$(replacement_kind)-owner",
                    )
                Mycelia._with_metamdbg_output_lock(snapshot_outdir) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        snapshot_reservation,
                        snapshot_outdir,
                    )
                end
                _test_metamdbg_error(
                    () -> Mycelia._with_metamdbg_output_domain_lock(
                        snapshot_outdir;
                        allowed_same_root_locks = (
                            snapshot_reservation.output_root_reservation_marker,
                        ),
                    ) do
                        Mycelia._bind_metamdbg_submission_job!(
                            snapshot_reservation,
                            "778";
                            pre_job_record_publication_hook = (
                                _reservation::NamedTuple,
                                _path::AbstractString,
                            ) -> error("synthetic snapshot bind failure"),
                        )
                    end,
                    ErrorException,
                    r"synthetic snapshot bind failure",
                )
                held_snapshot_path = joinpath(
                    temporary_root,
                    "held-pending-snapshot-$(replacement_kind)",
                )
                snapshot_hook = function (
                        records::Vector{NamedTuple},
                )
                    Test.@test length(records) == 1
                    if replacement_kind == :owner
                        mv(snapshot_reservation.path, held_snapshot_path)
                        mkdir(snapshot_reservation.path)
                        chmod(snapshot_reservation.path, 0o700)
                        write(
                            snapshot_reservation.contract_marker,
                            snapshot_reservation.contents,
                        )
                        chmod(snapshot_reservation.contract_marker, 0o600)
                    else
                        mv(
                            snapshot_reservation.output_root_reservation_marker,
                            held_snapshot_path,
                        )
                        write(
                            snapshot_reservation.output_root_reservation_marker,
                            snapshot_reservation.output_root_reservation_contents,
                        )
                        chmod(
                            snapshot_reservation.output_root_reservation_marker,
                            0o600,
                        )
                    end
                    return nothing
                end
                replacement_pattern = replacement_kind == :owner ?
                                      r"reservation was replaced after recovery inspection" :
                                      r"shared output-root reservation was replaced after recovery"
                _test_metamdbg_error(
                    () -> Mycelia._recover_metamdbg_pending_submission_job_records!(
                        snapshot_outdir;
                        post_initial_snapshot_hook = snapshot_hook,
                    ),
                    ErrorException,
                    replacement_pattern,
                )
                Test.@test !ispath(snapshot_reservation.job_marker)
                if replacement_kind == :owner
                    rm(snapshot_reservation.path; recursive = true)
                    mv(held_snapshot_path, snapshot_reservation.path)
                else
                    rm(snapshot_reservation.output_root_reservation_marker)
                    mv(
                        held_snapshot_path,
                        snapshot_reservation.output_root_reservation_marker,
                    )
                end
                snapshot_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        snapshot_outdir,
                    ),
                )
                Test.@test snapshot_metadata.pending_job_id == "778"
                snapshot_bound =
                    Mycelia.bind_metamdbg_submission_reservation_job!(
                        snapshot_metadata;
                        owner_token = snapshot_metadata.owner_token,
                        job_id = "778",
                        confirm_submitted = true,
                    )
                snapshot_reclaimed =
                    Mycelia.reclaim_metamdbg_submission_reservation!(
                        snapshot_bound;
                        owner_token = snapshot_bound.owner_token,
                        job_id = snapshot_bound.job_id,
                        confirm_cancelled = true,
                    )
                Test.@test snapshot_reclaimed.recovery_reason == :cancelled
            end

            symlink_pending_outdir = joinpath(
                temporary_root,
                "pending-symlink-race-output",
            )
            symlink_pending_outputs =
                Mycelia._metamdbg_output_paths(symlink_pending_outdir, 21)
            symlink_pending_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            symlink_pending_reservation =
                Mycelia._metamdbg_submission_reservation(
                    symlink_pending_outputs,
                    symlink_pending_contract,
                    21;
                    owner_token = "pending-symlink-owner",
                )
            Mycelia._with_metamdbg_output_lock(symlink_pending_outdir) do
                Mycelia._create_metamdbg_submission_reservation!(
                    symlink_pending_reservation,
                    symlink_pending_outdir,
                )
            end
            symlink_pending_path =
                Mycelia._metamdbg_pending_submission_job_path(
                    symlink_pending_outdir,
                    symlink_pending_reservation.output_root_reservation_capability,
                    "775",
                )
            write(symlink_pending_path, "")
            chmod(symlink_pending_path, 0o600)
            captured_empty_pending =
                Mycelia._metamdbg_pending_submission_job_record(
                    symlink_pending_reservation,
                    symlink_pending_path;
                    allow_incomplete = true,
                )
            held_empty_pending = joinpath(
                temporary_root,
                "held-empty-pending-job-record",
            )
            symlink_target = joinpath(
                temporary_root,
                "pending-symlink-do-not-truncate",
            )
            write(symlink_target, "do-not-truncate")
            chmod(symlink_target, 0o600)
            swap_pending_for_symlink = function (_pending::NamedTuple)
                mv(symlink_pending_path, held_empty_pending)
                symlink(symlink_target, symlink_pending_path)
                return nothing
            end
            _test_metamdbg_error(
                () -> Mycelia._complete_metamdbg_pending_submission_job_record!(
                    symlink_pending_reservation,
                    captured_empty_pending;
                    pre_descriptor_open_hook = swap_pending_for_symlink,
                ),
                ErrorException,
                r"descriptor was replaced before prefix completion",
            )
            Test.@test read(symlink_target, String) == "do-not-truncate"
            Test.@test islink(symlink_pending_path)
            Test.@test !ispath(symlink_pending_reservation.job_marker)
            rm(symlink_pending_path)
            mv(held_empty_pending, symlink_pending_path)
            symlink_pending_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    symlink_pending_outdir,
                ),
            )
            symlink_pending_bound =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    symlink_pending_metadata;
                    owner_token = symlink_pending_metadata.owner_token,
                    job_id = "775",
                    confirm_submitted = true,
                )
            symlink_pending_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    symlink_pending_bound;
                    owner_token = symlink_pending_bound.owner_token,
                    job_id = symlink_pending_bound.job_id,
                    confirm_cancelled = true,
                )
            Test.@test symlink_pending_reclaimed.recovery_reason == :cancelled

            pending_name_outdir = joinpath(
                temporary_root,
                "pending-name-job-bind-sigkill-output",
            )
            pending_name_reservation =
                _kill_metamdbg_job_record_publication!(
                    pending_name_outdir,
                    valid_reads,
                    :pending_name;
                    live_check = function (
                            process::Base.Process,
                            reservation::NamedTuple,
                    )
                        Test.@test Base.process_running(process)
                        pending_path =
                            Mycelia._metamdbg_pending_submission_job_path(
                                pending_name_outdir,
                                reservation.output_root_reservation_capability,
                                "771",
                            )
                        Test.@test isfile(pending_path)
                        Test.@test filesize(pending_path) == 0
                        Test.@test !ispath(reservation.job_marker)
                        live_metadata = only(
                            Mycelia.inspect_metamdbg_submission_reservations(
                                pending_name_outdir,
                            ),
                        )
                        Test.@test live_metadata.job_id === nothing
                        Test.@test live_metadata.pending_job_id == "771"
                        Test.@test !live_metadata.pending_job_complete
                        Test.@test live_metadata.submission_state ==
                                   :submission_pending
                        _test_metamdbg_error(
                            () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                                live_metadata;
                                owner_token = live_metadata.owner_token,
                                confirm_not_submitted = true,
                            ),
                            ArgumentError,
                            r"pending scheduler job evidence.*cannot be reclaimed",
                        )
                        _test_metamdbg_error(
                            () -> Mycelia.inspect_metamdbg_submission_reservations(
                                pending_name_outdir;
                                confirm_process_dead = true,
                            ),
                            ErrorException,
                            r"still names a live or remotely unverifiable process",
                        )
                        Test.@test isfile(pending_path)
                        Test.@test filesize(pending_path) == 0
                        return nothing
                    end,
                )
            pending_name_path =
                Mycelia._metamdbg_pending_submission_job_path(
                    pending_name_outdir,
                    pending_name_reservation.output_root_reservation_capability,
                    "771",
                )
            failing_pending_recovery = function (outdir::AbstractString)
                return Mycelia._recover_metamdbg_pending_submission_job_records!(
                    outdir;
                    pre_promotion_hook = (
                        _reservation::NamedTuple,
                        _pending::NamedTuple,
                    ) -> error("synthetic pending promotion failure"),
                )
            end
            _test_metamdbg_error(
                () -> Mycelia._inspect_metamdbg_submission_reservations(
                    pending_name_outdir;
                    confirm_process_dead = true,
                    pending_recovery_function = failing_pending_recovery,
                ),
                ErrorException,
                r"synthetic pending promotion failure",
            )
            pending_name_pid_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    pending_name_outdir,
                )
            pending_name_private_lock =
                Mycelia._metamdbg_output_lock_path(pending_name_outdir)
            pending_name_cleanup_reservation =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    pending_name_outdir,
                )
            Test.@test !ispath(pending_name_pid_lock)
            Test.@test !ispath(pending_name_private_lock)
            Test.@test !ispath(pending_name_cleanup_reservation)
            Test.@test isfile(pending_name_path)
            recovered_pending_name = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    pending_name_outdir,
                ),
            )
            Test.@test recovered_pending_name.job_id === nothing
            Test.@test recovered_pending_name.pending_job_id == "771"
            Test.@test recovered_pending_name.pending_job_complete
            recovered_pending_name =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    recovered_pending_name;
                    owner_token = recovered_pending_name.owner_token,
                    job_id = "771",
                    confirm_submitted = true,
                )
            Test.@test recovered_pending_name.job_id == "771"
            Test.@test recovered_pending_name.pending_job_id === nothing
            Test.@test isfile(pending_name_reservation.job_marker)
            Test.@test !ispath(pending_name_path)
            reclaimed_pending_name =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    recovered_pending_name;
                    owner_token = recovered_pending_name.owner_token,
                    job_id = recovered_pending_name.job_id,
                    confirm_cancelled = true,
                )
            Test.@test reclaimed_pending_name.recovery_reason == :cancelled
            Test.@test isempty(
                Mycelia.inspect_metamdbg_submission_reservations(
                    pending_name_outdir,
                ),
            )

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

            post_rename_outdir =
                joinpath(temporary_root, "post-rename-interruption")
            post_rename_outputs =
                Mycelia._metamdbg_output_paths(post_rename_outdir, 21)
            post_rename_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(valid_reads, nothing),
                3,
            )
            post_rename_reservation =
                Mycelia._metamdbg_submission_reservation(
                    post_rename_outputs,
                    post_rename_contract,
                    21;
                    owner_token = "post-rename-interruption-fixture",
                )
            post_rename_interrupt = InterruptException()
            post_rename_hook = (_reservation::NamedTuple) ->
                throw(post_rename_interrupt)
            observed_post_rename_interrupt = try
                Mycelia._with_metamdbg_output_lock(post_rename_outdir) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        post_rename_reservation,
                        post_rename_outdir;
                        post_rename_hook,
                    )
                end
                nothing
            catch caught
                caught
            end
            Test.@test observed_post_rename_interrupt === post_rename_interrupt
            Test.@test !ispath(post_rename_reservation.path)
            Test.@test !ispath(
                post_rename_reservation.output_root_reservation_marker,
            )
            Test.@test isempty(
                Mycelia._metamdbg_submission_reservation_paths(
                    post_rename_outdir,
                ),
            )
            post_rename_prefix =
                Mycelia._metamdbg_submission_reservation_prefix(
                    post_rename_outdir,
                )
            Test.@test isempty(filter(
                entry -> startswith(entry, post_rename_prefix),
                readdir(dirname(post_rename_outdir)),
            ))

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
            filesystem_root = dirname(abspath("/"))
            filesystem_root_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    filesystem_root,
                )
            root_probe_ancestor_locks =
                Mycelia._ancestor_output_root_reservation_lock_paths(
                    joinpath(temporary_root, "root-ancestor-probe"),
                )
            Test.@test filesystem_root_lock in root_probe_ancestor_locks

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
                    "active same-root output-root reservation",
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

            case_alias_root =
                joinpath(temporary_root, "Filesystem-Case-Alias")
            case_alias_spelling =
                joinpath(temporary_root, "filesystem-case-alias")
            mkpath(case_alias_root)
            if isdir(case_alias_spelling) && Base.Filesystem.samefile(
                    case_alias_root,
                    case_alias_spelling,
                )
                Test.@test Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(case_alias_root),
                ) == Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(
                        case_alias_spelling,
                    ),
                )
                Mycelia._with_metamdbg_output_domain_lock(
                    case_alias_root,
                ) do
                    case_alias_conflict = fetch(@async try
                        Mycelia._with_unicycler_output_lock(
                            case_alias_spelling,
                        ) do _reserved
                            :unexpected
                        end
                    catch caught
                        caught
                    end)
                    Test.@test case_alias_conflict isa ArgumentError
                    Test.@test occursin(
                        "active same-root output-root reservation",
                        sprint(showerror, case_alias_conflict),
                    )
                end
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

            runtime_fence_root = joinpath(
                temporary_root,
                "shared-domain-runtime-fence",
            )
            mkpath(runtime_fence_root)
            runtime_fence_outputs =
                Mycelia._metamdbg_output_paths(runtime_fence_root, 21)
            runtime_fence_reservation =
                Mycelia._metamdbg_submission_reservation(
                    runtime_fence_outputs,
                    input_contract,
                    21;
                    owner_token = "shared-domain-runtime-fence-owner",
                )
            runtime_fence_marker =
                runtime_fence_reservation.runtime_output_root_reservation_marker
            mkdir(runtime_fence_marker; mode = 0o700)
            Mycelia._fsync_metamdbg_directory(runtime_fence_marker)
            Mycelia._fsync_metamdbg_directory(
                dirname(runtime_fence_marker),
            )
            missing_fenced_parent = joinpath(
                runtime_fence_root,
                "missing-parent",
            )
            nested_fenced_output = joinpath(
                missing_fenced_parent,
                "output",
            )
            nested_fenced_pid_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    nested_fenced_output,
                )
            nested_fenced_private_lock =
                Mycelia._metamdbg_output_lock_path(nested_fenced_output)
            function require_runtime_fence_preflight!(
                    action::Function,
                    action_entered::Base.RefValue{Bool},
            )::Nothing
                observed_error = try
                    action()
                    nothing
                catch caught
                    caught isa InterruptException && rethrow()
                    caught
                end
                Test.@test observed_error isa ArgumentError
                if observed_error isa ArgumentError
                    Test.@test occursin(
                        runtime_fence_marker,
                        sprint(showerror, observed_error),
                    )
                end
                Test.@test !action_entered[]
                Test.@test !ispath(missing_fenced_parent)
                Test.@test !ispath(nested_fenced_pid_lock)
                Test.@test !ispath(nested_fenced_private_lock)
                return nothing
            end

            nested_meta_entered = Ref(false)
            require_runtime_fence_preflight!(nested_meta_entered) do
                Mycelia._with_metamdbg_output_domain_lock(
                    nested_fenced_output,
                ) do
                    nested_meta_entered[] = true
                end
            end
            nested_unicycler_entered = Ref(false)
            require_runtime_fence_preflight!(nested_unicycler_entered) do
                Mycelia._with_unicycler_output_lock(
                    nested_fenced_output,
                ) do _reserved
                    nested_unicycler_entered[] = true
                end
            end
            nested_autocycler_entered = Ref(false)
            require_runtime_fence_preflight!(nested_autocycler_entered) do
                Mycelia._with_autocycler_output_lock(
                    nested_fenced_output,
                ) do _reserved
                    nested_autocycler_entered[] = true
                end
            end
            rm(runtime_fence_marker; recursive = true)
            Mycelia._fsync_metamdbg_directory(
                dirname(runtime_fence_marker),
            )

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

            normal_cleanup_root =
                joinpath(temporary_root, "shared-domain-normal-cleanup")
            normal_cleanup_lock =
                Mycelia._metamdbg_output_lock_path(normal_cleanup_root)
            normal_cleanup_generic_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(
                        normal_cleanup_root,
                    ),
                )
            normal_cleanup_sentinel =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    normal_cleanup_root,
                )
            Test.@test Mycelia._with_metamdbg_output_domain_lock(
                normal_cleanup_root,
            ) do
                Test.@test isdir(normal_cleanup_sentinel)
                :normal_cleanup
            end == :normal_cleanup
            Test.@test !ispath(normal_cleanup_lock)
            Test.@test !ispath(normal_cleanup_generic_lock)
            Test.@test !ispath(normal_cleanup_sentinel)

            normal_exception_root =
                joinpath(temporary_root, "shared-domain-normal-exception")
            normal_exception_lock =
                Mycelia._metamdbg_output_lock_path(normal_exception_root)
            normal_exception_generic_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(
                        normal_exception_root,
                    ),
                )
            normal_exception_sentinel =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    normal_exception_root,
                )
            normal_primary_error = ErrorException("normal primary failure")
            observed_normal_primary = try
                Mycelia._with_metamdbg_output_domain_lock(
                    normal_exception_root,
                ) do
                    throw(normal_primary_error)
                end
                nothing
            catch caught
                caught
            end
            Test.@test observed_normal_primary === normal_primary_error
            Test.@test !ispath(normal_exception_lock)
            Test.@test !ispath(normal_exception_generic_lock)
            Test.@test !ispath(normal_exception_sentinel)

            function block_output_domain_cleanup!(
                    root::AbstractString,
            )::String
                sentinel =
                    Mycelia._metamdbg_lifecycle_cleanup_reservation_path(root)
                Test.@test isdir(sentinel)
                write(joinpath(sentinel, "blocker"), "retain cleanup evidence\n")
                return sentinel
            end

            for (case_label, primary_error) in (
                    (
                        :ordinary,
                        ErrorException("synthetic output-domain primary"),
                    ),
                    (:interrupt, InterruptException()),
                )
                primary_root = joinpath(
                    temporary_root,
                    "output-domain-primary-$(case_label)",
                )
                primary_generic_lock =
                    Mycelia._output_root_reservation_lock_path_from_canonical(
                        Mycelia._metamdbg_canonical_output_path(primary_root),
                    )
                retained_sentinel = Ref("")
                observed_primary = Test.@test_logs (
                    :warn,
                    r"output-domain cleanup failed while preserving",
                ) min_level=Logging.Warn match_mode=:any begin
                    try
                        Mycelia._with_metamdbg_output_domain_lock(
                            primary_root,
                        ) do
                            retained_sentinel[] =
                                block_output_domain_cleanup!(primary_root)
                            throw(primary_error)
                        end
                        nothing
                    catch caught
                        caught
                    end
                end
                Test.@test observed_primary === primary_error
                Test.@test isdir(retained_sentinel[])
                Test.@test !ispath(primary_generic_lock)
                rm(retained_sentinel[]; recursive = true, force = true)
            end

            success_cleanup_root = joinpath(
                temporary_root,
                "output-domain-success-cleanup-failure",
            )
            success_generic_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(
                        success_cleanup_root,
                    ),
                )
            success_sentinel = Ref("")
            success_cleanup_failure = try
                Mycelia._with_metamdbg_output_domain_lock(
                    success_cleanup_root,
                ) do
                    success_sentinel[] =
                        block_output_domain_cleanup!(success_cleanup_root)
                    return :success_before_cleanup
                end
                nothing
            catch caught
                caught
            end
            Test.@test success_cleanup_failure isa Exception
            Test.@test isdir(success_sentinel[])
            Test.@test !ispath(success_generic_lock)
            rm(success_sentinel[]; recursive = true, force = true)

            cleanup_parent =
                joinpath(temporary_root, "shared-domain-cleanup-failure")
            cleanup_root = joinpath(cleanup_parent, "target")
            cleanup_lock = Mycelia._metamdbg_output_lock_path(cleanup_root)
            cleanup_generic_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(cleanup_root),
                )
            cleanup_sentinel =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    cleanup_root,
                )
            cleanup_failure = try
                Mycelia._with_metamdbg_output_domain_lock(cleanup_root) do
                    Test.@test isdir(cleanup_sentinel)
                    write(joinpath(cleanup_lock, "blocker"), "block cleanup\n")
                    return :cleanup_should_fail
                end
                nothing
            catch caught
                caught
            end
            Test.@test cleanup_failure isa Exception
            Test.@test isdir(cleanup_lock)
            Test.@test isdir(cleanup_sentinel)
            Test.@test (stat(cleanup_sentinel).mode & 0o777) == 0o700
            Test.@test !ispath(cleanup_generic_lock)
            GC.gc(true)
            Test.@test isdir(cleanup_sentinel)
            _test_metamdbg_error(
                () -> Mycelia._with_unicycler_output_lock(cleanup_root) do _
                    :unexpected
                end,
                ArgumentError,
                r"active same-root output-root reservation",
            )
            _test_metamdbg_error(
                () -> Mycelia._with_autocycler_output_lock(
                    joinpath(cleanup_root, "child"),
                ) do _
                    :unexpected
                end,
                ArgumentError,
                r"active ancestor output-root reservation",
            )
            _test_metamdbg_error(
                () -> Mycelia._with_unicycler_output_lock(cleanup_parent) do _
                    :unexpected
                end,
                ArgumentError,
                r"active descendant output-root reservation",
            )
            cleanup_sibling = joinpath(cleanup_parent, "sibling")
            Test.@test Mycelia._with_unicycler_output_lock(
                cleanup_sibling,
            ) do _
                :cleanup_sibling_reserved
            end == :cleanup_sibling_reserved
            rm(cleanup_lock; recursive = true)
            rm(cleanup_sentinel)

            exceptional_cleanup_root = joinpath(
                temporary_root,
                "shared-domain-exceptional-cleanup-failure",
            )
            exceptional_cleanup_lock =
                Mycelia._metamdbg_output_lock_path(exceptional_cleanup_root)
            exceptional_cleanup_sentinel =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    exceptional_cleanup_root,
                )
            exceptional_primary_error =
                ErrorException("primary action failure wins")
            observed_exceptional_primary = try
                Mycelia._with_metamdbg_output_domain_lock(
                    exceptional_cleanup_root,
                ) do
                    write(
                        joinpath(exceptional_cleanup_lock, "blocker"),
                        "block cleanup\n",
                    )
                    throw(exceptional_primary_error)
                end
                nothing
            catch caught
                caught
            end
            Test.@test observed_exceptional_primary ===
                       exceptional_primary_error
            Test.@test isdir(exceptional_cleanup_lock)
            Test.@test isdir(exceptional_cleanup_sentinel)
            rm(exceptional_cleanup_lock; recursive = true)
            rm(exceptional_cleanup_sentinel)

            replaced_cleanup_root = joinpath(
                temporary_root,
                "shared-domain-replaced-cleanup-lock",
            )
            replaced_cleanup_lock =
                Mycelia._metamdbg_output_lock_path(replaced_cleanup_root)
            replaced_original_lock = replaced_cleanup_lock * ".original"
            replaced_cleanup_sentinel =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    replaced_cleanup_root,
                )
            replaced_cleanup_error = try
                Mycelia._with_metamdbg_output_domain_lock(
                    replaced_cleanup_root,
                ) do
                    mv(replaced_cleanup_lock, replaced_original_lock)
                    mkdir(replaced_cleanup_lock)
                    return :replacement_should_fail
                end
                nothing
            catch caught
                caught
            end
            Test.@test replaced_cleanup_error isa Exception
            Test.@test occursin(
                "was replaced before cleanup",
                sprint(showerror, replaced_cleanup_error),
            )
            Test.@test isdir(replaced_cleanup_lock)
            Test.@test isdir(replaced_cleanup_sentinel)
            rm(replaced_cleanup_lock)
            rm(replaced_original_lock)
            rm(replaced_cleanup_sentinel)

            long_component = repeat("l", 240)
            long_root = joinpath(temporary_root, long_component)
            long_sibling_parent = joinpath(temporary_root, "long-sibling")
            mkpath(long_sibling_parent)
            long_same_basename_root =
                joinpath(long_sibling_parent, long_component)
            long_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(long_root),
                )
            long_private_lock =
                Mycelia._metamdbg_output_lock_path(long_root)
            long_cleanup_sentinel =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(long_root)
            long_same_basename_lock =
                Mycelia._output_root_reservation_lock_path_from_canonical(
                    Mycelia._metamdbg_canonical_output_path(
                        long_same_basename_root,
                    ),
                )
            Test.@test basename(long_lock) !=
                       basename(long_same_basename_lock)
            for hashed_path in (
                    long_lock,
                    long_private_lock,
                    long_cleanup_sentinel,
                )
                Test.@test ncodeunits(basename(hashed_path)) <= 255
                Test.@test !occursin(long_component, basename(hashed_path))
            end
            long_reservation = reserve_output_root!(
                long_root,
                "long-basename-owner",
            )
            for hashed_path in (
                    long_reservation.path,
                    long_reservation.output_root_reservation_marker,
                    long_reservation.runtime_output_root_reservation_marker,
                )
                Test.@test ncodeunits(basename(hashed_path)) <= 255
                Test.@test !occursin(long_component, basename(hashed_path))
            end
            Test.@test isdir(long_reservation.path)
            Test.@test isfile(
                long_reservation.output_root_reservation_marker,
            )
            Test.@test length(Set((
                long_cleanup_sentinel,
                long_reservation.output_root_reservation_marker,
                long_reservation.runtime_output_root_reservation_marker,
            ))) == 3
            Mycelia._with_metamdbg_output_lock(long_root) do
                Mycelia._remove_metamdbg_submission_reservation!(
                    long_reservation,
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
            Test.@test occursin(
                "runtime refuses pending scheduler job evidence",
                script,
            )
            Test.@test occursin(
                "runtime refuses ambiguous pending scheduler job evidence",
                script,
            )
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
            Test.@test !occursin("reservation_tombstone=", script)
            Test.@test occursin(
                "mv -n -- \"\$submission_reservation_dir\" " *
                "\"\$runtime_submission_reservation_dir\"",
                script,
            )
            Test.@test occursin(
                "mv -n -- \"\$runtime_submission_reservation_dir\" " *
                "\"\$consumed_submission_reservation_dir\"",
                script,
            )
            Test.@test occursin(
                "mv -n -- \"\$contigs_new\" \"\$contigs_gz\"",
                script,
            )
            Test.@test occursin(
                "mv -n -- \"\$completion_new\" \"\$completion_marker\"",
                script,
            )
            Test.@test !occursin("mv -f", script)
            Test.@test !occursin("rm -f", script)
            critical_mutation_lines = filter(
                line -> occursin("mv -n --", line) ||
                        occursin("rm --", line) ||
                        occursin("rm -rf --", line),
                split(script, '\n'),
            )
            Test.@test !isempty(critical_mutation_lines)
            Test.@test all(critical_mutation_lines) do line
                normalized_line = strip(line)
                return startswith(normalized_line, "if ! ") ||
                       startswith(normalized_line, "elif ! ") ||
                       (
                           startswith(normalized_line, "elif mv -n --") &&
                           occursin("&&", normalized_line)
                       )
            end
            Test.@test occursin("fsync_file_and_parent", script)
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
            Test.@test occursin("sha256_stream < \"\$path\"", script)
            Test.@test !occursin("sha256sum -- \"\$path\"", script)
            Test.@test occursin(
                only(input_contract.contract.inputs).sha256,
                script,
            )
            script_lines = split(script, '\n')
            lock_position = findfirst(==("lock_acquired=1"), script_lines)
            runtime_marker_position = findfirst(
                ==("output_root_runtime_reservation_acquired=1"),
                script_lines,
            )
            validation_positions = findall(
                ==("validate_metamdbg_inputs"),
                script_lines,
            )
            reuse_position = findfirst(==("contract_exists=0"), script_lines)
            finalization_position = findfirst(
                ==("if [ \"\$contract_exists\" -eq 1 ]; then"),
                script_lines,
            )
            Test.@test length(validation_positions) == 5
            Test.@test runtime_marker_position < lock_position
            Test.@test lock_position < first(validation_positions) <
                       reuse_position < validation_positions[2] <
                       validation_positions[3] < finalization_position <
                       validation_positions[4] < validation_positions[5]
            capture_position = findfirst(
                ==("capture_metamdbg_artifact_snapshot() {"),
                script_lines,
            )
            unchanged_position = findfirst(
                ==("require_unchanged_metamdbg_artifacts() {"),
                script_lines,
            )
            contract_validation_position = findfirst(
                ==("require_unchanged_metamdbg_contract() {"),
                script_lines,
            )
            capture_lines = script_lines[
                capture_position:(unchanged_position - 1)
            ]
            unchanged_lines = script_lines[
                unchanged_position:(contract_validation_position - 1)
            ]
            Test.@test count(
                ==("  observe_metamdbg_artifact_snapshot || return 1"),
                capture_lines,
            ) == 2
            Test.@test count(
                ==("  validate_contigs \"\$contigs_gz\" || return 1"),
                capture_lines,
            ) == 1
            Test.@test count(
                ==("  validate_gfa \"\$graph_source\" || return 1"),
                capture_lines,
            ) == 1
            Test.@test count(
                ==("  observe_metamdbg_artifact_snapshot || return 1"),
                unchanged_lines,
            ) == 1
            Test.@test !any(
                line -> startswith(strip(line), "validate_"),
                unchanged_lines,
            )
            Test.@test occursin("reservation_inventory_index=0", script)
            Test.@test occursin(
                "reservation_inventory_index=\$((" *
                "reservation_inventory_index + 1))",
                script,
            )
            Test.@test occursin(
                "submission-reservation.\${reservation_inventory_index}.paths",
                script,
            )
            Test.@test !occursin(
                "reservation_entries_file=\"\$secure_tmpdir/" *
                "submission-reservation.paths\"",
                script,
            )
            Test.@test occursin(
                "contigs_size=\$artifact_contigs_size",
                script,
            )
            Test.@test occursin(
                "graph_sha256=\$artifact_graph_sha256",
                script,
            )
            Test.@test !occursin("\$\$", script)
            Test.@test !occursin("\${contigs_gz}.tmp", script)
            Test.@test !occursin("rmdir -- \"\$lock_dir\" || true", script)
            Test.@test !occursin("< <(find", script)
            Test.@test !occursin("find -P", script)
            Test.@test occursin("os.scandir(directory_descriptor)", script)
            Test.@test occursin("O_NOFOLLOW", script)
            Test.@test occursin("reservations-one-level", script)
            Test.@test occursin("reservations-recursive", script)
            Test.@test occursin("MAX_DEPTH = 256", script)
            Test.@test occursin("MAX_ENTRIES = 1000000", script)
            Test.@test occursin(
                "MAX_MANIFEST_BYTES = 268435456",
                script,
            )
            Test.@test occursin("while stack:", script)
            Test.@test occursin(
                "if second_manifest != first_manifest:",
                script,
            )
            Test.@test occursin(
                "successful inventory contains an unterminated NUL record",
                script,
            )
            Test.@test occursin(
                "could not completely enumerate \$inventory_label after " *
                "3 attempts",
                script,
            )
            Test.@test occursin("metamdbg_output_root_identity", script)
            Test.@test occursin(
                Mycelia._output_root_reservation_identity(dirname(abspath("/"))),
                script,
            )
            Test.@test occursin("  while :; do", script)
            Test.@test occursin(
                "require_owned_runtime_output_root_reservation",
                script,
            )
            Test.@test occursin(
                "bound_runtime_output_root_reservation_identity",
                script,
            )
            Test.@test !occursin(
                "root_identity=\$(sha256_text",
                script,
            )
            Test.@test occursin(
                ".mycelia-metamdbg-tmp.\${outdir_identity}.",
                script,
            )
            private_cleanup_position = findfirst(
                line -> occursin("! rmdir -- \"\$lock_dir\"", line),
                script_lines,
            )
            shared_cleanup_position = findfirst(
                line -> occursin(
                    "rmdir -- \"\$runtime_output_root_reservation\"",
                    line,
                ),
                script_lines,
            )
            Test.@test private_cleanup_position < shared_cleanup_position

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

            function create_bound_runtime_fixture(
                    runtime_outdir::String,
                    owner_token::String,
                    job_id::String;
                    fixture_input_contract::NamedTuple = input_contract,
            )::Tuple{NamedTuple, NamedTuple}
                runtime_outputs =
                    Mycelia._metamdbg_output_paths(runtime_outdir, 21)
                runtime_reservation =
                    Mycelia._metamdbg_submission_reservation(
                        runtime_outputs,
                        fixture_input_contract,
                        21;
                        owner_token,
                    )
                Mycelia._with_metamdbg_output_domain_lock(runtime_outdir) do
                    Mycelia._create_metamdbg_submission_reservation!(
                        runtime_reservation,
                        runtime_outdir,
                    )
                end
                bound_reservation =
                    Mycelia._bind_metamdbg_submission_job!(
                        runtime_reservation,
                        job_id,
                    )
                return runtime_outputs, bound_reservation
            end

            function remove_queued_runtime_fixture!(
                    runtime_outdir::String,
                    runtime_reservation::NamedTuple,
            )::Nothing
                Mycelia._with_metamdbg_output_lock(runtime_outdir) do
                    Mycelia._remove_metamdbg_submission_reservation!(
                        runtime_reservation,
                    )
                end
                return nothing
            end

            function write_side_effect_failure_wrapper!(
                    wrapper_path::String,
                    real_command::String,
                    target_path::String,
                    exit_code::Int,
            )::Nothing
                write(
                    wrapper_path,
                    "#!/usr/bin/env bash\n" *
                    "set -euo pipefail\n" *
                    "if [ \"\$#\" -gt 0 ] && [ \"\${!#}\" = " *
                    Base.shell_escape(target_path) * " ]; then\n" *
                    "  " * Base.shell_escape(real_command) *
                    " \"\$@\"\n" *
                    "  exit $(exit_code)\n" *
                    "fi\n" *
                    "exec " * Base.shell_escape(real_command) *
                    " \"\$@\"\n",
                )
                chmod(wrapper_path, 0o700)
                Test.@test success(Cmd(["bash", "-n", wrapper_path]))
                return nothing
            end

            function write_no_clobber_failure_wrapper!(
                    wrapper_path::String,
                    real_mv::String,
                    target_path::String,
                    occupant_command::String,
            )::Nothing
                write(
                    wrapper_path,
                    "#!/usr/bin/env bash\n" *
                    "set -euo pipefail\n" *
                    "if [ \"\$#\" -gt 0 ] && [ \"\${!#}\" = " *
                    Base.shell_escape(target_path) * " ]; then\n" *
                    "  $(occupant_command)\n" *
                    "  exit 73\n" *
                    "fi\n" *
                    "exec " * Base.shell_escape(real_mv) *
                    " \"\$@\"\n",
                )
                chmod(wrapper_path, 0o700)
                Test.@test success(Cmd(["bash", "-n", wrapper_path]))
                return nothing
            end

            contigs_race_outdir = joinpath(
                temporary_root,
                "executor-contigs-publication-race",
            )
            contigs_race_outputs, contigs_race_reservation =
                create_bound_runtime_fixture(
                    contigs_race_outdir,
                    "executor-contigs-publication-race-owner",
                    "881",
                )
            contigs_race_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                contigs_race_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = contigs_race_reservation,
                pre_contigs_publication_hook =
                    "printf '>replacement-contig\\nTGCA\\n' | gzip -c > " *
                    "\"\$contigs_gz\"",
            )
            contigs_race_script_path = joinpath(
                temporary_root,
                "executor-contigs-publication-race.sh",
            )
            write(contigs_race_script_path, contigs_race_script)
            real_mv_path = something(Sys.which("mv"), "/bin/mv")
            contigs_race_mv_directory = joinpath(
                temporary_root,
                "executor-contigs-publication-race-bin",
            )
            mkpath(contigs_race_mv_directory)
            write_no_clobber_failure_wrapper!(
                joinpath(contigs_race_mv_directory, "mv"),
                real_mv_path,
                contigs_race_outputs.contigs_gz,
                ":",
            )
            contigs_race_log = IOBuffer()
            contigs_race_process = Base.withenv(
                "SLURM_JOB_ID" => contigs_race_reservation.job_id,
                "PATH" =>
                    "$(contigs_race_mv_directory):$(get(ENV, "PATH", ""))",
            ) do
                return run(pipeline(
                    ignorestatus(`bash $(contigs_race_script_path)`),
                    stdout = contigs_race_log,
                    stderr = contigs_race_log,
                ))
            end
            contigs_race_output = String(take!(contigs_race_log))
            Test.@test !success(contigs_race_process)
            Test.@test occursin(
                "refused to overwrite a concurrent compressed-contigs " *
                "artifact",
                contigs_race_output,
            )
            Test.@test read(
                `gzip -cd -- $(contigs_race_outputs.contigs_gz)`,
                String,
            ) == ">replacement-contig\nTGCA\n"
            Test.@test !ispath(contigs_race_outputs.completion_marker)
            Test.@test !ispath(contigs_race_reservation.path)
            Test.@test isdir(contigs_race_reservation.runtime_path)
            Test.@test !ispath(
                contigs_race_reservation.output_root_reservation_marker,
            )
            Test.@test isdir(
                contigs_race_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test isdir(
                Mycelia._metamdbg_output_lock_path(contigs_race_outdir),
            )
            Test.@test any(
                name -> startswith(name, ".mycelia-metamdbg-tmp."),
                readdir(temporary_root),
            )

            for (
                    publication_name,
                    target_field,
                    job_id,
                    expected_failure,
                    generic_failure,
                ) in (
                    (
                        "graph-alias",
                        :graph_alias,
                        "892",
                        "refused to overwrite a concurrent graph alias",
                        "graph-alias publication command failed",
                    ),
                    (
                        "provenance-marker",
                        :contract_marker,
                        "893",
                        "refused to overwrite a concurrent provenance marker",
                        "provenance-marker publication command failed",
                    ),
                    (
                        "completion-manifest",
                        :completion_marker,
                        "894",
                        "refused to overwrite a concurrent completion manifest",
                        "completion-manifest publication command failed",
                    ),
                )
                no_clobber_outdir = joinpath(
                    temporary_root,
                    "executor-$(publication_name)-no-clobber-failure",
                )
                no_clobber_outputs, no_clobber_reservation =
                    create_bound_runtime_fixture(
                        no_clobber_outdir,
                        "executor-$(publication_name)-no-clobber-owner",
                        job_id,
                    )
                target_path = getproperty(no_clobber_outputs, target_field)
                replacement =
                    "concurrent $(publication_name) publication must survive\n"
                no_clobber_mv_directory = joinpath(
                    temporary_root,
                    "executor-$(publication_name)-no-clobber-bin",
                )
                mkpath(no_clobber_mv_directory)
                write_no_clobber_failure_wrapper!(
                    joinpath(no_clobber_mv_directory, "mv"),
                    real_mv_path,
                    target_path,
                    "printf '%s' " * Base.shell_escape(replacement) *
                    " > " * Base.shell_escape(target_path),
                )
                no_clobber_script = Mycelia._metamdbg_executor_script(
                    fake_asm,
                    fake_gfa,
                    no_clobber_outputs,
                    21,
                    input_contract;
                    conda_runner = fake_conda,
                    submission_reservation = no_clobber_reservation,
                )
                no_clobber_script_path = joinpath(
                    temporary_root,
                    "executor-$(publication_name)-no-clobber.sh",
                )
                write(no_clobber_script_path, no_clobber_script)
                no_clobber_log = IOBuffer()
                no_clobber_process = Base.withenv(
                    "SLURM_JOB_ID" => no_clobber_reservation.job_id,
                    "PATH" =>
                        "$(no_clobber_mv_directory):$(get(ENV, "PATH", ""))",
                ) do
                    return run(pipeline(
                        ignorestatus(`bash $(no_clobber_script_path)`),
                        stdout = no_clobber_log,
                        stderr = no_clobber_log,
                    ))
                end
                no_clobber_output = String(take!(no_clobber_log))
                Test.@test !success(no_clobber_process)
                Test.@test occursin(expected_failure, no_clobber_output)
                Test.@test !occursin(generic_failure, no_clobber_output)
                Test.@test read(target_path, String) == replacement
                Test.@test !ispath(no_clobber_reservation.path)
            end

            completion_race_outdir = joinpath(
                temporary_root,
                "executor-completion-publication-replacement",
            )
            completion_race_outputs, completion_race_reservation =
                create_bound_runtime_fixture(
                    completion_race_outdir,
                    "executor-completion-publication-replacement-owner",
                    "882",
                )
            completion_original =
                completion_race_outputs.completion_marker * ".original"
            completion_replacement = "replacement completion must survive\n"
            completion_race_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                completion_race_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = completion_race_reservation,
                post_completion_publication_hook =
                    "mv -- \"\$completion_marker\" " *
                    Base.shell_escape(completion_original) * "; " *
                    "printf '%s' " *
                    Base.shell_escape(completion_replacement) *
                    " > \"\$completion_marker\"",
            )
            completion_race_script_path = joinpath(
                temporary_root,
                "executor-completion-publication-replacement.sh",
            )
            write(completion_race_script_path, completion_race_script)
            completion_race_log = IOBuffer()
            completion_race_process = Base.withenv(
                "SLURM_JOB_ID" => completion_race_reservation.job_id,
            ) do
                return run(pipeline(
                    ignorestatus(`bash $(completion_race_script_path)`),
                    stdout = completion_race_log,
                    stderr = completion_race_log,
                ))
            end
            completion_race_output = String(take!(completion_race_log))
            Test.@test !success(completion_race_process)
            Test.@test occursin(
                "published completion manifest changed unexpectedly; " *
                "retained fail-closed evidence",
                completion_race_output,
            )
            Test.@test read(
                completion_race_outputs.completion_marker,
                String,
            ) == completion_replacement
            Test.@test isfile(completion_original)
            Test.@test !ispath(completion_race_outputs.contract_marker)
            Test.@test !ispath(completion_race_reservation.path)
            Test.@test isdir(completion_race_reservation.runtime_path)
            Test.@test isdir(
                completion_race_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test isdir(
                Mycelia._metamdbg_output_lock_path(completion_race_outdir),
            )

            pending_runtime_outdir = joinpath(
                temporary_root,
                "executor-pending-bind-cleanup",
            )
            pending_runtime_outputs, pending_runtime_reservation =
                create_bound_runtime_fixture(
                    pending_runtime_outdir,
                    "executor-pending-bind-cleanup-owner",
                    "883",
                )
            pending_runtime_path =
                Mycelia._metamdbg_pending_submission_job_path(
                    pending_runtime_outdir,
                    pending_runtime_reservation.output_root_reservation_capability,
                    pending_runtime_reservation.job_id,
                )
            Base.Filesystem.hardlink(
                pending_runtime_reservation.job_marker,
                pending_runtime_path,
            )
            Mycelia._fsync_metamdbg_directory(dirname(pending_runtime_path))
            pending_runtime_assembly_marker = joinpath(
                temporary_root,
                "executor-pending-bind-cleanup-assembly-ran",
            )
            pending_runtime_script = Mycelia._metamdbg_executor_script(
                "touch " * Base.shell_escape(pending_runtime_assembly_marker) *
                "; " * fake_asm,
                fake_gfa,
                pending_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = pending_runtime_reservation,
            )
            pending_runtime_script_path = joinpath(
                temporary_root,
                "executor-pending-bind-cleanup.sh",
            )
            write(pending_runtime_script_path, pending_runtime_script)
            pending_runtime_log = IOBuffer()
            pending_runtime_process = Base.withenv(
                "SLURM_JOB_ID" => pending_runtime_reservation.job_id,
            ) do
                return run(pipeline(
                    ignorestatus(`bash $(pending_runtime_script_path)`),
                    stdout = pending_runtime_log,
                    stderr = pending_runtime_log,
                ))
            end
            Test.@test !success(pending_runtime_process)
            Test.@test occursin(
                "runtime refuses pending scheduler job evidence",
                String(take!(pending_runtime_log)),
            )
            Test.@test isdir(pending_runtime_reservation.path)
            Test.@test isfile(
                pending_runtime_reservation.output_root_reservation_marker,
            )
            Test.@test isfile(pending_runtime_path)
            Test.@test Base.Filesystem.samefile(
                pending_runtime_reservation.job_marker,
                pending_runtime_path,
            )
            Test.@test !ispath(pending_runtime_reservation.runtime_path)
            Test.@test !ispath(
                pending_runtime_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(pending_runtime_outdir),
            )
            Test.@test !ispath(pending_runtime_assembly_marker)
            pending_runtime_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    pending_runtime_outdir,
                ),
            )
            Test.@test pending_runtime_metadata.submission_state ==
                       :submission_commit_cleanup_pending
            pending_runtime_bound =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    pending_runtime_metadata;
                    owner_token = pending_runtime_metadata.owner_token,
                    job_id = pending_runtime_metadata.job_id,
                    confirm_submitted = true,
                )
            pending_runtime_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    pending_runtime_bound;
                    owner_token = pending_runtime_bound.owner_token,
                    job_id = pending_runtime_bound.job_id,
                    confirm_cancelled = true,
                )
            Test.@test pending_runtime_reclaimed.recovery_reason == :cancelled

            foreign_pending_outdir = joinpath(
                temporary_root,
                "executor-foreign-pending-same-output",
            )
            foreign_pending_outputs, foreign_pending_reservation =
                create_bound_runtime_fixture(
                    foreign_pending_outdir,
                    "executor-foreign-pending-same-output-owner",
                    "884",
                )
            foreign_capability = repeat("a", 64)
            if foreign_capability ==
               foreign_pending_reservation.output_root_reservation_capability
                foreign_capability = repeat("b", 64)
            end
            foreign_pending_path =
                Mycelia._metamdbg_pending_submission_job_path(
                    foreign_pending_outdir,
                    foreign_capability,
                    "999",
                )
            open(foreign_pending_path, "w") do pending_io
                write(pending_io, "foreign accepted-job evidence\n")
                chmod(foreign_pending_path, 0o600)
                Mycelia._fsync_metamdbg_file(
                    pending_io,
                    foreign_pending_path,
                )
            end
            Mycelia._fsync_metamdbg_directory(dirname(foreign_pending_path))
            foreign_pending_assembly_marker = joinpath(
                temporary_root,
                "executor-foreign-pending-same-output-assembly-ran",
            )
            foreign_pending_script = Mycelia._metamdbg_executor_script(
                "touch " * Base.shell_escape(foreign_pending_assembly_marker) *
                "; " * fake_asm,
                fake_gfa,
                foreign_pending_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = foreign_pending_reservation,
            )
            foreign_pending_script_path = joinpath(
                temporary_root,
                "executor-foreign-pending-same-output.sh",
            )
            write(foreign_pending_script_path, foreign_pending_script)
            foreign_pending_log = IOBuffer()
            foreign_pending_process = Base.withenv(
                "SLURM_JOB_ID" => foreign_pending_reservation.job_id,
            ) do
                return run(pipeline(
                    ignorestatus(`bash $(foreign_pending_script_path)`),
                    stdout = foreign_pending_log,
                    stderr = foreign_pending_log,
                ))
            end
            Test.@test !success(foreign_pending_process)
            Test.@test occursin(
                "runtime refuses ambiguous pending scheduler job evidence",
                String(take!(foreign_pending_log)),
            )
            Test.@test isdir(foreign_pending_reservation.path)
            Test.@test isfile(
                foreign_pending_reservation.output_root_reservation_marker,
            )
            Test.@test isfile(foreign_pending_path)
            Test.@test !ispath(foreign_pending_reservation.runtime_path)
            Test.@test !ispath(
                foreign_pending_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(foreign_pending_outdir),
            )
            Test.@test !ispath(foreign_pending_assembly_marker)
            rm(foreign_pending_path)
            Mycelia._fsync_metamdbg_directory(dirname(foreign_pending_path))
            foreign_pending_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    foreign_pending_outdir,
                ),
            )
            foreign_pending_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    foreign_pending_metadata;
                    owner_token = foreign_pending_metadata.owner_token,
                    job_id = foreign_pending_metadata.job_id,
                    confirm_cancelled = true,
                )
            Test.@test foreign_pending_reclaimed.recovery_reason == :cancelled

            ignored_pending_parent = joinpath(
                temporary_root,
                "executor-ignores-other-output-pending-parent",
            )
            mkpath(ignored_pending_parent)
            ignored_pending_target_outdir = joinpath(
                ignored_pending_parent,
                "executor-ignores-other-output-pending-target",
            )
            ignored_pending_outputs, ignored_pending_reservation =
                create_bound_runtime_fixture(
                    ignored_pending_target_outdir,
                    "executor-ignores-other-output-pending-target-owner",
                    "885",
                )
            ignored_pending_other_outdir = joinpath(
                ignored_pending_parent,
                "executor-ignores-other-output-pending-other",
            )
            ignored_pending_other_outputs = Mycelia._metamdbg_output_paths(
                ignored_pending_other_outdir,
                21,
            )
            ignored_pending_other_reservation =
                Mycelia._metamdbg_submission_reservation(
                    ignored_pending_other_outputs,
                    input_contract,
                    21;
                    owner_token =
                        "executor-ignores-other-output-pending-other-owner",
                )
            Mycelia._with_metamdbg_output_lock(
                ignored_pending_other_outdir,
            ) do
                Mycelia._create_metamdbg_submission_reservation!(
                    ignored_pending_other_reservation,
                    ignored_pending_other_outdir,
                )
            end
            _test_metamdbg_error(
                () -> Mycelia._with_metamdbg_output_domain_lock(
                    ignored_pending_other_outdir;
                    allowed_same_root_locks = (
                        ignored_pending_other_reservation.output_root_reservation_marker,
                    ),
                ) do
                    Mycelia._bind_metamdbg_submission_job!(
                        ignored_pending_other_reservation,
                        "886";
                        pre_job_record_publication_hook = (
                            _reservation::NamedTuple,
                            _path::AbstractString,
                        ) -> error("synthetic other-output bind failure"),
                    )
                end,
                ErrorException,
                r"synthetic other-output bind failure",
            )
            ignored_pending_other_path =
                Mycelia._metamdbg_pending_submission_job_path(
                    ignored_pending_other_outdir,
                    ignored_pending_other_reservation.output_root_reservation_capability,
                    "886",
                )
            Test.@test isfile(ignored_pending_other_path)
            ignored_pending_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                ignored_pending_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = ignored_pending_reservation,
            )
            ignored_pending_script_path = joinpath(
                temporary_root,
                "executor-ignores-other-output-pending-target.sh",
            )
            write(ignored_pending_script_path, ignored_pending_script)
            Test.@test Base.withenv(
                "SLURM_JOB_ID" => ignored_pending_reservation.job_id,
            ) do
                success(`bash $(ignored_pending_script_path)`)
            end
            Test.@test isfile(ignored_pending_outputs.completion_marker)
            Test.@test !ispath(ignored_pending_reservation.path)
            Test.@test isdir(ignored_pending_reservation.consumed_path)
            Test.@test isfile(ignored_pending_other_path)
            ignored_pending_target_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    ignored_pending_target_outdir,
                ),
            )
            ignored_pending_target_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    ignored_pending_target_metadata;
                    owner_token = ignored_pending_target_metadata.owner_token,
                    job_id = ignored_pending_target_metadata.job_id,
                    confirm_terminal = :completed,
                )
            Test.@test ignored_pending_target_reclaimed.recovery_reason ==
                       :completed
            ignored_pending_other_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    ignored_pending_other_outdir,
                ),
            )
            ignored_pending_other_bound =
                Mycelia.bind_metamdbg_submission_reservation_job!(
                    ignored_pending_other_metadata;
                    owner_token = ignored_pending_other_metadata.owner_token,
                    job_id = "886",
                    confirm_submitted = true,
                )
            Test.@test !ispath(ignored_pending_other_path)
            ignored_pending_other_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    ignored_pending_other_bound;
                    owner_token = ignored_pending_other_bound.owner_token,
                    job_id = ignored_pending_other_bound.job_id,
                    confirm_cancelled = true,
                )
            Test.@test ignored_pending_other_reclaimed.recovery_reason ==
                       :cancelled
            Test.@test !ispath(ignored_pending_other_reservation.path)
            Test.@test !ispath(
                ignored_pending_other_reservation.output_root_reservation_marker,
            )

            multi_snapshot_outdir =
                joinpath(temporary_root, "multi-record-inspection-race")
            _, multi_consumed_reservation = create_bound_runtime_fixture(
                multi_snapshot_outdir,
                "multi-record-consumed-owner",
                "881",
            )
            Mycelia._with_metamdbg_output_domain_lock(
                multi_snapshot_outdir;
                allowed_same_root_locks = (
                    multi_consumed_reservation.output_root_reservation_marker,
                ),
            ) do
                Mycelia._require_metamdbg_submission_reservation!(
                    multi_consumed_reservation,
                )
                rm(multi_consumed_reservation.output_root_reservation_marker)
                Mycelia._fsync_metamdbg_directory(
                    dirname(
                        multi_consumed_reservation.output_root_reservation_marker,
                    ),
                )
                mv(
                    multi_consumed_reservation.path,
                    multi_consumed_reservation.consumed_path,
                )
                Mycelia._fsync_metamdbg_directory(
                    dirname(multi_consumed_reservation.path),
                )
            end
            multi_queued_outputs =
                Mycelia._metamdbg_output_paths(multi_snapshot_outdir, 21)
            multi_queued_reservation =
                Mycelia._metamdbg_submission_reservation(
                    multi_queued_outputs,
                    input_contract,
                    21;
                    owner_token = "multi-record-queued-owner",
                )
            Mycelia._with_metamdbg_output_domain_lock(
                multi_snapshot_outdir,
            ) do
                Mycelia._create_metamdbg_submission_reservation!(
                    multi_queued_reservation,
                    multi_snapshot_outdir,
                )
            end
            multi_snapshot_hook_calls = Ref(0)
            multi_snapshot_hook = function (snapshot::NamedTuple)
                multi_snapshot_hook_calls[] += 1
                Test.@test length(snapshot.reservations) == 2
                Test.@test !ispath(
                    Mycelia._output_root_reservation_lock_path_from_canonical(
                        multi_snapshot_outdir,
                    ),
                )
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(multi_snapshot_outdir),
                )
                Test.@test !ispath(
                    Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                        multi_snapshot_outdir,
                    ),
                )
                Mycelia._bind_metamdbg_submission_job_after_submit!(
                    multi_queued_reservation,
                    "882",
                )
                return nothing
            end
            _test_metamdbg_error(
                () -> Mycelia._inspect_metamdbg_submission_reservations(
                    multi_snapshot_outdir;
                    post_initial_snapshot_hook = multi_snapshot_hook,
                ),
                ErrorException,
                r"reservation job_id changed during recovery inspection",
            )
            Test.@test multi_snapshot_hook_calls[] == 1
            multi_snapshot_records =
                Mycelia.inspect_metamdbg_submission_reservations(
                    multi_snapshot_outdir,
                )
            Test.@test length(multi_snapshot_records) == 2
            multi_queued_metadata = only(filter(
                record -> record.reservation_state == :queued,
                multi_snapshot_records,
            ))
            multi_consumed_metadata = only(filter(
                record -> record.reservation_state == :consumed,
                multi_snapshot_records,
            ))
            Test.@test multi_queued_metadata.job_id == "882"
            Mycelia.reclaim_metamdbg_submission_reservation!(
                multi_queued_metadata;
                owner_token = multi_queued_metadata.owner_token,
                job_id = multi_queued_metadata.job_id,
                confirm_cancelled = true,
            )
            Mycelia.reclaim_metamdbg_submission_reservation!(
                multi_consumed_metadata;
                owner_token = multi_consumed_metadata.owner_token,
                job_id = multi_consumed_metadata.job_id,
                confirm_terminal = :completed,
            )

            inspection_race_outdir =
                joinpath(temporary_root, "inspection-runtime-race")
            inspection_race_outputs, inspection_race_reservation =
                create_bound_runtime_fixture(
                    inspection_race_outdir,
                    "inspection-runtime-race-owner",
                    "880",
                )
            inspection_race_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                inspection_race_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = inspection_race_reservation,
            )
            inspection_race_script_path = joinpath(
                temporary_root,
                "inspection-runtime-race.sh",
            )
            write(inspection_race_script_path, inspection_race_script)
            inspection_hook_calls = Ref(0)
            inspection_runtime_process = Ref{Union{Nothing, Base.Process}}(
                nothing,
            )
            inspection_race_hook = function (_reservation::NamedTuple)
                inspection_hook_calls[] += 1
                pid_lock =
                    Mycelia._output_root_reservation_lock_path_from_canonical(
                        inspection_race_outdir,
                    )
                private_lock =
                    Mycelia._metamdbg_output_lock_path(inspection_race_outdir)
                cleanup_reservation =
                    Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                        inspection_race_outdir,
                    )
                Test.@test !ispath(pid_lock)
                Test.@test !ispath(private_lock)
                Test.@test !ispath(cleanup_reservation)
                process = run(
                    addenv(
                        `bash $(inspection_race_script_path)`,
                        "SLURM_JOB_ID" => inspection_race_reservation.job_id,
                    );
                    wait = false,
                )
                inspection_runtime_process[] = process
                wait(process)
                Test.@test success(process)
                Test.@test !ispath(pid_lock)
                Test.@test !ispath(private_lock)
                Test.@test !ispath(cleanup_reservation)
                return nothing
            end
            _test_metamdbg_error(
                () -> Mycelia._inspect_metamdbg_submission_reservations(
                    inspection_race_outdir;
                    post_initial_snapshot_hook = inspection_race_hook,
                ),
                ErrorException,
                Regex(
                    "changed during recovery inspection|" *
                    "presence changed after recovery inspection|" *
                    "must be a regular, non-symlink directory",
                ),
            )
            Test.@test inspection_hook_calls[] == 1
            Test.@test inspection_runtime_process[] isa Base.Process
            if inspection_runtime_process[] isa Base.Process
                Test.@test success(inspection_runtime_process[])
            end
            Test.@test !ispath(inspection_race_reservation.path)
            Test.@test !ispath(
                inspection_race_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                inspection_race_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test isdir(inspection_race_reservation.consumed_path)
            inspection_race_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    inspection_race_outdir,
                ),
            )
            Test.@test inspection_race_metadata.reservation_state == :consumed
            reclaimed_inspection_race =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    inspection_race_metadata;
                    owner_token = inspection_race_metadata.owner_token,
                    job_id = inspection_race_metadata.job_id,
                    confirm_terminal = :completed,
                )
            Test.@test reclaimed_inspection_race.recovery_reason == :completed

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

            for (hierarchy_index, hierarchy_case) in
                enumerate((:ancestor, :descendant))
                hierarchy_parent = joinpath(
                    temporary_root,
                    "executor-runtime-$(hierarchy_case)-domain",
                )
                hierarchy_outdir = if hierarchy_case == :ancestor
                    joinpath(
                        hierarchy_parent,
                        "reserved-ancestor",
                        "nested",
                        "output",
                    )
                else
                    joinpath(hierarchy_parent, "output")
                end
                hierarchy_foreign_root = if hierarchy_case == :ancestor
                    joinpath(hierarchy_parent, "reserved-ancestor")
                else
                    joinpath(hierarchy_outdir, "reserved-descendant")
                end
                hierarchy_outputs, hierarchy_reservation =
                    create_bound_runtime_fixture(
                        hierarchy_outdir,
                        "runtime-$(hierarchy_case)-owner",
                        string(870 + hierarchy_index),
                    )
                mkpath(dirname(hierarchy_foreign_root))
                hierarchy_case == :descendant &&
                    mkpath(hierarchy_foreign_root)
                hierarchy_foreign_marker =
                    Mycelia._output_root_durable_reservation_path_from_canonical(
                        hierarchy_foreign_root,
                        "runtime-$(hierarchy_case)-foreign-owner",
                    )
                write(hierarchy_foreign_marker, "foreign owner\n")
                chmod(hierarchy_foreign_marker, 0o600)
                hierarchy_assembly_marker = joinpath(
                    temporary_root,
                    "runtime-$(hierarchy_case)-assembly-ran",
                )
                hierarchy_script = Mycelia._metamdbg_executor_script(
                    "touch $(Base.shell_escape(hierarchy_assembly_marker)); " *
                    fake_asm,
                    fake_gfa,
                    hierarchy_outputs,
                    21,
                    input_contract;
                    conda_runner = fake_conda,
                    submission_reservation = hierarchy_reservation,
                )
                hierarchy_script_path = joinpath(
                    temporary_root,
                    "executor-runtime-$(hierarchy_case).sh",
                )
                write(hierarchy_script_path, hierarchy_script)
                Test.@test !Base.withenv(
                    "SLURM_JOB_ID" => hierarchy_reservation.job_id,
                ) do
                    success(`bash $(hierarchy_script_path)`)
                end
                Test.@test !ispath(hierarchy_assembly_marker)
                Test.@test isdir(hierarchy_reservation.path)
                Test.@test isfile(
                    hierarchy_reservation.output_root_reservation_marker,
                )
                Test.@test !ispath(
                    hierarchy_reservation.runtime_output_root_reservation_marker,
                )
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(hierarchy_outdir),
                )
                Test.@test !ispath(
                    Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                        hierarchy_outdir,
                    ),
                )
                rm(hierarchy_foreign_marker)
                remove_queued_runtime_fixture!(
                    hierarchy_outdir,
                    hierarchy_reservation,
                )
            end

            sibling_runtime_parent =
                joinpath(temporary_root, "executor-runtime-sibling-domain")
            sibling_runtime_outdir =
                joinpath(sibling_runtime_parent, "output")
            sibling_foreign_root =
                joinpath(sibling_runtime_parent, "foreign-sibling")
            sibling_runtime_outputs, sibling_runtime_reservation =
                create_bound_runtime_fixture(
                    sibling_runtime_outdir,
                    "runtime-sibling-owner",
                    "873",
                )
            sibling_foreign_marker =
                Mycelia._output_root_durable_reservation_path_from_canonical(
                    sibling_foreign_root,
                    "runtime-sibling-foreign-owner",
                )
            write(sibling_foreign_marker, "foreign sibling owner\n")
            chmod(sibling_foreign_marker, 0o600)
            sibling_assembly_marker = joinpath(
                temporary_root,
                "runtime-sibling-assembly-ran",
            )
            sibling_runtime_script = Mycelia._metamdbg_executor_script(
                "touch $(Base.shell_escape(sibling_assembly_marker)); " *
                fake_asm,
                fake_gfa,
                sibling_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = sibling_runtime_reservation,
            )
            sibling_runtime_script_path =
                joinpath(temporary_root, "executor-runtime-sibling.sh")
            write(sibling_runtime_script_path, sibling_runtime_script)
            Test.@test Base.withenv(
                "SLURM_JOB_ID" => sibling_runtime_reservation.job_id,
            ) do
                success(`bash $(sibling_runtime_script_path)`)
            end
            Test.@test ispath(sibling_assembly_marker)
            Test.@test !ispath(sibling_runtime_reservation.path)
            Test.@test !ispath(
                sibling_runtime_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                sibling_runtime_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(sibling_runtime_outdir),
            )
            Test.@test isfile(sibling_foreign_marker)
            rm(sibling_foreign_marker)

            unreadable_runtime_parent =
                joinpath(temporary_root, "executor-unreadable-ancestor")
            unreadable_runtime_outdir =
                joinpath(unreadable_runtime_parent, "output")
            unreadable_outputs, unreadable_reservation =
                create_bound_runtime_fixture(
                    unreadable_runtime_outdir,
                    "runtime-unreadable-ancestor-owner",
                    "874",
                )
            unreadable_assembly_marker = joinpath(
                temporary_root,
                "runtime-unreadable-ancestor-assembly-ran",
            )
            unreadable_script = Mycelia._metamdbg_executor_script(
                "touch $(Base.shell_escape(unreadable_assembly_marker)); " *
                fake_asm,
                fake_gfa,
                unreadable_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = unreadable_reservation,
            )
            unreadable_script_path = joinpath(
                temporary_root,
                "executor-unreadable-ancestor.sh",
            )
            write(unreadable_script_path, unreadable_script)
            unreadable_runtime_log = IOBuffer()
            unreadable_runtime_process = try
                chmod(unreadable_runtime_parent, 0o333)
                Base.withenv(
                    "SLURM_JOB_ID" => unreadable_reservation.job_id,
                ) do
                    run(pipeline(
                        ignorestatus(`bash $(unreadable_script_path)`);
                        stdout = unreadable_runtime_log,
                        stderr = unreadable_runtime_log,
                    ))
                end
            finally
                chmod(unreadable_runtime_parent, 0o700)
            end
            unreadable_runtime_output =
                String(take!(unreadable_runtime_log))
            Test.@test !success(unreadable_runtime_process)
            Test.@test occursin(
                "could not durably fsync lifecycle path and containing " *
                "directory",
                unreadable_runtime_output,
            )
            Test.@test !ispath(unreadable_assembly_marker)
            Test.@test isdir(unreadable_reservation.path)
            Test.@test isfile(
                unreadable_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                unreadable_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(unreadable_runtime_outdir),
            )
            remove_queued_runtime_fixture!(
                unreadable_runtime_outdir,
                unreadable_reservation,
            )

            function enumeration_retry_remnants(
                    root::AbstractString,
            )::Vector{String}
                remnants = String[]
                for (directory, directories, files) in walkdir(root)
                    for name in vcat(directories, files)
                        if startswith(name, ".mycelia-metamdbg-tmp.") ||
                           endswith(name, ".attempt") ||
                           endswith(name, ".errors")
                            push!(remnants, joinpath(directory, name))
                        end
                    end
                end
                return remnants
            end

            aba_root = joinpath(temporary_root, "enumerator-aba-root")
            aba_replacement = joinpath(
                temporary_root,
                "enumerator-aba-replacement",
            )
            aba_parked = joinpath(temporary_root, "enumerator-aba-parked")
            mkpath(aba_root)
            mkpath(aba_replacement)
            aba_original_entry = joinpath(aba_root, "original-entry")
            aba_replacement_entry = joinpath(
                aba_replacement,
                "replacement-entry",
            )
            write(aba_original_entry, "original\n")
            aba_newline_entry = joinpath(aba_root, "original\nnewline-entry")
            write(aba_newline_entry, "newline\n")
            write(aba_replacement_entry, "replacement\n")
            aba_open_ready = joinpath(temporary_root, "aba-open.ready")
            aba_open_continue = joinpath(temporary_root, "aba-open.continue")
            aba_enumeration_ready = joinpath(
                temporary_root,
                "aba-enumeration.ready",
            )
            aba_enumeration_continue = joinpath(
                temporary_root,
                "aba-enumeration.continue",
            )
            aba_inventory = joinpath(temporary_root, "aba-inventory.paths")
            aba_open_hook = """
            with open($(JSON.json(aba_open_ready)), "w", encoding="utf-8") as stream:
                stream.write("ready\\n")
            while not os.path.exists($(JSON.json(aba_open_continue))):
                time.sleep(0.01)
            """
            aba_enumeration_hook = """
            with open($(JSON.json(aba_enumeration_ready)), "w", encoding="utf-8") as stream:
                stream.write("ready\\n")
            while not os.path.exists($(JSON.json(aba_enumeration_continue))):
                time.sleep(0.01)
            """
            aba_enumerator =
                Mycelia._metamdbg_runtime_directory_enumerator_python(
                    ;
                    post_open_hook = aba_open_hook,
                    post_enumeration_hook = aba_enumeration_hook,
                )
            python_path = Sys.which("python3")
            python_path === nothing && error(
                "python3 is required for the metaMDBG ABA enumerator test.",
            )
            depth_root = joinpath(
                temporary_root,
                "enumerator-depth-bound-root",
            )
            deep_directory = depth_root
            for _depth in 1:Mycelia._METAMDBG_ENUMERATION_MAX_DEPTH
                deep_directory = joinpath(deep_directory, "d")
                mkpath(deep_directory)
            end
            depth_marker = joinpath(
                deep_directory,
                Mycelia._OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX * "too-deep",
            )
            mkpath(depth_marker)
            depth_inventory = joinpath(
                temporary_root,
                "enumerator-depth-bound.paths",
            )
            depth_enumerator =
                Mycelia._metamdbg_runtime_directory_enumerator_python()
            depth_log = IOBuffer()
            depth_process = run(pipeline(
                ignorestatus(
                    `$(python_path) -c $(depth_enumerator) $(depth_root) $(depth_inventory) reservations-recursive $(Mycelia._OUTPUT_ROOT_RESERVATION_LOCK_PREFIX) $(Mycelia._OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX)`,
                );
                stdout = depth_log,
                stderr = depth_log,
            ))
            Test.@test !success(depth_process)
            Test.@test occursin(
                "directory manifest exceeds depth bound " *
                string(Mycelia._METAMDBG_ENUMERATION_MAX_DEPTH),
                String(take!(depth_log)),
            )
            Test.@test !ispath(depth_inventory)
            aba_process = run(
                `$(python_path) -c $(aba_enumerator) $(aba_root) $(aba_inventory) one-level "" ""`;
                wait = false,
            )
            Test.@test Base.timedwait(
                () -> isfile(aba_open_ready) || !Base.process_running(aba_process),
                30.0;
                pollint = 0.01,
            ) == :ok
            Test.@test Base.process_running(aba_process)
            mv(aba_root, aba_parked)
            mv(aba_replacement, aba_root)
            write(aba_open_continue, "continue\n")
            Test.@test Base.timedwait(
                () -> isfile(aba_enumeration_ready) ||
                      !Base.process_running(aba_process),
                30.0;
                pollint = 0.01,
            ) == :ok
            Test.@test Base.process_running(aba_process)
            mv(aba_root, aba_replacement)
            mv(aba_parked, aba_root)
            write(aba_enumeration_continue, "continue\n")
            wait(aba_process)
            Test.@test success(aba_process)
            aba_entries = filter(
                !isempty,
                split(String(read(aba_inventory)), '\0'),
            )
            Test.@test length(aba_entries) == 2
            Test.@test Set(aba_entries) == Set((
                aba_original_entry,
                aba_newline_entry,
            ))
            Test.@test isfile(aba_original_entry)
            Test.@test isfile(aba_newline_entry)
            Test.@test isfile(aba_replacement_entry)

            before_inventory_root = joinpath(
                temporary_root,
                "enumerator-before-inventory-root",
            )
            before_inventory_replacement = joinpath(
                temporary_root,
                "enumerator-before-inventory-replacement",
            )
            before_inventory_parked = joinpath(
                temporary_root,
                "enumerator-before-inventory-parked",
            )
            mkpath(before_inventory_root)
            mkpath(before_inventory_replacement)
            write(joinpath(before_inventory_root, "original"), "original\n")
            write(
                joinpath(before_inventory_replacement, "replacement"),
                "replacement\n",
            )
            before_inventory_ready = joinpath(
                temporary_root,
                "before-inventory.ready",
            )
            before_inventory_continue = joinpath(
                temporary_root,
                "before-inventory.continue",
            )
            before_inventory_path = joinpath(
                temporary_root,
                "before-inventory.paths",
            )
            before_inventory_hook = """
            with open($(JSON.json(before_inventory_ready)), "w", encoding="utf-8") as stream:
                stream.write("ready\\n")
            while not os.path.exists($(JSON.json(before_inventory_continue))):
                time.sleep(0.01)
            """
            before_inventory_enumerator =
                Mycelia._metamdbg_runtime_directory_enumerator_python(
                    ;
                    post_enumeration_hook = before_inventory_hook,
                )
            before_inventory_process = run(
                ignorestatus(
                    `$(python_path) -c $(before_inventory_enumerator) $(before_inventory_root) $(before_inventory_path) one-level "" ""`,
                );
                wait = false,
            )
            Test.@test Base.timedwait(
                () -> isfile(before_inventory_ready) ||
                      !Base.process_running(before_inventory_process),
                30.0;
                pollint = 0.01,
            ) == :ok
            Test.@test Base.process_running(before_inventory_process)
            mv(before_inventory_root, before_inventory_parked)
            mv(before_inventory_replacement, before_inventory_root)
            write(before_inventory_continue, "continue\n")
            wait(before_inventory_process)
            Test.@test !success(before_inventory_process)
            Test.@test !ispath(before_inventory_path)
            mv(before_inventory_root, before_inventory_replacement)
            mv(before_inventory_parked, before_inventory_root)
            Test.@test isfile(joinpath(before_inventory_root, "original"))
            Test.@test isfile(joinpath(
                before_inventory_replacement,
                "replacement",
            ))

            after_inventory_root = joinpath(
                temporary_root,
                "enumerator-after-inventory-root",
            )
            after_inventory_replacement = joinpath(
                temporary_root,
                "enumerator-after-inventory-replacement",
            )
            after_inventory_parked = joinpath(
                temporary_root,
                "enumerator-after-inventory-parked",
            )
            mkpath(after_inventory_root)
            mkpath(after_inventory_replacement)
            write(joinpath(after_inventory_root, "original"), "original\n")
            write(
                joinpath(after_inventory_replacement, "replacement"),
                "replacement\n",
            )
            after_inventory_ready = joinpath(
                temporary_root,
                "after-inventory.ready",
            )
            after_inventory_continue = joinpath(
                temporary_root,
                "after-inventory.continue",
            )
            after_inventory_path = joinpath(
                temporary_root,
                "after-inventory.paths",
            )
            after_inventory_hook = """
            with open($(JSON.json(after_inventory_ready)), "w", encoding="utf-8") as stream:
                stream.write("ready\\n")
            while not os.path.exists($(JSON.json(after_inventory_continue))):
                time.sleep(0.01)
            """
            after_inventory_enumerator =
                Mycelia._metamdbg_runtime_directory_enumerator_python(
                    ;
                    post_inventory_hook = after_inventory_hook,
                )
            after_inventory_process = run(
                ignorestatus(
                    `$(python_path) -c $(after_inventory_enumerator) $(after_inventory_root) $(after_inventory_path) one-level "" ""`,
                );
                wait = false,
            )
            Test.@test Base.timedwait(
                () -> isfile(after_inventory_ready) ||
                      !Base.process_running(after_inventory_process),
                30.0;
                pollint = 0.01,
            ) == :ok
            Test.@test Base.process_running(after_inventory_process)
            mv(after_inventory_root, after_inventory_parked)
            mv(after_inventory_replacement, after_inventory_root)
            write(after_inventory_continue, "continue\n")
            wait(after_inventory_process)
            Test.@test !success(after_inventory_process)
            Test.@test !ispath(after_inventory_path)
            mv(after_inventory_root, after_inventory_replacement)
            mv(after_inventory_parked, after_inventory_root)
            Test.@test isfile(joinpath(after_inventory_root, "original"))
            Test.@test isfile(joinpath(
                after_inventory_replacement,
                "replacement",
            ))

            function require_reservation_manifest_mutation_result!(
                    label::String,
                    setup::Function,
                    mutation::Function,
                    ;
                    expect_rejection::Bool,
                    enumeration_mode::String = "reservations-recursive",
            )::Nothing
                mutation_root = joinpath(
                    temporary_root,
                    "manifest-mutation-$(label)",
                )
                mkpath(mutation_root)
                setup(mutation_root)
                mutation_ready = joinpath(
                    temporary_root,
                    "manifest-mutation-$(label).ready",
                )
                mutation_continue = joinpath(
                    temporary_root,
                    "manifest-mutation-$(label).continue",
                )
                mutation_inventory = joinpath(
                    temporary_root,
                    "manifest-mutation-$(label).paths",
                )
                mutation_log = joinpath(
                    temporary_root,
                    "manifest-mutation-$(label).log",
                )
                mutation_hook = """
                with open($(JSON.json(mutation_ready)), "w", encoding="utf-8") as stream:
                    stream.write("ready\\n")
                while not os.path.exists($(JSON.json(mutation_continue))):
                    time.sleep(0.01)
                """
                mutation_enumerator =
                    Mycelia._metamdbg_runtime_directory_enumerator_python(
                        ;
                        post_inventory_hook = mutation_hook,
                    )
                mutation_command = if enumeration_mode == "reservations-recursive"
                    `$(python_path) -c $(mutation_enumerator) $(mutation_root) $(mutation_inventory) reservations-recursive $(Mycelia._OUTPUT_ROOT_RESERVATION_LOCK_PREFIX) $(Mycelia._OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX)`
                elseif enumeration_mode == "reservations-one-level"
                    `$(python_path) -c $(mutation_enumerator) $(mutation_root) $(mutation_inventory) reservations-one-level $(Mycelia._OUTPUT_ROOT_RESERVATION_LOCK_PREFIX)exact $(Mycelia._OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX) .mycelia-metamdbg-pending.fixture.`
                else
                    error(
                        "unsupported reservation manifest test mode: " *
                        enumeration_mode,
                    )
                end
                mutation_log_io = open(mutation_log, "w")
                mutation_process = run(
                    pipeline(
                        ignorestatus(mutation_command);
                        stdout = mutation_log_io,
                        stderr = mutation_log_io,
                    );
                    wait = false,
                )
                Test.@test Base.timedwait(
                    () -> isfile(mutation_ready) ||
                          !Base.process_running(mutation_process),
                    30.0;
                    pollint = 0.01,
                ) == :ok
                Test.@test Base.process_running(mutation_process)
                mutation(mutation_root)
                write(mutation_continue, "continue\n")
                wait(mutation_process)
                close(mutation_log_io)
                if expect_rejection
                    Test.@test !success(mutation_process)
                    Test.@test occursin(
                        "directory manifest changed before inventory " *
                        "publication",
                        read(mutation_log, String),
                    )
                    Test.@test !ispath(mutation_inventory)
                else
                    Test.@test success(mutation_process)
                    Test.@test isfile(mutation_inventory)
                end
                return nothing
            end

            require_reservation_manifest_mutation_result!(
                "unrelated-insert",
                root -> begin
                    mkpath(joinpath(root, "nested"))
                    write(joinpath(root, "nested", "stable"), "stable\n")
                end,
                root -> write(
                    joinpath(root, "nested", "inserted"),
                    "inserted\n",
                ),
                ;
                expect_rejection = false,
            )
            require_reservation_manifest_mutation_result!(
                "unrelated-remove",
                root -> begin
                    mkpath(joinpath(root, "nested"))
                    write(joinpath(root, "nested", "removed"), "removed\n")
                end,
                root -> rm(joinpath(root, "nested", "removed")),
                ;
                expect_rejection = false,
            )
            require_reservation_manifest_mutation_result!(
                "unrelated-file-to-directory",
                root -> write(joinpath(root, "changing"), "file\n"),
                root -> begin
                    rm(joinpath(root, "changing"))
                    mkpath(joinpath(root, "changing"))
                    write(
                        joinpath(root, "changing", "descendant"),
                        "descendant\n",
                    )
                end,
                ;
                expect_rejection = false,
            )
            require_reservation_manifest_mutation_result!(
                "unrelated-symlink-to-directory",
                root -> begin
                    write(joinpath(root, "symlink-target"), "target\n")
                    symlink("symlink-target", joinpath(root, "changing"))
                end,
                root -> begin
                    rm(joinpath(root, "changing"))
                    mkpath(joinpath(root, "changing"))
                end,
                ;
                expect_rejection = false,
            )
            require_reservation_manifest_mutation_result!(
                "unrelated-directory-identity",
                root -> begin
                    mkpath(joinpath(root, "changing"))
                    write(joinpath(root, "changing", "stable"), "stable\n")
                end,
                root -> begin
                    mv(
                        joinpath(root, "changing"),
                        joinpath(temporary_root, "parked-changing-directory"),
                    )
                    mkpath(joinpath(root, "changing"))
                    write(joinpath(root, "changing", "stable"), "stable\n")
                end,
                ;
                expect_rejection = false,
            )

            durable_reservation_name(label::String)::String =
                Mycelia._OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX * label
            require_reservation_manifest_mutation_result!(
                "one-level-unrelated-insert",
                root -> write(joinpath(root, "stable"), "stable\n"),
                root -> write(joinpath(root, "inserted"), "inserted\n"),
                ;
                expect_rejection = false,
                enumeration_mode = "reservations-one-level",
            )
            require_reservation_manifest_mutation_result!(
                "one-level-reservation-insert",
                root -> write(joinpath(root, "stable"), "stable\n"),
                root -> mkpath(joinpath(
                    root,
                    durable_reservation_name("one-level-inserted"),
                )),
                ;
                expect_rejection = true,
                enumeration_mode = "reservations-one-level",
            )
            require_reservation_manifest_mutation_result!(
                "reservation-insert",
                root -> mkpath(joinpath(root, "nested")),
                root -> mkpath(joinpath(
                    root,
                    "nested",
                    durable_reservation_name("inserted"),
                )),
                ;
                expect_rejection = true,
            )
            require_reservation_manifest_mutation_result!(
                "reservation-remove",
                root -> begin
                    mkpath(joinpath(root, "nested"))
                    mkpath(joinpath(
                        root,
                        "nested",
                        durable_reservation_name("removed"),
                    ))
                end,
                root -> rm(joinpath(
                    root,
                    "nested",
                    durable_reservation_name("removed"),
                )),
                ;
                expect_rejection = true,
            )
            require_reservation_manifest_mutation_result!(
                "reservation-file-to-directory",
                root -> write(
                    joinpath(root, durable_reservation_name("changing")),
                    "file\n",
                ),
                root -> begin
                    changing = joinpath(
                        root,
                        durable_reservation_name("changing"),
                    )
                    rm(changing)
                    mkpath(changing)
                end,
                ;
                expect_rejection = true,
            )
            require_reservation_manifest_mutation_result!(
                "reservation-symlink-to-directory",
                root -> begin
                    write(joinpath(root, "symlink-target"), "target\n")
                    symlink(
                        "symlink-target",
                        joinpath(
                            root,
                            durable_reservation_name("changing"),
                        ),
                    )
                end,
                root -> begin
                    changing = joinpath(
                        root,
                        durable_reservation_name("changing"),
                    )
                    rm(changing)
                    mkpath(changing)
                end,
                ;
                expect_rejection = true,
            )
            require_reservation_manifest_mutation_result!(
                "reservation-directory-identity",
                root -> begin
                    changing = joinpath(
                        root,
                        durable_reservation_name("changing"),
                    )
                    mkpath(changing)
                    write(joinpath(changing, "stable"), "stable\n")
                end,
                root -> begin
                    changing = joinpath(
                        root,
                        durable_reservation_name("changing"),
                    )
                    mv(
                        changing,
                        joinpath(
                            temporary_root,
                            "parked-reservation-changing-directory",
                        ),
                    )
                    mkpath(changing)
                    write(joinpath(changing, "stable"), "stable\n")
                end,
                ;
                expect_rejection = true,
            )

            nested_aba_root = joinpath(
                temporary_root,
                "enumerator-nested-aba-root",
            )
            nested_aba_directory = joinpath(nested_aba_root, "nested")
            nested_aba_replacement = joinpath(
                temporary_root,
                "enumerator-nested-aba-replacement",
            )
            nested_aba_parked = joinpath(
                temporary_root,
                "enumerator-nested-aba-parked",
            )
            mkpath(nested_aba_directory)
            mkpath(nested_aba_replacement)
            nested_aba_marker = joinpath(
                nested_aba_directory,
                Mycelia._OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX * "nested",
            )
            write(nested_aba_marker, "nested marker\n")
            nested_aba_open_ready = joinpath(
                temporary_root,
                "nested-aba-open.ready",
            )
            nested_aba_open_continue = joinpath(
                temporary_root,
                "nested-aba-open.continue",
            )
            nested_aba_enumeration_ready = joinpath(
                temporary_root,
                "nested-aba-enumeration.ready",
            )
            nested_aba_enumeration_continue = joinpath(
                temporary_root,
                "nested-aba-enumeration.continue",
            )
            nested_aba_inventory = joinpath(
                temporary_root,
                "nested-aba.paths",
            )
            nested_aba_open_hook = """
            if manifest_generation == 1 and relative_parts == (b"nested",):
                with open($(JSON.json(nested_aba_open_ready)), "w", encoding="utf-8") as stream:
                    stream.write("ready\\n")
                while not os.path.exists($(JSON.json(nested_aba_open_continue))):
                    time.sleep(0.01)
            """
            nested_aba_enumeration_hook = """
            if manifest_generation == 1 and relative_parts == (b"nested",):
                with open($(JSON.json(nested_aba_enumeration_ready)), "w", encoding="utf-8") as stream:
                    stream.write("ready\\n")
                while not os.path.exists($(JSON.json(nested_aba_enumeration_continue))):
                    time.sleep(0.01)
            """
            nested_aba_enumerator =
                Mycelia._metamdbg_runtime_directory_enumerator_python(
                    ;
                    post_directory_open_hook = nested_aba_open_hook,
                    post_directory_enumeration_hook =
                        nested_aba_enumeration_hook,
                )
            nested_aba_process = run(
                `$(python_path) -c $(nested_aba_enumerator) $(nested_aba_root) $(nested_aba_inventory) reservations-recursive $(Mycelia._OUTPUT_ROOT_RESERVATION_LOCK_PREFIX) $(Mycelia._OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX)`;
                wait = false,
            )
            Test.@test Base.timedwait(
                () -> isfile(nested_aba_open_ready) ||
                      !Base.process_running(nested_aba_process),
                30.0;
                pollint = 0.01,
            ) == :ok
            Test.@test Base.process_running(nested_aba_process)
            mv(nested_aba_directory, nested_aba_parked)
            mv(nested_aba_replacement, nested_aba_directory)
            write(nested_aba_open_continue, "continue\n")
            Test.@test Base.timedwait(
                () -> isfile(nested_aba_enumeration_ready) ||
                      !Base.process_running(nested_aba_process),
                30.0;
                pollint = 0.01,
            ) == :ok
            Test.@test Base.process_running(nested_aba_process)
            mv(nested_aba_directory, nested_aba_replacement)
            mv(nested_aba_parked, nested_aba_directory)
            write(nested_aba_enumeration_continue, "continue\n")
            wait(nested_aba_process)
            Test.@test success(nested_aba_process)
            nested_aba_entries = filter(
                !isempty,
                split(String(read(nested_aba_inventory)), '\0'),
            )
            Test.@test nested_aba_entries == [nested_aba_marker]
            Test.@test isfile(nested_aba_marker)

            generation_two_fence_parent = joinpath(
                temporary_root,
                "executor-generation-two-fence",
            )
            generation_two_fence_outdir = joinpath(
                generation_two_fence_parent,
                "output",
            )
            generation_two_fence_outputs,
            generation_two_fence_reservation =
                create_bound_runtime_fixture(
                    generation_two_fence_outdir,
                    "runtime-generation-two-fence-owner",
                    "891",
                )
            generation_two_fence_assembly_marker = joinpath(
                temporary_root,
                "runtime-generation-two-fence-assembly-ran",
            )
            generation_two_fence_ready = joinpath(
                temporary_root,
                "generation-two-fence.ready",
            )
            generation_two_fence_continue = joinpath(
                temporary_root,
                "generation-two-fence.continue",
            )
            generation_two_fence_hook = """
            expected_root = os.fsencode(
                $(JSON.json(generation_two_fence_parent))
            )
            ready_path = $(JSON.json(generation_two_fence_ready))
            continue_path = $(JSON.json(generation_two_fence_continue))
            if (
                manifest_generation == 2
                and relative_parts == ()
                and root == expected_root
            ):
                with open(ready_path, "w", encoding="utf-8") as stream:
                    stream.write("ready\\n")
                while not os.path.exists(continue_path):
                    time.sleep(0.01)
            """
            generation_two_fence_script =
                Mycelia._metamdbg_executor_script(
                    "touch " *
                    Base.shell_escape(
                        generation_two_fence_assembly_marker,
                    ) *
                    "; " * fake_asm,
                    fake_gfa,
                    generation_two_fence_outputs,
                    21,
                    input_contract;
                    conda_runner = fake_conda,
                    submission_reservation =
                        generation_two_fence_reservation,
                    directory_enumerator_post_directory_open_hook =
                        generation_two_fence_hook,
                )
            generation_two_fence_script_path = joinpath(
                temporary_root,
                "executor-generation-two-fence.sh",
            )
            write(
                generation_two_fence_script_path,
                generation_two_fence_script,
            )
            generation_two_fence_process = run(
                addenv(
                    `bash $(generation_two_fence_script_path)`,
                    "SLURM_JOB_ID" =>
                        generation_two_fence_reservation.job_id,
                );
                wait = false,
            )
            generation_two_fence_released = false
            generation_two_runtime_marker =
                generation_two_fence_reservation.runtime_output_root_reservation_marker
            generation_two_queued_marker =
                generation_two_fence_reservation.output_root_reservation_marker
            try
                Test.@test Base.timedwait(
                    () -> isfile(generation_two_fence_ready) ||
                          !Base.process_running(
                              generation_two_fence_process,
                          ),
                    30.0;
                    pollint = 0.01,
                ) == :ok
                Test.@test Base.process_running(
                    generation_two_fence_process,
                )
                Test.@test isdir(generation_two_runtime_marker)
                competing_parent = joinpath(
                    generation_two_fence_outdir,
                    "missing-competing-parent",
                )
                competing_outdir = joinpath(
                    competing_parent,
                    "output",
                )
                competing_pid_lock =
                    Mycelia._output_root_reservation_lock_path_from_canonical(
                        competing_outdir,
                    )
                competing_outputs = Mycelia._metamdbg_output_paths(
                    competing_outdir,
                    31,
                )
                competing_reservation =
                    Mycelia._metamdbg_submission_reservation(
                        competing_outputs,
                        input_contract,
                        31;
                        owner_token =
                            "generation-two-competing-owner",
                    )
                competing_action_entered = Ref(false)
                competing_error = try
                    Mycelia._with_allowed_output_root_ancestor_locks(
                        (generation_two_queued_marker,),
                    ) do
                        Mycelia._with_metamdbg_output_domain_lock(
                            competing_outdir,
                        ) do
                            competing_action_entered[] = true
                            Mycelia._create_metamdbg_submission_reservation!(
                                competing_reservation,
                                competing_outdir,
                            )
                        end
                    end
                    nothing
                catch caught
                    caught isa InterruptException && rethrow()
                    caught
                end
                Test.@test competing_error isa ArgumentError
                if competing_error isa ArgumentError
                    Test.@test occursin(
                        generation_two_runtime_marker,
                        sprint(showerror, competing_error),
                    )
                end
                Test.@test !competing_action_entered[]
                Test.@test !ispath(competing_parent)
                Test.@test !ispath(competing_pid_lock)
                Test.@test !ispath(competing_reservation.path)
                Test.@test !ispath(
                    competing_reservation.output_root_reservation_marker,
                )
                write(generation_two_fence_continue, "continue\n")
                generation_two_fence_released = true
            finally
                if !generation_two_fence_released
                    write(generation_two_fence_continue, "continue\n")
                end
                wait(generation_two_fence_process)
            end
            Test.@test success(generation_two_fence_process)
            Test.@test ispath(generation_two_fence_assembly_marker)
            Test.@test !ispath(generation_two_fence_reservation.path)
            Test.@test !ispath(generation_two_queued_marker)
            Test.@test !ispath(generation_two_runtime_marker)
            Test.@test isdir(
                generation_two_fence_reservation.consumed_path,
            )

            transient_find_parent = joinpath(
                temporary_root,
                "executor-transient-ancestor-find",
            )
            transient_find_outdir = joinpath(transient_find_parent, "output")
            transient_find_outputs, transient_find_reservation =
                create_bound_runtime_fixture(
                    transient_find_outdir,
                    "runtime-transient-ancestor-find-owner",
                    "887",
                )
            transient_find_assembly_marker = joinpath(
                temporary_root,
                "runtime-transient-ancestor-find-assembly-ran",
            )
            transient_conda_runner = joinpath(
                temporary_root,
                "transient-enumerator-conda",
            )
            transient_find_counter = joinpath(
                temporary_root,
                "transient-ancestor-find.count",
            )
            write(
                transient_conda_runner,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "original_args=(\"\$@\")\n" *
                "while [ \"\$#\" -gt 0 ] && [ \"\$1\" != python ]; do shift; done\n" *
                "if [ \"\$#\" -ge 5 ] && [ \"\$4\" = " *
                Base.shell_escape(transient_find_parent) * " ]; then\n" *
                "  count=0\n" *
                "  [ ! -f " *
                Base.shell_escape(transient_find_counter) *
                " ] || count=\$(cat " *
                Base.shell_escape(transient_find_counter) * ")\n" *
                "  count=\$((count + 1))\n" *
                "  printf '%s\\n' \"\$count\" > " *
                Base.shell_escape(transient_find_counter) * "\n" *
                "  if [ \"\$count\" -eq 1 ]; then\n" *
                "    printf '%s\\0' \"\$4/.mycelia-output-root-reservation.partial\" > \"\$5\"\n" *
                "    exit 97\n" *
                "  fi\n" *
                "fi\n" *
                "exec " * Base.shell_escape(fake_conda) *
                " \"\${original_args[@]}\"\n",
            )
            chmod(transient_conda_runner, 0o700)
            transient_find_script = Mycelia._metamdbg_executor_script(
                "touch " *
                Base.shell_escape(transient_find_assembly_marker) *
                "; " * fake_asm,
                fake_gfa,
                transient_find_outputs,
                21,
                input_contract;
                conda_runner = transient_conda_runner,
                submission_reservation = transient_find_reservation,
            )
            transient_find_script_path = joinpath(
                temporary_root,
                "executor-transient-ancestor-find.sh",
            )
            write(transient_find_script_path, transient_find_script)
            transient_find_process = Base.withenv(
                "SLURM_JOB_ID" => transient_find_reservation.job_id,
            ) do
                run(`bash $(transient_find_script_path)`)
            end
            Test.@test success(transient_find_process)
            Test.@test parse(Int, strip(read(
                transient_find_counter,
                String,
            ))) == 3
            Test.@test isempty(
                enumeration_retry_remnants(transient_find_parent),
            )
            Test.@test ispath(transient_find_assembly_marker)
            Test.@test !ispath(transient_find_reservation.path)
            Test.@test isdir(transient_find_reservation.consumed_path)
            Test.@test !ispath(
                transient_find_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                transient_find_reservation.runtime_output_root_reservation_marker,
            )
            transient_find_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    transient_find_outdir,
                ),
            )
            transient_find_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    transient_find_metadata;
                    owner_token = transient_find_metadata.owner_token,
                    job_id = transient_find_metadata.job_id,
                    confirm_terminal = :completed,
                )
            Test.@test transient_find_reclaimed.recovery_reason == :completed

            persistent_find_parent = joinpath(
                temporary_root,
                "executor-persistent-ancestor-find",
            )
            persistent_find_outdir = joinpath(
                persistent_find_parent,
                "output",
            )
            persistent_find_outputs, persistent_find_reservation =
                create_bound_runtime_fixture(
                    persistent_find_outdir,
                    "runtime-persistent-ancestor-find-owner",
                    "888",
                )
            persistent_find_assembly_marker = joinpath(
                temporary_root,
                "runtime-persistent-ancestor-find-assembly-ran",
            )
            persistent_conda_runner = joinpath(
                temporary_root,
                "persistent-enumerator-conda",
            )
            persistent_find_counter = joinpath(
                temporary_root,
                "persistent-ancestor-find.count",
            )
            write(
                persistent_conda_runner,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "original_args=(\"\$@\")\n" *
                "while [ \"\$#\" -gt 0 ] && [ \"\$1\" != python ]; do shift; done\n" *
                "if [ \"\$#\" -ge 5 ] && [ \"\$4\" = " *
                Base.shell_escape(persistent_find_parent) * " ]; then\n" *
                "  count=0\n" *
                "  [ ! -f " *
                Base.shell_escape(persistent_find_counter) *
                " ] || count=\$(cat " *
                Base.shell_escape(persistent_find_counter) * ")\n" *
                "  count=\$((count + 1))\n" *
                "  printf '%s\\n' \"\$count\" > " *
                Base.shell_escape(persistent_find_counter) * "\n" *
                "  printf '%s\\0' \"\$4/.mycelia-output-root-reservation.partial\" > \"\$5\"\n" *
                "  exit 97\n" *
                "fi\n" *
                "exec " * Base.shell_escape(fake_conda) *
                " \"\${original_args[@]}\"\n",
            )
            chmod(persistent_conda_runner, 0o700)
            persistent_find_script = Mycelia._metamdbg_executor_script(
                "touch " *
                Base.shell_escape(persistent_find_assembly_marker) *
                "; " * fake_asm,
                fake_gfa,
                persistent_find_outputs,
                21,
                input_contract;
                conda_runner = persistent_conda_runner,
                submission_reservation = persistent_find_reservation,
            )
            persistent_find_script_path = joinpath(
                temporary_root,
                "executor-persistent-ancestor-find.sh",
            )
            write(persistent_find_script_path, persistent_find_script)
            persistent_find_log = IOBuffer()
            persistent_find_process = Base.withenv(
                "SLURM_JOB_ID" => persistent_find_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(`bash $(persistent_find_script_path)`);
                    stdout = persistent_find_log,
                    stderr = persistent_find_log,
                ))
            end
            persistent_find_output = String(take!(persistent_find_log))
            Test.@test !success(persistent_find_process)
            Test.@test occursin(
                "could not completely enumerate same-root or ancestor " *
                "output reservations after 3 attempts",
                persistent_find_output,
            )
            Test.@test parse(Int, strip(read(
                persistent_find_counter,
                String,
            ))) == 3
            Test.@test !isempty(
                enumeration_retry_remnants(persistent_find_parent),
            )
            Test.@test !ispath(persistent_find_assembly_marker)
            Test.@test isdir(persistent_find_reservation.path)
            Test.@test isfile(
                persistent_find_reservation.output_root_reservation_marker,
            )
            Test.@test isdir(
                persistent_find_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(persistent_find_outdir),
            )
            rm(
                persistent_find_reservation.runtime_output_root_reservation_marker,
            )
            remove_queued_runtime_fixture!(
                persistent_find_outdir,
                persistent_find_reservation,
            )

            unterminated_find_parent = joinpath(
                temporary_root,
                "executor-unterminated-ancestor-find",
            )
            unterminated_find_outdir = joinpath(
                unterminated_find_parent,
                "output",
            )
            unterminated_find_outputs, unterminated_find_reservation =
                create_bound_runtime_fixture(
                    unterminated_find_outdir,
                    "runtime-unterminated-ancestor-find-owner",
                    "890",
                )
            unterminated_find_assembly_marker = joinpath(
                temporary_root,
                "runtime-unterminated-ancestor-find-assembly-ran",
            )
            unterminated_conda_runner = joinpath(
                temporary_root,
                "unterminated-enumerator-conda",
            )
            unterminated_find_counter = joinpath(
                temporary_root,
                "unterminated-ancestor-find.count",
            )
            write(
                unterminated_conda_runner,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "original_args=(\"\$@\")\n" *
                "while [ \"\$#\" -gt 0 ] && [ \"\$1\" != python ]; do shift; done\n" *
                "if [ \"\$#\" -ge 5 ] && [ \"\$4\" = " *
                Base.shell_escape(unterminated_find_parent) * " ]; then\n" *
                "  count=0\n" *
                "  [ ! -f " *
                Base.shell_escape(unterminated_find_counter) *
                " ] || count=\$(cat " *
                Base.shell_escape(unterminated_find_counter) * ")\n" *
                "  count=\$((count + 1))\n" *
                "  printf '%s\\n' \"\$count\" > " *
                Base.shell_escape(unterminated_find_counter) * "\n" *
                "  printf '%s' \"\$4/unterminated-record\" > \"\$5\"\n" *
                "  exit 0\n" *
                "fi\n" *
                "exec " * Base.shell_escape(fake_conda) *
                " \"\${original_args[@]}\"\n",
            )
            chmod(unterminated_conda_runner, 0o700)
            unterminated_find_script = Mycelia._metamdbg_executor_script(
                "touch " *
                Base.shell_escape(unterminated_find_assembly_marker) *
                "; " * fake_asm,
                fake_gfa,
                unterminated_find_outputs,
                21,
                input_contract;
                conda_runner = unterminated_conda_runner,
                submission_reservation = unterminated_find_reservation,
            )
            unterminated_find_script_path = joinpath(
                temporary_root,
                "executor-unterminated-ancestor-find.sh",
            )
            write(unterminated_find_script_path, unterminated_find_script)
            unterminated_find_log = IOBuffer()
            unterminated_find_process = Base.withenv(
                "SLURM_JOB_ID" => unterminated_find_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(`bash $(unterminated_find_script_path)`);
                    stdout = unterminated_find_log,
                    stderr = unterminated_find_log,
                ))
            end
            unterminated_find_output = String(take!(unterminated_find_log))
            Test.@test !success(unterminated_find_process)
            Test.@test occursin(
                "successful inventory contains an unterminated NUL record",
                unterminated_find_output,
            )
            Test.@test occursin(
                "could not completely enumerate same-root or ancestor " *
                "output reservations after 3 attempts",
                unterminated_find_output,
            )
            Test.@test parse(Int, strip(read(
                unterminated_find_counter,
                String,
            ))) == 3
            Test.@test !isempty(
                enumeration_retry_remnants(unterminated_find_parent),
            )
            Test.@test !ispath(unterminated_find_assembly_marker)
            Test.@test isdir(unterminated_find_reservation.path)
            Test.@test isfile(
                unterminated_find_reservation.output_root_reservation_marker,
            )
            Test.@test isdir(
                unterminated_find_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(unterminated_find_outdir),
            )
            rm(
                unterminated_find_reservation.runtime_output_root_reservation_marker,
            )
            remove_queued_runtime_fixture!(
                unterminated_find_outdir,
                unterminated_find_reservation,
            )

            descendant_find_outdir =
                joinpath(temporary_root, "executor-descendant-find", "output")
            descendant_find_outputs, descendant_find_reservation =
                create_bound_runtime_fixture(
                    descendant_find_outdir,
                    "runtime-descendant-find-owner",
                    "875",
                )
            mkpath(descendant_find_outdir)
            descendant_find_assembly_marker = joinpath(
                temporary_root,
                "runtime-descendant-find-assembly-ran",
            )
            descendant_conda_runner = joinpath(
                temporary_root,
                "descendant-enumerator-conda",
            )
            descendant_find_counter = joinpath(
                temporary_root,
                "descendant-enumerator.count",
            )
            write(
                descendant_conda_runner,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "original_args=(\"\$@\")\n" *
                "while [ \"\$#\" -gt 0 ] && [ \"\$1\" != python ]; do shift; done\n" *
                "if [ \"\$#\" -ge 5 ] && [ \"\$4\" = " *
                Base.shell_escape(descendant_find_outdir) * " ]; then\n" *
                "  count=0\n" *
                "  [ ! -f " *
                Base.shell_escape(descendant_find_counter) *
                " ] || count=\$(cat " *
                Base.shell_escape(descendant_find_counter) * ")\n" *
                "  count=\$((count + 1))\n" *
                "  printf '%s\\n' \"\$count\" > " *
                Base.shell_escape(descendant_find_counter) * "\n" *
                "  printf '%s\\0' \"\$4/$(Mycelia._OUTPUT_ROOT_DURABLE_RESERVATION_PREFIX)partial\" > \"\$5\"\n" *
                "  echo 'synthetic descendant enumerator failure' >&2\n" *
                "  exit 97\n" *
                "fi\n" *
                "exec " * Base.shell_escape(fake_conda) *
                " \"\${original_args[@]}\"\n",
            )
            chmod(descendant_conda_runner, 0o700)
            descendant_find_script = Mycelia._metamdbg_executor_script(
                "touch $(Base.shell_escape(descendant_find_assembly_marker)); " *
                fake_asm,
                fake_gfa,
                descendant_find_outputs,
                21,
                input_contract;
                conda_runner = descendant_conda_runner,
                submission_reservation = descendant_find_reservation,
            )
            descendant_find_script_path = joinpath(
                temporary_root,
                "executor-descendant-find.sh",
            )
            write(descendant_find_script_path, descendant_find_script)
            descendant_find_log = IOBuffer()
            descendant_find_process = Base.withenv(
                "SLURM_JOB_ID" => descendant_find_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(`bash $(descendant_find_script_path)`);
                    stdout = descendant_find_log,
                    stderr = descendant_find_log,
                ))
            end
            descendant_find_output = String(take!(descendant_find_log))
            Test.@test !success(descendant_find_process)
            Test.@test occursin(
                "could not completely enumerate descendant output " *
                "reservations",
                descendant_find_output,
            )
            Test.@test occursin(
                "synthetic descendant enumerator failure",
                descendant_find_output,
            )
            Test.@test parse(Int, strip(read(
                descendant_find_counter,
                String,
            ))) == 3
            Test.@test !isempty(enumeration_retry_remnants(
                dirname(descendant_find_outdir),
            ))
            Test.@test !ispath(descendant_find_assembly_marker)
            Test.@test isdir(descendant_find_reservation.path)
            Test.@test isfile(
                descendant_find_reservation.output_root_reservation_marker,
            )
            Test.@test isdir(
                descendant_find_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(descendant_find_outdir),
            )
            rm(
                descendant_find_reservation.runtime_output_root_reservation_marker,
            )
            remove_queued_runtime_fixture!(
                descendant_find_outdir,
                descendant_find_reservation,
            )

            descendant_aba_outdir = joinpath(
                temporary_root,
                "executor-descendant-aba",
                "output",
            )
            descendant_aba_outputs, descendant_aba_reservation =
                create_bound_runtime_fixture(
                    descendant_aba_outdir,
                    "runtime-descendant-aba-owner",
                    "889",
                )
            mkpath(descendant_aba_outdir)
            descendant_aba_foreign_root = joinpath(
                descendant_aba_outdir,
                "nested",
                "foreign-output",
            )
            descendant_aba_foreign_marker =
                Mycelia._output_root_durable_reservation_path_from_canonical(
                    descendant_aba_foreign_root,
                    "runtime-descendant-aba-foreign-owner",
                )
            mkpath(dirname(descendant_aba_foreign_marker))
            write(descendant_aba_foreign_marker, "foreign descendant\n")
            chmod(descendant_aba_foreign_marker, 0o600)
            descendant_aba_replacement = joinpath(
                dirname(descendant_aba_outdir),
                "replacement-output",
            )
            descendant_aba_parked = joinpath(
                dirname(descendant_aba_outdir),
                "parked-output",
            )
            mkpath(descendant_aba_replacement)
            descendant_aba_open_ready = joinpath(
                temporary_root,
                "descendant-aba-open.ready",
            )
            descendant_aba_open_continue = joinpath(
                temporary_root,
                "descendant-aba-open.continue",
            )
            descendant_aba_enumeration_ready = joinpath(
                temporary_root,
                "descendant-aba-enumeration.ready",
            )
            descendant_aba_enumeration_continue = joinpath(
                temporary_root,
                "descendant-aba-enumeration.continue",
            )
            descendant_aba_open_hook = """
            if mode == "reservations-recursive" and os.fsdecode(root) == $(JSON.json(descendant_aba_outdir)):
                with open($(JSON.json(descendant_aba_open_ready)), "w", encoding="utf-8") as stream:
                    stream.write("ready\\n")
                while not os.path.exists($(JSON.json(descendant_aba_open_continue))):
                    time.sleep(0.01)
            """
            descendant_aba_enumeration_hook = """
            if mode == "reservations-recursive" and os.fsdecode(root) == $(JSON.json(descendant_aba_outdir)):
                with open($(JSON.json(descendant_aba_enumeration_ready)), "w", encoding="utf-8") as stream:
                    stream.write("ready\\n")
                while not os.path.exists($(JSON.json(descendant_aba_enumeration_continue))):
                    time.sleep(0.01)
            """
            descendant_aba_assembly_marker = joinpath(
                temporary_root,
                "runtime-descendant-aba-assembly-ran",
            )
            descendant_aba_script = Mycelia._metamdbg_executor_script(
                "touch " *
                Base.shell_escape(descendant_aba_assembly_marker) *
                "; " * fake_asm,
                fake_gfa,
                descendant_aba_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = descendant_aba_reservation,
                directory_enumerator_post_open_hook =
                    descendant_aba_open_hook,
                directory_enumerator_post_enumeration_hook =
                    descendant_aba_enumeration_hook,
            )
            descendant_aba_script_path = joinpath(
                temporary_root,
                "executor-descendant-aba.sh",
            )
            write(descendant_aba_script_path, descendant_aba_script)
            descendant_aba_log_path = joinpath(
                temporary_root,
                "executor-descendant-aba.log",
            )
            descendant_aba_log_io = open(descendant_aba_log_path, "w")
            descendant_aba_process = run(
                pipeline(
                    ignorestatus(addenv(
                        `bash $(descendant_aba_script_path)`,
                        "SLURM_JOB_ID" => descendant_aba_reservation.job_id,
                    ));
                    stdout = descendant_aba_log_io,
                    stderr = descendant_aba_log_io,
                );
                wait = false,
            )
            Test.@test Base.timedwait(
                () -> isfile(descendant_aba_open_ready) ||
                      !Base.process_running(descendant_aba_process),
                30.0;
                pollint = 0.01,
            ) == :ok
            Test.@test Base.process_running(descendant_aba_process)
            mv(descendant_aba_outdir, descendant_aba_parked)
            mv(descendant_aba_replacement, descendant_aba_outdir)
            write(descendant_aba_open_continue, "continue\n")
            Test.@test Base.timedwait(
                () -> isfile(descendant_aba_enumeration_ready) ||
                      !Base.process_running(descendant_aba_process),
                30.0;
                pollint = 0.01,
            ) == :ok
            Test.@test Base.process_running(descendant_aba_process)
            mv(descendant_aba_outdir, descendant_aba_replacement)
            mv(descendant_aba_parked, descendant_aba_outdir)
            write(descendant_aba_enumeration_continue, "continue\n")
            wait(descendant_aba_process)
            close(descendant_aba_log_io)
            descendant_aba_output = read(descendant_aba_log_path, String)
            Test.@test !success(descendant_aba_process)
            Test.@test occursin(
                "active descendant shared output-root reservation",
                descendant_aba_output,
            )
            Test.@test !ispath(descendant_aba_assembly_marker)
            Test.@test isdir(descendant_aba_reservation.path)
            Test.@test isfile(
                descendant_aba_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                descendant_aba_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                Mycelia._metamdbg_output_lock_path(descendant_aba_outdir),
            )
            Test.@test isfile(descendant_aba_foreign_marker)
            Test.@test isempty(enumeration_retry_remnants(
                dirname(descendant_aba_outdir),
            ))
            rm(descendant_aba_foreign_marker)
            remove_queued_runtime_fixture!(
                descendant_aba_outdir,
                descendant_aba_reservation,
            )

            long_runtime_outdir = joinpath(
                temporary_root,
                "executor-long-runtime",
                repeat("r", 240),
            )
            long_runtime_outputs, long_runtime_reservation =
                create_bound_runtime_fixture(
                    long_runtime_outdir,
                    "runtime-long-basename-owner",
                    "876",
                )
            long_runtime_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                long_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = long_runtime_reservation,
            )
            long_runtime_script_path =
                joinpath(temporary_root, "executor-long-runtime.sh")
            write(long_runtime_script_path, long_runtime_script)
            Test.@test Base.withenv(
                "SLURM_JOB_ID" => long_runtime_reservation.job_id,
            ) do
                success(`bash $(long_runtime_script_path)`)
            end
            Test.@test isfile(long_runtime_outputs.completion_marker)
            Test.@test !ispath(long_runtime_reservation.path)
            Test.@test !ispath(
                long_runtime_reservation.runtime_output_root_reservation_marker,
            )

            near_runtime_outdir =
                joinpath(temporary_root, "executor-near-name-max-input")
            near_runtime_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(near_name_max_reads, nothing),
                3,
            )
            near_runtime_outputs, near_runtime_reservation =
                create_bound_runtime_fixture(
                    near_runtime_outdir,
                    "runtime-near-name-max-owner",
                    "910";
                    fixture_input_contract = near_runtime_contract,
                )
            near_runtime_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                near_runtime_outputs,
                21,
                near_runtime_contract;
                conda_runner = fake_conda,
                submission_reservation = near_runtime_reservation,
            )
            Test.@test occursin(
                Mycelia._metamdbg_staged_input_name(
                    1,
                    near_name_max_reads,
                ),
                near_runtime_script,
            )
            near_runtime_script_path =
                joinpath(temporary_root, "executor-near-name-max-input.sh")
            write(near_runtime_script_path, near_runtime_script)
            Test.@test Base.withenv(
                "SLURM_JOB_ID" => near_runtime_reservation.job_id,
            ) do
                success(`bash $(near_runtime_script_path)`)
            end
            Test.@test isfile(near_runtime_outputs.completion_marker)
            Test.@test !ispath(near_runtime_reservation.path)
            Test.@test !ispath(
                near_runtime_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                near_runtime_reservation.runtime_output_root_reservation_marker,
            )

            backslash_runtime_outdir =
                joinpath(temporary_root, "executor-backslash-input")
            backslash_runtime_contract = Mycelia._metamdbg_input_contract(
                Mycelia._metamdbg_selected_input(backslash_reads, nothing),
                3,
            )
            backslash_runtime_outputs, backslash_runtime_reservation =
                create_bound_runtime_fixture(
                    backslash_runtime_outdir,
                    "runtime-backslash-owner",
                    "911";
                    fixture_input_contract = backslash_runtime_contract,
                )
            backslash_runtime_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                backslash_runtime_outputs,
                21,
                backslash_runtime_contract;
                conda_runner = fake_conda,
                submission_reservation = backslash_runtime_reservation,
            )
            backslash_runtime_script_path =
                joinpath(temporary_root, "executor-backslash-input.sh")
            write(backslash_runtime_script_path, backslash_runtime_script)
            fake_sha256_directory =
                joinpath(temporary_root, "backslash-sha256-bin")
            mkpath(fake_sha256_directory)
            fake_sha256_path =
                joinpath(fake_sha256_directory, "sha256sum")
            real_sha256_path = Sys.which("sha256sum")
            real_sha256_path === nothing &&
                error("sha256sum is required for this executor fixture")
            write(
                fake_sha256_path,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "if [ \"\$#\" -gt 0 ] && " *
                "printf '%s' \"\${!#}\" | grep -Fq '\\'; then\n" *
                "  printf '\\\\'\n" *
                "fi\n" *
                "exec " * Base.shell_escape(real_sha256_path) * " \"\$@\"\n",
            )
            chmod(fake_sha256_path, 0o700)
            Test.@test success(`bash -n $(fake_sha256_path)`)
            Test.@test Base.withenv(
                "SLURM_JOB_ID" => backslash_runtime_reservation.job_id,
                "PATH" =>
                    "$(fake_sha256_directory):$(get(ENV, "PATH", ""))",
            ) do
                success(`bash $(backslash_runtime_script_path)`)
            end
            Test.@test isfile(backslash_runtime_outputs.completion_marker)
            Test.@test !ispath(backslash_runtime_reservation.path)
            Test.@test !ispath(
                backslash_runtime_reservation.output_root_reservation_marker,
            )
            Test.@test !ispath(
                backslash_runtime_reservation.runtime_output_root_reservation_marker,
            )

            move_fail_parent =
                joinpath(temporary_root, "executor-queued-move-failure")
            move_fail_outdir = joinpath(move_fail_parent, "target")
            move_fail_outputs, move_fail_reservation =
                create_bound_runtime_fixture(
                    move_fail_outdir,
                    "runtime-queued-move-failure-owner",
                    "913",
                )
            move_fail_assembly_marker = joinpath(
                temporary_root,
                "runtime-queued-move-failure-assembly-ran",
            )
            move_fail_script = Mycelia._metamdbg_executor_script(
                "touch $(Base.shell_escape(move_fail_assembly_marker)); " *
                fake_asm,
                fake_gfa,
                move_fail_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = move_fail_reservation,
            )
            move_fail_script_path = joinpath(
                temporary_root,
                "executor-queued-move-failure.sh",
            )
            write(move_fail_script_path, move_fail_script)
            fake_mv_directory =
                joinpath(temporary_root, "queued-move-failure-bin")
            mkpath(fake_mv_directory)
            fake_mv_path = joinpath(fake_mv_directory, "mv")
            real_mv_path = something(Sys.which("mv"), "/bin/mv")
            move_fail_quarantine =
                move_fail_reservation.output_root_reservation_marker *
                ".removing." * move_fail_reservation.job_id
            write_side_effect_failure_wrapper!(
                fake_mv_path,
                real_mv_path,
                move_fail_quarantine,
                97,
            )
            move_fail_log = IOBuffer()
            move_fail_process = Base.withenv(
                "SLURM_JOB_ID" => move_fail_reservation.job_id,
                "PATH" =>
                    "$(fake_mv_directory):$(get(ENV, "PATH", ""))",
            ) do
                run(pipeline(
                    ignorestatus(Cmd(["bash", move_fail_script_path]));
                    stdout = move_fail_log,
                    stderr = move_fail_log,
                ))
            end
            move_fail_output = String(take!(move_fail_log))
            Test.@test !success(move_fail_process)
            Test.@test occursin(
                "queued shared-reservation quarantine command failed",
                move_fail_output,
            )
            Test.@test occursin(
                "failed to consume its queued shared output-root reservation",
                move_fail_output,
            )
            Test.@test !ispath(move_fail_assembly_marker)
            Test.@test !ispath(
                move_fail_reservation.output_root_reservation_marker,
            )
            Test.@test isfile(move_fail_quarantine)
            Test.@test read(move_fail_quarantine, String) ==
                       move_fail_reservation.output_root_reservation_contents
            Test.@test !ispath(move_fail_reservation.path)
            Test.@test isdir(move_fail_reservation.runtime_path)
            Test.@test isdir(
                move_fail_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test isdir(
                Mycelia._metamdbg_output_lock_path(move_fail_outdir),
            )

            delete_fail_parent =
                joinpath(temporary_root, "executor-queued-delete-failure")
            delete_fail_outdir = joinpath(delete_fail_parent, "target")
            delete_fail_outputs, delete_fail_reservation =
                create_bound_runtime_fixture(
                    delete_fail_outdir,
                    "runtime-queued-delete-failure-owner",
                    "912",
                )
            delete_fail_assembly_marker = joinpath(
                temporary_root,
                "runtime-queued-delete-failure-assembly-ran",
            )
            delete_fail_script = Mycelia._metamdbg_executor_script(
                "touch $(Base.shell_escape(delete_fail_assembly_marker)); " *
                fake_asm,
                fake_gfa,
                delete_fail_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = delete_fail_reservation,
            )
            delete_fail_script_path = joinpath(
                temporary_root,
                "executor-queued-delete-failure.sh",
            )
            write(delete_fail_script_path, delete_fail_script)
            fake_rm_directory =
                joinpath(temporary_root, "queued-delete-failure-bin")
            mkpath(fake_rm_directory)
            fake_rm_path = joinpath(fake_rm_directory, "rm")
            real_rm_path = something(Sys.which("rm"), "/bin/rm")
            delete_fail_quarantine =
                delete_fail_reservation.output_root_reservation_marker *
                ".removing." * delete_fail_reservation.job_id
            write(
                fake_rm_path,
                "#!/usr/bin/env bash\n" *
                "set -euo pipefail\n" *
                "if [ \"\$#\" -gt 0 ] && [ \"\${!#}\" = " *
                Base.shell_escape(delete_fail_quarantine) * " ]; then\n" *
                "  " * Base.shell_escape(real_rm_path) * " -- " *
                Base.shell_escape(delete_fail_quarantine) * "\n" *
                "  exit 97\n" *
                "fi\n" *
                "exec " * Base.shell_escape(real_rm_path) * " \"\$@\"\n",
            )
            chmod(fake_rm_path, 0o700)
            Test.@test success(`bash -n $(fake_rm_path)`)
            delete_fail_log = IOBuffer()
            delete_fail_process = Base.withenv(
                "SLURM_JOB_ID" => delete_fail_reservation.job_id,
                "PATH" => "$(fake_rm_directory):$(get(ENV, "PATH", ""))",
            ) do
                run(pipeline(
                    ignorestatus(`bash $(delete_fail_script_path)`);
                    stdout = delete_fail_log,
                    stderr = delete_fail_log,
                ))
            end
            delete_fail_output = String(take!(delete_fail_log))
            Test.@test !success(delete_fail_process)
            Test.@test occursin(
                "failed to consume its queued shared output-root reservation",
                delete_fail_output,
            )
            Test.@test !ispath(delete_fail_assembly_marker)
            Test.@test !ispath(delete_fail_reservation.path)
            Test.@test isdir(delete_fail_reservation.runtime_path)
            Test.@test !ispath(
                delete_fail_reservation.output_root_reservation_marker,
            )
            Test.@test isdir(
                delete_fail_reservation.runtime_output_root_reservation_marker,
            )
            delete_fail_lock =
                Mycelia._metamdbg_output_lock_path(delete_fail_outdir)
            Test.@test isdir(delete_fail_lock)
            _test_metamdbg_error(
                () -> Mycelia._with_unicycler_output_lock(
                    delete_fail_outdir,
                ) do _
                    :unexpected
                end,
                ArgumentError,
                r"active same-root output-root reservation",
            )
            _test_metamdbg_error(
                () -> Mycelia._with_autocycler_output_lock(
                    joinpath(delete_fail_outdir, "child"),
                ) do _
                    :unexpected
                end,
                ArgumentError,
                r"active ancestor output-root reservation",
            )
            _test_metamdbg_error(
                () -> Mycelia._with_unicycler_output_lock(
                    delete_fail_parent,
                ) do _
                    :unexpected
                end,
                ArgumentError,
                r"active descendant output-root reservation",
            )
            rm(delete_fail_lock)
            rm(
                delete_fail_reservation.runtime_output_root_reservation_marker,
            )
            rm(delete_fail_reservation.runtime_path; recursive = true)

            runtime_cleanup_parent =
                joinpath(temporary_root, "executor-runtime-cleanup-domain")
            runtime_cleanup_outdir =
                joinpath(runtime_cleanup_parent, "target")
            runtime_cleanup_outputs, runtime_cleanup_reservation =
                create_bound_runtime_fixture(
                    runtime_cleanup_outdir,
                    "runtime-cleanup-blocker-owner",
                    "877",
                )
            runtime_cleanup_lock =
                Mycelia._metamdbg_output_lock_path(runtime_cleanup_outdir)
            runtime_cleanup_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                runtime_cleanup_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = runtime_cleanup_reservation,
                post_completion_publication_hook =
                    "touch \"\$lock_dir/cleanup-blocker\"",
            )
            runtime_cleanup_script_path =
                joinpath(temporary_root, "executor-runtime-cleanup.sh")
            write(runtime_cleanup_script_path, runtime_cleanup_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => runtime_cleanup_reservation.job_id,
            ) do
                success(`bash $(runtime_cleanup_script_path)`)
            end
            Test.@test isfile(runtime_cleanup_outputs.completion_marker)
            Test.@test !ispath(runtime_cleanup_reservation.path)
            Test.@test isdir(runtime_cleanup_reservation.runtime_path)
            Test.@test !ispath(
                runtime_cleanup_reservation.output_root_reservation_marker,
            )
            Test.@test isdir(runtime_cleanup_lock)
            Test.@test isfile(
                joinpath(runtime_cleanup_lock, "cleanup-blocker"),
            )
            Test.@test isdir(
                runtime_cleanup_reservation.runtime_output_root_reservation_marker,
            )
            runtime_cleanup_reservation_mode = stat(
                runtime_cleanup_reservation.runtime_output_root_reservation_marker,
            ).mode & 0o777
            Test.@test runtime_cleanup_reservation_mode == 0o700
            _test_metamdbg_error(
                () -> Mycelia._with_unicycler_output_lock(
                    runtime_cleanup_outdir,
                ) do _
                    :unexpected
                end,
                ArgumentError,
                r"active same-root output-root reservation",
            )
            _test_metamdbg_error(
                () -> Mycelia._with_autocycler_output_lock(
                    joinpath(runtime_cleanup_outdir, "child"),
                ) do _
                    :unexpected
                end,
                ArgumentError,
                r"active ancestor output-root reservation",
            )
            _test_metamdbg_error(
                () -> Mycelia._with_unicycler_output_lock(
                    runtime_cleanup_parent,
                ) do _
                    :unexpected
                end,
                ArgumentError,
                r"active descendant output-root reservation",
            )
            rm(runtime_cleanup_lock; recursive = true)
            rm(
                runtime_cleanup_reservation.runtime_output_root_reservation_marker,
            )
            rm(runtime_cleanup_reservation.runtime_path; recursive = true)

            exceptional_runtime_outdir = joinpath(
                temporary_root,
                "executor-runtime-exceptional-cleanup",
            )
            exceptional_runtime_outputs, exceptional_runtime_reservation =
                create_bound_runtime_fixture(
                    exceptional_runtime_outdir,
                    "runtime-exceptional-cleanup-owner",
                    "878",
                )
            exceptional_runtime_lock =
                Mycelia._metamdbg_output_lock_path(exceptional_runtime_outdir)
            exceptional_runtime_script = Mycelia._metamdbg_executor_script(
                "touch \"\$lock_dir/cleanup-blocker\"; exit 44",
                fake_gfa,
                exceptional_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = exceptional_runtime_reservation,
            )
            exceptional_runtime_script_path = joinpath(
                temporary_root,
                "executor-runtime-exceptional-cleanup.sh",
            )
            write(
                exceptional_runtime_script_path,
                exceptional_runtime_script,
            )
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => exceptional_runtime_reservation.job_id,
            ) do
                success(`bash $(exceptional_runtime_script_path)`)
            end
            Test.@test !ispath(exceptional_runtime_reservation.path)
            Test.@test isdir(exceptional_runtime_reservation.runtime_path)
            Test.@test !ispath(
                exceptional_runtime_reservation.output_root_reservation_marker,
            )
            Test.@test isdir(exceptional_runtime_lock)
            Test.@test isdir(
                exceptional_runtime_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test !ispath(
                exceptional_runtime_outputs.completion_marker,
            )
            rm(exceptional_runtime_lock; recursive = true)
            rm(
                exceptional_runtime_reservation.runtime_output_root_reservation_marker,
            )
            rm(exceptional_runtime_reservation.runtime_path; recursive = true)

            replaced_runtime_outdir = joinpath(
                temporary_root,
                "executor-runtime-replaced-cleanup-lock",
            )
            replaced_runtime_outputs, replaced_runtime_reservation =
                create_bound_runtime_fixture(
                    replaced_runtime_outdir,
                    "runtime-replaced-cleanup-owner",
                    "879",
                )
            replaced_runtime_lock =
                Mycelia._metamdbg_output_lock_path(replaced_runtime_outdir)
            replaced_runtime_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                replaced_runtime_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = replaced_runtime_reservation,
                post_completion_publication_hook =
                    "mv -- \"\$lock_dir\" \"\${lock_dir}.original\"; " *
                    "mkdir -m 700 -- \"\$lock_dir\"",
            )
            replaced_runtime_script_path = joinpath(
                temporary_root,
                "executor-runtime-replaced-cleanup-lock.sh",
            )
            write(replaced_runtime_script_path, replaced_runtime_script)
            Test.@test !Base.withenv(
                "SLURM_JOB_ID" => replaced_runtime_reservation.job_id,
            ) do
                success(`bash $(replaced_runtime_script_path)`)
            end
            Test.@test !ispath(replaced_runtime_reservation.path)
            Test.@test isdir(replaced_runtime_reservation.runtime_path)
            Test.@test !ispath(
                replaced_runtime_reservation.output_root_reservation_marker,
            )
            Test.@test isdir(replaced_runtime_lock)
            Test.@test isempty(readdir(replaced_runtime_lock))
            Test.@test isdir(
                replaced_runtime_reservation.runtime_output_root_reservation_marker,
            )
            rm(replaced_runtime_lock)
            rm(replaced_runtime_lock * ".original")
            rm(
                replaced_runtime_reservation.runtime_output_root_reservation_marker,
            )
            rm(replaced_runtime_reservation.runtime_path; recursive = true)

            replaced_shared_runtime_outdir = joinpath(
                temporary_root,
                "executor-runtime-replaced-shared-reservation",
            )
            replaced_shared_outputs, replaced_shared_reservation =
                create_bound_runtime_fixture(
                    replaced_shared_runtime_outdir,
                    "runtime-replaced-shared-owner",
                    "880",
                )
            replaced_shared_private_lock =
                Mycelia._metamdbg_output_lock_path(
                    replaced_shared_runtime_outdir,
                )
            replaced_shared_marker =
                replaced_shared_reservation.runtime_output_root_reservation_marker
            replaced_shared_original = replaced_shared_marker * ".original"
            replaced_shared_cleanup_sentinel =
                Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                    replaced_shared_runtime_outdir,
                )
            replaced_shared_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                replaced_shared_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = replaced_shared_reservation,
                post_completion_publication_hook =
                    "mv -- \"\$runtime_output_root_reservation\" " *
                    "\"\${runtime_output_root_reservation}.original\"; " *
                    "mkdir -m 700 -- \"\$runtime_output_root_reservation\"",
            )
            replaced_shared_script_path = joinpath(
                temporary_root,
                "executor-runtime-replaced-shared-reservation.sh",
            )
            write(replaced_shared_script_path, replaced_shared_script)
            Test.@test isempty(
                Mycelia._metamdbg_pending_submission_job_paths(
                    replaced_shared_runtime_outdir,
                ),
            )
            replaced_shared_log = IOBuffer()
            replaced_shared_process = Base.withenv(
                "SLURM_JOB_ID" => replaced_shared_reservation.job_id,
            ) do
                return run(pipeline(
                    ignorestatus(`bash $(replaced_shared_script_path)`),
                    stdout = replaced_shared_log,
                    stderr = replaced_shared_log,
                ))
            end
            replaced_shared_output = String(take!(replaced_shared_log))
            Test.@test !success(replaced_shared_process)
            Test.@test occursin(
                "runtime output-root reservation identity changed after " *
                "acquisition",
                replaced_shared_output,
            )
            Test.@test !ispath(replaced_shared_outputs.completion_marker)
            Test.@test !ispath(replaced_shared_reservation.path)
            Test.@test isdir(replaced_shared_reservation.runtime_path)
            Test.@test isdir(replaced_shared_private_lock)
            Test.@test isdir(replaced_shared_marker)
            Test.@test isdir(replaced_shared_original)
            Test.@test isdir(replaced_shared_cleanup_sentinel)
            _test_metamdbg_error(
                () -> Mycelia._with_unicycler_output_lock(
                    replaced_shared_runtime_outdir,
                ) do _
                    :unexpected
                end,
                ArgumentError,
                r"active same-root output-root reservation",
            )
            ispath(replaced_shared_private_lock) &&
                rm(replaced_shared_private_lock; recursive = true)
            ispath(replaced_shared_marker) &&
                rm(replaced_shared_marker; recursive = true)
            ispath(replaced_shared_original) &&
                rm(replaced_shared_original; recursive = true)
            ispath(replaced_shared_cleanup_sentinel) &&
                rm(replaced_shared_cleanup_sentinel; recursive = true)
            ispath(replaced_shared_reservation.runtime_path) && rm(
                replaced_shared_reservation.runtime_path;
                recursive = true,
            )
            if ispath(replaced_shared_reservation.path)
                remove_queued_runtime_fixture!(
                    replaced_shared_runtime_outdir,
                    replaced_shared_reservation,
                )
            end

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
                matching = (;
                    job_id = "900",
                    bound_cluster = nothing,
                    runtime_cluster = nothing,
                    succeeds = true,
                ),
                missing = (;
                    job_id = nothing,
                    bound_cluster = nothing,
                    runtime_cluster = nothing,
                    succeeds = false,
                ),
                mismatched = (;
                    job_id = "999",
                    bound_cluster = nothing,
                    runtime_cluster = nothing,
                    succeeds = false,
                ),
                federated_matching = (;
                    job_id = "900",
                    bound_cluster = "cluster-a",
                    runtime_cluster = "cluster-a",
                    succeeds = true,
                ),
                federated_mismatched_cluster = (;
                    job_id = "900",
                    bound_cluster = "cluster-a",
                    runtime_cluster = "cluster-b",
                    succeeds = false,
                ),
                federated_missing_cluster = (;
                    job_id = "900",
                    bound_cluster = "cluster-a",
                    runtime_cluster = nothing,
                    succeeds = false,
                ),
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
                        runtime_case.bound_cluster,
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
                    "SLURM_CLUSTER_NAME" => runtime_case.runtime_cluster,
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
            Test.@test executable_completion.manifest.toolchain.package_inventory ==
                       _test_metamdbg_toolchain().package_inventory
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

            artifact_drift_outdir =
                joinpath(temporary_root, "executor-artifact-drift")
            artifact_drift_outputs, artifact_drift_reservation =
                create_bound_runtime_fixture(
                    artifact_drift_outdir,
                    "artifact-drift-owner",
                    "916",
                )
            artifact_drift_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                artifact_drift_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = artifact_drift_reservation,
                post_package_inventory_hook =
                    "printf 'H\\tVN:Z:1.0\\nS\\tcontig-1\\tTGCA\\n' > " *
                    "\"\$graph_source\"",
            )
            artifact_drift_script_path =
                joinpath(temporary_root, "artifact-drift.sh")
            write(artifact_drift_script_path, artifact_drift_script)
            artifact_drift_log = IOBuffer()
            artifact_drift_process = Base.withenv(
                "SLURM_JOB_ID" => artifact_drift_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(Cmd([
                        "bash",
                        artifact_drift_script_path,
                    ]));
                    stdout = artifact_drift_log,
                    stderr = artifact_drift_log,
                ))
            end
            Test.@test !success(artifact_drift_process)
            Test.@test occursin(
                "artifacts changed after package-inventory validation",
                String(take!(artifact_drift_log)),
            )
            Test.@test !ispath(artifact_drift_outputs.contigs_gz)
            Test.@test !ispath(artifact_drift_outputs.graph_alias)
            Test.@test !ispath(artifact_drift_outputs.contract_marker)
            Test.@test !ispath(artifact_drift_outputs.completion_marker)
            Test.@test isfile(artifact_drift_outputs.contigs_plain)
            Test.@test isfile(joinpath(
                artifact_drift_outdir,
                "assemblyGraph_k21_4bps.gfa",
            ))

            contract_drift_outdir =
                joinpath(temporary_root, "executor-contract-drift")
            contract_drift_outputs, contract_drift_reservation =
                create_bound_runtime_fixture(
                    contract_drift_outdir,
                    "contract-drift-owner",
                    "917",
                )
            contract_drift_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                contract_drift_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = contract_drift_reservation,
                post_completion_publication_hook =
                    "printf 'mutated-contract\\n' > \"\$contract_marker\"",
            )
            contract_drift_script_path =
                joinpath(temporary_root, "contract-drift.sh")
            write(contract_drift_script_path, contract_drift_script)
            contract_drift_log = IOBuffer()
            contract_drift_process = Base.withenv(
                "SLURM_JOB_ID" => contract_drift_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(Cmd([
                        "bash",
                        contract_drift_script_path,
                    ]));
                    stdout = contract_drift_log,
                    stderr = contract_drift_log,
                ))
            end
            Test.@test !success(contract_drift_process)
            Test.@test occursin(
                "provenance contract changed after exact binding",
                String(take!(contract_drift_log)),
            )
            Test.@test !ispath(contract_drift_outputs.completion_marker)
            Test.@test isfile(contract_drift_outputs.contract_marker)
            Test.@test read(
                contract_drift_outputs.contract_marker,
                String,
            ) == "mutated-contract\n"
            Test.@test !ispath(contract_drift_outputs.graph_alias)

            for (binding_kind, binding_job_id, expected_binding_error) in (
                    (
                        :owner,
                        "920",
                        "queued owner-record identity changed after binding",
                    ),
                    (
                        :shared_marker,
                        "921",
                        "queued shared-reservation identity changed after " *
                        "binding",
                    ),
                    (
                        :contract_child,
                        "922",
                        "queued owner children changed after binding",
                    ),
                    (
                        :job_child,
                        "923",
                        "queued owner children changed after binding",
                    ),
                )
                binding_outdir = joinpath(
                    temporary_root,
                    "executor-write-once-$(binding_kind)",
                )
                binding_outputs, binding_reservation =
                    create_bound_runtime_fixture(
                        binding_outdir,
                        "write-once-$(binding_kind)-owner",
                        binding_job_id,
                    )
                binding_path_variable = if binding_kind == :owner
                    "submission_reservation_dir"
                elseif binding_kind == :shared_marker
                    "submission_output_root_reservation"
                elseif binding_kind == :contract_child
                    "submission_reservation_contract"
                else
                    "submission_reservation_job"
                end
                held_binding = "\${secure_tmpdir}/held-$(binding_kind)"
                replace_binding = if binding_kind == :owner
                    "mv -- \"\$$(binding_path_variable)\" " *
                    "\"$(held_binding)\"; " *
                    "cp -R -- \"$(held_binding)\" " *
                    "\"\$$(binding_path_variable)\"; " *
                    "chmod 700 \"\$$(binding_path_variable)\"; " *
                    "chmod 600 \"\$$(binding_path_variable)\"/*.json"
                else
                    "mv -- \"\$$(binding_path_variable)\" " *
                    "\"$(held_binding)\"; " *
                    "cp -- \"$(held_binding)\" " *
                    "\"\$$(binding_path_variable)\"; " *
                    "chmod 600 \"\$$(binding_path_variable)\""
                end
                binding_script = Mycelia._metamdbg_executor_script(
                    fake_asm,
                    fake_gfa,
                    binding_outputs,
                    21,
                    input_contract;
                    conda_runner = fake_conda,
                    submission_reservation = binding_reservation,
                    post_runtime_marker_publication_hook = replace_binding,
                )
                binding_script_path = joinpath(
                    temporary_root,
                    "write-once-$(binding_kind).sh",
                )
                write(binding_script_path, binding_script)
                binding_log = IOBuffer()
                binding_process = Base.withenv(
                    "SLURM_JOB_ID" => binding_reservation.job_id,
                ) do
                    run(pipeline(
                        ignorestatus(Cmd([
                            "bash",
                            binding_script_path,
                        ]));
                        stdout = binding_log,
                        stderr = binding_log,
                    ))
                end
                Test.@test !success(binding_process)
                Test.@test occursin(
                    expected_binding_error,
                    String(take!(binding_log)),
                )
                for cleanup_path in (
                        binding_reservation.path,
                        binding_reservation.path * ".held",
                        binding_reservation.path *
                        ".held-contract_child",
                        binding_reservation.path * ".held-job_child",
                        binding_reservation.runtime_path,
                        binding_reservation.consumed_path,
                        binding_reservation.output_root_reservation_marker,
                        binding_reservation.output_root_reservation_marker *
                        ".held",
                        binding_reservation.runtime_output_root_reservation_marker,
                        Mycelia._metamdbg_output_lock_path(binding_outdir),
                    )
                    if ispath(cleanup_path) || islink(cleanup_path)
                        rm(cleanup_path; force = true, recursive = true)
                    end
                end
            end

            for (release_phase, release_job_id) in (
                    (:runtime_marker, "918"),
                    (:consumed_owner, "919"),
                )
                release_outdir = joinpath(
                    temporary_root,
                    "executor-$(release_phase)-late-mutation",
                )
                release_outputs, release_reservation =
                    create_bound_runtime_fixture(
                        release_outdir,
                        "$(release_phase)-late-mutation-owner",
                        release_job_id,
                    )
                mutation_hook = if release_phase == :runtime_marker
                    "rm -- \"\$contigs_gz\""
                else
                    "printf 'late-lifecycle-mutation\\n' > " *
                    "\"\$contract_marker\""
                end
                release_script = if release_phase == :runtime_marker
                    Mycelia._metamdbg_executor_script(
                        fake_asm,
                        fake_gfa,
                        release_outputs,
                        21,
                        input_contract;
                        conda_runner = fake_conda,
                        submission_reservation = release_reservation,
                        post_runtime_marker_release_hook = mutation_hook,
                    )
                else
                    Mycelia._metamdbg_executor_script(
                        fake_asm,
                        fake_gfa,
                        release_outputs,
                        21,
                        input_contract;
                        conda_runner = fake_conda,
                        submission_reservation = release_reservation,
                        post_consumed_owner_rename_hook = mutation_hook,
                    )
                end
                release_script_path = joinpath(
                    temporary_root,
                    "$(release_phase)-late-mutation.sh",
                )
                write(release_script_path, release_script)
                release_log = IOBuffer()
                release_process = Base.withenv(
                    "SLURM_JOB_ID" => release_reservation.job_id,
                ) do
                    run(pipeline(
                        ignorestatus(Cmd([
                            "bash",
                            release_script_path,
                        ]));
                        stdout = release_log,
                        stderr = release_log,
                    ))
                end
                Test.@test !success(release_process)
                release_output = String(take!(release_log))
                expected_release_error = if release_phase == :runtime_marker
                    "compressed contigs changed type during artifact snapshot"
                else
                    "provenance contract changed during final lifecycle " *
                    "cleanup"
                end
                Test.@test occursin(
                    expected_release_error,
                    release_output,
                )
                Test.@test isfile(release_outputs.completion_marker)
                if release_phase == :runtime_marker
                    Test.@test !ispath(release_outputs.contigs_gz)
                else
                    Test.@test read(
                        release_outputs.contract_marker,
                        String,
                    ) == "late-lifecycle-mutation\n"
                end
                Test.@test isdir(release_reservation.consumed_path)
            end

            resurrected_outdir = joinpath(
                temporary_root,
                "executor-released-marker-resurrection",
            )
            resurrected_outputs, resurrected_reservation =
                create_bound_runtime_fixture(
                    resurrected_outdir,
                    "released-marker-resurrection-owner",
                    "920",
                )
            resurrected_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                resurrected_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = resurrected_reservation,
                post_runtime_marker_release_hook =
                    "mkdir -m 700 -- " *
                    "\"\$runtime_output_root_reservation\"",
            )
            resurrected_script_path = joinpath(
                temporary_root,
                "released-marker-resurrection.sh",
            )
            write(resurrected_script_path, resurrected_script)
            resurrected_log = IOBuffer()
            resurrected_process = Base.withenv(
                "SLURM_JOB_ID" => resurrected_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(Cmd(["bash", resurrected_script_path]));
                    stdout = resurrected_log,
                    stderr = resurrected_log,
                ))
            end
            resurrected_output = String(take!(resurrected_log))
            Test.@test !success(resurrected_process)
            Test.@test occursin(
                "released lifecycle path reappeared before successful exit",
                resurrected_output,
            )
            Test.@test isdir(
                resurrected_reservation.runtime_output_root_reservation_marker,
            )
            Test.@test isdir(resurrected_reservation.consumed_path)
            Test.@test isfile(resurrected_outputs.completion_marker)
            rm(
                resurrected_reservation.runtime_output_root_reservation_marker,
            )
            rm(resurrected_reservation.consumed_path; recursive = true)

            replaced_child_outdir = joinpath(
                temporary_root,
                "executor-consumed-child-replacement",
            )
            replaced_child_outputs, replaced_child_reservation =
                create_bound_runtime_fixture(
                    replaced_child_outdir,
                    "consumed-child-replacement-owner",
                    "921",
                )
            replaced_child_inode = stat(
                replaced_child_reservation.contract_marker,
            ).inode
            replaced_child_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                replaced_child_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = replaced_child_reservation,
                post_consumed_owner_rename_hook =
                    "cp -- \"\$consumed_submission_reservation_contract\" " *
                    "\"\${consumed_submission_reservation_contract}.new\"; " *
                    "chmod 600 " *
                    "\"\${consumed_submission_reservation_contract}.new\"; " *
                    "rm -- \"\$consumed_submission_reservation_contract\"; " *
                    "mv -- " *
                    "\"\${consumed_submission_reservation_contract}.new\" " *
                    "\"\$consumed_submission_reservation_contract\"",
            )
            replaced_child_script_path = joinpath(
                temporary_root,
                "consumed-child-replacement.sh",
            )
            write(replaced_child_script_path, replaced_child_script)
            replaced_child_log = IOBuffer()
            replaced_child_process = Base.withenv(
                "SLURM_JOB_ID" => replaced_child_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(Cmd(["bash", replaced_child_script_path]));
                    stdout = replaced_child_log,
                    stderr = replaced_child_log,
                ))
            end
            replaced_child_output = String(take!(replaced_child_log))
            replaced_consumed_contract = joinpath(
                replaced_child_reservation.consumed_path,
                "contract.json",
            )
            Test.@test !success(replaced_child_process)
            Test.@test occursin(
                "consumed lifecycle contract changed before successful exit",
                replaced_child_output,
            )
            Test.@test isfile(replaced_child_outputs.completion_marker)
            Test.@test read(replaced_consumed_contract, String) ==
                       replaced_child_reservation.contents
            Test.@test stat(replaced_consumed_contract).inode !=
                       replaced_child_inode
            rm(replaced_child_reservation.consumed_path; recursive = true)

            third_child_outdir = joinpath(
                temporary_root,
                "executor-consumed-third-child",
            )
            third_child_outputs, third_child_reservation =
                create_bound_runtime_fixture(
                    third_child_outdir,
                    "consumed-third-child-owner",
                    "922",
                )
            third_child_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                third_child_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = third_child_reservation,
                post_consumed_owner_rename_hook =
                    "printf 'unexpected\\n' > " *
                    "\"\$consumed_submission_reservation_dir/third.json\"; " *
                    "chmod 600 " *
                    "\"\$consumed_submission_reservation_dir/third.json\"",
            )
            third_child_script_path = joinpath(
                temporary_root,
                "consumed-third-child.sh",
            )
            write(third_child_script_path, third_child_script)
            third_child_log = IOBuffer()
            third_child_process = Base.withenv(
                "SLURM_JOB_ID" => third_child_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(Cmd(["bash", third_child_script_path]));
                    stdout = third_child_log,
                    stderr = third_child_log,
                ))
            end
            third_child_output = String(take!(third_child_log))
            third_child_path = joinpath(
                third_child_reservation.consumed_path,
                "third.json",
            )
            Test.@test !success(third_child_process)
            Test.@test occursin(
                "consumed lifecycle owner has an unexpected entry set",
                third_child_output,
            )
            Test.@test isfile(third_child_path)
            Test.@test isfile(third_child_outputs.completion_marker)
            rm(third_child_reservation.consumed_path; recursive = true)

            completion_symlink_outdir = joinpath(
                temporary_root,
                "executor-final-completion-symlink",
            )
            completion_symlink_outputs, completion_symlink_reservation =
                create_bound_runtime_fixture(
                    completion_symlink_outdir,
                    "final-completion-symlink-owner",
                    "923",
                )
            completion_symlink_original =
                completion_symlink_outputs.completion_marker * ".original"
            completion_symlink_script = Mycelia._metamdbg_executor_script(
                fake_asm,
                fake_gfa,
                completion_symlink_outputs,
                21,
                input_contract;
                conda_runner = fake_conda,
                submission_reservation = completion_symlink_reservation,
                post_runtime_marker_release_hook =
                    "mv -- \"\$completion_marker\" " *
                    Base.shell_escape(completion_symlink_original) * "; " *
                    "ln -s -- " *
                    Base.shell_escape(basename(completion_symlink_original)) *
                    " \"\$completion_marker\"",
            )
            completion_symlink_script_path = joinpath(
                temporary_root,
                "final-completion-symlink.sh",
            )
            write(
                completion_symlink_script_path,
                completion_symlink_script,
            )
            completion_symlink_log = IOBuffer()
            completion_symlink_process = Base.withenv(
                "SLURM_JOB_ID" => completion_symlink_reservation.job_id,
            ) do
                run(pipeline(
                    ignorestatus(Cmd([
                        "bash",
                        completion_symlink_script_path,
                    ]));
                    stdout = completion_symlink_log,
                    stderr = completion_symlink_log,
                ))
            end
            completion_symlink_output =
                String(take!(completion_symlink_log))
            Test.@test !success(completion_symlink_process)
            Test.@test occursin(
                "completion manifest changed during final lifecycle cleanup",
                completion_symlink_output,
            )
            Test.@test islink(completion_symlink_outputs.completion_marker)
            Test.@test isfile(completion_symlink_original)
            Test.@test isdir(completion_symlink_reservation.consumed_path)
            rm(completion_symlink_reservation.consumed_path; recursive = true)

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
            Test.@test !ispath(publication_race_outputs.contract_marker)
            Test.@test !ispath(publication_race_outputs.graph_alias)
            Test.@test isfile(publication_race_outputs.contigs_gz)
            Test.@test !ispath(publication_race_reservation.path)

            for (failure_command, job_id, exit_code) in (
                    ("mv", "914", 97),
                    ("rm", "915", 98),
                )
                rollback_outdir = joinpath(
                    temporary_root,
                    "executor-completion-$(failure_command)-failure",
                )
                rollback_outputs, rollback_reservation =
                    create_bound_runtime_fixture(
                        rollback_outdir,
                        "completion-$(failure_command)-failure-owner",
                        job_id,
                    )
                rollback_quarantine =
                    rollback_outputs.completion_marker *
                    ".removing." *
                    Mycelia._output_root_reservation_identity(
                        rollback_outdir,
                    ) *
                    "." * job_id * ".completion-marker"
                rollback_bin = joinpath(
                    temporary_root,
                    "completion-$(failure_command)-failure-bin",
                )
                mkpath(rollback_bin)
                rollback_wrapper = joinpath(
                    rollback_bin,
                    failure_command,
                )
                real_command = something(
                    Sys.which(failure_command),
                    "/bin/$(failure_command)",
                )
                write_side_effect_failure_wrapper!(
                    rollback_wrapper,
                    real_command,
                    rollback_quarantine,
                    exit_code,
                )
                rollback_script = Mycelia._metamdbg_executor_script(
                    fake_asm,
                    fake_gfa,
                    rollback_outputs,
                    21,
                    input_contract;
                    conda_runner = fake_conda,
                    submission_reservation = rollback_reservation,
                    post_completion_publication_hook =
                        "printf '>contig-1\\nTGCA\\n' | gzip -c > " *
                        "\"\$contigs_gz\"",
                )
                rollback_script_path = joinpath(
                    temporary_root,
                    "completion-$(failure_command)-failure.sh",
                )
                write(rollback_script_path, rollback_script)
                rollback_log = IOBuffer()
                rollback_process = Base.withenv(
                    "SLURM_JOB_ID" => rollback_reservation.job_id,
                    "PATH" =>
                        "$(rollback_bin):$(get(ENV, "PATH", ""))",
                ) do
                    run(pipeline(
                        ignorestatus(Cmd(["bash", rollback_script_path]));
                        stdout = rollback_log,
                        stderr = rollback_log,
                    ))
                end
                rollback_output = String(take!(rollback_log))
                Test.@test !success(rollback_process)
                expected_failure = failure_command == "mv" ?
                                   "completion-marker rollback quarantine " *
                                   "command failed" :
                                   "exact completion-marker rollback " *
                                   "unlink failed"
                Test.@test occursin(expected_failure, rollback_output)
                Test.@test occursin(
                    "retained nonexact completion evidence",
                    rollback_output,
                )
                Test.@test !ispath(rollback_outputs.completion_marker)
                Test.@test ispath(rollback_quarantine) ==
                           (failure_command == "mv")
                Test.@test !ispath(rollback_outputs.contract_marker)
                Test.@test !ispath(rollback_outputs.graph_alias)
                Test.@test isfile(rollback_outputs.contigs_gz)
                Test.@test !ispath(rollback_reservation.path)
                Test.@test isdir(rollback_reservation.runtime_path)
                Test.@test isdir(
                    rollback_reservation.runtime_output_root_reservation_marker,
                )
                Test.@test isdir(
                    Mycelia._metamdbg_output_lock_path(rollback_outdir),
                )
            end

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
                    stdout = "123\n",
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
            Test.@test submitted_job_record["job_cluster"] === nothing
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
                    confirm_terminal = :running,
                ),
                ArgumentError,
                r"accepts only :failed or :completed",
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
                ArgumentError,
                r"active same-root output-root reservation",
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
                    stdout = "456\n",
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
            Test.@test immediate_complete.provenance.package_inventory ==
                       _test_metamdbg_toolchain().package_inventory
            consumed_metadata = only(
                Mycelia.inspect_metamdbg_submission_reservations(
                    immediate_outdir,
                ),
            )
            Test.@test consumed_metadata.reservation_state == :consumed
            Test.@test consumed_metadata.submission_state == :consumed
            Test.@test consumed_metadata.lifecycle_owner == :runtime
            Test.@test consumed_metadata.job_id == "456"
            Test.@test consumed_metadata.queued_reservation_identity === nothing
            Test.@test consumed_metadata.runtime_reservation_identity === nothing
            consumed_reclaimed =
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    consumed_metadata;
                    owner_token = consumed_metadata.owner_token,
                    job_id = consumed_metadata.job_id,
                    confirm_terminal = :completed,
                )
            Test.@test consumed_reclaimed.recovery_reason == :completed
            Test.@test isempty(
                Mycelia.inspect_metamdbg_submission_reservations(
                    immediate_outdir,
                ),
            )
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
                missing_stdout = (
                    result = Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = "12344",
                    ),
                    message = r"returned no exact stdout evidence",
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
                    message =
                        r"stdout did not contain the exact returned job reference",
                ),
                multiple_stdout_references = (
                    result = Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = "12346",
                        stdout = "12346\n12347\n",
                    ),
                    message =
                        r"stdout did not contain the exact returned job reference",
                ),
                noncanonical_cluster = (
                    result = Mycelia.SubmitResult(
                        ok = true,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = "12345",
                        job_cluster = " cluster-a ",
                        stdout = "12345;cluster-a\n",
                    ),
                    message = r"returned a non-exact job cluster",
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

            Test.@testset "accepted audit failures retain exact handles" begin
                accepted_audit_submission = function (
                        job_id::AbstractString,
                        job_cluster::AbstractString,
                        audit_failure::AbstractString,
                )
                    reference = "$(job_id);$(job_cluster)"
                    return Mycelia.SubmitResult(
                        ok = false,
                        dry_run = false,
                        held = true,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        job_id = String(job_id),
                        job_cluster = String(job_cluster),
                        stdout = "$(reference)\n",
                        errors = [String(audit_failure)],
                    )
                end
                release_calls = Ref(0)
                forbidden_release = function (
                        _job_id::AbstractString;
                        job_cluster::Union{Nothing, AbstractString} = nothing,
                )
                    release_calls[] += 1
                    error(
                        "accepted audit failure must not release " *
                        "$(repr(job_cluster))",
                    )
                end

                bound_outdir = joinpath(
                    temporary_root,
                    "accepted-audit-bound-output",
                )
                bound_failure = try
                    Mycelia._run_metamdbg(;
                        hifi_reads = valid_reads,
                        outdir = bound_outdir,
                        executor = Mycelia.SlurmExecutor(dry_run = false),
                        site = :scg,
                        dependency_checker = _test_metamdbg_toolchain,
                        local_runner = forbidden_runner,
                        submission_runner = (
                            _job::Mycelia.JobSpec,
                            _executor::Mycelia.AbstractExecutor,
                        ) -> accepted_audit_submission(
                            "333",
                            "cluster-a",
                            "synthetic post-acceptance audit drift",
                        ),
                        submission_release_runner = forbidden_release,
                    )
                    nothing
                catch caught
                    caught
                end
                Test.@test bound_failure isa ErrorException
                bound_failure_message = sprint(showerror, bound_failure)
                Test.@test occursin(
                    "synthetic post-acceptance audit drift",
                    bound_failure_message,
                )
                Test.@test occursin(
                    "exact scheduler reference was durably bound",
                    bound_failure_message,
                )
                Test.@test occursin(
                    "job was not released",
                    bound_failure_message,
                )
                Test.@test release_calls[] == 0
                bound_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        bound_outdir,
                    ),
                )
                Test.@test bound_metadata.submission_state == :submitted
                Test.@test bound_metadata.job_id == "333"
                Test.@test bound_metadata.job_cluster == "cluster-a"
                bound_reclaimed =
                    Mycelia.reclaim_metamdbg_submission_reservation!(
                        bound_metadata;
                        owner_token = bound_metadata.owner_token,
                        job_id = bound_metadata.job_id,
                        job_cluster = bound_metadata.job_cluster,
                        confirm_cancelled = true,
                    )
                Test.@test bound_reclaimed.recovery_reason == :cancelled

                bind_failure_outdir = joinpath(
                    temporary_root,
                    "accepted-audit-bind-failure-output",
                )
                bind_failure = ErrorException(
                    "synthetic accepted-audit durable bind failure",
                )
                observed_bind_failure = try
                    Mycelia._run_metamdbg(;
                        hifi_reads = valid_reads,
                        outdir = bind_failure_outdir,
                        executor = Mycelia.SlurmExecutor(dry_run = false),
                        site = :scg,
                        dependency_checker = _test_metamdbg_toolchain,
                        local_runner = forbidden_runner,
                        submission_runner = (
                            _job::Mycelia.JobSpec,
                            _executor::Mycelia.AbstractExecutor,
                        ) -> accepted_audit_submission(
                            "334",
                            "cluster-b",
                            "synthetic bind-path audit drift",
                        ),
                        submission_job_binder = function (
                                reservation::NamedTuple,
                                job_id::AbstractString,
                                job_cluster::AbstractString,
                        )
                            return Mycelia._bind_metamdbg_submission_job!(
                                reservation,
                                job_id,
                                job_cluster;
                                post_pending_job_record_creation_hook =
                                    (_reservation, _path) ->
                                        throw(bind_failure),
                            )
                        end,
                        submission_release_runner = forbidden_release,
                    )
                    nothing
                catch caught
                    caught
                end
                Test.@test observed_bind_failure isa ErrorException
                bind_failure_message = sprint(
                    showerror,
                    observed_bind_failure,
                )
                Test.@test occursin(
                    "synthetic bind-path audit drift",
                    bind_failure_message,
                )
                Test.@test occursin(
                    "synthetic accepted-audit durable bind failure",
                    bind_failure_message,
                )
                Test.@test release_calls[] == 0
                pending_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        bind_failure_outdir,
                    ),
                )
                Test.@test pending_metadata.submission_state ==
                           :submission_pending
                Test.@test pending_metadata.job_id === nothing
                Test.@test pending_metadata.pending_job_id == "334"
                Test.@test pending_metadata.pending_job_cluster == "cluster-b"
                recovered_bound =
                    Mycelia.bind_metamdbg_submission_reservation_job!(
                        pending_metadata;
                        owner_token = pending_metadata.owner_token,
                        job_id = pending_metadata.pending_job_id,
                        job_cluster = pending_metadata.pending_job_cluster,
                        confirm_submitted = true,
                    )
                recovered_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        bind_failure_outdir,
                    ),
                )
                Test.@test recovered_metadata.job_id == recovered_bound.job_id
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    recovered_metadata;
                    owner_token = recovered_metadata.owner_token,
                    job_id = recovered_metadata.job_id,
                    job_cluster = recovered_metadata.job_cluster,
                    confirm_cancelled = true,
                )

                interrupt_outdir = joinpath(
                    temporary_root,
                    "accepted-audit-bind-interrupt-output",
                )
                bind_interrupt = InterruptException()
                observed_interrupt = Test.@test_logs (
                    :warn,
                    r"backend audit failed after scheduler acceptance.*synthetic interrupt-path audit drift.*SLURM job 335;cluster-c",
                ) min_level=Logging.Warn match_mode=:any begin
                    try
                        Mycelia._run_metamdbg(;
                            hifi_reads = valid_reads,
                            outdir = interrupt_outdir,
                            executor =
                                Mycelia.SlurmExecutor(dry_run = false),
                            site = :scg,
                            dependency_checker = _test_metamdbg_toolchain,
                            local_runner = forbidden_runner,
                            submission_runner = (
                                _job::Mycelia.JobSpec,
                                _executor::Mycelia.AbstractExecutor,
                            ) -> accepted_audit_submission(
                                "335",
                                "cluster-c",
                                "synthetic interrupt-path audit drift",
                            ),
                            submission_job_binder = function (
                                    reservation::NamedTuple,
                                    job_id::AbstractString,
                                    job_cluster::AbstractString,
                            )
                                return Mycelia._bind_metamdbg_submission_job!(
                                    reservation,
                                    job_id,
                                    job_cluster;
                                    post_pending_job_record_creation_hook =
                                        (_reservation, _path) ->
                                            throw(bind_interrupt),
                                )
                            end,
                            submission_release_runner = forbidden_release,
                        )
                        nothing
                    catch caught
                        caught
                    end
                end
                Test.@test observed_interrupt === bind_interrupt
                Test.@test release_calls[] == 0
                interrupt_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        interrupt_outdir,
                    ),
                )
                Test.@test interrupt_metadata.submission_state ==
                           :submission_pending
                Test.@test interrupt_metadata.pending_job_id == "335"
                Test.@test interrupt_metadata.pending_job_cluster ==
                           "cluster-c"
                interrupt_bound =
                    Mycelia.bind_metamdbg_submission_reservation_job!(
                        interrupt_metadata;
                        owner_token = interrupt_metadata.owner_token,
                        job_id = interrupt_metadata.pending_job_id,
                        job_cluster = interrupt_metadata.pending_job_cluster,
                        confirm_submitted = true,
                    )
                interrupt_recovered = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        interrupt_outdir,
                    ),
                )
                Test.@test interrupt_recovered.job_id == interrupt_bound.job_id
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    interrupt_recovered;
                    owner_token = interrupt_recovered.owner_token,
                    job_id = interrupt_recovered.job_id,
                    job_cluster = interrupt_recovered.job_cluster,
                    confirm_cancelled = true,
                )
                Test.@test release_calls[] == 0
            end

            Test.@testset "accepted callback binds before postcheck" begin
                callback_result = function (
                        record::NamedTuple,
                        ok::Bool,
                        errors::Vector{String},
                )
                    return Mycelia.SubmitResult(
                        ok = ok,
                        dry_run = false,
                        held = record.held,
                        scheduler_acceptance = :accepted,
                        site = :scg,
                        backend = :sbatch,
                        artifact_path = record.artifact_path,
                        job_id = record.job_id,
                        job_cluster = record.job_cluster,
                        stdout = record.stdout,
                        errors = errors,
                    )
                end
                accepted_record = function (
                        job_id::AbstractString,
                        job_cluster::AbstractString,
                        artifact_path::AbstractString,
                )
                    return (;
                        held = true,
                        job_id = String(job_id),
                        job_cluster = String(job_cluster),
                        stdout = "$(job_id);$(job_cluster)\n",
                        artifact_path = String(artifact_path),
                    )
                end

                success_outdir = joinpath(
                    temporary_root,
                    "accepted-callback-success-output",
                )
                success_record = accepted_record(
                    "601",
                    "cluster-callback-a",
                    joinpath(temporary_root, "callback-success.sbatch"),
                )
                success_callback_calls = Ref(0)
                success_bind_calls = Ref(0)
                success_release_calls = Ref(0)
                success_runner = function (
                        _job::Mycelia.JobSpec,
                        executor::Mycelia.AbstractExecutor;
                        accepted_sbatch_callback::Function,
                )
                    Test.@test executor isa Mycelia.SlurmExecutor
                    success_callback_calls[] += 1
                    accepted_sbatch_callback(success_record)
                    return callback_result(success_record, true, String[])
                end
                success = Mycelia._run_metamdbg(;
                    hifi_reads = valid_reads,
                    outdir = success_outdir,
                    executor = Mycelia.SlurmExecutor(dry_run = false),
                    site = :scg,
                    dependency_checker = _test_metamdbg_toolchain,
                    local_runner = forbidden_runner,
                    submission_runner = success_runner,
                    submission_runner_accepts_callback = true,
                    submission_job_binder = function (
                            reservation::NamedTuple,
                            job_id::AbstractString,
                            job_cluster::AbstractString,
                    )
                        success_bind_calls[] += 1
                        return Mycelia._bind_metamdbg_submission_job_after_submit!(
                            reservation,
                            job_id,
                            job_cluster,
                        )
                    end,
                    submission_release_runner = function (
                            job_id::AbstractString;
                            job_cluster::Union{Nothing, AbstractString} =
                                nothing,
                    )
                        success_release_calls[] += 1
                        Test.@test job_cluster == "cluster-callback-a"
                        return String(job_id)
                    end,
                )
                Test.@test success.status == :submitted
                Test.@test success_callback_calls[] == 1
                Test.@test success_bind_calls[] == 1
                Test.@test success_release_calls[] == 1
                Test.@test success.submission_reservation.job_id == "601"
                Test.@test success.submission_reservation.job_cluster ==
                           "cluster-callback-a"
                success_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        success_outdir,
                    ),
                )
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    success_metadata;
                    owner_token = success_metadata.owner_token,
                    job_id = success_metadata.job_id,
                    job_cluster = success_metadata.job_cluster,
                    confirm_cancelled = true,
                )

                mismatch_outdir = joinpath(
                    temporary_root,
                    "accepted-callback-mismatch-output",
                )
                mismatch_record = accepted_record(
                    "605",
                    "cluster-callback-e",
                    joinpath(temporary_root, "callback-mismatch.sbatch"),
                )
                mismatch_result_record = merge(mismatch_record, (;
                    job_id = "606",
                    stdout = "606;cluster-callback-e\n",
                ))
                mismatch_bind_calls = Ref(0)
                mismatch_release_calls = Ref(0)
                mismatch_runner = function (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor;
                        accepted_sbatch_callback::Function,
                )
                    accepted_sbatch_callback(mismatch_record)
                    return callback_result(
                        mismatch_result_record,
                        true,
                        String[],
                    )
                end
                mismatch_failure = try
                    Mycelia._run_metamdbg(;
                        hifi_reads = valid_reads,
                        outdir = mismatch_outdir,
                        executor = Mycelia.SlurmExecutor(dry_run = false),
                        site = :scg,
                        dependency_checker = _test_metamdbg_toolchain,
                        local_runner = forbidden_runner,
                        submission_runner = mismatch_runner,
                        submission_runner_accepts_callback = true,
                        submission_job_binder = function (
                                reservation::NamedTuple,
                                job_id::AbstractString,
                                job_cluster::AbstractString,
                        )
                            mismatch_bind_calls[] += 1
                            return Mycelia._bind_metamdbg_submission_job_after_submit!(
                                reservation,
                                job_id,
                                job_cluster,
                            )
                        end,
                        submission_release_runner = function (
                                _job_id::AbstractString;
                                job_cluster::Union{
                                    Nothing,
                                    AbstractString,
                                } = nothing,
                        )
                            mismatch_release_calls[] += 1
                            error("release forbidden $(repr(job_cluster))")
                        end,
                    )
                    nothing
                catch caught
                    caught
                end
                Test.@test mismatch_failure isa ErrorException
                Test.@test occursin(
                    "disagree about the exact job id",
                    sprint(showerror, mismatch_failure),
                )
                Test.@test mismatch_bind_calls[] == 1
                Test.@test mismatch_release_calls[] == 0
                mismatch_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        mismatch_outdir,
                    ),
                )
                Test.@test mismatch_metadata.job_id == "605"
                Test.@test mismatch_metadata.job_cluster ==
                           "cluster-callback-e"
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    mismatch_metadata;
                    owner_token = mismatch_metadata.owner_token,
                    job_id = mismatch_metadata.job_id,
                    job_cluster = mismatch_metadata.job_cluster,
                    confirm_cancelled = true,
                )

                failure_outdir = joinpath(
                    temporary_root,
                    "accepted-callback-bind-failure-output",
                )
                failure_record = accepted_record(
                    "602",
                    "cluster-callback-b",
                    joinpath(temporary_root, "callback-failure.sbatch"),
                )
                callback_bind_failure = ErrorException(
                    "synthetic synchronous callback bind failure",
                )
                failure_bind_calls = Ref(0)
                failure_release_calls = Ref(0)
                failure_runner = function (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor;
                        accepted_sbatch_callback::Function,
                )
                    callback_error = try
                        accepted_sbatch_callback(failure_record)
                        nothing
                    catch caught
                        caught isa InterruptException && rethrow()
                        caught
                    end
                    Test.@test callback_error === callback_bind_failure
                    return callback_result(
                        failure_record,
                        false,
                        [
                            "Accepted sbatch callback failed after scheduler " *
                            "acceptance: " *
                            sprint(showerror, callback_error),
                        ],
                    )
                end
                observed_callback_failure = try
                    Mycelia._run_metamdbg(;
                        hifi_reads = valid_reads,
                        outdir = failure_outdir,
                        executor =
                            Mycelia.SlurmExecutor(dry_run = false),
                        site = :scg,
                        dependency_checker = _test_metamdbg_toolchain,
                        local_runner = forbidden_runner,
                        submission_runner = failure_runner,
                        submission_runner_accepts_callback = true,
                        submission_job_binder = function (
                                reservation::NamedTuple,
                                job_id::AbstractString,
                                job_cluster::AbstractString,
                        )
                            failure_bind_calls[] += 1
                            return Mycelia._bind_metamdbg_submission_job!(
                                reservation,
                                job_id,
                                job_cluster;
                                post_pending_job_record_creation_hook =
                                    (_reservation, _path) ->
                                        throw(callback_bind_failure),
                            )
                        end,
                        submission_release_runner = function (
                                _job_id::AbstractString;
                                job_cluster::Union{
                                    Nothing,
                                    AbstractString,
                                } = nothing,
                        )
                            failure_release_calls[] += 1
                            error("release forbidden $(repr(job_cluster))")
                        end,
                    )
                    nothing
                catch caught
                    caught
                end
                Test.@test observed_callback_failure isa ErrorException
                callback_failure_message = sprint(
                    showerror,
                    observed_callback_failure,
                )
                Test.@test occursin(
                    "synthetic synchronous callback bind failure",
                    callback_failure_message,
                )
                Test.@test failure_bind_calls[] == 1
                Test.@test failure_release_calls[] == 0
                callback_pending = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        failure_outdir,
                    ),
                )
                Test.@test callback_pending.pending_job_id == "602"
                Test.@test callback_pending.pending_job_cluster ==
                           "cluster-callback-b"
                callback_bound =
                    Mycelia.bind_metamdbg_submission_reservation_job!(
                        callback_pending;
                        owner_token = callback_pending.owner_token,
                        job_id = callback_pending.pending_job_id,
                        job_cluster = callback_pending.pending_job_cluster,
                        confirm_submitted = true,
                    )
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    only(Mycelia.inspect_metamdbg_submission_reservations(
                        failure_outdir,
                    ));
                    owner_token = callback_bound.owner_token,
                    job_id = callback_bound.job_id,
                    job_cluster = callback_bound.job_cluster,
                    confirm_cancelled = true,
                )

                interrupt_outdir = joinpath(
                    temporary_root,
                    "accepted-callback-bind-interrupt-output",
                )
                interrupt_record = accepted_record(
                    "603",
                    "cluster-callback-c",
                    joinpath(temporary_root, "callback-interrupt.sbatch"),
                )
                callback_bind_interrupt = InterruptException()
                interrupt_bind_calls = Ref(0)
                interrupt_release_calls = Ref(0)
                interrupt_result_record = accepted_record(
                    "606",
                    "cluster-result",
                    joinpath(temporary_root, "callback-interrupt-result.sbatch"),
                )
                interrupt_runner = function (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor;
                        accepted_sbatch_callback::Function,
                )
                    retained_interrupt = try
                        accepted_sbatch_callback(interrupt_record)
                        nothing
                    catch caught
                        caught
                    end
                    Test.@test retained_interrupt === callback_bind_interrupt
                    return callback_result(
                        interrupt_result_record,
                        false,
                        ["callback interrupt retained after acceptance"],
                    )
                end
                observed_callback_interrupt = Test.@test_logs (
                    :warn,
                    r"SLURM job 603;cluster-callback-c.*pending scheduler-job evidence were preserved",
                ) min_level=Logging.Warn match_mode=:any begin
                    try
                        Mycelia._run_metamdbg(;
                            hifi_reads = valid_reads,
                            outdir = interrupt_outdir,
                            executor =
                                Mycelia.SlurmExecutor(dry_run = false),
                            site = :scg,
                            dependency_checker = _test_metamdbg_toolchain,
                            local_runner = forbidden_runner,
                            submission_runner = interrupt_runner,
                            submission_runner_accepts_callback = true,
                            submission_job_binder = function (
                                    reservation::NamedTuple,
                                    job_id::AbstractString,
                                    job_cluster::AbstractString,
                            )
                                interrupt_bind_calls[] += 1
                                return Mycelia._bind_metamdbg_submission_job!(
                                    reservation,
                                    job_id,
                                    job_cluster;
                                    post_pending_job_record_creation_hook =
                                        (_reservation, _path) ->
                                            throw(callback_bind_interrupt),
                                )
                            end,
                            submission_release_runner = function (
                                    _job_id::AbstractString;
                                    job_cluster::Union{
                                        Nothing,
                                        AbstractString,
                                    } = nothing,
                            )
                                interrupt_release_calls[] += 1
                                error("release forbidden $(repr(job_cluster))")
                            end,
                        )
                        nothing
                    catch caught
                        caught
                    end
                end
                Test.@test observed_callback_interrupt ===
                           callback_bind_interrupt
                Test.@test interrupt_bind_calls[] == 1
                Test.@test interrupt_release_calls[] == 0
                callback_interrupted = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        interrupt_outdir,
                    ),
                )
                Test.@test callback_interrupted.pending_job_id == "603"
                Test.@test callback_interrupted.pending_job_cluster ==
                           "cluster-callback-c"
                interrupt_bound =
                    Mycelia.bind_metamdbg_submission_reservation_job!(
                        callback_interrupted;
                        owner_token = callback_interrupted.owner_token,
                        job_id = callback_interrupted.pending_job_id,
                        job_cluster = callback_interrupted.pending_job_cluster,
                        confirm_submitted = true,
                    )
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    only(Mycelia.inspect_metamdbg_submission_reservations(
                        interrupt_outdir,
                    ));
                    owner_token = interrupt_bound.owner_token,
                    job_id = interrupt_bound.job_id,
                    job_cluster = interrupt_bound.job_cluster,
                    confirm_cancelled = true,
                )

                duplicate_outdir = joinpath(
                    temporary_root,
                    "accepted-callback-duplicate-after-interrupt-output",
                )
                duplicate_record = accepted_record(
                    "607",
                    "cluster-callback-duplicate",
                    joinpath(temporary_root, "callback-duplicate.sbatch"),
                )
                duplicate_interrupt = InterruptException()
                duplicate_bind_calls = Ref(0)
                duplicate_release_calls = Ref(0)
                duplicate_runner = function (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor;
                        accepted_sbatch_callback::Function,
                )
                    first_error = try
                        accepted_sbatch_callback(duplicate_record)
                        nothing
                    catch caught
                        caught
                    end
                    Test.@test first_error === duplicate_interrupt
                    second_error = try
                        accepted_sbatch_callback(duplicate_record)
                        nothing
                    catch caught
                        caught
                    end
                    Test.@test second_error === duplicate_interrupt
                    return callback_result(
                        duplicate_record,
                        false,
                        ["duplicate callback retained original interrupt"],
                    )
                end
                observed_duplicate_interrupt = Test.@test_logs (
                    :warn,
                    r"SLURM job 607;cluster-callback-duplicate.*pending scheduler-job evidence were preserved",
                ) min_level=Logging.Warn match_mode=:any begin
                    try
                        Mycelia._run_metamdbg(;
                            hifi_reads = valid_reads,
                            outdir = duplicate_outdir,
                            executor =
                                Mycelia.SlurmExecutor(dry_run = false),
                            site = :scg,
                            dependency_checker = _test_metamdbg_toolchain,
                            local_runner = forbidden_runner,
                            submission_runner = duplicate_runner,
                            submission_runner_accepts_callback = true,
                            submission_job_binder = function (
                                    reservation::NamedTuple,
                                    job_id::AbstractString,
                                    job_cluster::AbstractString,
                            )
                                duplicate_bind_calls[] += 1
                                return Mycelia._bind_metamdbg_submission_job!(
                                    reservation,
                                    job_id,
                                    job_cluster;
                                    post_pending_job_record_creation_hook =
                                        (_reservation, _path) ->
                                            throw(duplicate_interrupt),
                                )
                            end,
                            submission_release_runner = function (
                                    _job_id::AbstractString;
                                    job_cluster::Union{
                                        Nothing,
                                        AbstractString,
                                    } = nothing,
                            )
                                duplicate_release_calls[] += 1
                                error("release forbidden $(repr(job_cluster))")
                            end,
                        )
                        nothing
                    catch caught
                        caught
                    end
                end
                Test.@test observed_duplicate_interrupt === duplicate_interrupt
                Test.@test duplicate_bind_calls[] == 1
                Test.@test duplicate_release_calls[] == 0
                duplicate_pending = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        duplicate_outdir,
                    ),
                )
                Test.@test duplicate_pending.pending_job_id == "607"
                Test.@test duplicate_pending.pending_job_cluster ==
                           "cluster-callback-duplicate"
                duplicate_bound =
                    Mycelia.bind_metamdbg_submission_reservation_job!(
                        duplicate_pending;
                        owner_token = duplicate_pending.owner_token,
                        job_id = duplicate_pending.pending_job_id,
                        job_cluster = duplicate_pending.pending_job_cluster,
                        confirm_submitted = true,
                    )
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    only(Mycelia.inspect_metamdbg_submission_reservations(
                        duplicate_outdir,
                    ));
                    owner_token = duplicate_bound.owner_token,
                    job_id = duplicate_bound.job_id,
                    job_cluster = duplicate_bound.job_cluster,
                    confirm_cancelled = true,
                )

                validation_outdir = joinpath(
                    temporary_root,
                    "accepted-callback-validation-interrupt-output",
                )
                validation_interrupt = InterruptException()
                validation_record = (;
                    held = true,
                    job_id = _TestInterruptingMetamdbgString(
                        validation_interrupt,
                    ),
                    job_cluster = nothing,
                    stdout = "608\n",
                    artifact_path = joinpath(
                        temporary_root,
                        "callback-validation-interrupt.sbatch",
                    ),
                )
                Test.@test Base.ncodeunits(validation_record.job_id) == 1
                Test.@test Base.codeunit(validation_record.job_id, 1) ==
                           Base.codeunit("x", 1)
                Test.@test Base.isvalid(validation_record.job_id, 1)
                validation_release_calls = Ref(0)
                validation_runner_failure = ErrorException(
                    "synthetic runner failure after validation interrupt",
                )
                validation_runner = function (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor;
                        accepted_sbatch_callback::Function,
                )
                    callback_error = try
                        accepted_sbatch_callback(validation_record)
                        nothing
                    catch caught
                        caught
                    end
                    Test.@test callback_error === validation_interrupt
                    throw(validation_runner_failure)
                end
                observed_validation_interrupt = Test.@test_logs (
                    :warn,
                    r"before it retained a canonical exact scheduler reference.*unbound reservation.*scheduler name.*squeue --name",
                ) min_level=Logging.Warn match_mode=:any begin
                    try
                        Mycelia._run_metamdbg(;
                            hifi_reads = valid_reads,
                            outdir = validation_outdir,
                            executor =
                                Mycelia.SlurmExecutor(dry_run = false),
                            site = :scg,
                            dependency_checker = _test_metamdbg_toolchain,
                            local_runner = forbidden_runner,
                            submission_runner = validation_runner,
                            submission_runner_accepts_callback = true,
                            submission_job_binder = (
                                _reservation::NamedTuple,
                                _job_id::AbstractString,
                            ) -> error("validation binder must not run"),
                            submission_release_runner = function (
                                    _job_id::AbstractString;
                                    job_cluster::Union{
                                        Nothing,
                                        AbstractString,
                                    } = nothing,
                            )
                                validation_release_calls[] += 1
                                error("release forbidden $(repr(job_cluster))")
                            end,
                        )
                        nothing
                    catch caught
                        caught
                    end
                end
                Test.@test observed_validation_interrupt ===
                           validation_interrupt
                Test.@test validation_release_calls[] == 0
                validation_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        validation_outdir,
                    ),
                )
                Test.@test validation_metadata.job_id === nothing
                Test.@test validation_metadata.pending_job_id === nothing
                Test.@test validation_metadata.submission_state == :reserved
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    validation_metadata;
                    owner_token = validation_metadata.owner_token,
                    confirm_not_submitted = true,
                )

                postcheck_outdir = joinpath(
                    temporary_root,
                    "accepted-callback-postcheck-interrupt-output",
                )
                postcheck_record = accepted_record(
                    "604",
                    "cluster-callback-d",
                    joinpath(temporary_root, "callback-postcheck.sbatch"),
                )
                postcheck_interrupt = InterruptException()
                postcheck_bind_calls = Ref(0)
                postcheck_release_calls = Ref(0)
                postcheck_runner = function (
                        _job::Mycelia.JobSpec,
                        _executor::Mycelia.AbstractExecutor;
                        accepted_sbatch_callback::Function,
                )
                    accepted_sbatch_callback(postcheck_record)
                    throw(postcheck_interrupt)
                end
                observed_postcheck_interrupt = Test.@test_logs (
                    :warn,
                    r"SLURM job 604;cluster-callback-d.*exact scheduler reference was durably bound",
                ) min_level=Logging.Warn match_mode=:any begin
                    try
                        Mycelia._run_metamdbg(;
                            hifi_reads = valid_reads,
                            outdir = postcheck_outdir,
                            executor =
                                Mycelia.SlurmExecutor(dry_run = false),
                            site = :scg,
                            dependency_checker = _test_metamdbg_toolchain,
                            local_runner = forbidden_runner,
                            submission_runner = postcheck_runner,
                            submission_runner_accepts_callback = true,
                            submission_job_binder = function (
                                    reservation::NamedTuple,
                                    job_id::AbstractString,
                                    job_cluster::AbstractString,
                            )
                                postcheck_bind_calls[] += 1
                                return Mycelia._bind_metamdbg_submission_job_after_submit!(
                                    reservation,
                                    job_id,
                                    job_cluster,
                                )
                            end,
                            submission_release_runner = function (
                                    _job_id::AbstractString;
                                    job_cluster::Union{
                                        Nothing,
                                        AbstractString,
                                    } = nothing,
                            )
                                postcheck_release_calls[] += 1
                                error("release forbidden $(repr(job_cluster))")
                            end,
                        )
                        nothing
                    catch caught
                        caught
                    end
                end
                Test.@test observed_postcheck_interrupt === postcheck_interrupt
                Test.@test postcheck_bind_calls[] == 1
                Test.@test postcheck_release_calls[] == 0
                postcheck_metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        postcheck_outdir,
                    ),
                )
                Test.@test postcheck_metadata.job_id == "604"
                Test.@test postcheck_metadata.job_cluster ==
                           "cluster-callback-d"
                Mycelia.reclaim_metamdbg_submission_reservation!(
                    postcheck_metadata;
                    owner_token = postcheck_metadata.owner_token,
                    job_id = postcheck_metadata.job_id,
                    job_cluster = postcheck_metadata.job_cluster,
                    confirm_cancelled = true,
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
                        stdout = "111\n",
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
                        stdout = "222\n",
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

        Test.@testset "runtime transition crash recovery is identity-bound" begin
            for phase in (
                    :post_runtime_marker,
                    :post_private_rename,
                    :post_queued_release,
                    :post_runtime_release,
                    :post_consumed_rename,
            )
                runtime_outdir = joinpath(
                    temporary_root,
                    "$(phase)-runtime-transition-output",
                )
                live_metadata = Ref{Any}(nothing)
                bound = _kill_metamdbg_runtime_transition!(
                    runtime_outdir,
                    valid_reads,
                    phase;
                    live_check = function (
                            process::Base.Process,
                            reservation::NamedTuple,
                    )
                        Test.@test Base.process_running(process)
                        owner_path = if phase == :post_runtime_marker
                            reservation.path
                        elseif phase == :post_consumed_rename
                            reservation.consumed_path
                        else
                            reservation.runtime_path
                        end
                        Test.@test isdir(owner_path)
                        runtime_marker_exists = isdir(
                            reservation.runtime_output_root_reservation_marker,
                        )
                        Test.@test runtime_marker_exists == !(phase in (
                            :post_runtime_release,
                            :post_consumed_rename,
                        ))
                        private_lock_path =
                            Mycelia._metamdbg_output_lock_path(runtime_outdir)
                        private_lock_exists = isdir(private_lock_path)
                        Test.@test private_lock_exists == !(phase in (
                            :post_runtime_marker,
                            :post_runtime_release,
                            :post_consumed_rename,
                        ))
                        queued_exists = isfile(
                            reservation.output_root_reservation_marker,
                        )
                        Test.@test queued_exists ==
                                   (phase in (
                            :post_runtime_marker,
                            :post_private_rename,
                        ))
                        before = (;
                            owner = Mycelia._metamdbg_output_lock_identity(
                                owner_path,
                            ),
                            runtime = runtime_marker_exists ?
                                      Mycelia._metamdbg_output_lock_identity(
                                reservation.runtime_output_root_reservation_marker,
                            ) : nothing,
                            private = private_lock_exists ?
                                      Mycelia._metamdbg_output_lock_identity(
                                private_lock_path,
                            ) : nothing,
                        )
                        if phase == :post_consumed_rename
                            _test_metamdbg_error(
                                () -> Mycelia.inspect_metamdbg_submission_reservations(
                                    runtime_outdir;
                                    confirm_process_dead = true,
                                ),
                                ErrorException,
                                r"requires a pre-existing canonical local PID",
                            )
                            Test.@test Mycelia._metamdbg_output_lock_identity(
                                owner_path,
                            ) == before.owner
                            confirmed = only(
                                Mycelia.inspect_metamdbg_submission_reservations(
                                    runtime_outdir,
                                ),
                            )
                            Test.@test confirmed.reservation_state == :consumed
                        else
                            _test_metamdbg_error(
                                () -> Mycelia.inspect_metamdbg_submission_reservations(
                                runtime_outdir;
                                confirm_process_dead = true,
                                ),
                                ErrorException,
                                r"cannot recover scheduler-owned runtime state",
                            )
                        end
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            owner_path,
                        ) == before.owner
                        if runtime_marker_exists
                            Test.@test Mycelia._metamdbg_output_lock_identity(
                                reservation.runtime_output_root_reservation_marker,
                            ) == before.runtime
                        end
                        if private_lock_exists
                            Test.@test Mycelia._metamdbg_output_lock_identity(
                                private_lock_path,
                            ) == before.private
                        end
                        Test.@test isfile(
                            reservation.output_root_reservation_marker,
                        ) == queued_exists
                        live_metadata[] = only(
                            Mycelia.inspect_metamdbg_submission_reservations(
                                runtime_outdir,
                            ),
                        )
                        Test.@test live_metadata[].lifecycle_owner == :runtime
                        Test.@test live_metadata[].private_lock_identity ==
                                   before.private
                        return nothing
                    end,
                )
                metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        runtime_outdir,
                    ),
                )
                expected_owner_path = if phase == :post_runtime_marker
                    bound.path
                elseif phase == :post_consumed_rename
                    bound.consumed_path
                else
                    bound.runtime_path
                end
                Test.@test metadata.path == expected_owner_path
                Test.@test metadata.job_id == bound.job_id
                Test.@test metadata.owner_token == bound.owner_token
                Test.@test metadata.reservation_state == if phase ==
                                                            :post_runtime_marker
                    :runtime_claiming
                elseif phase == :post_private_rename
                    :runtime_transition
                elseif phase == :post_runtime_release
                    :runtime_release_pending
                elseif phase == :post_consumed_rename
                    :consumed
                else
                    :runtime
                end
                Test.@test metadata.reservation_identity ==
                           live_metadata[].reservation_identity
                _test_metamdbg_error(
                    () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                        metadata;
                        owner_token = metadata.owner_token,
                        confirm_not_submitted = true,
                    ),
                    ArgumentError,
                    r"cannot be reclaimed as not submitted",
                )
                _test_metamdbg_error(
                    () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                        metadata;
                        owner_token = metadata.owner_token,
                        job_id = "9999",
                        confirm_terminal = :failed,
                    ),
                    ErrorException,
                    r"job id does not match",
                )
                Test.@test isdir(expected_owner_path)
                recovered =
                    Mycelia.reclaim_metamdbg_submission_reservation!(
                        metadata;
                        owner_token = metadata.owner_token,
                        job_id = metadata.job_id,
                        confirm_terminal = :failed,
                    )
                Test.@test recovered.status == :reclaimed
                Test.@test recovered.recovery_reason == :failed
                Test.@test !ispath(expected_owner_path)
                Test.@test !ispath(
                    bound.runtime_output_root_reservation_marker,
                )
                Test.@test !ispath(bound.output_root_reservation_marker)
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(runtime_outdir),
                )
                Test.@test isempty(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        runtime_outdir,
                    ),
                )
            end
        end

        Test.@testset "explicit reclaim crash recovery stays discoverable" begin
            for phase in (:post_private_rename, :post_shared_release)
                reclaim_outdir = joinpath(
                    temporary_root,
                    "$(phase)-reclaim-transition-output",
                )
                live_identity = Ref{Any}(nothing)
                live_pid_identity = Ref{Any}(nothing)
                live_private_identity = Ref{Any}(nothing)
                live_cleanup_identity = Ref{Any}(nothing)
                reservation = _kill_metamdbg_reclaim_transition!(
                    reclaim_outdir,
                    valid_reads,
                    phase;
                    live_check = function (
                            process::Base.Process,
                            queued_reservation::NamedTuple,
                    )
                        Test.@test Base.process_running(process)
                        Test.@test isdir(queued_reservation.reclaiming_path)
                        Test.@test isdir(
                            Mycelia._metamdbg_output_lock_path(reclaim_outdir),
                        )
                        pid_lock_path =
                            Mycelia._output_root_reservation_lock_path_from_canonical(
                                reclaim_outdir,
                            )
                        cleanup_path =
                            Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                                reclaim_outdir,
                            )
                        Test.@test isfile(pid_lock_path)
                        Test.@test isdir(cleanup_path)
                        queued_marker_exists = isfile(
                            queued_reservation.output_root_reservation_marker,
                        )
                        Test.@test queued_marker_exists ==
                                   (phase == :post_private_rename)
                        owner_identity =
                            Mycelia._metamdbg_output_lock_identity(
                                queued_reservation.reclaiming_path,
                            )
                        private_identity =
                            Mycelia._metamdbg_output_lock_identity(
                                Mycelia._metamdbg_output_lock_path(
                                    reclaim_outdir,
                                ),
                            )
                        pid_identity =
                            _test_metamdbg_output_root_pid_lock_identity(
                                pid_lock_path,
                            )
                        cleanup_identity =
                            Mycelia._metamdbg_output_lock_identity(cleanup_path)
                        _test_metamdbg_error(
                            () -> Mycelia.inspect_metamdbg_submission_reservations(
                                reclaim_outdir;
                                confirm_process_dead = true,
                            ),
                            ErrorException,
                            r"still names a live or remotely unverifiable process",
                        )
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            queued_reservation.reclaiming_path,
                        ) == owner_identity
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            Mycelia._metamdbg_output_lock_path(reclaim_outdir),
                        ) == private_identity
                        Test.@test _test_metamdbg_output_root_pid_lock_identity(
                            pid_lock_path,
                        ) == pid_identity
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            cleanup_path,
                        ) == cleanup_identity
                        inspected = only(
                            Mycelia.inspect_metamdbg_submission_reservations(
                                reclaim_outdir,
                            ),
                        )
                        Test.@test inspected.lifecycle_owner == :recovery
                        Test.@test inspected.private_lock_identity ==
                                   private_identity
                        Test.@test inspected.cleanup_reservation_identity ==
                                   cleanup_identity
                        _test_metamdbg_error(
                            () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                                inspected;
                                owner_token = inspected.owner_token,
                                confirm_not_submitted = true,
                            ),
                            ArgumentError,
                            r"reinspect with confirm_process_dead=true",
                        )
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            queued_reservation.reclaiming_path,
                        ) == owner_identity
                        Test.@test _test_metamdbg_output_root_pid_lock_identity(
                            pid_lock_path,
                        ) == pid_identity
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            Mycelia._metamdbg_output_lock_path(reclaim_outdir),
                        ) == private_identity
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            cleanup_path,
                        ) == cleanup_identity
                        live_identity[] = inspected.reservation_identity
                        live_pid_identity[] = pid_identity
                        live_private_identity[] = private_identity
                        live_cleanup_identity[] = cleanup_identity
                        return nothing
                    end,
                )
                metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        reclaim_outdir;
                        confirm_process_dead = true,
                    ),
                )
                Test.@test live_pid_identity[] !== nothing
                Test.@test live_private_identity[] !== nothing
                Test.@test live_cleanup_identity[] !== nothing
                Test.@test metadata.path == reservation.reclaiming_path
                Test.@test metadata.reservation_identity == live_identity[]
                Test.@test metadata.private_lock_identity === nothing
                Test.@test metadata.cleanup_reservation_identity === nothing
                Test.@test metadata.reservation_state == if phase ==
                                                            :post_private_rename
                    :reclaiming
                else
                    :reclaim_release_pending
                end
                first_recovery_identity = metadata.reservation_identity
                _kill_metamdbg_reclaim_transition!(
                    reclaim_outdir,
                    valid_reads,
                    phase;
                    resume_reclaiming = true,
                    live_check = function (
                            process::Base.Process,
                            reclaiming_reservation::NamedTuple,
                    )
                        Test.@test Base.process_running(process)
                        owner_identity =
                            Mycelia._metamdbg_output_lock_identity(
                                reclaiming_reservation.reclaiming_path,
                            )
                        Test.@test owner_identity == first_recovery_identity
                        pid_lock_path =
                            Mycelia._output_root_reservation_lock_path_from_canonical(
                                reclaim_outdir,
                            )
                        private_lock_path =
                            Mycelia._metamdbg_output_lock_path(reclaim_outdir)
                        cleanup_path =
                            Mycelia._metamdbg_lifecycle_cleanup_reservation_path(
                                reclaim_outdir,
                            )
                        pid_identity =
                            _test_metamdbg_output_root_pid_lock_identity(
                                pid_lock_path,
                            )
                        private_identity =
                            Mycelia._metamdbg_output_lock_identity(
                                private_lock_path,
                            )
                        cleanup_identity =
                            Mycelia._metamdbg_output_lock_identity(cleanup_path)
                        live_metadata = only(
                            Mycelia.inspect_metamdbg_submission_reservations(
                                reclaim_outdir,
                            ),
                        )
                        _test_metamdbg_error(
                            () -> Mycelia.reclaim_metamdbg_submission_reservation!(
                                live_metadata;
                                owner_token = live_metadata.owner_token,
                                confirm_not_submitted = true,
                            ),
                            ArgumentError,
                            r"reinspect with confirm_process_dead=true",
                        )
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            reclaiming_reservation.reclaiming_path,
                        ) == owner_identity
                        Test.@test _test_metamdbg_output_root_pid_lock_identity(
                            pid_lock_path,
                        ) == pid_identity
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            private_lock_path,
                        ) == private_identity
                        Test.@test Mycelia._metamdbg_output_lock_identity(
                            cleanup_path,
                        ) == cleanup_identity
                        return nothing
                    end,
                )
                metadata = only(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        reclaim_outdir;
                        confirm_process_dead = true,
                    ),
                )
                Test.@test metadata.reservation_identity ==
                           first_recovery_identity
                Test.@test metadata.private_lock_identity === nothing
                Test.@test metadata.cleanup_reservation_identity === nothing
                recovered =
                    Mycelia.reclaim_metamdbg_submission_reservation!(
                        metadata;
                        owner_token = metadata.owner_token,
                        confirm_not_submitted = true,
                    )
                Test.@test recovered.recovery_reason == :not_submitted
                Test.@test !ispath(reservation.reclaiming_path)
                Test.@test !ispath(
                    reservation.output_root_reservation_marker,
                )
                Test.@test !ispath(
                    Mycelia._metamdbg_output_lock_path(reclaim_outdir),
                )
                Test.@test isempty(
                    Mycelia.inspect_metamdbg_submission_reservations(
                        reclaim_outdir,
                    ),
                )
            end
        end
    end
end
