import Test
import Mycelia

function _has_message(messages::Vector{String}, needle::String)
    return any(msg -> occursin(needle, msg), messages)
end

function _with_env(f::Function, overrides::AbstractDict)
    previous = Dict{String, Union{Nothing, String}}()
    for (key, value) in overrides
        key_string = String(key)
        previous[key_string] = get(ENV, key_string, nothing)
        if value === nothing
            pop!(ENV, key_string, nothing)
        else
            ENV[key_string] = String(value)
        end
    end

    try
        return f()
    finally
        for (key, value) in previous
            if value === nothing
                pop!(ENV, key, nothing)
            else
                ENV[key] = value
            end
        end
    end
end

function _normalize_fixture_logdir(text::AbstractString)
    legacy_logdir = "/Users/cameronprybol/workspace/slurmlogs"
    return replace(String(text), legacy_logdir => Mycelia.DEFAULT_SLURM_LOGDIR)
end

function _with_git_email(f::Function, email::AbstractString)
    mktempdir() do repo_dir
        run(`git -C $repo_dir init --quiet`)
        run(`git -C $repo_dir config user.name TestUser`)
        run(`git -C $repo_dir config user.email $email`)
        cd(repo_dir) do
            return f()
        end
    end
end

Test.@testset "JobSpec JSON roundtrip" begin
    job = Mycelia.JobSpec(
        job_name = "json-roundtrip",
        cmd = "echo hello",
        site = :scg,
        partition = "batch",
        account = "PI_EXAMPLE",
        time_limit = "12:00:00",
        nodes = 1,
        ntasks = 1,
        cpus_per_task = 4,
        mem_gb = 16,
        mail_type = "END,FAIL"
    )

    json_text = Mycelia.job_spec_to_json(job)
    loaded = Mycelia.job_spec_from_json(json_text)

    Test.@test loaded.job_name == job.job_name
    Test.@test loaded.cmd == job.cmd
    Test.@test loaded.site == :scg
    Test.@test loaded.partition == "batch"
    Test.@test loaded.account == "PI_EXAMPLE"
    Test.@test loaded.mail_type == ["END", "FAIL"]
end

Test.@testset "Time parsing and site validation" begin
    Test.@test Mycelia.parse_time_limit_seconds("12:34:56") == 45_296
    Test.@test Mycelia.parse_time_limit_seconds("2-01:02:03") == 176_523
    Test.@test isnothing(Mycelia.parse_time_limit_seconds("99:99"))

    bad_scg = Mycelia.JobSpec(
        job_name = "bad-scg",
        cmd = "echo fail",
        site = :scg,
        partition = "nih_s10",
        account = "PI_EXAMPLE",
        time_limit = "10:00:00",
        nodes = 2,
        cpus_per_task = 2,
        mem_gb = 4
    )
    report_scg = Mycelia.validate(bad_scg)
    Test.@test _has_message(report_scg.errors, "nodes=1")

    missing_scg_account = Mycelia.JobSpec(
        job_name = "missing-scg-account",
        cmd = "echo fail",
        site = :scg,
        partition = "nih_s10",
        time_limit = "10:00:00",
        nodes = 1,
        cpus_per_task = 2,
        mem_gb = 4
    )
    report_missing = Mycelia.validate(missing_scg_account)
    Test.@test _has_message(report_missing.errors, "requires account")

    bad_lawrencium = Mycelia.JobSpec(
        job_name = "bad-lr",
        cmd = "echo fail",
        site = :lawrencium,
        time_limit = "04:00:00",
        nodes = 1,
        cpus_per_task = 2,
        mem_gb = 8
    )
    report_l = Mycelia.validate(bad_lawrencium; check_lawrencium_associations = false)
    Test.@test _has_message(report_l.errors, "require partition")
    Test.@test _has_message(report_l.errors, "require account")
end

Test.@testset "NERSC policy validation and salloc quoting" begin
    shared_gpu_bad = Mycelia.JobSpec(
        job_name = "nersc-shared-bad",
        cmd = "python run.py",
        site = :nersc,
        account = "m1234_g",
        qos = "shared",
        constraint = "gpu",
        gpus_per_node = 4,
        time_limit = "04:00:00",
        cpus_per_task = 8,
        ntasks = 1,
        mem_gb = 64
    )
    shared_report = Mycelia.validate(shared_gpu_bad)
    Test.@test _has_message(shared_report.errors, "1 or 2 GPUs")

    debug_too_long = Mycelia.JobSpec(
        job_name = "nersc-debug-too-long",
        cmd = "python debug.py",
        site = :nersc,
        account = "m1234",
        qos = "debug",
        constraint = "cpu",
        time_limit = "01:00:00",
        cpus_per_task = 2,
        mem_gb = 4
    )
    debug_report = Mycelia.validate(debug_too_long)
    Test.@test _has_message(debug_report.errors, "debug QoS max walltime")

    preempt_too_short = Mycelia.JobSpec(
        job_name = "nersc-preempt-too-short",
        cmd = "python preempt.py",
        site = :nersc,
        account = "m1234",
        qos = "preempt",
        constraint = "cpu",
        time_limit = "01:30:00",
        cpus_per_task = 2,
        mem_gb = 4
    )
    preempt_report = Mycelia.validate(preempt_too_short)
    Test.@test _has_message(preempt_report.errors, "minimum walltime")

    interactive_hbm80 = Mycelia.JobSpec(
        job_name = "nersc-interactive-hbm80",
        cmd = "python shell.py",
        site = :nersc,
        account = "m1234_g",
        qos = "interactive",
        constraint = "gpu&hbm80g",
        gpus_per_node = 1,
        time_limit = "01:00:00",
        cpus_per_task = 8,
        mem_gb = 64
    )
    interactive_cmd = Mycelia.render_salloc(interactive_hbm80)
    Test.@test occursin("-C \"gpu&hbm80g\"", interactive_cmd)
end

Test.@testset "Template rendering and fixture checks" begin
    nersc_job = Mycelia.JobSpec(
        job_name = "fixture-nersc",
        cmd = "python train.py --epochs 1",
        site = :nersc,
        account = "m1234_g",
        qos = "regular",
        constraint = "gpu",
        gpus_per_node = 1,
        ntasks = 1,
        cpus_per_task = 16,
        mem_gb = 64,
        time_limit = "04:00:00",
        mail_user = "user@example.org",
        mail_type = "END,FAIL"
    )

    scg_job = Mycelia.JobSpec(
        job_name = "fixture-scg",
        cmd = "julia --project=. workflow.jl",
        site = :scg,
        partition = "nih_s10",
        account = "PI_EXAMPLE",
        time_limit = "24:00:00",
        nodes = 1,
        ntasks = 1,
        cpus_per_task = 12,
        mem_gb = 96,
        mail_user = "user@example.org",
        mail_type = "END,FAIL"
    )

    docker_job = Mycelia.JobSpec(
        job_name = "fixture-docker",
        cmd = "python main.py --input /workspace/input.tsv",
        site = :local,
        time_limit = "01:00:00",
        container_image = "ghcr.io/acme/example:1.0.0",
        container_mounts = ["/tmp/data:/workspace:ro"],
        container_workdir = "/workspace",
        env = Dict("OMP_NUM_THREADS" => "4")
    )

    cloud_job = Mycelia.JobSpec(
        job_name = "fixture-cloudbuild",
        cmd = "python batch.py --quick",
        site = :cloudbuild,
        time_limit = "01:00:00",
        container_image = "gcr.io/my-project/example:sha-deadbeef",
        env = Dict("MY_MODE" => "ci")
    )

    rendered_nersc = Mycelia.render_sbatch(nersc_job)
    rendered_scg = Mycelia.render_sbatch(scg_job)
    rendered_docker = Mycelia.render_docker_run(docker_job)
    rendered_cloud = Mycelia.render_cloudbuild(cloud_job)

    Test.@test !occursin("{{", rendered_nersc)
    Test.@test !occursin("{{", rendered_scg)
    Test.@test !occursin("{{", rendered_docker)
    Test.@test !occursin("{{", rendered_cloud)

    fixture_root = joinpath(dirname(@__DIR__), "fixtures", "job_templates")

    expected_nersc = _normalize_fixture_logdir(
        read(joinpath(fixture_root, "nersc_gpu_regular_1gpu_1task.sbatch"), String))
    expected_scg = _normalize_fixture_logdir(
        read(joinpath(fixture_root, "scg_nih_s10.sbatch"), String))
    expected_docker = read(joinpath(fixture_root, "docker_run.sh"), String)
    expected_cloud = read(joinpath(fixture_root, "cloudbuild_prebuilt.yaml"), String)

    Test.@test rendered_nersc == expected_nersc
    Test.@test rendered_scg == expected_scg
    Test.@test rendered_docker == expected_docker
    Test.@test rendered_cloud == expected_cloud
end

Test.@testset "SLURM wrapper entrypoints" begin
    Test.@testset "lawrencium_sbatch uses env fallbacks with collect executor" begin
        mktempdir() do logdir
            _with_env(Dict(
                "LRC_ACCOUNT" => "pc_example",
                "SLURM_MAIL_USER" => "slurm@example.org"
            )) do
                collector = Mycelia.CollectExecutor()
                outcome = Mycelia.lawrencium_sbatch(
                    job_name = "lawrencium-wrapper",
                    cmd = "echo hello",
                    logdir = logdir,
                    executor = collector
                )

                Test.@test outcome == 1
                Test.@test length(collector.jobs) == 1

                job = only(collector.jobs)
                Test.@test job.site == :lawrencium
                Test.@test job.partition == "lr6"
                Test.@test job.qos == "lr_normal"
                Test.@test job.account == "pc_example"
                Test.@test job.mail_user == "slurm@example.org"
                Test.@test job.output_path == joinpath(logdir, "%j.%x.out")
                Test.@test job.error_path == joinpath(logdir, "%j.%x.err")
            end
        end
    end

    Test.@testset "lawrencium_sbatch falls back to git email" begin
        mktempdir() do logdir
            _with_env(Dict(
                "LRC_ACCOUNT" => "pc_git",
                "SLURM_MAIL_USER" => nothing
            )) do
                _with_git_email("git-lawrencium@example.org") do
                    collector = Mycelia.CollectExecutor()
                    outcome = Mycelia.lawrencium_sbatch(
                        job_name = "lawrencium-git-fallback",
                        cmd = "echo hello",
                        logdir = logdir,
                        executor = collector
                    )

                    Test.@test outcome == 1
                    Test.@test only(collector.jobs).mail_user == "git-lawrencium@example.org"
                end
            end
        end
    end

    Test.@testset "lawrencium_sbatch dry_run delegates to slurm submit executor" begin
        mktempdir() do logdir
            result = Mycelia.lawrencium_sbatch(
                job_name = "lawrencium-dry-run",
                mail_user = "dryrun@example.org",
                account = "pc_dryrun",
                cmd = "echo dry-run",
                logdir = logdir,
                dry_run = true,
                executor = :slurm
            )

            Test.@test result isa Mycelia.SubmitResult
            Test.@test result.ok
            Test.@test result.dry_run
            Test.@test result.backend == :sbatch
            Test.@test result.site == :lawrencium
            Test.@test occursin("#SBATCH --partition=lr6", something(result.artifact_text, ""))
        end
    end

    Test.@testset "scg_sbatch keeps explicit overrides" begin
        mktempdir() do logdir
            collector = Mycelia.CollectExecutor()
            outcome = Mycelia.scg_sbatch(
                job_name = "scg-wrapper",
                mail_user = "override@example.org",
                account = "PI_OVERRIDE",
                partition = "owners",
                time = "12:00:00",
                cpus_per_task = 8,
                mem_gb = 32,
                cmd = "echo hi",
                logdir = logdir,
                executor = collector
            )

            Test.@test outcome == 1
            job = only(collector.jobs)
            Test.@test job.site == :scg
            Test.@test job.partition == "owners"
            Test.@test job.account == "PI_OVERRIDE"
            Test.@test job.mail_user == "override@example.org"
            Test.@test job.time_limit == "12:00:00"
            Test.@test job.cpus_per_task == 8
            Test.@test job.mem_gb == 32.0
        end
    end

    Test.@testset "scg_sbatch falls back to git email" begin
        mktempdir() do logdir
            _with_env(Dict("SLURM_MAIL_USER" => nothing)) do
                _with_git_email("git-scg@example.org") do
                    collector = Mycelia.CollectExecutor()
                    outcome = Mycelia.scg_sbatch(
                        job_name = "scg-git-fallback",
                        account = "PI_GIT",
                        cmd = "echo hi",
                        logdir = logdir,
                        executor = collector
                    )

                    Test.@test outcome == 1
                    Test.@test only(collector.jobs).mail_user == "git-scg@example.org"
                end
            end
        end
    end

    Test.@testset "scg_sbatch dry_run delegates to slurm submit executor" begin
        mktempdir() do logdir
            result = Mycelia.scg_sbatch(
                job_name = "scg-dry-run",
                mail_user = "dryrun@example.org",
                account = "PI_DRYRUN",
                cmd = "echo dry-run",
                logdir = logdir,
                dry_run = true,
                executor = :slurm
            )

            Test.@test result isa Mycelia.SubmitResult
            Test.@test result.ok
            Test.@test result.dry_run
            Test.@test result.backend == :sbatch
            Test.@test result.site == :scg
            Test.@test occursin("sbatch", something(result.submit_command, ""))
            Test.@test occursin("#SBATCH --partition=nih_s10", something(result.artifact_text, ""))
        end
    end

    Test.@testset "nersc wrappers validate account and join command vectors" begin
        _with_env(Dict("NERSC_ACCOUNT" => nothing)) do
            Test.@test_throws ErrorException Mycelia.nersc_sbatch_shared(
                job_name = "missing-account",
                mail_user = "user@example.org",
                cmd = "echo fail"
            )
            Test.@test_throws ErrorException Mycelia.nersc_sbatch(
                job_name = "missing-direct-account",
                mail_user = "user@example.org",
                cmd = "echo fail"
            )
        end

        mktempdir() do logdir
            _with_env(Dict("NERSC_ACCOUNT" => "m1234")) do
                collector = Mycelia.CollectExecutor()
                outcome = Mycelia.nersc_sbatch_shared(
                    job_name = "nersc-shared-wrapper",
                    mail_user = "user@example.org",
                    cmd = "echo shared",
                    logdir = logdir,
                    executor = collector
                )

                Test.@test outcome == 1
                shared_job = only(collector.jobs)
                Test.@test shared_job.site == :nersc
                Test.@test shared_job.qos == "shared"
                Test.@test shared_job.account == "m1234"
                Test.@test shared_job.constraint == "cpu"
                Test.@test shared_job.mem_gb == 2.0
            end
        end

        mktempdir() do logdir
            mktempdir() do scriptdir
                collector = Mycelia.CollectExecutor()
                outcome = Mycelia.nersc_sbatch(
                    job_name = "nersc-command-block",
                    mail_user = "user@example.org",
                    account = "m5678",
                    cmd = ["module load myenv", "python workflow.py --epochs 1"],
                    logdir = logdir,
                    scriptdir = scriptdir,
                    constraint = "gpu",
                    gpus_per_node = 1,
                    executor = collector
                )

                Test.@test outcome == 1
                job = only(collector.jobs)
                Test.@test job.cmd == "module load myenv\npython workflow.py --epochs 1"
                Test.@test job.site == :nersc
                Test.@test job.constraint == "gpu"
                Test.@test job.gpus_per_node == 1
                Test.@test job.output_path == joinpath(logdir, "%j.%x.out")
                Test.@test job.error_path == joinpath(logdir, "%j.%x.err")
            end
        end

        mktempdir() do logdir
            result = Mycelia.nersc_sbatch(
                job_name = "nersc-dry-run",
                mail_user = "user@example.org",
                account = "m3456",
                cmd = "echo dry-run",
                logdir = logdir,
                dry_run = true,
                executor = :slurm
            )

            Test.@test result isa Mycelia.SubmitResult
            Test.@test result.ok
            Test.@test result.dry_run
            Test.@test result.backend == :sbatch
            Test.@test result.site == :nersc
            Test.@test occursin("#SBATCH --constraint=cpu", something(result.artifact_text, ""))
        end
    end

    Test.@testset "lovelace_run writes logs and reports failures" begin
        mktempdir() do logdir
            Test.@test Mycelia.lovelace_run(
                cmd = "printf 'hello from lovelace'",
                logdir = logdir,
                job_name = "lovelace-success"
            )

            out_files = filter(name -> endswith(name, ".out"), readdir(logdir))
            err_files = filter(name -> endswith(name, ".err"), readdir(logdir))
            Test.@test length(out_files) == 1
            Test.@test length(err_files) == 1
            Test.@test read(joinpath(logdir, only(out_files)), String) == "hello from lovelace"
            Test.@test isempty(read(joinpath(logdir, only(err_files)), String))
        end

        mktempdir() do logdir
            Test.@test !Mycelia.lovelace_run(
                cmd = "echo boom >&2; exit 7",
                logdir = logdir,
                job_name = "lovelace-failure"
            )

            err_files = filter(name -> endswith(name, ".err"), readdir(logdir))
            Test.@test length(err_files) == 1
            Test.@test occursin("boom", read(joinpath(logdir, only(err_files)), String))
        end
    end

    Test.@testset "submit_job routes to configured backend" begin
        mktempdir() do logdir
            _with_env(Dict(
                "LRC_ACCOUNT" => "pc_submit",
                "SLURM_MAIL_USER" => "route@example.org"
            )) do
                collector = Mycelia.CollectExecutor()
                outcome = Mycelia.submit_job(
                    site = "lawrencium",
                    job_name = "submit-lawrencium",
                    cmd = "echo routed lawrencium",
                    logdir = logdir,
                    executor = collector
                )

                Test.@test outcome == 1
                Test.@test only(collector.jobs).site == :lawrencium
            end
        end

        mktempdir() do logdir
            collector = Mycelia.CollectExecutor()
            outcome = Mycelia.submit_job(
                site = "SCG",
                job_name = "submit-route",
                mail_user = "user@example.org",
                account = "PI_ROUTED",
                cmd = "echo routed",
                logdir = logdir,
                executor = collector
            )

            Test.@test outcome == 1
            Test.@test only(collector.jobs).site == :scg
        end

        mktempdir() do logdir
            collector = Mycelia.CollectExecutor()
            outcome = Mycelia.submit_job(
                site = "nersc",
                job_name = "submit-nersc",
                mail_user = "user@example.org",
                account = "m9012",
                cmd = "echo routed nersc",
                logdir = logdir,
                executor = collector
            )

            Test.@test outcome == 1
            Test.@test only(collector.jobs).site == :nersc
        end

        mktempdir() do logdir
            Test.@test Mycelia.submit_job(
                site = "lovelace",
                job_name = "submit-lovelace",
                cmd = "printf routed",
                logdir = logdir
            )
        end

        Test.@test_throws ErrorException Mycelia.submit_job(
            site = "unknown-site",
            job_name = "bad-route",
            cmd = "echo nope"
        )
    end
end

Test.@testset "Lawrencium sacctmgr parser" begin
    sample = """
Account|User|Partition|QOS|
pc_project|alice|lr6|lr_normal,lr_debug|
pc_gpu|alice|es1|es_normal|
"""

    table_rows = Mycelia.parse_sacctmgr_pipe_table(sample)
    Test.@test length(table_rows) == 2

    parsed = Mycelia.parse_lawrencium_associations(sample)
    Test.@test length(parsed) == 2
    Test.@test parsed[1]["account"] == "pc_project"
    Test.@test parsed[1]["partition"] == "lr6"
    Test.@test parsed[1]["qos"] == ["lr_normal", "lr_debug"]
    Test.@test parsed[2]["qos"] == ["es_normal"]
end
