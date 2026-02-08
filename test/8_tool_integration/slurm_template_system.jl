import Test
import Mycelia

function _has_message(messages::Vector{String}, needle::String)
    return any(msg -> occursin(needle, msg), messages)
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

    expected_nersc = read(joinpath(fixture_root, "nersc_gpu_regular_1gpu_1task.sbatch"), String)
    expected_scg = read(joinpath(fixture_root, "scg_nih_s10.sbatch"), String)
    expected_docker = read(joinpath(fixture_root, "docker_run.sh"), String)
    expected_cloud = read(joinpath(fixture_root, "cloudbuild_prebuilt.yaml"), String)

    Test.@test rendered_nersc == expected_nersc
    Test.@test rendered_scg == expected_scg
    Test.@test rendered_docker == expected_docker
    Test.@test rendered_cloud == expected_cloud
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
