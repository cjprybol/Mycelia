import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function _executor_test_write_file(path::String, content::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

function _executor_test_read(path::String)
    return read(path, String)
end

Test.@testset "Execution backend helpers" begin
    Test.@test Mycelia.resolve_executor(nothing) isa Mycelia.LocalExecutor
    Test.@test Mycelia.resolve_executor(:local) isa Mycelia.LocalExecutor
    Test.@test Mycelia.resolve_executor(:slurm) isa Mycelia.SlurmExecutor
    Test.@test Mycelia.resolve_executor(:collect) isa Mycelia.CollectExecutor
    Test.@test Mycelia.resolve_executor(:dry_run) isa Mycelia.DryRunExecutor
    Test.@test Mycelia.resolve_executor(:dryrun) isa Mycelia.DryRunExecutor
    Test.@test_throws ErrorException Mycelia.resolve_executor(:unknown_backend)
    Test.@test_throws ErrorException Mycelia.resolve_executor(1)

    Test.@test Mycelia.command_string(`echo hello`) == "echo hello"
    Test.@test Mycelia.command_string("echo hello") == "echo hello"

    job = Mycelia.build_execution_job(
        cmd = "echo hello",
        job_name = " helper-job ",
        site = "scg",
        time_limit = "00:10:00",
        cpus_per_task = 2,
        mem_gb = 4,
        partition = "nih_s10",
        account = "PI_EXAMPLE",
        mail_user = "user@example.org",
        output_path = " stdout.log ",
        error_path = " stderr.log ",
        workdir = " run-dir ",
        template_name = " scg/cpu "
    )

    Test.@test job isa Mycelia.JobSpec
    Test.@test job.job_name == " helper-job "
    Test.@test job.site == :scg
    Test.@test job.cpus_per_task == 2
    Test.@test job.mem_gb == 4.0
    Test.@test job.output_path == "stdout.log"
    Test.@test job.error_path == "stderr.log"
    Test.@test job.workdir == "run-dir"
    Test.@test job.template_name == "scg/cpu"
end

Test.@testset "Local executor writes redirected output" begin
    mktempdir() do tmpdir
        out_path = joinpath(tmpdir, "logs", "stdout.txt")
        err_path = joinpath(tmpdir, "logs", "stderr.txt")
        workdir = joinpath(tmpdir, "workspace")
        marker_path = joinpath(workdir, "marker.txt")

        job = Mycelia.JobSpec(
            job_name = "local-exec",
            cmd = "pwd > marker.txt; printf 'stdout-line'; printf 'stderr-line' >&2",
            site = :local,
            time_limit = "00:05:00",
            workdir = workdir,
            output_path = out_path,
            error_path = err_path,
            container_image = "ghcr.io/acme/example:1.0.0"
        )

        Test.@test Mycelia.execute(job, Mycelia.LocalExecutor()) == true
        Test.@test isfile(out_path)
        Test.@test isfile(err_path)
        Test.@test isfile(marker_path)
        Test.@test strip(_executor_test_read(marker_path)) == workdir
        Test.@test _executor_test_read(out_path) == "stdout-line"
        Test.@test _executor_test_read(err_path) == "stderr-line"
    end
end

Test.@testset "Executor collection and dry-run submission" begin
    job = Mycelia.JobSpec(
        job_name = "collector-test",
        cmd = "echo hi",
        site = :scg,
        time_limit = "00:10:00",
        partition = "nih_s10",
        account = "PI_EXAMPLE",
        nodes = 1,
        ntasks = 1,
        cpus_per_task = 1,
        mem_gb = 4
    )

    collector = Mycelia.CollectExecutor()
    idx = Mycelia.execute(job, collector)
    Test.@test idx == 1
    Test.@test length(collector.jobs) == 1
    Test.@test collector.jobs[1] isa Mycelia.JobSpec

    dry_executor = Mycelia.DryRunExecutor()
    result = Mycelia.execute(job, dry_executor)
    Test.@test result isa Mycelia.SubmitResult
    Test.@test result.dry_run
    Test.@test result.ok
    Test.@test length(dry_executor.jobs) == 1
    Test.@test length(dry_executor.results) == 1
end

Test.@testset "Wrapper collection paths remain routable through executor infrastructure" begin
    mktempdir() do tmpdir
        fastq1 = _executor_test_write_file(joinpath(tmpdir, "reads_1.fastq"), "@r1\nACGT\n+\n!!!!\n")
        fastq2 = _executor_test_write_file(joinpath(tmpdir, "reads_2.fastq"), "@r2\nACGT\n+\n!!!!\n")
        query_fasta = _executor_test_write_file(joinpath(tmpdir, "query.faa"), ">q1\nMKT\n")
        reference_fasta = _executor_test_write_file(joinpath(tmpdir, "ref.faa"), ">r1\nMKT\n")
        query_sig = _executor_test_write_file(joinpath(tmpdir, "query.sig"), "signature")
        db_sig = _executor_test_write_file(joinpath(tmpdir, "db.sig"), "signature")
        cluster_fasta = _executor_test_write_file(joinpath(tmpdir, "cluster.fa"), ">a\nACGT\n")
        genome_a = _executor_test_write_file(joinpath(tmpdir, "genome_a.fa"), ">chr1\nACGT\n")
        genome_b = _executor_test_write_file(joinpath(tmpdir, "genome_b.fa"), ">chr1\nAGGT\n")

        assembly_exec = Mycelia.CollectExecutor()
        megahit_out = Mycelia.run_megahit(
            fastq1 = fastq1,
            fastq2 = fastq2,
            outdir = joinpath(tmpdir, "megahit"),
            executor = assembly_exec,
            site = :scg
        )
        Test.@test length(assembly_exec.jobs) == 1
        Test.@test assembly_exec.jobs[1].site == :scg
        Test.@test occursin("megahit", assembly_exec.jobs[1].cmd)
        Test.@test megahit_out.contigs == joinpath(tmpdir, "megahit", "final.contigs.fa")

        classification_exec = Mycelia.CollectExecutor()
        sourmash_out = Mycelia.run_sourmash_search(
            query_sig = query_sig,
            database_sig = db_sig,
            outdir = joinpath(tmpdir, "sourmash"),
            executor = classification_exec,
            site = :scg
        )
        Test.@test length(classification_exec.jobs) == 1
        Test.@test occursin("sourmash search", classification_exec.jobs[1].cmd)
        Test.@test endswith(sourmash_out.results_csv, "_search_results.csv")

        align_exec = Mycelia.CollectExecutor()
        blast_out = Mycelia.run_blastp_search(
            query_fasta = query_fasta,
            reference_fasta = reference_fasta,
            output_dir = joinpath(tmpdir, "blast"),
            executor = align_exec,
            site = :scg
        )
        Test.@test length(align_exec.jobs) == 1
        Test.@test occursin("blastp", align_exec.jobs[1].cmd)
        Test.@test endswith(blast_out, "_blastp_results.tsv")

        cluster_exec = Mycelia.CollectExecutor()
        cluster_out = Mycelia.mmseqs_easy_cluster(
            fasta = cluster_fasta,
            output = joinpath(tmpdir, "mmseqs_cluster"),
            tmp = joinpath(tmpdir, "mmseqs_tmp"),
            executor = cluster_exec,
            site = :scg
        )
        Test.@test length(cluster_exec.jobs) == 1
        Test.@test occursin("easy-cluster", cluster_exec.jobs[1].cmd)
        Test.@test endswith(cluster_out, "_cluster.tsv")

        pangenome_exec = Mycelia.CollectExecutor()
        cactus_out = Mycelia.construct_pangenome_cactus(
            [genome_a, genome_b],
            ["REF", "ALT"],
            joinpath(tmpdir, "cactus"),
            "REF";
            output_formats = ["gfa", "vcf"],
            executor = pangenome_exec,
            site = :scg
        )
        Test.@test length(pangenome_exec.jobs) == 1
        Test.@test occursin("cactus-pangenome", pangenome_exec.jobs[1].cmd)
        Test.@test haskey(cactus_out, "gfa")
        Test.@test haskey(cactus_out, "vcf")
    end
end

Test.@testset "SLURM wrapper executor support and environment fallbacks" begin
    collector = Mycelia.CollectExecutor()
    outcome = Mycelia.scg_sbatch(
        job_name = "executor-test",
        mail_user = "user@example.org",
        account = "PI_EXAMPLE",
        cmd = "echo hello",
        executor = collector
    )
    Test.@test outcome == 1
    Test.@test length(collector.jobs) == 1
    Test.@test collector.jobs[1].site == :scg

    env_collector = Mycelia.CollectExecutor()
    env_outcome = Base.withenv(
        "SLURM_MAIL_USER" => "slurm@example.org",
        "LRC_ACCOUNT" => "LRC_PROJECT"
    ) do
        Mycelia.lawrencium_sbatch(
            job_name = "lawrencium-env",
            cmd = "echo lr",
            executor = env_collector
        )
    end
    Test.@test env_outcome == 1
    Test.@test env_collector.jobs[1].site == :lawrencium
    Test.@test env_collector.jobs[1].account == "LRC_PROJECT"
    Test.@test env_collector.jobs[1].mail_user == "slurm@example.org"

    dry_run_result = Base.withenv(
        "SLURM_MAIL_USER" => "slurm@example.org",
        "SCG_ACCOUNT" => "PI_ENV"
    ) do
        Mycelia.scg_sbatch(
            job_name = "scg-dry-run",
            cmd = "echo scg",
            dry_run = true,
            executor = :slurm
        )
    end
    Test.@test dry_run_result isa Mycelia.SubmitResult
    Test.@test dry_run_result.ok
    Test.@test dry_run_result.dry_run
    Test.@test dry_run_result.backend == :sbatch

    nersc_collector = Mycelia.CollectExecutor()
    shared_outcome = Base.withenv("NERSC_ACCOUNT" => "m1234") do
        Mycelia.nersc_sbatch_shared(
            job_name = "shared-job",
            mail_user = "user@example.org",
            cmd = "echo nersc shared",
            executor = nersc_collector
        )
    end
    Test.@test shared_outcome == 1
    Test.@test nersc_collector.jobs[1].site == :nersc
    Test.@test nersc_collector.jobs[1].qos == "shared"
    Test.@test nersc_collector.jobs[1].account == "m1234"

    vector_cmd_collector = Mycelia.CollectExecutor()
    nersc_outcome = Base.withenv("NERSC_ACCOUNT" => "m5678") do
        Mycelia.nersc_sbatch(
            job_name = "nersc-vector-cmd",
            mail_user = "user@example.org",
            cmd = ["echo first", "echo second"],
            executor = vector_cmd_collector
        )
    end
    Test.@test nersc_outcome == 1
    Test.@test vector_cmd_collector.jobs[1].cmd == "echo first\necho second"
    Test.@test_throws ErrorException Mycelia.nersc_sbatch_shared(
        job_name = "missing-account",
        mail_user = "user@example.org",
        cmd = "echo fail",
        executor = Mycelia.CollectExecutor()
    )
end

Test.@testset "Dispatch helpers and local run wrappers" begin
    mktempdir() do tmpdir
        Test.@test Mycelia.lovelace_run(
            cmd = "printf 'lovelace-ok'",
            logdir = tmpdir,
            job_name = "lovelace-pass"
        )
        Test.@test any(name -> occursin("lovelace-pass", name), readdir(tmpdir))

        Test.@test !Mycelia.lovelace_run(
            cmd = "exit 1",
            logdir = tmpdir,
            job_name = "lovelace-fail"
        )

        routed = Base.withenv(
            "SLURM_MAIL_USER" => "router@example.org",
            "SCG_ACCOUNT" => "PI_ROUTED"
        ) do
            Mycelia.submit_job(
                site = "scg",
                job_name = "submit-job-scg",
                cmd = "echo routed",
                dry_run = true,
                executor = :slurm
            )
        end
        Test.@test routed isa Mycelia.SubmitResult
        Test.@test routed.ok
        Test.@test routed.backend == :sbatch

        Test.@test Mycelia.submit_job(
            site = "lovelace",
            cmd = "printf 'dispatch-ok'",
            logdir = tmpdir,
            job_name = "submit-job-lovelace"
        )
        Test.@test_throws ErrorException Mycelia.submit_job(site = "unknown", cmd = "echo no")
    end
end

Test.@testset "Submission helper branches" begin
    docker_job = Mycelia.JobSpec(
        job_name = "docker-dry-run",
        cmd = "python main.py --input /workspace/input.tsv",
        site = :local,
        time_limit = "01:00:00",
        container_image = "ghcr.io/acme/example:1.0.0",
        container_mounts = ["/tmp/data:/workspace:ro"],
        container_workdir = "/workspace",
        env = Dict("OMP_NUM_THREADS" => "4")
    )

    cloud_job = Mycelia.JobSpec(
        job_name = "cloudbuild-write-only",
        cmd = "python batch.py --quick",
        site = :cloudbuild,
        time_limit = "01:00:00",
        container_image = "gcr.io/my-project/example:sha-deadbeef",
        env = Dict("MY_MODE" => "ci")
    )

    interactive_job = Mycelia.JobSpec(
        job_name = "interactive-scg",
        cmd = "bash",
        site = :scg,
        partition = "interactive",
        account = "PI_INTERACTIVE",
        time_limit = "00:30:00",
        nodes = 1,
        ntasks = 1,
        cpus_per_task = 4,
        mem_gb = 8
    )

    bad_job = Mycelia.JobSpec(
        job_name = "bad-scg",
        cmd = "echo fail",
        site = :scg,
        partition = "nih_s10",
        time_limit = "00:10:00",
        nodes = 2,
        cpus_per_task = 1,
        mem_gb = 1
    )

    docker_dry = Mycelia.submit(docker_job)
    Test.@test docker_dry.ok
    Test.@test docker_dry.dry_run
    Test.@test docker_dry.backend == :docker
    Test.@test occursin("docker run --rm", something(docker_dry.artifact_text, ""))

    mktempdir() do tmpdir
        cloud_path = joinpath(tmpdir, "cloudbuild.yaml")
        cloud_result = Mycelia.submit(cloud_job; dry_run = false, path = cloud_path)
        Test.@test cloud_result.ok
        Test.@test !cloud_result.dry_run
        Test.@test cloud_result.backend == :cloudbuild
        Test.@test cloud_result.artifact_path == cloud_path
        Test.@test isfile(cloud_path)
        Test.@test any(warning -> occursin("execute_cloudbuild=true", warning), cloud_result.warnings)
    end

    interactive_dry = Mycelia.submit(interactive_job)
    Test.@test interactive_dry.ok
    Test.@test interactive_dry.backend == :salloc
    Test.@test interactive_dry.dry_run
    Test.@test occursin("salloc", something(interactive_dry.submit_command, ""))

    interactive_write_only = Mycelia.submit(interactive_job; dry_run = false)
    Test.@test interactive_write_only.ok
    Test.@test !interactive_write_only.dry_run
    Test.@test interactive_write_only.backend == :salloc
    Test.@test any(
        warning -> occursin("execute_interactive=true", warning),
        interactive_write_only.warnings
    )

    invalid_result = Mycelia.submit(bad_job)
    Test.@test !invalid_result.ok
    Test.@test invalid_result.backend == :validation
    Test.@test !isempty(invalid_result.errors)
end

Test.@testset "Scheduler utility helpers" begin
    Test.@test Mycelia.list_lawrencium_associations(execute = false) ==
               "sacctmgr show association -p user=\$USER"
    Test.@test Mycelia.list_lawrencium_qos_limits(execute = false) ==
               "sacctmgr show qos -p format=name,maxtres,maxwall,mintres"

    malformed = "Account|Partition|\nvalid|lr6|\nbad|extra|field|\n"
    parsed_rows = Mycelia.parse_sacctmgr_pipe_table(malformed)
    Test.@test length(parsed_rows) == 1
    Test.@test parsed_rows[1]["account"] == "valid"

    assoc_sample = """
Acct|PartitionName|DefaultQOS|
pc_project|lr6|lr_normal,lr_debug|
pc_gpu|es1||
"""
    parsed_assoc = Mycelia.parse_lawrencium_associations(assoc_sample)
    Test.@test parsed_assoc[1]["account"] == "pc_project"
    Test.@test parsed_assoc[1]["partition"] == "lr6"
    Test.@test parsed_assoc[1]["qos"] == ["lr_normal", "lr_debug"]
    Test.@test parsed_assoc[2]["qos"] == String[]

    summary_buffer = IOBuffer()
    summary = Mycelia.summarize_job("123"; io = summary_buffer)
    summary_output = String(take!(summary_buffer))
    Test.@test isa(summary, Dict{String, String})
    Test.@test all(key -> haskey(summary, key), ["scontrol", "sacct", "sstat", "seff"])
    Test.@test occursin("=== Job Summary: 123 ===", summary_output)
    Test.@test_throws ErrorException Mycelia.summarize_job("invalid-job-id")

    warning_buffer = IOBuffer()
    Test.@test isnothing(Mycelia.nersc_login_node_limits_warning(io = warning_buffer))
    Test.@test occursin("56 GB RAM", String(take!(warning_buffer)))
end
