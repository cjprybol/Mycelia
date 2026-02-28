import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function _write_file(path::String, content::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

Test.@testset "Executor Resolution" begin
    Test.@test Mycelia.resolve_executor(nothing) isa Mycelia.LocalExecutor
    Test.@test Mycelia.resolve_executor(:local) isa Mycelia.LocalExecutor
    Test.@test Mycelia.resolve_executor(:slurm) isa Mycelia.SlurmExecutor
    Test.@test Mycelia.resolve_executor(:collect) isa Mycelia.CollectExecutor
    Test.@test Mycelia.resolve_executor(:dry_run) isa Mycelia.DryRunExecutor
    Test.@test_throws ErrorException Mycelia.resolve_executor(:unknown_backend)
    Test.@test !isdefined(Mycelia, :Execution)
end

Test.@testset "Executor Dispatch" begin
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
    Test.@test length(dry_executor.jobs) == 1
    Test.@test length(dry_executor.results) == 1
end

Test.@testset "Wrapper Collection Paths" begin
    tmpdir = mktempdir()

    fastq1 = _write_file(joinpath(tmpdir, "reads_1.fastq"), "@r1\nACGT\n+\n!!!!\n")
    fastq2 = _write_file(joinpath(tmpdir, "reads_2.fastq"), "@r2\nACGT\n+\n!!!!\n")
    query_fasta = _write_file(joinpath(tmpdir, "query.faa"), ">q1\nMKT\n")
    reference_fasta = _write_file(joinpath(tmpdir, "ref.faa"), ">r1\nMKT\n")
    query_sig = _write_file(joinpath(tmpdir, "query.sig"), "signature")
    db_sig = _write_file(joinpath(tmpdir, "db.sig"), "signature")
    cluster_fasta = _write_file(joinpath(tmpdir, "cluster.fa"), ">a\nACGT\n")
    genome_a = _write_file(joinpath(tmpdir, "genome_a.fa"), ">chr1\nACGT\n")
    genome_b = _write_file(joinpath(tmpdir, "genome_b.fa"), ">chr1\nAGGT\n")

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

Test.@testset "Backward Compatibility Default Local Path" begin
    tmpdir = mktempdir()
    fastq = _write_file(joinpath(tmpdir, "reads.fastq"), "@r1\nACGT\n+\n!!!!\n")
    outdir = joinpath(tmpdir, "megahit_local")
    _write_file(joinpath(outdir, "final.contigs.fa"), ">k141_1\nACGT\n")
    _write_file(joinpath(outdir, "final.contigs.fastg"), ">FASTG:begin\n")
    _write_file(joinpath(outdir, "final.contigs.fastg.gfa"), "H\tVN:Z:1.0\n")

    result = Mycelia.run_megahit(fastq1 = fastq, outdir = outdir)
    Test.@test result.contigs == joinpath(outdir, "final.contigs.fa")
    Test.@test result.fastg == joinpath(outdir, "final.contigs.fastg")
    Test.@test result.gfa == joinpath(outdir, "final.contigs.fastg.gfa")
end

Test.@testset "SLURM Wrapper Executor Support" begin
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
end
