import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function with_fake_conda_runner(f::Function)
    runner_path = Mycelia.CONDA_RUNNER
    runner_dir = dirname(runner_path)
    backup_dir = mktempdir()
    backup_path = joinpath(backup_dir, "conda")

    if isfile(runner_path)
        mv(runner_path, backup_path; force = true)
    end

    mkpath(runner_dir)
    open(runner_path, "w") do io
        write(io, raw"""
        #!/bin/sh
        set -eu

        if [ "${5:-}" != "claMLST" ]; then
          echo "unexpected command: $*" >&2
          exit 1
        fi

        case "${6:-}" in
          import)
            mkdir -p "$(dirname "$7")"
            : > "$7"
            ;;
          search)
            printf 'scheme\tst\nfake\t1\n'
            ;;
          *)
            echo "unexpected subcommand: ${6:-}" >&2
            exit 1
            ;;
        esac
        """)
    end
    chmod(runner_path, 0o755)

    try
        return f()
    finally
        rm(runner_path; force = true)
        if isfile(backup_path)
            mv(backup_path, runner_path; force = true)
        end
        rm(backup_dir; recursive = true, force = true)
    end
end

Test.@testset "Classification wrapper executor coverage" begin
    temp_dir = mktempdir()
    try
        genome_file = joinpath(temp_dir, "genome.fasta")
        open(genome_file, "w") do io
            println(io, ">contig1")
            println(io, "ATGCATGC")
        end

        db_path = joinpath(temp_dir, "custom", "claMLSTDB")
        outdir = joinpath(temp_dir, "clamlst")
        executor = Mycelia.CollectExecutor()

        output_tsv = Mycelia.run_clamlst(
            genome_file;
            db_path = db_path,
            outdir = outdir,
            threads = 3,
            executor = executor,
            site = :scg,
            job_name = "clamlst-collect"
        )

        Test.@test output_tsv == joinpath(outdir, "genome_clamlst.tsv")
        Test.@test length(executor.jobs) == 1
        Test.@test executor.jobs[1].job_name == "clamlst-collect"
        Test.@test executor.jobs[1].site == :scg
        Test.@test executor.jobs[1].cpus_per_task == 3
        Test.@test occursin("claMLST import", executor.jobs[1].cmd)
        Test.@test occursin("claMLST search", executor.jobs[1].cmd)
        Test.@test occursin(db_path, executor.jobs[1].cmd)
        Test.@test occursin(genome_file, executor.jobs[1].cmd)
    finally
        rm(temp_dir; recursive = true, force = true)
    end
end

Test.@testset "Classification wrapper cache and compatibility" begin
    temp_dir = mktempdir()
    try
        missing_genome = joinpath(temp_dir, "missing.fasta")
        cached_outdir = joinpath(temp_dir, "cached")
        mkpath(cached_outdir)
        cached_output = joinpath(cached_outdir, "missing_clamlst.tsv")
        write(cached_output, "scheme\tst\n")

        Test.@test Mycelia.run_clamlst(
            missing_genome;
            outdir = cached_outdir
        ) == cached_output

        alias_genome = joinpath(temp_dir, "alias.fasta")
        open(alias_genome, "w") do io
            println(io, ">contig1")
            println(io, "ATGCATGC")
        end

        alias_db_path = joinpath(temp_dir, "alias", "claMLSTDB")
        alias_executor = Mycelia.CollectExecutor()
        alias_output = Test.@test_deprecated Mycelia.run_clamlst(
            alias_genome;
            db_dir = alias_db_path,
            outdir = joinpath(temp_dir, "alias_out"),
            executor = alias_executor
        )

        Test.@test alias_output == joinpath(temp_dir, "alias_out", "alias_clamlst.tsv")
        Test.@test length(alias_executor.jobs) == 1
        Test.@test occursin(alias_db_path, alias_executor.jobs[1].cmd)
    finally
        rm(temp_dir; recursive = true, force = true)
    end
end

Test.@testset "Classification wrapper uncached local execution" begin
    temp_dir = mktempdir()
    try
        genome_file = joinpath(temp_dir, "local.fasta")
        open(genome_file, "w") do io
            println(io, ">contig1")
            println(io, "ATGCATGC")
        end

        db_path = joinpath(temp_dir, "local", "claMLSTDB")
        outdir = joinpath(temp_dir, "local_out")
        output_tsv = with_fake_conda_runner() do
            Mycelia.run_clamlst(
                genome_file;
                db_path = db_path,
                outdir = outdir
            )
        end

        Test.@test output_tsv == joinpath(outdir, "local_clamlst.tsv")
        Test.@test isfile(output_tsv)
        Test.@test isfile(db_path)
        Test.@test read(output_tsv, String) == "scheme\tst\nfake\t1\n"
    finally
        rm(temp_dir; recursive = true, force = true)
    end
end

Test.@testset "Classification wrapper executor force-update and validation" begin
    temp_dir = mktempdir()
    try
        plain_input = joinpath(temp_dir, "compressed.fasta")
        open(plain_input, "w") do io
            println(io, ">contig1")
            println(io, "ATGCATGC")
        end
        gz_input = joinpath(temp_dir, "compressed.fasta.gz")
        open(gz_input, "w") do io
            run(pipeline(`gzip -c $(plain_input)`, stdout = io))
        end

        forced_db_path = joinpath(temp_dir, "forced", "claMLSTDB")
        mkpath(forced_db_path)
        force_executor = Mycelia.CollectExecutor()
        forced_output = Mycelia.run_clamlst(
            gz_input;
            db_path = forced_db_path,
            outdir = joinpath(temp_dir, "forced_out"),
            force_db_update = true,
            executor = force_executor
        )

        Test.@test forced_output == joinpath(temp_dir, "forced_out", "compressed_clamlst.tsv")
        Test.@test length(force_executor.jobs) == 1
        Test.@test occursin("claMLST import", force_executor.jobs[1].cmd)
        Test.@test occursin("--force", force_executor.jobs[1].cmd)
        Test.@test occursin("gunzip -c", force_executor.jobs[1].cmd)
        Test.@test occursin("mktemp", force_executor.jobs[1].cmd)

        mismatch_executor = Mycelia.CollectExecutor()
        mismatch_genome = joinpath(temp_dir, "mismatch.fasta")
        open(mismatch_genome, "w") do io
            println(io, ">contig1")
            println(io, "ATGCATGC")
        end

        Test.@test_throws ErrorException Mycelia.run_clamlst(
            mismatch_genome;
            db_path = joinpath(temp_dir, "path_a"),
            db_dir = joinpath(temp_dir, "path_b"),
            outdir = joinpath(temp_dir, "mismatch_out"),
            executor = mismatch_executor
        )
    finally
        rm(temp_dir; recursive = true, force = true)
    end
end
