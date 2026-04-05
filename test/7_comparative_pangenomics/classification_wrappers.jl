import Test
import Mycelia

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
