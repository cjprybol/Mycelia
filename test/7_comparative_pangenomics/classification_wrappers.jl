import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
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
