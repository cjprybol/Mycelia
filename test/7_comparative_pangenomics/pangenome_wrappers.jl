# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/7_comparative_pangenomics/pangenome_wrappers.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/7_comparative_pangenomics/pangenome_wrappers.jl", "test/7_comparative_pangenomics", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

Test.@testset "Pangenome wrapper validation" begin
    # PGGB should fail fast on missing genome files
    genomes = ["missing1.fasta", "missing2.fasta"]
    outdir = mktempdir()
    Test.@test_throws ErrorException Mycelia.construct_pangenome_pggb(genomes, outdir)

    # Cactus should catch mismatched names before touching the filesystem
    genome_files = ["g1.fasta", "g2.fasta"]
    genome_names = ["sample1"]
    Test.@test_throws ErrorException Mycelia.construct_pangenome_cactus(
        genome_files, genome_names, outdir, "sample1")

    # vg deconstruct validation
    Test.@test_throws ErrorException Mycelia.call_variants_from_pggb_graph("missing.gfa", "ref")
end

Test.@testset "Pangenome wrapper executor collection" begin
    temp_dir = mktempdir()
    try
        genome_a = joinpath(temp_dir, "genome_a.fasta")
        genome_b = joinpath(temp_dir, "genome_b.fasta")
        open(genome_a, "w") do io
            println(io, ">genome_a")
            println(io, "ATGCGATGCA")
        end
        open(genome_b, "w") do io
            println(io, ">genome_b")
            println(io, "ATGCGTTTAA")
        end

        executor = Mycelia.CollectExecutor()
        output_dir = joinpath(temp_dir, "pggb")
        gfa_path = Mycelia.construct_pangenome_pggb(
            [genome_a, genome_b],
            output_dir;
            threads = 4,
            additional_args = ["--poa-length-target", "1000"],
            executor = executor,
            site = :scg,
            job_name = "pggb-collect"
        )

        Test.@test gfa_path == joinpath(output_dir, "pangenome.gfa")
        Test.@test length(executor.jobs) == 1
        Test.@test executor.jobs[1].job_name == "pggb-collect"
        Test.@test executor.jobs[1].site == :scg
        Test.@test executor.jobs[1].cpus_per_task == 4
        Test.@test occursin("samtools faidx", executor.jobs[1].cmd)
        Test.@test occursin("pggb", executor.jobs[1].cmd)
        Test.@test occursin("--poa-length-target 1000", executor.jobs[1].cmd)
    finally
        isdir(temp_dir) && rm(temp_dir, recursive = true, force = true)
    end
end

Test.@testset "Cactus wrapper executor collection" begin
    temp_dir = mktempdir()
    try
        reference = joinpath(temp_dir, "reference.fasta")
        sample = joinpath(temp_dir, "sample.fasta")
        open(reference, "w") do io
            println(io, ">reference")
            println(io, "ATGCGATGCA")
        end
        open(sample, "w") do io
            println(io, ">sample")
            println(io, "ATGCGTTTAA")
        end

        executor = Mycelia.CollectExecutor()
        output_dir = joinpath(temp_dir, "cactus")
        output_files = Mycelia.construct_pangenome_cactus(
            [reference, sample],
            ["reference", "sample"],
            output_dir,
            "reference";
            max_cores = 6,
            max_memory_gb = 24,
            output_formats = ["gfa", "odgi"],
            executor = executor,
            site = :nersc,
            job_name = "cactus-collect"
        )

        Test.@test output_files == Dict(
            "gfa" => joinpath(output_dir, "pangenome.gfa"),
            "odgi" => joinpath(output_dir, "pangenome.odgi")
        )
        Test.@test isfile(joinpath(output_dir, "cactus_config.txt"))
        Test.@test length(executor.jobs) == 1
        Test.@test executor.jobs[1].job_name == "cactus-collect"
        Test.@test executor.jobs[1].site == :nersc
        Test.@test executor.jobs[1].cpus_per_task == 6
        Test.@test executor.jobs[1].mem_gb == 24
        Test.@test occursin("podman-hpc", executor.jobs[1].cmd)
        Test.@test occursin("--reference reference", executor.jobs[1].cmd)
        Test.@test occursin("--gfa", executor.jobs[1].cmd)
        Test.@test occursin("--odgi", executor.jobs[1].cmd)
    finally
        isdir(temp_dir) && rm(temp_dir, recursive = true, force = true)
    end
end

Test.@testset "Pangenome wrapper index validation" begin
    temp_dir = mktempdir()
    try
        graph_file = joinpath(temp_dir, "graph.gfa")
        open(graph_file, "w") do io
            println(io, "H\tVN:Z:1.0")
        end

        Test.@test isempty(Mycelia.index_pangenome_graph(graph_file; index_types = ["mystery"]))
        Test.@test_throws ErrorException Mycelia.convert_gfa_to_vg_format(joinpath(temp_dir, "missing.gfa"))
        Test.@test_throws ErrorException Mycelia.index_pangenome_graph(joinpath(temp_dir, "missing.vg"))
    finally
        isdir(temp_dir) && rm(temp_dir, recursive = true, force = true)
    end
end
