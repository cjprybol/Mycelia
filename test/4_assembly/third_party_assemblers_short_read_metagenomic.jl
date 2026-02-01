# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/third_party_assemblers_short_read_metagenomic.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/third_party_assemblers_short_read_metagenomic.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

import Test
import Mycelia
import StableRNGs
import BioSequences
import FASTX

Test.@testset "Short Read Metagenomic Assembly" begin
    threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

    mktempdir() do dir
        meta_ref_fasta = joinpath(dir, "metagenomic_ref.fasta")

        rng1 = StableRNGs.StableRNG(456)
        rng2 = StableRNGs.StableRNG(457)
        genome1 = BioSequences.randdnaseq(rng1, 4_000)
        genome2 = BioSequences.randdnaseq(rng2, 3_000)
        Mycelia.write_fasta(
            outfile = meta_ref_fasta,
            records = [
                FASTX.FASTA.Record("test_metagenomic_genome_1", genome1),
                FASTX.FASTA.Record("test_metagenomic_genome_2", genome2)
            ]
        )

        simulated_reads = Mycelia.simulate_illumina_reads(
            fasta = meta_ref_fasta, coverage = 10, rndSeed = 456, quiet = true)
        fastq1 = simulated_reads.forward_reads
        fastq2 = simulated_reads.reverse_reads

        Test.@testset "MEGAHIT" begin
            outdir = joinpath(dir, "megahit_assembly")
            isdir(outdir) && rm(outdir, recursive = true, force = true)
            try
                result = Mycelia.run_megahit(
                    fastq1 = fastq1,
                    fastq2 = fastq2,
                    outdir = outdir,
                    k_list = "21",
                    threads = threads
                )
                Test.@test result.outdir == outdir
                Test.@test isfile(result.contigs)
            finally
                rm(outdir, recursive = true, force = true)
            end
        end

        Test.@testset "metaSPAdes" begin
            outdir = joinpath(dir, "metaspades_assembly")
            isdir(outdir) && rm(outdir, recursive = true, force = true)
            try
                result = Mycelia.run_metaspades(
                    fastq1 = fastq1,
                    fastq2 = fastq2,
                    outdir = outdir,
                    k_list = "21",
                    threads = threads
                )
                Test.@test result.outdir == outdir
                Test.@test isfile(result.contigs)
                Test.@test isfile(result.scaffolds)
            finally
                rm(outdir, recursive = true, force = true)
            end
        end
    end
end
