# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/third_party_assemblers_short_read_isolate.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/third_party_assemblers_short_read_isolate.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

import Test
import Mycelia
import StableRNGs
import BioSequences
import FASTX

Test.@testset "Short Read Isolate Assembly" begin
    threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

    mktempdir() do dir
        ref_fasta = joinpath(dir, "isolate_ref.fasta")
        rng = StableRNGs.StableRNG(123)
        isolate_genome = BioSequences.randdnaseq(rng, 5_000)
        Mycelia.write_fasta(outfile=ref_fasta, records=[FASTX.FASTA.Record("test_isolate_genome", isolate_genome)])

        simulated_reads = Mycelia.simulate_illumina_reads(fasta=ref_fasta, coverage=10, rndSeed=123, quiet=true)
        fastq1 = simulated_reads.forward_reads
        fastq2 = simulated_reads.reverse_reads

        Test.@testset "SPAdes" begin
            spades_outdir = joinpath(dir, "spades_assembly")
            isdir(spades_outdir) && rm(spades_outdir, recursive=true, force=true)
            try
                result = Mycelia.run_spades(
                    fastq1=fastq1,
                    fastq2=fastq2,
                    outdir=spades_outdir,
                    k_list="21",
                    threads=threads
                )
                Test.@test result.outdir == spades_outdir
                Test.@test isfile(result.contigs)
                Test.@test isfile(result.scaffolds)
            finally
                rm(spades_outdir, recursive=true, force=true)
            end
        end

        Test.@testset "SKESA" begin
            skesa_outdir = joinpath(dir, "skesa_assembly")
            isdir(skesa_outdir) && rm(skesa_outdir, recursive=true, force=true)
            try
                result = Mycelia.run_skesa(
                    fastq1=fastq1,
                    fastq2=fastq2,
                    outdir=skesa_outdir,
                    threads=threads
                )
                Test.@test result.outdir == skesa_outdir
                Test.@test isfile(result.contigs)
            finally
                rm(skesa_outdir, recursive=true, force=true)
            end
        end

        Test.@testset "Velvet" begin
            velvet_outdir = joinpath(dir, "velvet_assembly")
            isdir(velvet_outdir) && rm(velvet_outdir, recursive=true, force=true)
            try
                result = Mycelia.run_velvet(fastq1; fastq2=fastq2, outdir=velvet_outdir, k=21)
                Test.@test result.outdir == velvet_outdir
                Test.@test isfile(result.contigs)
            finally
                rm(velvet_outdir, recursive=true, force=true)
            end
        end
    end
end

