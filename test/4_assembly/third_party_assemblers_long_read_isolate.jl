# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/third_party_assemblers_long_read_isolate.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/third_party_assemblers_long_read_isolate.jl", "test/4_assembly", execute=false)'
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

Test.@testset "Long Read Isolate Assembly" begin
    threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

    mktempdir() do dir
        Test.@testset "Flye" begin
            ref_fasta = joinpath(dir, "flye_isolate_ref.fasta")
            rng = StableRNGs.StableRNG(456)
            genome = BioSequences.randdnaseq(rng, 10_000)
            Mycelia.write_fasta(outfile=ref_fasta, records=[FASTX.FASTA.Record("flye_isolate_genome", genome)])

            reads_gz = Mycelia.simulate_nanopore_reads(fasta=ref_fasta, quantity="15x", quiet=true)
            reads_fastq = joinpath(dir, "flye_reads.fq")
            run(pipeline(`gunzip -c $(reads_gz)`, reads_fastq))

            outdir = joinpath(dir, "flye_assembly")
            isdir(outdir) && rm(outdir, recursive=true, force=true)
            try
                result = Mycelia.run_flye(
                    fastq=reads_fastq,
                    outdir=outdir,
                    genome_size="10k",
                    read_type="nano-hq",
                    min_overlap=1000,
                    threads=threads
                )
                Test.@test result.outdir == outdir
                Test.@test isfile(result.assembly)
            finally
                rm(outdir, recursive=true, force=true)
            end
        end

        Test.@testset "Hifiasm" begin
            ref_fasta = joinpath(dir, "hifiasm_reference.fasta")
            rng = StableRNGs.StableRNG(333)
            genome = BioSequences.randdnaseq(rng, 15_000)
            Mycelia.write_fasta(outfile=ref_fasta, records=[FASTX.FASTA.Record("hifiasm_test_genome", genome)])

            reads_gz = Mycelia.simulate_pacbio_reads(fasta=ref_fasta, quantity="10x", quiet=true)
            reads_fastq = joinpath(dir, "hifiasm_reads.fq")
            run(pipeline(`gunzip -c $(reads_gz)`, reads_fastq))

            outdir = joinpath(dir, "hifiasm_assembly")
            isdir(outdir) && rm(outdir, recursive=true, force=true)
            try
                result = Mycelia.run_hifiasm(
                    fastq=reads_fastq,
                    outdir=outdir,
                    bloom_filter=0,
                    threads=threads
                )
                Test.@test result.outdir == outdir
                Test.@test result.hifiasm_outprefix == joinpath(outdir, basename(reads_fastq) * ".hifiasm")
            finally
                rm(outdir, recursive=true, force=true)
            end
        end
    end
end
