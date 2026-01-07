# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/third_party_assemblers_long_read_metagenomic.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/third_party_assemblers_long_read_metagenomic.jl", "test/4_assembly", execute=false)'
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

Test.@testset "Long Read Metagenomic Assembly" begin
    threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

    mktempdir() do dir
        ref_fasta = joinpath(dir, "metaflye_specific_ref.fasta")
        rng1 = StableRNGs.StableRNG(800)
        rng2 = StableRNGs.StableRNG(801)
        genome1 = BioSequences.randdnaseq(rng1, 4_000)
        genome2 = BioSequences.randdnaseq(rng2, 3_000)
        Mycelia.write_fasta(
            outfile=ref_fasta,
            records=[
                FASTX.FASTA.Record("meta_flye_genome_1", genome1),
                FASTX.FASTA.Record("meta_flye_genome_2", genome2),
            ]
        )

        reads_gz = Mycelia.simulate_pacbio_reads(fasta=ref_fasta, quantity="10x", quiet=true)
        reads_fastq = joinpath(dir, "meta_flye_reads.fq")
        run(pipeline(`gunzip -c $(reads_gz)`, reads_fastq))

        outdir = joinpath(dir, "metaflye_assembly")
        isdir(outdir) && rm(outdir, recursive=true, force=true)
        try
            result = Mycelia.run_metaflye(
                fastq=reads_fastq,
                outdir=outdir,
                genome_size="7k",
                read_type="pacbio-hifi",
                threads=threads
            )
            Test.@test result.outdir == outdir
            Test.@test isfile(result.assembly)
        finally
            rm(outdir, recursive=true, force=true)
        end
    end
end

