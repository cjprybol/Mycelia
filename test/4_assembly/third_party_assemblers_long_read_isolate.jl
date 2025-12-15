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

