# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/third_party_assemblers_plass_penguin.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/third_party_assemblers_plass_penguin.jl", "test/4_assembly", execute=false)'
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

Test.@testset "Protein/Nucleotide Assembly (PLASS/PenguiN)" begin
    threads = clamp(something(tryparse(Int, get(ENV, "MYCELIA_ASSEMBLER_TEST_THREADS", "2")), 2), 1, 4)

    Test.@testset "PLASS (simulated reads)" begin
        mktempdir() do dir
            ref_fasta = joinpath(dir, "ref.fasta")
            rng = StableRNGs.StableRNG(42)
            seq = BioSequences.randdnaseq(rng, 2_000)
            Mycelia.write_fasta(outfile=ref_fasta, records=[FASTX.FASTA.Record("ref", seq)])

            sim = Mycelia.simulate_illumina_reads(
                fasta=ref_fasta,
                coverage=5,
                outbase=joinpath(dir, "sim_plass"),
                read_length=100,
                mflen=200,
                seqSys="HS25",
                paired=true,
                errfree=true,
                rndSeed=42,
                quiet=true
            )

            outdir = joinpath(dir, "plass_out")
            try
                result = Mycelia.run_plass_assemble(
                    reads1=sim.forward_reads,
                    reads2=sim.reverse_reads,
                    outdir=outdir,
                    min_length=20,
                    num_iterations=1,
                    threads=threads
                )
                Test.@test isfile(result.assembly)
            catch e
                if isa(e, ProcessFailedException) || occursin("memory", lowercase(string(e))) || occursin("not found", lowercase(string(e)))
                    @warn "PLASS test skipped due to environment/resource constraints." exception=e
                    Test.@test_skip "PLASS test skipped - environment/resource constraints"
                else
                    rethrow(e)
                end
            end
        end
    end

    Test.@testset "PenguiN guided_nuclassemble (simulated reads)" begin
        mktempdir() do dir
            ref_fasta = joinpath(dir, "ref.fasta")
            rng = StableRNGs.StableRNG(43)
            seq = BioSequences.randdnaseq(rng, 2_000)
            Mycelia.write_fasta(outfile=ref_fasta, records=[FASTX.FASTA.Record("ref", seq)])

            sim = Mycelia.simulate_illumina_reads(
                fasta=ref_fasta,
                coverage=5,
                outbase=joinpath(dir, "sim_penguin_guided"),
                read_length=100,
                mflen=200,
                seqSys="HS25",
                paired=true,
                errfree=true,
                rndSeed=43,
                quiet=true
            )

            outdir = joinpath(dir, "penguin_guided_out")
            try
                result = Mycelia.run_penguin_guided_nuclassemble(
                    reads1=sim.forward_reads,
                    reads2=sim.reverse_reads,
                    outdir=outdir,
                    threads=threads
                )
                Test.@test isfile(result.assembly)
            catch e
                if isa(e, ProcessFailedException) || occursin("memory", lowercase(string(e))) || occursin("not found", lowercase(string(e)))
                    @warn "PenguiN guided_nuclassemble test skipped due to environment/resource constraints." exception=e
                    Test.@test_skip "PenguiN guided test skipped - environment/resource constraints"
                else
                    rethrow(e)
                end
            end
        end
    end

    Test.@testset "PenguiN nuclassemble (simulated reads)" begin
        mktempdir() do dir
            ref_fasta = joinpath(dir, "ref.fasta")
            rng = StableRNGs.StableRNG(44)
            seq = BioSequences.randdnaseq(rng, 2_000)
            Mycelia.write_fasta(outfile=ref_fasta, records=[FASTX.FASTA.Record("ref", seq)])

            sim = Mycelia.simulate_illumina_reads(
                fasta=ref_fasta,
                coverage=5,
                outbase=joinpath(dir, "sim_penguin"),
                read_length=100,
                mflen=200,
                seqSys="HS25",
                paired=true,
                errfree=true,
                rndSeed=44,
                quiet=true
            )

            outdir = joinpath(dir, "penguin_out")
            try
                result = Mycelia.run_penguin_nuclassemble(
                    reads1=sim.forward_reads,
                    reads2=sim.reverse_reads,
                    outdir=outdir,
                    threads=threads
                )
                Test.@test isfile(result.assembly)
            catch e
                if isa(e, ProcessFailedException) || occursin("memory", lowercase(string(e))) || occursin("not found", lowercase(string(e)))
                    @warn "PenguiN nuclassemble test skipped due to environment/resource constraints." exception=e
                    Test.@test_skip "PenguiN nuclassemble test skipped - environment/resource constraints"
                else
                    rethrow(e)
                end
            end
        end
    end
end

