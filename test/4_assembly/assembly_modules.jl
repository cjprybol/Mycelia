# From the Mycelia base directory, run the tests with:
# 
# ```bash
# julia --project=. -e 'include("test/4_assembly/assembly_modules.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/assembly_modules.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise
import Test
import Mycelia
Test.@testset "assembly modules" begin
    Test.@testset "1. Pre‐processing & Read QC" begin
    end
    Test.@testset "2. k‑mer Analysis" begin
    end
    Test.@testset "3. Short Read Assembly" begin
        mktempdir() do dir
            fastq1 = joinpath(dir, "reads_1.fq")
            fastq2 = joinpath(dir, "reads_2.fq")
            open(fastq1, "w") do io
                println(io, "@r1/1")
                println(io, "ACGTACGTACGTACGTACGT")
                println(io, "+")
                println(io, "IIIIIIIIIIIIIIIIIIII")
            end
            open(fastq2, "w") do io
                println(io, "@r1/2")
                println(io, "TGCATGCATGCATGCATGCA")
                println(io, "+")
                println(io, "IIIIIIIIIIIIIIIIIIII")
            end
            
            # Test MEGAHIT
            outdir = joinpath(dir, "megahit")
            mkpath(outdir)
            touch(joinpath(outdir, "final.contigs.fa"))
            result = Mycelia.run_megahit(fastq1=fastq1, fastq2=fastq2, outdir=outdir)
            Test.@test result.outdir == outdir
            Test.@test result.contigs == joinpath(outdir, "final.contigs.fa")
            
            # Test metaSPAdes
            outdir = joinpath(dir, "metaspades")
            mkpath(outdir)
            touch(joinpath(outdir, "contigs.fasta"))
            touch(joinpath(outdir, "scaffolds.fasta"))
            result = Mycelia.run_metaspades(fastq1=fastq1, fastq2=fastq2, outdir=outdir)
            Test.@test result.outdir == outdir
            Test.@test result.contigs == joinpath(outdir, "contigs.fasta")
            Test.@test result.scaffolds == joinpath(outdir, "scaffolds.fasta")
        end
    end
    
    Test.@testset "4. Long Read Assembly" begin
        mktempdir() do dir
            fastq = joinpath(dir, "long_reads.fq")
            open(fastq, "w") do io
                println(io, "@r1")
                println(io, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
                println(io, "+")
                println(io, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII")
            end
            
            # Test Flye
            outdir = joinpath(dir, "flye")
            mkpath(outdir)
            touch(joinpath(outdir, "assembly.fasta"))
            result = Mycelia.run_flye(fastq=fastq, outdir=outdir, genome_size="5m")
            Test.@test result.outdir == outdir
            Test.@test result.assembly == joinpath(outdir, "assembly.fasta")
            
            # Test Canu
            outdir = joinpath(dir, "canu")
            mkpath(outdir)
            touch(joinpath(outdir, "reads.contigs.fasta"))
            result = Mycelia.run_canu(fastq=fastq, outdir=outdir, genome_size="5m")
            Test.@test result.outdir == outdir
            Test.@test result.assembly == joinpath(outdir, "reads.contigs.fasta")
            
            # Test hifiasm
            outdir = joinpath(dir, "hifiasm")
            mkpath(outdir)
            prefix = joinpath(outdir, basename(fastq) * ".hifiasm")
            touch(prefix * ".dummy")
            result = Mycelia.run_hifiasm(fastq=fastq, outdir=outdir)
            Test.@test result.outdir == outdir
            Test.@test result.hifiasm_outprefix == prefix
        end
    end
    
    Test.@testset "5. Hybrid Assembly" begin
        mktempdir() do dir
            short1 = joinpath(dir, "short_1.fq")
            short2 = joinpath(dir, "short_2.fq")
            long_reads = joinpath(dir, "long_reads.fq")
            
            open(short1, "w") do io
                println(io, "@r1/1")
                println(io, "ACGTACGTACGTACGTACGT")
                println(io, "+")
                println(io, "IIIIIIIIIIIIIIIIIIII")
            end
            open(short2, "w") do io
                println(io, "@r1/2")
                println(io, "TGCATGCATGCATGCATGCA")
                println(io, "+")
                println(io, "IIIIIIIIIIIIIIIIIIII")
            end
            open(long_reads, "w") do io
                println(io, "@r1")
                println(io, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
                println(io, "+")
                println(io, "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII")
            end
            
            # Test Unicycler
            outdir = joinpath(dir, "unicycler")
            mkpath(outdir)
            touch(joinpath(outdir, "assembly.fasta"))
            result = Mycelia.run_unicycler(short_1=short1, short_2=short2, long_reads=long_reads, outdir=outdir)
            Test.@test result.outdir == outdir
            Test.@test result.assembly == joinpath(outdir, "assembly.fasta")
        end
    end
    Test.@testset "6. Probabilistic Assembly (Mycelia)" begin
        mktempdir() do dir
            fastq = joinpath(dir, "reads.fq")
            open(fastq, "w") do io
                println(io, "@r1")
                println(io, "ACGTACGTACGTACGTACGT")
                println(io, "+")
                println(io, "IIIIIIIIIIIIIIIIIIII")
                println(io, "@r2")
                println(io, "CGTACGTACGTACGTACGTA")
                println(io, "+")
                println(io, "IIIIIIIIIIIIIIIIIIII")
            end
            
            # Test string graph building
            graph = Mycelia.string_to_ngram_graph(s="ACGTACGTACGTACGTACGT", n=5)
            Test.@test Graphs.nv(graph) > 0
            
            # Test Viterbi error correction functions exist
            Test.@test hasmethod(Mycelia.viterbi_maximum_likelihood_traversals, (Any,))
            Test.@test hasmethod(Mycelia.polish_fastq, (Any,))
        end
    end
    
    Test.@testset "7. Assembly merging" begin
    end
    Test.@testset "8. Polishing & Error Correction" begin
    end
    Test.@testset "6. Strain resolution" begin
    end
    Test.@testset "7. Validation & Quality Control" begin
    end
end
