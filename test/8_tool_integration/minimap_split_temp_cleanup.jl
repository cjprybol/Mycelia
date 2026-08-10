# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/8_tool_integration/minimap_split_temp_cleanup.jl")'
# ```
#
# minimap2 writes one `PREFIX.NNNN.tmp` chunk per index part when run against a
# multi-part index, and unlinks them itself on clean exit. They survive only when
# the process is killed first (walltime, OOM, scancel) -- and nothing can consume
# them afterwards, since minimap2 cannot resume from them. A 2026-07 CAMI_I_LOW
# run aborted this way and stranded 783 GiB on NERSC scratch.
#
# These tests need no external tools: they exercise the prefix derivation and the
# reclaim helper against synthetic files, and check that the command builders
# hand the caller a prefix to clean up with (run_mapping is never invoked).

import Test
import Mycelia

Test.@testset "minimap2 split-index temp cleanup" begin
    Test.@testset "split prefix derives from the output path" begin
        Test.@test Mycelia.minimap_split_prefix("/x/y.sorted.bam") ==
                   "/x/y.sorted.bam.tmp"
    end

    Test.@testset "reclaims chunks, spares siblings" begin
        mktempdir() do dir
            outfile = joinpath(dir, "sample.ref.mmi.minimap2.sorted.bam")
            split_prefix = Mycelia.minimap_split_prefix(outfile)

            # Three minimap2 split chunks, 1 KiB each.
            chunks = [string(split_prefix, ".", lpad(i, 4, '0'), ".tmp") for i in 0:2]
            for c in chunks
                write(c, rand(UInt8, 1024))
            end

            # Siblings that share a stem and MUST survive: the real output BAM,
            # its index, and samtools' own sort temps (a different infix).
            survivors = [
                outfile,
                outfile * ".bai",
                string(outfile, ".sort.tmp.0000.bam"),
                joinpath(dir, "unrelated.sorted.bam")
            ]
            for s in survivors
                write(s, "keep")
            end

            result = Mycelia.cleanup_minimap_split_temps(split_prefix; verbose = false)

            Test.@test result.removed == 3
            Test.@test result.bytes == 3 * 1024
            for c in chunks
                Test.@test !isfile(c)
            end
            for s in survivors
                Test.@test isfile(s)
            end
        end
    end

    Test.@testset "is idempotent and safe on a missing directory" begin
        mktempdir() do dir
            prefix = joinpath(dir, "nothing-here.bam.tmp")
            r = Mycelia.cleanup_minimap_split_temps(prefix; verbose = false)
            Test.@test r.removed == 0
            Test.@test r.bytes == 0
        end
        r = Mycelia.cleanup_minimap_split_temps(
            joinpath("/nonexistent-dir-for-mycelia-test", "x.bam.tmp"); verbose = false)
        Test.@test r.removed == 0
    end

    Test.@testset "builders return a split_prefix the caller can clean up" begin
        mktempdir() do dir
            ref = joinpath(dir, "ref.fa")
            write(ref, ">ref\n" * "ACGT"^100 * "\n")
            fq = joinpath(dir, "reads.fq")
            write(fq, "@r1\nACGT\n+\nIIII\n")
            outfile = joinpath(dir, "out.sorted.bam")

            res = Mycelia.minimap_map(;
                fasta = ref, fastq = fq, mapping_type = "sr",
                outfile = outfile, as_string = true)

            Test.@test haskey(res, :split_prefix)
            Test.@test res.split_prefix == Mycelia.minimap_split_prefix(res.outfile)
            # The prefix handed back must be the one minimap2 is actually told
            # to use -- otherwise cleanup would target the wrong files.
            Test.@test occursin("--split-prefix=$(res.split_prefix)", res.cmd)
        end
    end
end
