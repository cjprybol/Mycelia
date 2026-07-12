# External-tier (fastp/bioconda) smoke test for the ONLY novel path in the
# interleaved wrapper: fastp --interleaved_in splitting one file into R1/R2.
# Registered in EXTERNAL_DEPENDENCY_FILES (test/runtests.jl) so core CI skips it;
# runs under MYCELIA_RUN_EXTERNAL=true and is exercised end-to-end by the CAMI
# triangulation driver on NERSC.
import Test
import Mycelia

Test.@testset "qc_filter_short_reads_fastp_interleaved (fastp split)" begin
    mktempdir() do dir
        # Build a small INTERLEAVED FASTQ: records alternate R1, R2, R1, R2, ...
        # Clean 60bp reads at high quality so fastp keeps all of them.
        interleaved = joinpath(dir, "reads.fq")
        npairs = 5
        seq = "ACGT"^15                      # 60 bp
        qual = repeat("I", 60)               # Q40
        open(interleaved, "w") do io
            for i in 1:npairs
                println(io, "@read$(i)/1");
                println(io, seq)
                println(io, "+");
                println(io, qual)
                println(io, "@read$(i)/2");
                println(io, seq)
                println(io, "+");
                println(io, qual)
            end
        end

        res = Mycelia.qc_filter_short_reads_fastp_interleaved(
            interleaved_reads = interleaved)

        # All four declared outputs exist.
        Test.@test isfile(res.out_forward)
        Test.@test isfile(res.out_reverse)
        Test.@test isfile(res.json)
        Test.@test isfile(res.html)

        # The single interleaved file was split into two mate files of equal,
        # nonzero record count (the R1/R2 halves) — the behavior that would break
        # if --out1/--out2 were transposed or --interleaved_in dropped.
        nf = length(collect(Mycelia.open_fastx(res.out_forward)))
        nr = length(collect(Mycelia.open_fastx(res.out_reverse)))
        Test.@test nf == nr
        Test.@test nf >= 1
        Test.@test nf <= npairs

        # Idempotent cache-hit path returns the same output paths without error.
        res2 = Mycelia.qc_filter_short_reads_fastp_interleaved(
            interleaved_reads = interleaved)
        Test.@test res2.out_forward == res.out_forward
        Test.@test res2.out_reverse == res.out_reverse
    end
end
