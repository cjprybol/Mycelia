# OPT4 PARALLEL FREEZE-STREAK BYTE-IDENTITY (td-jbjd, review fix C1)
# =======================================================================
#
# corrector_opt4_freeze_state_test.jl and corrector_opt4_frozen_read_skip_test.jl
# exercise the freeze-streak bookkeeping only under SERIAL decode (no test in
# either file passes `enable_parallel=true`), so the parallel collect-results
# freeze-streak update sites in `_improve_read_set_likelihood_impl` (the
# `Threads.@threads` batch loop's per-read streak/diagnostic update, guarded on
# `skip_flags[i]` / `_frozen_read_at`) were UNTESTED under real concurrency —
# review finding C1. CI runs single-threaded by default, so `@threads` degrades
# to sequential execution and proves nothing; this test relaunches itself under
# `julia --threads=4` to get genuine concurrent batch decode (mirrors
# parallel_soft_em_byte_identity_test.jl's self-hoist harness, opt1).
#
# WHAT IS ASSERTED
#   * The parallel path is ACTUALLY taken (parallel_decode_batches > 0) and
#     freezing ACTUALLY engages (frozen_reads_skipped > 0) in both arms — two
#     positive controls, so a silent revert to serial or a freeze mechanism
#     that never fires cannot pass vacuously.
#   * skip_frozen_reads=true + enable_parallel=true produces BYTE-IDENTICAL
#     corrected FASTQ output (identifier, sequence, quality) to
#     skip_frozen_reads=true + enable_parallel=false (serial), on a fixture
#     with batch_size < n_reads (multiple batches per pass) and enough passes
#     (multi-k, multi-iteration, aggressive threshold) that freezing reliably
#     engages within the first rung.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_opt4_parallel_freeze_byte_identity_test.jl")'

import Test
import Mycelia
import FASTX
import Random

if Threads.nthreads() == 1 && get(ENV, "MYCELIA_OPT4_PFB_CHILD", "0") != "1"
    # Self-hoist to real threads, same pattern as opt1's
    # parallel_soft_em_byte_identity_test.jl.
    proj = Base.active_project()
    thisfile = @__FILE__
    cmd = `$(Base.julia_cmd()) --project=$(proj) --threads=4 -e "include(\"$(thisfile)\")"`
    Test.@testset "opt4 parallel freeze-streak byte-identity (hoisted to -t4)" begin
        Test.@test success(pipeline(
            setenv(cmd, "MYCELIA_OPT4_PFB_CHILD" => "1"; dir = pwd());
            stdout = stdout, stderr = stderr))
    end
else
    const _OPT4PFB_BASES = ['A', 'C', 'G', 'T']

    function _opt4pfb_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
        reflen = length(ref)
        records = FASTX.FASTQ.Record[]
        for i in 1:n_reads
            s = rand(rng, 1:(reflen - readlen + 1))
            seq = collect(ref[s:(s + readlen - 1)])
            if i <= n_err
                p = rand(rng, 1:readlen)
                seq[p] = rand(rng, filter(!=(seq[p]), _OPT4PFB_BASES))
            end
            push!(records,
                FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
        end
        return records
    end

    # Full-record fingerprint (identifier, sequence, quality) — mirrors
    # corrector_opt4_frozen_read_skip_test.jl's `_opt4_recs`.
    function _opt4pfb_recs(records)
        [(String(FASTX.identifier(r)), FASTX.sequence(String, r),
             String(FASTX.quality(r))) for r in records]
    end

    function _opt4pfb_final_recs(result)
        open(FASTX.FASTQ.Reader, result[:metadata][:final_fastq_file]) do rd
            _opt4pfb_recs(collect(rd))
        end
    end

    Test.@testset "opt4 parallel freeze-streak byte-identity (td-jbjd, C1)" begin
        Test.@test Threads.nthreads() > 1     # guard: real threads

        rng = Random.MersenneTwister(606)
        ref = join(rand(rng, _OPT4PFB_BASES, 900))
        reads = _opt4pfb_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
        tmp = mktempdir()
        fq = joinpath(tmp, "in.fastq")
        Mycelia.write_fastq(records = reads, filename = fq)

        # batch_size=50 < n_reads=150 => 3 batches/pass, exercising the
        # parallel collect-results loop across multiple `Threads.@threads`
        # dispatches. freeze_streak_threshold=1 is maximally aggressive (any
        # single no-improvement pass makes a read frozen-eligible on the very
        # next pass), and max_iterations_per_k=3 gives each of the 3 k-rungs
        # enough passes for freezing to engage reliably — the same
        # config shown to produce frozen_reads_skipped > 0 in
        # corrector_opt4_freeze_state_test.jl's "scope" testset.
        run = parallel -> Mycelia.mycelia_iterative_assemble(fq;
            max_k = 17, n_k_rungs = 3, max_iterations_per_k = 3,
            graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            hard_window = true, soft_em = false, batch_size = 50,
            verbose = false, enable_checkpointing = false,
            enable_parallel = parallel,
            skip_frozen_reads = true, freeze_streak_threshold = 1,
            output_dir = joinpath(tmp, parallel ? "parallel" : "serial"))

        r_parallel = Base.redirect_stderr(devnull) do
            run(true)
        end
        r_serial = Base.redirect_stderr(devnull) do
            run(false)
        end

        parallel_errors = r_parallel[:metadata][:corrector_errors]
        serial_errors = r_serial[:metadata][:corrector_errors]

        # Positive control (a): the parallel path was ACTUALLY taken. A silent
        # revert-to-serial regression would leave this at 0 even though
        # byte-identity still (trivially) holds.
        Test.@test parallel_errors[:parallel_decode_batches] > 0
        Test.@test serial_errors[:parallel_decode_batches] == 0

        # Positive control (b): freezing ACTUALLY engaged in BOTH arms. A
        # config that never crosses the freeze threshold would make the
        # byte-identity comparison below vacuous (both arms trivially equal
        # because neither ever exercises the freeze-skip code path).
        Test.@test parallel_errors[:frozen_reads_skipped] > 0
        Test.@test serial_errors[:frozen_reads_skipped] > 0

        # The actual lock: serial and parallel decode must agree exactly, read
        # for read, byte for byte, under freezing.
        Test.@test _opt4pfb_final_recs(r_parallel) == _opt4pfb_final_recs(r_serial)
    end
end
