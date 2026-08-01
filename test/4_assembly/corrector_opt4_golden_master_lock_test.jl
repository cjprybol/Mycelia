# OPT4 FROZEN-READ SKIP — GOLDEN LOCK vs PRE-OPT4 MASTER (td-jbjd, pass 3)
# =============================================================================
#
# Pass 1's disabled-path lock (corrector_opt4_frozen_read_skip_test.jl) proves
# `skip_frozen_reads=false` (explicit) is byte-identical to the kwarg being
# OMITTED — both on THIS branch. That is a real but strictly weaker claim than
# "opt4 is a true no-op when off": a bug shared by both call sites on this
# branch (e.g. a stray side effect introduced by the freeze-streak plumbing
# that fires regardless of the flag) would not be caught by comparing the
# branch to itself.
#
# This test closes that gap by comparing THIS branch's disabled-path output
# directly against pre-opt4 MASTER (merge b032fab3, the commit
# cp/opt4-frozen-read-design forked from) — a completely independent
# checkout that has never seen the opt4 diff at all. If this test ever
# diverges, that is a CRITICAL finding: it means opt4-off is not a true no-op,
# contradicting the pass 1/2 acceptance criterion ("skip_frozen_reads=false
# (default) is a bit-identical no-op vs pre-opt4 master (locked by test)").
#
# PROVENANCE OF THE HARDCODED HASH
#   Generated 2026-07-27 by running the fixture below (verbatim) via:
#     LD_LIBRARY_PATH='' julia --project=<checkout> \
#       /path/to/opt4_golden_fixture.jl [--opt4-off]
#   against THREE checkouts, all producing the IDENTICAL SHA256 below:
#     1. pre-opt4 master, git commit b032fab3 (no opt4 kwargs — they don't
#        exist on that commit; ran with no opt4 flags at all)
#     2. cp/opt4-frozen-read-design (commit 62f0b3d1, this branch's pass-2
#        HEAD at write time) with `skip_frozen_reads=false` explicit
#     3. cp/opt4-frozen-read-design (same commit) with all opt4 kwargs OMITTED
#   Fixture params: MersenneTwister(20260726) seed, 900bp random reference,
#   150 reads/80bp/25 injected single-base substitution errors, max_k=17,
#   n_k_rungs=3, max_iterations_per_k=3, graph_mode=:canonical,
#   skip_solid=true, cheap_correct=true, hard_window=true, soft_em=false,
#   batch_size=50, enable_checkpointing=false.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_opt4_golden_master_lock_test.jl")'

import Test
import Mycelia
import FASTX
import Random
import SHA

# Hash captured from pre-opt4 master (b032fab3) — see provenance note above.
# DO NOT update this constant to make a failing test pass; a mismatch means
# the disabled path regressed and must be investigated, not silenced.
const _OPT4_GOLDEN_MASTER_SHA256 = "e0c249c092399d98afb40137ef4bb07ae3d9c1d1c94b3aafc01d95a6516b6c54"

const _OPT4GM_BASES = ['A', 'C', 'G', 'T']

# Verbatim copy of the fixture generator used to produce the golden hash above
# (kept identical on purpose — do not "clean up" the parameters without
# regenerating the golden hash against master again).
function _opt4gm_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        if i <= n_err
            p = rand(rng, 1:readlen)
            seq[p] = rand(rng, filter(!=(seq[p]), _OPT4GM_BASES))
        end
        push!(records,
            FASTX.FASTQ.Record("r$i", String(seq), String(fill('I', readlen))))
    end
    return records
end

function _opt4gm_run(fq, out_dir; skip_frozen_reads_explicit::Bool)
    common = (; max_k = 17, n_k_rungs = 3, max_iterations_per_k = 3,
        graph_mode = :canonical, skip_solid = true, cheap_correct = true,
        hard_window = true, soft_em = false, batch_size = 50,
        verbose = false, enable_checkpointing = false, output_dir = out_dir)
    if skip_frozen_reads_explicit
        return Mycelia.mycelia_iterative_assemble(fq; common...,
            skip_frozen_reads = false, freeze_streak_threshold = 2,
            freeze_across_rungs = false)
    end
    return Mycelia.mycelia_iterative_assemble(fq; common...)
end

Test.@testset "opt4 golden lock vs pre-opt4 master (td-jbjd pass 3)" begin
    rng = Random.MersenneTwister(20260726)
    ref = join(rand(rng, _OPT4GM_BASES, 900))
    reads = _opt4gm_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 25)
    tmp = mktempdir()
    fq = joinpath(tmp, "in.fastq")
    Mycelia.write_fastq(records = reads, filename = fq)

    Test.@testset "skip_frozen_reads=false (explicit) matches golden master hash" begin
        result = _opt4gm_run(
            fq, joinpath(tmp, "explicit"); skip_frozen_reads_explicit = true)
        digest = bytes2hex(SHA.sha256(read(result[:metadata][:final_fastq_file])))
        Test.@test digest == _OPT4_GOLDEN_MASTER_SHA256
    end

    Test.@testset "opt4 kwargs omitted matches golden master hash" begin
        result = _opt4gm_run(
            fq, joinpath(tmp, "omitted"); skip_frozen_reads_explicit = false)
        digest = bytes2hex(SHA.sha256(read(result[:metadata][:final_fastq_file])))
        Test.@test digest == _OPT4_GOLDEN_MASTER_SHA256
    end
end
