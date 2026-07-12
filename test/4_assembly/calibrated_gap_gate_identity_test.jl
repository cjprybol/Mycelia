# Calibrated decode gate (td-4osf Stage 3): byte-identity when OFF, and correct
# gating direction when ON.
#
# The gate reverts a decoded base to the OBSERVED base wherever the per-position
# Viterbi gap is below `calibrated_gap_threshold`. Two guarantees:
#   1. OFF (threshold === nothing) and the permissive limit (threshold = -Inf,
#      which engages gap recording + the gate but reverts nothing) produce
#      BYTE-IDENTICAL output — the telemetry + no-op gate never perturb the decode.
#   2. A large threshold reverts every finite-gap correction back toward the
#      observed read, so the gated result can only move AWAY from the corrector's
#      output and TOWARD the raw read (the gate suppresses, never fabricates).
#
# Run:
#   julia --project=. test/4_assembly/calibrated_gap_gate_identity_test.jl

import Test
import Mycelia
import FASTX
import Random

rec(id, seq) = FASTX.FASTQ.Record(id, seq, String(fill('I', length(seq))))
seqs(reads) = [FASTX.sequence(String, r) for r in reads]
mismatches(a, b) = count(x -> x[1] != x[2], zip(collect(a), collect(b)))

Test.@testset "calibrated gap gate: OFF byte-identity + ON reverts to observed" begin
    k = 7
    clean_seq = "ATGCGTACGTACGTTAGCCGATACAGGTCA"
    reads = [rec("clean$i", clean_seq) for i in 1:10]
    err_seq = collect(clean_seq)
    err_seq[15] = err_seq[15] == 'A' ? 'C' : 'A'          # single substitution
    err_seq = String(err_seq)
    push!(reads, rec("err1", err_seq))
    graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k;
        dataset_id = "gate", mode = :canonical, memory_profile = :full)

    run_at(thresh) = begin
        Random.seed!(11)   # deterministic: pin any global-RNG draws
        out,
        _ = Mycelia.improve_read_set_likelihood(reads, graph, k;
            graph_mode = :canonical, calibrated_gap_threshold = thresh)
        out
    end

    off = run_at(nothing)
    permissive = run_at(-Inf)
    strict = run_at(100.0)

    # (1) BYTE-IDENTITY (the core guarantee): gate OFF === gate engaged-but-permissive.
    # Proves the gap-recording telemetry + the no-op gate never perturb the decode.
    Test.@test seqs(off) == seqs(permissive)
    Test.@test length(off) == length(reads)

    # (2) MONOTONE REVERT-TOWARD-OBSERVED: gating can only move a corrected base back
    # to the observed base, so under any threshold each read is at least as close to
    # its OWN observed input as the ungated correction is (the gate suppresses, never
    # fabricates). Guaranteed regardless of whether a correction actually fired.
    # (Behavioral proof that the gate changes real corrections lives in the
    # full-pipeline frontier smoke, where err=0.10 reliably triggers corrections.)
    obs = seqs(reads)
    for i in eachindex(reads)
        Test.@test mismatches(seqs(strict)[i], obs[i]) <= mismatches(seqs(off)[i], obs[i])
    end

    # Solid clean reads are untouched under every setting (no spurious edits).
    for out in (off, permissive, strict)
        for i in 1:10
            Test.@test seqs(out)[i] == clean_seq
        end
    end
end
