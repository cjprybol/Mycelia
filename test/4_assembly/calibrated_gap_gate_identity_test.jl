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


# td-ufsu item 1: the same OFF vs -Inf guarantee on the windowed path, which
# applies the gate per sub-read. Residual risk is a sub-read coordinate-frame
# bug; this smoke needs a >1024-observation read with hard windows that decode.
Test.@testset "windowed path: OFF byte-identity + gate engaged" begin
    k = 13
    rng = Random.MersenneTwister(42)
    reflen = 2000
    readlen = 1200
    n_clean = 20
    n_err = 3
    err_locus = 1000
    bases = ['A', 'C', 'G', 'T']
    ref = join(rand(rng, bases, reflen))
    reads = FASTX.FASTQ.Record[]
    for i in 1:n_clean
        s = rand(rng, 1:(reflen - readlen + 1))
        push!(reads, rec("c$i", ref[s:(s + readlen - 1)]))
    end
    err_reads = FASTX.FASTQ.Record[]
    for i in 1:n_err
        s = err_locus - div(readlen, 2)
        off = div(readlen, 2) + 1
        seq = collect(ref[s:(s + readlen - 1)])
        seq[off] = first(filter(!=(seq[off]), bases))
        er = rec("e$i", String(seq))
        push!(err_reads, er)
        push!(reads, er)
    end

    graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k;
        dataset_id = "win-gate", mode = :canonical, memory_profile = :full)
    hard = Mycelia._hard_vertex_set(graph, k)
    Test.@test !isempty(hard)

    r = err_reads[1]
    n_obs = length(FASTX.sequence(r)) - k + 1
    Test.@test n_obs > 1024
    wins = Mycelia._hard_window_ranges(r, k, hard; pad = k, max_window = 500)
    Test.@test !isempty(wins)

    function run_at(thresh)
        diag = Mycelia.CorrectorDiagnostics()
        record, _improved, decoded, _divergent = Mycelia.improve_read_likelihood_windowed_detail(
            r, graph, k, hard;
            graph_mode = :canonical, calibrated_gap_threshold = thresh,
            diagnostics = diag)
        return (record = record, decoded = decoded, diag = diag)
    end

    off = run_at(nothing)
    permissive = run_at(-Inf)
    strict = run_at(100.0)

    # Hard windows actually decoded (otherwise this is a passthrough, not a gate smoke).
    Test.@test off.decoded > 0
    Test.@test permissive.decoded > 0
    Test.@test permissive.decoded == off.decoded

    # BYTE-IDENTITY: gate OFF === gate engaged-but-permissive.
    Test.@test FASTX.sequence(String, off.record) == FASTX.sequence(String, permissive.record)

    # Gate engaged and did not fall open: the contract (path + gaps) held.
    Test.@test permissive.diag.gate_skipped[] == 0

    # Gating can only move a corrected base back toward the observed read.
    obs = FASTX.sequence(String, r)
    Test.@test mismatches(FASTX.sequence(String, strict.record), obs) <=
               mismatches(FASTX.sequence(String, off.record), obs)
end

# td-ufsu item 3: AbstractString overloads must forward calibrated_gap_threshold
# so the string API is not a silent completeness hole off the production path.
Test.@testset "AbstractString overloads forward calibrated_gap_threshold" begin
    k = 7
    clean_seq = "ATGCGTACGTACGTTAGCCGATACAGGTCA"
    reads = [rec("clean$i", clean_seq) for i in 1:10]
    err_seq = collect(clean_seq)
    err_seq[15] = err_seq[15] == 'A' ? 'C' : 'A'
    err_seq = String(err_seq)
    push!(reads, rec("err1", err_seq))
    graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k;
        dataset_id = "str-gate", mode = :canonical, memory_profile = :full)
    qual = String(fill('I', length(err_seq)))

    off_seq, _ = Mycelia.find_optimal_sequence_path(err_seq, qual, graph, k;
        graph_mode = :canonical, calibrated_gap_threshold = nothing)
    perm_seq, _ = Mycelia.find_optimal_sequence_path(err_seq, qual, graph, k;
        graph_mode = :canonical, calibrated_gap_threshold = -Inf)
    Test.@test off_seq == perm_seq

    off_v = Mycelia.try_viterbi_path_improvement(err_seq, qual, graph, k;
        graph_mode = :canonical, calibrated_gap_threshold = nothing)
    perm_v = Mycelia.try_viterbi_path_improvement(err_seq, qual, graph, k;
        graph_mode = :canonical, calibrated_gap_threshold = -Inf)
    Test.@test off_v == perm_v
end
