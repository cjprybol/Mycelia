# Unit test for the Rhizomorph correction ACCURACY benchmark's SIMULATOR.
#
# Validates the truth-retention + length-preservation invariants of
# `simulate_substitution_reads` (which needs Mycelia for
# mutate_dna_substitution_fraction). The pure per-base metric math is tested
# separately in rhizomorph_correction_accuracy_metrics_test.jl (Mycelia-free,
# millisecond runtime). The corrector integration (correct_reads_scalable,
# run_accuracy_benchmark) is exercised by the MYCELIA_RCA_SMOKE end-to-end run,
# whose scale + scored-fraction guard is its safety net; note the FASTQ read-ID
# round-trip (write read_$i -> corrector -> read back) is ONLY observable there,
# never in these unit tests, so the smoke run is a required verification step.
#
# Run:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/rhizomorph_correction_accuracy_benchmark_test.jl

import Test
import Random
import BioSequences
import CSV
import DataFrames
import Statistics

include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_benchmark.jl"))

Test.@testset "simulate_substitution_reads — truth retention + length preservation" begin
    rng = Random.MersenneTwister(7)
    refseq = BioSequences.randdnaseq(rng, 500)
    err = 0.05
    readlen = 100
    records, truth_by_id,
    observed_by_id,
    injected_total,
    sampled_bases = simulate_substitution_reads(
        refseq, readlen, 5.0, err, rng; assigned_q = 20)

    Test.@test !isempty(records)
    Test.@test length(truth_by_id) == length(records)

    # Every read: truth and observed same length (substitution-only), and observed
    # differs from truth at exactly the injected positions.
    per_read_injected = 0
    for rec in records
        rid = FASTX.identifier(rec)
        T = truth_by_id[rid]
        O = observed_by_id[rid]
        Test.@test length(T) == length(O)
        Test.@test length(O) == readlen
        Tc = collect(T)
        Oc = collect(O)
        diffs = count(p -> Oc[p] != Tc[p], eachindex(Tc))
        per_read_injected += diffs
        # ceil(err * readlen) substitutions injected per read (distinct positions,
        # always to a different base) => exactly that many diffs.
        Test.@test diffs == ceil(Int, err * readlen)
    end
    Test.@test injected_total == per_read_injected
    Test.@test sampled_bases == length(records) * readlen
end

Test.@testset "simulate_substitution_reads — readlen > genome clamps to genome length" begin
    # The effective_readlen = min(readlen, glen) clamp: a 60 bp genome with a
    # 150 bp requested readlen must emit 60 bp reads, not error or pad.
    rng = Random.MersenneTwister(11)
    glen = 60
    refseq = BioSequences.randdnaseq(rng, glen)
    records, truth_by_id,
    observed_by_id,
    injected_total,
    sampled_bases = simulate_substitution_reads(
        refseq, 150, 5.0, 0.05, rng; assigned_q = 20)
    Test.@test !isempty(records)
    for rec in records
        Test.@test length(FASTX.sequence(rec)) == glen   # clamped to genome length
    end
    Test.@test sampled_bases == length(records) * glen
end

Test.@testset "scramble_reads — Control B per-read permutation invariants" begin
    # scramble_reads must (a) preserve each read's length, (b) preserve the
    # per-read injected-error COUNT (a permutation applied identically to observed
    # and truth is a bijection, so mismatches move but do not appear/disappear —
    # this is what keeps real and null recall comparable), and (c) actually
    # scramble base ORDER (so shared k-mer coverage is destroyed). Pure + fast.
    rng = Random.MersenneTwister(3)
    refseq = BioSequences.randdnaseq(rng, 400)
    err = 0.05
    readlen = 100
    records, truth_by_id,
    observed_by_id, _injected_total,
    _sampled_bases = simulate_substitution_reads(
        refseq, readlen, 5.0, err, rng; assigned_q = 20)

    scr_records, scr_truth_by_id,
    scr_observed_by_id = scramble_reads(records, truth_by_id, observed_by_id, 20)

    Test.@test length(scr_records) == length(records)
    Test.@test Set(keys(scr_truth_by_id)) == Set(keys(truth_by_id))
    order_changed = false
    for rec in records
        rid = FASTX.identifier(rec)
        T = truth_by_id[rid]
        O = observed_by_id[rid]
        sT = scr_truth_by_id[rid]
        sO = scr_observed_by_id[rid]
        # (a) length preserved
        Test.@test length(sT) == length(T)
        Test.@test length(sO) == length(O)
        # base multiset preserved (permutation only)
        Test.@test sort(collect(sT)) == sort(collect(T))
        Test.@test sort(collect(sO)) == sort(collect(O))
        # (b) injected-error count preserved
        real_diffs = count(p -> collect(O)[p] != collect(T)[p], eachindex(collect(T)))
        scr_diffs = count(p -> collect(sO)[p] != collect(sT)[p], eachindex(collect(sT)))
        Test.@test scr_diffs == real_diffs
        # (c) order genuinely changed for at least some read (100 bases => ~certain)
        order_changed |= (sT != T)
    end
    Test.@test order_changed

    # Determinism: same seed => identical scramble.
    scr_records2, scr_truth2, _ = scramble_reads(records, truth_by_id, observed_by_id, 20)
    Test.@test scr_truth2 == scr_truth_by_id
    Test.@test [FASTX.sequence(String, r) for r in scr_records2] ==
               [FASTX.sequence(String, r) for r in scr_records]
end

Test.@testset "run_accuracy_benchmark — Control A (err=0) + Control B (scramble null)" begin
    # End-to-end smoke: runs the WIRED corrector on a small synthetic genome
    # across an err=0.0 (Control A) cell and one real-error cell, plus the
    # read-scramble null (Control B). Asserts robust invariants only — no
    # hard-coded collapse magnitude — and reports the observed numbers via the
    # benchmark's own stdout. Heavy (invokes mycelia_iterative_assemble) but
    # bounded by the small ENV config below.
    # err=0.0 is Control A; err=0.01 is the real-recovery cell (low enough error
    # that the real corrector genuinely fixes most injected errors, giving a wide
    # real>>null gap for the C1 positive control below).
    withenv(
        "MYCELIA_RCA_SMOKE" => "true",
        "MYCELIA_RCA_SMOKE_GENOME_LEN" => "1500",
        "MYCELIA_RCA_COVERAGE" => "8",
        "MYCELIA_RCA_ERR" => "0.0,0.01",
        "MYCELIA_RCA_READLEN" => "120",
        "MYCELIA_RCA_K" => "21",
        "MYCELIA_RCA_SEED" => "42"
    ) do
        result = run_accuracy_benchmark()
        rows = CSV.read(result.csv_path, DataFrames.DataFrame)

        # --- Control A: the err=0.0 over-correction cell exists and is bounded ---
        control_a = rows[rows.error_rate .== 0.0, :]
        Test.@test DataFrames.nrow(control_a) >= 1
        Test.@test all(control_a.injected_errors .== 0)
        Test.@test all(control_a.reads_scored .> 0)
        for v in control_a.over_correction_rate
            Test.@test !isnan(v)
            Test.@test 0.0 <= v <= 1.0
        end
        # recall/precision are NaN at err=0.0 (injected==0) — must not crash the
        # reader, and that NaN is the correct "undefined here" signal.
        Test.@test all(isnan, control_a.recall)

        # --- Control B (C1 positive control): real recall must EXCEED null by a
        # strict margin. A weaker `null <= real` check would pass a silent real-
        # corrector regression (runs but fixes nothing => real recall 0.0, which
        # is finite and still >= null 0.0), collapsing the "real >> null" claim to
        # "real == null == 0" under green CI. Assert a STRICT gap safely below the
        # observed real recall (~1.0 at err=0.01) so a regression to ~0 fails loud.
        Test.@test isfinite(result.real_pooled_recall)   # err=0.01 cell has injected errors
        Test.@test isfinite(result.null_pooled_recall)
        Test.@test result.real_pooled_recall >= result.null_pooled_recall + 0.05
        # Both arms must have RUN (no caught-but-recorded corrector exception): a
        # crashed cell recorded as ok=false must not be scored as a legitimate
        # recall — otherwise a real-arm crash (recall NaN/0) could satisfy the gap
        # against an also-crashed null arm.
        Test.@test all(rows.ok)
        Test.@test all(rows.null_ok)
        # Null columns present, seeded, and in-range where defined.
        Test.@test all(rows.null_seed .== 20260711)
        for v in rows.null_recall
            Test.@test isnan(v) || (0.0 <= v <= 1.0)
        end
    end
end
