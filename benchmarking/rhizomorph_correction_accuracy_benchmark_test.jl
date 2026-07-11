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
