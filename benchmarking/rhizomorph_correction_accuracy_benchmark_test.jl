# Unit test for the Rhizomorph correction ACCURACY benchmark.
#
# Validates the pure per-base metric classifier `per_base_metrics` on tiny
# hand-constructed cases with a KNOWN correct answer, and the truth-retention +
# length-preservation invariants of `simulate_substitution_reads`. This does not
# run the corrector (that is the slow, network-optional smoke run); it pins the
# metric MATH, which is what makes the benchmark's numbers trustworthy.
#
# Including the benchmark file loads Mycelia (needed for the substitution
# simulator's mutate_dna_substitution_fraction); the metric assertions
# themselves are pure.
#
# Run:
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/rhizomorph_correction_accuracy_benchmark_test.jl

import Test
import Random
import BioSequences

include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_benchmark.jl"))

Test.@testset "rhizomorph correction accuracy — per_base_metrics" begin
    # --- Case A: a clean fix -------------------------------------------------
    # 1 injected substitution (pos 2 C->G), corrected back to truth.
    truth = Dict("r1" => "ACGT")
    observed = Dict("r1" => "AGGT")
    corrected = Dict("r1" => "ACGT")
    m = per_base_metrics(truth, observed, corrected)
    Test.@test m.injected == 1
    Test.@test m.total_edits == 1
    Test.@test m.tp == 1
    Test.@test m.over == 0
    Test.@test m.recall == 1.0
    Test.@test m.precision == 1.0
    Test.@test m.over_rate == 0.0
    Test.@test m.correction_rate == 1 / 4
    Test.@test m.reads_scored == 1

    # --- Case B: an over-correction (a correct base made wrong) ---------------
    # No injected error; corrector changes pos 4 T->A. precision must drop, and
    # recall is NaN because there were zero errors to recall.
    truth = Dict("r1" => "ACGT")
    observed = Dict("r1" => "ACGT")
    corrected = Dict("r1" => "ACGA")
    m = per_base_metrics(truth, observed, corrected)
    Test.@test m.injected == 0
    Test.@test m.total_edits == 1
    Test.@test m.tp == 0
    Test.@test m.over == 1
    Test.@test isnan(m.recall)          # 0 injected -> undefined recall
    Test.@test m.precision == 0.0       # 1 edit, 0 good -> precision 0
    Test.@test m.over_rate == 1 / 4

    # --- Case C: a missed error (no edit at the error) ------------------------
    truth = Dict("r1" => "ACGT")
    observed = Dict("r1" => "AGGT")
    corrected = Dict("r1" => "AGGT")
    m = per_base_metrics(truth, observed, corrected)
    Test.@test m.injected == 1
    Test.@test m.total_edits == 0
    Test.@test m.tp == 0
    Test.@test m.recall == 0.0
    Test.@test isnan(m.precision)       # 0 edits -> undefined precision
    Test.@test m.over_rate == 0.0

    # --- Case D: a dropped read (absent from corrected set) -------------------
    truth = Dict("r1" => "ACGT", "r2" => "TTTT")
    observed = Dict("r1" => "AGGT", "r2" => "TTTT")
    corrected = Dict("r1" => "ACGT")   # r2 dropped by the corrector
    m = per_base_metrics(truth, observed, corrected)
    Test.@test m.reads_dropped == 1
    Test.@test m.reads_scored == 1     # only r1 scored
    Test.@test m.tp == 1               # r1's fix still counted

    # --- Case E: a length-mismatched corrected read (Viterbi indel path) ------
    # Excluded from the positional metric, and reported.
    truth = Dict("r1" => "ACGT")
    observed = Dict("r1" => "AGGT")
    corrected = Dict("r1" => "ACGTA")  # corrector emitted a longer read
    m = per_base_metrics(truth, observed, corrected)
    Test.@test m.reads_excluded_len == 1
    Test.@test m.reads_scored == 0
    Test.@test m.total_bases == 0

    # --- Aggregation across reads --------------------------------------------
    # r1 clean fix (1 injected, 1 fixed); r2 over-correction (0 injected, 1 spurious edit).
    truth = Dict("r1" => "ACGT", "r2" => "GGGG")
    observed = Dict("r1" => "AGGT", "r2" => "GGGG")
    corrected = Dict("r1" => "ACGT", "r2" => "GGCG")
    m = per_base_metrics(truth, observed, corrected)
    Test.@test m.injected == 1
    Test.@test m.total_edits == 2      # 1 true fix + 1 over-correction
    Test.@test m.tp == 1
    Test.@test m.over == 1
    Test.@test m.recall == 1.0         # the single injected error was fixed
    Test.@test m.precision == 1 / 2    # half the edits were spurious
end

Test.@testset "rhizomorph correction accuracy — simulate_substitution_reads" begin
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
        diffs = count(p -> collect(O)[p] != collect(T)[p], eachindex(collect(T)))
        per_read_injected += diffs
        # ceil(err * readlen) substitutions injected per read (distinct positions,
        # always to a different base) => exactly that many diffs.
        Test.@test diffs == ceil(Int, err * readlen)
    end
    Test.@test injected_total == per_read_injected
    Test.@test sampled_bases == length(records) * readlen
end
