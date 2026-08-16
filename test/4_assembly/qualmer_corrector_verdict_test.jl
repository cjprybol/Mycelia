# Failure semantics of the corrector quality-sensitivity verdicts (bead td-4e19d.2).
#
# WHAT THIS FILE EXISTS TO PREVENT: a run in which the measurements FAILED reporting a
# confident scientific verdict. `benchmarking/qualmer_corrector_quality_sensitivity.jl`
# writes the sentinel strings "ERROR"/"NA" into a row whose condition threw, and those
# sentinels are then read back by the verdict logic. Two ways that has already gone
# wrong, both caught in PR #453 review:
#
#   1. `length(unique(digests)) == 1` over `["ERROR", "ERROR", "ERROR", "ERROR"]` is
#      `true` — an all-failed run reported PERFECT AGREEMENT between conditions, which
#      is precisely the study's headline finding.
#   2. A single `"NA"` among real quality values makes `unique(...)` length 2 — a
#      failure manufactured the quality difference that the identity result needs in
#      order to mean anything.
#
# Excluding the sentinels fixes the arithmetic but not the class: reporting `false`
# instead is still a fabricated verdict in the other direction ("the assemblies DIFFER",
# "the decoder NEVER RAN" — the latter being the exact string the script's own
# INTERPRETATION GUIDE keys on to reclassify the whole study as a gating result). So the
# contract pinned here is stronger than "no sentinels in the arithmetic":
#
#   A NEGATIVE finding requires COMPLETE evidence. Where it is missing, the field is
#   "NA" — never a boolean.
#
# Run:
#   julia --project=. test/4_assembly/qualmer_corrector_verdict_test.jl

# Wrapped in a module for the same reason as track_a_baseline_benchmark_test.jl:
# runtests.jl includes every test file into one shared `Main`, and the benchmark script
# (plus the probe it includes) defines top-level consts — K, DNA_BASES, CHEMISTRIES —
# that collide with other benchmark scripts already included there.
module QualmerCorrectorVerdictTest

import Test

# PROBE_INCLUDE_ONLY suppresses the probe's own driver; the corrector script's driver is
# guarded by `abspath(PROGRAM_FILE) == @__FILE__`. Including this file therefore defines
# `corrector_verdict` and nothing else runs — no reference download, no read simulation,
# no assemblies.
PROBE_INCLUDE_ONLY = true
include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "qualmer_corrector_quality_sensitivity.jl"))

# Row builders that mirror the exact schemas stage_b_decision_probe / stage_c_end_to_end
# emit. Kept minimal but field-complete, so a schema change in the harness surfaces here
# as a build error rather than as a silently skipped assertion.
function stage_b_row(condition; differing, failures = 0, n_reads = 100)
    (stage = "B_decision", chemistry = "illumina", coverage = 10, seed = 42,
        condition = condition, n_reads = n_reads,
        reads_changed_by_decoder = differing, reads_differing_from_oracle = differing,
        n_decode_failures = failures, n_empty_decodes = 0,
        stage_b_trustworthy = failures == 0)
end

function stage_c_ok_row(condition; digest, joint_quality, decode = "[0.0, 0.4]",
        skip = "[1.0, 0.6]")
    (stage = "C_end_to_end", chemistry = "illumina", coverage = 10, seed = 42,
        condition = condition, n_reads = 100, digest = digest,
        n_contigs = 3, largest_contig = 900,
        mean_joint_quality = joint_quality, skip_fraction_per_pass = skip,
        decode_fraction_per_pass = decode, cheap_corrections_total = "12",
        wall_seconds = 1.0)
end

function stage_c_error_row(condition)
    (stage = "C_end_to_end", chemistry = "illumina", coverage = 10, seed = 42,
        condition = condition, n_reads = 100, digest = "ERROR",
        n_contigs = -1, largest_contig = -1,
        mean_joint_quality = "NA", skip_fraction_per_pass = "NA",
        decode_fraction_per_pass = "NA", cheap_corrections_total = "NA",
        wall_seconds = -1.0)
end

healthy_b_rows() = [stage_b_row(c; differing = c == "oracle" ? 0 : 7)
                    for c in CORRECTOR_CONDITIONS]

Test.@testset "corrector verdict failure semantics (td-4e19d.2)" begin

    # --- THE FINDING: an all-errors run must not produce a confident verdict --------
    Test.@testset "every condition failed => no confident verdict" begin
        c_rows = [stage_c_error_row(c) for c in CORRECTOR_CONDITIONS]
        digests = fill("ERROR", length(CORRECTOR_CONDITIONS))
        verdict = corrector_verdict("illumina", 100, healthy_b_rows(), c_rows, digests)

        Test.@test verdict.n_conditions_failed == length(CORRECTOR_CONDITIONS)
        Test.@test verdict.verdict_is_interpretable == false

        # The original defect, pinned directly: `unique(["ERROR", ...])` has one element,
        # so this field used to be `true` — an all-failed run publishing the study's
        # headline agreement result.
        Test.@test verdict.end_to_end_assemblies_identical != true
        # And the half-fix, pinned too: `false` asserts the assemblies DIFFER, which is
        # equally unsupported by zero successful assemblies.
        Test.@test verdict.end_to_end_assemblies_identical != false
        Test.@test verdict.end_to_end_assemblies_identical == "NA"

        # `decoder_ever_ran = false` is not inert either — the script's INTERPRETATION
        # GUIDE keys on it to reclassify the end-to-end result as a GATING result.
        Test.@test verdict.end_to_end_decoder_ever_ran == "NA"
        Test.@test verdict.end_to_end_quality_reached_graph == "NA"
        Test.@test verdict.end_to_end_skip_fraction == "NA"

        # Nothing on the row may be a boolean that a reader could quote as a result.
        end_to_end_fields = (verdict.end_to_end_assemblies_identical,
            verdict.end_to_end_quality_reached_graph,
            verdict.end_to_end_decoder_ever_ran)
        Test.@test !any(f isa Bool for f in end_to_end_fields)
    end

    # --- A single failure is enough to withdraw a claim about ALL conditions --------
    Test.@testset "one failed condition withdraws the identity verdict" begin
        c_rows = [stage_c_ok_row("oracle"; digest = "aaa", joint_quality = "30.1"),
            stage_c_ok_row("constant40"; digest = "aaa", joint_quality = "40.0"),
            stage_c_ok_row("constant02"; digest = "aaa", joint_quality = "2.0"),
            stage_c_error_row("antioracle")]
        digests = ["aaa", "aaa", "aaa", "ERROR"]
        verdict = corrector_verdict("illumina", 100, healthy_b_rows(), c_rows, digests)

        Test.@test verdict.n_conditions_failed == 1
        Test.@test verdict.verdict_is_interpretable == false
        # Three surviving conditions agree, but "the assemblies are identical" is a
        # claim about four.
        Test.@test verdict.end_to_end_assemblies_identical == "NA"
        # A POSITIVE finding may still rest on the survivors: two distinct joint
        # qualities were observed, and no failure can un-observe them.
        Test.@test verdict.end_to_end_quality_reached_graph == true
        Test.@test verdict.end_to_end_decoder_ever_ran == true
    end

    # --- The "NA" sentinel must not manufacture a quality difference ----------------
    Test.@testset "a failed row's NA cannot create a quality difference" begin
        # Both surviving conditions reached the graph at the SAME joint quality, so
        # there is no observed difference. The failed row carries mean_joint_quality
        # "NA", which under the old `unique(...)` over all rows made the set size 2 and
        # reported a difference that never happened.
        c_rows = [stage_c_ok_row("oracle"; digest = "aaa", joint_quality = "30.0"),
            stage_c_ok_row("constant40"; digest = "aaa", joint_quality = "30.0"),
            stage_c_error_row("constant02"), stage_c_error_row("antioracle")]
        digests = ["aaa", "aaa", "ERROR", "ERROR"]
        verdict = corrector_verdict("illumina", 100, healthy_b_rows(), c_rows, digests)

        Test.@test verdict.end_to_end_quality_reached_graph != true
        Test.@test verdict.end_to_end_quality_reached_graph == "NA"
    end

    # --- Stage B's own failures gate Stage B's negative -----------------------------
    Test.@testset "a decoder that threw on every read is not a quality-blind decoder" begin
        # `reads_differing_from_oracle = 0` in every condition is bit-for-bit the
        # signature of the quality-blindness finding this script exists to establish.
        # It is also what a decoder that threw on every read produces, because the
        # failure path passes the read through unchanged.
        b_rows = [stage_b_row(c; differing = 0, failures = 100) for c in CORRECTOR_CONDITIONS]
        c_rows = [stage_c_ok_row(c; digest = "aaa", joint_quality = "30.0")
                  for c in CORRECTOR_CONDITIONS]
        verdict = corrector_verdict("illumina", 100, b_rows, c_rows,
            fill("aaa", length(CORRECTOR_CONDITIONS)))

        Test.@test verdict.n_stage_b_decode_failures == 100 * length(CORRECTOR_CONDITIONS)
        Test.@test verdict.decoder_decisions_depend_on_quality == "NA"
        Test.@test verdict.verdict_is_interpretable == false
    end

    # --- A clean run must still produce real, quotable booleans ---------------------
    # Without this, "return NA everywhere" would pass every assertion above — the
    # detector has to be shown answering, not just refusing.
    Test.@testset "a complete run still yields booleans" begin
        identical = corrector_verdict("illumina", 100, healthy_b_rows(),
            [stage_c_ok_row("oracle"; digest = "aaa", joint_quality = "30.0"),
                stage_c_ok_row("constant40"; digest = "aaa", joint_quality = "40.0"),
                stage_c_ok_row("constant02"; digest = "aaa", joint_quality = "2.0"),
                stage_c_ok_row("antioracle"; digest = "aaa", joint_quality = "5.0")],
            fill("aaa", length(CORRECTOR_CONDITIONS)))

        Test.@test identical.n_conditions_failed == 0
        Test.@test identical.n_stage_b_decode_failures == 0
        Test.@test identical.verdict_is_interpretable == true
        Test.@test identical.end_to_end_assemblies_identical === true
        Test.@test identical.end_to_end_quality_reached_graph === true
        Test.@test identical.end_to_end_decoder_ever_ran === true
        Test.@test identical.decoder_decisions_depend_on_quality === true
        Test.@test identical.end_to_end_skip_fraction == "[1.0, 0.6]"

        # ...and a complete run whose assemblies genuinely differ reports `false`, so
        # the `!= false` assertions above are pinning the ALL-ERRORS case specifically
        # and not a field that can never be false.
        differing = corrector_verdict("illumina", 100, healthy_b_rows(),
            [stage_c_ok_row("oracle"; digest = "aaa", joint_quality = "30.0"),
                stage_c_ok_row("constant40"; digest = "bbb", joint_quality = "40.0"),
                stage_c_ok_row("constant02"; digest = "ccc", joint_quality = "2.0"),
                stage_c_ok_row("antioracle"; digest = "ddd", joint_quality = "5.0")],
            ["aaa", "bbb", "ccc", "ddd"])
        Test.@test differing.verdict_is_interpretable == true
        Test.@test differing.end_to_end_assemblies_identical === false

        # A complete run in which the decoder truly never fired reports `false`, which
        # is a real negative result and must remain expressible.
        never_ran = corrector_verdict("illumina", 100,
            [stage_b_row(c; differing = 0) for c in CORRECTOR_CONDITIONS],
            [stage_c_ok_row(c; digest = "aaa", joint_quality = "30.0",
                 decode = "[0.0, 0.0]") for c in CORRECTOR_CONDITIONS],
            fill("aaa", length(CORRECTOR_CONDITIONS)))
        Test.@test never_ran.verdict_is_interpretable == true
        Test.@test never_ran.end_to_end_decoder_ever_ran === false
        Test.@test never_ran.decoder_decisions_depend_on_quality === false
        Test.@test never_ran.end_to_end_quality_reached_graph === false
    end
end

end # module
