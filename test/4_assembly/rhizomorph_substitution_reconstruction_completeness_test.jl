# Substitution-decode reconstruction COMPLETENESS (td-pw7p.1). Regression guard
# for the defect where the ungated full-read substitution corrector emitted
# TRUNCATED / length-divergent reads, which the length-preserving scoring contract
# then silently EXCLUDED — so only a small fraction of reads were scored (e.g.
# 5/200) and no frontier verdict was valid.
#
# The contract this test pins: in substitution mode EVERY corrected record is
# length-preserving OR the read fails open to its original sequence — either way
# it keeps its identifier and is scored at input length. So `reads_scored` must
# equal the read count and every corrected length must equal its input length.
#
# TRUNCATION IS NOW ABSENT AT THIS CONFIG: until the td-tg05 decoder fix (weighted
# decode graph resolves both-strand edges to Forward), the ungated decode dead-ended
# on mixed-strand edges and the fail-open guard fired here at both error rates. It
# now fires zero times (Perlmutter job 60917657 and CI at 410695c8), so the
# real-decode case asserts `substitution_length_divergences == 0`: a return of
# truncation (td-5vak) fails loudly instead of being silently absorbed by the guard.
# The guard itself is pinned by the forced testset at the bottom of this file, which
# injects a short-path kernel and checks the read fails open and is counted.
# `injected > 0` is the positive control (errors exist to correct).
#
# WHAT THIS TEST DOES NOT PIN (verified elsewhere / deferred):
#   * WHICH reads failed open, and that a failed-open read equals its ORIGINAL
#     sequence byte-for-byte, plus quality preservation — `correct_reads_at_point`
#     returns sequences only, not records, and exposes counts not per-read flags.
#     The fail-open path returning the original record verbatim (ID + qualities +
#     sequence) is verified by code review of `try_viterbi_path_improvement` ->
#     `find_optimal_sequence_path` (returns the original `read` on `nothing`).
#   * That completeness came from PRESERVED corrections vs TOTAL pass-through — we
#     assert completeness and the absence of truncation, not a correction rate.
#   * Indel mode (deliberately exempt from the guard) — this harness is
#     substitution-only; indel-mode length-change survival is a separate test.
#   * One geometry (500 bp genome, 150 bp reads, 10x, k=21, single pass, two error
#     rates) — matches the bead's acceptance config; broader geometries not swept.
#
# Run:
#   julia --project=. test/4_assembly/rhizomorph_substitution_reconstruction_completeness_test.jl

import Test
import Random
import BioSequences
import FASTX
import MetaGraphsNext
import Mycelia

include(joinpath(@__DIR__, "..", "..", "benchmarking", "rhizomorph_correction_pr_curve.jl"))

# Substitution-only, ungated full-read decode — the exact operating point whose
# truncation the fix addresses. Single pass keeps the test CI-fast; the per-read
# Viterbi decode (and thus the truncation the fix guards) runs every pass.
const _COMPLETENESS_POINT = (
    name = "noskip+nogate", skip_solid = false, cheap_correct = true,
    hard_window = false, soft_em = false, n_k_rungs = 1, max_iterations_per_k = 1,
    beam_width = nothing, calibrated_gap_threshold = nothing)

function _run_completeness_case(; err::Float64, seed::Int)
    k = 21
    rng = Random.MersenneTwister(seed)
    genome = BioSequences.randdnaseq(rng, 500)
    records, truth_by_id, observed_by_id,
    injected, _ = simulate_substitution_reads(genome, 150, 10, err, rng)

    Test.@test injected > 0                     # positive control: errors exist to correct
    Test.@test !isempty(records)

    input_len_by_id = Dict(FASTX.identifier(r) => length(FASTX.sequence(String, r))
    for r in records)

    Random.seed!(seed)
    corrected,
    corrector_errors = correct_reads_at_point(
        records, k, _COMPLETENESS_POINT; return_corrector_errors = true)

    # (1) IDENTIFIER PRESERVATION: every input read is present in the output.
    for id in keys(input_len_by_id)
        Test.@test haskey(corrected, id)
    end

    # (2) LENGTH-PRESERVING CONTRACT (the non-vacuous core): every corrected
    # substitution record equals its input length — no truncated/length-divergent
    # reads survive to the output. Violated pre-fix for the truncated reads.
    n_length_matched = 0
    for (id, corrected_seq) in corrected
        if haskey(input_len_by_id, id)
            Test.@test length(corrected_seq) == input_len_by_id[id]
            n_length_matched += length(corrected_seq) == input_len_by_id[id] ? 1 : 0
        end
    end
    Test.@test n_length_matched == length(input_len_by_id)

    # (3) COMPLETENESS: with lengths preserved, per_base_metrics scores every read
    # (no length exclusions) — the property the frontier verdict depends on.
    m = per_base_metrics(truth_by_id, observed_by_id, corrected)
    Test.@test m.reads_scored == length(records)

    # (4) NO TRUNCATION: after td-tg05 every decode at this config spans the read,
    # so the fail-open guard never fires. (2)/(3) would still pass if truncation
    # came back (the guard absorbs it), so this is what surfaces a regression.
    # Guard coverage lives in the forced testset below.
    Test.@test get(corrector_errors, :substitution_length_divergences, 0) == 0

    # (5) TELEMETRY-DRIFT GUARD: the exported corrector-error dict must expose one
    # key per CorrectorDiagnostics field. If a future field is added to the struct
    # but not to the finalize export dict (a real, previously-latent omission —
    # gate_skipped), this count mismatch catches it.
    Test.@test length(corrector_errors) == fieldcount(Mycelia.CorrectorDiagnostics)
end

Test.@testset "substitution reconstruction completeness (fail-open length contract)" begin
    # Both operating points named in the bead acceptance criteria.
    _run_completeness_case(err = 0.05, seed = 20260713)
    _run_completeness_case(err = 0.10, seed = 20260714)
end

# A substitution kernel that returns a one-unit path for a multi-unit read: the
# shape a dead-ended (truncated) decode produced before td-tg05.
function _short_path_substitution_kernel(
        graph,
        observations::AbstractVector;
        config::Mycelia.ViterbiCorrectionConfig,
        weighted_graph::Any = nothing,
)::Mycelia.ViterbiCorrectionResult
    @assert length(only(observations)) > 1
    path = Mycelia.Rhizomorph.ViterbiDecodingResult(
        nothing, 0.0, Dict{Symbol, Any}())
    return Mycelia.ViterbiCorrectionResult(
        Any[Any[first(MetaGraphsNext.labels(graph))]],
        Mycelia.Rhizomorph.ViterbiDecodingResult[path],
        Dict{Symbol, Any}(),
    )
end

Test.@testset "substitution fail-open guard fires on a short decode (forced)" begin
    sequence = "ACGTTGCAAGGCTTACG"
    read = FASTX.FASTQ.Record("r", sequence, repeat("I", length(sequence)))
    graph = Mycelia.Rhizomorph.build_qualmer_graph([read], 5; mode = :doublestrand)

    # Control: the real kernel on a read the graph was built from spans the read,
    # so the guard does not fire.
    control = Mycelia.CorrectorDiagnostics()
    Mycelia.try_viterbi_path_improvement(
        read, graph, 5; graph_mode = :doublestrand, diagnostics = control)
    Test.@test control.substitution_length_divergences[] == 0

    # Forced truncation: the read fails open (`nothing` -> caller keeps the
    # original record) and the divergence is counted exactly once.
    forced = Mycelia.CorrectorDiagnostics()
    result = Mycelia.try_viterbi_path_improvement(
        read, graph, 5; graph_mode = :doublestrand, diagnostics = forced,
        correction_kernel = _short_path_substitution_kernel)
    Test.@test result === nothing
    Test.@test forced.substitution_length_divergences[] == 1
    Test.@test forced.structural_errors[] == 0
end
