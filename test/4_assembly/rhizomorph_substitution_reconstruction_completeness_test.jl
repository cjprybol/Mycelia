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
# NON-VACUITY IS A DEDICATED DEAD-END, not these seeds. The two acceptance seeds
# used to truncate because `first_evidence_strand` flipped along a read and the
# substitution walk died. Labeling mixed edges Forward removes that dead-end, so
# those seeds now complete and `substitution_length_divergences` is 0 there. The
# fail-open path is still required when a walk genuinely cannot consume the read
# (the suffix is absent from the graph). That case is its own testset: the
# counter fires, and the returned record is the original read, not a truncated
# correction. `injected > 0` remains the positive control that errors exist.
#
# WHAT THIS TEST DOES NOT PIN (verified elsewhere / deferred):
#   * WHICH of the acceptance-seed reads would have failed open — those seeds no
#     longer truncate. Byte-for-byte fail-open (ID + qualities + sequence) is the
#     dead-end testset, which calls `find_optimal_sequence_path` directly.
#   * That completeness came from PRESERVED corrections vs TOTAL pass-through — at
#     this operating point a `recall > 0` floor would be flaky; we assert
#     completeness, not a correction rate.
#   * Indel mode (deliberately exempt from the guard) — this harness is
#     substitution-only; indel-mode length-change survival is covered by
#     `indel_length_change_survival_test.jl` (td-3qjd).
#   * One geometry (500 bp genome, 150 bp reads, 10x, k=21, single pass, two error
#     rates) — matches the bead's acceptance config; broader geometries not swept.
#
# Run:
#   julia --project=. test/4_assembly/rhizomorph_substitution_reconstruction_completeness_test.jl

import Test
import Random
import BioSequences
import FASTX
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

    # These seeds no longer truncate: mixed double-strand edges are walked as
    # Forward, so the substitution path consumes the read. The fail-open counter
    # is covered by the dead-end testset below, not by requiring a truncation
    # that this geometry no longer produces.
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

# Positive control for the fail-open path. The acceptance seeds above no longer
# die mid-read, so they cannot certify that a short substitution path is
# discarded. This graph is one non-repetitive path, so its last k-mer has no
# outgoing edge. The read continues past that end; the walk stops, the counter
# fires, and the caller returns the original record instead of a truncated one.
# A repetitive reference does not work: the decoder can stay on the cycle for
# the whole read and never truncate.
Test.@testset "dead-end substitution walk fails open to the original read" begin
    k = 7
    reference = "ATGCTAGCTAGTCTAGCATG"
    reference_quality = repeat("I", length(reference))
    graph = Mycelia.Rhizomorph.build_qualmer_graph(
        [FASTX.FASTQ.Record("graph_read", reference, reference_quality)],
        k;
        mode = :doublestrand,
    )
    read_sequence = reference * "AAAAAAA"
    read = FASTX.FASTQ.Record("deadend", read_sequence, repeat("I", length(read_sequence)))
    diagnostics = Mycelia.CorrectorDiagnostics()
    improved, improvement = Mycelia.find_optimal_sequence_path(
        read, graph, k; graph_mode = :doublestrand, diagnostics = diagnostics)
    Test.@test diagnostics.substitution_length_divergences[] > 0
    Test.@test improvement == 0.0
    Test.@test FASTX.identifier(improved) == FASTX.identifier(read)
    Test.@test FASTX.sequence(String, improved) == read_sequence
    Test.@test length(FASTX.sequence(String, improved)) == length(read_sequence)
    Test.@test String(FASTX.quality(improved)) == String(FASTX.quality(read))
end
