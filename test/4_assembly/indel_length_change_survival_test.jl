# Indel-mode exemption from the substitution fail-open length guard (td-3qjd).
#
# The td-pw7p.1 guard (`indel_params === nothing && length diverges`) is
# deliberately skipped in indel mode, because length-changing corrections are
# the point. That exemption was review-verified but untested; a regression that
# applied the guard to indel mode would fail-open every length-changing
# correction and still pass the substitution-only completeness suite.
#
# This test calls `try_viterbi_path_improvement` with `indel_params` and the
# C1 single-base insertion/deletion k-mer fixtures, and asserts the
# length-changed record is kept (not failed open).
#
# Run:
#   julia --project=. -e 'include("test/4_assembly/indel_length_change_survival_test.jl")'

import Test
import FASTX
import Mycelia

function _survival_quality_record(identifier::AbstractString, sequence::AbstractString)
    return FASTX.FASTQ.Record(identifier, sequence, repeat("I", length(sequence)))
end

Test.@testset "indel-mode length-changing corrections survive the substitution fail-open guard (td-3qjd)" begin
    k = 7
    reference = Mycelia.random_fasta_record(moltype = :DNA, seed = 1, L = 40)
    reference_string = FASTX.sequence(String, reference)
    graph = Mycelia.Rhizomorph.build_kmer_graph_singlestrand([reference], k)
    params = Mycelia.IndelDecodeParams(
        0.10, 0.30, 0.30, 0.10, 0.10, 3, 3, nothing)

    indel_position = 20
    cases = (
        "insertion" => reference_string[1:indel_position] * "A" *
                       reference_string[(indel_position + 1):end],
        "deletion" => reference_string[1:(indel_position - 1)] *
                      reference_string[(indel_position + 1):end],
    )

    for (label, observed_string) in cases
        Test.@testset "$(label) is kept at a new length" begin
            Test.@test length(observed_string) != length(reference_string)
            observed = _survival_quality_record("obs_$(label)", observed_string)
            diagnostics = Mycelia.CorrectorDiagnostics()
            result = Mycelia.try_viterbi_path_improvement(
                observed, graph, k;
                graph_mode = :singlestrand,
                beam_width = typemax(Int),
                diagnostics = diagnostics,
                indel_params = params,
            )
            Test.@test result !== nothing
            improved_record, _likelihood = result
            improved_seq = FASTX.sequence(String, improved_record)
            Test.@test FASTX.identifier(improved_record) == FASTX.identifier(observed)
            Test.@test length(improved_seq) != length(observed_string)
            Test.@test length(improved_seq) == length(reference_string)
            Test.@test improved_seq == reference_string
            Test.@test diagnostics.substitution_length_divergences[] == 0
            Test.@test diagnostics.indel_attempts[] == 1
            Test.@test diagnostics.indel_decodes[] == 1
            Test.@test diagnostics.truncated_decodes[] == 0
            Test.@test diagnostics.trace_contract_errors[] == 0
        end
    end
end

println("✓ Indel-mode length-change survival tests completed")
