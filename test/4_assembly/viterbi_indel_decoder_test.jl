# Indel-aware pair-HMM extension to the graph-Viterbi read corrector.
#
# Milestone A: the config carries indel/gap fields that DEFAULT to the
#   substitution-collapse values (indel_moves=false, all gap masses 0), so the
#   existing decoder is byte-identical and the whole existing suite passes with
#   no edits.
# Milestone B: an inserted-unit and a deleted-unit observation each correct back
#   to the reference under an indel-enabled config (beam_width=typemax), and the
#   substitution fixture stays byte-identical with indel_moves=false.
# Milestone C: bounding knobs (adaptive band, D_max, insertion-run cap) and a
#   small nanopore-style smoke that provably routes through the indel kernel.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/viterbi_indel_decoder_test.jl")'

import BioSequences
import FASTX
import MetaGraphsNext
import Mycelia
import Test

function indel_decoded_label_strings(
        result::Mycelia.ViterbiCorrectionResult
)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

Test.@testset "Indel-aware pair-HMM Viterbi correction" begin
    Test.@testset "Milestone A: config indel fields default to substitution-collapse" begin
        cfg = Mycelia.ViterbiCorrectionConfig()
        # The gap machinery is OFF by default, so the decoder is byte-identical.
        Test.@test cfg.indel_moves == false
        Test.@test cfg.insertion_fraction == 0.0
        Test.@test cfg.deletion_fraction == 0.0
        Test.@test cfg.insertion_extend_probability == 0.0
        Test.@test cfg.deletion_extend_probability == 0.0
        Test.@test cfg.deletion_max_run == 0
        Test.@test cfg.max_insertion_run == 0
        Test.@test cfg.band_width === nothing

        # Fields are settable.
        cfg2 = Mycelia.ViterbiCorrectionConfig(
            indel_moves = true,
            insertion_fraction = 0.3,
            deletion_fraction = 0.3,
            insertion_extend_probability = 0.1,
            deletion_extend_probability = 0.1,
            deletion_max_run = 3,
            max_insertion_run = 3,
            band_width = 12
        )
        Test.@test cfg2.indel_moves == true
        Test.@test cfg2.insertion_fraction == 0.3
        Test.@test cfg2.deletion_fraction == 0.3
        Test.@test cfg2.deletion_max_run == 3
        Test.@test cfg2.max_insertion_run == 3
        Test.@test cfg2.band_width == 12
    end
end
