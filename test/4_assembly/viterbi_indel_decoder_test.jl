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

    # Shared OLC fixture (mirrors viterbi_variable_length_graph_correction_test.jl):
    # a 5-vertex overlap-3 chain ATGCG -> GCGTA -> GTACC -> ACCGT -> CGTAA. The true
    # read under test covers only the FIRST THREE vertices; the graph deliberately
    # CONTINUES past that end (ACCGT, CGTAA) so a frame-shifted substitution-only
    # decode over-walks onto a wrong downstream vertex instead of coincidentally
    # stopping at the true terminal — making the insertion case discriminating.
    # The indel-enabled profile: err 0.10, indels split 30/30 of the error mass,
    # short affine gaps, exact (unbounded) beam so the assertion tests the true
    # argmax rather than the beam approximation.
    function _indel_olc_graph()
        records = [
            FASTX.FASTA.Record("read_1", BioSequences.dna"ATGCG"),
            FASTX.FASTA.Record("read_2", BioSequences.dna"GCGTA"),
            FASTX.FASTA.Record("read_3", BioSequences.dna"GTACC"),
            FASTX.FASTA.Record("read_4", BioSequences.dna"ACCGT"),
            FASTX.FASTA.Record("read_5", BioSequences.dna"CGTAA")
        ]
        return Mycelia.Rhizomorph.build_fasta_graph_olc(records; min_overlap = 3)
    end

    function _indel_config(; indel_moves = true)
        return Mycelia.ViterbiCorrectionConfig(
            error_rate = 0.10,
            indel_moves = indel_moves,
            insertion_fraction = 0.30,
            deletion_fraction = 0.30,
            insertion_extend_probability = 0.10,
            deletion_extend_probability = 0.10,
            deletion_max_run = 3,
            max_insertion_run = 3,
            beam_width = typemax(Int)
        )
    end

    Test.@testset "Milestone B: deletion-shifted read corrects back to reference" begin
        graph = _indel_olc_graph()
        # The middle read (GCGTA) is missing entirely — a unit-level deletion. The
        # terminal (GTACC) matches a graph vertex exactly, anchoring the endpoint
        # two vertices past the start. The substitution-only lock-step decoder
        # cannot reach that target in a single read step (it under-walks and never
        # reaches GTACC); the pair-HMM must take a D-move (advance the graph edge
        # ATGCG->GCGTA without consuming a read unit) and then match GTACC, landing
        # the 3-vertex reference walk.
        observed = [BioSequences.dna"ATGCG", BioSequences.dna"GTACC"]
        result = Mycelia.correct_observations(graph, [observed]; config = _indel_config())
        Test.@test indel_decoded_label_strings(result) == ["ATGCG", "GCGTA", "GTACC"]
    end

    Test.@testset "Milestone B: insertion-shifted read corrects back to reference" begin
        graph = _indel_olc_graph()
        # An extra spurious read unit (AAAAA) is inserted after the first — a
        # unit-level insertion. The last unit (GTACA) is a corrupted terminal (not a
        # graph vertex) so the decode is free-endpoint. The lock-step decoder treats
        # the junk as a real step and frame-shifts one vertex too far, over-walking
        # onto the spurious downstream vertex ACCGT ([ATGCG, GCGTA, GTACC, ACCGT]);
        # the pair-HMM must take an I-move (consume AAAAA while staying on ATGCG) and
        # recover the true 3-vertex walk, correcting GTACA back to GTACC.
        observed = [
            BioSequences.dna"ATGCG",
            BioSequences.dna"AAAAA",
            BioSequences.dna"GCGTA",
            BioSequences.dna"GTACA"
        ]
        result = Mycelia.correct_observations(graph, [observed]; config = _indel_config())
        Test.@test indel_decoded_label_strings(result) == ["ATGCG", "GCGTA", "GTACC"]
    end

    Test.@testset "Milestone B: substitution fixture stays byte-identical (indel_moves=false)" begin
        graph = _indel_olc_graph()
        # read_2 carries a substitution (GCGTT for GCGTA). Under indel_moves=false
        # the decoder is the untouched substitution path and must reproduce the
        # reference walk exactly, as the existing variable-length suite asserts.
        observed = [
            BioSequences.dna"ATGCG",
            BioSequences.dna"GCGTT",
            BioSequences.dna"GTACC"
        ]
        result = Mycelia.correct_observations(
            graph, [observed]; config = _indel_config(indel_moves = false))
        Test.@test indel_decoded_label_strings(result) == ["ATGCG", "GCGTA", "GTACC"]

        # And with indel_moves=true the same equal-length substitution fixture still
        # corrects to the reference (gap moves never beat the pure-M substitution
        # path when the read length matches the walk length).
        result_indel = Mycelia.correct_observations(
            graph, [observed]; config = _indel_config(indel_moves = true))
        Test.@test indel_decoded_label_strings(result_indel) ==
                   ["ATGCG", "GCGTA", "GTACC"]
    end
end
