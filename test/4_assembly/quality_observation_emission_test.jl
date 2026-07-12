# CHECKPOINT test for the graph-as-HMM correction core (td-jqdf / td-tps5):
# the Viterbi emission must respond to the READ's per-base Phred quality, carried
# by `Mycelia.QualityObservation`, not just the graph's population-average quality.
#
# Two levels of evidence:
#   1. Emission logp: the SAME observed k-mer against the SAME graph node yields a
#      DIFFERENT emission log-probability when the read base at the mismatch is
#      low-quality vs high-quality (low quality => smaller mismatch penalty).
#   2. Correction decision: the SAME observed error k-mer is corrected toward the
#      well-supported graph vertex when the read base is low-quality, but kept when
#      the read base is high-quality — a decision driven purely by per-read quality.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/quality_observation_emission_test.jl")'

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

function _quality_record(
        identifier::AbstractString,
        sequence::AbstractString,
        phred::AbstractVector{<:Integer}
)::FASTX.FASTQ.Record
    quality = String([Char(score + 33) for score in phred])
    return FASTX.FASTQ.Record(identifier, sequence, quality)
end

function _decoded_strings(result::Mycelia.ViterbiCorrectionResult)::Vector{String}
    return [string(label) for label in only(result.corrected_observations)]
end

Test.@testset "QualityObservation per-read emission (td-jqdf / td-tps5)" begin
    Test.@testset "wrapper accessors delegate to the k-mer and expose raw Phred" begin
        qo = Mycelia.QualityObservation(Kmers.DNAKmer{3}("TGA"), UInt8[40, 30, 2])
        # string + alphabet resolution see through the wrapper to the k-mer
        Test.@test Mycelia._viterbi_unit_string(qo) == "TGA"
        Test.@test Mycelia._infer_viterbi_alphabet(qo) == :DNA
        # raw Phred is returned verbatim (no ASCII decode heuristic)
        Test.@test Mycelia._viterbi_direct_quality_scores(qo) == Float64[40.0, 30.0, 2.0]
        # label-comparison unwraps to the k-mer; bare units pass through unchanged
        Test.@test Mycelia._viterbi_label_unit(qo) == Kmers.DNAKmer{3}("TGA")
        Test.@test Mycelia._viterbi_label_unit(Kmers.DNAKmer{3}("TGA")) ==
                   Kmers.DNAKmer{3}("TGA")
    end

    Test.@testset "same k-mer + node: low-Q base gives a higher emission logp than high-Q" begin
        observed = Kmers.DNAKmer{3}("TGA")
        node = Kmers.DNAKmer{3}("TGC")  # mismatch only at the last base

        qo_low = Mycelia.QualityObservation(observed, UInt8[40, 40, 2])
        qo_high = Mycelia.QualityObservation(observed, UInt8[40, 40, 40])

        logp_low = Mycelia.default_viterbi_emission_logp(
            qo_low,
            node,
            :DNA;
            quality = Mycelia._viterbi_direct_quality_scores(qo_low),
            error_rate = 0.01
        )
        logp_high = Mycelia.default_viterbi_emission_logp(
            qo_high,
            node,
            :DNA;
            quality = Mycelia._viterbi_direct_quality_scores(qo_high),
            error_rate = 0.01
        )

        # A low-quality read base at the mismatch means "I don't trust this base",
        # so observing TGA when the truth is TGC is cheap; a high-quality base makes
        # the same mismatch expensive. The emission therefore MUST respond to the
        # read's per-base quality.
        Test.@test logp_low > logp_high
    end

    Test.@testset "correction decision flips on per-read quality alone" begin
        # A graph with a well-supported truth vertex (TGC) and a weakly-supported
        # error branch (TGA) reachable from the same predecessor (ATG).
        records = FASTX.FASTQ.Record[]
        for index in 1:30
            push!(records, _quality_record("truth_$index", "ATGC", [40, 40, 40, 40]))
        end
        for index in 1:3
            push!(records, _quality_record("variant_$index", "ATGA", [40, 40, 40, 40]))
        end
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records,
            3;
            dataset_id = "quality_observation_decision",
            mode = :singlestrand
        )

        anchor = Mycelia.QualityObservation(Kmers.DNAKmer{3}("ATG"), UInt8[40, 40, 40])

        # SAME observed error k-mer TGA; only the read's per-base quality differs.
        low_quality = [
            anchor,
            Mycelia.QualityObservation(Kmers.DNAKmer{3}("TGA"), UInt8[2, 2, 2])
        ]
        high_quality = [
            anchor,
            Mycelia.QualityObservation(Kmers.DNAKmer{3}("TGA"), UInt8[40, 40, 40])
        ]

        low_result = Mycelia.correct_observations(graph, [low_quality])
        high_result = Mycelia.correct_observations(graph, [high_quality])

        # Low-quality error base => trust the well-supported graph vertex => correct
        # TGA -> TGC. High-quality base => keep the read's TGA. The decision is
        # driven entirely by per-read quality on an identical observed k-mer.
        Test.@test _decoded_strings(low_result) == ["ATG", "TGC"]
        Test.@test _decoded_strings(high_result) == ["ATG", "TGA"]
        Test.@test _decoded_strings(low_result) != _decoded_strings(high_result)
        Test.@test only(low_result.paths).diagnostics[:emission_scoring] == :quality_aware
    end
end
