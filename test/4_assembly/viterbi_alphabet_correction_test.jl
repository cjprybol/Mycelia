# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_alphabet_correction_test.jl")'
# ```

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

function _decoded_labels(result::Mycelia.ViterbiCorrectionResult)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

function _path_result(
        result::Mycelia.ViterbiCorrectionResult
)::Mycelia.Rhizomorph.ViterbiDecodingResult
    return only(result.paths)
end

Test.@testset "Alphabet-parameterized Viterbi correction" begin
    Test.@testset "DNA k-mer emissions recover substitution errors" begin
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records,
            3;
            dataset_id = "viterbi_dna",
            mode = :singlestrand
        )
        observed = [
            Kmers.DNAKmer{3}("ATG"),
            Kmers.DNAKmer{3}("TGA"),
            Kmers.DNAKmer{3}("GCG"),
            Kmers.DNAKmer{3}("CGT")
        ]

        result = Mycelia.correct_observations(graph, [observed])
        path = _path_result(result)

        Test.@test result.diagnostics[:alphabet] == :DNA
        Test.@test result.diagnostics[:emission_model] == :alphabet_parameterized
        Test.@test path.diagnostics[:emission_scoring] == :alphabet_parameterized
        Test.@test path.diagnostics[:alphabet] == :DNA
        Test.@test _decoded_labels(result) == ["ATG", "TGC", "GCG", "CGT"]
    end

    Test.@testset "RNA k-mer emissions recover substitution errors" begin
        records = [FASTX.FASTA.Record("rna", BioSequences.rna"AUGCGU")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records,
            3;
            dataset_id = "viterbi_rna",
            mode = :singlestrand,
            type_hint = :RNA
        )
        observed = [
            Kmers.RNAKmer{3}("AUG"),
            Kmers.RNAKmer{3}("UGA"),
            Kmers.RNAKmer{3}("GCG"),
            Kmers.RNAKmer{3}("CGU")
        ]

        result = Mycelia.correct_observations(graph, [observed])
        path = _path_result(result)

        Test.@test result.diagnostics[:alphabet] == :RNA
        Test.@test path.diagnostics[:alphabet] == :RNA
        Test.@test _decoded_labels(result) == ["AUG", "UGC", "GCG", "CGU"]
    end

    Test.@testset "AA k-mer emissions recover substitution errors" begin
        records = [FASTX.FASTA.Record("aa", BioSequences.aa"MKWVTF")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records,
            3;
            dataset_id = "viterbi_aa",
            mode = :singlestrand,
            type_hint = :AA
        )
        observed = [
            Kmers.AAKmer{3}("MKW"),
            Kmers.AAKmer{3}("KWA"),
            Kmers.AAKmer{3}("WVT"),
            Kmers.AAKmer{3}("VTF")
        ]

        result = Mycelia.correct_observations(graph, [observed])
        path = _path_result(result)

        Test.@test result.diagnostics[:alphabet] == :AA
        Test.@test path.diagnostics[:alphabet] == :AA
        Test.@test _decoded_labels(result) == ["MKW", "KWV", "WVT", "VTF"]
    end

    Test.@testset "alphabet-specific substitution costs" begin
        dna_mismatch = Mycelia.default_viterbi_emission_logp(
            Kmers.DNAKmer{3}("TGA"),
            Kmers.DNAKmer{3}("TGC"),
            :DNA;
            error_rate = 0.03
        )
        aa_mismatch = Mycelia.default_viterbi_emission_logp(
            Kmers.AAKmer{3}("KWA"),
            Kmers.AAKmer{3}("KWV"),
            :AA;
            error_rate = 0.03
        )

        Test.@test dna_mismatch > aa_mismatch
        Test.@test_throws ArgumentError Mycelia.default_viterbi_emission_logp(
            Kmers.DNAKmer{3}("TGA"),
            Kmers.DNAKmer{3}("TGC"),
            :RNA;
            error_rate = 0.03
        )
    end
end
