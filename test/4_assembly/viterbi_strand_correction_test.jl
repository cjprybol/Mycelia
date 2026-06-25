# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_strand_correction_test.jl")'
# ```

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

function _strand_decoded_labels(result::Mycelia.ViterbiCorrectionResult)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

function _strand_path_result(
        result::Mycelia.ViterbiCorrectionResult
)::Mycelia.Rhizomorph.ViterbiDecodingResult
    return only(result.paths)
end

function _canonical_label_strings(labels::AbstractVector)::Vector{String}
    return [string(BioSequences.canonical(label)) for label in labels]
end

Test.@testset "Strand-aware Viterbi correction" begin
    Test.@testset "BioSequences RNA reverse-complement behavior is documented" begin
        rna = BioSequences.rna"AUGC"
        rna_rc = BioSequences.reverse_complement(rna)
        rna_kmer = Kmers.RNAKmer{3}("AUG")
        rna_kmer_rc = BioSequences.reverse_complement(rna_kmer)

        Test.@test string(rna_rc) == "GCAU"
        Test.@test typeof(rna_rc) == typeof(rna)
        Test.@test string(rna_kmer_rc) == "CAU"
        Test.@test typeof(rna_kmer_rc) == typeof(rna_kmer)
        Test.@test string(BioSequences.canonical(rna_kmer)) == "AUG"
    end

    Test.@testset "DNA doublestrand correction can start on reverse-complement strand" begin
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records,
            3;
            dataset_id = "viterbi_dna_ds",
            mode = :doublestrand
        )
        expected = [
            Kmers.DNAKmer{3}("ACG"),
            Kmers.DNAKmer{3}("CGC"),
            Kmers.DNAKmer{3}("GCA"),
            Kmers.DNAKmer{3}("CAT")
        ]
        observed = [expected[1], Kmers.DNAKmer{3}("CCC"), expected[3], expected[4]]

        result = Mycelia.correct_observations(graph, [observed])
        path = _strand_path_result(result)

        Test.@test result.diagnostics[:alphabet] == :DNA
        Test.@test result.diagnostics[:strand_mode] == :doublestrand
        Test.@test path.diagnostics[:strand_mode] == :doublestrand
        Test.@test _strand_decoded_labels(result) == [string(label) for label in expected]
        Test.@test all(step.strand == Mycelia.Rhizomorph.Reverse for step in path.path.steps)
    end

    Test.@testset "RNA doublestrand correction uses BioSequences U-aware RC" begin
        records = [FASTX.FASTA.Record("rna", BioSequences.rna"AUGCGU")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records,
            3;
            dataset_id = "viterbi_rna_ds",
            mode = :doublestrand,
            type_hint = :RNA
        )
        expected = [
            Kmers.RNAKmer{3}("ACG"),
            Kmers.RNAKmer{3}("CGC"),
            Kmers.RNAKmer{3}("GCA"),
            Kmers.RNAKmer{3}("CAU")
        ]
        observed = [expected[1], Kmers.RNAKmer{3}("CCC"), expected[3], expected[4]]

        result = Mycelia.correct_observations(graph, [observed])
        path = _strand_path_result(result)

        Test.@test result.diagnostics[:alphabet] == :RNA
        Test.@test result.diagnostics[:strand_mode] == :doublestrand
        Test.@test path.diagnostics[:strand_mode] == :doublestrand
        Test.@test _strand_decoded_labels(result) == [string(label) for label in expected]
        Test.@test all(step.strand == Mycelia.Rhizomorph.Reverse for step in path.path.steps)
    end

    Test.@testset "DNA/RNA canonical correction scores observed RC against canonical labels" begin
        dna_records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
        dna_graph = Mycelia.Rhizomorph.build_kmer_graph(
            dna_records,
            3;
            dataset_id = "viterbi_dna_canonical",
            mode = :canonical
        )
        dna_expected_rc = [
            Kmers.DNAKmer{3}("ACG"),
            Kmers.DNAKmer{3}("CGC"),
            Kmers.DNAKmer{3}("GCA"),
            Kmers.DNAKmer{3}("CAT")
        ]
        dna_observed = [
            dna_expected_rc[1],
            Kmers.DNAKmer{3}("CCC"),
            dna_expected_rc[3],
            dna_expected_rc[4]
        ]
        dna_result = Mycelia.correct_observations(dna_graph, [dna_observed])

        Test.@test dna_result.diagnostics[:strand_mode] == :canonical
        Test.@test _strand_path_result(dna_result).diagnostics[:strand_mode] == :canonical
        Test.@test _strand_decoded_labels(dna_result) == _canonical_label_strings(dna_expected_rc)

        rna_records = [FASTX.FASTA.Record("rna", BioSequences.rna"AUGCGU")]
        rna_graph = Mycelia.Rhizomorph.build_kmer_graph(
            rna_records,
            3;
            dataset_id = "viterbi_rna_canonical",
            mode = :canonical,
            type_hint = :RNA
        )
        rna_expected_rc = [
            Kmers.RNAKmer{3}("ACG"),
            Kmers.RNAKmer{3}("CGC"),
            Kmers.RNAKmer{3}("GCA"),
            Kmers.RNAKmer{3}("CAU")
        ]
        rna_observed = [
            rna_expected_rc[1],
            Kmers.RNAKmer{3}("CCC"),
            rna_expected_rc[3],
            rna_expected_rc[4]
        ]
        rna_result = Mycelia.correct_observations(rna_graph, [rna_observed])

        Test.@test rna_result.diagnostics[:alphabet] == :RNA
        Test.@test rna_result.diagnostics[:strand_mode] == :canonical
        Test.@test _strand_decoded_labels(rna_result) == _canonical_label_strings(rna_expected_rc)
    end

    Test.@testset "AA and text remain reverse-complement naive singlestrand" begin
        aa_records = [FASTX.FASTA.Record("aa", BioSequences.aa"MKWVTF")]
        aa_graph = Mycelia.Rhizomorph.build_kmer_graph(
            aa_records,
            3;
            dataset_id = "viterbi_aa_single",
            mode = :singlestrand,
            type_hint = :AA
        )
        aa_observed = [
            Kmers.AAKmer{3}("MKW"),
            Kmers.AAKmer{3}("KWA"),
            Kmers.AAKmer{3}("WVT"),
            Kmers.AAKmer{3}("VTF")
        ]
        aa_result = Mycelia.correct_observations(aa_graph, [aa_observed])

        Test.@test aa_result.diagnostics[:alphabet] == :AA
        Test.@test aa_result.diagnostics[:strand_mode] == :singlestrand
        Test.@test aa_result.diagnostics[:reverse_complement_support] == false
        Test.@test _strand_decoded_labels(aa_result) == ["MKW", "KWV", "WVT", "VTF"]
        Test.@test_throws ArgumentError Mycelia.correct_observations(
            aa_graph,
            [aa_observed];
            config = Mycelia.ViterbiCorrectionConfig(strand_mode = :doublestrand)
        )

        text_graph = Mycelia.Rhizomorph.build_ngram_graph(["a1b2c3"], 2)
        text_observed = ["a1", "1x", "b2", "2c", "c3"]
        text_result = Mycelia.correct_observations(text_graph, [text_observed])

        Test.@test text_result.diagnostics[:alphabet] == :TEXT
        Test.@test text_result.diagnostics[:strand_mode] == :singlestrand
        Test.@test text_result.diagnostics[:reverse_complement_support] == false
        Test.@test _strand_decoded_labels(text_result) == ["a1", "1b", "b2", "2c", "c3"]
        Test.@test_throws ArgumentError Mycelia.correct_observations(
            text_graph,
            [text_observed];
            config = Mycelia.ViterbiCorrectionConfig(strand_mode = :canonical)
        )
    end
end
