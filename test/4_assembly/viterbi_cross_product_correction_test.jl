# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_cross_product_correction_test.jl")'
# ```

import BioSequences
import FASTX
import Kmers
import Mycelia
import Test

struct CrossProductQualityKmerObservation{K}
    Kmer::K
    quality_scores::Vector{UInt8}
end

struct CrossProductQualitySequenceObservation{S}
    sequence::S
    quality_scores::Vector{UInt8}
end

function cross_product_quality_record(
    identifier::AbstractString, sequence::AbstractString, phred::AbstractVector{<:Integer}
)::FASTX.FASTQ.Record
    quality = String([Char(score + 33) for score in phred])
    return FASTX.FASTQ.Record(identifier, sequence, quality)
end

function cross_product_decoded_label_strings(
    result::Mycelia.ViterbiCorrectionResult
)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

function cross_product_assert_recovery(
    result::Mycelia.ViterbiCorrectionResult,
    expected::Vector{String};
    alphabet::Symbol,
    strand_mode::Symbol,
    emission_model::Symbol,
    transition_model::Symbol=:normalized_edge_weight,
)::Nothing
    path = only(result.paths)

    Test.@test result.diagnostics[:alphabet] == alphabet
    Test.@test result.diagnostics[:strand_mode] == strand_mode
    Test.@test result.diagnostics[:emission_model] == emission_model
    Test.@test result.diagnostics[:transition_model] == transition_model
    Test.@test path.diagnostics[:emission_scoring] == emission_model
    Test.@test path.diagnostics[:transition_scoring] == transition_model
    Test.@test cross_product_decoded_label_strings(result) == expected
    return nothing
end

function cross_product_error_message(callable::Function)::String
    try
        callable()
    catch error
        Test.@test error isa ArgumentError
        return sprint(showerror, error)
    end
    Test.@test false
    return ""
end

function cross_product_kmer_case(alphabet::Symbol, strand_mode::Symbol)::NamedTuple
    if alphabet == :DNA && strand_mode == :singlestrand
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id="b7_kmer_dna_single", mode=:singlestrand
        )
        observed = [
            Kmers.DNAKmer{3}("ATG"),
            Kmers.DNAKmer{3}("TGA"),
            Kmers.DNAKmer{3}("GCG"),
            Kmers.DNAKmer{3}("CGT"),
        ]
        return (graph=graph, observed=observed, expected=["ATG", "TGC", "GCG", "CGT"])
    elseif alphabet == :DNA && strand_mode in (:doublestrand, :canonical)
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id="b7_kmer_dna_$(strand_mode)", mode=strand_mode
        )
        expected_labels = [
            Kmers.DNAKmer{3}("ACG"),
            Kmers.DNAKmer{3}("CGC"),
            Kmers.DNAKmer{3}("GCA"),
            Kmers.DNAKmer{3}("CAT"),
        ]
        observed = [
            expected_labels[1],
            Kmers.DNAKmer{3}("CCC"),
            expected_labels[3],
            expected_labels[4],
        ]
        expected = if strand_mode == :canonical
            [string(BioSequences.canonical(label)) for label in expected_labels]
        else
            [string(label) for label in expected_labels]
        end
        return (graph=graph, observed=observed, expected=expected)
    elseif alphabet == :RNA && strand_mode == :singlestrand
        records = [FASTX.FASTA.Record("rna", BioSequences.rna"AUGCGU")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id="b7_kmer_rna_single", mode=:singlestrand, type_hint=:RNA
        )
        observed = [
            Kmers.RNAKmer{3}("AUG"),
            Kmers.RNAKmer{3}("UGA"),
            Kmers.RNAKmer{3}("GCG"),
            Kmers.RNAKmer{3}("CGU"),
        ]
        return (graph=graph, observed=observed, expected=["AUG", "UGC", "GCG", "CGU"])
    elseif alphabet == :RNA && strand_mode in (:doublestrand, :canonical)
        records = [FASTX.FASTA.Record("rna", BioSequences.rna"AUGCGU")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records,
            3;
            dataset_id="b7_kmer_rna_$(strand_mode)",
            mode=strand_mode,
            type_hint=:RNA,
        )
        expected_labels = [
            Kmers.RNAKmer{3}("ACG"),
            Kmers.RNAKmer{3}("CGC"),
            Kmers.RNAKmer{3}("GCA"),
            Kmers.RNAKmer{3}("CAU"),
        ]
        observed = [
            expected_labels[1],
            Kmers.RNAKmer{3}("CCC"),
            expected_labels[3],
            expected_labels[4],
        ]
        expected = if strand_mode == :canonical
            [string(BioSequences.canonical(label)) for label in expected_labels]
        else
            [string(label) for label in expected_labels]
        end
        return (graph=graph, observed=observed, expected=expected)
    elseif alphabet == :AA && strand_mode == :singlestrand
        records = [FASTX.FASTA.Record("aa", BioSequences.aa"MKWVTF")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id="b7_kmer_aa_single", mode=:singlestrand, type_hint=:AA
        )
        observed = [
            Kmers.AAKmer{3}("MKW"),
            Kmers.AAKmer{3}("KWA"),
            Kmers.AAKmer{3}("WVT"),
            Kmers.AAKmer{3}("VTF"),
        ]
        return (graph=graph, observed=observed, expected=["MKW", "KWV", "WVT", "VTF"])
    end
    throw(ArgumentError("unsupported k-mer B7 case: $alphabet/$strand_mode"))
end

function cross_product_qualmer_case(alphabet::Symbol, strand_mode::Symbol)::NamedTuple
    if alphabet == :DNA
        records = [
            cross_product_quality_record("dna_$index", "ATGCGT", [40, 40, 40, 40, 40, 40])
            for index in 1:8
        ]
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; dataset_id="b7_qualmer_dna_$(strand_mode)", mode=strand_mode
        )
        if strand_mode == :singlestrand
            observed = [
                CrossProductQualityKmerObservation(
                    Kmers.DNAKmer{3}("ATG"), UInt8[40, 40, 40]
                ),
                CrossProductQualityKmerObservation(Kmers.DNAKmer{3}("TGA"), UInt8[2, 2, 2]),
                CrossProductQualityKmerObservation(
                    Kmers.DNAKmer{3}("GCG"), UInt8[40, 40, 40]
                ),
                CrossProductQualityKmerObservation(
                    Kmers.DNAKmer{3}("CGT"), UInt8[40, 40, 40]
                ),
            ]
            expected = ["ATG", "TGC", "GCG", "CGT"]
        else
            expected_labels = [
                Kmers.DNAKmer{3}("ACG"),
                Kmers.DNAKmer{3}("CGC"),
                Kmers.DNAKmer{3}("GCA"),
                Kmers.DNAKmer{3}("CAT"),
            ]
            observed = [
                CrossProductQualityKmerObservation(expected_labels[1], UInt8[40, 40, 40]),
                CrossProductQualityKmerObservation(Kmers.DNAKmer{3}("CCC"), UInt8[2, 2, 2]),
                CrossProductQualityKmerObservation(expected_labels[3], UInt8[40, 40, 40]),
                CrossProductQualityKmerObservation(expected_labels[4], UInt8[40, 40, 40]),
            ]
            expected = if strand_mode == :canonical
                [string(BioSequences.canonical(label)) for label in expected_labels]
            else
                [string(label) for label in expected_labels]
            end
        end
        return (graph=graph, observed=observed, expected=expected)
    elseif alphabet == :RNA
        records = [
            cross_product_quality_record("rna_$index", "AUGCGU", [40, 40, 40, 40, 40, 40])
            for index in 1:8
        ]
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records,
            3;
            dataset_id="b7_qualmer_rna_$(strand_mode)",
            mode=strand_mode,
            type_hint=:RNA,
        )
        if strand_mode == :singlestrand
            observed = [
                CrossProductQualityKmerObservation(
                    Kmers.RNAKmer{3}("AUG"), UInt8[40, 40, 40]
                ),
                CrossProductQualityKmerObservation(Kmers.RNAKmer{3}("UGA"), UInt8[2, 2, 2]),
                CrossProductQualityKmerObservation(
                    Kmers.RNAKmer{3}("GCG"), UInt8[40, 40, 40]
                ),
                CrossProductQualityKmerObservation(
                    Kmers.RNAKmer{3}("CGU"), UInt8[40, 40, 40]
                ),
            ]
            expected = ["AUG", "UGC", "GCG", "CGU"]
        else
            expected_labels = [
                Kmers.RNAKmer{3}("ACG"),
                Kmers.RNAKmer{3}("CGC"),
                Kmers.RNAKmer{3}("GCA"),
                Kmers.RNAKmer{3}("CAU"),
            ]
            observed = [
                CrossProductQualityKmerObservation(expected_labels[1], UInt8[40, 40, 40]),
                CrossProductQualityKmerObservation(Kmers.RNAKmer{3}("CCC"), UInt8[2, 2, 2]),
                CrossProductQualityKmerObservation(expected_labels[3], UInt8[40, 40, 40]),
                CrossProductQualityKmerObservation(expected_labels[4], UInt8[40, 40, 40]),
            ]
            expected = if strand_mode == :canonical
                [string(BioSequences.canonical(label)) for label in expected_labels]
            else
                [string(label) for label in expected_labels]
            end
        end
        return (graph=graph, observed=observed, expected=expected)
    elseif alphabet == :AA && strand_mode == :singlestrand
        records = [
            cross_product_quality_record("aa_$index", "MKWVTF", [40, 40, 40, 40, 40, 40])
            for index in 1:8
        ]
        graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; dataset_id="b7_qualmer_aa_single", mode=:singlestrand, type_hint=:AA
        )
        observed = [
            CrossProductQualityKmerObservation(Kmers.AAKmer{3}("MKW"), UInt8[40, 40, 40]),
            CrossProductQualityKmerObservation(Kmers.AAKmer{3}("KWA"), UInt8[2, 2, 2]),
            CrossProductQualityKmerObservation(Kmers.AAKmer{3}("WVT"), UInt8[40, 40, 40]),
            CrossProductQualityKmerObservation(Kmers.AAKmer{3}("VTF"), UInt8[40, 40, 40]),
        ]
        return (graph=graph, observed=observed, expected=["MKW", "KWV", "WVT", "VTF"])
    end
    throw(ArgumentError("unsupported qualmer B7 case: $alphabet/$strand_mode"))
end

function cross_product_fasta_case(alphabet::Symbol)::NamedTuple
    if alphabet == :DNA
        records = [
            FASTX.FASTA.Record("read_1", BioSequences.dna"ATGCG"),
            FASTX.FASTA.Record("read_2", BioSequences.dna"GCGTA"),
            FASTX.FASTA.Record("read_3", BioSequences.dna"GTACC"),
        ]
        graph = Mycelia.Rhizomorph.build_fasta_graph_olc(records; min_overlap=3)
        observed = [
            BioSequences.dna"ATGCG", BioSequences.dna"GCGTT", BioSequences.dna"GTACC"
        ]
        return (graph=graph, observed=observed, expected=["ATGCG", "GCGTA", "GTACC"])
    elseif alphabet == :RNA
        records = [
            FASTX.FASTA.Record("read_1", BioSequences.rna"AUGCG"),
            FASTX.FASTA.Record("read_2", BioSequences.rna"GCGUA"),
            FASTX.FASTA.Record("read_3", BioSequences.rna"GUACC"),
        ]
        graph = Mycelia.Rhizomorph.build_fasta_graph_olc(
            records; min_overlap=3, type_hint=:RNA
        )
        observed = [
            BioSequences.rna"AUGCG", BioSequences.rna"GCGUU", BioSequences.rna"GUACC"
        ]
        return (graph=graph, observed=observed, expected=["AUGCG", "GCGUA", "GUACC"])
    elseif alphabet == :AA
        records = [
            FASTX.FASTA.Record("read_1", BioSequences.aa"MKWVT"),
            FASTX.FASTA.Record("read_2", BioSequences.aa"WVTFG"),
            FASTX.FASTA.Record("read_3", BioSequences.aa"TFGHI"),
        ]
        graph = Mycelia.Rhizomorph.build_fasta_graph_olc(
            records; min_overlap=3, type_hint=:AA
        )
        observed = [BioSequences.aa"MKWVT", BioSequences.aa"WVTAG", BioSequences.aa"TFGHI"]
        return (graph=graph, observed=observed, expected=["MKWVT", "WVTFG", "TFGHI"])
    end
    throw(ArgumentError("unsupported FASTA B7 alphabet: $alphabet"))
end

function cross_product_fastq_case(alphabet::Symbol)::NamedTuple
    if alphabet == :DNA
        records = [
            cross_product_quality_record("read_1", "ATGCG", [40, 40, 40, 40, 40]),
            cross_product_quality_record("read_2", "GCGTA", [40, 40, 40, 40, 40]),
            cross_product_quality_record("read_3", "GTACC", [40, 40, 40, 40, 40]),
        ]
        graph = Mycelia.Rhizomorph.build_fastq_graph_olc(records; min_overlap=3)
        observed = [
            BioSequences.dna"ATGCG",
            CrossProductQualitySequenceObservation(
                BioSequences.dna"GCGTT", UInt8[2, 2, 2, 2, 2]
            ),
            BioSequences.dna"GTACC"
        ]
        return (graph=graph, observed=observed, expected=["ATGCG", "GCGTA", "GTACC"])
    elseif alphabet == :RNA
        records = [
            cross_product_quality_record("read_1", "AUGCG", [40, 40, 40, 40, 40]),
            cross_product_quality_record("read_2", "GCGUA", [40, 40, 40, 40, 40]),
            cross_product_quality_record("read_3", "GUACC", [40, 40, 40, 40, 40]),
        ]
        graph = Mycelia.Rhizomorph.build_fastq_graph_olc(records; min_overlap=3)
        observed = [
            BioSequences.rna"AUGCG",
            CrossProductQualitySequenceObservation(
                BioSequences.rna"GCGUU", UInt8[2, 2, 2, 2, 2]
            ),
            BioSequences.rna"GUACC"
        ]
        return (graph=graph, observed=observed, expected=["AUGCG", "GCGUA", "GUACC"])
    elseif alphabet == :AA
        records = [
            cross_product_quality_record("read_1", "MKWVT", [40, 40, 40, 40, 40]),
            cross_product_quality_record("read_2", "WVTFG", [40, 40, 40, 40, 40]),
            cross_product_quality_record("read_3", "TFGHI", [40, 40, 40, 40, 40]),
        ]
        graph = Mycelia.Rhizomorph.build_fastq_graph_olc(records; min_overlap=3)
        observed = [
            BioSequences.aa"MKWVT",
            CrossProductQualitySequenceObservation(
                BioSequences.aa"WVTAG", UInt8[2, 2, 2, 2, 2]
            ),
            BioSequences.aa"TFGHI"
        ]
        return (graph=graph, observed=observed, expected=["MKWVT", "WVTFG", "TFGHI"])
    end
    throw(ArgumentError("unsupported FASTQ B7 alphabet: $alphabet"))
end

Test.@testset "Viterbi cross-product correction matrix" begin
    fixed_sequence_cases = [
        (graph_type=:kmer, alphabet=:DNA, strand_mode=:singlestrand),
        (graph_type=:kmer, alphabet=:DNA, strand_mode=:doublestrand),
        (graph_type=:kmer, alphabet=:DNA, strand_mode=:canonical),
        (graph_type=:kmer, alphabet=:RNA, strand_mode=:singlestrand),
        (graph_type=:kmer, alphabet=:RNA, strand_mode=:doublestrand),
        (graph_type=:kmer, alphabet=:RNA, strand_mode=:canonical),
        (graph_type=:kmer, alphabet=:AA, strand_mode=:singlestrand),
    ]

    quality_kmer_cases = [
        (graph_type=:qualmer, alphabet=:DNA, strand_mode=:singlestrand),
        (graph_type=:qualmer, alphabet=:DNA, strand_mode=:doublestrand),
        (graph_type=:qualmer, alphabet=:DNA, strand_mode=:canonical),
        (graph_type=:qualmer, alphabet=:RNA, strand_mode=:singlestrand),
        (graph_type=:qualmer, alphabet=:RNA, strand_mode=:doublestrand),
        (graph_type=:qualmer, alphabet=:RNA, strand_mode=:canonical),
        (graph_type=:qualmer, alphabet=:AA, strand_mode=:singlestrand),
    ]

    fasta_cases = [
        (graph_type=:fasta, alphabet=:DNA, strand_mode=:singlestrand),
        (graph_type=:fasta, alphabet=:RNA, strand_mode=:singlestrand),
        (graph_type=:fasta, alphabet=:AA, strand_mode=:singlestrand),
    ]

    fastq_cases = [
        (graph_type=:fastq, alphabet=:DNA, strand_mode=:singlestrand),
        (graph_type=:fastq, alphabet=:RNA, strand_mode=:singlestrand),
        (graph_type=:fastq, alphabet=:AA, strand_mode=:singlestrand),
    ]

    text_cases = [
        (graph_type=:ngram, alphabet=:TEXT, strand_mode=:singlestrand),
        (graph_type=:string, alphabet=:TEXT, strand_mode=:singlestrand),
    ]

    valid_cases = vcat(
        fixed_sequence_cases, quality_kmer_cases, fasta_cases, fastq_cases, text_cases
    )
    Test.@test length(valid_cases) == 22

    for case in fixed_sequence_cases
        Test.@testset "$(case.graph_type) $(case.alphabet) $(case.strand_mode)" begin
            fixture = cross_product_kmer_case(case.alphabet, case.strand_mode)
            result = Mycelia.correct_observations(fixture.graph, [fixture.observed])
            cross_product_assert_recovery(
                result,
                fixture.expected;
                alphabet=case.alphabet,
                strand_mode=case.strand_mode,
                emission_model=:alphabet_parameterized,
            )
        end
    end

    for case in quality_kmer_cases
        Test.@testset "$(case.graph_type) $(case.alphabet) $(case.strand_mode)" begin
            fixture = cross_product_qualmer_case(case.alphabet, case.strand_mode)
            result = Mycelia.correct_observations(fixture.graph, [fixture.observed])
            cross_product_assert_recovery(
                result,
                fixture.expected;
                alphabet=case.alphabet,
                strand_mode=case.strand_mode,
                emission_model=:quality_aware,
            )
        end
    end

    for case in fasta_cases
        Test.@testset "$(case.graph_type) $(case.alphabet) $(case.strand_mode)" begin
            fixture = cross_product_fasta_case(case.alphabet)
            result = Mycelia.correct_observations(fixture.graph, [fixture.observed])
            cross_product_assert_recovery(
                result,
                fixture.expected;
                alphabet=case.alphabet,
                strand_mode=:singlestrand,
                emission_model=:alphabet_parameterized,
                transition_model=:normalized_overlap_length,
            )
        end
    end

    for case in fastq_cases
        Test.@testset "$(case.graph_type) $(case.alphabet) $(case.strand_mode)" begin
            fixture = cross_product_fastq_case(case.alphabet)
            result = Mycelia.correct_observations(fixture.graph, [fixture.observed])
            cross_product_assert_recovery(
                result,
                fixture.expected;
                alphabet=case.alphabet,
                strand_mode=:singlestrand,
                emission_model=:quality_aware,
                transition_model=:normalized_overlap_length,
            )
        end
    end

    Test.@testset "ngram TEXT singlestrand" begin
        graph = Mycelia.Rhizomorph.build_ngram_graph(
            ["ABCDEF"], 3; dataset_id="b7_ngram_text"
        )
        observed = ["ABC", "BCX", "CDE", "DEF"]
        result = Mycelia.correct_observations(graph, [observed])
        cross_product_assert_recovery(
            result,
            ["ABC", "BCD", "CDE", "DEF"];
            alphabet=:TEXT,
            strand_mode=:singlestrand,
            emission_model=:alphabet_parameterized,
        )
    end

    Test.@testset "string TEXT singlestrand" begin
        graph = Mycelia.Rhizomorph.build_string_graph(
            ["ABCDE", "CDEFG", "EFGHI"]; dataset_id="b7_string_text", min_overlap=3
        )
        observed = ["ABCDE", "CDXFG", "EFGHI"]
        result = Mycelia.correct_observations(graph, [observed])
        cross_product_assert_recovery(
            result,
            ["ABCDE", "CDEFG", "EFGHI"];
            alphabet=:TEXT,
            strand_mode=:singlestrand,
            emission_model=:alphabet_parameterized,
            transition_model=:normalized_overlap_length,
        )
    end

    Test.@testset "invalid reverse-complement combinations are rejected clearly" begin
        aa_fixture = cross_product_kmer_case(:AA, :singlestrand)
        aa_double_message = cross_product_error_message() do
            Mycelia.correct_observations(
                aa_fixture.graph,
                [aa_fixture.observed];
                config=Mycelia.ViterbiCorrectionConfig(strand_mode=:doublestrand),
            )
        end
        aa_canonical_message = cross_product_error_message() do
            Mycelia.correct_observations(
                aa_fixture.graph,
                [aa_fixture.observed];
                config=Mycelia.ViterbiCorrectionConfig(strand_mode=:canonical),
            )
        end

        text_graph = Mycelia.Rhizomorph.build_ngram_graph(["ABCDEF"], 3)
        text_observed = ["ABC", "BCX", "CDE", "DEF"]
        text_double_message = cross_product_error_message() do
            Mycelia.correct_observations(
                text_graph,
                [text_observed];
                config=Mycelia.ViterbiCorrectionConfig(strand_mode=:doublestrand),
            )
        end
        text_canonical_message = cross_product_error_message() do
            Mycelia.correct_observations(
                text_graph,
                [text_observed];
                config=Mycelia.ViterbiCorrectionConfig(strand_mode=:canonical),
            )
        end

        for message in (
            aa_double_message,
            aa_canonical_message,
            text_double_message,
            text_canonical_message,
        )
            Test.@test occursin("requires a DNA/RNA alphabet", message)
            Test.@test occursin("reverse-complement naive", message)
        end
        Test.@test occursin("alphabet AA", aa_double_message)
        Test.@test occursin("alphabet AA", aa_canonical_message)
        Test.@test occursin("alphabet TEXT", text_double_message)
        Test.@test occursin("alphabet TEXT", text_canonical_message)
    end
end
