# From the Mycelia base directory, run this test with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/viterbi_modality_validation_test.jl")'
# ```
#
# Cross-modality validation for the universal-sequence corrector. Extends the
# tiny-fixture pattern in viterbi_alphabet_correction_test.jl with:
#
#   1. A single-stranded RNA "molecule" fixture (~300 nt) carrying simulated
#      substitution errors, corrected under mode=:singlestrand with RNA k-mers.
#   2. A reverse-complement guard test locking the RC-naive contract: an AA
#      (protein) alphabet asked for :doublestrand or :canonical must throw.
#   3. A single-strand-vs-canonical discrimination fixture proving the two
#      strand modes produce DIFFERENT corrected paths (guards against the modes
#      silently collapsing).

import BioSequences
import FASTX
import Kmers
import Mycelia
import StableRNGs
import Test

function _mv_decoded_labels(result::Mycelia.ViterbiCorrectionResult)::Vector{String}
    decoded = only(result.corrected_observations)
    return [string(label) for label in decoded]
end

Test.@testset "Cross-modality Viterbi correction validation" begin
    Test.@testset "single-stranded RNA molecule (~300nt) recovers substitutions" begin
        # A single-stranded RNA molecule: one strand, no reverse complement.
        rng = StableRNGs.StableRNG(1234)
        bases = ('A', 'C', 'G', 'U')
        genome = String([bases[rand(rng, 1:4)] for _ in 1:300])
        rna = BioSequences.LongRNA{4}(genome)
        records = [FASTX.FASTA.Record("ssRNA_genome", rna)]

        k = 7
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, k;
            dataset_id = "ssrna_molecule",
            mode = :singlestrand,
            type_hint = :RNA
        )

        true_kmers = [
            Kmers.RNAKmer{k}(genome[i:(i + k - 1)])
            for i in 1:(length(genome) - k + 1)
        ]
        expected = [string(kmer) for kmer in true_kmers]

        # Simulate substitution errors at well-separated interior k-mers. Each
        # error is a single-base substitution at the k-mer's middle position,
        # surrounded by correct k-mers so the graph forces the true path.
        observed = copy(true_kmers)
        error_positions = (40, 120, 200)
        mid = cld(k, 2)
        for pos in error_positions
            original = string(observed[pos])
            chars = collect(original)
            chars[mid] = chars[mid] == 'A' ? 'C' : 'A'
            observed[pos] = Kmers.RNAKmer{k}(String(chars))
        end
        # Confirm the corruption actually perturbed the observations.
        Test.@test any(string(observed[pos]) != expected[pos] for pos in error_positions)

        result = Mycelia.correct_observations(graph, [observed])

        Test.@test result.diagnostics[:alphabet] == :RNA
        Test.@test result.diagnostics[:strand_mode] == :singlestrand
        Test.@test _mv_decoded_labels(result) == expected
    end

    Test.@testset "RC guard: AA alphabet rejects doublestrand and canonical" begin
        # Amino-acid k-mers are reverse-complement naive; requesting a
        # double-stranded or canonical decode must throw an ArgumentError.
        records = [FASTX.FASTA.Record("aa", BioSequences.aa"MKWVTF")]
        graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3;
            dataset_id = "aa_rc_guard",
            mode = :singlestrand,
            type_hint = :AA
        )
        observed = [
            Kmers.AAKmer{3}("MKW"),
            Kmers.AAKmer{3}("KWV"),
            Kmers.AAKmer{3}("WVT"),
            Kmers.AAKmer{3}("VTF"),
        ]

        # Sanity: the singlestrand decode is well-formed for the same fixture.
        ok = Mycelia.correct_observations(graph, [observed])
        Test.@test ok.diagnostics[:alphabet] == :AA
        Test.@test ok.diagnostics[:strand_mode] == :singlestrand

        for bad_mode in (:doublestrand, :canonical)
            local threw = false
            local message = ""
            try
                Mycelia.correct_observations(
                    graph,
                    [observed];
                    config = Mycelia.ViterbiCorrectionConfig(strand_mode = bad_mode)
                )
            catch err
                threw = true
                Test.@test err isa ArgumentError
                message = sprint(showerror, err)
            end
            Test.@test threw
            Test.@test occursin("requires a DNA/RNA alphabet", message)
            Test.@test occursin("reverse-complement naive", message)
            Test.@test occursin("alphabet AA", message)
        end
    end

    Test.@testset "ssDNA vs canonical produce different corrected paths" begin
        # The same sequence's k-mers, corrected under :singlestrand vs
        # :canonical, must yield different label paths. If the two modes
        # collapsed into one, these paths would be identical.
        records = [FASTX.FASTA.Record("dna", BioSequences.dna"ATGCGT")]

        # Single-strand: forward-oriented k-mers; one corrupted (TGA<-TGC).
        graph_ss = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id = "disc_ss", mode = :singlestrand
        )
        observed_ss = [
            Kmers.DNAKmer{3}("ATG"),
            Kmers.DNAKmer{3}("TGA"),
            Kmers.DNAKmer{3}("GCG"),
            Kmers.DNAKmer{3}("CGT"),
        ]
        result_ss = Mycelia.correct_observations(graph_ss, [observed_ss])
        labels_ss = _mv_decoded_labels(result_ss)

        # Canonical: the canonical k-mers of the same sequence; one corrupted
        # (CCC<-CGC). Expected labels are the canonical forms.
        graph_cn = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; dataset_id = "disc_cn", mode = :canonical
        )
        canonical_labels = [
            Kmers.DNAKmer{3}("ACG"),
            Kmers.DNAKmer{3}("CGC"),
            Kmers.DNAKmer{3}("GCA"),
            Kmers.DNAKmer{3}("CAT"),
        ]
        observed_cn = [
            canonical_labels[1],
            Kmers.DNAKmer{3}("CCC"),
            canonical_labels[3],
            canonical_labels[4],
        ]
        result_cn = Mycelia.correct_observations(graph_cn, [observed_cn])
        labels_cn = _mv_decoded_labels(result_cn)
        expected_cn = [string(BioSequences.canonical(label)) for label in canonical_labels]

        Test.@test result_ss.diagnostics[:strand_mode] == :singlestrand
        Test.@test result_cn.diagnostics[:strand_mode] == :canonical
        Test.@test labels_ss == ["ATG", "TGC", "GCG", "CGT"]
        Test.@test labels_cn == expected_cn
        # The load-bearing assertion: the two modes do NOT collapse.
        Test.@test labels_ss != labels_cn
    end
end
