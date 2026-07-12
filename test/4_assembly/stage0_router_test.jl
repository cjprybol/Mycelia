# Unit tests for the Stage 0 input→workflow router (td-26tt).
# Exercises the empty / single / sub-k / all-ambiguous read quartet and the
# fail-loud degenerate-input guard (review F4/F5) + quality-evidence gate (I1).

import Test
import Mycelia
import FASTX

const R = Mycelia.Rhizomorph

# helper: a FASTQ record with all-Q40 bases
rec(id, seq) = FASTX.FASTQ.Record(id, seq, String(fill('I', length(seq))))

Test.@testset "Stage 0 assembly router" begin

    Test.@testset "extract_input_features on a normal library" begin
        recs = [rec("r$i", "ATGCGTACGTACGTACGTACGT") for i in 1:10]
        f = R.extract_input_features(recs; k = 7)
        Test.@test f.num_reads == 10
        Test.@test f.n_informative_kmers > 0
        Test.@test f.has_quality_evidence
        Test.@test isfinite(f.estimated_coverage_depth)
        Test.@test f.estimated_coverage_depth > 0
        Test.@test 0.0 <= f.unique_to_total_kmer_ratio <= 1.0
        Test.@test f.mean_read_length == 22.0
    end

    Test.@testset "degenerate inputs signal zero informative k-mers" begin
        # empty read set
        f_empty = R.extract_input_features(FASTX.FASTQ.Record[]; k = 7)
        Test.@test f_empty.num_reads == 0
        Test.@test f_empty.n_informative_kmers == 0
        # reads shorter than k → no k-mers
        f_short = R.extract_input_features([rec("s", "ACG")]; k = 7)
        Test.@test f_short.num_reads == 1
        Test.@test f_short.n_informative_kmers == 0
        # all-ambiguous reads → UnambiguousDNAMers yields nothing
        f_ambig = R.extract_input_features([rec("n", "NNNNNNNNNN")]; k = 7)
        Test.@test f_ambig.n_informative_kmers == 0
    end

    Test.@testset "F4: select_classifier fails loud on degenerate input" begin
        f_empty = R.extract_input_features(FASTX.FASTQ.Record[]; k = 7)
        Test.@test_throws ArgumentError R.select_classifier(f_empty)
        f_short = R.extract_input_features([rec("s", "ACG")]; k = 7)
        Test.@test_throws ArgumentError R.select_classifier(f_short)
    end

    Test.@testset "select_classifier returns a valid arm for a real library" begin
        recs = [rec("r$i", "ATGCGTACGTTTGCACATCGGA") for i in 1:30]
        f = R.extract_input_features(recs; k = 7)
        sel = R.select_classifier(f)
        Test.@test sel.classifier isa R.KmerClassifier
        Test.@test sel.rationale isa AbstractString && !isempty(sel.rationale)
    end
end
