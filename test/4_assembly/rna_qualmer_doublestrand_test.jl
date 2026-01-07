# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rna_qualmer_doublestrand_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/rna_qualmer_doublestrand_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate("../..")
## end
## using Revise

# RNA Qualmer Doublestrand Graph Test
# Tests for doublestrand RNA qualmer graph construction where forward and reverse
# complement k-mers with quality scores are preserved as separate vertices.

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "RNA Qualmer Doublestrand Graph" begin

    Test.@testset "Doublestrand Qualmer - Basic" begin
        # RNA sequence with quality scores
        seq = "AUGCAU"
        qual = [30, 35, 32, 28, 31, 29]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_doublestrand([record], 3; dataset_id="test")

        # In "AUGCAU": AUG, UGC, GCA, CAU
        # Doublestrand preserves forward + RC vertices (no canonical merge here)
        Test.@test Mycelia.Rhizomorph.vertex_count(graph) == 4
    end

    Test.@testset "Quality Scores Preserved in Canonical Form" begin
        seq = "AUGC"
        qual = [30, 35, 32, 28]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_doublestrand([record], 3; dataset_id="test")

        canon_aug = BioSequences.canonical(Kmers.RNAKmer{3}("AUG"))
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, canon_aug)

        # Should have quality evidence
        evidence_set = vertex_data.evidence["test"]["read_001"]
        Test.@test length(evidence_set) >= 1

        # At least one entry should be QualityEvidenceEntry
        has_quality = any(e -> e isa Mycelia.Rhizomorph.QualityEvidenceEntry, evidence_set)
        Test.@test has_quality
    end

    Test.@testset "Evidence Merging from Both Strands with Quality" begin
        # Forward read
        record1 = FASTX.FASTQ.Record("read_fwd", "AUGC", "IIII")  # Q40

        # Reverse read (RC of AUGC is GCAU)
        record2 = FASTX.FASTQ.Record("read_rev", "GCAU", "JJJJ")  # Q41

        graph = Mycelia.Rhizomorph.build_qualmer_graph_doublestrand([record1, record2], 3; dataset_id="test")

        # Both reads should contribute to canonical k-mers
        canon_aug = BioSequences.canonical(Kmers.RNAKmer{3}("AUG"))
        count = Mycelia.Rhizomorph.get_vertex_observation_count(graph, canon_aug)

        # Should have 2 observations (one from each read)
        Test.@test count == 2
    end

    Test.@testset "Quality Functions Work on Doublestrand" begin
        seq = "AUGC"
        qual1 = [30, 30, 30, 30]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [30, 30, 30, 30]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_doublestrand([record1, record2], 3; dataset_id="test")

        canon_aug = BioSequences.canonical(Kmers.RNAKmer{3}("AUG"))
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, canon_aug)

        # Joint quality should combine both observations
        joint_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "test")
        Test.@test !isnothing(joint_qual)

        # Two Q30 observations should combine to Q60
        Test.@test all(q == 60 for q in joint_qual)
    end

    Test.@testset "AA Sequences Rejected for Doublestrand" begin
        record = FASTX.FASTQ.Record("protein_001", "MKVLW", "IIIII")

        # AA sequences don't have reverse complement
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_qualmer_graph_doublestrand([record], 3)
    end
end

println("✓ RNA doublestrand qualmer graph tests completed")
