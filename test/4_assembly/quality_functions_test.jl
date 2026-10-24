# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/quality_functions_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/quality_functions_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Quality Functions Test
# Tests for joint quality calculation and quality-aware operations.

import Test
import Mycelia
import BioSequences
import FASTX
import Kmers

Test.@testset "Quality Functions" begin
    Test.@testset "combine_phred_scores - Basic" begin
        # Two Q10 observations combine to Q20
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[10, 10])
        Test.@test combined == UInt8(20)

        # Three Q20 observations combine to Q60
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[20, 20, 20])
        Test.@test combined == UInt8(60)

        # Single observation returns same score
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[30])
        Test.@test combined == UInt8(30)
    end

    Test.@testset "combine_phred_scores - Can Exceed 60" begin
        # Three Q30 observations combine to Q90
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[30, 30, 30])
        Test.@test combined == UInt8(90)

        # Ten Q20 observations combine to Q200
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[20, 20, 20, 20, 20, 20, 20, 20, 20, 20])
        Test.@test combined == UInt8(200)
    end

    Test.@testset "combine_phred_scores - Clamps at 255" begin
        # Many high-quality observations clamp at 255
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[60, 60, 60, 60, 60])
        Test.@test combined == UInt8(255)
    end

    Test.@testset "combine_phred_scores - Empty" begin
        combined = Mycelia.Rhizomorph.combine_phred_scores(UInt8[])
        Test.@test combined == UInt8(0)
    end

    Test.@testset "get_vertex_joint_quality - Single Observation" begin
        seq = "ATGC"
        qual = [30, 30, 30, 30]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record], 3; dataset_id = "test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        joint_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "test")

        # Single observation: joint quality = original quality
        Test.@test joint_qual == UInt8[30, 30, 30]
    end

    Test.@testset "get_vertex_joint_quality - Multiple Observations" begin
        # Two observations with Q10
        seq = "ATGC"
        qual1 = [10, 10, 10, 10]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [10, 10, 10, 10]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id = "test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        joint_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "test")

        # Two Q10 observations combine to Q20
        Test.@test joint_qual == UInt8[20, 20, 20]
    end

    Test.@testset "get_vertex_joint_quality - Different Quality Profiles" begin
        # Two observations with different qualities
        seq = "ATGC"
        qual1 = [20, 30, 25, 15]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [30, 20, 25, 35]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id = "test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        joint_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "test")

        # Position-wise combination: [20+30, 30+20, 25+25]
        Test.@test joint_qual == UInt8[50, 50, 50]
    end

    Test.@testset "get_vertex_mean_quality" begin
        seq = "ATGC"
        qual1 = [20, 30, 40, 50]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [30, 30, 30, 30]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id = "test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        mean_qual = Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, "test")

        # Mean of [20,30] [30,30] [40,30]
        Test.@test mean_qual ≈ [25.0, 30.0, 35.0]
    end

    Test.@testset "get_vertex_min_quality" begin
        seq = "ATGC"
        qual1 = [20, 30, 40, 50]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [30, 20, 30, 40]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id = "test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        min_qual = Mycelia.Rhizomorph.get_vertex_min_quality(vertex_data, "test")

        # Min of [20,30] [30,20] [40,30]
        Test.@test min_qual == UInt8[20, 20, 30]
    end

    Test.@testset "filter_by_quality" begin
        # Create graph with k-mers of varying quality
        seq1 = "ATGC"  # High quality
        qual1 = [40, 40, 40, 40]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq1, qual_str1)

        seq2 = "CGTA"  # Low quality
        qual2 = [10, 10, 10, 10]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq2, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record1, record2], 3; dataset_id = "test")

        # Filter for Q30+
        high_qual = Mycelia.Rhizomorph.filter_by_quality(graph, 30, "test")

        # Should include ATG and TGC (from high-quality read)
        Test.@test Kmers.DNAKmer{3}("ATG") in high_qual
        Test.@test Kmers.DNAKmer{3}("TGC") in high_qual

        # Should NOT include CGT and GTA (from low-quality read)
        Test.@test !(Kmers.DNAKmer{3}("CGT") in high_qual)
        Test.@test !(Kmers.DNAKmer{3}("GTA") in high_qual)
    end

    Test.@testset "Quality Functions - Nonexistent Dataset" begin
        seq = "ATGC"
        qual = [30, 30, 30, 30]
        qual_str = String([Char(q + 33) for q in qual])
        record = FASTX.FASTQ.Record("read_001", seq, qual_str)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand([record], 3; dataset_id = "test")
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        # Query nonexistent dataset
        Test.@test isnothing(Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "nonexistent"))
        Test.@test isnothing(Mycelia.Rhizomorph.get_vertex_mean_quality(vertex_data, "nonexistent"))
        Test.@test isnothing(Mycelia.Rhizomorph.get_vertex_min_quality(vertex_data, "nonexistent"))
        Test.@test isnothing(Mycelia.Rhizomorph.get_vertex_posterior_log_odds(vertex_data, "nonexistent"))
    end

    Test.@testset "fragment_observation_id" begin
        Test.@test Mycelia.Rhizomorph.fragment_observation_id("fragA/1") == "fragA"
        Test.@test Mycelia.Rhizomorph.fragment_observation_id("fragA/2") == "fragA"
        Test.@test Mycelia.Rhizomorph.fragment_observation_id("fragB.1") == "fragB"
        Test.@test Mycelia.Rhizomorph.fragment_observation_id("fragB.2") == "fragB"
        Test.@test Mycelia.Rhizomorph.fragment_observation_id("fragC_R1") == "fragC"
        Test.@test Mycelia.Rhizomorph.fragment_observation_id("fragC_R2") == "fragC"
        Test.@test Mycelia.Rhizomorph.fragment_observation_id("plain_read") == "plain_read"
        # Only trailing suffixes strip; embedded `/1` mid-string is preserved.
        Test.@test Mycelia.Rhizomorph.fragment_observation_id("x/1y") == "x/1y"
    end

    Test.@testset "collapse_molecule_qualities - mean of error probabilities" begin
        # Two entries, raw Q10 and Q20 at each position (Phred+33 encoded).
        # mean_p = (10^(-1) + 10^(-2)) / 2 = 0.055
        q10 = UInt8[10 + 33, 10 + 33, 10 + 33]
        q20 = UInt8[20 + 33, 20 + 33, 20 + 33]
        mean_p = Mycelia.Rhizomorph.collapse_molecule_qualities([q10, q20])
        Test.@test length(mean_p) == 3
        Test.@test all(p -> isapprox(p, 0.055; rtol = 1e-12), mean_p)

        # Identical entries: mean_p equals single-entry p
        single = Mycelia.Rhizomorph.collapse_molecule_qualities([q10])
        cloned = Mycelia.Rhizomorph.collapse_molecule_qualities([q10, q10])
        Test.@test single == cloned
    end

    Test.@testset ":quality joint path stays bit-identical (UInt8 sum)" begin
        # Pin the pre-P1 UInt8 independence path. Additive posterior APIs must not
        # change get_vertex_joint_quality / combine_phred_scores outputs.
        seq = "ATGC"
        qual1 = [20, 30, 25, 15]
        qual_str1 = String([Char(q + 33) for q in qual1])
        record1 = FASTX.FASTQ.Record("read_001", seq, qual_str1)

        qual2 = [30, 20, 25, 35]
        qual_str2 = String([Char(q + 33) for q in qual2])
        record2 = FASTX.FASTQ.Record("read_002", seq, qual_str2)

        graph = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(
            [record1, record2], 3; dataset_id = "test"
        )
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.get_vertex_data(graph, kmer)

        joint_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, "test")
        Test.@test joint_qual isa Vector{UInt8}
        Test.@test joint_qual == UInt8[50, 50, 50]
        Test.@test Mycelia.Rhizomorph.combine_phred_scores(UInt8[60, 60, 60, 60, 60]) ==
                  UInt8(255)
    end

    Test.@testset "get_vertex_posterior_log_odds - finite, monotonic under constant Q" begin
        seq = "ATGC"
        qual = [40, 40, 40, 40]
        qual_str = String([Char(q + 33) for q in qual])

        records_1 = [FASTX.FASTQ.Record("mol_001", seq, qual_str)]
        records_2 = [
            FASTX.FASTQ.Record("mol_001", seq, qual_str),
            FASTX.FASTQ.Record("mol_002", seq, qual_str),
        ]
        records_3 = [
            FASTX.FASTQ.Record("mol_001", seq, qual_str),
            FASTX.FASTQ.Record("mol_002", seq, qual_str),
            FASTX.FASTQ.Record("mol_003", seq, qual_str),
        ]

        function lambda_for(records)
            g = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(
                records, 3; dataset_id = "test"
            )
            vd = Mycelia.Rhizomorph.get_vertex_data(g, Kmers.DNAKmer{3}("ATG"))
            return Mycelia.Rhizomorph.get_vertex_posterior_log_odds(vd, "test")
        end

        Λ1 = lambda_for(records_1)
        Λ2 = lambda_for(records_2)
        Λ3 = lambda_for(records_3)

        Test.@test Λ1 isa Float64
        Test.@test isfinite(Λ1) && isfinite(Λ2) && isfinite(Λ3)
        Test.@test !isnan(Λ1) && !isnan(Λ2) && !isnan(Λ3)
        Test.@test Λ1 < Λ2 < Λ3

        scores = let
            g = Mycelia.Rhizomorph.build_qualmer_graph_singlestrand(
                records_2, 3; dataset_id = "test"
            )
            vd = Mycelia.Rhizomorph.get_vertex_data(g, Kmers.DNAKmer{3}("ATG"))
            Mycelia.Rhizomorph.get_vertex_posterior_log_odds_scores(vd, "test")
        end
        Test.@test scores isa Vector{Float64}
        Test.@test length(scores) == 3
        Test.@test all(isfinite, scores)
        Test.@test minimum(scores) == Λ2
    end

    Test.@testset "duplicate-read invariance (same observation_id clone)" begin
        # Cloning a second entry under the SAME observation_id (distinct position so
        # Set keeps both) with identical qualities must leave Λ bit-identical.
        kmer = Kmers.DNAKmer{3}("ATG")
        qual = UInt8[40 + 33, 40 + 33, 40 + 33]

        base = Mycelia.Rhizomorph.QualmerVertexData(kmer)
        Mycelia.Rhizomorph.add_evidence!(
            base,
            "test",
            "read_dup",
            Mycelia.Rhizomorph.QualityEvidenceEntry(1, Mycelia.Rhizomorph.Forward, copy(qual)),
        )
        Λ_base = Mycelia.Rhizomorph.get_vertex_posterior_log_odds(base, "test")

        cloned = Mycelia.Rhizomorph.QualmerVertexData(kmer)
        Mycelia.Rhizomorph.add_evidence!(
            cloned,
            "test",
            "read_dup",
            Mycelia.Rhizomorph.QualityEvidenceEntry(1, Mycelia.Rhizomorph.Forward, copy(qual)),
        )
        Mycelia.Rhizomorph.add_evidence!(
            cloned,
            "test",
            "read_dup",
            Mycelia.Rhizomorph.QualityEvidenceEntry(2, Mycelia.Rhizomorph.Forward, copy(qual)),
        )
        Λ_cloned = Mycelia.Rhizomorph.get_vertex_posterior_log_odds(cloned, "test")

        Test.@test Λ_base isa Float64
        Test.@test Λ_cloned isa Float64
        Test.@test Λ_base === Λ_cloned  # bit-identical Float64

        # :quality path still compounds every entry (unchanged semantics).
        joint_base = Mycelia.Rhizomorph.get_vertex_joint_quality(base, "test")
        joint_cloned = Mycelia.Rhizomorph.get_vertex_joint_quality(cloned, "test")
        Test.@test joint_base == UInt8[40, 40, 40]
        Test.@test joint_cloned == UInt8[80, 80, 80]
    end

    Test.@testset "PE fragment key: /1 and /2 collapse to one molecule" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        qual = UInt8[30 + 33, 30 + 33, 30 + 33]

        pe = Mycelia.Rhizomorph.QualmerVertexData(kmer)
        Mycelia.Rhizomorph.add_evidence!(
            pe,
            "test",
            "fragX/1",
            Mycelia.Rhizomorph.QualityEvidenceEntry(1, Mycelia.Rhizomorph.Forward, copy(qual)),
        )
        Mycelia.Rhizomorph.add_evidence!(
            pe,
            "test",
            "fragX/2",
            Mycelia.Rhizomorph.QualityEvidenceEntry(1, Mycelia.Rhizomorph.Forward, copy(qual)),
        )

        single = Mycelia.Rhizomorph.QualmerVertexData(kmer)
        Mycelia.Rhizomorph.add_evidence!(
            single,
            "test",
            "fragX",
            Mycelia.Rhizomorph.QualityEvidenceEntry(1, Mycelia.Rhizomorph.Forward, copy(qual)),
        )

        Λ_pe = Mycelia.Rhizomorph.get_vertex_posterior_log_odds(pe, "test")
        Λ_single = Mycelia.Rhizomorph.get_vertex_posterior_log_odds(single, "test")
        Test.@test Λ_pe === Λ_single

        # Distinct fragment ids still count as two molecules → higher Λ.
        two_mol = Mycelia.Rhizomorph.QualmerVertexData(kmer)
        Mycelia.Rhizomorph.add_evidence!(
            two_mol,
            "test",
            "fragA/1",
            Mycelia.Rhizomorph.QualityEvidenceEntry(1, Mycelia.Rhizomorph.Forward, copy(qual)),
        )
        Mycelia.Rhizomorph.add_evidence!(
            two_mol,
            "test",
            "fragB/2",
            Mycelia.Rhizomorph.QualityEvidenceEntry(1, Mycelia.Rhizomorph.Forward, copy(qual)),
        )
        Λ_two = Mycelia.Rhizomorph.get_vertex_posterior_log_odds(two_mol, "test")
        Test.@test Λ_two > Λ_pe
    end
end

println("✓ Quality functions tests completed")
