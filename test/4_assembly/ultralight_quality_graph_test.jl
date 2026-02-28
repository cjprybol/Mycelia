import Test
import Mycelia

Test.@testset "Ultralight Quality Graph Mode" begin
    # Helper to create FASTQ records with known quality scores
    function make_fastq_record(id, seq, quals)
        Mycelia.FASTX.FASTQ.Record(id, seq, quals)
    end

    Test.@testset "Ultralight Quality K-mer Graph - basic construction" begin
        # Create FASTQ records with quality scores
        seq = "ATGATGATG"
        quals = UInt8[63, 63, 63, 63, 63, 63, 63, 63, 63]  # all Q30 (Phred+33)
        records = [make_fastq_record("seq1", seq, quals)]

        full_graph = Mycelia.Rhizomorph.build_qualmer_graph(records, 3)
        ulq_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight_quality)

        # Same topology
        Test.@test Mycelia.Graphs.nv(full_graph.graph) == Mycelia.Graphs.nv(ulq_graph.graph)
        Test.@test Mycelia.Graphs.ne(full_graph.graph) == Mycelia.Graphs.ne(ulq_graph.graph)

        # Same labels
        full_labels = sort(collect(Mycelia.MetaGraphsNext.labels(full_graph)); by = string)
        ulq_labels = sort(collect(Mycelia.MetaGraphsNext.labels(ulq_graph)); by = string)
        Test.@test length(full_labels) == length(ulq_labels)
        for (fl, ul) in zip(full_labels, ulq_labels)
            Test.@test string(fl) == string(ul)
        end

        # count_evidence_entries matches
        for label in full_labels
            full_entries = Mycelia.Rhizomorph.count_evidence_entries(full_graph[label])
            ulq_entries = Mycelia.Rhizomorph.count_evidence_entries(ulq_graph[label])
            Test.@test full_entries == ulq_entries
        end
    end

    Test.@testset "Ultralight Quality - joint quality matches full mode" begin
        # Two reads covering the same k-mer with known quality
        seq1 = "ATGATG"
        seq2 = "ATGATG"
        quals1 = UInt8[43, 48, 45, 43, 48, 45]  # Q10, Q15, Q12, Q10, Q15, Q12
        quals2 = UInt8[53, 53, 53, 53, 53, 53]  # Q20, Q20, Q20, Q20, Q20, Q20
        records = [
            make_fastq_record("read1", seq1, quals1),
            make_fastq_record("read2", seq2, quals2)
        ]

        full_graph = Mycelia.Rhizomorph.build_qualmer_graph(records, 3)
        ulq_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight_quality)

        # Compare joint quality scores for each k-mer
        for label in Mycelia.MetaGraphsNext.labels(full_graph)
            full_jq = Mycelia.Rhizomorph.get_vertex_joint_quality(
                full_graph[label], "dataset_01")
            ulq_jq = Mycelia.Rhizomorph.get_vertex_joint_quality(
                ulq_graph[label], "dataset_01")

            if !isnothing(full_jq)
                Test.@test !isnothing(ulq_jq)
                Test.@test length(full_jq) == length(ulq_jq)
                # Joint quality should match (incremental == batch)
                for i in eachindex(full_jq)
                    Test.@test full_jq[i] == ulq_jq[i]
                end
            end
        end
    end

    Test.@testset "Ultralight Quality - no observation ID tracking" begin
        seq = "ATGATG"
        quals = UInt8[63, 63, 63, 63, 63, 63]
        records = [make_fastq_record("read1", seq, quals)]

        ulq_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight_quality)

        first_label = first(Mycelia.MetaGraphsNext.labels(ulq_graph))
        vdata = ulq_graph[first_label]

        # No observation ID tracking
        Test.@test isnothing(Mycelia.Rhizomorph.get_all_observation_ids(
            vdata, "dataset_01"))

        # But dataset IDs are tracked
        Test.@test "dataset_01" in Mycelia.Rhizomorph.get_all_dataset_ids(vdata)
    end

    Test.@testset "Ultralight Quality - edge data" begin
        e = Mycelia.Rhizomorph.UltralightQualityEdgeData()
        Test.@test e.overlap_length == 0
        Test.@test Mycelia.Rhizomorph.count_evidence(e) == 0
        Test.@test isempty(e.from_joint_quality)
        Test.@test isempty(e.to_joint_quality)
    end

    Test.@testset "Ultralight Quality - vertex data direct" begin
        kmer = Mycelia.Kmers.DNAKmer{3}("ATG")
        v = Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer)

        Test.@test v.total_count == 0
        Test.@test isempty(v.joint_quality)
        Test.@test isempty(v.dataset_joint_quality)

        # Add quality evidence
        Mycelia.Rhizomorph.add_evidence!(
            v, "ds1", "obs1",
            Mycelia.Rhizomorph.QualityEvidenceEntry(
                1, Mycelia.Rhizomorph.Forward, UInt8[43, 48, 45]))  # Q10, Q15, Q12

        Test.@test v.total_count == 1
        Test.@test length(v.joint_quality) == 3
        Test.@test v.joint_quality == UInt8[10, 15, 12]

        # Add second observation — quality should accumulate
        Mycelia.Rhizomorph.add_evidence!(
            v, "ds1", "obs2",
            Mycelia.Rhizomorph.QualityEvidenceEntry(
                1, Mycelia.Rhizomorph.Forward, UInt8[53, 53, 53]))  # Q20, Q20, Q20

        Test.@test v.total_count == 2
        Test.@test v.joint_quality == UInt8[30, 35, 32]  # Q10+Q20=Q30, Q15+Q20=Q35, etc.

        # Dataset-specific quality
        ds_qual = Mycelia.Rhizomorph.get_vertex_joint_quality(v, "ds1")
        Test.@test !isnothing(ds_qual)
        Test.@test ds_qual == UInt8[30, 35, 32]
    end

    Test.@testset "Ultralight Quality via build_qualmer_graph" begin
        seq = "ATGATG"
        quals = UInt8[63, 63, 63, 63, 63, 63]
        records = [make_fastq_record("read1", seq, quals)]

        ulq_graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; memory_profile = :ultralight_quality)

        Test.@test Mycelia.Graphs.nv(ulq_graph.graph) > 0
    end
end

Test.@testset "Lightweight Quality Graph Mode" begin
    function make_fastq_record(id, seq, quals)
        Mycelia.FASTX.FASTQ.Record(id, seq, quals)
    end

    Test.@testset "Lightweight Quality K-mer Graph - basic construction" begin
        seq = "ATGATGATG"
        quals = UInt8[63, 63, 63, 63, 63, 63, 63, 63, 63]
        records = [make_fastq_record("seq1", seq, quals)]

        full_graph = Mycelia.Rhizomorph.build_qualmer_graph(records, 3)
        lwq_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :lightweight_quality)

        # Same topology
        Test.@test Mycelia.Graphs.nv(full_graph.graph) == Mycelia.Graphs.nv(lwq_graph.graph)
        Test.@test Mycelia.Graphs.ne(full_graph.graph) == Mycelia.Graphs.ne(lwq_graph.graph)

        # count_evidence_entries matches
        for label in Mycelia.MetaGraphsNext.labels(full_graph)
            full_entries = Mycelia.Rhizomorph.count_evidence_entries(full_graph[label])
            lwq_entries = Mycelia.Rhizomorph.count_evidence_entries(lwq_graph[label])
            Test.@test full_entries == lwq_entries
        end
    end

    Test.@testset "Lightweight Quality - has observation ID tracking" begin
        seq = "ATGATG"
        quals = UInt8[63, 63, 63, 63, 63, 63]
        records = [make_fastq_record("read1", seq, quals)]

        lwq_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :lightweight_quality)

        first_label = first(Mycelia.MetaGraphsNext.labels(lwq_graph))
        vdata = lwq_graph[first_label]

        # Observation IDs ARE tracked (unlike ultralight quality)
        obs_ids = Mycelia.Rhizomorph.get_all_observation_ids(vdata, "dataset_01")
        Test.@test !isnothing(obs_ids)
        Test.@test "read1" in obs_ids
    end

    Test.@testset "Lightweight Quality - joint quality matches ultralight quality" begin
        seq1 = "ATGATG"
        seq2 = "ATGATG"
        quals1 = UInt8[43, 48, 45, 43, 48, 45]
        quals2 = UInt8[53, 53, 53, 53, 53, 53]
        records = [
            make_fastq_record("read1", seq1, quals1),
            make_fastq_record("read2", seq2, quals2)
        ]

        ulq_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight_quality)
        lwq_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :lightweight_quality)

        # Joint quality should be identical between ultralight_quality and lightweight_quality
        for label in Mycelia.MetaGraphsNext.labels(ulq_graph)
            ulq_jq = Mycelia.Rhizomorph.get_vertex_joint_quality(
                ulq_graph[label], "dataset_01")
            lwq_jq = Mycelia.Rhizomorph.get_vertex_joint_quality(
                lwq_graph[label], "dataset_01")

            if !isnothing(ulq_jq)
                Test.@test !isnothing(lwq_jq)
                Test.@test ulq_jq == lwq_jq
            end
        end
    end

    Test.@testset "Lightweight Quality - count_total_observations uses unique obs" begin
        # Same read visiting ATG twice (positions 1 and 4 in ATGATG)
        seq = "ATGATG"
        quals = UInt8[63, 63, 63, 63, 63, 63]
        records = [make_fastq_record("read1", seq, quals)]

        lwq_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :lightweight_quality)

        # ATG appears at positions 1 and 4
        # count_evidence_entries = 2 (two add_evidence! calls)
        # count_total_observations = 1 (unique obs ID "read1")
        for label in Mycelia.MetaGraphsNext.labels(lwq_graph)
            vdata = lwq_graph[label]
            entries = Mycelia.Rhizomorph.count_evidence_entries(vdata)
            obs = Mycelia.Rhizomorph.count_total_observations(vdata)
            Test.@test obs <= entries
        end
    end

    Test.@testset "Lightweight Quality via build_qualmer_graph" begin
        seq = "ATGATG"
        quals = UInt8[63, 63, 63, 63, 63, 63]
        records = [make_fastq_record("read1", seq, quals)]

        lwq_graph = Mycelia.Rhizomorph.build_qualmer_graph(
            records, 3; memory_profile = :lightweight_quality)

        Test.@test Mycelia.Graphs.nv(lwq_graph.graph) > 0
    end

    Test.@testset "Lightweight Quality vertex data direct" begin
        kmer = Mycelia.Kmers.DNAKmer{3}("ATG")
        v = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)

        Test.@test v.total_count == 0
        Test.@test isempty(v.joint_quality)
        Test.@test isempty(v.dataset_observations)

        Mycelia.Rhizomorph.add_evidence!(
            v, "ds1", "obs1",
            Mycelia.Rhizomorph.QualityEvidenceEntry(
                1, Mycelia.Rhizomorph.Forward, UInt8[43, 48, 45]))

        Test.@test v.total_count == 1
        Test.@test v.joint_quality == UInt8[10, 15, 12]
        Test.@test "obs1" in v.dataset_observations["ds1"]
    end

    Test.@testset "All four profiles produce same topology" begin
        seq = "ATGATGATG"
        quals = UInt8[63, 63, 63, 63, 63, 63, 63, 63, 63]
        records = [make_fastq_record("seq1", seq, quals)]

        full = Mycelia.Rhizomorph.build_qualmer_graph(records, 3)
        lw = Mycelia.Rhizomorph.build_kmer_graph(records, 3; memory_profile = :lightweight)
        ul = Mycelia.Rhizomorph.build_kmer_graph(records, 3; memory_profile = :ultralight)
        ulq = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight_quality)
        lwq = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :lightweight_quality)

        nv = Mycelia.Graphs.nv(full.graph)
        ne = Mycelia.Graphs.ne(full.graph)

        Test.@test Mycelia.Graphs.nv(lw.graph) == nv
        Test.@test Mycelia.Graphs.nv(ul.graph) == nv
        Test.@test Mycelia.Graphs.nv(ulq.graph) == nv
        Test.@test Mycelia.Graphs.nv(lwq.graph) == nv

        Test.@test Mycelia.Graphs.ne(lw.graph) == ne
        Test.@test Mycelia.Graphs.ne(ul.graph) == ne
        Test.@test Mycelia.Graphs.ne(ulq.graph) == ne
        Test.@test Mycelia.Graphs.ne(lwq.graph) == ne
    end
end
