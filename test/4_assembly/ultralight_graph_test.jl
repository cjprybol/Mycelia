import Test
import Mycelia

Test.@testset "Ultralight Graph Mode" begin
    Test.@testset "Ultralight N-gram Graph - basic construction" begin
        texts = ["ABCDE", "BCDEF", "CDEFG"]

        full_graph = Mycelia.Rhizomorph.build_ngram_graph(texts, 3)
        ul_graph = Mycelia.Rhizomorph.build_ngram_graph(
            texts, 3; memory_profile = :ultralight)

        # Same topology
        Test.@test Mycelia.Graphs.nv(full_graph.graph) == Mycelia.Graphs.nv(ul_graph.graph)
        Test.@test Mycelia.Graphs.ne(full_graph.graph) == Mycelia.Graphs.ne(ul_graph.graph)

        # Same vertex labels
        full_labels = sort(collect(Mycelia.MetaGraphsNext.labels(full_graph)))
        ul_labels = sort(collect(Mycelia.MetaGraphsNext.labels(ul_graph)))
        Test.@test full_labels == ul_labels

        # count_evidence_entries should match
        for label in full_labels
            full_count = Mycelia.Rhizomorph.count_evidence_entries(full_graph[label])
            ul_count = Mycelia.Rhizomorph.count_evidence_entries(ul_graph[label])
            Test.@test full_count == ul_count
        end

        # Edge counts should match
        for (src, dst) in Mycelia.MetaGraphsNext.edge_labels(full_graph)
            full_count = Mycelia.Rhizomorph.count_evidence(full_graph[src, dst])
            ul_count = Mycelia.Rhizomorph.count_evidence(ul_graph[src, dst])
            Test.@test full_count == ul_count
        end
    end

    Test.@testset "Ultralight N-gram Graph - multi-dataset" begin
        texts1 = ["ABCDE"]
        texts2 = ["BCDEF"]

        ul_graph = Mycelia.Rhizomorph.build_ngram_graph(
            texts1, 3; dataset_id = "ds1", memory_profile = :ultralight)
        Mycelia.Rhizomorph.add_observations_to_ngram_graph_ultralight!(
            ul_graph, texts2, 3; dataset_id = "ds2")

        # BCD vertex should have evidence from both datasets
        ul_bcd_ds = sort(Mycelia.Rhizomorph.get_all_dataset_ids(ul_graph["BCD"]))
        Test.@test "ds1" in ul_bcd_ds
        Test.@test "ds2" in ul_bcd_ds
    end

    Test.@testset "Ultralight N-gram Graph - statistics compatibility" begin
        texts = ["hello world", "world peace"]
        ul_graph = Mycelia.Rhizomorph.build_ngram_graph(
            texts, 3; memory_profile = :ultralight)

        stats = Mycelia.Rhizomorph.get_ngram_statistics(ul_graph)
        Test.@test stats[:num_vertices] > 0
        Test.@test stats[:num_edges] > 0
        Test.@test stats[:n] == 3
        Test.@test !isnothing(stats[:most_common_ngram])
    end

    Test.@testset "Ultralight N-gram Graph - analysis functions" begin
        texts = ["ABCABC", "BCABCA"]
        ul_graph = Mycelia.Rhizomorph.build_ngram_graph(
            texts, 3; memory_profile = :ultralight)

        common = Mycelia.Rhizomorph.find_high_coverage_ngrams(ul_graph, 2)
        Test.@test !isempty(common)

        unique_ngrams = Mycelia.Rhizomorph.find_unique_ngrams(ul_graph, "dataset_01")
        Test.@test isa(unique_ngrams, Vector{String})

        shared = Mycelia.Rhizomorph.find_shared_ngrams(ul_graph, ["dataset_01"])
        Test.@test !isempty(shared)
    end

    Test.@testset "Ultralight evidence function interface" begin
        v = Mycelia.Rhizomorph.UltralightStringVertexData("ABC")
        Test.@test Mycelia.Rhizomorph.count_total_observations(v) == 0
        Test.@test Mycelia.Rhizomorph.count_evidence_entries(v) == 0
        Test.@test isempty(Mycelia.Rhizomorph.get_all_dataset_ids(v))

        # Add evidence and verify counts
        Mycelia.Rhizomorph.add_evidence!(
            v, "ds1", "obs1",
            Mycelia.Rhizomorph.EvidenceEntry(1, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(
            v, "ds1", "obs2",
            Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(
            v, "ds2", "obs3",
            Mycelia.Rhizomorph.EvidenceEntry(3, Mycelia.Rhizomorph.Forward))

        # count_total_observations returns total_count (NOT unique obs)
        Test.@test Mycelia.Rhizomorph.count_total_observations(v) == 3
        Test.@test Mycelia.Rhizomorph.count_evidence_entries(v) == 3
        Test.@test Mycelia.Rhizomorph.count_dataset_observations(v, "ds1") == 2
        Test.@test Mycelia.Rhizomorph.count_dataset_observations(v, "ds2") == 1
        Test.@test sort(Mycelia.Rhizomorph.get_all_dataset_ids(v)) == ["ds1", "ds2"]

        # Observation IDs NOT tracked in ultralight mode
        Test.@test isnothing(Mycelia.Rhizomorph.get_all_observation_ids(v, "ds1"))

        # Positional evidence queries return nothing
        Test.@test isnothing(Mycelia.Rhizomorph.get_dataset_evidence(v, "ds1"))
        Test.@test isnothing(Mycelia.Rhizomorph.get_observation_evidence(v, "ds1", "obs1"))

        # Strand defaults
        Test.@test Mycelia.Rhizomorph.first_evidence_strand(v) == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "Ultralight edge data interface" begin
        e = Mycelia.Rhizomorph.UltralightEdgeData(2)
        Test.@test e.overlap_length == 2
        Test.@test Mycelia.Rhizomorph.count_evidence(e) == 0

        Mycelia.Rhizomorph.add_evidence!(
            e, "ds1", "obs1",
            Mycelia.Rhizomorph.EdgeEvidenceEntry(1, 2, Mycelia.Rhizomorph.Forward))

        Test.@test Mycelia.Rhizomorph.count_evidence(e) == 1
        Test.@test Mycelia.Rhizomorph.compute_edge_weight(e) == 1.0
        Test.@test Mycelia.Rhizomorph.compute_edge_coverage(e) == 1
    end

    Test.@testset "Ultralight K-mer Graph - basic construction" begin
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG"),
            Mycelia.FASTX.FASTA.Record("seq2", "ATGATCATG")
        ]

        full_graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3)
        ul_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight)

        # Same topology
        Test.@test Mycelia.Graphs.nv(full_graph.graph) == Mycelia.Graphs.nv(ul_graph.graph)
        Test.@test Mycelia.Graphs.ne(full_graph.graph) == Mycelia.Graphs.ne(ul_graph.graph)

        # Same labels
        full_labels = sort(collect(Mycelia.MetaGraphsNext.labels(full_graph));
            by = string)
        ul_labels = sort(collect(Mycelia.MetaGraphsNext.labels(ul_graph));
            by = string)
        Test.@test length(full_labels) == length(ul_labels)
        for (fl, ul) in zip(full_labels, ul_labels)
            Test.@test string(fl) == string(ul)
        end

        # count_evidence_entries matches
        for label in full_labels
            full_entries = Mycelia.Rhizomorph.count_evidence_entries(full_graph[label])
            ul_entries = Mycelia.Rhizomorph.count_evidence_entries(ul_graph[label])
            Test.@test full_entries == ul_entries
        end

        # Edge evidence counts match
        for (src, dst) in Mycelia.MetaGraphsNext.edge_labels(full_graph)
            full_count = Mycelia.Rhizomorph.count_evidence(full_graph[src, dst])
            ul_count = Mycelia.Rhizomorph.count_evidence(ul_graph[src, dst])
            Test.@test full_count == ul_count
        end
    end

    Test.@testset "Ultralight K-mer Graph - statistics" begin
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG")
        ]

        ul_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight)
        stats = Mycelia.Rhizomorph.get_kmer_statistics(ul_graph)

        Test.@test stats[:num_vertices] > 0
        Test.@test stats[:num_edges] > 0
        Test.@test stats[:k] == 3
    end

    Test.@testset "Ultralight mode errors on non-singlestrand" begin
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG")
        ]

        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight, mode = :doublestrand)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight, mode = :canonical)
    end

    Test.@testset "Ultralight backward compat - lightweight=true still works" begin
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG")
        ]

        lw_graph = Mycelia.Rhizomorph.build_kmer_graph(records, 3; lightweight = true)
        Test.@test Mycelia.Graphs.nv(lw_graph.graph) > 0

        # Verify it's actually a lightweight graph (has observation tracking)
        first_label = first(Mycelia.MetaGraphsNext.labels(lw_graph))
        vdata = lw_graph[first_label]
        Test.@test vdata isa Mycelia.Rhizomorph.LightweightKmerVertexData

        obs_ids = Mycelia.Rhizomorph.get_all_observation_ids(vdata, "dataset_01")
        Test.@test !isnothing(obs_ids)
    end

    Test.@testset "Ultralight path-finding compatibility" begin
        texts = ["ABCDE", "BCDEF"]
        ul_graph = Mycelia.Rhizomorph.build_ngram_graph(
            texts, 3; memory_profile = :ultralight)

        weighted = Mycelia.Rhizomorph.weighted_graph_from_rhizomorph(ul_graph)
        Test.@test Mycelia.Graphs.nv(weighted.graph) == Mycelia.Graphs.nv(ul_graph.graph)
        Test.@test Mycelia.Graphs.ne(weighted.graph) >= Mycelia.Graphs.ne(ul_graph.graph)
    end

    Test.@testset "Ultralight type mutable struct behavior" begin
        v = Mycelia.Rhizomorph.UltralightStringVertexData("test")
        Test.@test v.total_count == 0

        v.total_count = 42
        Test.@test v.total_count == 42

        e = Mycelia.Rhizomorph.UltralightEdgeData()
        Test.@test e.total_count == 0
        e.total_count = 10
        Test.@test e.total_count == 10
    end

    Test.@testset "Ultralight N-gram Graph - from file" begin
        mktempdir() do dir
            path = joinpath(dir, "doc.txt")
            open(path, "w") do io
                println(io, "ABCD")
                println(io, "BCDA")
            end

            full_graph = Mycelia.Rhizomorph.build_ngram_graph_from_file(path, 2)
            ul_graph = Mycelia.Rhizomorph.build_ngram_graph_from_file(
                path, 2; memory_profile = :ultralight)

            Test.@test Mycelia.Graphs.nv(full_graph.graph) ==
                       Mycelia.Graphs.nv(ul_graph.graph)
            Test.@test Mycelia.Graphs.ne(full_graph.graph) ==
                       Mycelia.Graphs.ne(ul_graph.graph)
        end
    end

    Test.@testset "Ultralight N-gram Graph - from multiple files" begin
        mktempdir() do dir
            path1 = joinpath(dir, "doc1.txt")
            path2 = joinpath(dir, "doc2.txt")
            open(path1, "w") do io
                println(io, "ABCD")
            end
            open(path2, "w") do io
                println(io, "BCXY")
            end

            full_graph = Mycelia.Rhizomorph.build_ngram_graph_from_files(
                [path1, path2], 2)
            ul_graph = Mycelia.Rhizomorph.build_ngram_graph_from_files(
                [path1, path2], 2; memory_profile = :ultralight)

            Test.@test Mycelia.Graphs.nv(full_graph.graph) ==
                       Mycelia.Graphs.nv(ul_graph.graph)
            Test.@test Mycelia.Graphs.ne(full_graph.graph) ==
                       Mycelia.Graphs.ne(ul_graph.graph)
        end
    end

    Test.@testset "Ultralight K-mer Graph - from file" begin
        mktempdir() do dir
            path = joinpath(dir, "test.fasta")
            open(path, "w") do io
                println(io, ">seq1")
                println(io, "ATGATGATG")
                println(io, ">seq2")
                println(io, "ATGATCATG")
            end

            full_graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(path, 3)
            ul_graph = Mycelia.Rhizomorph.build_kmer_graph_from_file(
                path, 3; memory_profile = :ultralight)

            Test.@test Mycelia.Graphs.nv(full_graph.graph) ==
                       Mycelia.Graphs.nv(ul_graph.graph)
            Test.@test Mycelia.Graphs.ne(full_graph.graph) ==
                       Mycelia.Graphs.ne(ul_graph.graph)
        end
    end

    Test.@testset "Ultralight K-mer Graph - from multiple files" begin
        mktempdir() do dir
            path1 = joinpath(dir, "sample1.fasta")
            path2 = joinpath(dir, "sample2.fasta")
            open(path1, "w") do io
                println(io, ">seq1")
                println(io, "ATGATGATG")
            end
            open(path2, "w") do io
                println(io, ">seq2")
                println(io, "ATGATCATG")
            end

            full_graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [path1, path2], 3)
            ul_graph = Mycelia.Rhizomorph.build_kmer_graph_from_files(
                [path1, path2], 3; memory_profile = :ultralight)

            Test.@test Mycelia.Graphs.nv(full_graph.graph) ==
                       Mycelia.Graphs.nv(ul_graph.graph)
            Test.@test Mycelia.Graphs.ne(full_graph.graph) ==
                       Mycelia.Graphs.ne(ul_graph.graph)
        end
    end

    Test.@testset "Ultralight vs Lightweight semantic difference" begin
        # When a k-mer appears twice in the same read, count_total_observations
        # diverges: lightweight counts unique obs IDs, ultralight counts total entries
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATG")
        ]

        lw_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :lightweight)
        ul_graph = Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight)

        # ATG appears at positions 1 and 4 in ATGATG (from same read seq1)
        # Lightweight: count_total_observations = 1 (unique obs ID "seq1")
        # Ultralight: count_total_observations = 2 (total add_evidence! calls)
        # count_evidence_entries should be 2 in both
        for label in Mycelia.MetaGraphsNext.labels(lw_graph)
            lw_entries = Mycelia.Rhizomorph.count_evidence_entries(lw_graph[label])
            ul_entries = Mycelia.Rhizomorph.count_evidence_entries(ul_graph[label])
            Test.@test lw_entries == ul_entries
        end
    end

    Test.@testset "Quality profiles error on n-gram graphs" begin
        texts = ["ABCDE"]
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_ngram_graph(
            texts, 3; memory_profile = :ultralight_quality)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_ngram_graph(
            texts, 3; memory_profile = :lightweight_quality)
    end

    Test.@testset "Quality profiles error on FASTA input" begin
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG")
        ]
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :ultralight_quality)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :lightweight_quality)
    end

    Test.@testset "Invalid memory_profile errors" begin
        records = [
            Mycelia.FASTX.FASTA.Record("seq1", "ATGATGATG")
        ]
        Test.@test_throws ErrorException Mycelia.Rhizomorph.build_kmer_graph(
            records, 3; memory_profile = :invalid)
    end
end
