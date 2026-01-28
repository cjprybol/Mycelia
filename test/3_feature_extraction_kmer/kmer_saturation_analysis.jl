import Test
import Mycelia
import BioSequences
import Kmers
import FASTX

Test.@testset "K-mer saturation analysis" begin
    sequences = [
        BioSequences.LongDNA{4}("ATGCA"),
        BioSequences.LongDNA{4}("ATGGA"),
        BioSequences.LongDNA{4}("ATGTA"),
        BioSequences.LongDNA{4}("ATGAA"),
    ]

    result = Mycelia.analyze_kmer_saturation(
        sequences,
        [3, 4];
        max_sampling_points=3,
        min_sequences=2,
    )

    Test.@test result.optimal_k in [3, 4]
    Test.@test length(result.saturation_levels) == 2
    Test.@test length(result.sampling_points) == length(result.unique_kmer_counts[3])
    Test.@test length(result.sampling_points) == length(result.unique_kmer_counts[4])
    Test.@test all(k -> haskey(result.unique_kmer_counts, k), [3, 4])

    Test.@test_throws ArgumentError Mycelia.analyze_kmer_saturation(sequences, Int[]; min_sequences=2)
    Test.@test_throws ArgumentError Mycelia.analyze_kmer_saturation(sequences, [3]; min_sequences=10)

    kmer_counts = Dict(
        Kmers.DNAKmer{3}(BioSequences.LongDNA{4}("ATG")) => 5,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{4}("TGC")) => 4,
        Kmers.DNAKmer{3}(BioSequences.LongDNA{4}("GCA")) => 1,
    )

    threshold_result = Mycelia.find_optimal_assembly_threshold(
        kmer_counts,
        3;
        min_threshold=1,
        max_threshold=5,
    )

    Test.@test threshold_result.optimal_threshold in threshold_result.thresholds_tested
    Test.@test length(threshold_result.thresholds_tested) == length(threshold_result.connectivity_scores)
    Test.@test length(threshold_result.thresholds_tested) == length(threshold_result.component_counts)
    Test.@test length(threshold_result.thresholds_tested) == length(threshold_result.average_connectivity)
end

Test.@testset "K-mer saturation analysis - sequence normalization" begin
    string_sequences = [
        "ATGCA",
        "ATGGA",
        "ATGTA",
    ]

    string_result = Mycelia.analyze_kmer_saturation(
        string_sequences,
        [3];
        max_sampling_points=2,
        min_sequences=3,
    )

    Test.@test string_result.optimal_k == 3
    Test.@test length(string_result.saturation_levels) == 1
    Test.@test length(string_result.sampling_points) == length(string_result.unique_kmer_counts[3])

    fasta_records = [
        FASTX.FASTA.Record("seq1", "ATGCA"),
        FASTX.FASTA.Record("seq2", "ATGGA"),
        FASTX.FASTA.Record("seq3", "ATGTA"),
    ]

    record_result = Mycelia.analyze_kmer_saturation(
        fasta_records,
        [3];
        max_sampling_points=2,
        min_sequences=3,
    )

    Test.@test record_result.optimal_k == 3
    Test.@test length(record_result.saturation_levels) == 1
    Test.@test length(record_result.sampling_points) == length(record_result.unique_kmer_counts[3])

    mixed_sequences = [
        BioSequences.LongDNA{4}("ATGC"),
        BioSequences.LongRNA{4}("AUGC"),
    ]
    Test.@test_throws ArgumentError Mycelia.analyze_kmer_saturation(
        mixed_sequences,
        [3];
        min_sequences=2,
    )

    Test.@test_throws ArgumentError Mycelia.analyze_kmer_saturation(
        [1, 2, 3],
        [3];
        min_sequences=2,
    )
end

Test.@testset "K-mer saturation analysis - helpers" begin
    x = [1, 2, 3, 4, 5]
    vmax = 50.0
    km = 2.0
    y = [round(Int, (vmax * xi) / (km + xi)) for xi in x]

    params, r_squared = Mycelia._fit_michaelis_menten(x, y)
    Test.@test length(params) == 2
    Test.@test all(isfinite, params)
    Test.@test r_squared > 0.8

    aa_kmers = [
        Kmers.AAKmer{3}(BioSequences.LongAA("ACD")),
        Kmers.AAKmer{3}(BioSequences.LongAA("CDE")),
        Kmers.AAKmer{3}(BioSequences.LongAA("DEF")),
    ]

    adjacency = Mycelia._build_kmer_overlap_graph(aa_kmers, 3)
    Test.@test length(adjacency) == 3
    Test.@test 2 in adjacency[1]
    Test.@test 3 in adjacency[2]
    Test.@test isempty(adjacency[3])
    Test.@test_throws ArgumentError Mycelia._build_kmer_overlap_graph(aa_kmers, 4)
    Test.@test isempty(Mycelia._build_kmer_overlap_graph(Kmers.AAKmer{3}[], 3))

    component_adjacency = Vector{Vector{Int}}([[2], [1], Int[], [5], [4]])
    Test.@test Mycelia._count_connected_components(component_adjacency) == 3
    Test.@test isapprox(Mycelia._calculate_average_connectivity(component_adjacency), 0.8; atol=1e-12)

    Test.@test_throws ArgumentError Mycelia.find_optimal_assembly_threshold(
        Dict{Kmers.DNAKmer{3}, Int}(),
        3,
    )

    kmer_sequences = [
        "AAA", "AAC", "AAG", "AAT", "ACA",
        "ACC", "ACG", "ACT", "AGA", "AGC",
    ]
    kmer_counts = Dict(
        Kmers.DNAKmer{3}(BioSequences.LongDNA{4}(seq)) => 5
        for seq in kmer_sequences
    )
    Test.@test_throws ArgumentError Mycelia.find_optimal_assembly_threshold(
        kmer_counts,
        4;
        min_threshold=1,
        max_threshold=5,
    )
end
