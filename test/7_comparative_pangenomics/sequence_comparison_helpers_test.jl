import Test
import Mycelia

Test.@testset "Sequence Comparison Helpers" begin
    Test.@testset "Entropy and richness" begin
        seq = "ATAT"
        Test.@test Mycelia.shannon_entropy(seq) > 0.0
        Test.@test Mycelia.renyi_entropy(seq; α=2) > 0.0
        Test.@test Mycelia.kmer_richness(seq, 2; normalize=false) == 2

        profile, summary = Mycelia.linguistic_complexity(seq; kmax=2)
        Test.@test length(profile) == 2
        Test.@test summary <= 1.0
    end

    Test.@testset "Mash distance helper" begin
        distance = Mycelia.mash_distance_from_jaccard(0.5, 21)
        Test.@test distance > 0.0
    end

    Test.@testset "Hash and file utilities" begin
        Test.@test Mycelia.seq2sha256("ATCG") == Mycelia.seq2sha256("ATCG")

        mktempdir() do dir
            fasta_path = joinpath(dir, "test.fasta")
            open(fasta_path, "w") do io
                println(io, ">seq1")
                println(io, "ATCG")
                println(io, ">seq2")
                println(io, "ATGC")
            end
            Test.@test Mycelia.count_fasta_records(fasta_path) == 2
            Test.@test !isempty(Mycelia.genome_pair_id(fasta_path, fasta_path))
        end
    end

    Test.@testset "Circular sequence helpers" begin
        Test.@test Mycelia.reverse_complement_ascii("ACGT") == "ACGT"
        Test.@test Mycelia.booth_min_rotation_index(Vector{UInt8}(codeunits("baba"))) == 2
        Test.@test Mycelia.rotate_bytes(Vector{UInt8}(codeunits("abcd")), 3) == "cdab"

        canonical, orientation, _ = Mycelia.canonical_circular_sequence("AAAA")
        Test.@test canonical == "AAAA"
        Test.@test orientation == :forward

        Test.@test Mycelia.header_says_circular("plasmid circular genome")
        Test.@test Mycelia.has_terminal_overlap("ATGCATGC"; overlap_bp=4)
    end
end
