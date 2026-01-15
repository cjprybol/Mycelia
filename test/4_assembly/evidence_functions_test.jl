# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/evidence_functions_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/evidence_functions_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Evidence Helper Functions Test
# Tests for evidence manipulation and query helper functions as specified
# in the rhizomorph graph ecosystem planning document section 1.3.1.
# These functions provide efficient evidence queries and manipulation.

import Test
import Mycelia
import Kmers

Test.@testset "Evidence Helper Functions" begin

    Test.@testset "add_evidence! - Vertex" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Add evidence
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))

        Test.@test haskey(vertex.evidence, "dataset_01")
        Test.@test haskey(vertex.evidence["dataset_01"], "read_001")
    end

    Test.@testset "add_evidence! - Edge" begin
        edge = Mycelia.Rhizomorph.KmerEdgeData()

        # Add evidence
        Mycelia.Rhizomorph.add_evidence!(edge, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))

        Test.@test haskey(edge.evidence, "dataset_01")
        Test.@test haskey(edge.evidence["dataset_01"], "read_001")
    end

    Test.@testset "get_dataset_evidence - Vertex" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Add evidence from multiple datasets
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_02", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))

        # Get evidence from specific dataset
        dataset_01_evidence = Mycelia.Rhizomorph.get_dataset_evidence(vertex, "dataset_01")
        dataset_02_evidence = Mycelia.Rhizomorph.get_dataset_evidence(vertex, "dataset_02")

        Test.@test !isnothing(dataset_01_evidence)
        Test.@test !isnothing(dataset_02_evidence)
        Test.@test haskey(dataset_01_evidence, "read_001")
        Test.@test haskey(dataset_02_evidence, "read_001")
    end

    Test.@testset "get_dataset_evidence - Nonexistent Dataset" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # No evidence added
        result = Mycelia.Rhizomorph.get_dataset_evidence(vertex, "nonexistent")

        Test.@test isnothing(result)
    end

    Test.@testset "get_observation_evidence - Vertex" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Add evidence from multiple observations
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward))

        # Get evidence from specific observation
        read_001_evidence = Mycelia.Rhizomorph.get_observation_evidence(vertex, "dataset_01", "read_001")
        read_002_evidence = Mycelia.Rhizomorph.get_observation_evidence(vertex, "dataset_01", "read_002")

        Test.@test !isnothing(read_001_evidence)
        Test.@test !isnothing(read_002_evidence)
        Test.@test length(read_001_evidence) >= 1
        Test.@test length(read_002_evidence) >= 1
    end

    Test.@testset "get_observation_evidence - Nonexistent Observation" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        result = Mycelia.Rhizomorph.get_observation_evidence(vertex, "dataset_01", "nonexistent")

        Test.@test isnothing(result)
    end

    Test.@testset "get_position_evidence - All Positions" begin
        # Get all evidence at a specific position across all observations

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Add evidence at same position from different observations
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_003",
                             Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward))

        # Get all evidence at position 5
        position_5_evidence = Mycelia.Rhizomorph.get_position_evidence(vertex, "dataset_01", 5)

        Test.@test length(position_5_evidence) == 2  # read_001 and read_002
    end

    Test.@testset "merge_evidence_sets - Union" begin
        # Merge evidence from two observations

        evidence1 = Set{Mycelia.Rhizomorph.EvidenceEntry}([
            Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward),
            Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward)
        ])

        evidence2 = Set{Mycelia.Rhizomorph.EvidenceEntry}([
            Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward),  # Duplicate
            Mycelia.Rhizomorph.EvidenceEntry(15, Mycelia.Rhizomorph.Forward)
        ])

        merged = Mycelia.Rhizomorph.merge_evidence_sets(evidence1, evidence2)

        Test.@test length(merged) == 3  # Union removes duplicate
        Test.@test Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward) in merged
        Test.@test Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward) in merged
        Test.@test Mycelia.Rhizomorph.EvidenceEntry(15, Mycelia.Rhizomorph.Forward) in merged
    end

    Test.@testset "flip_evidence_strand - Reverse Complement" begin
        # Flip strand orientation of evidence (for RC operations)

        entry = Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward)
        flipped = Mycelia.Rhizomorph.flip_evidence_strand(entry)

        Test.@test flipped.position == entry.position
        Test.@test flipped.strand == Mycelia.Rhizomorph.Reverse

        # Flip back
        double_flipped = Mycelia.Rhizomorph.flip_evidence_strand(flipped)
        Test.@test double_flipped.strand == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "flip_evidence_strand - Quality Evidence" begin
        # Flipping quality evidence should reverse quality scores

        entry = Mycelia.Rhizomorph.QualityEvidenceEntry(5, Mycelia.Rhizomorph.Forward, UInt8[30, 35, 40])
        flipped = Mycelia.Rhizomorph.flip_evidence_strand(entry)

        Test.@test flipped.position == entry.position
        Test.@test flipped.strand == Mycelia.Rhizomorph.Reverse
        Test.@test flipped.quality_scores == reverse(entry.quality_scores)
        Test.@test flipped.quality_scores == UInt8[40, 35, 30]
    end

    Test.@testset "count_total_observations - Vertex" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Add evidence from multiple observations
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_02", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(3, Mycelia.Rhizomorph.Forward))

        # Total unique observations across all datasets
        total = Mycelia.Rhizomorph.count_total_observations(vertex)

        Test.@test total == 3  # 3 unique observation IDs
    end

    Test.@testset "count_dataset_observations" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_02", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(3, Mycelia.Rhizomorph.Forward))

        # Count observations in specific dataset
        dataset_01_count = Mycelia.Rhizomorph.count_dataset_observations(vertex, "dataset_01")
        dataset_02_count = Mycelia.Rhizomorph.count_dataset_observations(vertex, "dataset_02")

        Test.@test dataset_01_count == 2
        Test.@test dataset_02_count == 1
    end

    Test.@testset "count_evidence_entries - Total Evidence" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Same observation can have multiple evidence entries
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(50, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))

        # Total evidence entries (not unique observations)
        total_entries = Mycelia.Rhizomorph.count_evidence_entries(vertex)

        Test.@test total_entries == 3  # 3 total evidence entries
    end

    Test.@testset "collect_evidence_entries and strands" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Reverse))

        entries = Mycelia.Rhizomorph.collect_evidence_entries(vertex.evidence)
        Test.@test length(entries) == 2

        strands = Mycelia.Rhizomorph.collect_evidence_strands(vertex.evidence)
        Test.@test Set(strands) == Set([Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse])

        first_strand = Mycelia.Rhizomorph.first_evidence_strand(vertex.evidence)
        Test.@test first_strand in (Mycelia.Rhizomorph.Forward, Mycelia.Rhizomorph.Reverse)
    end

    Test.@testset "filter_evidence_by_strand" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Add evidence from both strands
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Reverse))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_003",
                             Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward))

        # Filter to only Forward strand
        forward_only = Mycelia.Rhizomorph.filter_evidence_by_strand(vertex, Mycelia.Rhizomorph.Forward)

        # Should have 2 observations with Forward evidence
        total_forward = sum(length(obs_evidence)
                           for dataset_evidence in values(forward_only.evidence)
                           for obs_evidence in values(dataset_evidence))

        Test.@test total_forward == 2
    end

    Test.@testset "get_all_dataset_ids" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_02", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_03", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(3, Mycelia.Rhizomorph.Forward))

        dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex)

        Test.@test length(dataset_ids) == 3
        Test.@test "dataset_01" in dataset_ids
        Test.@test "dataset_02" in dataset_ids
        Test.@test "dataset_03" in dataset_ids
    end

    Test.@testset "get_all_observation_ids - For Dataset" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex = Mycelia.Rhizomorph.KmerVertexData(kmer)

        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex, "dataset_01", "read_003",
                             Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward))

        observation_ids = Mycelia.Rhizomorph.get_all_observation_ids(vertex, "dataset_01")

        Test.@test length(observation_ids) == 3
        Test.@test "read_001" in observation_ids
        Test.@test "read_002" in observation_ids
        Test.@test "read_003" in observation_ids
    end
end

Test.@testset "Evidence Functions - Edge Support" begin

    Test.@testset "get_dataset_evidence - Edge" begin
        edge = Mycelia.Rhizomorph.KmerEdgeData()

        Mycelia.Rhizomorph.add_evidence!(edge, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge, "dataset_02", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(8, 9, Mycelia.Rhizomorph.Forward))

        dataset_01_evidence = Mycelia.Rhizomorph.get_dataset_evidence(edge, "dataset_01")

        Test.@test !isnothing(dataset_01_evidence)
        Test.@test haskey(dataset_01_evidence, "read_001")
    end

    Test.@testset "count_total_observations - Edge" begin
        edge = Mycelia.Rhizomorph.KmerEdgeData()

        Mycelia.Rhizomorph.add_evidence!(edge, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(10, 11, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge, "dataset_02", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(3, 4, Mycelia.Rhizomorph.Forward))

        total = Mycelia.Rhizomorph.count_total_observations(edge)

        Test.@test total == 3
    end
end

println("✓ Evidence helper functions tests defined (will fail until implementation)")
