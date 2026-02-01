# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/evidence_structures_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/evidence_structures_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Evidence Structures Test
# Tests for the fundamental evidence tracking data structures defined in the
# rhizomorph graph ecosystem planning document.
# These tests define the CORRECT behavior before implementation.

import Test
import Mycelia

Test.@testset "Evidence Entry Structures" begin
    Test.@testset "EvidenceEntry - Basic Construction" begin
        # EvidenceEntry should store position and strand only
        # observation_id is stored as Dict key, not in the entry itself

        entry = Mycelia.Rhizomorph.EvidenceEntry(42, Mycelia.Rhizomorph.Forward)

        Test.@test entry.position == 42
        Test.@test entry.strand == Mycelia.Rhizomorph.Forward
        Test.@test typeof(entry.position) == Int
        Test.@test typeof(entry.strand) == Mycelia.Rhizomorph.StrandOrientation
    end

    Test.@testset "EvidenceEntry - Different Strands" begin
        forward_entry = Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward)
        reverse_entry = Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Reverse)

        Test.@test forward_entry.strand == Mycelia.Rhizomorph.Forward
        Test.@test reverse_entry.strand == Mycelia.Rhizomorph.Reverse
        Test.@test forward_entry.position == reverse_entry.position
        Test.@test forward_entry != reverse_entry  # Different strands means different entries
    end

    Test.@testset "EvidenceEntry - Equality and Hashing" begin
        entry1 = Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward)
        entry2 = Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward)
        entry3 = Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Reverse)
        entry4 = Mycelia.Rhizomorph.EvidenceEntry(6, Mycelia.Rhizomorph.Forward)

        # Same position and strand should be equal
        Test.@test entry1 == entry2
        Test.@test hash(entry1) == hash(entry2)

        # Different strand should not be equal
        Test.@test entry1 != entry3

        # Different position should not be equal
        Test.@test entry1 != entry4
    end

    Test.@testset "QualityEvidenceEntry - With Quality Scores" begin
        # QualityEvidenceEntry adds quality scores to basic evidence
        quality_scores = UInt8[30, 35, 32, 28, 40]

        entry = Mycelia.Rhizomorph.QualityEvidenceEntry(10, Mycelia.Rhizomorph.Forward, quality_scores)

        Test.@test entry.position == 10
        Test.@test entry.strand == Mycelia.Rhizomorph.Forward
        Test.@test entry.quality_scores == quality_scores
        Test.@test typeof(entry.quality_scores) == Vector{UInt8}
        Test.@test all(q -> 0 <= q <= 60, entry.quality_scores)  # Valid Phred range
    end

    Test.@testset "QualityEvidenceEntry - Quality Score Validation" begin
        # Quality scores should be UInt8 in range 0-60 for standard Phred scores
        # (though joint scores can go up to 255)

        valid_scores = UInt8[20, 30, 40]
        entry = Mycelia.Rhizomorph.QualityEvidenceEntry(5, Mycelia.Rhizomorph.Forward, valid_scores)
        Test.@test !isempty(entry.quality_scores)
        Test.@test all(q -> q isa UInt8, entry.quality_scores)
    end

    Test.@testset "EdgeEvidenceEntry - Basic Construction" begin
        # EdgeEvidenceEntry tracks transitions between positions

        entry = Mycelia.Rhizomorph.EdgeEvidenceEntry(10, 11, Mycelia.Rhizomorph.Forward)

        Test.@test entry.from_position == 10
        Test.@test entry.to_position == 11
        Test.@test entry.strand == Mycelia.Rhizomorph.Forward
        Test.@test typeof(entry.from_position) == Int
        Test.@test typeof(entry.to_position) == Int
    end

    Test.@testset "EdgeEvidenceEntry - Adjacent Positions" begin
        # For k-mer graphs, edges typically connect adjacent positions
        entry = Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward)

        Test.@test entry.to_position == entry.from_position + 1
    end

    Test.@testset "EdgeQualityEvidenceEntry - With Qualities" begin
        # EdgeQualityEvidenceEntry includes quality at both ends

        from_quality = UInt8[30, 35, 32]
        to_quality = UInt8[28, 40, 38]

        entry = Mycelia.Rhizomorph.EdgeQualityEvidenceEntry(
            10, 11, Mycelia.Rhizomorph.Forward, from_quality, to_quality
        )

        Test.@test entry.from_position == 10
        Test.@test entry.to_position == 11
        Test.@test entry.strand == Mycelia.Rhizomorph.Forward
        Test.@test entry.from_quality == from_quality
        Test.@test entry.to_quality == to_quality
        Test.@test typeof(entry.from_quality) == Vector{UInt8}
        Test.@test typeof(entry.to_quality) == Vector{UInt8}
    end
end

Test.@testset "Evidence Storage Structure - Nested Dictionaries" begin
    Test.@testset "Evidence Dictionary Structure" begin
        # Evidence is stored as: dataset_id -> observation_id -> Set{EvidenceEntry}
        # This provides O(1) lookup by dataset or observation

        evidence = Dict{String, Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}}()

        # Add evidence from dataset "dataset_01", observation "read_001"
        dataset_id = "dataset_01"
        observation_id = "read_001"

        if !haskey(evidence, dataset_id)
            evidence[dataset_id] = Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}()
        end

        if !haskey(evidence[dataset_id], observation_id)
            evidence[dataset_id][observation_id] = Set{Mycelia.Rhizomorph.EvidenceEntry}()
        end

        push!(evidence[dataset_id][observation_id],
            Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        push!(evidence[dataset_id][observation_id],
            Mycelia.Rhizomorph.EvidenceEntry(10, Mycelia.Rhizomorph.Forward))

        # Test structure
        Test.@test haskey(evidence, dataset_id)
        Test.@test haskey(evidence[dataset_id], observation_id)
        Test.@test length(evidence[dataset_id][observation_id]) == 2
    end

    Test.@testset "Multiple Datasets and Observations" begin
        evidence = Dict{String, Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}}()

        # Dataset 1, Observation 1
        evidence["dataset_01"] = Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}()
        evidence["dataset_01"]["read_001"] = Set{Mycelia.Rhizomorph.EvidenceEntry}()
        push!(evidence["dataset_01"]["read_001"],
            Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))

        # Dataset 1, Observation 2
        evidence["dataset_01"]["read_002"] = Set{Mycelia.Rhizomorph.EvidenceEntry}()
        push!(evidence["dataset_01"]["read_002"],
            Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))

        # Dataset 2, Observation 1
        evidence["dataset_02"] = Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}()
        evidence["dataset_02"]["read_001"] = Set{Mycelia.Rhizomorph.EvidenceEntry}()
        push!(evidence["dataset_02"]["read_001"],
            Mycelia.Rhizomorph.EvidenceEntry(3, Mycelia.Rhizomorph.Reverse))

        # Verify structure
        Test.@test length(keys(evidence)) == 2  # Two datasets
        Test.@test length(keys(evidence["dataset_01"])) == 2  # Two observations in dataset 1
        Test.@test length(keys(evidence["dataset_02"])) == 1  # One observation in dataset 2
    end

    Test.@testset "Evidence Efficiency - O(1) Dataset Query" begin
        evidence = Dict{String, Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}}()

        # Add evidence from multiple datasets
        for dataset_num in 1:100
            dataset_id = "dataset_$(lpad(dataset_num, 3, '0'))"
            evidence[dataset_id] = Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}()
            for obs_num in 1:10
                obs_id = "read_$(lpad(obs_num, 3, '0'))"
                evidence[dataset_id][obs_id] = Set{Mycelia.Rhizomorph.EvidenceEntry}()
                push!(evidence[dataset_id][obs_id],
                    Mycelia.Rhizomorph.EvidenceEntry(obs_num, Mycelia.Rhizomorph.Forward))
            end
        end

        # O(1) access to specific dataset
        target_dataset = "dataset_050"
        Test.@test haskey(evidence, target_dataset)
        Test.@test length(evidence[target_dataset]) == 10

        # O(1) access to specific observation within dataset
        target_obs = "read_005"
        Test.@test haskey(evidence[target_dataset], target_obs)
    end

    Test.@testset "Memory Efficiency - Observation ID as Key" begin
        # Key insight: observation_id is stored ONCE as dict key,
        # not repeated in every EvidenceEntry

        evidence = Dict{String, Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}}()
        observation_id = "very_long_observation_identifier_12345"
        dataset_id = "dataset_01"

        # Add multiple evidence entries for same observation
        evidence[dataset_id] = Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}()
        evidence[dataset_id][observation_id] = Set{Mycelia.Rhizomorph.EvidenceEntry}()
        for pos in 1:100
            push!(evidence[dataset_id][observation_id],
                Mycelia.Rhizomorph.EvidenceEntry(pos, Mycelia.Rhizomorph.Forward))
        end

        # observation_id is stored once as key, not 100 times in entries
        Test.@test haskey(evidence[dataset_id], observation_id)
        Test.@test length(evidence[dataset_id][observation_id]) == 100

        # Each entry only stores position and strand (much smaller)
        for entry in evidence[dataset_id][observation_id]
            Test.@test !hasfield(typeof(entry), :observation_id)
        end
    end
end

Test.@testset "Evidence Structure Types" begin
    Test.@testset "Type Stability" begin
        # All evidence structures should be concretely typed

        entry = Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward)
        Test.@test isconcretetype(typeof(entry))

        qual_entry = Mycelia.Rhizomorph.QualityEvidenceEntry(5, Mycelia.Rhizomorph.Forward, UInt8[30, 35])
        Test.@test isconcretetype(typeof(qual_entry))

        edge_entry = Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward)
        Test.@test isconcretetype(typeof(edge_entry))

        edge_qual_entry = Mycelia.Rhizomorph.EdgeQualityEvidenceEntry(
            5, 6, Mycelia.Rhizomorph.Forward, UInt8[30], UInt8[35]
        )
        Test.@test isconcretetype(typeof(edge_qual_entry))
    end
end

println("✓ Evidence structures tests executed")
