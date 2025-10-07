# KmerEdgeData Test
#
# Tests for the correct KmerEdgeData structure as specified in the
# rhizomorph graph ecosystem planning document.
#
# KEY REQUIREMENTS:
# 1. Use nested Dict evidence structure (not flat vectors)
# 2. Use EdgeEvidenceEntry structs (not tuples)
# 3. NO premature weight calculation
# 4. Support efficient dataset/observation queries
#
# Run with: julia --project=. test/4_assembly/kmer_edge_data_test.jl

import Test
import Mycelia

Test.@testset "KmerEdgeData - Correct Structure" begin

    Test.@testset "Basic Construction" begin
        # Edge data should use nested dict evidence structure

        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        # Should have evidence field
        Test.@test hasfield(typeof(edge_data), :evidence)

        # Evidence should be nested dict structure
        Test.@test typeof(edge_data.evidence) ==
            Dict{String, Dict{String, Set{Mycelia.Rhizomorph.EdgeEvidenceEntry}}}

        # Should start empty
        Test.@test isempty(edge_data.evidence)
    end

    Test.@testset "Adding Edge Evidence" begin
        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        dataset_id = "dataset_01"
        observation_id = "read_001"

        # Add evidence for edge transition
        Mycelia.Rhizomorph.add_evidence!(edge_data, dataset_id, observation_id,
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))

        # Verify structure
        Test.@test haskey(edge_data.evidence, dataset_id)
        Test.@test haskey(edge_data.evidence[dataset_id], observation_id)
        Test.@test length(edge_data.evidence[dataset_id][observation_id]) == 1

        # Get the entry and verify it's a struct
        entry = first(edge_data.evidence[dataset_id][observation_id])
        Test.@test isa(entry, Mycelia.Rhizomorph.EdgeEvidenceEntry)
        Test.@test entry.from_position == 5
        Test.@test entry.to_position == 6
        Test.@test entry.strand == Mycelia.Rhizomorph.Forward
    end

    Test.@testset "Multiple Edge Observations" begin
        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        dataset_id = "dataset_01"

        # Same edge seen multiple times in different reads
        Mycelia.Rhizomorph.add_evidence!(edge_data, dataset_id, "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge_data, dataset_id, "read_002",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(10, 11, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge_data, dataset_id, "read_003",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(15, 16, Mycelia.Rhizomorph.Reverse))

        # Should have 3 observations
        Test.@test length(edge_data.evidence[dataset_id]) == 3
    end

    Test.@testset "Same Read Multiple Edge Observations" begin
        # A long read might traverse the same edge multiple times
        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        dataset_id = "dataset_01"
        observation_id = "read_001"

        Mycelia.Rhizomorph.add_evidence!(edge_data, dataset_id, observation_id,
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge_data, dataset_id, observation_id,
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(120, 121, Mycelia.Rhizomorph.Forward))

        # Should have 2 evidence entries from same read
        Test.@test length(edge_data.evidence[dataset_id][observation_id]) == 2
    end

    Test.@testset "Evidence from Multiple Datasets" begin
        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        # Edge seen in multiple datasets
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_02", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(8, 9, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_03", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(3, 4, Mycelia.Rhizomorph.Forward))

        # Should have 3 datasets
        Test.@test length(keys(edge_data.evidence)) == 3
    end

    Test.@testset "Efficient Dataset Query" begin
        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        # Add evidence from many datasets
        for i in 1:100
            dataset_id = "dataset_$(lpad(i, 3, '0'))"
            Mycelia.Rhizomorph.add_evidence!(edge_data, dataset_id, "read_001",
                                 Mycelia.Rhizomorph.EdgeEvidenceEntry(i, i+1, Mycelia.Rhizomorph.Forward))
        end

        # O(1) access to specific dataset
        target_dataset = "dataset_050"
        dataset_evidence = Mycelia.Rhizomorph.get_dataset_evidence(edge_data, target_dataset)

        Test.@test !isnothing(dataset_evidence)
        Test.@test haskey(dataset_evidence, "read_001")
    end

    Test.@testset "Efficient Observation Query" begin
        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        dataset_id = "dataset_01"

        # Add evidence from many observations
        for i in 1:100
            obs_id = "read_$(lpad(i, 3, '0'))"
            Mycelia.Rhizomorph.add_evidence!(edge_data, dataset_id, obs_id,
                                 Mycelia.Rhizomorph.EdgeEvidenceEntry(i, i+1, Mycelia.Rhizomorph.Forward))
        end

        # O(1) access to specific observation
        target_obs = "read_050"
        obs_evidence = Mycelia.Rhizomorph.get_observation_evidence(edge_data, dataset_id, target_obs)

        Test.@test !isnothing(obs_evidence)
        Test.@test length(obs_evidence) >= 1
    end
end

Test.@testset "KmerEdgeData - NO Premature Computation" begin

    Test.@testset "No Stored Weight" begin
        # OLD WRONG: Edge stored computed weight
        # NEW CORRECT: Compute weight on-demand from evidence

        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        # Add evidence
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))

        # Should NOT have pre-computed weight field
        # (Or if it exists, should be computed lazily)
        Test.@test !hasfield(typeof(edge_data), :weight) ||
                   hasmethod(Mycelia.Rhizomorph.compute_edge_weight, (typeof(edge_data),))
    end

    Test.@testset "No Stored Strand Orientations" begin
        # OLD WRONG: Stored src_strand, dst_strand in edge data
        # NEW CORRECT: Strand info is in EdgeEvidenceEntry, varies per observation

        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        # Same edge can be traversed in different strand orientations
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(10, 11, Mycelia.Rhizomorph.Reverse))

        # Edge data should NOT have single src_strand/dst_strand fields
        # (Different observations may have different strand orientations)
        if hasfield(typeof(edge_data), :src_strand)
            # If these fields exist, they should be computed/dominant, not stored
            @warn "Edge has strand fields - ensure they're computed, not stored"
        end

        # Verify both strand orientations are in evidence
        all_entries = vcat([collect(obs_evidence)
                           for obs_evidence in values(edge_data.evidence["dataset_01"])]...)
        strands = [entry.strand for entry in all_entries]
        Test.@test Mycelia.Rhizomorph.Forward in strands
        Test.@test Mycelia.Rhizomorph.Reverse in strands
    end
end

Test.@testset "KmerEdgeData - Critical Anti-Patterns" begin

    Test.@testset "NEVER Use Nested Tuples" begin
        # OLD WRONG: Vector{Tuple{Tuple{Int,Int,Strand}, Tuple{Int,Int,Strand}}}
        # NEW CORRECT: Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}

        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))

        # Get an evidence entry
        evidence_set = edge_data.evidence["dataset_01"]["read_001"]
        entry = first(evidence_set)

        # Should be EdgeEvidenceEntry struct, not nested tuples
        Test.@test !isa(entry, Tuple)
        Test.@test isa(entry, Mycelia.Rhizomorph.EdgeEvidenceEntry)
        Test.@test hasfield(typeof(entry), :from_position)
        Test.@test hasfield(typeof(entry), :to_position)
        Test.@test hasfield(typeof(entry), :strand)
    end

    Test.@testset "NEVER Use Flat Coverage Vector" begin
        # OLD WRONG: Vector{Tuple{...}}
        # NEW CORRECT: Dict{String, Dict{String, Set{EdgeEvidenceEntry}}}

        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        # Should NOT have flat coverage field
        Test.@test !hasfield(typeof(edge_data), :coverage) ||
                   (typeof(edge_data.coverage) != Vector)

        # Should have nested evidence dict
        Test.@test hasfield(typeof(edge_data), :evidence)
        Test.@test typeof(edge_data.evidence) <: Dict
    end

    Test.@testset "Evidence Struct Has Correct Fields" begin
        # EdgeEvidenceEntry should have exactly these fields:
        # - from_position::Int
        # - to_position::Int
        # - strand::StrandOrientation

        entry = Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward)

        Test.@test hasfield(typeof(entry), :from_position)
        Test.@test hasfield(typeof(entry), :to_position)
        Test.@test hasfield(typeof(entry), :strand)

        # Should NOT have observation_id (that's the dict key)
        Test.@test !hasfield(typeof(entry), :observation_id)
    end
end

Test.@testset "KmerEdgeData - Weight Calculation" begin

    Test.@testset "Compute Weight from Evidence On-Demand" begin
        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        # Add evidence
        for i in 1:10
            Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_$(lpad(i, 3, '0'))",
                                 Mycelia.Rhizomorph.EdgeEvidenceEntry(i, i+1, Mycelia.Rhizomorph.Forward))
        end

        # Compute weight from evidence count
        weight = Mycelia.Rhizomorph.compute_edge_weight(edge_data)

        # Weight should be based on number of observations
        # (e.g., count of observations, or sum of evidence entries)
        total_observations = sum(length(obs_evidence)
                               for obs_evidence in values(edge_data.evidence["dataset_01"]))
        Test.@test weight > 0
        Test.@test weight isa Float64
    end

    Test.@testset "Coverage from Evidence" begin
        edge_data = Mycelia.Rhizomorph.KmerEdgeData()

        # Add varying amounts of evidence per observation
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(5, 6, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(50, 51, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(edge_data, "dataset_01", "read_002",
                             Mycelia.Rhizomorph.EdgeEvidenceEntry(10, 11, Mycelia.Rhizomorph.Forward))

        # Coverage should be computable from evidence
        coverage = Mycelia.Rhizomorph.compute_edge_coverage(edge_data)

        # Different ways to define coverage:
        # - Number of unique observations: 2
        # - Total evidence entries: 3
        Test.@test coverage >= 2  # At least 2 observations
    end
end

println("✓ KmerEdgeData tests defined (will fail until implementation)")
