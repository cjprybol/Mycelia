# KmerVertexData Test
#
# Tests for the correct KmerVertexData structure as specified in the
# rhizomorph graph ecosystem planning document.
#
# KEY REQUIREMENTS:
# 1. Store observed k-mer (NOT canonical)
# 2. Use nested Dict evidence structure
# 3. Support efficient dataset/observation queries
#
# Run with: julia --project=. test/4_assembly/kmer_vertex_data_test.jl

import Test
import Mycelia
import Kmers
import BioSequences

Test.@testset "KmerVertexData - Correct Structure" begin

    Test.@testset "Basic Construction - Stores OBSERVED K-mer" begin
        # CRITICAL: Vertex should store the k-mer AS OBSERVED, not canonical

        # Create a DNA k-mer
        kmer = Kmers.DNAKmer{3}("ATG")

        # Create vertex data with empty evidence
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Test that it stores the EXACT k-mer we gave it
        Test.@test vertex_data.Kmer == kmer
        Test.@test string(vertex_data.Kmer) == "ATG"

        # Should NOT have a field called "canonical_kmer"
        Test.@test !hasfield(typeof(vertex_data), :canonical_kmer)

        # Should have field called "Kmer" (capital K)
        Test.@test hasfield(typeof(vertex_data), :Kmer)
    end

    Test.@testset "Stores Forward vs Reverse K-mers Separately" begin
        # Two observations of reverse complement k-mers should create
        # TWO DIFFERENT vertices (strand-specific construction)

        forward_kmer = Kmers.DNAKmer{3}("ATG")
        reverse_kmer = BioSequences.reverse_complement(forward_kmer)  # CAT

        vertex_forward = Mycelia.Rhizomorph.KmerVertexData(forward_kmer)
        vertex_reverse = Mycelia.Rhizomorph.KmerVertexData(reverse_kmer)

        # These should store DIFFERENT k-mers
        Test.@test vertex_forward.Kmer == forward_kmer
        Test.@test vertex_reverse.Kmer == reverse_kmer
        Test.@test vertex_forward.Kmer != vertex_reverse.Kmer

        # Verify the actual sequences
        Test.@test string(vertex_forward.Kmer) == "ATG"
        Test.@test string(vertex_reverse.Kmer) == "CAT"
    end

    Test.@testset "Evidence Structure - Nested Dictionaries" begin
        # Evidence should be: Dict{String, Dict{String, Set{EvidenceEntry}}}

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Should have evidence field
        Test.@test hasfield(typeof(vertex_data), :evidence)

        # Evidence should be nested dict structure
        Test.@test typeof(vertex_data.evidence) ==
            Dict{String, Dict{String, Set{Mycelia.Rhizomorph.EvidenceEntry}}}

        # Should start empty
        Test.@test isempty(vertex_data.evidence)
    end

    Test.@testset "Adding Evidence from Single Observation" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Add evidence from one observation
        dataset_id = "dataset_01"
        observation_id = "read_001"

        Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, observation_id,
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))

        # Verify structure
        Test.@test haskey(vertex_data.evidence, dataset_id)
        Test.@test haskey(vertex_data.evidence[dataset_id], observation_id)
        Test.@test length(vertex_data.evidence[dataset_id][observation_id]) == 1
    end

    Test.@testset "Adding Multiple Evidence from Same Observation" begin
        # A long read might visit the same k-mer multiple times
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        dataset_id = "dataset_01"
        observation_id = "read_001"

        # Add evidence at multiple positions
        Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, observation_id,
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, observation_id,
                             Mycelia.Rhizomorph.EvidenceEntry(120, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, observation_id,
                             Mycelia.Rhizomorph.EvidenceEntry(245, Mycelia.Rhizomorph.Forward))

        # Should have 3 evidence entries from same observation
        Test.@test length(vertex_data.evidence[dataset_id][observation_id]) == 3
    end

    Test.@testset "Evidence from Multiple Observations" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        dataset_id = "dataset_01"

        # Evidence from multiple reads
        Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, "read_002",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, "read_003",
                             Mycelia.Rhizomorph.EvidenceEntry(12, Mycelia.Rhizomorph.Reverse))

        # Should have 3 observations
        Test.@test length(vertex_data.evidence[dataset_id]) == 3

        # Each observation has its own evidence
        Test.@test haskey(vertex_data.evidence[dataset_id], "read_001")
        Test.@test haskey(vertex_data.evidence[dataset_id], "read_002")
        Test.@test haskey(vertex_data.evidence[dataset_id], "read_003")
    end

    Test.@testset "Evidence from Multiple Datasets" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Evidence from different datasets
        Mycelia.Rhizomorph.add_evidence!(vertex_data, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex_data, "dataset_02", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(8, Mycelia.Rhizomorph.Forward))
        Mycelia.Rhizomorph.add_evidence!(vertex_data, "dataset_03", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(3, Mycelia.Rhizomorph.Forward))

        # Should have 3 datasets
        Test.@test length(keys(vertex_data.evidence)) == 3
        Test.@test haskey(vertex_data.evidence, "dataset_01")
        Test.@test haskey(vertex_data.evidence, "dataset_02")
        Test.@test haskey(vertex_data.evidence, "dataset_03")
    end

    Test.@testset "Efficient Dataset Query - O(1)" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Add evidence from many datasets
        for i in 1:100
            dataset_id = "dataset_$(lpad(i, 3, '0'))"
            Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, "read_001",
                                 Mycelia.Rhizomorph.EvidenceEntry(i, Mycelia.Rhizomorph.Forward))
        end

        # O(1) access to specific dataset
        target_dataset = "dataset_050"
        dataset_evidence = Mycelia.Rhizomorph.get_dataset_evidence(vertex_data, target_dataset)

        Test.@test !isnothing(dataset_evidence)
        Test.@test haskey(dataset_evidence, "read_001")
    end

    Test.@testset "Efficient Observation Query - O(1)" begin
        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        dataset_id = "dataset_01"

        # Add evidence from many observations
        for i in 1:100
            obs_id = "read_$(lpad(i, 3, '0'))"
            Mycelia.Rhizomorph.add_evidence!(vertex_data, dataset_id, obs_id,
                                 Mycelia.Rhizomorph.EvidenceEntry(i, Mycelia.Rhizomorph.Forward))
        end

        # O(1) access to specific observation
        target_obs = "read_050"
        obs_evidence = Mycelia.Rhizomorph.get_observation_evidence(vertex_data, dataset_id, target_obs)

        Test.@test !isnothing(obs_evidence)
        Test.@test length(obs_evidence) >= 1
    end

    Test.@testset "Type Parameterization" begin
        # KmerVertexData should be parameterized by k-mer type

        dna_kmer = Kmers.DNAKmer{3}("ATG")
        dna_vertex = Mycelia.Rhizomorph.KmerVertexData(dna_kmer)

        rna_kmer = Kmers.RNAKmer{3}("AUG")
        rna_vertex = Mycelia.Rhizomorph.KmerVertexData(rna_kmer)

        aa_kmer = Kmers.AAKmer{3}("MET")
        aa_vertex = Mycelia.Rhizomorph.KmerVertexData(aa_kmer)

        # Should have different parameterized types
        Test.@test typeof(dna_vertex) != typeof(rna_vertex)
        Test.@test typeof(dna_vertex) != typeof(aa_vertex)
        Test.@test typeof(rna_vertex) != typeof(aa_vertex)

        # But all should be KmerVertexData
        Test.@test typeof(dna_vertex) <: Mycelia.Rhizomorph.KmerVertexData
        Test.@test typeof(rna_vertex) <: Mycelia.Rhizomorph.KmerVertexData
        Test.@test typeof(aa_vertex) <: Mycelia.Rhizomorph.KmerVertexData
    end

    Test.@testset "Different K Sizes" begin
        # Should work with different k-mer sizes

        k3_kmer = Kmers.DNAKmer{3}("ATG")
        k3_vertex = Mycelia.Rhizomorph.KmerVertexData(k3_kmer)
        Test.@test length(k3_vertex.Kmer) == 3

        k31_kmer = Kmers.DNAKmer{31}("ATGATGATGATGATGATGATGATGATGATGA")
        k31_vertex = Mycelia.Rhizomorph.KmerVertexData(k31_kmer)
        Test.@test length(k31_vertex.Kmer) == 31

        k101_kmer = Kmers.DNAKmer{101}("A" ^ 101)
        k101_vertex = Mycelia.Rhizomorph.KmerVertexData(k101_kmer)
        Test.@test length(k101_vertex.Kmer) == 101
    end
end

Test.@testset "KmerVertexData - Critical Anti-Patterns" begin

    Test.@testset "NEVER Store Canonical K-mer" begin
        # The OLD WRONG implementation stored canonical_kmer
        # The NEW CORRECT implementation stores observed Kmer

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Should NOT have canonical_kmer field
        Test.@test !hasfield(typeof(vertex_data), :canonical_kmer)

        # If someone tries to access canonical representation, they should
        # compute it from the stored k-mer, not store it
        Test.@test hasmethod(BioSequences.canonical, (typeof(kmer),))
    end

    Test.@testset "NEVER Use Flat Coverage Vector" begin
        # The OLD WRONG implementation used: Vector{Tuple{Int, Int, StrandOrientation}}
        # The NEW CORRECT implementation uses: Dict{String, Dict{String, Set{EvidenceEntry}}}

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        # Should NOT have a flat coverage vector
        Test.@test !hasfield(typeof(vertex_data), :coverage) ||
                   (typeof(vertex_data.coverage) != Vector{Tuple{Int, Int, Mycelia.Rhizomorph.StrandOrientation}})

        # Should have nested evidence dict
        Test.@test hasfield(typeof(vertex_data), :evidence)
        Test.@test typeof(vertex_data.evidence) <: Dict
    end

    Test.@testset "Evidence NOT Tuples - Use Structs" begin
        # OLD WRONG: Tuple{Int, Int, StrandOrientation}
        # NEW CORRECT: EvidenceEntry struct

        kmer = Kmers.DNAKmer{3}("ATG")
        vertex_data = Mycelia.Rhizomorph.KmerVertexData(kmer)

        Mycelia.Rhizomorph.add_evidence!(vertex_data, "dataset_01", "read_001",
                             Mycelia.Rhizomorph.EvidenceEntry(5, Mycelia.Rhizomorph.Forward))

        # Get an evidence entry
        evidence_set = vertex_data.evidence["dataset_01"]["read_001"]
        entry = first(evidence_set)

        # Should be a struct, not a tuple
        Test.@test !isa(entry, Tuple)
        Test.@test isa(entry, Mycelia.Rhizomorph.EvidenceEntry)
        Test.@test hasfield(typeof(entry), :position)
        Test.@test hasfield(typeof(entry), :strand)
    end
end

println("✓ KmerVertexData tests defined (will fail until implementation)")
