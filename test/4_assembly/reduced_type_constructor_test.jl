# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/reduced_type_constructor_test.jl")'
# ```
#
# And to turn this file into a jupyter notebook, run:
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("test/4_assembly/reduced_type_constructor_test.jl", "test/4_assembly", execute=false)'
# ```

## If running Literate notebook, ensure the package is activated:
## import Pkg
## if isinteractive()
##     Pkg.activate(joinpath(@__DIR__, "..", ".."))
## end
## using Revise

# Reduced Type Constructor and Field Tests
#
# Tests for all 10 new reduced vertex/edge data types in the Rhizomorph module:
#
# Ultralight (4): UltralightKmerVertexData, UltralightBioSequenceVertexData,
#                 UltralightStringVertexData, UltralightEdgeData
#
# Ultralight Quality (3): UltralightQualityKmerVertexData,
#                          UltralightQualityBioSequenceVertexData,
#                          UltralightQualityEdgeData
#
# Lightweight Quality (3): LightweightQualityKmerVertexData,
#                           LightweightQualityBioSequenceVertexData,
#                           LightweightQualityEdgeData
#
# KEY REQUIREMENTS VERIFIED:
# 1. Default constructors work and produce valid instances
# 2. All expected fields exist with correct types
# 3. Mutable fields (total_count) can be modified
# 4. Initial values are correct (total_count=0, empty dicts, empty vectors)
# 5. Fields absent in lighter variants are confirmed absent:
#    - ultralight types do NOT have :dataset_observations
#    - non-quality types do NOT have :joint_quality
#    - edge data does NOT have :evidence (reduced types)

import Test
import Mycelia

Test.@testset "Reduced Type Constructor and Field Tests" begin

    # ------------------------------------------------------------------
    # Shared test inputs
    # ------------------------------------------------------------------
    kmer = Mycelia.Kmers.DNAKmer{3}("ATG")
    biosequence = Mycelia.BioSequences.LongDNA{4}("ATGATG")

    # ==================================================================
    # ULTRALIGHT VERTEX TYPES
    # ==================================================================

    Test.@testset "UltralightKmerVertexData{T}" begin
        Test.@testset "Construction" begin
            v = Mycelia.Rhizomorph.UltralightKmerVertexData(kmer)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightKmerVertexData)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightKmerVertexData{typeof(kmer)})
        end

        Test.@testset "Type parameterization" begin
            v = Mycelia.Rhizomorph.UltralightKmerVertexData(kmer)
            Test.@test typeof(v) <: Mycelia.Rhizomorph.UltralightKmerVertexData
        end

        Test.@testset "Expected fields present" begin
            v = Mycelia.Rhizomorph.UltralightKmerVertexData(kmer)
            Test.@test hasfield(typeof(v), :Kmer)
            Test.@test hasfield(typeof(v), :total_count)
            Test.@test hasfield(typeof(v), :dataset_counts)
        end

        Test.@testset "Initial field values" begin
            v = Mycelia.Rhizomorph.UltralightKmerVertexData(kmer)
            Test.@test v.Kmer == kmer
            Test.@test string(v.Kmer) == "ATG"
            Test.@test v.total_count == 0
            Test.@test isempty(v.dataset_counts)
            Test.@test v.dataset_counts isa Dict{String, Int}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            v = Mycelia.Rhizomorph.UltralightKmerVertexData(kmer)
            v.total_count = 42
            Test.@test v.total_count == 42
        end

        Test.@testset "Absent fields (ultralight has no observation tracking)" begin
            v = Mycelia.Rhizomorph.UltralightKmerVertexData(kmer)
            # No observation ID tracking (that lives only in Lightweight variants)
            Test.@test !hasfield(typeof(v), :dataset_observations)
            # No quality (that lives only in quality variants)
            Test.@test !hasfield(typeof(v), :joint_quality)
            Test.@test !hasfield(typeof(v), :dataset_joint_quality)
            # No full evidence dict (that lives only in full, non-reduced types)
            Test.@test !hasfield(typeof(v), :evidence)
        end

        Test.@testset "Explicit type parameter constructor" begin
            v = Mycelia.Rhizomorph.UltralightKmerVertexData{typeof(kmer)}(kmer)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightKmerVertexData{typeof(kmer)})
            Test.@test v.Kmer == kmer
        end
    end

    Test.@testset "UltralightBioSequenceVertexData{T}" begin
        Test.@testset "Construction" begin
            v = Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightBioSequenceVertexData)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightBioSequenceVertexData{typeof(biosequence)})
        end

        Test.@testset "Expected fields present" begin
            v = Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence)
            Test.@test hasfield(typeof(v), :sequence)
            Test.@test hasfield(typeof(v), :total_count)
            Test.@test hasfield(typeof(v), :dataset_counts)
        end

        Test.@testset "Initial field values" begin
            v = Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence)
            Test.@test v.sequence == biosequence
            Test.@test v.total_count == 0
            Test.@test isempty(v.dataset_counts)
            Test.@test v.dataset_counts isa Dict{String, Int}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            v = Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence)
            v.total_count = 7
            Test.@test v.total_count == 7
        end

        Test.@testset "Absent fields" begin
            v = Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence)
            Test.@test !hasfield(typeof(v), :dataset_observations)
            Test.@test !hasfield(typeof(v), :joint_quality)
            Test.@test !hasfield(typeof(v), :dataset_joint_quality)
            Test.@test !hasfield(typeof(v), :evidence)
            # BioSequence variant uses :sequence, not :Kmer
            Test.@test !hasfield(typeof(v), :Kmer)
        end

        Test.@testset "Explicit type parameter constructor" begin
            v = Mycelia.Rhizomorph.UltralightBioSequenceVertexData{typeof(biosequence)}(biosequence)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightBioSequenceVertexData{typeof(biosequence)})
            Test.@test v.sequence == biosequence
        end
    end

    Test.@testset "UltralightStringVertexData" begin
        Test.@testset "Construction" begin
            v = Mycelia.Rhizomorph.UltralightStringVertexData("ATG")
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightStringVertexData)
        end

        Test.@testset "Expected fields present" begin
            v = Mycelia.Rhizomorph.UltralightStringVertexData("ATG")
            Test.@test hasfield(typeof(v), :string_value)
            Test.@test hasfield(typeof(v), :total_count)
            Test.@test hasfield(typeof(v), :dataset_counts)
        end

        Test.@testset "Initial field values" begin
            v = Mycelia.Rhizomorph.UltralightStringVertexData("ATG")
            Test.@test v.string_value == "ATG"
            Test.@test v.total_count == 0
            Test.@test isempty(v.dataset_counts)
            Test.@test v.dataset_counts isa Dict{String, Int}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            v = Mycelia.Rhizomorph.UltralightStringVertexData("ATG")
            v.total_count = 100
            Test.@test v.total_count == 100
        end

        Test.@testset "Absent fields" begin
            v = Mycelia.Rhizomorph.UltralightStringVertexData("ATG")
            Test.@test !hasfield(typeof(v), :dataset_observations)
            Test.@test !hasfield(typeof(v), :joint_quality)
            Test.@test !hasfield(typeof(v), :dataset_joint_quality)
            Test.@test !hasfield(typeof(v), :evidence)
            # String variant uses :string_value, not :Kmer or :sequence
            Test.@test !hasfield(typeof(v), :Kmer)
            Test.@test !hasfield(typeof(v), :sequence)
        end
    end

    # ==================================================================
    # ULTRALIGHT QUALITY VERTEX TYPES
    # ==================================================================

    Test.@testset "UltralightQualityKmerVertexData{T}" begin
        Test.@testset "Construction" begin
            v = Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightQualityKmerVertexData)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightQualityKmerVertexData{typeof(kmer)})
        end

        Test.@testset "Expected fields present" begin
            v = Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer)
            Test.@test hasfield(typeof(v), :Kmer)
            Test.@test hasfield(typeof(v), :total_count)
            Test.@test hasfield(typeof(v), :dataset_counts)
            Test.@test hasfield(typeof(v), :joint_quality)
            Test.@test hasfield(typeof(v), :dataset_joint_quality)
        end

        Test.@testset "Initial field values" begin
            v = Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer)
            Test.@test v.Kmer == kmer
            Test.@test v.total_count == 0
            Test.@test isempty(v.dataset_counts)
            Test.@test v.dataset_counts isa Dict{String, Int}
            Test.@test isempty(v.joint_quality)
            Test.@test v.joint_quality isa Vector{UInt8}
            Test.@test isempty(v.dataset_joint_quality)
            Test.@test v.dataset_joint_quality isa Dict{String, Vector{UInt8}}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            v = Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer)
            v.total_count = 5
            Test.@test v.total_count == 5
        end

        Test.@testset "Absent fields (ultralight-quality has no obs tracking)" begin
            v = Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer)
            # No observation ID tracking — that is added in lightweight_quality
            Test.@test !hasfield(typeof(v), :dataset_observations)
            # No full evidence dict
            Test.@test !hasfield(typeof(v), :evidence)
        end

        Test.@testset "Explicit type parameter constructor" begin
            v = Mycelia.Rhizomorph.UltralightQualityKmerVertexData{typeof(kmer)}(kmer)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightQualityKmerVertexData{typeof(kmer)})
        end
    end

    Test.@testset "UltralightQualityBioSequenceVertexData{T}" begin
        Test.@testset "Construction" begin
            v = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData{typeof(biosequence)})
        end

        Test.@testset "Expected fields present" begin
            v = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence)
            Test.@test hasfield(typeof(v), :sequence)
            Test.@test hasfield(typeof(v), :total_count)
            Test.@test hasfield(typeof(v), :dataset_counts)
            Test.@test hasfield(typeof(v), :joint_quality)
            Test.@test hasfield(typeof(v), :dataset_joint_quality)
        end

        Test.@testset "Initial field values" begin
            v = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence)
            Test.@test v.sequence == biosequence
            Test.@test v.total_count == 0
            Test.@test isempty(v.dataset_counts)
            Test.@test v.dataset_counts isa Dict{String, Int}
            Test.@test isempty(v.joint_quality)
            Test.@test v.joint_quality isa Vector{UInt8}
            Test.@test isempty(v.dataset_joint_quality)
            Test.@test v.dataset_joint_quality isa Dict{String, Vector{UInt8}}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            v = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence)
            v.total_count = 3
            Test.@test v.total_count == 3
        end

        Test.@testset "Absent fields" begin
            v = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence)
            Test.@test !hasfield(typeof(v), :dataset_observations)
            Test.@test !hasfield(typeof(v), :evidence)
            Test.@test !hasfield(typeof(v), :Kmer)
        end

        Test.@testset "Explicit type parameter constructor" begin
            v = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData{typeof(biosequence)}(biosequence)
            Test.@test isa(v, Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData{typeof(biosequence)})
            Test.@test v.sequence == biosequence
        end
    end

    # ==================================================================
    # LIGHTWEIGHT QUALITY VERTEX TYPES
    # ==================================================================

    Test.@testset "LightweightQualityKmerVertexData{T}" begin
        Test.@testset "Construction" begin
            v = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)
            Test.@test isa(v, Mycelia.Rhizomorph.LightweightQualityKmerVertexData)
            Test.@test isa(v, Mycelia.Rhizomorph.LightweightQualityKmerVertexData{typeof(kmer)})
        end

        Test.@testset "Expected fields present" begin
            v = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)
            Test.@test hasfield(typeof(v), :Kmer)
            Test.@test hasfield(typeof(v), :total_count)
            Test.@test hasfield(typeof(v), :dataset_counts)
            Test.@test hasfield(typeof(v), :dataset_observations)
            Test.@test hasfield(typeof(v), :joint_quality)
            Test.@test hasfield(typeof(v), :dataset_joint_quality)
        end

        Test.@testset "Initial field values" begin
            v = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)
            Test.@test v.Kmer == kmer
            Test.@test v.total_count == 0
            Test.@test isempty(v.dataset_counts)
            Test.@test v.dataset_counts isa Dict{String, Int}
            Test.@test isempty(v.dataset_observations)
            Test.@test v.dataset_observations isa Dict{String, Set{String}}
            Test.@test isempty(v.joint_quality)
            Test.@test v.joint_quality isa Vector{UInt8}
            Test.@test isempty(v.dataset_joint_quality)
            Test.@test v.dataset_joint_quality isa Dict{String, Vector{UInt8}}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            v = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)
            v.total_count = 11
            Test.@test v.total_count == 11
        end

        Test.@testset "Absent fields (no full evidence dict)" begin
            v = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)
            Test.@test !hasfield(typeof(v), :evidence)
        end

        Test.@testset "Has observation tracking (unlike ultralight)" begin
            v = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)
            # This is the key differentiator from ultralight_quality
            Test.@test hasfield(typeof(v), :dataset_observations)
        end

        Test.@testset "Explicit type parameter constructor" begin
            v = Mycelia.Rhizomorph.LightweightQualityKmerVertexData{typeof(kmer)}(kmer)
            Test.@test isa(v, Mycelia.Rhizomorph.LightweightQualityKmerVertexData{typeof(kmer)})
        end
    end

    Test.@testset "LightweightQualityBioSequenceVertexData{T}" begin
        Test.@testset "Construction" begin
            v = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence)
            Test.@test isa(v, Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData)
            Test.@test isa(v, Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData{typeof(biosequence)})
        end

        Test.@testset "Expected fields present" begin
            v = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence)
            Test.@test hasfield(typeof(v), :sequence)
            Test.@test hasfield(typeof(v), :total_count)
            Test.@test hasfield(typeof(v), :dataset_counts)
            Test.@test hasfield(typeof(v), :dataset_observations)
            Test.@test hasfield(typeof(v), :joint_quality)
            Test.@test hasfield(typeof(v), :dataset_joint_quality)
        end

        Test.@testset "Initial field values" begin
            v = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence)
            Test.@test v.sequence == biosequence
            Test.@test v.total_count == 0
            Test.@test isempty(v.dataset_counts)
            Test.@test v.dataset_counts isa Dict{String, Int}
            Test.@test isempty(v.dataset_observations)
            Test.@test v.dataset_observations isa Dict{String, Set{String}}
            Test.@test isempty(v.joint_quality)
            Test.@test v.joint_quality isa Vector{UInt8}
            Test.@test isempty(v.dataset_joint_quality)
            Test.@test v.dataset_joint_quality isa Dict{String, Vector{UInt8}}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            v = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence)
            v.total_count = 99
            Test.@test v.total_count == 99
        end

        Test.@testset "Absent fields" begin
            v = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence)
            Test.@test !hasfield(typeof(v), :evidence)
            Test.@test !hasfield(typeof(v), :Kmer)
        end

        Test.@testset "Has observation tracking (unlike ultralight)" begin
            v = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence)
            Test.@test hasfield(typeof(v), :dataset_observations)
        end

        Test.@testset "Explicit type parameter constructor" begin
            v = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData{typeof(biosequence)}(biosequence)
            Test.@test isa(v, Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData{typeof(biosequence)})
            Test.@test v.sequence == biosequence
        end
    end

    # ==================================================================
    # ULTRALIGHT EDGE TYPE
    # ==================================================================

    Test.@testset "UltralightEdgeData" begin
        Test.@testset "Default construction (overlap_length=0)" begin
            e = Mycelia.Rhizomorph.UltralightEdgeData()
            Test.@test isa(e, Mycelia.Rhizomorph.UltralightEdgeData)
        end

        Test.@testset "Construction with explicit overlap_length" begin
            e = Mycelia.Rhizomorph.UltralightEdgeData(5)
            Test.@test isa(e, Mycelia.Rhizomorph.UltralightEdgeData)
            Test.@test e.overlap_length == 5
        end

        Test.@testset "Expected fields present" begin
            e = Mycelia.Rhizomorph.UltralightEdgeData()
            Test.@test hasfield(typeof(e), :overlap_length)
            Test.@test hasfield(typeof(e), :total_count)
            Test.@test hasfield(typeof(e), :dataset_counts)
        end

        Test.@testset "Initial field values" begin
            e = Mycelia.Rhizomorph.UltralightEdgeData()
            Test.@test e.overlap_length == 0
            Test.@test e.total_count == 0
            Test.@test isempty(e.dataset_counts)
            Test.@test e.dataset_counts isa Dict{String, Int}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            e = Mycelia.Rhizomorph.UltralightEdgeData()
            e.total_count = 15
            Test.@test e.total_count == 15
        end

        Test.@testset "Absent fields (ultralight has no observation tracking)" begin
            e = Mycelia.Rhizomorph.UltralightEdgeData()
            Test.@test !hasfield(typeof(e), :dataset_observations)
            Test.@test !hasfield(typeof(e), :from_joint_quality)
            Test.@test !hasfield(typeof(e), :to_joint_quality)
            Test.@test !hasfield(typeof(e), :evidence)
        end
    end

    # ==================================================================
    # ULTRALIGHT QUALITY EDGE TYPE
    # ==================================================================

    Test.@testset "UltralightQualityEdgeData" begin
        Test.@testset "Default construction (overlap_length=0)" begin
            e = Mycelia.Rhizomorph.UltralightQualityEdgeData()
            Test.@test isa(e, Mycelia.Rhizomorph.UltralightQualityEdgeData)
        end

        Test.@testset "Construction with explicit overlap_length" begin
            e = Mycelia.Rhizomorph.UltralightQualityEdgeData(10)
            Test.@test isa(e, Mycelia.Rhizomorph.UltralightQualityEdgeData)
            Test.@test e.overlap_length == 10
        end

        Test.@testset "Expected fields present" begin
            e = Mycelia.Rhizomorph.UltralightQualityEdgeData()
            Test.@test hasfield(typeof(e), :overlap_length)
            Test.@test hasfield(typeof(e), :total_count)
            Test.@test hasfield(typeof(e), :dataset_counts)
            Test.@test hasfield(typeof(e), :from_joint_quality)
            Test.@test hasfield(typeof(e), :to_joint_quality)
        end

        Test.@testset "Initial field values" begin
            e = Mycelia.Rhizomorph.UltralightQualityEdgeData()
            Test.@test e.overlap_length == 0
            Test.@test e.total_count == 0
            Test.@test isempty(e.dataset_counts)
            Test.@test e.dataset_counts isa Dict{String, Int}
            Test.@test isempty(e.from_joint_quality)
            Test.@test e.from_joint_quality isa Vector{UInt8}
            Test.@test isempty(e.to_joint_quality)
            Test.@test e.to_joint_quality isa Vector{UInt8}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            e = Mycelia.Rhizomorph.UltralightQualityEdgeData()
            e.total_count = 8
            Test.@test e.total_count == 8
        end

        Test.@testset "Absent fields (no obs tracking in ultralight variants)" begin
            e = Mycelia.Rhizomorph.UltralightQualityEdgeData()
            Test.@test !hasfield(typeof(e), :dataset_observations)
            Test.@test !hasfield(typeof(e), :evidence)
        end
    end

    # ==================================================================
    # LIGHTWEIGHT QUALITY EDGE TYPE
    # ==================================================================

    Test.@testset "LightweightQualityEdgeData" begin
        Test.@testset "Default construction (overlap_length=0)" begin
            e = Mycelia.Rhizomorph.LightweightQualityEdgeData()
            Test.@test isa(e, Mycelia.Rhizomorph.LightweightQualityEdgeData)
        end

        Test.@testset "Construction with explicit overlap_length" begin
            e = Mycelia.Rhizomorph.LightweightQualityEdgeData(3)
            Test.@test isa(e, Mycelia.Rhizomorph.LightweightQualityEdgeData)
            Test.@test e.overlap_length == 3
        end

        Test.@testset "Expected fields present" begin
            e = Mycelia.Rhizomorph.LightweightQualityEdgeData()
            Test.@test hasfield(typeof(e), :overlap_length)
            Test.@test hasfield(typeof(e), :total_count)
            Test.@test hasfield(typeof(e), :dataset_counts)
            Test.@test hasfield(typeof(e), :dataset_observations)
            Test.@test hasfield(typeof(e), :from_joint_quality)
            Test.@test hasfield(typeof(e), :to_joint_quality)
        end

        Test.@testset "Initial field values" begin
            e = Mycelia.Rhizomorph.LightweightQualityEdgeData()
            Test.@test e.overlap_length == 0
            Test.@test e.total_count == 0
            Test.@test isempty(e.dataset_counts)
            Test.@test e.dataset_counts isa Dict{String, Int}
            Test.@test isempty(e.dataset_observations)
            Test.@test e.dataset_observations isa Dict{String, Set{String}}
            Test.@test isempty(e.from_joint_quality)
            Test.@test e.from_joint_quality isa Vector{UInt8}
            Test.@test isempty(e.to_joint_quality)
            Test.@test e.to_joint_quality isa Vector{UInt8}
        end

        Test.@testset "Mutable field total_count can be modified" begin
            e = Mycelia.Rhizomorph.LightweightQualityEdgeData()
            e.total_count = 200
            Test.@test e.total_count == 200
        end

        Test.@testset "Has observation tracking (unlike ultralight variants)" begin
            e = Mycelia.Rhizomorph.LightweightQualityEdgeData()
            # Key differentiator from UltralightQualityEdgeData
            Test.@test hasfield(typeof(e), :dataset_observations)
        end

        Test.@testset "Absent fields (no full evidence dict)" begin
            e = Mycelia.Rhizomorph.LightweightQualityEdgeData()
            Test.@test !hasfield(typeof(e), :evidence)
        end
    end

    # ==================================================================
    # CROSS-TYPE DIFFERENTIATION TESTS
    # ==================================================================

    Test.@testset "Ultralight vs Lightweight Quality field differences" begin
        Test.@testset "Ultralight kmer vertex has no dataset_observations" begin
            uv = Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer)
            Test.@test !hasfield(typeof(uv), :dataset_observations)
        end

        Test.@testset "LightweightQuality kmer vertex has dataset_observations" begin
            lv = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)
            Test.@test hasfield(typeof(lv), :dataset_observations)
        end

        Test.@testset "Ultralight edge has no dataset_observations" begin
            ue = Mycelia.Rhizomorph.UltralightQualityEdgeData()
            Test.@test !hasfield(typeof(ue), :dataset_observations)
        end

        Test.@testset "LightweightQuality edge has dataset_observations" begin
            le = Mycelia.Rhizomorph.LightweightQualityEdgeData()
            Test.@test hasfield(typeof(le), :dataset_observations)
        end
    end

    Test.@testset "Non-quality ultralight types have no quality fields" begin
        vk = Mycelia.Rhizomorph.UltralightKmerVertexData(kmer)
        vb = Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence)
        vs = Mycelia.Rhizomorph.UltralightStringVertexData("ATG")
        ee = Mycelia.Rhizomorph.UltralightEdgeData()

        for v in (vk, vb, vs, ee)
            Test.@test !hasfield(typeof(v), :joint_quality)
            Test.@test !hasfield(typeof(v), :dataset_joint_quality)
        end
        Test.@test !hasfield(typeof(ee), :from_joint_quality)
        Test.@test !hasfield(typeof(ee), :to_joint_quality)
    end

    Test.@testset "Quality types (both ultralight and lightweight) have quality fields" begin
        uvk = Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer)
        uvb = Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence)
        lvk = Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer)
        lvb = Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence)

        for v in (uvk, uvb, lvk, lvb)
            Test.@test hasfield(typeof(v), :joint_quality)
            Test.@test hasfield(typeof(v), :dataset_joint_quality)
        end

        ue = Mycelia.Rhizomorph.UltralightQualityEdgeData()
        le = Mycelia.Rhizomorph.LightweightQualityEdgeData()
        for e in (ue, le)
            Test.@test hasfield(typeof(e), :from_joint_quality)
            Test.@test hasfield(typeof(e), :to_joint_quality)
        end
    end

    Test.@testset "All reduced types start with total_count == 0" begin
        types_and_instances = [
            Mycelia.Rhizomorph.UltralightKmerVertexData(kmer),
            Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.UltralightStringVertexData("ATG"),
            Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer),
            Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer),
            Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.UltralightEdgeData(),
            Mycelia.Rhizomorph.UltralightQualityEdgeData(),
            Mycelia.Rhizomorph.LightweightQualityEdgeData()
        ]
        for instance in types_and_instances
            Test.@test instance.total_count == 0
        end
    end

    Test.@testset "All reduced types start with empty dataset_counts" begin
        types_and_instances = [
            Mycelia.Rhizomorph.UltralightKmerVertexData(kmer),
            Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.UltralightStringVertexData("ATG"),
            Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer),
            Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer),
            Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.UltralightEdgeData(),
            Mycelia.Rhizomorph.UltralightQualityEdgeData(),
            Mycelia.Rhizomorph.LightweightQualityEdgeData()
        ]
        for instance in types_and_instances
            Test.@test isempty(instance.dataset_counts)
        end
    end

    Test.@testset "No reduced type carries a full evidence dict" begin
        types_and_instances = [
            Mycelia.Rhizomorph.UltralightKmerVertexData(kmer),
            Mycelia.Rhizomorph.UltralightBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.UltralightStringVertexData("ATG"),
            Mycelia.Rhizomorph.UltralightQualityKmerVertexData(kmer),
            Mycelia.Rhizomorph.UltralightQualityBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.LightweightQualityKmerVertexData(kmer),
            Mycelia.Rhizomorph.LightweightQualityBioSequenceVertexData(biosequence),
            Mycelia.Rhizomorph.UltralightEdgeData(),
            Mycelia.Rhizomorph.UltralightQualityEdgeData(),
            Mycelia.Rhizomorph.LightweightQualityEdgeData()
        ]
        for instance in types_and_instances
            Test.@test !hasfield(typeof(instance), :evidence)
        end
    end
end

println("Reduced type constructor and field tests executed")
