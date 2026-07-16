# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_qualmer_aggregate_quality_test.jl")'
# ```
#
# Tests for corrector-compatible AGGREGATE qualmer quality storage (td-n8ax):
# a coverage counter + running per-position MEAN Phred per DISTINCT k-mer,
# O(distinct k-mers), composed with the opt-in coverage prefilter (PR #425), so
# E. coli-scale assembly stays under a SPAdes-class memory ceiling.
#
# Every behavior is OPT-IN: the default corrector path stays memory_profile=:full
# (per-observation evidence) and byte-identical. These testsets fail until the
# aggregate storage is made corrector-compatible (B1 accessor + B2 exact-mean
# storage + B4 mode/prefilter + wiring).

import Test
import Mycelia
import FASTX
import BioSequences
import MetaGraphsNext
import Graphs

# ---------------------------------------------------------------------------
# Shared fixtures / helpers
# ---------------------------------------------------------------------------

"A short reference with no exact repeats at k=7, so k-mers are distinct."
const AGG_REF = "ATCGGCTAATGCCGATTGCACGTACGTTAGCTAGGCATG"

"""
Build `n` identical FASTQ reads of `seq` with a strictly-ascending per-position
Phred so each k-mer's correct per-position mean is non-uniform (catches position
ordering) AND high enough that a running SUM over `n` observations saturates a
UInt8 accumulator (raw Phred ~40 * n>=7 > 255). This is what discriminates an
exact unclamped mean (correct) from the clamped joint-quality sum (wrong).
"""
function agg_identical_fastq_reads(seq::AbstractString; n::Int = 10)
    qual = join(Char(min(35 + j, 60) + 33) for j in 1:length(seq))  # ~Q36..Q60
    reads = FASTX.FASTQ.Record[]
    for i in 1:n
        push!(reads, FASTX.FASTQ.Record("r$(i)", seq, qual))
    end
    return reads
end

Test.@testset "Rhizomorph aggregate qualmer quality (td-n8ax)" begin

    # -----------------------------------------------------------------------
    # A. EXACT-MEAN EQUIVALENCE: get_vertex_mean_quality on the aggregate
    #    (:lightweight_quality) graph == the :full graph, per position, even at
    #    coverage that saturates a UInt8 running sum. Proves B1 (accessor) + B2
    #    (unclamped storage).
    # -----------------------------------------------------------------------
    Test.@testset "A: exact per-position mean quality matches :full" begin
        k = 7
        reads = agg_identical_fastq_reads(AGG_REF; n = 10)
        fq = Mycelia.Rhizomorph._prepare_fastq_observations(reads)

        full_graph = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :singlestrand, memory_profile = :full)
        agg_graph = Mycelia.Rhizomorph.build_qualmer_graph(
            fq, k; mode = :singlestrand, memory_profile = :lightweight_quality)

        full_labels = Set(MetaGraphsNext.labels(full_graph))
        agg_labels = Set(MetaGraphsNext.labels(agg_graph))
        # Same distinct k-mers (as-observed keys) in both storage modes.
        Test.@test full_labels == agg_labels
        Test.@test !isempty(agg_labels)

        ds = "dataset_01"
        compared = 0
        for label in full_labels
            full_mean = Mycelia.Rhizomorph.get_vertex_mean_quality(full_graph[label], ds)
            agg_mean = Mycelia.Rhizomorph.get_vertex_mean_quality(agg_graph[label], ds)
            # The aggregate MUST expose a real mean, not the nothing/UInt8(2) filler.
            Test.@test agg_mean !== nothing
            Test.@test full_mean !== nothing
            Test.@test length(agg_mean) == k
            # Numerically equal to the per-observation mean (exact, not clamped).
            Test.@test all(isapprox.(agg_mean, full_mean; atol = 1e-9))
            compared += 1
        end
        Test.@test compared == length(full_labels)
    end
end

println("✓ Rhizomorph aggregate qualmer quality tests completed")
