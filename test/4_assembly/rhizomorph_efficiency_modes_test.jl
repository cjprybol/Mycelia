# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/rhizomorph_efficiency_modes_test.jl")'
# ```
#
# Tests for the three additive assembly-efficiency modes:
#   Mode 1: reverse-complement contig dedup   (config.dedup_revcomp)
#   Mode 2: canonical graph mode              (graph_mode = Canonical)
#   Mode 3a: memory_profile threading         (config.memory_profile)
#   Mode 3b: unitig compaction                (config.compact_unitigs)
#
# Every mode is opt-in; the default AssemblyConfig must reproduce today's
# DoubleStrand / :full / no-dedup / no-compaction behavior exactly.

import Test
import Mycelia
import FASTX
import BioSequences
import MetaGraphsNext

const R = Mycelia.Rhizomorph

# ---------------------------------------------------------------------------
# Shared fixtures / helpers
# ---------------------------------------------------------------------------

"A small reference with no exact repeats at k=7, so assembly is well-defined."
const EFF_REF = "ATCGGCTAATGCCGATTGCACGTACGTTAGCTAGGCATG"

"Tile the reference into overlapping FASTA reads (no quality => k-mer graph arm)."
function eff_tiling_reads(ref::AbstractString; read_len::Int = 15)
    reads = FASTX.FASTA.Record[]
    for i in 1:(length(ref) - read_len + 1)
        push!(reads, FASTX.FASTA.Record("r$(i)", ref[i:(i + read_len - 1)]))
    end
    return reads
end

function eff_rc(seq::AbstractString)::String
    return String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(String(seq))))
end

"Canonical set of k-mers (as strings) for a collection of sequences."
function eff_canonical_kmers(seqs, k::Int)
    kmers = Set{String}()
    for s in seqs
        str = String(s)
        for i in 1:(length(str) - k + 1)
            sub = str[i:(i + k - 1)]
            push!(kmers, min(sub, eff_rc(sub)))
        end
    end
    return kmers
end

"True if any two distinct contigs in the set are reverse complements of each other."
function eff_has_rc_pair(contigs)
    s = Set(String.(contigs))
    for c in contigs
        rc = eff_rc(c)
        if rc != c && rc in s
            return true
        end
    end
    return false
end

Test.@testset "Rhizomorph efficiency modes" begin

    Test.@testset "Mode 1: reverse-complement dedup" begin
        reads = eff_tiling_reads(EFF_REF)
        k = 7

        # (a) Default (dedup off) still emits RC pairs from the DoubleStrand arm.
        cfg_default = R.AssemblyConfig(;
            k = k, graph_mode = R.DoubleStrand, use_quality_scores = false)
        res_default = R.assemble_genome(reads, cfg_default)
        Test.@test !isempty(res_default.contigs)
        Test.@test eff_has_rc_pair(res_default.contigs)

        # (b) With dedup on, contig count drops and no two remaining contigs are RCs.
        cfg_dedup = R.AssemblyConfig(;
            k = k, graph_mode = R.DoubleStrand, use_quality_scores = false,
            dedup_revcomp = true)
        res_dedup = R.assemble_genome(reads, cfg_dedup)
        Test.@test !isempty(res_dedup.contigs)
        Test.@test length(res_dedup.contigs) < length(res_default.contigs)
        Test.@test !eff_has_rc_pair(res_dedup.contigs)
        # No exact duplicates remain either.
        Test.@test length(unique(res_dedup.contigs)) == length(res_dedup.contigs)

        # (c) The deduped set still covers the reference genome content: every
        # canonical reference k-mer is present among the deduped contigs.
        ref_kmers = eff_canonical_kmers([EFF_REF], k)
        contig_kmers = eff_canonical_kmers(res_dedup.contigs, k)
        Test.@test issubset(ref_kmers, contig_kmers)

        # And dedup preserves canonical coverage relative to the full set (no loss).
        default_kmers = eff_canonical_kmers(res_default.contigs, k)
        Test.@test contig_kmers == default_kmers
    end

end

println("✓ Rhizomorph efficiency mode tests completed")
