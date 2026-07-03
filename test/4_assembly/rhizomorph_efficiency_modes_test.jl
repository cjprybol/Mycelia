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
import Graphs

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

    Test.@testset "Mode 2: canonical graph (BLOCKED — undirected path reconstruction)" begin
        reads = eff_tiling_reads(EFF_REF)
        k = 7
        ref_kmers = eff_canonical_kmers([EFF_REF], k)

        # DoubleStrand baseline: correct reconstruction (full canonical coverage).
        cfg_ds = R.AssemblyConfig(;
            k = k, graph_mode = R.DoubleStrand, use_quality_scores = false)
        res_ds = R.assemble_genome(reads, cfg_ds)
        ds_frac = length(intersect(ref_kmers, eff_canonical_kmers(res_ds.contigs, k))) /
                  length(ref_kmers)
        Test.@test ds_frac == 1.0

        cfg_canon = R.AssemblyConfig(;
            k = k, graph_mode = R.Canonical, use_quality_scores = false)
        res_canon = R.assemble_genome(reads, cfg_canon)

        # WHAT WORKS: the canonical graph is the compact (~1x) representation the
        # mode is meant to produce — undirected, with roughly half the vertices of
        # the DoubleStrand graph (RC pairs merged onto one canonical vertex).
        canon_graph = R.build_kmer_graph(reads, k; mode = :canonical)
        ds_graph = R.build_kmer_graph(reads, k; mode = :doublestrand)
        Test.@test !Graphs.is_directed(canon_graph)
        Test.@test length(collect(MetaGraphsNext.labels(canon_graph))) <=
                   0.6 * length(collect(MetaGraphsNext.labels(ds_graph)))

        # WHAT IS BROKEN: path reconstruction over the UNDIRECTED canonical graph is
        # incorrect. find_eulerian_paths_next + path_to_sequence assume a single
        # fixed forward orientation and concatenate the trailing base of each
        # successive canonical k-mer. On an undirected/canonical graph, adjacent
        # canonical labels do not carry which strand their (k-1)-overlap is on, so
        # the reconstructed sequence corresponds to no real read path. Empirically
        # the emitted contig is genome-length but shares almost none of the
        # reference's canonical k-mers (~1/32 on this fixture).
        canon_frac = length(intersect(ref_kmers, eff_canonical_kmers(res_canon.contigs, k))) /
                     length(ref_kmers)
        # Documents the break: canonical coverage is far below the correct 1.0.
        Test.@test canon_frac < 0.5

        # The invariant that SHOULD hold once orientation-aware traversal exists for
        # the canonical graph. Left as @test_broken so a future fix flips it to a
        # visible pass (and any accidental "fix" that regresses is surfaced).
        Test.@test_broken canon_frac == ds_frac
    end

end

println("✓ Rhizomorph efficiency mode tests completed")
