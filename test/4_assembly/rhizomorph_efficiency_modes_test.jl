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
import Logging

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

"Tile the reference into overlapping FASTQ reads (quality => qualmer graph arm)."
function eff_tiling_fastq_reads(ref::AbstractString; read_len::Int = 15)
    reads = FASTX.FASTQ.Record[]
    for i in 1:(length(ref) - read_len + 1)
        sub = ref[i:(i + read_len - 1)]
        qual = repeat("I", length(sub))  # Phred+33 'I' == Q40
        push!(reads, FASTX.FASTQ.Record("r$(i)", sub, qual))
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
        cfg_default = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false)
        res_default = Mycelia.Rhizomorph.assemble_genome(reads, cfg_default)
        Test.@test !isempty(res_default.contigs)
        Test.@test eff_has_rc_pair(res_default.contigs)

        # (b) With dedup on, contig count drops and no two remaining contigs are RCs.
        cfg_dedup = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false,
            dedup_revcomp = true)
        res_dedup = Mycelia.Rhizomorph.assemble_genome(reads, cfg_dedup)
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
        cfg_ds = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false)
        res_ds = Mycelia.Rhizomorph.assemble_genome(reads, cfg_ds)
        ds_frac = length(intersect(ref_kmers, eff_canonical_kmers(res_ds.contigs, k))) /
                  length(ref_kmers)
        Test.@test ds_frac == 1.0

        # DoubleStrand reconstruction is flagged valid.
        Test.@test res_ds.assembly_stats["reconstruction_valid"] == true

        cfg_canon = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.Canonical, use_quality_scores = false)
        # FIX 1: Canonical assembly must warn UNCONDITIONALLY (not verbose-gated)
        # that contig reconstruction is not yet correct.
        res_canon = Test.@test_logs (:warn,) match_mode = :any Mycelia.Rhizomorph.assemble_genome(reads, cfg_canon)

        # FIX 1: and it must flag the result as an invalid reconstruction.
        Test.@test res_canon.assembly_stats["reconstruction_valid"] == false

        # WHAT WORKS: the canonical graph is the compact (~1x) representation the
        # mode is meant to produce — undirected, with roughly half the vertices of
        # the DoubleStrand graph (RC pairs merged onto one canonical vertex).
        canon_graph = Mycelia.Rhizomorph.build_kmer_graph(reads, k; mode = :canonical)
        ds_graph = Mycelia.Rhizomorph.build_kmer_graph(reads, k; mode = :doublestrand)
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

    Test.@testset "Mode 3a: memory_profile threading" begin
        reads = eff_tiling_reads(EFF_REF)
        k = 7

        cfg_full = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false,
            memory_profile = :full)
        res_full = Mycelia.Rhizomorph.assemble_genome(reads, cfg_full)

        cfg_light = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false,
            memory_profile = :lightweight)
        res_light = Mycelia.Rhizomorph.assemble_genome(reads, cfg_light)

        # memory_profile is an internal footprint change, not an output change:
        # the assembled contigs must be identical (order-independent).
        Test.@test Set(res_light.contigs) == Set(res_full.contigs)
        Test.@test length(res_light.contigs) == length(res_full.contigs)

        # Same expectation for the most compact profile.
        cfg_ultra = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false,
            memory_profile = :ultralight)
        res_ultra = Mycelia.Rhizomorph.assemble_genome(reads, cfg_ultra)
        Test.@test Set(res_ultra.contigs) == Set(res_full.contigs)
    end

    Test.@testset "Mode 3b: unitig compaction (partial — k-mer graph caveat)" begin
        reads = eff_tiling_reads(EFF_REF)
        k = 7

        cfg_plain = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false)
        res_plain = Mycelia.Rhizomorph.assemble_genome(reads, cfg_plain)

        cfg_compact = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false,
            compact_unitigs = true)
        # FIX 4: compaction on a fixed-length k-mer graph is a no-op and must warn.
        res_compact = Test.@test_logs (:warn,) match_mode = :any Mycelia.Rhizomorph.assemble_genome(reads, cfg_compact)

        # Contract: compaction must not change the assembled contig sequences.
        Test.@test Set(res_compact.contigs) == Set(res_plain.contigs)

        # When requested, simplified_graph is populated; default leaves it nothing.
        Test.@test res_compact.simplified_graph !== nothing
        Test.@test res_plain.simplified_graph === nothing

        # Documented caveat: collapse_linear_chains! is a no-op on fixed-length
        # k-mer graphs (collapsing would change the label type from Kmer to
        # BioSequence). So the simplified graph currently has the SAME vertex count
        # as the full graph — no real compaction happens yet.
        n_full = length(collect(MetaGraphsNext.labels(res_compact.graph)))
        n_simpl = length(collect(MetaGraphsNext.labels(res_compact.simplified_graph)))
        Test.@test n_simpl == n_full

        # FIX 4: the stats disclose that compaction was requested but ineffective.
        Test.@test res_compact.assembly_stats["unitig_compaction_requested"] == true
        Test.@test res_compact.assembly_stats["unitig_compaction_effective"] == false

        # The invariant that SHOULD hold once fixed->variable conversion is wired in
        # (a linear tiling should compact to strictly fewer vertices). Left as
        # @test_broken to flag the known gap without red-failing the suite.
        Test.@test_broken n_simpl < n_full
    end

    Test.@testset "Config validation guards" begin
        # (a) An unrecognized memory_profile is rejected at construction.
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(;
            k = 7, use_quality_scores = false, memory_profile = :bogus)

        # (b) Canonical (like DoubleStrand) requires a defined reverse complement,
        # so it is rejected for amino-acid and general-string sequence types.
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(;
            k = 7, sequence_type = BioSequences.LongAA,
            graph_mode = Mycelia.Rhizomorph.Canonical)
        Test.@test_throws ErrorException Mycelia.Rhizomorph.AssemblyConfig(;
            k = 7, sequence_type = String,
            graph_mode = Mycelia.Rhizomorph.Canonical)
    end

    Test.@testset "Mode 1b: structural RC-dedup in find_contigs_next" begin
        # find_contigs_next(; rc_aware=true) marks each walked vertex's
        # reverse-complement partner visited, so the DoubleStrand graph's reverse
        # strand is never independently traversed. This is a STRUCTURAL fix for the
        # QUAST-duplication-~2.0 pathology: unlike post-hoc whole-contig string
        # dedup, it removes RC twins even when their fragment breakpoints are
        # OFFSET between strands (the empirically-observed case where contig-level
        # string dedup halves the count yet leaves duplication at 2.0).
        reads = eff_tiling_reads(EFF_REF)
        k = 7
        g = Mycelia.Rhizomorph.build_kmer_graph(reads, k; mode = :doublestrand)

        c_off = [string(c.sequence)
                 for c in Mycelia.Rhizomorph.find_contigs_next(g; min_contig_length = 1)]
        c_on = [string(c.sequence)
                for c in Mycelia.Rhizomorph.find_contigs_next(
                    g; min_contig_length = 1, rc_aware = true)]

        # (a) Default (rc_aware=false) is unchanged and still emits RC pairs.
        Test.@test eff_has_rc_pair(c_off)
        # (b) rc_aware strictly reduces the contig count (removes the RC strand).
        Test.@test length(c_on) < length(c_off)
        # (c) No RC twins remain after structural dedup.
        Test.@test !eff_has_rc_pair(c_on)
        # (d) Genome content is PRESERVED: canonical k-mer coverage is unchanged
        #     (accuracy is not traded for the contig-count reduction).
        Test.@test eff_canonical_kmers(c_on, k) == eff_canonical_kmers(c_off, k)
        Test.@test issubset(eff_canonical_kmers([EFF_REF], k), eff_canonical_kmers(c_on, k))

        # (e) The helper is a no-op on labels with no defined reverse complement.
        Test.@test Mycelia.Rhizomorph._rc_partner_label("ACGT") === nothing
    end

    Test.@testset "FIX 2: efficiency flags on the quality/qualmer path warn" begin
        # The efficiency modes are only honored on the non-quality k-mer path.
        # When use_quality_scores=true (the default, and what FASTQ input auto-sets)
        # AND an efficiency flag is requested, construction must warn UNCONDITIONALLY
        # that the flag is a no-op on the quality path.
        Test.@test_logs (:warn,) match_mode = :any Mycelia.Rhizomorph.AssemblyConfig(;
            k = 7, use_quality_scores = true, dedup_revcomp = true)
        Test.@test_logs (:warn,) match_mode = :any Mycelia.Rhizomorph.AssemblyConfig(;
            k = 7, use_quality_scores = true, compact_unitigs = true)
        Test.@test_logs (:warn,) match_mode = :any Mycelia.Rhizomorph.AssemblyConfig(;
            k = 7, use_quality_scores = true, memory_profile = :lightweight)

        # And it must NOT warn for the default (no efficiency flag) quality config —
        # existing quality assemblies stay byte-for-byte unchanged with no noise.
        Test.@test_logs min_level = Logging.Warn Mycelia.Rhizomorph.AssemblyConfig(;
            k = 7, use_quality_scores = true)
    end

end

println("✓ Rhizomorph efficiency mode tests completed")
