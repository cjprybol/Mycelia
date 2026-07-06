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

    Test.@testset "Mode 2: canonical graph (orientation-aware reconstruction)" begin
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
        # Canonical assembly now reconstructs correctly and must NOT warn about
        # invalid reconstruction (the orientation-aware traversal is in place).
        res_canon = Test.@test_logs min_level = Logging.Warn Mycelia.Rhizomorph.assemble_genome(reads, cfg_canon)

        # And it must flag the result as a VALID reconstruction.
        Test.@test res_canon.assembly_stats["reconstruction_valid"] == true

        # The canonical graph is the compact (~1x) representation the mode is meant
        # to produce — undirected, with roughly half the vertices of the
        # DoubleStrand graph (RC pairs merged onto one canonical vertex).
        canon_graph = Mycelia.Rhizomorph.build_kmer_graph(reads, k; mode = :canonical)
        ds_graph = Mycelia.Rhizomorph.build_kmer_graph(reads, k; mode = :doublestrand)
        Test.@test !Graphs.is_directed(canon_graph)
        Test.@test length(collect(MetaGraphsNext.labels(canon_graph))) <=
                   0.6 * length(collect(MetaGraphsNext.labels(ds_graph)))

        # FIXED: path reconstruction over the UNDIRECTED canonical graph is now
        # orientation-aware. find_eulerian_paths_next handles undirected graphs and
        # path_to_sequence recovers each canonical k-mer's orientation from the
        # (k-1) overlap with its predecessor, reverse-complementing where the
        # overlap is on the reverse strand. Canonical coverage therefore reaches
        # the correct 1.0, matching DoubleStrand.
        canon_frac = length(intersect(ref_kmers, eff_canonical_kmers(res_canon.contigs, k))) /
                     length(ref_kmers)
        Test.@test canon_frac == 1.0

        # The invariant that orientation-aware canonical traversal must satisfy:
        # canonical genome-fraction matches the DoubleStrand baseline.
        Test.@test canon_frac == ds_frac
    end

    Test.@testset "Mode 2b: canonical repeat / RC-overlap (orientation disambiguation)" begin
        # EFF_REF above has no reverse-complement-symmetric junctions at k=7, so a
        # naive sequence-derived reconstruction (pick the first orientation whose
        # (k-1) prefix overlaps the predecessor) happens to succeed there. This
        # fixture DOES contain such junctions: at k=6 the dinucleotide run and the
        # TTGGCAAT / palindromic region create reverse-complement-symmetric overlaps
        # where BOTH orientations of the next canonical k-mer overlap the previous
        # oriented k-mer. Greedy sequence-derivation silently picks the wrong strand
        # there and corrupts the contig (empirically it recovers only ~0.75 of the
        # reference canonical k-mers). Orientation-aware reconstruction consumes the
        # STORED strand evidence the canonical builder recorded to break the tie, so
        # the emitted contig recovers the full reference (canon_frac == 1.0).
        rc_ref = "TGTGTGTTGGCAATAAACCGCCATCCGACTGAGCGCC"
        k = 6
        reads = eff_tiling_reads(rc_ref)
        ref_kmers = eff_canonical_kmers([rc_ref], k)

        cfg_ds = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = false)
        res_ds = Mycelia.Rhizomorph.assemble_genome(reads, cfg_ds)
        ds_frac = length(intersect(ref_kmers, eff_canonical_kmers(res_ds.contigs, k))) /
                  length(ref_kmers)
        Test.@test ds_frac == 1.0

        cfg_canon = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.Canonical, use_quality_scores = false)
        # Orientation-aware canonical reconstruction must succeed WITHOUT any warning
        # (no invalid-reconstruction warning, and no unresolved-overlap fallback).
        res_canon = Test.@test_logs min_level = Logging.Warn Mycelia.Rhizomorph.assemble_genome(reads, cfg_canon)
        Test.@test res_canon.assembly_stats["reconstruction_valid"] == true

        # The single canonical eulerian path exists (path-finding is not the failure
        # mode being exercised here) and every vertex on it carries an unambiguous
        # stored strand, so the disambiguation has a clean authoritative signal.
        canon_graph = Mycelia.Rhizomorph.build_kmer_graph(reads, k; mode = :canonical)
        canon_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(canon_graph)
        Test.@test length(canon_paths) == 1

        # The payload assertion: canonical reconstruction recovers the FULL reference
        # (which greedy sequence-derivation cannot), matching DoubleStrand.
        canon_frac = length(intersect(ref_kmers, eff_canonical_kmers(res_canon.contigs, k))) /
                     length(ref_kmers)
        Test.@test canon_frac == 1.0
        Test.@test canon_frac == ds_frac
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

    Test.@testset "Mode 3a: reduced profiles resolve branches (regression)" begin
        # The clean-tiling fixtures above never hit a branch vertex, so they cannot
        # detect the reduced-profile fragmentation bug (contigs.jl `_edge_support`
        # returned a constant 1 for Lightweight/Ultralight edges, which lack an
        # `evidence` field — collapsing the dominance test so `find_linear_path`
        # broke at every branch). Inject low-coverage "error" tips into a
        # high-depth tiling to create interior branches with UNEQUAL edge support,
        # then require the reduced profiles to reconstruct as well as :full.
        ref = EFF_REF
        k = 7
        reads = FASTX.FASTA.Record[]
        idx = 0
        for _ in 1:8, i in 1:(length(ref) - 15 + 1)
            idx += 1
            push!(reads, FASTX.FASTA.Record("r$(idx)", ref[i:(i + 15 - 1)]))
        end
        for i in (8, 15, 22)
            sub = collect(ref[i:(i + 15 - 1)])
            sub[8] = sub[8] == 'A' ? 'C' : 'A'
            idx += 1
            push!(reads, FASTX.FASTA.Record("e$(idx)", String(sub)))
        end

        cfg_full = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand,
            use_quality_scores = false, memory_profile = :full)
        res_full = Mycelia.Rhizomorph.assemble_genome(reads, cfg_full)
        full_longest = maximum(length.(res_full.contigs))

        for profile in (:lightweight, :ultralight)
            cfg_r = Mycelia.Rhizomorph.AssemblyConfig(;
                k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand,
                use_quality_scores = false, memory_profile = profile)
            res_r = Mycelia.Rhizomorph.assemble_genome(reads, cfg_r)
            r_longest = maximum(length.(res_r.contigs))
            # Reduced profiles must resolve branches, not shatter into k-mer-sized
            # fragments. Pre-fix the longest reduced contig was ~1/3 of :full.
            Test.@test r_longest >= 0.6 * full_longest
        end
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
        # Compaction is now effective (convert_fixed_to_variable +
        # collapse_linear_chains!), so no ineffective-compaction warning is emitted.
        res_compact = Test.@test_logs min_level = Logging.Warn Mycelia.Rhizomorph.assemble_genome(reads, cfg_compact)

        # Contract: compaction must not change the assembled contig sequences (the
        # contigs come from the k-mer traversal, not from simplified_graph).
        Test.@test Set(res_compact.contigs) == Set(res_plain.contigs)

        # When requested, simplified_graph is populated; default leaves it nothing.
        Test.@test res_compact.simplified_graph !== nothing
        Test.@test res_plain.simplified_graph === nothing

        # Keystone: convert_fixed_to_variable() relabels each k-mer vertex as a
        # variable-length BioSequence vertex, so collapse_linear_chains! can now
        # merge non-branching runs into unitigs. A linear tiling collapses to
        # strictly fewer vertices than the full k-mer graph.
        n_full = length(collect(MetaGraphsNext.labels(res_compact.graph)))
        n_simpl = length(collect(MetaGraphsNext.labels(res_compact.simplified_graph)))
        Test.@test n_simpl < n_full

        # The stats disclose that compaction was requested AND effective.
        Test.@test res_compact.assembly_stats["unitig_compaction_requested"] == true
        Test.@test res_compact.assembly_stats["unitig_compaction_effective"] == true
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

    Test.@testset "Canonical qualmer: orientation-aware sequence AND quality" begin
        # The quality (qualmer) arm reconstructs canonical contigs through the SAME
        # orientation-aware path as the k-mer arm, and per-base quality is oriented
        # to match. Exercise it on the RC-overlap fixture with a STRICTLY MONOTONIC
        # per-position quality gradient: because each read encodes the same quality
        # for a given reference position, the correct per-position mean quality is
        # itself monotonic along the reference. If any reverse-oriented k-mer emitted
        # the wrong-end quality (the pre-fix bug used vertex_quality[end] for a
        # reverse k-mer, whose emitted base is actually complement(first base) ->
        # vertex_quality[1]), the output quality would be scrambled and NON-monotonic.
        rc_ref = "TGTGTGTTGGCAATAAACCGCCATCCGACTGAGCGCC"
        k = 6
        read_len = 15
        qual_for(pos) = Char(min(20 + pos, 60) + 33)  # strictly increasing along ref
        fastq_reads = FASTX.FASTQ.Record[]
        for i in 1:(length(rc_ref) - read_len + 1)
            sub = rc_ref[i:(i + read_len - 1)]
            qs = join(qual_for(i + j - 1) for j in 1:read_len)
            push!(fastq_reads, FASTX.FASTQ.Record("r$(i)", sub, qs))
        end
        ref_kmers = eff_canonical_kmers([rc_ref], k)

        cfg_ds = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.DoubleStrand, use_quality_scores = true)
        res_ds = Mycelia.Rhizomorph.assemble_genome(fastq_reads, cfg_ds)
        ds_frac = length(intersect(ref_kmers, eff_canonical_kmers(res_ds.contigs, k))) /
                  length(ref_kmers)
        Test.@test ds_frac == 1.0

        cfg_canon = Mycelia.Rhizomorph.AssemblyConfig(;
            k = k, graph_mode = Mycelia.Rhizomorph.Canonical, use_quality_scores = true)
        # Canonical qualmer must reconstruct correctly and emit NO warning (no
        # invalid-reconstruction warning, no unresolved-overlap or length-clamp warn).
        res_canon = Test.@test_logs min_level = Logging.Warn Mycelia.Rhizomorph.assemble_genome(fastq_reads, cfg_canon)
        Test.@test res_canon.assembly_stats["reconstruction_valid"] == true

        # Sequence correctness: full canonical coverage, matching DoubleStrand.
        canon_frac = length(intersect(ref_kmers, eff_canonical_kmers(res_canon.contigs, k))) /
                     length(ref_kmers)
        Test.@test canon_frac == 1.0
        Test.@test canon_frac == ds_frac

        # Rebuild the canonical qualmer contig record to inspect sequence/quality.
        fq = Mycelia.Rhizomorph._prepare_fastq_observations(fastq_reads)
        gc = Mycelia.Rhizomorph.build_qualmer_graph(fq, k; mode = :canonical)
        canon_paths = Mycelia.Rhizomorph.find_eulerian_paths_next(gc)
        Test.@test length(canon_paths) == 1
        path = canon_paths[argmax(length.(canon_paths))]

        # The orientation-aware path MUST actually reverse-complement some k-mers,
        # otherwise the quality-orientation logic would be untested.
        flags = Mycelia.Rhizomorph._resolve_path_orientations(path, gc)
        Test.@test count(!, flags) > 0

        rec = Mycelia.Rhizomorph._qualmer_path_to_consensus_fastq(path, gc, "canon")
        seq = String(FASTX.sequence(rec))
        qual = collect(FASTX.quality_scores(rec))

        # Sequence == reference (or its reverse complement) and quality/sequence stay
        # length-consistent (the length clamp must never shorten a correct contig).
        Test.@test seq == rc_ref || seq == eff_rc(rc_ref)
        Test.@test length(qual) == length(seq)

        # FIX 3 payload: the per-base quality tracks the emitted base's orientation,
        # so the monotonic reference-quality gradient survives (ascending if the
        # contig came out forward, descending if reverse-complemented). A misoriented
        # quality vector would break this.
        Test.@test issorted(qual) || issorted(qual; rev = true)
    end

end

println("✓ Rhizomorph efficiency mode tests completed")
