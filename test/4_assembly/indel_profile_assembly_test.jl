# Step 4a (indel-wiring): drive the merged indel pair-HMM decoder from an ERROR
# PROFILE (sequencing technology), not read length or a hardcoded tech. The
# indel kernel (`_viterbi_correct_observation_indel`) is reachable only via
# `correct_observations`; this suite proves `assemble_genome` now enables it for
# indel-prone technologies (nanopore/pacbio) while keeping the DEFAULT
# (:illumina) on the substitution-only oracle path byte-for-byte.
#
# The profile -> fractions -> indel_moves mapping and the threshold gate are pure
# functions (asserted without running the slow corrector); two end-to-end
# assemblies prove the nanopore path enables indels and the illumina path stays
# byte-identical to a run with no new parameter.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/indel_profile_assembly_test.jl")'

import Test
import Mycelia
import FASTX
import BioSequences
import Random
import Statistics

# Nanopore/illumina read set sampled from a random reference through
# `Mycelia.observe(...)` at the requested tech (returns a (seq, quals) TUPLE).
# Returns `(records, refseq)` so callers can score assemblies against the truth.
function _profile_reads(; tech::Symbol = :nanopore, genome_len = 1000, readlen = 100,
        coverage = 20, err = 0.05, seed = 42)
    rec = Mycelia.random_fasta_record(moltype = :DNA, seed = seed, L = genome_len)
    refseq = FASTX.sequence(BioSequences.LongDNA{4}, rec)
    rng = Random.MersenneTwister(seed)
    glen = length(refseq)
    n_reads = max(1, ceil(Int, coverage * glen / readlen))
    records = FASTX.FASTQ.Record[]
    for i in 1:n_reads
        start = rand(rng, 1:(glen - readlen + 1))
        frag = refseq[start:(start + readlen - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        obs = Mycelia.observe(frag; error_rate = err, tech = tech)
        obs_seq = obs[1]
        quals = obs[2]
        isempty(obs_seq) && continue
        qstr = String([Char(q + 33) for q in quals])
        push!(records, FASTX.FASTQ.Record("read_$(i)", string(obs_seq), qstr))
    end
    return records, refseq
end

# Best Levenshtein identity of any assembled contig against the reference, in
# either orientation (assembly is strand-agnostic). 0.0 when there are no contigs.
function _best_identity_to_ref(contigs::Vector{String}, refseq)
    isempty(contigs) && return 0.0
    ref = string(refseq)
    rc = string(BioSequences.reverse_complement(BioSequences.LongDNA{4}(ref)))
    best = 0.0
    for c in contigs
        for target in (ref, rc)
            r = Mycelia.assess_alignment(c, target)
            best = max(best, r.total_matches / (r.total_matches + r.total_edits))
        end
    end
    return best
end

Test.@testset "indel-aware correction wired via sequencing-tech error profile" begin
    R = Mycelia.Rhizomorph

    Test.@testset "error-profile presets mirror simulation rates + fractions (pure)" begin
        np = Mycelia.indel_error_profile(:nanopore)
        # Nanopore ABSOLUTE base rate + conditional insertion/deletion fractions
        # mirror observe()'s per-tech values (base 0.10, split 0.30 / 0.30 near
        # simulation.jl:1928).
        Test.@test np.base_error_rate ≈ 0.10
        Test.@test np.insertion_fraction ≈ 0.30
        Test.@test np.deletion_fraction ≈ 0.30
        Test.@test np.insertion_extend_probability > 0.0
        Test.@test np.deletion_extend_probability > 0.0

        il = Mycelia.indel_error_profile(:illumina)
        # Illumina carries its TRUE conditional fractions (mirroring observe():
        # base 0.005, split 0.05 / 0.05) — NOT hand-zeroed. Its ABSOLUTE indel rate
        # (base × summed fractions) is what is negligible.
        Test.@test il.base_error_rate ≈ 0.005
        Test.@test il.insertion_fraction ≈ 0.05
        Test.@test il.deletion_fraction ≈ 0.05
        Test.@test il.base_error_rate * (il.insertion_fraction + il.deletion_fraction) <
                   0.02

        pb = Mycelia.indel_error_profile(:pacbio)
        # PacBio: base 0.11, split 0.40 / 0.40 ⇒ absolute indel rate ≈ 0.088.
        Test.@test pb.base_error_rate ≈ 0.11
        Test.@test pb.base_error_rate * (pb.insertion_fraction + pb.deletion_fraction) >
                   0.02

        ul = Mycelia.indel_error_profile(:ultima)
        # Ultima carries its TRUE conditional fractions (observe(): base 1e-6,
        # split 0.025 / 0.025); the absolute indel rate is what is negligible.
        Test.@test ul.insertion_fraction ≈ 0.025
        Test.@test ul.deletion_fraction ≈ 0.025
        Test.@test ul.base_error_rate * (ul.insertion_fraction + ul.deletion_fraction) <
                   0.02

        Test.@test_throws Exception Mycelia.indel_error_profile(:bogus)
    end

    Test.@testset "threshold gate: profile -> indel_moves" begin
        # Above threshold => indel-aware; negligible => substitution-only.
        Test.@test Mycelia.profile_enables_indels(:nanopore) == true
        Test.@test Mycelia.profile_enables_indels(:pacbio) == true
        Test.@test Mycelia.profile_enables_indels(:illumina) == false
        Test.@test Mycelia.profile_enables_indels(:ultima) == false
        # The gate is a strict threshold on the ABSOLUTE indel rate
        # (base_error_rate × summed fractions), not the conditional fractions.
        # Nanopore's absolute rate ≈ 0.06 sits below a 0.7 threshold.
        Test.@test Mycelia.profile_enables_indels(:nanopore; threshold = 0.7) == false
        # Illumina's absolute rate ≈ 5e-4: OFF at the default 0.02, but ON once the
        # threshold drops below it (1e-5) — proving the gate keys on the absolute
        # product, not the (nonzero) conditional fractions.
        Test.@test Mycelia.profile_enables_indels(:illumina; threshold = 1e-5) == true
    end

    Test.@testset "AssemblyConfig threads + validates sequencing_tech" begin
        # Default is :illumina (unchanged behavior).
        Test.@test R.AssemblyConfig(k = 13).sequencing_tech == :illumina
        Test.@test R.AssemblyConfig(k = 13, sequencing_tech = :nanopore).sequencing_tech ==
                   :nanopore
        Test.@test_throws Exception R.AssemblyConfig(k = 13, sequencing_tech = :bogus)
    end

    Test.@testset "tier knobs carry indel bounds" begin
        sc = R._corrector_strategy_knobs(:scalable)
        # :scalable => bounded band + small run caps.
        Test.@test sc.band_width isa Int
        Test.@test sc.band_width > 0
        Test.@test sc.deletion_max_run == 3
        Test.@test sc.max_insertion_run >= 1

        ex = R._corrector_strategy_knobs(:exhaustive)
        # :exhaustive => unbounded band + larger caps.
        Test.@test ex.band_width === nothing
        Test.@test ex.deletion_max_run > sc.deletion_max_run
    end

    Test.@testset "illumina default is the substitution-only oracle (byte-identical)" begin
        reads, _ = _profile_reads(tech = :illumina, err = 0.01, seed = 7)
        base = R.assemble_genome(reads; k = 13, corrector = :iterative)
        il = R.assemble_genome(reads; k = 13, corrector = :iterative,
            sequencing_tech = :illumina)
        # Oracle guard (PR #408 review, FIX 4): the byte-identity check below is
        # only meaningful if the corrector actually produced contigs — otherwise
        # `[] == []` passes vacuously and hides a corrector that ate every read.
        Test.@test !isempty(base.contigs)
        # Contigs byte-identical: the illumina profile threads NO indel params, so
        # the corrector runs the substitution decode unchanged.
        Test.@test il.contigs == base.contigs
        Test.@test base.assembly_stats["indel_moves"] == false
        Test.@test il.assembly_stats["indel_moves"] == false
        Test.@test il.assembly_stats["sequencing_tech"] == "illumina"
    end

    Test.@testset "nanopore enables indel-aware correction" begin
        reads, _ = _profile_reads(tech = :nanopore, err = 0.05, seed = 7)
        np = R.assemble_genome(reads; k = 13, corrector = :iterative,
            sequencing_tech = :nanopore)
        Test.@test np isa R.AssemblyResult
        Test.@test !isempty(np.contigs)
        # The provenance stamp proves the indel decode path was engaged.
        Test.@test np.assembly_stats["indel_moves"] == true
        Test.@test np.assembly_stats["sequencing_tech"] == "nanopore"
        Test.@test np.assembly_stats["indel_decodes"] >= 0
        Test.@test np.assembly_stats["window_divergences"] >= 0
    end

    # -- END-TO-END correction proof (PR #408 review, FIX 2) -------------------
    # The stamp checks above are config-derived tautologies: `indel_moves==true`
    # only proves the config enabled the gap moves, not that assemble_genome
    # actually CORRECTED an indel. This testset proves the value of the wiring end
    # to end and would FAIL if any of the threading forwards (profile -> gate ->
    # IndelDecodeParams -> base_error_rate -> ViterbiCorrectionConfig.error_rate)
    # were dropped: it correct-assembles indel-prone nanopore reads and shows the
    # indel-aware arm is materially CLOSER to the truth than the substitution-only
    # (:illumina) arm run on the SAME reads.
    #
    # Regime: a short (60 bp) non-repetitive reference, PARTIAL nanopore reads
    # (40 bp) at 30x, err=0.05. Partial reads must be tiled through the graph, so
    # uncorrected indels frameshift the overlaps and fragment / degrade the
    # substitution-only assembly — exactly where indel-aware correction pays off.
    # The outcome is stochastic per seed, so we aggregate over several seeds and
    # assert the DIFFERENTIAL on the MEAN identity (self-calibrating — no magic
    # threshold) plus a majority-of-seeds win and a recovery floor.
    Test.@testset "nanopore indel-aware correction beats substitution-only (E2E)" begin
        seeds = 1:5
        id_np = Float64[]
        id_il = Float64[]
        np_indel_moves = Bool[]
        np_indel_decodes = Int[]
        np_nonempty = Bool[]
        for seed in seeds
            reads,
            refseq = _profile_reads(tech = :nanopore, genome_len = 60,
                readlen = 40, coverage = 30, err = 0.05, seed = seed)
            np = R.assemble_genome(reads; k = 17, corrector = :iterative,
                sequencing_tech = :nanopore)
            il = R.assemble_genome(reads; k = 17, corrector = :iterative,
                sequencing_tech = :illumina)
            push!(id_np, _best_identity_to_ref(np.contigs, refseq))
            push!(id_il, _best_identity_to_ref(il.contigs, refseq))
            push!(np_indel_moves, np.assembly_stats["indel_moves"] == true)
            push!(np_indel_decodes, np.assembly_stats["indel_decodes"])
            push!(np_nonempty, !isempty(np.contigs))
        end
        mean_np = Statistics.mean(id_np)
        mean_il = Statistics.mean(id_il)
        wins = count(id_np .>= id_il)
        @info "E2E nanopore-vs-illumina identity" seeds = collect(seeds) id_np id_il mean_np mean_il wins

        # The nanopore arm actually engaged the indel decode and produced contigs.
        Test.@test all(np_indel_moves)
        Test.@test all(>(0), np_indel_decodes)
        Test.@test all(np_nonempty)
        # DIFFERENTIAL (the PR's value): indel-aware correction lands closer to the
        # truth than substitution-only on the SAME nanopore reads, on average.
        Test.@test mean_np > mean_il
        # Nanopore wins or ties the majority of seeds (corroborates the mean).
        Test.@test wins >= 3
        # Nanopore recovers most of the reference (proves indels were corrected, not
        # merely that the arm ran) — a substantial-identity floor with margin.
        Test.@test mean_np > 0.80
    end
end
