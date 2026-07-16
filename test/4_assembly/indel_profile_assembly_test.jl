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
import Kmers
import Random
import Statistics

if !isdefined(Main, :test_throws_message)
    include(joinpath(dirname(@__DIR__), "test_helpers.jl"))
end

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
        # Legacy :pacbio remains the CLR-like compatibility alias.
        Test.@test pb.base_error_rate ≈ 0.11
        Test.@test pb.base_error_rate * (pb.insertion_fraction + pb.deletion_fraction) >
                   0.02
        Test.@test Mycelia.indel_error_profile(:pacbio_clr) == pb

        hifi = Mycelia.indel_error_profile(:pacbio_hifi)
        # HiFi has an exact high-accuracy profile and deliberately stays off the
        # CLR-like indel pair-HMM until a validated HiFi indel split is modeled.
        Test.@test hifi.base_error_rate ≈ 0.001
        Test.@test hifi.insertion_fraction == 0.0
        Test.@test hifi.deletion_fraction == 0.0

        ul = Mycelia.indel_error_profile(:ultima)
        # Ultima carries its TRUE conditional fractions (observe(): base 1e-6,
        # split 0.025 / 0.025); the absolute indel rate is what is negligible.
        Test.@test ul.insertion_fraction ≈ 0.025
        Test.@test ul.deletion_fraction ≈ 0.025
        Test.@test ul.base_error_rate * (ul.insertion_fraction + ul.deletion_fraction) <
                   0.02

        test_throws_message(ErrorException, "unknown sequencing technology") do
            Mycelia.indel_error_profile(:bogus)
        end
    end

    Test.@testset "threshold gate: profile -> indel_moves" begin
        # Above threshold => indel-aware; negligible => substitution-only.
        Test.@test Mycelia.profile_enables_indels(:nanopore) == true
        Test.@test Mycelia.profile_enables_indels(:pacbio) == true
        Test.@test Mycelia.profile_enables_indels(:pacbio_clr) == true
        Test.@test Mycelia.profile_enables_indels(:pacbio_hifi) == false
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

    Test.@testset "substitution fallback-rate config + forwarding" begin
        config_options = (;
            alphabet = :DNA,
            graph_mode = :canonical,
            max_steps = 7,
            beam_width = typemax(Int),
            max_successors_per_state = typemax(Int),
            beam_score_margin = Inf,
            record_position_gaps = false,
        )

        # The default-nothing path deliberately omits `error_rate` from the
        # ViterbiCorrectionConfig constructor. Illumina, Ultima, and legacy direct
        # callers therefore retain the historical substitution-only 0.01 fallback.
        legacy = Mycelia._iterative_viterbi_correction_config(; config_options...)
        Test.@test legacy.error_rate == 0.01
        Test.@test legacy.indel_moves == false

        # HiFi remains substitution-only, but its quality-free fallback is the exact
        # technology-profile rate rather than the legacy generic default.
        hifi = Mycelia._iterative_viterbi_correction_config(
            ; config_options..., substitution_error_rate = 0.001)
        Test.@test hifi.error_rate == 0.001
        Test.@test hifi.indel_moves == false
        Test.@test R._substitution_error_rate(:pacbio_hifi) == 0.001
        Test.@test R._substitution_error_rate(:illumina) === nothing
        Test.@test R._substitution_error_rate(:ultima) === nothing

        profile = Mycelia.indel_error_profile(:nanopore)
        indel_params = Mycelia.IndelDecodeParams(
            profile.base_error_rate,
            profile.insertion_fraction,
            profile.deletion_fraction,
            profile.insertion_extend_probability,
            profile.deletion_extend_probability,
            3,
            3,
            16,
        )
        test_throws_message(ArgumentError, "mutually exclusive") do
            Mycelia._iterative_viterbi_correction_config(
                ; config_options..., indel_params,
                substitution_error_rate = 0.001)
        end
        test_throws_message(ArgumentError, "error_rate must be in (0, 0.5)") do
            Mycelia._iterative_viterbi_correction_config(
                ; config_options..., substitution_error_rate = 0.0)
        end

        # Exercise the real non-windowed and windowed forwarding chains. The hard
        # set contains every canonical k-mer in the probe, guaranteeing that the
        # windowed arm reaches its per-window Viterbi config builder.
        sequences = [
            "ACGTACGTACGT",
            "ACGTACGTACGT",
            "ACGTTCGTACGT",
        ]
        records = FASTX.FASTQ.Record[
            FASTX.FASTQ.Record("forward_$(index)", sequence, repeat("I", length(sequence)))
            for (index, sequence) in enumerate(sequences)
        ]
        k = 5
        graph = R.build_qualmer_graph(records, k; mode = :canonical)
        probe = last(records)
        probe_sequence = FASTX.sequence(BioSequences.LongDNA{4}, probe)
        hard_vertices = Set(
            BioSequences.canonical(kmer) for
            (kmer, _) in Kmers.UnambiguousDNAMers{k}(probe_sequence)
        )

        direct = Mycelia.try_viterbi_path_improvement(
            probe,
            graph,
            k;
            graph_mode = :canonical,
            substitution_error_rate = 0.001,
        )
        Test.@test direct isa Union{Nothing, Tuple{FASTX.FASTQ.Record, Float64}}
        direct_diagnostics = Mycelia.CorrectorDiagnostics()
        Test.@test Mycelia.try_viterbi_path_improvement(
            probe,
            graph,
            k;
            graph_mode = :canonical,
            diagnostics = direct_diagnostics,
            substitution_error_rate = 0.0,
        ) === nothing
        Test.@test direct_diagnostics.structural_errors[] > 0

        windowed, _improved, decoded_windows,
        divergent_windows = Mycelia.improve_read_likelihood_windowed_detail(
            probe,
            graph,
            k,
            hard_vertices;
            graph_mode = :canonical,
            substitution_error_rate = 0.001,
        )
        Test.@test decoded_windows >= 1
        Test.@test divergent_windows == 0
        Test.@test length(FASTX.sequence(windowed)) == length(FASTX.sequence(probe))
        window_diagnostics = Mycelia.CorrectorDiagnostics()
        Mycelia.improve_read_likelihood_windowed_detail(
            probe,
            graph,
            k,
            hard_vertices;
            graph_mode = :canonical,
            diagnostics = window_diagnostics,
            substitution_error_rate = 0.0,
        )
        Test.@test window_diagnostics.structural_errors[] > 0

        batch, _improvements, _skip_fraction, _cheap_corrections,
        _decode_gated = Mycelia.improve_read_set_likelihood(
            records,
            graph,
            k;
            graph_mode = :canonical,
            hard_vertices,
            windowed_decode = true,
            substitution_error_rate = 0.001,
        )
        Test.@test length(batch) == length(records)
        batch_diagnostics = Mycelia.CorrectorDiagnostics()
        Mycelia.improve_read_set_likelihood(
            records,
            graph,
            k;
            graph_mode = :canonical,
            hard_vertices,
            windowed_decode = true,
            diagnostics = batch_diagnostics,
            substitution_error_rate = 0.0,
        )
        Test.@test batch_diagnostics.structural_errors[] > 0
    end

    Test.@testset "AssemblyConfig threads + validates sequencing_tech" begin
        # Default is :illumina (unchanged behavior).
        Test.@test R.AssemblyConfig(k = 13).sequencing_tech == :illumina
        Test.@test R.AssemblyConfig(k = 13, sequencing_tech = :nanopore).sequencing_tech ==
                   :nanopore
        Test.@test R.AssemblyConfig(
            k = 13,
            sequencing_tech = :pacbio_hifi,
        ).sequencing_tech == :pacbio_hifi
        test_throws_message(ErrorException, "sequencing_tech must be one of") do
            R.AssemblyConfig(k = 13, sequencing_tech = :bogus)
        end
    end

    Test.@testset "tier knobs carry indel bounds" begin
        sc = R._corrector_strategy_knobs(:scalable)
        # :scalable => bounded band + small run caps.
        Test.@test sc.band_width isa Int
        Test.@test sc.band_width > 0
        Test.@test sc.deletion_max_run == 3
        Test.@test sc.max_insertion_run >= 1
        Test.@test sc.indel_schedule == :frontier_budgeted

        ex = R._corrector_strategy_knobs(:exhaustive)
        # :exhaustive => unbounded band + larger caps and no runtime classifier.
        Test.@test ex.band_width === nothing
        Test.@test ex.deletion_max_run > sc.deletion_max_run
        Test.@test ex.indel_schedule == :unrestricted
    end

    Test.@testset "illumina + ultima preserve substitution oracle bytes" begin
        reads, _ = _profile_reads(tech = :illumina, err = 0.01, seed = 7)
        base = R.assemble_genome(reads; k = 13, corrector = :iterative)
        il = R.assemble_genome(reads; k = 13, corrector = :iterative,
            sequencing_tech = :illumina)
        ul = R.assemble_genome(reads; k = 13, corrector = :iterative,
            sequencing_tech = :ultima)
        hifi = R.assemble_genome(reads; k = 13, corrector = :iterative,
            sequencing_tech = :pacbio_hifi)
        # Oracle guard (PR #408 review, FIX 4): the byte-identity check below is
        # only meaningful if the corrector actually produced contigs — otherwise
        # `[] == []` passes vacuously and hides a corrector that ate every read.
        Test.@test !isempty(base.contigs)
        # Contigs byte-identical: the illumina profile threads NO indel params, so
        # the corrector runs the substitution decode unchanged.
        Test.@test il.contigs == base.contigs
        Test.@test ul.contigs == base.contigs
        Test.@test base.assembly_stats["indel_moves"] == false
        Test.@test il.assembly_stats["indel_moves"] == false
        Test.@test ul.assembly_stats["indel_moves"] == false
        Test.@test hifi.assembly_stats["indel_moves"] == false
        Test.@test il.assembly_stats["sequencing_tech"] == "illumina"
        Test.@test ul.assembly_stats["sequencing_tech"] == "ultima"
        Test.@test hifi.assembly_stats["sequencing_tech"] == "pacbio_hifi"
        Test.@test hifi.assembly_stats["substitution_error_rate"] == 0.001
        Test.@test il.assembly_stats["substitution_error_rate"] === nothing

        il_telemetry = il.assembly_stats["indel_rung_telemetry"]
        Test.@test il_telemetry isa AbstractVector
        Test.@test all(row -> !get(row, :profile_requested, false), il_telemetry)
        Test.@test all(row -> get(row, :requested, 0) == 0, il_telemetry)
        Test.@test il.assembly_stats["indel_requested"] == 0
        Test.@test il.assembly_stats["indel_attempted"] == 0
        Test.@test il.assembly_stats["indel_completed"] == 0
        Test.@test il.assembly_stats["indel_truncated"] == 0
        Test.@test il.assembly_stats["indel_engaged"] == 0
        Test.@test il.assembly_stats["indel_decodes"] ==
                   il.assembly_stats["indel_completed"]
        Test.@test il.assembly_stats["truncated_decodes"] ==
                   il.assembly_stats["indel_truncated"]
    end

    Test.@testset "nanopore enables indel-aware correction" begin
        reads, _ = _profile_reads(tech = :nanopore, err = 0.05, seed = 7)
        np = R.assemble_genome(reads; k = 13, corrector = :iterative,
            sequencing_tech = :nanopore)
        Test.@test np isa R.AssemblyResult
        Test.@test !isempty(np.contigs)
        # The legacy provenance stamp records profile intent. Runtime engagement
        # is surfaced separately by the aggregate and per-rung telemetry.
        Test.@test np.assembly_stats["indel_moves"] == true
        Test.@test np.assembly_stats["sequencing_tech"] == "nanopore"
        np_telemetry = np.assembly_stats["indel_rung_telemetry"]
        Test.@test np_telemetry isa AbstractVector
        Test.@test np.assembly_stats["indel_requested"] >= 0
        Test.@test np.assembly_stats["indel_attempted"] >= 0
        Test.@test np.assembly_stats["indel_completed"] >= 0
        Test.@test np.assembly_stats["indel_truncated"] >= 0
        Test.@test np.assembly_stats["indel_engaged"] >= 0
        for field in (:requested, :attempted, :completed, :truncated, :engaged)
            key = "indel_$(field)"
            Test.@test np.assembly_stats[key] ==
                       sum(get(row, field, 0) for row in np_telemetry)
        end
        Test.@test np.assembly_stats["indel_requested"] >=
                   np.assembly_stats["indel_attempted"]
        Test.@test np.assembly_stats["indel_attempted"] ==
                   np.assembly_stats["indel_completed"] +
                   np.assembly_stats["indel_truncated"]
        Test.@test np.assembly_stats["indel_completed"] >=
                   np.assembly_stats["indel_engaged"]
        Test.@test np.assembly_stats["indel_decodes"] ==
                   np.assembly_stats["indel_completed"]
        Test.@test np.assembly_stats["truncated_decodes"] ==
                   np.assembly_stats["indel_truncated"]
        Test.@test np.assembly_stats["trace_contract_errors"] == 0
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
        np_indel_requested = Int[]
        np_indel_decodes = Int[]
        np_indel_engaged = Int[]
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
            push!(np_indel_requested, np.assembly_stats["indel_requested"])
            push!(np_indel_decodes, np.assembly_stats["indel_decodes"])
            push!(np_indel_engaged, np.assembly_stats["indel_engaged"])
            push!(np_nonempty, !isempty(np.contigs))
        end
        mean_np = Statistics.mean(id_np)
        mean_il = Statistics.mean(id_il)
        wins = count(id_np .>= id_il)
        @info "E2E nanopore-vs-illumina identity" seeds = collect(seeds) id_np id_il mean_np mean_il wins

        # Every nanopore arm requests indel service. The branching/frontier
        # classifier may correctly reject individual high-branching seeds, so the
        # matrix requires actual completed and gap-engaged decodes across the
        # fixed seeds rather than forcing every seed through an unaffordable rung.
        Test.@test all(np_indel_moves)
        Test.@test all(>(0), np_indel_requested)
        Test.@test any(>(0), np_indel_decodes)
        Test.@test any(>(0), np_indel_engaged)
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
