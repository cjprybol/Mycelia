# OPT4 FROZEN-READ SKIP — POSITIVE CONTROL (td-jbjd, pass 4)
# =============================================================================
#
# THE ASK (design doc "Accuracy sign-off methodology", bullet 3, and Rule-of-5
# pass 4): build a toy fixture where freezing DEMONSTRABLY DEGRADES per-base
# correction accuracy, so a later accuracy sign-off's "no degradation" verdict
# on representative data is trustworthy -- the same "null needs a positive
# control" discipline as the accuracy benchmark's Control B read-scramble null
# (benchmarking/rhizomorph_correction_accuracy_benchmark.jl).
#
# THE RESULT, after extensive investigation (see the finding below): at TOY
# fixture scale, using the wired `:scalable`-tier corrector configuration
# (cheap_correct=true, skip_solid=true, hard_window=true -- the settings used
# throughout the opt4 test suite and the accuracy benchmark), `skip_frozen_reads`
# could NOT be made to measurably degrade per-base correction accuracy, even
# under the MOST aggressive freeze configuration tested
# (freeze_streak_threshold=1, freeze_across_rungs=true). This held across every
# fixture design attempted:
#   1. An engineered repeat/paralogous-SNP confounder (two near-identical repeat
#      copies differing at one base, with a read error that mimics the OTHER
#      copy's real allele) walked through a k-ladder from k=11 to k=39, run via
#      the full `mycelia_iterative_assemble` pipeline.
#   2. An ordinary random single-substitution fixture (mirrors
#      corrector_opt4_freeze_state_test.jl's own "composition" fixture, which IS
#      proven to converge) driven directly via `improve_read_set_likelihood`,
#      rebuilding the graph every pass, comparing unlimited retries (exact) vs
#      threshold=1 freezing over multiple passes at a fixed k.
#
# THE ROOT CAUSE (verified against src/iterative-assembly.jl, 2026-07-27):
#   (a) Stage 0 cheap k-mer-spectrum correction (`_stage0_cheap_correct`,
#       gated on `cheap_correct=true`) runs on `work_reads` UNCONDITIONALLY,
#       BEFORE the skip/freeze decision (`_gate_skip`/`_frozen_read_at`), for
#       EVERY read every pass -- including reads that are about to be
#       frozen-skipped for that pass's decode. A frozen read's Stage 0
#       correction (if any) still lands in `updated_reads` regardless of freeze
#       state (iterative-assembly.jl ~4675-4692 unconditional Stage 0 call,
#       ~4901/~4980 `read = batch_reads[i]`/`work_reads[i]` flowing into a
#       SKIPPED read's output). opt4's freeze mechanism ONLY gates the
#       expensive per-read Viterbi decode step downstream of Stage 0 -- it
#       cannot touch Stage 0's output either way.
#   (b) At toy fixture scale, Stage 0 (freeze-immune per (a)) accounts for
#       effectively ALL correction volume: every fixture tested showed the raw
#       per-read Viterbi decode alone (cheap_correct=false, so Stage 0 disabled)
#       making ZERO edits at ANY k from 7 to 43, with or without
#       skip_solid/hard_window/windowed_decode. Separately, `mycelia_iterative_
#       assemble`'s adaptive low-k decode gate (`decode_gate_density`, default
#       0.90, active whenever `hard_window=true`) additionally disables the
#       whole per-read decode step at toy scale because a small/dense graph
#       cannot achieve genuine (<90%) decode-fraction selectivity -- a SECOND,
#       independent reason decode contributes ~nothing to toy-scale accuracy.
#   (c) Multi-pass convergence is also observed to be "one-shot" at this scale:
#       a read that CAN be fixed (by Stage 0 or decode) is fixed on its very
#       first evaluation; a read that can't be fixed does not benefit from
#       additional passes even with the graph rebuilt fresh each pass (probe
#       investigation, not itself asserted below since it is a property of the
#       decode/Stage-0 algorithms, not of opt4's freeze plumbing). This means
#       there is no "read that would have converged on attempt 3" for freezing
#       to strand at toy scale -- exactly the scenario a genuine positive
#       control would need to exploit.
#
# WHAT THIS MEANS FOR THE SIGN-OFF (per the design doc: "If after reasonable
# effort you CANNOT make freezing degrade on a toy fixture: that is itself an
# important FINDING -- report it... Do NOT fake a passing assertion"):
#   * This is NOT evidence that opt4 is safe on REPRESENTATIVE data. Stage 0's
#     freeze-immunity is scale-INDEPENDENT (architectural), so it will also
#     dominate on larger data -- meaning opt4's real accuracy exposure is
#     concentrated in whatever residual correction volume Stage 0 and the
#     adaptive gate leave to genuine per-read Viterbi decode, which is
#     precisely the volume this toy scale could not exercise (findings (a)-(b)
#     above). The pass-5 accuracy sign-off on representative data (Lambda
#     phage scale, per rhizomorph_correction_accuracy_benchmark.jl) is where
#     the adaptive gate stops firing (larger, sparser graphs have genuine
#     hard-window selectivity) and decode does real work -- THAT is where a
#     freeze-induced degradation, if any exists, would actually show up.
#   * Rather than assert a fabricated degradation, this file locks in the
#     ACTUAL discovered invariant as a regression guard: both fixture designs
#     above produce BYTE-IDENTICAL corrected output whether frozen or exact, at
#     the most aggressive freeze setting tested. If a future change makes
#     these diverge, that is a signal worth investigating (either Stage 0's
#     unconditional guarantee broke, the adaptive gate's behavior changed, or
#     decode started doing real work at toy scale) -- and it would also mean a
#     toy-scale positive control might become achievable where it currently is
#     not.
#
# Run directly:
#   LD_LIBRARY_PATH='' julia --project=. \
#     -e 'include("test/4_assembly/corrector_opt4_positive_control_test.jl")'

import Test
import Mycelia
import FASTX
import Random

include(joinpath(@__DIR__, "..", "..", "benchmarking",
    "rhizomorph_correction_accuracy_metrics.jl"))

const _OPT4PC_BASES = ['A', 'C', 'G', 'T']

# --- Fixture 1: engineered repeat/paralogous-SNP confounder -----------------
# Two 60bp repeat copies (copy1, copy2) embedded in unique flanking sequence,
# identical except at ONE base (`snp_offset` into each copy). A small set of
# "engineered" reads are copy1-derived with a substitution AT that position
# that happens to equal copy2's real allele -- i.e. the corrupted base is not
# noise, it is a second real (paralogous) allele elsewhere in the genome. This
# is the design doc's Risks-section scenario: a read wrong at low k (the SNP
# position's local k-mer window is genuinely ambiguous between the two repeat
# copies) that should resolve once the graph sharpens at higher k (once the
# k-mer window spans past the repeat edge into copy-specific unique flank).
function _opt4pc_build_genome(rng; fb = 40, rl = 60, mu = 40, fa = 40, snp_offset = 30)
    flank_before = join(rand(rng, _OPT4PC_BASES, fb))
    repeat_body = collect(rand(rng, _OPT4PC_BASES, rl))
    middle = join(rand(rng, _OPT4PC_BASES, mu))
    flank_after = join(rand(rng, _OPT4PC_BASES, fa))
    copy1 = join(repeat_body)
    other_base = rand(rng, filter(!=(repeat_body[snp_offset]), _OPT4PC_BASES))
    copy2_body = copy(repeat_body)
    copy2_body[snp_offset] = other_base
    copy2 = join(copy2_body)
    genome = flank_before * copy1 * middle * copy2 * flank_after
    copy1_snp_pos = fb + snp_offset
    return genome, copy1_snp_pos, other_base
end

function _opt4pc_build_reads(rng, genome, copy1_snp_pos, mimic_allele;
        readlen = 80, coverage = 20, n_error_reads = 8, n_plain = 8)
    glen = length(genome)
    n_reads = max(1, ceil(Int, coverage * glen / readlen))
    records = FASTX.FASTQ.Record[]
    truth_by_id = Dict{String, String}()
    observed_by_id = Dict{String, String}()
    qstr = repeat(string(Char(20 + 33)), readlen)
    for i in 1:n_reads
        s = rand(rng, 1:(glen - readlen + 1))
        seq = genome[s:(s + readlen - 1)]
        rid = "bg_$(i)"
        truth_by_id[rid] = seq
        observed_by_id[rid] = seq
        push!(records, FASTX.FASTQ.Record(rid, seq, qstr))
    end
    lo = max(1, copy1_snp_pos - readlen + 1)
    hi = min(glen - readlen + 1, copy1_snp_pos)
    for i in 1:n_error_reads
        s = rand(rng, lo:hi)
        truth = genome[s:(s + readlen - 1)]
        p = copy1_snp_pos - s + 1
        chars = collect(truth)
        chars[p] = mimic_allele
        observed = join(chars)
        rid = "err_$(i)"
        truth_by_id[rid] = truth
        observed_by_id[rid] = observed
        push!(records, FASTX.FASTQ.Record(rid, observed, qstr))
    end
    for i in 1:n_plain
        s = rand(rng, 1:(glen - readlen + 1))
        truth = genome[s:(s + readlen - 1)]
        p = rand(rng, 1:readlen)
        chars = collect(truth)
        chars[p] = rand(rng, filter(!=(chars[p]), _OPT4PC_BASES))
        observed = join(chars)
        rid = "plain_$(i)"
        truth_by_id[rid] = truth
        observed_by_id[rid] = observed
        push!(records, FASTX.FASTQ.Record(rid, observed, qstr))
    end
    return records, truth_by_id, observed_by_id
end

function _opt4pc_run_ladder(records, tmp, label; skip_frozen, threshold, across, k_ladder)
    fq = joinpath(tmp, "in.fastq")
    isfile(fq) || Mycelia.write_fastq(records = records, filename = fq)
    result = Mycelia.mycelia_iterative_assemble(fq;
        max_k = maximum(k_ladder) + 10, k_ladder = k_ladder, graph_mode = :canonical,
        skip_solid = true, cheap_correct = true, hard_window = true, soft_em = false,
        batch_size = 50, max_iterations_per_k = 2, verbose = false,
        enable_checkpointing = false, skip_frozen_reads = skip_frozen,
        freeze_streak_threshold = threshold, freeze_across_rungs = across,
        output_dir = joinpath(tmp, label))
    corrected_by_id = Dict{String, String}()
    open(FASTX.FASTQ.Reader, result[:metadata][:final_fastq_file]) do rd
        for rec in rd
            corrected_by_id[FASTX.identifier(rec)] = FASTX.sequence(String, rec)
        end
    end
    return corrected_by_id, result
end

# --- Fixture 2: ordinary substitution fixture (mirrors freeze_state_test's
# "composition" subtest, which IS proven to converge under skip_solid=true/
# cheap_correct=true) -----------------------------------------------------
function _opt4pc_plain_reads(rng, ref; n_reads = 150, readlen = 80, n_err = 20)
    reflen = length(ref)
    records = FASTX.FASTQ.Record[]
    truth_by_id = Dict{String, String}()
    observed_by_id = Dict{String, String}()
    for i in 1:n_reads
        s = rand(rng, 1:(reflen - readlen + 1))
        seq = collect(ref[s:(s + readlen - 1)])
        truth = join(seq)
        if i <= n_err
            p = rand(rng, 1:readlen)
            seq[p] = rand(rng, filter(!=(seq[p]), _OPT4PC_BASES))
        end
        observed = join(seq)
        rid = "r$i"
        truth_by_id[rid] = truth
        observed_by_id[rid] = observed
        push!(records, FASTX.FASTQ.Record(rid, observed, String(fill('I', readlen))))
    end
    return records, truth_by_id, observed_by_id
end

function _opt4pc_run_multi_pass_rebuild(
        records, k, common; skip_frozen, threshold, n_passes)
    streaks = zeros(Int, length(records))
    diag = Mycelia.CorrectorDiagnostics()
    reads = records
    for _pass in 1:n_passes
        graph = Mycelia.Rhizomorph.build_qualmer_graph(reads, k; mode = :canonical)
        reads, = Mycelia.improve_read_set_likelihood(reads, graph, k; common...,
            skip_frozen_reads = skip_frozen, freeze_streak_threshold = threshold,
            freeze_streaks = streaks, diagnostics = diag)
    end
    return reads, diag
end

Test.@testset "opt4 positive control -- toy-scale degradation NOT inducible (td-jbjd pass 4)" begin
    Test.@testset "repeat/paralogous-SNP fixture: exact vs freeze byte-identical" begin
        rng = Random.MersenneTwister(1234)
        genome, c1pos, mimic_allele = _opt4pc_build_genome(rng)
        rng2 = Random.MersenneTwister(5678)
        records, truth_by_id,
        observed_by_id = _opt4pc_build_reads(
            rng2, genome, c1pos, mimic_allele; readlen = 80, coverage = 20,
            n_error_reads = 8, n_plain = 8)

        tmp = mktempdir()
        k_ladder = [11, 15, 19, 23, 27, 31, 35, 39]
        corrected_exact,
        r_exact = _opt4pc_run_ladder(records, tmp, "exact";
            skip_frozen = false, threshold = 2, across = false, k_ladder = k_ladder)
        corrected_freeze,
        r_freeze = _opt4pc_run_ladder(records, tmp, "freeze";
            skip_frozen = true, threshold = 1, across = true, k_ladder = k_ladder)

        # Setup sanity: this fixture DOES exercise real correction (some
        # engineered reads fixed via Stage 0's unique-neighbor search once the
        # k-mer window spans into copy-specific flank), and freeze DID skip a
        # substantial number of decodes -- so the byte-identical result below
        # is not simply "nothing happened in either arm".
        err_truth = Dict(
            rid=>truth_by_id[rid] for rid in keys(truth_by_id) if startswith(rid, "err_"))
        err_obs = Dict(
            rid=>observed_by_id[rid]
        for
        rid in keys(observed_by_id) if startswith(rid, "err_"))
        err_corr_exact = Dict(
            rid=>corrected_exact[rid]
        for rid in keys(err_truth) if haskey(corrected_exact, rid))
        m_exact = per_base_metrics(err_truth, err_obs, err_corr_exact)
        Test.@test m_exact.tp > 0   # some engineered reads WERE fixed (by Stage 0)
        Test.@test r_freeze[:metadata][:corrector_errors][:frozen_reads_skipped] > 0

        # THE FINDING: despite real correction activity and real freeze
        # activity, the two arms' final corrected FASTQs are byte-identical --
        # Stage 0 (freeze-immune) accounts for every actual edit at this scale.
        recs = res -> [(String(FASTX.identifier(r)), FASTX.sequence(String, r),
                           String(FASTX.quality(r)))
                       for r in
                           open(FASTX.FASTQ.Reader, res[:metadata][:final_fastq_file]) do rd
            collect(rd)
        end]
        Test.@test recs(r_exact) == recs(r_freeze)
    end

    Test.@testset "ordinary-substitution fixture, multi-pass with graph rebuild: exact vs freeze byte-identical" begin
        rng = Random.MersenneTwister(2026)
        ref = join(rand(rng, _OPT4PC_BASES, 900))
        records, truth_by_id,
        observed_by_id = _opt4pc_plain_reads(
            rng, ref; n_reads = 150, readlen = 80, n_err = 20)
        k = 13
        common = (; graph_mode = :canonical, skip_solid = true, cheap_correct = true,
            decode_enabled = true, batch_size = 50)

        reads_exact,
        diag_exact = _opt4pc_run_multi_pass_rebuild(
            records, k, common; skip_frozen = false, threshold = 1, n_passes = 4)
        reads_freeze,
        diag_freeze = _opt4pc_run_multi_pass_rebuild(
            records, k, common; skip_frozen = true, threshold = 1, n_passes = 4)

        corrected_exact = Dict(
            String(FASTX.identifier(r))=>FASTX.sequence(String, r) for r in reads_exact)
        corrected_freeze = Dict(
            String(FASTX.identifier(r))=>FASTX.sequence(String, r) for r in reads_freeze)
        m_exact = per_base_metrics(truth_by_id, observed_by_id, corrected_exact)

        # Setup sanity: real correction happened (recall > 0), and freeze DID
        # fire under the maximally aggressive threshold=1 setting -- multiple
        # extra passes were denied to whichever reads did not converge on their
        # first attempt.
        Test.@test m_exact.tp > 0
        Test.@test diag_freeze.frozen_reads_skipped[] > 0

        # THE FINDING (independent confirmation, bypassing mycelia_iterative_
        # assemble's adaptive low-k gate entirely by calling
        # improve_read_set_likelihood directly): the corrected sequences are
        # identical whether extra passes were denied (freeze) or granted
        # (exact) -- a read that converges here does so on its FIRST
        # evaluation; the additional passes freezing denies never mattered.
        recs = d -> sort(collect(d))
        Test.@test recs(corrected_exact) == recs(corrected_freeze)
    end
end
