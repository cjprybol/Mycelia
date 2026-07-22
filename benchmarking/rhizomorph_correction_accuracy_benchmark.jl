# Rhizomorph Correction ACCURACY Benchmark — precision / over-correction / recall
# ==============================================================================
#
# Measures PER-BASE correction quality of the WIRED :scalable corrector against
# injected-error ground truth. This closes the gap left by the two existing
# harnesses:
#
#   * benchmarking/viterbi_accuracy_benchmark.jl ("B8") measures recall and net
#     edit-distance reduction — but NOT per-base precision / over-correction rate
#     — and targets the bare `correct_observations` decode primitive, NOT the
#     shipped assemble_genome(corrector=:iterative, strategy=:scalable) pipeline.
#   * benchmarking/rhizomorph_correction_validation_sweep.jl measures ASSEMBLY
#     metrics (N50 / genome_fraction), not per-base correction accuracy.
#
# The central question (user direction 2026-07-11): does the wired corrector
# fix injected errors WITHOUT over-correcting (changing a correct base to a
# wrong one), and is its edit VOLUME ~ the injected error rate?
#
# TWO CONTROLS (the pair the manuscript flags as missing next to a recall table):
#   * Control A — over-correction on UN-CORRUPTED input. The default error-rate
#     sweep includes an err=0.0 cell: no substitutions are injected, so every edit
#     the corrector makes on that row is a false positive. recall/precision are
#     NaN there (injected==0, an "undefined here" signal, not a measured 0); the
#     load-bearing columns are over_correction_rate and correction_rate, which
#     directly measure whether the WIRED corrector damages already-correct reads.
#   * Control B — read-scramble null (the discriminating null). For each read a
#     fixed-seed PER-READ permutation is applied to its bases BEFORE correction
#     (and the SAME permutation to that read's truth for scoring). Distinct
#     per-read permutations destroy the shared k-mer coverage the corrector relies
#     on, so the internally-built graph has no consensus and CANNOT recover the
#     (still-present) injected errors: null recall should COLLAPSE relative to the
#     real arm. Null columns (null_recall, null_over_correction_rate, ...) sit
#     beside the real metrics in the same CSV row. The injected-error COUNT is
#     preserved by the permutation (a bijection applied identically to observed
#     and truth), so real and null recall are directly comparable.
#
# WHY per-base and not assembly metrics: a prior program error (td-wlfq
# retraction) quoted an assembly-level genome-fraction "improvement" that was
# actually a read-coverage artifact. Correction quality is a per-read / per-base
# property and must be measured there, against known injected errors.
#
# CORRECTOR ENTRY POINT: this harness calls `Mycelia.mycelia_iterative_assemble`
# directly with the EXACT knobs the `:scalable` branch of
# `_corrector_strategy_knobs` supplies to the wired path (both live in
# src/rhizomorph/assembly.jl; symbol-anchored, not line-anchored, since that file
# is actively edited). `mycelia_iterative_assemble` IS the read-corrector half of
# assemble_genome(corrector=:iterative, strategy=:scalable) — assemble_genome
# merely re-assembles its corrected reads into contigs (which would discard the
# per-read identity this harness needs). Corrected reads are read back from
# `result[:metadata][:final_fastq_file]` and joined to their truth fragments by
# read ID (unkmerizable/skipped reads pass through uncorrected; the ID-join
# surfaces any drop-outs explicitly rather than silently miscounting).
#
# GROUND TRUTH (substitution-only, this version): errors are injected with
# `Mycelia.mutate_dna_substitution_fraction` (length-preserving, no indels), so
# truth/observed/corrected are positionally comparable and per-base metrics are
# EXACT with no alignment. Injected positions are recovered by diffing observed
# vs truth (E = {p : observed[p] != truth[p]}). Indel-bearing regimes
# (nanopore/pacbio) need an alignment-based classifier and are a tracked
# follow-on; a length-guard here excludes (and reports) any read whose corrected
# length differs from truth so a Viterbi indel path can never corrupt the metric.
#
# READ QUALITY: reads are emitted with a UNIFORM Phred (MYCELIA_RCA_ASSIGNED_Q,
# default 20). This deliberately does NOT leak truth into the per-base quality
# string — the corrector must detect errors from GRAPH COVERAGE (a substitution
# spawns low-coverage k-mers), which is the mechanism under test. A truth-derived
# quality string would make the benchmark trivially easy and defensibility-void.
#
# METRICS (aggregated over all scored reads in a cell):
#   recall          = true_fixes / injected_errors
#   precision       = true_fixes / total_edits          (total_edits = C != O)
#   over_correction = over_corrections / correct_positions  (O==T but C!=T)
#   correction_rate = total_edits / total_bases         (compare to error_rate)
# where a true_fix is an injected position corrected to truth (C==T at O!=T). An
# edit is one of three classes — true_fix, mis_fix (error edited to a WRONG base),
# or over_correction (correct base made wrong) — and they partition exactly:
# total_edits == true_fixes + mis_fixes + over_corrections. mis_fixes vs
# over_corrections attributes precision loss to failed corrections vs collateral
# damage. Zero-denominator metrics are NaN ("undefined here", not a measured 0).
#
# GUARDS: a VERDICT requires BOTH the scale floor AND a minimum aggregate SCORED
# fraction (MYCELIA_RCA_MIN_SCORED_FRACTION) so clean-looking metrics on a tiny
# surviving subset (mass drops/exclusions) cannot be quoted as validation; a
# crashed cell is marked ok=false in the CSV; an empty corrector output errors.
#
# SCALE GUARD: reuses rhizomorph_scale_guard.jl (via the sweep include). Below
# the floor the run is labelled SMOKE-ONLY and MUST NOT be quoted as validation.
#
# Usage:
#   # Local smoke (synthetic ~2 kb genome, no network, SMOKE-ONLY):
#   MYCELIA_RCA_SMOKE=true LD_LIBRARY_PATH='' \
#     julia --project=. benchmarking/rhizomorph_correction_accuracy_benchmark.jl
#
#   # Full run (Lambda phage, real coverage, VERDICT) on Lovelace:
#   MYCELIA_RCA_COVERAGE=30 LD_LIBRARY_PATH='' \
#     julia --project=. benchmarking/rhizomorph_correction_accuracy_benchmark.jl
#
# Environment variables (mirror the validation-sweep interface, RCA prefix):
#   MYCELIA_RCA_SMOKE            truthy -> synthetic genome, no download (default false)
#   MYCELIA_RCA_ACCESSION       reference accession (default NC_001416, Lambda phage)
#   MYCELIA_RCA_SMOKE_GENOME_LEN synthetic genome length in smoke mode (default 2000)
#   MYCELIA_RCA_ERR             comma-separated substitution rates (default 0.0,0.01,0.05,0.10; the 0.0 cell is Control A)
#   MYCELIA_RCA_READLEN         comma-separated read lengths (default 150)
#   MYCELIA_RCA_COVERAGE        target fold coverage (default 30; smoke default 10)
#   MYCELIA_RCA_K               assembly k-mer size (default 21)
#   MYCELIA_RCA_BATCH_SIZE      corrector read-batch size (default 10000; set small, e.g. 50, to expose between-batch GC/parallel effects)
#   MYCELIA_RCA_ASSIGNED_Q      uniform Phred assigned to every base (default 20)
#   MYCELIA_RCA_SEED            RNG seed (default 42)
#   MYCELIA_RCA_SCALE_FLOOR     override the scale-guard floor in bases

import Pkg
if isinteractive()
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# Reuse acquire_reference, regime_for_readlen, the _parse_* helpers, and the
# scale guard (rhizomorph_scale_guard.jl, which the sweep itself includes) rather
# than duplicating them. The sweep's run_sweep() is guarded by
# `PROGRAM_FILE == @__FILE__`, so including it here only DEFINES functions.
include(joinpath(@__DIR__, "rhizomorph_correction_validation_sweep.jl"))

# === Substitution-only read simulation with retained truth ==================

"""
Simulate fixed-length substitution-only reads for one (error_rate, readlen) cell.
Each read samples a random fragment (either strand) of `refseq`, injects
`ceil(error_rate * readlen)` substitutions (so the realized rate is >= nominal;
e.g. err=0.01 at readlen=150 -> ceil(1.5)=2 -> 0.0133 realized) via
`Mycelia.mutate_dna_substitution_fraction` (length-preserving), and is emitted
with a uniform Phred `assigned_q`.

Returns `(records, truth_by_id, observed_by_id, injected_total, sampled_bases)`:
`truth_by_id[id]` / `observed_by_id[id]` are the pre-/post-error read sequences
(equal length), keyed by the FASTQ read id so corrected reads can be joined back.
"""
function simulate_substitution_reads(refseq::BioSequences.LongDNA{4}, readlen::Int,
        coverage::Real, error_rate::Float64, rng::Random.AbstractRNG; assigned_q::Int = 20)
    glen = length(refseq)
    effective_readlen = min(readlen, glen)
    n_reads = max(1, ceil(Int, coverage * glen / effective_readlen))
    records = Vector{FASTX.FASTQ.Record}()
    truth_by_id = Dict{String, String}()
    observed_by_id = Dict{String, String}()
    injected_total = 0
    sampled_bases = 0
    qstr_for = (n) -> repeat(string(Char(assigned_q + 33)), n)
    for i in 1:n_reads
        start = glen == effective_readlen ? 1 : rand(rng, 1:(glen - effective_readlen + 1))
        frag = refseq[start:(start + effective_readlen - 1)]
        rand(rng, Bool) && (frag = BioSequences.reverse_complement(frag))
        truth = String(frag)
        observed_seq = Mycelia.mutate_dna_substitution_fraction(frag; fraction = error_rate, rng = rng)
        observed = String(observed_seq)
        rid = "read_$(i)"
        # Substitution-only + length-preserving => positional diff is the exact
        # injected-error set for this read.
        injected_total += count(p -> observed[p] != truth[p], eachindex(truth))
        sampled_bases += effective_readlen
        truth_by_id[rid] = truth
        observed_by_id[rid] = observed
        push!(records, FASTX.FASTQ.Record(rid, observed, qstr_for(length(observed))))
    end
    return records, truth_by_id, observed_by_id, injected_total, sampled_bases
end

# === Run the wired :scalable corrector on one read set ======================

"""
Correct `records` with the WIRED :scalable knobs (exact replica of the
`:scalable` branch of `_corrector_strategy_knobs` as consumed by the
`mycelia_iterative_assemble` call in `_assemble_with_iterative_corrector`, both
in src/rhizomorph/assembly.jl). Returns `(corrected_by_id, result_dict)` where
`corrected_by_id[id]` is the corrected read sequence keyed by read id.
"""
function correct_reads_scalable(records::Vector{FASTX.FASTQ.Record}, k::Int)
    input_dir = mktempdir()
    output_dir = mktempdir()
    temp_fastq = joinpath(input_dir, "corrector_input.fastq")
    open(FASTX.FASTQ.Writer, temp_fastq) do w
        for r in records
            write(w, r)
        end
    end
    max_k = max(k, 13)
    batch_size = parse(Int, get(ENV, "MYCELIA_RCA_BATCH_SIZE", "10000"))
    batch_size > 0 ||
        error("MYCELIA_RCA_BATCH_SIZE must be a positive integer; got $batch_size")
    result = Mycelia.mycelia_iterative_assemble(
        temp_fastq;
        max_k = max_k,
        skip_solid = true,
        graph_mode = :doublestrand,
        n_k_rungs = 3,
        max_iterations_per_k = 2,
        batch_size = batch_size,
        hard_window = true,
        soft_em = true,
        cheap_correct = true,
        beam_width = nothing,
        verbose = false,
        enable_checkpointing = false,
        output_dir = output_dir
    )
    corrected_fastq = get(result[:metadata], :final_fastq_file, nothing)
    if corrected_fastq === nothing
        error("iterative corrector metadata is missing :final_fastq_file; " *
              "keys=$(collect(keys(result[:metadata])))")
    elseif !isfile(corrected_fastq)
        error("iterative corrector :final_fastq_file points at a nonexistent path: $corrected_fastq")
    end
    corrected_by_id = Dict{String, String}()
    open(FASTX.FASTQ.Reader, corrected_fastq) do reader
        for rec in reader
            corrected_by_id[FASTX.identifier(rec)] = String(FASTX.sequence(BioSequences.LongDNA{4}, rec))
        end
    end
    # A structurally-valid but EMPTY corrected FASTQ (zero records) is never a
    # legitimate result for a non-empty input — it would otherwise flow to
    # per_base_metrics as "every read dropped" (all-NaN metrics) with no error.
    # Fail loud so an empty-output corrector bug can't masquerade as a real run.
    if !isempty(records) && isempty(corrected_by_id)
        error("iterative corrector produced 0 corrected reads from $(length(records)) input reads " *
              "(corrected FASTQ $(corrected_fastq) was empty); refusing to score an empty result")
    end
    return corrected_by_id, result
end

# === Control B: read-scramble null ==========================================
# Fixed seed for the per-read scramble null so the null CSV columns are
# byte-reproducible across runs. Each read's permutation is drawn from
# MersenneTwister(RCA_NULL_SEED + read_index), so permutations DIFFER across reads
# (which is what destroys the shared k-mer coverage), yet the whole null arm is
# deterministic.
const RCA_NULL_SEED = 20260711

"""
    scramble_reads(records, truth_by_id, observed_by_id, assigned_q; seed=RCA_NULL_SEED)

Build the Control-B null read set. For each read draw a PER-READ base permutation
(seeded from `seed + read_index`) and apply it to BOTH the observed sequence (the
input the corrector sees) and that read's truth (used for scoring). Because each
read's permutation is different, the scrambled reads share no k-mers, so the
correction graph built internally from them has no coverage consensus and cannot
recover injected errors — null recall should collapse relative to the real arm.

The permutation is a bijection applied identically to observed and truth, so the
per-read injected-error COUNT (positions where observed != truth) is preserved:
real-arm and null-arm recall are therefore directly comparable. Returns
`(scrambled_records, scrambled_truth_by_id, scrambled_observed_by_id)` mirroring
`simulate_substitution_reads` so the SAME `correct_reads_scalable` +
`per_base_metrics` path scores the null arm.
"""
function scramble_reads(records::Vector{FASTX.FASTQ.Record},
        truth_by_id::Dict{String, String}, observed_by_id::Dict{String, String},
        assigned_q::Int; seed::Int = RCA_NULL_SEED)
    scrambled_records = Vector{FASTX.FASTQ.Record}()
    scrambled_truth_by_id = Dict{String, String}()
    scrambled_observed_by_id = Dict{String, String}()
    qstr_for = (n) -> repeat(string(Char(assigned_q + 33)), n)
    for (i, rec) in enumerate(records)
        rid = FASTX.identifier(rec)
        Tc = collect(truth_by_id[rid])
        Oc = collect(observed_by_id[rid])
        perm = Random.randperm(Random.MersenneTwister(seed + i), length(Oc))
        scrambled_truth = String(Tc[perm])
        scrambled_observed = String(Oc[perm])
        scrambled_truth_by_id[rid] = scrambled_truth
        scrambled_observed_by_id[rid] = scrambled_observed
        push!(scrambled_records,
            FASTX.FASTQ.Record(rid, scrambled_observed, qstr_for(length(scrambled_observed))))
    end
    return scrambled_records, scrambled_truth_by_id, scrambled_observed_by_id
end

# === Per-base metric computation ============================================
# `per_base_metrics` lives in a pure, Mycelia-free companion so its
# classification math is unit-testable in milliseconds. See that file's header
# for the full position classification and the
# `total_edits == true_fixes + mis_fixes + over_corrections` partition invariant.
include(joinpath(@__DIR__, "rhizomorph_correction_accuracy_metrics.jl"))

# === Main sweep =============================================================

function run_accuracy_benchmark()
    smoke = _truthy(get(ENV, "MYCELIA_RCA_SMOKE", "false"))
    accession = get(ENV, "MYCELIA_RCA_ACCESSION", "NC_001416")  # Lambda phage
    smoke_len = parse(Int, get(ENV, "MYCELIA_RCA_SMOKE_GENOME_LEN", "2000"))
    # Control A lives at err=0.0 (no injected errors -> every edit is over-
    # correction). Keep it FIRST so the sweep's first cell is the specificity
    # control.
    errs = _parse_float_list(get(ENV, "MYCELIA_RCA_ERR", ""), [0.0, 0.01, 0.05, 0.10])
    readlens = _parse_int_list(get(ENV, "MYCELIA_RCA_READLEN", ""), [150])
    coverage = parse(Float64, get(ENV, "MYCELIA_RCA_COVERAGE", smoke ? "10" : "30"))
    k = parse(Int, get(ENV, "MYCELIA_RCA_K", "21"))
    assigned_q = parse(Int, get(ENV, "MYCELIA_RCA_ASSIGNED_Q", "20"))
    seed = parse(Int, get(ENV, "MYCELIA_RCA_SEED", "42"))
    scale_floor = parse(Float64, get(ENV, "MYCELIA_RCA_SCALE_FLOOR", string(SCALE_FLOOR_BASES)))
    # Minimum aggregate fraction of reads that must actually be SCORED (joined +
    # positionally comparable) for a VERDICT — guards against clean-looking
    # metrics computed on a tiny surviving subset. SMOKE-ONLY otherwise.
    min_scored_fraction = parse(Float64, get(ENV, "MYCELIA_RCA_MIN_SCORED_FRACTION", "0.5"))

    println("=== Rhizomorph Correction ACCURACY Benchmark (wired :scalable corrector) ===")
    println("Start time     : $(Dates.now())")
    println("Mode           : $(smoke ? "SMOKE (synthetic genome, no network)" : "FULL (download $accession)")")
    println("Error rates    : $errs  (substitution-only)")
    println("Read lengths   : $readlens")
    println("Coverage       : $(coverage)x")
    println("k / assigned_Q : $k / Q$assigned_q")
    println("Corrector      : mycelia_iterative_assemble with :scalable knobs (skip_solid, hard_window, soft_em, cheap_correct, doublestrand)")
    println("Scale floor    : $(scale_floor) bases")

    workdir = mktempdir(prefix = "rca_bench_")
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)

    println("\n--- Acquiring reference ---")
    refseq, ref_path,
    ref_label = acquire_reference(
        smoke = smoke, accession = accession, smoke_len = smoke_len, seed = seed, workdir = workdir)
    glen = length(refseq)
    println("Reference: $ref_label ($glen bp)")

    # Harness-local guard (the wired pipeline does NOT clamp k — it uses config.k
    # directly): only triggers in degenerate tiny-genome / short-read configs. At
    # the default k=21, readlen>=150 it never fires, and max_k=max(k,13) matches
    # the authoritative floor either way.
    min_readlen = minimum(min.(readlens, glen))
    if k > min_readlen
        @warn "k=$k exceeds shortest effective read length $min_readlen; clamping (harness-local, not in the wired pipeline)"
        k = min_readlen
    end

    rng = Random.MersenneTwister(seed)

    rows = DataFrames.DataFrame(
        reference = String[], genome_len = Int[], error_rate = Float64[],
        regime = String[], readlen = Int[], target_coverage = Float64[],
        effective_coverage = Float64[], k = Int[], assigned_q = Int[],
        ok = Bool[], n_reads = Int[], reads_scored = Int[], scored_fraction = Float64[],
        reads_dropped = Int[], reads_excluded_len = Int[], corrected_unjoined = Int[],
        total_bases = Int[], injected_errors = Int[], total_edits = Int[],
        true_fixes = Int[], mis_fixes = Int[], over_corrections = Int[],
        correct_positions = Int[], recall = Float64[], precision = Float64[],
        over_correction_rate = Float64[], correction_rate = Float64[],
        corrector_runtime_s = Float64[],
        # Control B — read-scramble null (same reads, per-read base permutation).
        null_seed = Int[], null_ok = Bool[], null_reads_scored = Int[],
        null_injected_errors = Int[], null_true_fixes = Int[],
        null_recall = Float64[], null_precision = Float64[],
        null_over_correction_rate = Float64[], null_correction_rate = Float64[],
        null_runtime_s = Float64[]
    )

    min_effective_coverage = Inf
    # Aggregate scored-fraction across cells governs the VERDICT alongside the
    # scale floor: clean-looking recall/precision on a tiny scored subset (mass
    # drops/exclusions) must NOT be quotable as validation.
    total_reads_all = 0
    total_scored_all = 0
    # Control-B aggregation (over cells with injected errors, where recall is
    # defined): pooled true_fixes / injected for the real and null arms.
    real_tp_all = 0
    real_injected_all = 0
    null_tp_all = 0
    null_injected_all = 0

    println("\n--- Sweeping (error_rate x read-regime) on the wired :scalable corrector ---")
    for err in errs
        for readlen in readlens
            regime, _tech = regime_for_readlen(readlen)
            println("\n[cell] err=$err  regime=$regime  readlen=$readlen")

            records, truth_by_id,
            observed_by_id,
            injected_total,
            sampled_bases = simulate_substitution_reads(
                refseq, readlen, coverage, err, rng; assigned_q = assigned_q)
            eff_cov = glen == 0 ? 0.0 : round(sampled_bases / glen; digits = 2)
            min_effective_coverage = min(min_effective_coverage, eff_cov)
            println("  simulated $(length(records)) reads, effective coverage $(eff_cov)x, injected $(injected_total) substitutions")

            local corrected_by_id
            t0 = time()
            ok = true
            try
                corrected_by_id, _result = correct_reads_scalable(records, k)
            catch e
                # A user interrupt (Ctrl-C during a long run) must abort, not be
                # swallowed as a "corrector failure" that records NaN and marches on.
                e isa InterruptException && rethrow()
                @warn "Corrector failed for this cell — recording NaN metrics" err readlen exception = (
                    e, catch_backtrace())
                ok = false
                corrected_by_id = Dict{String, String}()
            end
            runtime = round(time() - t0; digits = 3)

            m = per_base_metrics(truth_by_id, observed_by_id, corrected_by_id)
            total_reads_all += length(records)
            total_scored_all += m.reads_scored

            # --- Control B: read-scramble null arm ---------------------------
            # Same reads + same injected errors, but each read's bases are
            # permuted (per-read seed) before correction so no shared k-mer
            # coverage survives. Score the corrected scramble against the SAME
            # permutation of truth. Null recall should collapse vs the real arm.
            scr_records, scr_truth_by_id,
            scr_observed_by_id = scramble_reads(
                records, truth_by_id, observed_by_id, assigned_q)
            local null_corrected_by_id
            tn0 = time()
            null_ok = true
            try
                null_corrected_by_id, _null_result = correct_reads_scalable(scr_records, k)
            catch e
                e isa InterruptException && rethrow()
                @warn "Null-arm corrector failed for this cell — recording NaN null metrics" err readlen exception = (
                    e, catch_backtrace())
                null_ok = false
                null_corrected_by_id = Dict{String, String}()
            end
            null_runtime = round(time() - tn0; digits = 3)
            mn = per_base_metrics(scr_truth_by_id, scr_observed_by_id, null_corrected_by_id)
            # Mirror the real-arm ID-mismatch diagnostic on the null arm: a null
            # recall of 0 is EXPECTED (collapse), but reads_scored==0 with
            # corrected reads that did not join to truth is a read-ID FORMAT
            # mismatch (renamed ids), byte-indistinguishable from a legitimate
            # collapse in the metrics — surface it rather than let a join bug
            # masquerade as "the null correctly recovered nothing".
            if null_ok && mn.reads_scored == 0 && mn.corrected_unjoined > 0
                @warn "Null arm: no corrected read joined to truth by id, yet the corrector emitted reads — " *
                      "likely a read-ID FORMAT mismatch (renamed ids), not a genuine null collapse" err readlen corrected_unjoined = mn.corrected_unjoined
            end
            # Pool recall numerator/denominator over injected-bearing cells only.
            if m.injected > 0
                real_tp_all += m.tp
                real_injected_all += m.injected
            end
            if mn.injected > 0
                null_tp_all += mn.tp
                null_injected_all += mn.injected
            end
            # ID-format mismatch diagnostic: corrected reads exist but none
            # joined to truth => the corrector renamed ids, NOT a genuine drop.
            # Without this, every read reports as "dropped" and the operator
            # misdiagnoses a join bug as a corrector that fixed nothing.
            if ok && m.reads_scored == 0 && m.corrected_unjoined > 0
                @warn "No corrected read joined to truth by id, yet the corrector emitted reads — " *
                      "likely a read-ID FORMAT mismatch (renamed ids), not genuine drops" err readlen corrected_unjoined = m.corrected_unjoined
            end
            push!(rows,
                (
                    reference = ref_label, genome_len = glen, error_rate = err,
                    regime = regime, readlen = readlen, target_coverage = coverage,
                    effective_coverage = eff_cov, k = k, assigned_q = assigned_q,
                    ok = ok, n_reads = length(records), reads_scored = m.reads_scored,
                    scored_fraction = m.scored_fraction, reads_dropped = m.reads_dropped,
                    reads_excluded_len = m.reads_excluded_len,
                    corrected_unjoined = m.corrected_unjoined,
                    total_bases = m.total_bases, injected_errors = m.injected,
                    total_edits = m.total_edits, true_fixes = m.tp, mis_fixes = m.mis_fixes,
                    over_corrections = m.over, correct_positions = m.correct_positions,
                    recall = m.recall, precision = m.precision,
                    over_correction_rate = m.over_rate, correction_rate = m.correction_rate,
                    corrector_runtime_s = runtime,
                    null_seed = RCA_NULL_SEED, null_ok = null_ok,
                    null_reads_scored = mn.reads_scored, null_injected_errors = mn.injected,
                    null_true_fixes = mn.tp, null_recall = mn.recall,
                    null_precision = mn.precision, null_over_correction_rate = mn.over_rate,
                    null_correction_rate = mn.correction_rate, null_runtime_s = null_runtime))
            println("    ok=$ok scored=$(m.reads_scored)/$(length(records)) " *
                    "(frac=$(round(m.scored_fraction; digits=3))) " *
                    "dropped=$(m.reads_dropped) excl_len=$(m.reads_excluded_len)")
            println("    recall=$(round(m.recall; digits=4)) precision=$(round(m.precision; digits=4)) " *
                    "mis_fixes=$(m.mis_fixes) over_corr_rate=$(round(m.over_rate; digits=6)) " *
                    "correction_rate=$(round(m.correction_rate; digits=4)) (vs err=$err) $(runtime)s")
            println("    [null] scored=$(mn.reads_scored)/$(length(scr_records)) " *
                    "recall=$(round(mn.recall; digits=4)) (real=$(round(m.recall; digits=4))) " *
                    "over_corr_rate=$(round(mn.over_rate; digits=6)) $(null_runtime)s")
        end
    end

    # === Summary table ===
    println("\n=== Per-(error_rate, regime) correction accuracy ===")
    _fmt = (v) -> rpad(string(v), 12)
    println(join(
        _fmt.(["err", "regime", "readlen", "recall", "null_recall", "precision",
            "over_rate", "corr_rate", "scored"]),
        ""))
    for r in eachrow(sort(rows, [:error_rate, :readlen]))
        println(join(
            _fmt.([r.error_rate, r.regime, r.readlen,
                round(r.recall; digits = 4), round(r.null_recall; digits = 4),
                round(r.precision; digits = 4), round(r.over_correction_rate; digits = 6),
                round(r.correction_rate; digits = 4), r.reads_scored]),
            ""))
    end

    # === Control A (over-correction on un-corrupted input) ===
    # The err=0.0 rows: no errors injected, so over_correction_rate is a direct
    # specificity readout (recall/precision are NaN there by construction).
    control_a = rows[rows.error_rate .== 0.0, :]
    println("\n=== Control A: over-correction on un-corrupted input (err=0.0) ===")
    if DataFrames.nrow(control_a) == 0
        println("  (no err=0.0 cell in this sweep — set MYCELIA_RCA_ERR to include 0.0)")
    else
        for r in eachrow(control_a)
            println("  regime=$(r.regime) readlen=$(r.readlen): " *
                    "over_correction_rate=$(round(r.over_correction_rate; digits=6)) " *
                    "correction_rate=$(round(r.correction_rate; digits=6)) " *
                    "(edits=$(r.total_edits)/$(r.total_bases) bases, scored=$(r.reads_scored))")
        end
    end

    # === Control B (read-scramble null) aggregate recall comparison ===
    real_pooled_recall = real_injected_all == 0 ? NaN : real_tp_all / real_injected_all
    null_pooled_recall = null_injected_all == 0 ? NaN : null_tp_all / null_injected_all
    println("\n=== Control B: read-scramble null vs real (pooled over injected-bearing cells) ===")
    println("  real recall = $(round(real_pooled_recall; digits=4)) " *
            "($(real_tp_all)/$(real_injected_all) injected errors fixed)")
    println("  null recall = $(round(null_pooled_recall; digits=4)) " *
            "($(null_tp_all)/$(null_injected_all) injected errors fixed)")
    if isfinite(real_pooled_recall) && isfinite(null_pooled_recall)
        println("  recall collapse (real - null) = $(round(real_pooled_recall - null_pooled_recall; digits=4))")
    end

    # === Scale + scored-fraction guard ===
    # Two independent gates must BOTH pass for a VERDICT:
    #   (1) scale floor  — enough sequenced bases (coverage x genome length)
    #   (2) scored floor — enough reads actually entered the positional metric,
    #       so recall/precision don't rest on a tiny surviving subset.
    eff_cov_for_guard = isfinite(min_effective_coverage) ? min_effective_coverage : 0.0
    metric = scale_metric_bases(eff_cov_for_guard, glen)
    scale_ok = scale_verdict_allowed(eff_cov_for_guard, glen; floor = scale_floor)
    agg_scored_fraction = total_reads_all == 0 ? 0.0 : total_scored_all / total_reads_all
    scored_ok = agg_scored_fraction >= min_scored_fraction
    verdict_allowed = scale_ok && scored_ok
    mode = verdict_allowed ? "VERDICT" : "SMOKE-ONLY"

    rows[!, :mode] = fill(mode, DataFrames.nrow(rows))
    rows[!, :scale_metric_bases] = fill(metric, DataFrames.nrow(rows))
    rows[!, :scale_floor_bases] = fill(scale_floor, DataFrames.nrow(rows))
    rows[!, :agg_scored_fraction] = fill(round(agg_scored_fraction; digits = 4), DataFrames.nrow(rows))
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    csv_path = joinpath(results_dir, "rhizomorph_correction_accuracy_$(timestamp).csv")
    CSV.write(csv_path, rows)
    println("\nCSV written: $csv_path")
    println("Aggregate scored fraction: $(round(agg_scored_fraction; digits=4)) (floor $(min_scored_fraction))")

    if verdict_allowed
        println("\n=== VERDICT (scale $(round(metric; digits=0)) bases >= $(scale_floor); scored $(round(agg_scored_fraction; digits=3)) >= $(min_scored_fraction)) ===")
        println("Interpretation per cell:")
        println("  recall high + precision high + over_correction_rate ~0 => corrector fixes errors cleanly")
        println("  precision < recall via mis_fixes => failed corrections; via over_corrections => collateral damage")
        println("  correction_rate should track error_rate (edits ~ true errors); >> error_rate => over-editing")
    else
        reasons = String[]
        scale_ok ||
            push!(reasons, "scale metric $(round(metric; digits=0)) < floor $(scale_floor)")
        scored_ok || push!(reasons,
            "scored fraction $(round(agg_scored_fraction; digits=3)) < floor $(min_scored_fraction)")
        println("\n=== SMOKE-ONLY ($(join(reasons, "; "))) ===")
        println("This run MUST NOT be quoted as validation. " *
                (scale_ok ? "" : "Raise coverage/genome size. ") *
                (scored_ok ? "" :
                 "Too few reads scored (mass drops/exclusions) — investigate before trusting metrics. "))
        println("Run the full Lambda-phage config on Lovelace for a VERDICT.")
    end

    println("\n=== Accuracy benchmark complete: $(Dates.now()) ===")
    return (csv_path = csv_path, mode = mode, rows = rows,
        real_pooled_recall = real_pooled_recall, null_pooled_recall = null_pooled_recall)
end

# Only run when invoked as a script; `include`-ing this file (from the unit test)
# just defines the functions.
if abspath(PROGRAM_FILE) == @__FILE__
    run_accuracy_benchmark()
end
