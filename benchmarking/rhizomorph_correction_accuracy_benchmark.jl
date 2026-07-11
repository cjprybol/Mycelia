# Rhizomorph Correction ACCURACY Benchmark — precision / over-correction / recall
# ==============================================================================
#
# Measures PER-BASE correction quality of the WIRED :scalable corrector against
# injected-error ground truth. This closes the gap left by the two existing
# harnesses:
#
#   * benchmarking/viterbi_accuracy_benchmark.jl ("B8") measures RECALL ONLY and
#     targets the bare `correct_observations` decode primitive — NOT the shipped
#     assemble_genome(corrector=:iterative, strategy=:scalable) pipeline.
#   * benchmarking/rhizomorph_correction_validation_sweep.jl measures ASSEMBLY
#     metrics (N50 / genome_fraction), not per-base correction accuracy.
#
# The central question (user direction 2026-07-11): does the wired corrector
# fix injected errors WITHOUT over-correcting (changing a correct base to a
# wrong one), and is its edit VOLUME ~ the injected error rate?
#
# WHY per-base and not assembly metrics: a prior program error (td-wlfq
# retraction) quoted an assembly-level genome-fraction "improvement" that was
# actually a read-coverage artifact. Correction quality is a per-read / per-base
# property and must be measured there, against known injected errors.
#
# CORRECTOR ENTRY POINT: this harness calls `Mycelia.mycelia_iterative_assemble`
# directly with the EXACT knobs `_corrector_strategy_knobs(:scalable)` supplies
# to the wired path (src/rhizomorph/assembly.jl:783). That function IS the
# read-corrector half of assemble_genome(corrector=:iterative, strategy=:scalable)
# — assemble_genome merely re-assembles its corrected reads into contigs (which
# would discard the per-read identity this harness needs). Corrected reads are
# read back from `result[:metadata][:final_fastq_file]` and joined to their truth
# fragments by read ID (unkmerizable/skipped reads pass through uncorrected; the
# ID-join surfaces any drop-outs explicitly rather than silently miscounting).
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
# where a true_fix is an injected position corrected to truth (C==T at O!=T).
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
#   MYCELIA_RCA_ERR             comma-separated substitution rates (default 0.01,0.05,0.10)
#   MYCELIA_RCA_READLEN         comma-separated read lengths (default 150)
#   MYCELIA_RCA_COVERAGE        target fold coverage (default 30; smoke default 10)
#   MYCELIA_RCA_K               assembly k-mer size (default 21)
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
substitutions at `error_rate` via `Mycelia.mutate_dna_substitution_fraction`
(length-preserving), and is emitted with a uniform Phred `assigned_q`.

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
Correct `records` with the WIRED :scalable knobs (exact replica of
`_corrector_strategy_knobs(:scalable)` as consumed at
src/rhizomorph/assembly.jl:879). Returns `(corrected_by_id, result_dict)` where
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
    result = Mycelia.mycelia_iterative_assemble(
        temp_fastq;
        max_k = max_k,
        skip_solid = true,
        graph_mode = :doublestrand,
        n_k_rungs = 3,
        max_iterations_per_k = 2,
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
    return corrected_by_id, result
end

# === Per-base metric computation ============================================

"""
Compute per-base correction metrics over all reads by joining `corrected_by_id`
to `truth_by_id` / `observed_by_id` on read id. Reads absent from the corrected
set (`reads_dropped`) and reads whose corrected length != truth length
(`reads_excluded_len`, a Viterbi indel path) are counted and EXCLUDED from the
positional metric rather than silently miscounted.

Position classification for a scored read (truth T, observed O, corrected C, |T|=|O|=|C|):
  injected error   : O[p] != T[p]
  edit             : C[p] != O[p]
  true_fix         : O[p] != T[p]  AND  C[p] == T[p]
  over_correction  : O[p] == T[p]  AND  C[p] != T[p]   (a correct base made wrong)
"""
function per_base_metrics(truth_by_id::Dict{String, String},
        observed_by_id::Dict{String, String}, corrected_by_id::Dict{String, String})
    tp = 0                 # true fixes
    total_edits = 0        # positions where C != O
    injected = 0           # positions where O != T
    over = 0               # over-corrections (O==T, C!=T)
    correct_positions = 0  # positions where O == T
    total_bases = 0
    reads_scored = 0
    reads_excluded_len = 0
    reads_dropped = 0
    for (rid, T) in truth_by_id
        O = observed_by_id[rid]
        if !haskey(corrected_by_id, rid)
            reads_dropped += 1
            continue
        end
        C = corrected_by_id[rid]
        # collect to Char vectors: ACGT are ASCII, but this is robust to any
        # codeunit subtlety and lets us index positionally.
        Tc = collect(T)
        Oc = collect(O)
        Cc = collect(C)
        if length(Cc) != length(Tc)
            reads_excluded_len += 1
            continue
        end
        reads_scored += 1
        L = length(Tc)
        total_bases += L
        for p in 1:L
            terr = Oc[p] != Tc[p]
            edited = Cc[p] != Oc[p]
            terr ? (injected += 1) : (correct_positions += 1)
            edited && (total_edits += 1)
            (terr && Cc[p] == Tc[p]) && (tp += 1)
            (!terr && Cc[p] != Tc[p]) && (over += 1)
        end
    end
    recall = injected == 0 ? NaN : tp / injected
    precision = total_edits == 0 ? NaN : tp / total_edits
    over_rate = correct_positions == 0 ? NaN : over / correct_positions
    correction_rate = total_bases == 0 ? NaN : total_edits / total_bases
    return (; tp, total_edits, injected, over, correct_positions, total_bases,
        recall, precision, over_rate, correction_rate,
        reads_scored, reads_excluded_len, reads_dropped)
end

# === Main sweep =============================================================

function run_accuracy_benchmark()
    smoke = _truthy(get(ENV, "MYCELIA_RCA_SMOKE", "false"))
    accession = get(ENV, "MYCELIA_RCA_ACCESSION", "NC_001416")  # Lambda phage
    smoke_len = parse(Int, get(ENV, "MYCELIA_RCA_SMOKE_GENOME_LEN", "2000"))
    errs = _parse_float_list(get(ENV, "MYCELIA_RCA_ERR", ""), [0.01, 0.05, 0.10])
    readlens = _parse_int_list(get(ENV, "MYCELIA_RCA_READLEN", ""), [150])
    coverage = parse(Float64, get(ENV, "MYCELIA_RCA_COVERAGE", smoke ? "10" : "30"))
    k = parse(Int, get(ENV, "MYCELIA_RCA_K", "21"))
    assigned_q = parse(Int, get(ENV, "MYCELIA_RCA_ASSIGNED_Q", "20"))
    seed = parse(Int, get(ENV, "MYCELIA_RCA_SEED", "42"))
    scale_floor = parse(Float64, get(ENV, "MYCELIA_RCA_SCALE_FLOOR", string(SCALE_FLOOR_BASES)))

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

    min_readlen = minimum(min.(readlens, glen))
    if k > min_readlen
        @warn "k=$k exceeds shortest effective read length $min_readlen; clamping"
        k = min_readlen
    end

    rng = Random.MersenneTwister(seed)

    rows = DataFrames.DataFrame(
        reference = String[], genome_len = Int[], error_rate = Float64[],
        regime = String[], readlen = Int[], target_coverage = Float64[],
        effective_coverage = Float64[], k = Int[], assigned_q = Int[],
        n_reads = Int[], reads_scored = Int[], reads_dropped = Int[],
        reads_excluded_len = Int[], total_bases = Int[],
        injected_errors = Int[], total_edits = Int[], true_fixes = Int[],
        over_corrections = Int[], correct_positions = Int[],
        recall = Float64[], precision = Float64[], over_correction_rate = Float64[],
        correction_rate = Float64[], corrector_runtime_s = Float64[]
    )

    min_effective_coverage = Inf

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
                @warn "Corrector failed for this cell — recording NaN metrics" err readlen exception = (
                    e, catch_backtrace())
                ok = false
                corrected_by_id = Dict{String, String}()
            end
            runtime = round(time() - t0; digits = 3)

            m = per_base_metrics(truth_by_id, observed_by_id, corrected_by_id)
            push!(rows,
                (
                    reference = ref_label, genome_len = glen, error_rate = err,
                    regime = regime, readlen = readlen, target_coverage = coverage,
                    effective_coverage = eff_cov, k = k, assigned_q = assigned_q,
                    n_reads = length(records), reads_scored = m.reads_scored,
                    reads_dropped = m.reads_dropped, reads_excluded_len = m.reads_excluded_len,
                    total_bases = m.total_bases, injected_errors = m.injected,
                    total_edits = m.total_edits, true_fixes = m.tp,
                    over_corrections = m.over, correct_positions = m.correct_positions,
                    recall = m.recall, precision = m.precision,
                    over_correction_rate = m.over_rate, correction_rate = m.correction_rate,
                    corrector_runtime_s = runtime))
            println("    ok=$ok scored=$(m.reads_scored)/$(length(records)) " *
                    "dropped=$(m.reads_dropped) excl_len=$(m.reads_excluded_len)")
            println("    recall=$(round(m.recall; digits=4)) precision=$(round(m.precision; digits=4)) " *
                    "over_corr_rate=$(round(m.over_rate; digits=6)) " *
                    "correction_rate=$(round(m.correction_rate; digits=4)) (vs err=$err) $(runtime)s")
        end
    end

    # === Summary table ===
    println("\n=== Per-(error_rate, regime) correction accuracy ===")
    _fmt = (v) -> rpad(string(v), 12)
    println(join(
        _fmt.(["err", "regime", "readlen", "recall", "precision",
            "over_rate", "corr_rate", "scored"]),
        ""))
    for r in eachrow(sort(rows, [:error_rate, :readlen]))
        println(join(
            _fmt.([r.error_rate, r.regime, r.readlen,
                round(r.recall; digits = 4), round(r.precision; digits = 4),
                round(r.over_correction_rate; digits = 6),
                round(r.correction_rate; digits = 4), r.reads_scored]),
            ""))
    end

    # === Scale guard ===
    eff_cov_for_guard = isfinite(min_effective_coverage) ? min_effective_coverage : 0.0
    metric = scale_metric_bases(eff_cov_for_guard, glen)
    verdict_allowed = scale_verdict_allowed(eff_cov_for_guard, glen; floor = scale_floor)
    mode = verdict_allowed ? "VERDICT" : "SMOKE-ONLY"

    rows[!, :mode] = fill(mode, DataFrames.nrow(rows))
    rows[!, :scale_metric_bases] = fill(metric, DataFrames.nrow(rows))
    rows[!, :scale_floor_bases] = fill(scale_floor, DataFrames.nrow(rows))
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    csv_path = joinpath(results_dir, "rhizomorph_correction_accuracy_$(timestamp).csv")
    CSV.write(csv_path, rows)
    println("\nCSV written: $csv_path")

    if verdict_allowed
        println("\n=== VERDICT (scale metric $(round(metric; digits=0)) bases >= floor $(scale_floor)) ===")
        println("Interpretation per cell:")
        println("  recall high + precision high + over_correction_rate ~0 => corrector fixes errors cleanly")
        println("  precision < recall or over_correction_rate elevated    => corrector over-corrects")
        println("  correction_rate should track error_rate (edits ~ true errors); >> error_rate => over-editing")
    else
        println("\n=== SMOKE-ONLY (scale metric $(round(metric; digits=0)) bases < floor $(scale_floor)) ===")
        println("This run is BELOW the scale floor and MUST NOT be quoted as validation.")
        println("It confirms the harness parses + runs end-to-end. Raise coverage/genome size")
        println("(or run the full Lambda-phage config on Lovelace) for a VERDICT.")
    end

    println("\n=== Accuracy benchmark complete: $(Dates.now()) ===")
    return (csv_path = csv_path, mode = mode, rows = rows)
end

# Only run when invoked as a script; `include`-ing this file (from the unit test)
# just defines the functions.
if abspath(PROGRAM_FILE) == @__FILE__
    run_accuracy_benchmark()
end
