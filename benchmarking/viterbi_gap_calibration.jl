# Tier-2 per-base Viterbi-gap calibration (td-4osf).
#
# Stage 2 (Route a, self-contained): drive `correct_observations` directly with
# `record_position_gaps=true` on a controlled singlestrand graph, align each
# finite per-position Viterbi gap to the corrected base, label it against truth,
# and fit + compare THREE calibration models head-to-head (isotonic, logistic,
# binned). AUROC (discrimination) is a property of the raw gap ranking — reported
# ONCE, identical across models; ECE/Brier (calibration) isolate where the models
# differ. Not a CI test (drives many whole-graph decodes); the CI-fast assertions
# live in test/4_assembly/viterbi_gap_calibration_test.jl.

import Mycelia
import Random
import FASTX

include(joinpath(@__DIR__, "calibration_metrics.jl"))
include(joinpath(@__DIR__, "benchmark_kmer_graph_fixture.jl"))

"""
    collect_gap_truth(fixture, observed_sequence, truth_sequence; config)

Decode one fixed-length observation and align `position_gaps[i]` with
`path.steps[i + 1]`. A k-mer transition contributes the newly emitted base at
read position `i + k`; its label is true iff that corrected base matches truth.
Collapsed-frontier `Inf` sentinels are excluded before calibration.
"""
function collect_gap_truth(fixture::NamedTuple,
        observed_sequence::AbstractString,
        truth_sequence::AbstractString;
        config::Mycelia.ViterbiCorrectionConfig =
        Mycelia.ViterbiCorrectionConfig(record_position_gaps = true))::NamedTuple
    length(observed_sequence) == length(truth_sequence) ||
        throw(ArgumentError("observed and truth sequences must have equal length"))
    observation = fixture.to_observation(observed_sequence)
    result = Mycelia.correct_observations(fixture.graph, [observation]; config = config)
    # NB: `something(x, throw(e))` would raise EAGERLY (throw is an argument,
    # evaluated before `something` runs) — use a lazy nothing-check so we only
    # raise when the decode genuinely failed.
    path = result.paths[1].path
    path === nothing && throw(ArgumentError("decode produced no path"))
    gaps = result.paths[1].diagnostics[:position_gaps]
    length(gaps) == length(path.steps) - 1 ||
        throw(ArgumentError("position-gap/path alignment contract violated"))
    k = length(String(path.steps[1].vertex_label))
    scores = Float64[]
    labels = Bool[]
    positions = Int[]
    for i in eachindex(gaps)
        isfinite(gaps[i]) || continue
        corrected_kmer = String(path.steps[i + 1].vertex_label)
        position = i + k
        position <= lastindex(truth_sequence) || continue
        push!(scores, gaps[i])
        push!(labels, corrected_kmer[end] == truth_sequence[position])
        push!(positions, position)
    end
    return (; scores, labels, positions, result)
end

# The three calibration models compared head-to-head. Each is a (fit, predict)
# pair over the same (gap, label) sample; the runtime gate consumes only a raw
# gap threshold, so these differ only in how they turn the gap RANKING into
# calibrated probabilities (the ECE/Brier axis).
const CALIBRATION_MODELS = (
    (name = "isotonic", fit = fit_isotonic_map, predict = predict_isotonic),
    (name = "logistic", fit = fit_logistic_map, predict = predict_logistic),
    (name = "binned", fit = fit_binned_map, predict = predict_binned)
)

"""Fit+evaluate one (fit, predict) calibration model on finite gap/truth pairs."""
function calibrate_gap_probability(fit::Function, predict::Function,
        scores::AbstractVector{<:Real}, labels::AbstractVector{Bool};
        nbins::Int = 10)::NamedTuple
    model = fit(scores, labels)
    probabilities = predict(model, scores)
    return (
        model = model,
        probabilities = probabilities,
        ece = expected_calibration_error(probabilities, labels; nbins = nbins),
        brier = brier_score(probabilities, labels),
        reliability = reliability_bins(probabilities, labels; nbins = nbins)
    )
end

"""
    compare_gap_calibrations(scores, labels; nbins=10) -> NamedTuple

Fit all three calibration models on the same finite (gap, label) sample. Returns
`(auroc, rows)` where `auroc` is the shared raw-gap discrimination (identical
across models) and each row is `(model, ece, brier, cal)`. ECE/Brier isolate
CALIBRATION quality — the only axis on which the models differ.
"""
function compare_gap_calibrations(scores::AbstractVector{<:Real},
        labels::AbstractVector{Bool}; nbins::Int = 10)::NamedTuple
    rows = [(model = m.name,
                cal = calibrate_gap_probability(
                    m.fit, m.predict, scores, labels; nbins = nbins))
            for m in CALIBRATION_MODELS]
    return (auroc = auroc(scores, labels),
        rows = [(model = r.model, ece = r.cal.ece, brier = r.cal.brier, cal = r.cal)
                for r in rows])
end

# Substitute each base with prob `rate` to a DIFFERENT base (length-preserving),
# so the observed read aligns 1:1 with its truth for per-position labeling.
function _inject_substitutions(seq::AbstractString, rate::Float64, rng)::String
    chars = collect(seq)
    alts = Dict('A' => "CGT", 'C' => "AGT", 'G' => "ACT", 'T' => "ACG")
    for i in eachindex(chars)
        if rand(rng) < rate
            opts = get(alts, chars[i], "ACGT")
            chars[i] = opts[rand(rng, 1:lastindex(opts))]
        end
    end
    return String(chars)
end

"""
    build_coverage_fixture(reads, k; moltype=:DNA) -> NamedTuple

Build a `collect_gap_truth`-compatible fixture (`.graph`, `.to_observation`) from
a READ SET rather than a single truth sequence. The read-coverage k-mer graph
carries the error- and repeat-induced BRANCHES the decoder must navigate — the
only place the per-position Viterbi margin is FINITE. A truth-only fixture graph
is linear (unique k-mers, no competitor), so every gap is `Inf` and there is no
calibration signal; the coverage graph is what the production gate actually
decodes against.
"""
function build_coverage_fixture(reads::Vector{<:FASTX.FASTQ.Record}, k::Int;
        moltype::Symbol = :DNA)::NamedTuple
    graph = Mycelia.Rhizomorph.build_kmer_graph(reads, k;
        dataset_id = "gapcal", mode = :singlestrand)
    to_observation = seq -> _bench_convert_kmers(_bench_sequence_kmers(seq, k), k, moltype)
    return (graph = graph, to_observation = to_observation)
end

"""
    run_gap_calibration(; kwargs...) -> Vector

At each error rate (through err=0.10): simulate a substitution read set from a
controlled genome, build the read-COVERAGE k-mer graph (where the Viterbi margin
is finite — see `build_coverage_fixture`), decode each read collecting finite
per-position gaps with truth labels, then fit + compare all three calibration
models. Prints a head-to-head table and writes
`results/rhizomorph_gap_calibration.csv`. Deterministic under `seed`.
"""
function run_gap_calibration(;
        k::Int = 11,
        genome_length::Int = 1200,
        readlen::Int = 120,
        coverage::Int = 30,
        error_rates::Vector{Float64} = [0.05, 0.08, 0.10],
        seed::Int = 42,
        nbins::Int = 10,
        results_dir::String = joinpath(@__DIR__, "results"))::Vector
    rng = Random.MersenneTwister(seed)
    truth_genome = String(rand(rng, ['A', 'C', 'G', 'T'], genome_length))
    n_reads = max(1, ceil(Int, coverage * genome_length / readlen))
    mkpath(results_dir)
    out_rows = NamedTuple[]
    println("=== Rhizomorph Tier-2 gap calibration (3-model head-to-head) ===")
    println("genome=$(genome_length)  k=$(k)  readlen=$(readlen)  coverage=$(coverage)x  " *
            "n_reads=$(n_reads)  seed=$(seed)")
    for er in error_rates
        # Errored read set + each read's own clean truth (substitution-only, so the
        # observed read aligns 1:1 with its truth for per-position labeling).
        reads = FASTX.FASTQ.Record[]
        truths = String[]
        for i in 1:n_reads
            st = rand(rng, 1:(genome_length - readlen + 1))
            clean = truth_genome[st:(st + readlen - 1)]
            observed = _inject_substitutions(clean, er, rng)
            push!(reads, FASTX.FASTQ.Record("r$(i)", observed, String(fill('I', readlen))))
            push!(truths, clean)
        end
        cov = build_coverage_fixture(reads, k)
        cfg = Mycelia.ViterbiCorrectionConfig(record_position_gaps = true, error_rate = er)
        scores = Float64[]
        labels = Bool[]
        for (rec, clean) in zip(reads, truths)
            gt = try
                collect_gap_truth(cov, FASTX.sequence(String, rec), clean; config = cfg)
            catch
                continue   # skip a failed/contract-violating decode
            end
            append!(scores, gt.scores)
            append!(labels, gt.labels)
        end
        if isempty(scores) || all(labels) || !any(labels)
            println("err=$(er): insufficient class balance (n=$(length(scores))) — skipped")
            continue
        end
        cmp = compare_gap_calibrations(scores, labels; nbins = nbins)
        println("err=$(er): AUROC=$(round(cmp.auroc; digits = 3))  n=$(length(scores))  " *
                "positive_frac=$(round(count(labels) / length(labels); digits = 3))")
        for r in cmp.rows
            println("   $(rpad(r.model, 9))  ECE=$(round(r.ece; digits = 4))  " *
                    "Brier=$(round(r.brier; digits = 4))")
            push!(out_rows,
                (error_rate = er, model = r.model, n = length(scores),
                    auroc = cmp.auroc, ece = r.ece, brier = r.brier))
        end
    end
    csv_path = joinpath(results_dir, "rhizomorph_gap_calibration.csv")
    open(csv_path, "w") do io
        println(io, "error_rate,model,n,auroc,ece,brier")
        for r in out_rows
            println(io, join((r.error_rate, r.model, r.n, r.auroc, r.ece, r.brier), ","))
        end
    end
    println("Wrote $(csv_path)")
    return out_rows
end

# Run when invoked directly (julia benchmarking/viterbi_gap_calibration.jl).
if abspath(PROGRAM_FILE) == @__FILE__
    run_gap_calibration()
end
