# Is the independent-observation assumption behind joint quality calibrated? (td-2tg8)
#
# THE ASSUMPTION UNDER TEST
#   Rhizomorph aggregates per-position quality across observations by SUMMING Phred
#   (`Q_combined = Q₁ + … + Qₙ`, clamped to 255) — the `aggregate_quality_scores_independence`
#   model of `planning-docs/rhizomorph-graph-ecosystem-plan.md:498-550`, implemented as
#   `combine_phred_scores` / `get_vertex_joint_quality`. Summing in log space is
#   MULTIPLYING error probabilities, which is valid only if the observations are
#   INDEPENDENT.
#
# WHY IT MIGHT NOT BE
#   Technical artifacts recur across reads by construction — that is what makes them
#   artifacts — so their errors are correlated, not independent. Summing over correlated
#   observations OVERSTATES confidence, and overstates it fastest exactly where artifacts
#   recur, which is the case quality-aware assembly exists to catch.
#
# WHAT THIS MEASURES
#   A k-mer's joint quality is a PREDICTION: Q implies P(this k-mer is wrong) = 10^(-Q/10).
#   Ground truth is knowable here — a k-mer either does or does not occur in the reference
#   genome — so the prediction can be scored directly.
#
#     predicted error rate  = mean over k-mers of 10^(-Q_joint/10)
#     observed error rate   = fraction of k-mers absent from the reference
#
#   Calibrated  => observed ≈ predicted.
#   Overconfident => observed >> predicted. That is the independence assumption failing.
#
#   Reported per OBSERVATION-COUNT bin, because the whole question is how confidence
#   should grow with repeated observation.
#
# THREE THINGS ARE SEPARATED, because conflating them would misattribute the cause:
#
#   1. CLAMPING vs MODEL. `combine_phred_scores` clamps to 255, so a well-covered k-mer's
#      predicted error rate bottoms out at 10^-25.5 regardless of the model. This script
#      recomputes the UNCLAMPED sum from the raw evidence, so "the clamp saturated" is
#      never mistaken for "the model is overconfident".
#
#   2. RANDOM vs CORRELATED error. Two regimes are run: random substitutions only (where
#      independence is closest to true), and the same reads plus a SYSTEMATIC artifact
#      recurring at a fixed locus in every read covering it (where it is plainly false).
#      The contrast between them is the actual answer.
#
#   3. INDEPENDENCE vs CONSERVATIVE aggregation. The design specifies a conservative
#      alternative — P(all correct) = ∏P_i, Q = -10log₁₀(1-P) — which has NO implementation
#      in `src/` (audit td-8zle). It is computed here for comparison only; this script
#      does not change any production behaviour.
#
# SATURATION
#   Also reports what fraction of TRUE k-mers sit at the 255 clamp. Above that ceiling the
#   score stops discriminating and a quality gate silently degrades toward a coverage gate.
#
# USAGE
#   julia --project=. benchmarking/qualmer_joint_quality_calibration.jl
#   julia --project=. benchmarking/qualmer_joint_quality_calibration.jl --smoke
#   julia --project=. benchmarking/qualmer_joint_quality_calibration.jl --output-dir DIR

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import BioSequences
import FASTX
import StableRNGs
import Statistics
import MetaGraphsNext

# Reuse the constant-quality probe's chemistry-stratified simulator so both scripts
# describe the same read model.
PROBE_INCLUDE_ONLY = true
include(joinpath(@__DIR__, "qualmer_quality_channel_probe.jl"))

const CAL_OUTPUT_DIR = something(arg_value("--output-dir"),
    joinpath(@__DIR__, "results", "qualmer_quality_channel_probe"))

const CAL_COVERAGES = SMOKE ? [10] : [5, 10, 20, 40]
const CAL_SEED = 42
const CAL_K = 31

# === Ground truth ===

"""
Every k-mer present in the reference, in both orientations. Membership here is the
ground truth against which a k-mer's predicted error probability is scored.
"""
function reference_kmer_set(reference, k::Int)
    kmers = Set{String}()
    for i in 1:(length(reference) - k + 1)
        kmer = reference[i:(i + k - 1)]
        push!(kmers, String(kmer))
        push!(kmers, String(BioSequences.reverse_complement(kmer)))
    end
    return kmers
end

# === Per-vertex quality aggregation, recomputed from raw evidence ===

"""
Collect each observation's per-position Phred vector for one vertex, decoded from the
Phred+33 storage used by `QualityEvidenceEntry`.

Recomputed here rather than read from `get_vertex_joint_quality` so the aggregation can
be done UNCLAMPED and under more than one model. The clamped API value is reported
alongside, so clamp effects stay separable from model effects.
"""
function observation_quality_vectors(vertex_data)
    dataset_ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
    isempty(dataset_ids) && return (vectors = Vector{Vector{Float64}}(), n_entries = 0,
        n_observation_ids = 0)
    evidence = Mycelia.Rhizomorph.get_dataset_evidence(vertex_data, first(dataset_ids))
    evidence === nothing && return (vectors = Vector{Vector{Float64}}(), n_entries = 0,
        n_observation_ids = 0)

    # The evidence map is observation_id => [entries], and ONE observation id can carry
    # SEVERAL quality entries. Both counts are returned because the independence model's
    # `n` is only meaningful if it counts independent things: production
    # `get_vertex_joint_quality` sums over ENTRIES, so if these two counts diverge the
    # model is already compounding within a single observation, before any question of
    # correlation between observations arises.
    vectors = Vector{Vector{Float64}}()
    n_observation_ids = 0
    for (_, observation) in evidence
        observation_had_quality = false
        for entry in observation
            if entry isa Mycelia.Rhizomorph.QualityEvidenceEntry
                push!(vectors, Float64.(Int.(entry.quality_scores) .- 33))
                observation_had_quality = true
            end
        end
        observation_had_quality && (n_observation_ids += 1)
    end
    return (vectors = vectors, n_entries = length(vectors),
        n_observation_ids = n_observation_ids)
end

phred_to_error(q) = 10.0^(-q / 10.0)

"""
Independence (the production model): Phred sums across observations, per position.
Returned UNCLAMPED, then reduced to the k-mer's weakest position — the same reduction
the traversal gate uses.
"""
function joint_independence_positions(vectors)
    isempty(vectors) && return nothing
    k = length(first(vectors))
    return [sum(v[pos] for v in vectors) for pos in 1:k]
end

"""
WEAKEST-POSITION reduction: the k-mer scores as its worst position. This is what the
traversal gate uses, and it is a defensible *heuristic* — one bad base invalidates a
k-mer.

It is NOT the probability that the k-mer is wrong, and must not be scored as though it
were. See `union_error_probability`.
"""
weakest_position_quality(joint) = joint === nothing ? nothing : minimum(joint)

"""
UNION reduction: the probability that the k-mer is wrong is the probability that ANY of
its k positions is wrong,

    P(k-mer wrong) = 1 - prod_j (1 - p_j)

NOT `max_j p_j`, which is only a lower bound. For k=31 the two can differ by up to a
factor of k (~15 Phred points), so scoring the weakest-position value against
whole-k-mer ground truth builds a systematic overconfidence into the measurement that
has nothing to do with correlated observations. Reporting both separates the estimand
error from the modelling error.
"""
function union_error_probability(joint)
    joint === nothing && return nothing
    # Computed as -expm1(sum(log1p(-p))) rather than the literal 1 - prod(1 - p).
    # The literal form catastrophically cancels: once p < ~1e-16, (1 - p) rounds to
    # exactly 1.0 in Float64, the product is exactly 1.0, and the union probability
    # collapses to a hard 0.0 — which then looks like INFINITE overconfidence rather
    # than the very small number it actually is. log1p/expm1 stay accurate all the way
    # down, which matters here precisely because the interesting bins are the ones with
    # astronomically small predicted error.
    log_all_correct = 0.0
    for q in joint
        log_all_correct += log1p(-phred_to_error(q))
    end
    return -expm1(log_all_correct)
end

"""
Conservative aggregation (design doc; no implementation in `src/`):
P(all correct) = ∏ P_i(correct), then Q = -10log₁₀(1 - P(all correct)).

For independent observations this is very close to the sum; it diverges when any single
observation is poor, because one bad read caps the combined confidence instead of being
outvoted by good ones.
"""
function joint_conservative(vectors)
    isempty(vectors) && return nothing
    k = length(first(vectors))
    worst = Inf
    for pos in 1:k
        p_all_correct = 1.0
        for v in vectors
            p_all_correct *= (1.0 - phred_to_error(v[pos]))
        end
        p_error = max(1.0 - p_all_correct, eps(Float64))
        worst = min(worst, -10.0 * log10(p_error))
    end
    return worst
end

# === Main ===

function calibration_main()
    mkpath(CAL_OUTPUT_DIR)
    reference = load_reference(REFERENCE_FASTA)
    truth = reference_kmer_set(reference, CAL_K)
    profiles = SMOKE ? CHEMISTRIES[1:1] : CHEMISTRIES

    vertex_rows = NamedTuple[]
    bin_rows = NamedTuple[]
    saturation_rows = NamedTuple[]

    for profile in profiles, coverage in CAL_COVERAGES, regime in ("random", "artifact")
        # The `artifact` regime injects a CORRELATED error: one fixed REFERENCE locus,
        # mis-called the same way in every read covering it, always reported at low
        # quality. Its evidence count therefore matches real sequence at the same depth
        # and only its quality distinguishes it — which is the exact case the
        # independence assumption is suspected to mishandle.
        artifact_loci = regime == "artifact" ? [length(reference) ÷ 2] : nothing
        reads = simulate_read_set(reference, profile, coverage, CAL_SEED;
            artifact_loci = artifact_loci)

        records = [FASTX.FASTQ.Record("read_$(i)", read.sequence,
                       UInt8[c ? profile.q_error : profile.q_correct
                             for c in read.corrupted])
                   for (i, read) in enumerate(reads)]

        graph = Mycelia.Rhizomorph.build_qualmer_graph(records, CAL_K; mode = :doublestrand)

        for label in collect(MetaGraphsNext.labels(graph))
            vertex_data = graph[label]
            observations = observation_quality_vectors(vertex_data)
            vectors = observations.vectors
            isempty(vectors) && continue

            sequence = String(label)
            rc = String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(sequence)))
            is_true_kmer = (sequence in truth) || (rc in truth)

            joint = joint_independence_positions(vectors)
            conservative = joint_conservative(vectors)
            clamped = let
                ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
                value = isempty(ids) ? nothing :
                        Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, first(ids))
                value === nothing ? nothing : Float64(minimum(value))
            end

            push!(vertex_rows,
                (chemistry = profile.name, coverage = coverage, regime = regime,
                    n_observations = observations.n_entries,
                    n_observation_ids = observations.n_observation_ids,
                    is_true_kmer = is_true_kmer,
                    joint_independence_unclamped = round(
                        weakest_position_quality(joint); digits = 2),
                    # The estimand-correct prediction for "is this k-mer wrong".
                    p_union = union_error_probability(joint),
                    p_weakest = phred_to_error(weakest_position_quality(joint)),
                    joint_conservative = round(conservative; digits = 2),
                    joint_clamped_api = clamped === nothing ? -1.0 : clamped))
        end

        cell = filter(
            r -> r.chemistry == profile.name && r.coverage == coverage &&
                 r.regime == regime,
            vertex_rows)

        # Saturation: how much of the TRUE k-mer set is pinned at the clamp, where the
        # score can no longer discriminate.
        true_kmers = filter(r -> r.is_true_kmer, cell)
        push!(saturation_rows,
            (chemistry = profile.name, coverage = coverage, regime = regime,
                n_true_kmers = length(true_kmers),
                frac_true_at_clamp = isempty(true_kmers) ? 0.0 :
                                     round(
                    count(r -> r.joint_clamped_api >= 255.0, true_kmers) /
                    length(true_kmers); digits = 4),
                median_unclamped_joint = isempty(true_kmers) ? 0.0 :
                                         Statistics.median(
                    [r.joint_independence_unclamped for r in true_kmers])))

        # Calibration by observation-count bin.
        for bin in ((1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, 1_000_000))
            members = filter(r -> bin[1] <= r.n_observations <= bin[2], cell)
            isempty(members) && continue
            observed = count(r -> !r.is_true_kmer, members) / length(members)
            # `p_union` is the estimand-correct prediction (P(any position wrong));
            # `p_weakest` is what the traversal gate actually uses. Both are scored, so
            # the estimand mismatch is measured rather than silently inflating the
            # apparent overconfidence.
            predicted_union = Statistics.mean(r.p_union for r in members)
            predicted_weakest = Statistics.mean(r.p_weakest for r in members)
            predicted_con = Statistics.mean(
                phred_to_error(r.joint_conservative) for r in members)
            # NO silent fallback to 1.0 here. A zero id count means the count is BROKEN,
            # and defaulting it to the "everything is fine" value reports a bug as a clean
            # result -- which is exactly what happened once already. Emit -1.0 so a broken
            # count is loud.
            entry_to_id_ratio = any(r -> r.n_observation_ids == 0, members) ? -1.0 :
                                Statistics.mean(
                r.n_observations / r.n_observation_ids for r in members)
            # A zero-error bin is NOT "no information" and must not be dropped.
            #
            # Emitting `-Inf` for `observed == 0` censored 49 of 125 bins out of every
            # summary statistic, and the censoring is ASYMMETRIC in consequence: a bin
            # with no observed errors cannot falsify a model predicting ~0 errors
            # (independence), but it is precisely where a model predicting a SUBSTANTIAL
            # error rate (conservative) can be falsified. Summarising only the surviving
            # bins therefore removed the evidence against the model the ADR preferred —
            # and `-Inf` reads as the SAFE end of a scale documented as
            # "positive = overconfident", so the censoring was invisible.
            #
            # Instead: with 0 errors in N trials the one-sided 95% upper bound on the
            # true rate is ~3/N (rule of three), so the model is overconfident by AT
            # LEAST log10(bound / predicted) when predicted exceeds the bound. That is a
            # real, reportable lower bound rather than a hole in the table.
            rule_of_three = 3.0 / length(members)
            over(p) =
                if observed > 0
                    round(log10(observed / max(p, 1e-300)); digits = 2)
                elseif p > rule_of_three
                    # Overpredicts beyond what zero observed errors permits.
                    round(log10(rule_of_three / max(p, 1e-300)); digits = 2)
                else
                    # Prediction is consistent with zero observed errors; no evidence
                    # either way. Distinct from both a measurement and a censored hole.
                    0.0
                end
            # Carried so a reader can tell a measured value from a bounded one.
            bin_is_measured = observed > 0
            push!(bin_rows,
                (chemistry = profile.name, coverage = coverage, regime = regime,
                    obs_bin = "$(bin[1])-$(bin[2] > 1000 ? "inf" : bin[2])",
                    n_kmers = length(members),
                    # >1 means evidence entries outnumber observation ids, i.e. the
                    # independence model is compounding within one observation.
                    mean_entries_per_observation_id = round(
                        entry_to_id_ratio; digits = 3),
                    observed_error_rate = round(observed; sigdigits = 4),
                    # TRUE when observed > 0 (a measurement); FALSE when the bin had
                    # zero observed errors and log10_over_* is a rule-of-three LOWER
                    # BOUND instead. Any summary over these columns must say which
                    # subset it used -- summarising only measured bins is what
                    # invalidated ADR section 4.2.
                    bin_is_measured = bin_is_measured,
                    rule_of_three_upper_bound = round(rule_of_three; sigdigits = 4),
                    predicted_union = round(predicted_union; sigdigits = 4),
                    predicted_weakest = round(predicted_weakest; sigdigits = 4),
                    predicted_conservative = round(predicted_con; sigdigits = 4),
                    # How many orders of magnitude the model is wrong by. Positive =
                    # OVERCONFIDENT (observed errors exceed what the model allows, or
                    # exceed the rule-of-three bound when nothing was observed).
                    log10_over_union = over(predicted_union),
                    log10_over_weakest = over(predicted_weakest),
                    log10_over_conservative = over(predicted_con)))
        end

        @info "cell done" chemistry=profile.name coverage regime n_vertices=length(cell)
    end

    write_tsv(joinpath(CAL_OUTPUT_DIR, "qualmer_joint_quality_calibration_bins.tsv"), bin_rows)
    write_tsv(joinpath(CAL_OUTPUT_DIR, "qualmer_joint_quality_saturation.tsv"), saturation_rows)

    report_calibration(bin_rows, saturation_rows)
    return bin_rows
end

function report_calibration(bin_rows, saturation_rows)
    println("\n=== SATURATION: fraction of TRUE k-mers pinned at the 255 clamp ===")
    println("(above the clamp the score stops discriminating; the gate degrades toward coverage)")
    for row in saturation_rows
        println("  [$(rpad(row.chemistry, 8))] cov=$(lpad(row.coverage, 3))x " *
                "$(rpad(row.regime, 9)) frac_at_clamp=$(rpad(row.frac_true_at_clamp, 7)) " *
                "median_unclamped_joint=$(row.median_unclamped_joint)")
    end

    println("\n=== CALIBRATION: observed vs predicted error rate, by observation count ===")
    println("(log10_over > 0 means OVERCONFIDENT — real errors exceed what the model allows)")
    for chem in unique(r.chemistry for r in bin_rows)
        for regime in unique(r.regime for r in bin_rows)
            rows = filter(r -> r.chemistry == chem && r.regime == regime, bin_rows)
            isempty(rows) && continue
            println("\n  [$(chem) / $(regime)]")
            for row in rows
                println("    cov=$(lpad(row.coverage, 3))x obs=$(rpad(row.obs_bin, 8)) " *
                        "n=$(rpad(row.n_kmers, 6)) " *
                        "entries/id=$(rpad(row.mean_entries_per_observation_id, 6)) " *
                        "observed=$(rpad(row.observed_error_rate, 10)) " *
                        "over_union=$(rpad(row.log10_over_union, 7)) " *
                        "over_weakest=$(rpad(row.log10_over_weakest, 7)) " *
                        "over_conservative=$(row.log10_over_conservative)")
            end
        end
    end

    println("\nINTERPRETATION")
    println("  log10_over ~ 0  : that model is calibrated in that bin.")
    println("  log10_over >> 0 : that model overstates confidence there.")
    println("  over_union is the ESTIMAND-CORRECT score for 'is this k-mer wrong'.")
    println("  over_weakest is what the traversal gate actually uses; the gap between")
    println("  union and weakest is measurement error, NOT evidence about correlation.")
    println("  Compare `random` vs `artifact` at the SAME coverage: THAT gap is the cost")
    println("  of assuming independence when errors are correlated.")
    println("  entries/id > 1 means the model compounds within one observation, before")
    println("  any question of correlation BETWEEN observations arises.")
    return nothing
end

calibration_main()
