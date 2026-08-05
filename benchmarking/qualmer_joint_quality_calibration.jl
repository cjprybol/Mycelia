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
    isempty(dataset_ids) && return Vector{Vector{Float64}}()
    evidence = Mycelia.Rhizomorph.get_dataset_evidence(vertex_data, first(dataset_ids))
    evidence === nothing && return Vector{Vector{Float64}}()

    vectors = Vector{Vector{Float64}}()
    for observation in values(evidence)
        for entry in observation
            if entry isa Mycelia.Rhizomorph.QualityEvidenceEntry
                push!(vectors, Float64.(Int.(entry.quality_scores) .- 33))
            end
        end
    end
    return vectors
end

phred_to_error(q) = 10.0^(-q / 10.0)

"""
Independence (the production model): Phred sums across observations, per position.
Returned UNCLAMPED, then reduced to the k-mer's weakest position — the same reduction
the traversal gate uses.
"""
function joint_independence(vectors)
    isempty(vectors) && return nothing
    k = length(first(vectors))
    return minimum(sum(v[pos] for v in vectors) for pos in 1:k)
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
            vectors = observation_quality_vectors(vertex_data)
            isempty(vectors) && continue

            sequence = String(label)
            rc = String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(sequence)))
            is_true_kmer = (sequence in truth) || (rc in truth)

            independence = joint_independence(vectors)
            conservative = joint_conservative(vectors)
            clamped = let
                ids = Mycelia.Rhizomorph.get_all_dataset_ids(vertex_data)
                value = isempty(ids) ? nothing :
                        Mycelia.Rhizomorph.get_vertex_joint_quality(vertex_data, first(ids))
                value === nothing ? nothing : Float64(minimum(value))
            end

            push!(vertex_rows,
                (chemistry = profile.name, coverage = coverage, regime = regime,
                    n_observations = length(vectors), is_true_kmer = is_true_kmer,
                    joint_independence_unclamped = round(independence; digits = 2),
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
            predicted_ind = Statistics.mean(
                phred_to_error(r.joint_independence_unclamped) for r in members)
            predicted_con = Statistics.mean(
                phred_to_error(r.joint_conservative) for r in members)
            push!(bin_rows,
                (chemistry = profile.name, coverage = coverage, regime = regime,
                    obs_bin = "$(bin[1])-$(bin[2] > 1000 ? "inf" : bin[2])",
                    n_kmers = length(members),
                    observed_error_rate = round(observed; sigdigits = 4),
                    predicted_independence = round(predicted_ind; sigdigits = 4),
                    predicted_conservative = round(predicted_con; sigdigits = 4),
                    # How many orders of magnitude the model is wrong by. Positive =
                    # OVERCONFIDENT (observed errors exceed what the model allows).
                    log10_overconfidence_independence = observed <= 0 ? -Inf :
                                                        round(
                        log10(observed / max(predicted_ind, 1e-300)); digits = 2)))
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
                        "observed=$(rpad(row.observed_error_rate, 10)) " *
                        "pred_indep=$(rpad(row.predicted_independence, 11)) " *
                        "log10_over=$(row.log10_overconfidence_independence)")
            end
        end
    end

    println("\nINTERPRETATION")
    println("  log10_over ~ 0        : independence is calibrated in that bin.")
    println("  log10_over >> 0       : summing Phred overstates confidence there.")
    println("  Compare `random` vs `artifact` at the SAME coverage: the gap between them")
    println("  is the cost of assuming independence when errors are correlated.")
    return nothing
end

calibration_main()
