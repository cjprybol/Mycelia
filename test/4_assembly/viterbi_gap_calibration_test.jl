# Multi-feature per-base correction-confidence calibration (td-21eg).
#
# The finite-gap signal only exists on a branching graph. This focused test
# drives the coverage-graph harness at small scale, verifies feature alignment,
# and enforces read-grouped train/holdout isolation.

import Test
import Mycelia
import FASTX

include(joinpath(@__DIR__, "..", "..", "benchmarking", "viterbi_gap_calibration.jl"))
if !isdefined(Main, :test_throws_message)
    include(joinpath(dirname(@__DIR__), "test_helpers.jl"))
end

Test.@testset "Grouped held-out correction-confidence calibration" begin
    synthetic_rows = [(error_rate = er, read_id = "r$(read_id)")
                      for er in (0.05, 0.10) for read_id in 1:10 for _ in 1:2]
    split_a = grouped_read_split(synthetic_rows; holdout_fraction = 0.2, seed = 9)
    split_b = grouped_read_split(synthetic_rows; holdout_fraction = 0.2, seed = 9)
    Test.@test split_a.training_indices == split_b.training_indices
    Test.@test split_a.heldout_indices == split_b.heldout_indices
    Test.@test isempty(intersect(split_a.training_groups, split_a.heldout_groups))
    Test.@test length(split_a.heldout_groups) == 4
    test_throws_message(ArgumentError, "holdout_fraction must be") do
        grouped_read_split(synthetic_rows; holdout_fraction = 1.0)
    end

    artifact = mktempdir() do dir
        result = run_gap_calibration(;
            k = 9,
            genome_length = 300,
            readlen = 80,
            coverage = 15,
            error_rates = [0.08, 0.10],
            seed = 1,
            replicates = 1,
            results_dir = dir,
            return_artifact = true
        )
        csv_path = joinpath(dir, "rhizomorph_gap_calibration.csv")
        model_path = joinpath(dir, "rhizomorph_gap_calibration_model.csv")
        manifest_path = joinpath(dir, "rhizomorph_gap_calibration_manifest.csv")
        Test.@test isfile(csv_path)
        Test.@test isfile(model_path)
        Test.@test isfile(manifest_path)
        Test.@test startswith(readline(csv_path),
            "scope,error_rate,model,n,positive_frac,auroc,ece,brier")
        Test.@test readline(model_path) == "model,term,coefficient"
        Test.@test readline(manifest_path) == "key,value"
        result
    end

    Test.@test artifact.feature_names == CORRECTION_CONFIDENCE_FEATURES
    Test.@test artifact.feature_names == (
        :raw_gap,
        :min_kmer_support,
        :all_stage0_solid,
        :competing_branch_support_ratio
    )
    Test.@test length(artifact.multifeature_model.b) == 4
    Test.@test artifact.gap_model.b isa Float64
    Test.@test artifact.metrics_path ==
               joinpath(dirname(artifact.model_path), "rhizomorph_gap_calibration.csv")
    Test.@test artifact.manifest_path == joinpath(
        dirname(artifact.model_path), "rhizomorph_gap_calibration_manifest.csv")
    Test.@test artifact.training.n > 0
    Test.@test artifact.heldout.n > 0
    Test.@test size(artifact.training.features) == (artifact.training.n, 4)
    Test.@test size(artifact.heldout.features) == (artifact.heldout.n, 4)
    Test.@test all(row.candidate_edit for row in artifact.training.rows)
    Test.@test all(row.candidate_edit for row in artifact.heldout.rows)
    Test.@test artifact.serving_config == (
        max_k = 13,
        skip_solid = false,
        graph_mode = :doublestrand,
        n_k_rungs = 3,
        max_iterations_per_k = 2,
        hard_window = false,
        soft_em = true,
        cheap_correct = true,
        beam_width = nothing,
        assigned_q = 20
    )
    Test.@test length(artifact.dataset) >= artifact.training.n + artifact.heldout.n
    Test.@test all(row.k >= 1 for row in artifact.dataset)
    Test.@test isempty(intersect(
        Set(artifact.training.groups), Set(artifact.heldout.groups)))
    for error_rate in (0.08, 0.10)
        Test.@test any(group -> group[1] == error_rate,
            artifact.split.training_groups)
        Test.@test any(group -> group[1] == error_rate,
            artifact.split.heldout_groups)
        Test.@test artifact.contract_read_passes[error_rate] > 0
        Test.@test artifact.contract_skips[error_rate] /
                   artifact.contract_read_passes[error_rate] <=
                   artifact.max_contract_skip_fraction
    end

    metric_scopes = Set((row.scope, row.error_rate) for row in artifact.rows)
    Test.@test ("pooled", nothing) in metric_scopes
    Test.@test ("error-rate", 0.08) in metric_scopes
    Test.@test ("error-rate", 0.10) in metric_scopes
    for scope in metric_scopes
        scope_rows = filter(
            row -> (row.scope, row.error_rate) == scope, artifact.rows)
        Test.@test Set(row.model for row in scope_rows) ==
                   Set(["multifeature-logistic", "gap-only-logistic"])
        for row in scope_rows
            Test.@test row.n > 0
            Test.@test isfinite(row.ece) && 0.0 <= row.ece <= 1.0
            Test.@test isfinite(row.brier) && 0.0 <= row.brier <= 1.0
            Test.@test all(bin.count > 0 for bin in row.reliability)
        end
    end
    pooled = filter(row -> row.scope == "pooled", artifact.rows)
    Test.@test all(!isnan(row.auroc) for row in pooled)
    Test.@test only(filter(row -> row.model == "gap-only-logistic", pooled)).auroc > 0.5

    test_throws_message(ArgumentError, "max_contract_skip_fraction") do
        run_gap_calibration(max_contract_skip_fraction = 1.1)
    end
end

Test.@testset "Per-base correction-confidence feature alignment" begin
    rng = Random.MersenneTwister(3)
    genome = String(rand(rng, ['A', 'C', 'G', 'T'], 500))
    k = 11
    reads = FASTX.FASTQ.Record[]
    truths = String[]
    for i in 1:40
        start = rand(rng, 1:(500 - 120 + 1))
        clean = genome[start:(start + 119)]
        observed = _inject_substitutions(clean, 0.08, rng)
        push!(reads, FASTX.FASTQ.Record(
            "r$(i)", observed, String(fill('I', 120))))
        push!(truths, clean)
    end
    fixture = build_coverage_fixture(reads, k)
    solid_kmers = Mycelia._solid_kmer_set(fixture.graph)
    config = Mycelia.ViterbiCorrectionConfig(
        record_position_gaps = true, error_rate = 0.08,
        strand_mode = :doublestrand)
    any_finite = false
    for (record, clean) in zip(reads, truths)
        observed = FASTX.sequence(String, record)
        gt = try
            collect_gap_truth(fixture, observed, clean;
                config = config, solid_kmers = solid_kmers)
        catch error
            if error isa ArgumentError && occursin("decode produced no path", error.msg)
                continue
            end
            rethrow()
        end
        n_rows = length(gt.scores)
        Test.@test length(gt.labels) == n_rows
        Test.@test length(gt.positions) == n_rows
        Test.@test length(gt.candidate_edits) == n_rows
        Test.@test length(gt.feature_rows) == n_rows
        Test.@test size(gt.features) == (n_rows, 4)
        Test.@test all(isfinite, gt.features)
        Test.@test all(gt.features[:, 2] .>= 1.0)
        Test.@test all(value -> value == 0.0 || value == 1.0, gt.features[:, 3])
        Test.@test all(value -> 0.0 <= value <= 1.0, gt.features[:, 4])

        path = only(gt.result.paths).path
        Test.@test path !== nothing
        corrected_sequence = String(Mycelia.Rhizomorph.path_to_sequence(path, fixture.graph))
        gaps = only(gt.result.paths).diagnostics[:position_gaps]
        for row_index in eachindex(gt.positions)
            gap_index = gt.positions[row_index] - k
            expected_features = Mycelia.correction_confidence_features(
                fixture.graph, path, gap_index, k, solid_kmers,
                Float64(gaps[gap_index]))
            Test.@test gt.feature_rows[row_index].raw_gap == gt.scores[row_index]
            Test.@test gt.feature_rows[row_index].min_kmer_support ==
                       expected_features.min_kmer_support
            Test.@test gt.feature_rows[row_index].all_stage0_solid ==
                       expected_features.all_stage0_solid
            Test.@test gt.feature_rows[row_index].competing_branch_support_ratio ==
                       expected_features.competing_branch_support_ratio
            Test.@test gt.candidate_edits[row_index] ==
                       (corrected_sequence[gt.positions[row_index]] !=
                        observed[gt.positions[row_index]])
            Test.@test gt.labels[row_index] ==
                       (corrected_sequence[gt.positions[row_index]] ==
                        clean[gt.positions[row_index]])
        end
        any_finite |= n_rows > 0
    end
    Test.@test any_finite
end
