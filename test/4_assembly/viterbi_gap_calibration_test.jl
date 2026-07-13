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
            assigned_q = 21,
            seed = 1,
            replicates = 1,
            results_dir = dir,
            return_artifact = true
        )
        csv_path = result.metrics_path
        model_path = result.model_path
        audit_path = result.replicate_audit_path
        manifest_path = result.manifest_path
        Test.@test dirname(result.bundle_path) == dir
        Test.@test startswith(
            basename(result.bundle_path), "rhizomorph_gap_calibration_")
        Test.@test isfile(csv_path)
        Test.@test isfile(model_path)
        Test.@test isfile(audit_path)
        Test.@test isfile(manifest_path)
        Test.@test startswith(readline(csv_path),
            "scope,error_rate,model,n,positive_frac,auroc,ece,brier")
        Test.@test readline(model_path) == "model,term,coefficient"
        Test.@test startswith(readline(audit_path),
            "error_rate,replicate,contract_skips,read_passes")
        Test.@test readline(manifest_path) == "key,value"
        manifest_lines = readlines(manifest_path)
        Test.@test "artifact_status,publishable" in manifest_lines
        Test.@test "assigned_q,21" in manifest_lines
        Test.@test "evaluation_scope,within_cohort_read_grouped_holdout" in
                   manifest_lines
        Test.@test "whole_graph_holdout,false" in manifest_lines
        Test.@test "error_generator_distribution,bernoulli_per_base" in
                   manifest_lines
        Test.@test "frontier_simulator_parity_claimed,false" in manifest_lines
        Test.@test "feature_schema_version,td-21eg-v1" in manifest_lines
        Test.@test any(line -> startswith(line, "dataset_sha256,"), manifest_lines)
        Test.@test any(line -> startswith(line, "model_sha256,"), manifest_lines)
        Test.@test any(line -> startswith(line, "source_sha,"), manifest_lines)
        Test.@test any(
            line -> startswith(
                line, "optimizer_multifeature_logistic_final_loss,"),
            manifest_lines)
        Test.@test "logistic_l2,0.01" in manifest_lines

        # Publishing the same coherent contents again must create a distinct
        # bundle rather than overwrite or mix with the first run.
        second = _publish_calibration_bundle(
            dir, result.rows, result.multifeature_model, result.gap_model,
            result.replicate_audit, ["test_bundle" => "true"])
        Test.@test second.bundle_path != result.bundle_path
        Test.@test all(isfile, (
            second.metrics_path,
            second.model_path,
            second.replicate_audit_path,
            second.manifest_path
        ))
        result
    end

    Test.@test artifact.feature_names == CORRECTION_CONFIDENCE_FEATURES
    Test.@test artifact.feature_names == (
        :raw_gap,
        :min_kmer_support,
        :all_stage0_solid,
        :competing_branch_support_ratio,
        :collapsed_frontier
    )
    Test.@test length(artifact.multifeature_model.b) == 5
    Test.@test artifact.gap_model.b isa Vector{Float64}
    Test.@test length(artifact.gap_model.b) == 2
    Test.@test artifact.metrics_path ==
               joinpath(dirname(artifact.model_path), "rhizomorph_gap_calibration.csv")
    Test.@test artifact.manifest_path == joinpath(
        dirname(artifact.model_path), "rhizomorph_gap_calibration_manifest.csv")
    Test.@test artifact.replicate_audit_path == joinpath(
        dirname(artifact.model_path),
        "rhizomorph_gap_calibration_replicate_audit.csv")
    Test.@test artifact.bundle_path == dirname(artifact.model_path)
    Test.@test !isempty(artifact.run_id)
    Test.@test length(artifact.dataset_sha256) == 64
    Test.@test length(artifact.model_sha256) == 64
    Test.@test artifact.feature_schema_version == "td-21eg-v1"
    Test.@test artifact.evaluation_scope == :within_cohort_read_grouped_holdout
    Test.@test artifact.multifeature_model.diagnostics.converged
    Test.@test artifact.gap_model.diagnostics.converged
    Test.@test artifact.training.n > 0
    Test.@test artifact.heldout.n > 0
    Test.@test size(artifact.training.features) == (artifact.training.n, 5)
    Test.@test size(artifact.heldout.features) == (artifact.heldout.n, 5)
    Test.@test size(artifact.training.gap_features) == (artifact.training.n, 2)
    Test.@test size(artifact.heldout.gap_features) == (artifact.heldout.n, 2)
    Test.@test all(row.candidate_edit for row in artifact.training.rows)
    Test.@test all(row.candidate_edit for row in artifact.heldout.rows)
    Test.@test any(artifact.training.collapsed_frontier)
    Test.@test any(artifact.heldout.collapsed_frontier)
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
        assigned_q = 21
    )
    Test.@test artifact.assigned_q == 21
    Test.@test length(artifact.dataset) >= artifact.training.n + artifact.heldout.n
    Test.@test all(row.replicate == 1 for row in artifact.dataset)
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
    Test.@test length(artifact.replicate_audit) == 2
    Test.@test all(row -> row.read_passes > 0, artifact.replicate_audit)
    Test.@test all(row -> row.candidate_events > 0, artifact.replicate_audit)
    Test.@test all(row -> row.candidate_bearing_reads > 0,
        artifact.replicate_audit)
    Test.@test _contract_skip_fraction(1, 4, 0.25) == 0.25
    test_throws_message(ArgumentError, "exceeds configured bound") do
        _contract_skip_fraction(1, 4, 0.2)
    end
    test_throws_message(ArgumentError, "zero usable candidate events") do
        _validate_replicate_audit([(
            error_rate = 0.10,
            replicate = 2,
            contract_skips = 0,
            read_passes = 10,
            candidate_events = 0,
            candidate_bearing_reads = 0
        )])
    end
    test_throws_message(ArgumentError, "did not converge") do
        _require_converged_logistic((
                a = 0.0,
                b = [0.0],
                diagnostics = (
                    converged = false,
                    iterations = 1,
                    gradient_norm = 1.0
                )
            ), "test")
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
            Test.@test isfinite(row.auroc) && 0.0 <= row.auroc <= 1.0
            Test.@test isfinite(row.ece) && 0.0 <= row.ece <= 1.0
            Test.@test isfinite(row.brier) && 0.0 <= row.brier <= 1.0
            Test.@test all(bin.count > 0 for bin in row.reliability)
        end
    end
    pooled = filter(row -> row.scope == "pooled", artifact.rows)
    Test.@test all(!isnan(row.auroc) for row in pooled)
    Test.@test Set(row.model for row in pooled) ==
               Set(["multifeature-logistic", "gap-only-logistic"])

    test_throws_message(ArgumentError, "max_contract_skip_fraction") do
        run_gap_calibration(max_contract_skip_fraction = 1.1)
    end
    Test.@test _uniform_phred_quality_string(3, 0) == "!!!"
    Test.@test _uniform_phred_quality_string(3, 20) == "555"
    Test.@test _uniform_phred_quality_string(3, 93) == "~~~"
    assigned_record = _calibration_fastq_record("assigned-q", "ACGT", 21)
    assigned_quality = _uniform_phred_quality_string(4, 21)
    Test.@test String(FASTX.quality(assigned_record)) == assigned_quality
    Test.@test all(character -> Int(character) - 33 == 21,
        FASTX.quality(assigned_record))
    test_throws_message(ArgumentError, "supported Phred range 0:93") do
        run_gap_calibration(assigned_q = -1)
    end
    test_throws_message(ArgumentError, "supported Phred range 0:93") do
        run_gap_calibration(assigned_q = 94)
    end
    test_throws_message(ArgumentError, "needs both correctness classes") do
        _heldout_metric_rows(
            (a = 0.0, b = zeros(5)),
            (a = 0.0, b = zeros(2)),
            (
                labels = [true],
                rows = [(error_rate = 0.10,)],
                collapsed_frontier = [false],
                features = zeros(1, 5),
                gap_features = zeros(1, 2),
                scores = [0.0]
            ); nbins = 2)
    end
    test_throws_message(ArgumentError, "error-rate scope") do
        _heldout_metric_rows(
            (a = 0.0, b = zeros(5)),
            (a = 0.0, b = zeros(2)),
            (
                labels = [true, false, true, true],
                rows = [
                    (error_rate = 0.05,),
                    (error_rate = 0.05,),
                    (error_rate = 0.10,),
                    (error_rate = 0.10,)
                ],
                collapsed_frontier = falses(4),
                features = zeros(4, 5),
                gap_features = zeros(4, 2),
                scores = zeros(4)
            ); nbins = 2)
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
    checked_boundary_windows = false
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
        Test.@test size(gt.features) == (n_rows, 5)
        Test.@test all(isfinite, gt.features)
        Test.@test all(gt.features[:, 2] .>= 1.0)
        Test.@test all(value -> value == 0.0 || value == 1.0, gt.features[:, 3])
        Test.@test all(value -> 0.0 <= value <= 1.0, gt.features[:, 4])
        Test.@test all(value -> value == 0.0 || value == 1.0, gt.features[:, 5])

        path = only(gt.result.paths).path
        Test.@test path !== nothing
        corrected_sequence = String(Mycelia.Rhizomorph.path_to_sequence(path, fixture.graph))
        gaps = only(gt.result.paths).diagnostics[:position_gaps]
        checked_boundary_windows |= n_rows >= 3
        for row_index in eachindex(gt.positions)
            gap_index = gt.positions[row_index] - k
            first_step = gap_index + 1
            last_step = min(gap_index + k, length(path.steps))
            overlapping_labels = [path.steps[index].vertex_label
                                  for index in first_step:last_step]
            expected_min_support = minimum(
                Float64(Mycelia._vertex_coverage(fixture.graph[label]))
                for label in overlapping_labels)
            expected_all_solid = all(label -> label in solid_kmers,
                overlapping_labels)
            expected_branch_ratio = Float64(path.steps[first_step].probability)
            Test.@test gt.feature_rows[row_index].raw_gap == gt.scores[row_index]
            Test.@test gt.feature_rows[row_index].min_kmer_support ==
                       expected_min_support
            Test.@test gt.feature_rows[row_index].all_stage0_solid ==
                       expected_all_solid
            Test.@test gt.feature_rows[row_index].competing_branch_support_ratio ==
                       expected_branch_ratio
            Test.@test gt.feature_rows[row_index].collapsed_frontier ==
                       (gt.scores[row_index] == Inf)
            Test.@test gt.candidate_edits[row_index] ==
                       (corrected_sequence[gt.positions[row_index]] !=
                        observed[gt.positions[row_index]])
            Test.@test gt.labels[row_index] ==
                       (corrected_sequence[gt.positions[row_index]] ==
                        clean[gt.positions[row_index]])
        end
        any_finite |= any(isfinite, gt.scores)
    end
    Test.@test any_finite
    Test.@test checked_boundary_windows
end
