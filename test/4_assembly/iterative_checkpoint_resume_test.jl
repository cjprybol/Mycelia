# Checkpoint JSON must round-trip typed history and resume at the saved NEXT-pass
# cursor without replaying a completed pass or duplicating indel telemetry.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/iterative_checkpoint_resume_test.jl")'

import FASTX
import JSON
import Mycelia
import Test

function checkpoint_test_error(
        operation::F,
        expected_type::Type{T},
        expected_fragment::AbstractString,
)::Nothing where {F <: Function, T <: Exception}
    observed_exception = try
        operation()
        nothing
    catch exception
        exception
    end
    Test.@test typeof(observed_exception) === expected_type
    if typeof(observed_exception) === expected_type
        Test.@test Base.occursin(
            expected_fragment,
            Base.sprint(Base.showerror, observed_exception),
        )
    end
    return nothing
end

function checkpoint_frontier_metric_error(
        mutation::F,
        valid_checkpoint::AbstractDict,
        checkpoint_file::String,
        input_fastq::String,
        output_directory::String,
        expected_fragment::AbstractString,
)::Nothing where {F <: Function}
    corrupted_checkpoint = Base.deepcopy(valid_checkpoint)
    telemetry_row = Base.first(
        corrupted_checkpoint["indel_rung_telemetry"])
    mutation(telemetry_row)
    Base.open(checkpoint_file, "w") do io
        JSON.print(io, corrupted_checkpoint, 2)
    end
    checkpoint_test_error(ArgumentError, expected_fragment) do
        checkpoint_frontier_invocation(
            input_fastq, output_directory, 1)
    end
    return nothing
end

function checkpoint_resume_reads()::Vector{FASTX.FASTQ.Record}
    sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    quality = repeat("I", length(sequence))
    return FASTX.FASTQ.Record[
        FASTX.FASTQ.Record("checkpoint_1", sequence, quality),
        FASTX.FASTQ.Record("checkpoint_2", sequence, quality),
    ]
end

function checkpoint_resume_invocation(
        input_fastq::String,
        output_directory::String;
        indel_params::Union{Nothing, Mycelia.IndelDecodeParams} = nothing,
        k_ladder::Union{Nothing, Vector{Int}} = nothing,
)::Dict{Symbol, Any}
    return Mycelia.mycelia_iterative_assemble(
        input_fastq;
        max_k = 3,
        max_iterations_per_k = 2,
        improvement_threshold = 0.0,
        stop_on_no_change = false,
        graph_mode = :singlestrand,
        verbose = false,
        enable_checkpointing = true,
        checkpoint_interval = 1,
        output_dir = output_directory,
        indel_params = indel_params,
        k_ladder = k_ladder,
    )
end

function checkpoint_frontier_params()::Mycelia.IndelDecodeParams
    return Mycelia.IndelDecodeParams(
        0.05,
        0.2,
        0.2,
        0.1,
        0.1,
        3,
        3,
        16,
    )
end

function checkpoint_frontier_invocation(
        input_fastq::String,
        output_directory::String,
        max_iterations_per_k::Int,
)::Dict{Symbol, Any}
    return Mycelia.mycelia_iterative_assemble(
        input_fastq;
        max_k = 5,
        max_iterations_per_k = max_iterations_per_k,
        improvement_threshold = 0.0,
        stop_on_no_change = false,
        graph_mode = :canonical,
        verbose = false,
        enable_checkpointing = true,
        checkpoint_interval = 1,
        output_dir = output_directory,
        k_ladder = [3, 5],
        hard_window = true,
        windowed_decode = true,
        indel_params = checkpoint_frontier_params(),
        indel_schedule = :frontier_budgeted,
    )
end

Test.@testset "Fresh input hash and FASTQ parse share one descriptor" begin
    Base.mktempdir() do temporary_directory
        input_fastq = Base.joinpath(temporary_directory, "input.fastq")
        replacement_fastq = Base.joinpath(
            temporary_directory, "replacement.fastq")
        original_reads = checkpoint_resume_reads()
        replacement_reads = FASTX.FASTQ.Record[
            FASTX.FASTQ.Record("replacement", "TGCATGCA", "IIIIIIII"),
        ]
        Mycelia.write_fastq(records = original_reads, filename = input_fastq)
        Mycelia.write_fastq(
            records = replacement_reads, filename = replacement_fastq)
        original_sha256 = Mycelia._stream_sha256_file(input_fastq)

        observed_sha256,
        observed_reads = Mycelia._load_hashed_input_fastq(
            input_fastq;
            after_initial_hash = () -> Base.mv(
                replacement_fastq,
                input_fastq;
                force = true,
            ),
        )
        Test.@test observed_sha256 == original_sha256
        Test.@test [
            FASTX.sequence(String, read) for read in observed_reads
        ] == [FASTX.sequence(String, read) for read in original_reads]
        current_reads = open(FASTX.FASTQ.Reader, input_fastq) do reader
            collect(reader)
        end
        Test.@test [
            FASTX.sequence(String, read) for read in current_reads
        ] == [FASTX.sequence(String, read) for read in replacement_reads]

        Mycelia.write_fastq(records = original_reads, filename = input_fastq)
        replacement_bytes_path = Base.joinpath(
            temporary_directory, "replacement_bytes.fastq")
        Mycelia.write_fastq(
            records = replacement_reads,
            filename = replacement_bytes_path,
        )
        replacement_bytes = Base.read(replacement_bytes_path)
        checkpoint_test_error(
            ArgumentError,
            "input FASTQ changed while it was being hashed and parsed",
        ) do
            Mycelia._load_hashed_input_fastq(
                input_fastq;
                after_initial_hash = () -> Base.open(input_fastq, "w") do io
                    Base.write(io, replacement_bytes)
                end,
            )
        end
    end
end

Test.@testset "Checkpoint FASTQ rejects in-place mutation during parse" begin
    Base.mktempdir() do temporary_directory
        checkpoint_basename = "reads_k3_iter1_20260714T000000.fastq"
        checkpoint_fastq = Base.joinpath(
            temporary_directory, checkpoint_basename)
        original_reads = checkpoint_resume_reads()
        replacement_fastq = Base.joinpath(
            temporary_directory, "replacement.fastq")
        replacement_reads = FASTX.FASTQ.Record[
            FASTX.FASTQ.Record("replacement", "TGCATGCA", "IIIIIIII"),
        ]
        Mycelia.write_fastq(
            records = original_reads, filename = checkpoint_fastq)
        Mycelia.write_fastq(
            records = replacement_reads, filename = replacement_fastq)
        expected_sha256 = Mycelia._stream_sha256_file(checkpoint_fastq)
        replacement_bytes = Base.read(replacement_fastq)

        checkpoint_test_error(
            ArgumentError,
            "checkpoint current_fastq_file changed while it was being parsed",
        ) do
            Mycelia._load_validated_checkpoint_fastq(
                checkpoint_fastq,
                expected_sha256,
                Base.realpath(temporary_directory),
                checkpoint_basename;
                after_initial_hash = () -> Base.open(
                    checkpoint_fastq, "w") do io
                    Base.write(io, replacement_bytes)
                end,
            )
        end
    end
end

Test.@testset "Legacy frontier telemetry remains permissive" begin
    legacy_metric = Mycelia._symbolize_indel_frontier_metric(
        Dict{String, Any}(
            "anchored" => true,
            "reason" => "complete",
        ),
    )
    Test.@test legacy_metric == Dict{Symbol, Any}(
        :anchored => true,
        :reason => :complete,
    )

    legacy_checkpoint = Dict{String, Any}(
        "indel_rung_telemetry" => Any[
            Dict{String, Any}(
                "profile_requested" => true,
                "requested" => 5,
                "attempted" => 3,
                "completed" => 2,
                "truncated" => 1,
                "engaged" => 1,
                "graph_source" => "mixed",
                "decision_reason" => "sparse_frontier_affordable",
                "cleaned_frontier_evaluated" => 1,
                "raw_frontier_metrics" => Any[Dict{String, Any}(
                    "anchored" => true,
                    "reason" => "complete",
                )],
            ),
        ],
    )
    restored_row = Base.only(
        Mycelia._restore_indel_rung_telemetry(legacy_checkpoint))
    Test.@test restored_row[:requested] == 5
    Test.@test restored_row[:raw_frontier_evaluated] == 0
    Test.@test restored_row[:cleaned_frontier_evaluated] == 1
    Test.@test restored_row[:admitted_windows] == 0
    Test.@test restored_row[:rejected_windows] == 0
    Test.@test restored_row[:decision_reason] == :sparse_frontier_affordable
    Test.@test Base.only(restored_row[:raw_frontier_metrics]) == legacy_metric
end

Test.@testset "Schema-v2 resumes nonzero frontier telemetry exactly" begin
    Base.mktempdir() do temporary_directory
        sequence = "AAAGCTTAGGGAGAGTAGAAATAATATAGA"
        input_fastq = Base.joinpath(temporary_directory, "input.fastq")
        output_directory = Base.joinpath(temporary_directory, "run")
        Mycelia.write_fastq(
            records = FASTX.FASTQ.Record[
                FASTX.FASTQ.Record(
                    "multirung",
                    sequence,
                    repeat("I", length(sequence)),
                ),
            ],
            filename = input_fastq,
        )

        first_result = checkpoint_frontier_invocation(
            input_fastq, output_directory, 1)
        first_telemetry = Base.deepcopy(
            first_result[:metadata][:indel_rung_telemetry])
        Test.@test length(first_telemetry) == 2
        rejected_row, completed_row = first_telemetry
        Test.@test rejected_row[:k] == 3
        Test.@test rejected_row[:requested] > 0
        Test.@test rejected_row[:cleaned_frontier_evaluated] > 0
        Test.@test !isempty(rejected_row[:raw_frontier_metrics])
        Test.@test !isempty(rejected_row[:cleaned_frontier_metrics])
        Test.@test !isempty(rejected_row[:graph_cleanup])
        Test.@test completed_row[:k] == 5
        Test.@test completed_row[:requested] > 0
        Test.@test completed_row[:attempted] > 0
        Test.@test completed_row[:completed] > 0
        Test.@test completed_row[:attempted] ==
                   completed_row[:completed] + completed_row[:truncated]

        checkpoint_file = Base.joinpath(
            output_directory, "checkpoints", "latest_checkpoint.json")
        checkpoint = JSON.parsefile(checkpoint_file)
        Test.@test checkpoint["resume_configuration"]["version"] == 2
        Test.@test !isempty(
            first(checkpoint["indel_rung_telemetry"])[
                "cleaned_frontier_metrics"])

        # Schema-v2 frontier samples are exact scientific evidence, not open-ended
        # metadata. Normal probes restore every generated field and the one partial
        # producer shape (`:encode_error`) restores exactly its six context fields.
        normal_roundtrip = checkpoint_frontier_invocation(
            input_fastq, output_directory, 1)
        normal_metric = Base.first(
            Base.first(normal_roundtrip[:metadata][:indel_rung_telemetry])[
                :raw_frontier_metrics])
        Test.@test Base.Set(Base.keys(normal_metric)) ==
                   Mycelia._CHECKPOINT_REQUIRED_INDEL_METRIC_KEYS

        template_row = Base.first(checkpoint["indel_rung_telemetry"])
        template_metric = Base.first(template_row["raw_frontier_metrics"])
        truncated_sample_row = Base.deepcopy(template_row)
        truncated_sample_row["requested"] = 65
        truncated_sample_row["attempted"] = 0
        truncated_sample_row["completed"] = 0
        truncated_sample_row["truncated"] = 0
        truncated_sample_row["engaged"] = 0
        truncated_sample_row["admitted"] = false
        truncated_sample_row["admitted_windows"] = 0
        truncated_sample_row["rejected_windows"] = 65
        truncated_sample_row["graph_source"] = "substitution"
        truncated_sample_row["decision_reason"] = "frontier_budget_exceeded"
        truncated_sample_row["raw_frontier_evaluated"] = 65
        truncated_sample_row["cleaned_frontier_evaluated"] = 0
        truncated_sample_row["raw_frontier_metrics"] = Any[
            Base.deepcopy(template_metric) for _ in 1:64
        ]
        truncated_sample_row["cleaned_frontier_metrics"] = Any[]
        normalized_truncated_sample = Mycelia._normalize_indel_rung_telemetry(
            truncated_sample_row, true)
        Test.@test normalized_truncated_sample[:raw_frontier_evaluated] == 65
        Test.@test Base.length(
            normalized_truncated_sample[:raw_frontier_metrics]) == 64

        independent_samples_row = Base.deepcopy(truncated_sample_row)
        independent_samples_row["requested"] = 64
        independent_samples_row["rejected_windows"] = 64
        independent_samples_row["raw_frontier_evaluated"] = 64
        independent_samples_row["cleaned_frontier_evaluated"] = 64
        independent_samples_row["cleaned_frontier_metrics"] = Any[
            Base.deepcopy(template_metric) for _ in 1:64
        ]
        normalized_independent_samples =
            Mycelia._normalize_indel_rung_telemetry(
                independent_samples_row, true)
        Test.@test Base.length(
            normalized_independent_samples[:raw_frontier_metrics]) == 64
        Test.@test Base.length(
            normalized_independent_samples[:cleaned_frontier_metrics]) == 64

        boundary_row = Base.deepcopy(template_row)
        boundary_k = boundary_row["k"]
        minimum_span_metric = Base.deepcopy(template_metric)
        minimum_span_metric["window_start"] = 1
        minimum_span_metric["window_stop"] = boundary_k
        minimum_span_metric["window_length"] = 1
        minimum_span_metric["reason"] = "complete"
        minimum_span_metric["anchored"] = true
        minimum_span_metric["admitted"] = true
        minimum_span_metric["completed_columns"] = 1
        minimum_span_metric["frontier_area"] = 1
        minimum_span_metric["edge_expansions"] = 0
        minimum_span_metric["peak_frontier"] = 1
        minimum_span_metric["frontier_work"] = 1
        maximum_span_metric = Base.deepcopy(minimum_span_metric)
        maximum_span_metric["window_stop"] =
            Mycelia._INDEL_FRONTIER_MAX_WINDOW
        maximum_span_metric["window_length"] =
            Mycelia._INDEL_FRONTIER_MAX_WINDOW - boundary_k + 1
        maximum_span_metric["completed_columns"] =
            maximum_span_metric["window_length"]
        maximum_span_metric["frontier_area"] =
            maximum_span_metric["window_length"]
        maximum_span_metric["frontier_work"] =
            maximum_span_metric["window_length"]
        boundary_row["requested"] = 2
        boundary_row["attempted"] = 2
        boundary_row["completed"] = 2
        boundary_row["truncated"] = 0
        boundary_row["engaged"] = 0
        boundary_row["admitted"] = true
        boundary_row["admitted_windows"] = 2
        boundary_row["rejected_windows"] = 0
        boundary_row["graph_source"] = "raw"
        boundary_row["decision_reason"] = "raw_frontier_affordable"
        boundary_row["raw_frontier_evaluated"] = 2
        boundary_row["cleaned_frontier_evaluated"] = 0
        boundary_row["raw_frontier_metrics"] = Any[
            minimum_span_metric,
            maximum_span_metric,
        ]
        boundary_row["cleaned_frontier_metrics"] = Any[]
        normalized_boundaries = Mycelia._normalize_indel_rung_telemetry(
            boundary_row, true)
        Test.@test Base.first(normalized_boundaries[:raw_frontier_metrics])[
            :window_stop] == boundary_k
        Test.@test Base.last(normalized_boundaries[:raw_frontier_metrics])[
            :window_stop] == Mycelia._INDEL_FRONTIER_MAX_WINDOW

        reduced_window_checkpoint = Base.deepcopy(checkpoint)
        reduced_window_row = Base.first(
            reduced_window_checkpoint["indel_rung_telemetry"])
        reduced_window_metric = Base.first(
            reduced_window_row["raw_frontier_metrics"])
        maximum_window_length =
            reduced_window_metric["window_stop"] -
            reduced_window_metric["window_start"] -
            reduced_window_row["k"] + 2
        Test.@test maximum_window_length > 1
        reduced_window_length = maximum_window_length - 1
        reduced_window_metric["window_length"] = reduced_window_length
        reduced_window_metric["reason"] = "no_start_state"
        reduced_window_metric["anchored"] = true
        reduced_window_metric["admitted"] = false
        reduced_window_metric["completed_columns"] = 0
        reduced_window_metric["frontier_area"] = 0
        reduced_window_metric["edge_expansions"] = 0
        reduced_window_metric["peak_frontier"] = 0
        reduced_window_metric["frontier_work"] = 0
        Base.open(checkpoint_file, "w") do io
            JSON.print(io, reduced_window_checkpoint, 2)
        end
        reduced_window_roundtrip = checkpoint_frontier_invocation(
            input_fastq, output_directory, 1)
        restored_reduced_window_metric = Base.first(
            Base.first(
                reduced_window_roundtrip[:metadata][:indel_rung_telemetry],
            )[:raw_frontier_metrics],
        )
        Test.@test restored_reduced_window_metric[:window_length] ==
                   reduced_window_length
        Test.@test restored_reduced_window_metric[:reason] == :no_start_state

        empty_observation_checkpoint = Base.deepcopy(checkpoint)
        empty_observation_row = Base.first(
            empty_observation_checkpoint["indel_rung_telemetry"])
        empty_observation_metric = Base.first(
            empty_observation_row["raw_frontier_metrics"])
        empty_observation_metric["window_length"] = 0
        empty_observation_metric["reason"] = "empty_observation"
        empty_observation_metric["anchored"] = false
        empty_observation_metric["admitted"] = false
        empty_observation_metric["completed_columns"] = 0
        empty_observation_metric["frontier_area"] = 0
        empty_observation_metric["edge_expansions"] = 0
        empty_observation_metric["peak_frontier"] = 0
        empty_observation_metric["frontier_work"] = 0
        Base.open(checkpoint_file, "w") do io
            JSON.print(io, empty_observation_checkpoint, 2)
        end
        empty_observation_roundtrip = checkpoint_frontier_invocation(
            input_fastq, output_directory, 1)
        restored_empty_observation_metric = Base.first(
            Base.first(
                empty_observation_roundtrip[:metadata][:indel_rung_telemetry],
            )[:raw_frontier_metrics],
        )
        Test.@test restored_empty_observation_metric[:window_length] == 0
        Test.@test restored_empty_observation_metric[:reason] == :empty_observation

        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "checkpoint indel frontier metric lacks mandatory field reason",
        ) do telemetry_row
            telemetry_row["raw_frontier_metrics"] = Any[Dict{String, Any}()]
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "does not match schema v2; missing=frontier_work",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            Base.delete!(metric, "frontier_work")
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "window bounds must be positive and ordered",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["window_start"] = metric["window_stop"] + 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "window span exceeds $(Mycelia._INDEL_FRONTIER_MAX_WINDOW) bases",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["window_stop"] = metric["window_start"] +
                                    Mycelia._INDEL_FRONTIER_MAX_WINDOW
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "window_length must be nonnegative and not exceed " *
            "window span - k + 1",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["window_length"] += 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "frontier_work must equal saturating frontier_area + edge_expansions",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["frontier_work"] += 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "admitted does not match frontier predicate",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["admitted"] = !metric["admitted"]
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "reason=complete is inconsistent with anchor and completion fields",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["reason"] = "complete"
            metric["anchored"] = false
            metric["admitted"] = false
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "completed_columns exceeds window_length",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["completed_columns"] = metric["window_length"] + 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "frontier_area is smaller than completed_columns",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["reason"] = "work_limit"
            metric["anchored"] = true
            metric["admitted"] = false
            metric["completed_columns"] = 1
            metric["frontier_area"] = 0
            metric["edge_expansions"] =
                telemetry_row["frontier_work_limit"] + 1
            metric["frontier_work"] = metric["edge_expansions"]
            metric["peak_frontier"] = 0
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "branch_vertices exceeds vertex_count",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["branch_vertices"] = metric["vertex_count"] + 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "branch_fraction does not match branch_vertices / vertex_count",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["branch_fraction"] = metric["branch_fraction"] == 0.0 ?
                                        0.5 : 0.0
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "peak_frontier exceeds frontier_area",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["peak_frontier"] = metric["frontier_area"] + 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "reason=complete is inconsistent with anchor and completion fields",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            over_limit = telemetry_row["frontier_work_limit"] + 1
            metric["reason"] = "complete"
            metric["anchored"] = true
            metric["admitted"] = false
            metric["completed_columns"] = metric["window_length"]
            metric["frontier_area"] = metric["window_length"]
            metric["edge_expansions"] =
                over_limit - metric["frontier_area"]
            metric["frontier_work"] = over_limit
            metric["peak_frontier"] = 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "reason=complete requires a positive vertex_count",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["reason"] = "complete"
            metric["anchored"] = true
            metric["admitted"] = true
            metric["completed_columns"] = metric["window_length"]
            metric["frontier_area"] = metric["window_length"]
            metric["edge_expansions"] = 0
            metric["frontier_work"] = metric["window_length"]
            metric["peak_frontier"] = 1
            metric["vertex_count"] = 0
            metric["edge_count"] = 0
            metric["branch_vertices"] = 0
            metric["join_vertices"] = 0
            metric["branch_fraction"] = 0.0
            metric["join_fraction"] = 0.0
            metric["max_out_degree"] = 0
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "zero completed_columns requires zero frontier_area and peak_frontier",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["reason"] = "unanchored_start"
            metric["anchored"] = false
            metric["admitted"] = false
            metric["completed_columns"] = 0
            metric["frontier_area"] = 1
            metric["edge_expansions"] = 0
            metric["frontier_work"] = 1
            metric["peak_frontier"] = 0
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "reason=frontier_exhausted is inconsistent with anchor and completion fields",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["reason"] = "frontier_exhausted"
            metric["anchored"] = true
            metric["admitted"] = false
            metric["completed_columns"] = 0
            metric["frontier_area"] = 0
            metric["edge_expansions"] = 0
            metric["frontier_work"] = 0
            metric["peak_frontier"] = 0
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "reason=frontier_exhausted is inconsistent with anchor and completion fields",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            over_limit = telemetry_row["frontier_work_limit"] + 1
            metric["reason"] = "frontier_exhausted"
            metric["anchored"] = true
            metric["admitted"] = false
            metric["completed_columns"] = 1
            metric["frontier_area"] = 1
            metric["edge_expansions"] = over_limit - 1
            metric["frontier_work"] = over_limit
            metric["peak_frontier"] = 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "raw_frontier_metrics sample count does not match " *
            "raw_frontier_evaluated",
        ) do telemetry_row
            Base.empty!(telemetry_row["raw_frontier_metrics"])
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "raw frontier evaluations must equal requested windows",
        ) do telemetry_row
            telemetry_row["requested"] += 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "admitted + rejected must equal requested for a " *
            "frontier-evaluated decision",
        ) do telemetry_row
            telemetry_row["admitted_windows"] = 0
            telemetry_row["rejected_windows"] = 0
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "admitted does not match admitted_windows for a " *
            "frontier-evaluated decision",
        ) do telemetry_row
            telemetry_row["admitted"] = !telemetry_row["admitted"]
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "cleaned frontier evaluations exceed raw evaluations",
        ) do telemetry_row
            original_metric = Base.first(
                telemetry_row["raw_frontier_metrics"])
            telemetry_row["decision_reason"] = "unrestricted_semantics"
            telemetry_row["raw_frontier_evaluated"] = 0
            Base.empty!(telemetry_row["raw_frontier_metrics"])
            telemetry_row["cleaned_frontier_evaluated"] = 1
            telemetry_row["cleaned_frontier_metrics"] = Any[
                Base.deepcopy(original_metric),
            ]
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "frontier_work_limit must equal " *
            "$(Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT)",
        ) do telemetry_row
            telemetry_row["frontier_work_limit"] =
                Mycelia._DEFAULT_INDEL_FRONTIER_WORK_LIMIT - 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "frontier_metric_sample_limit must equal " *
            "$(Mycelia._INDEL_FRONTIER_TELEMETRY_SAMPLE_LIMIT)",
        ) do telemetry_row
            telemetry_row["frontier_metric_sample_limit"] =
                Mycelia._INDEL_FRONTIER_TELEMETRY_SAMPLE_LIMIT - 1
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "profile-disabled checkpoint indel telemetry must have zero counters",
        ) do telemetry_row
            telemetry_row["profile_requested"] = false
        end
        checkpoint_frontier_metric_error(
            checkpoint,
            checkpoint_file,
            input_fastq,
            output_directory,
            "reason=encode_error does not match schema v2",
        ) do telemetry_row
            metric = Base.first(telemetry_row["raw_frontier_metrics"])
            metric["reason"] = "encode_error"
            metric["anchored"] = false
            metric["admitted"] = false
        end

        saturating_checkpoint = Base.deepcopy(checkpoint)
        saturating_row = Base.first(
            saturating_checkpoint["indel_rung_telemetry"])
        saturating_metric = Base.first(
            saturating_row["raw_frontier_metrics"])
        saturating_metric["reason"] = "work_limit"
        saturating_metric["anchored"] = true
        saturating_metric["admitted"] = false
        saturating_metric["completed_columns"] = 1
        saturating_metric["frontier_area"] = typemax(Int)
        saturating_metric["edge_expansions"] = 1
        saturating_metric["frontier_work"] = typemax(Int)
        saturating_metric["peak_frontier"] = typemax(Int)
        Base.open(checkpoint_file, "w") do io
            JSON.print(io, saturating_checkpoint, 2)
        end
        saturating_roundtrip = checkpoint_frontier_invocation(
            input_fastq, output_directory, 1)
        restored_saturating_metric = Base.first(
            Base.first(saturating_roundtrip[:metadata][:indel_rung_telemetry])[
                :raw_frontier_metrics])
        Test.@test restored_saturating_metric[:frontier_area] == typemax(Int)
        Test.@test restored_saturating_metric[:edge_expansions] == 1
        Test.@test restored_saturating_metric[:peak_frontier] == typemax(Int)
        Test.@test restored_saturating_metric[:frontier_work] == typemax(Int)

        encode_error_checkpoint = Base.deepcopy(checkpoint)
        encode_error_row = Base.first(
            encode_error_checkpoint["indel_rung_telemetry"])
        original_metric = Base.first(encode_error_row["raw_frontier_metrics"])
        encode_error_row["raw_frontier_metrics"][1] = Dict{String, Any}(
            "anchored" => false,
            "reason" => "encode_error",
            "read_index" => original_metric["read_index"],
            "window_start" => original_metric["window_start"],
            "window_stop" => original_metric["window_stop"],
            "admitted" => false,
        )
        Base.open(checkpoint_file, "w") do io
            JSON.print(io, encode_error_checkpoint, 2)
        end
        encode_error_roundtrip = checkpoint_frontier_invocation(
            input_fastq, output_directory, 1)
        restored_encode_error = Base.first(
            Base.first(encode_error_roundtrip[:metadata][:indel_rung_telemetry])[
                :raw_frontier_metrics])
        Test.@test Base.Set(Base.keys(restored_encode_error)) ==
                   Mycelia._CHECKPOINT_ENCODE_ERROR_INDEL_METRIC_KEYS
        Test.@test restored_encode_error[:reason] == :encode_error
        Test.@test !restored_encode_error[:anchored]
        Test.@test !restored_encode_error[:admitted]

        # Restore the producer-written checkpoint before exercising cursor resume.
        Base.open(checkpoint_file, "w") do io
            JSON.print(io, checkpoint, 2)
        end
        checkpoint["next_k"] = checkpoint["current_k"]
        checkpoint["next_iteration"] = 2
        checkpoint["run_complete"] = false
        checkpoint["resume_configuration"]["max_iterations_per_k"] = 2
        Base.open(checkpoint_file, "w") do io
            JSON.print(io, checkpoint, 2)
        end

        resumed_result = checkpoint_frontier_invocation(
            input_fastq, output_directory, 2)
        resumed_metadata = resumed_result[:metadata]
        resumed_telemetry = resumed_metadata[:indel_rung_telemetry]
        Test.@test length(resumed_telemetry) == 3
        Test.@test resumed_telemetry[1:2] == first_telemetry
        Test.@test resumed_telemetry[3][:k] == 5
        Test.@test resumed_telemetry[3][:iteration] == 2
        for row in resumed_telemetry
            Test.@test row[:attempted] == row[:completed] + row[:truncated]
            Test.@test row[:engaged] <= row[:completed]
        end
        Test.@test resumed_metadata[:indel_requested] ==
                   sum(Int(row[:requested]) for row in resumed_telemetry)
        Test.@test resumed_metadata[:indel_attempted] ==
                   sum(Int(row[:attempted]) for row in resumed_telemetry)
        Test.@test resumed_metadata[:indel_completed] ==
                   sum(Int(row[:completed]) for row in resumed_telemetry)
        Test.@test resumed_metadata[:indel_truncated] ==
                   sum(Int(row[:truncated]) for row in resumed_telemetry)
        Test.@test resumed_metadata[:indel_engaged] ==
                   sum(Int(row[:engaged]) for row in resumed_telemetry)

        completed_result = checkpoint_frontier_invocation(
            input_fastq, output_directory, 2)
        completed_metadata = completed_result[:metadata]
        Test.@test completed_metadata[:iteration_history] ==
                   resumed_metadata[:iteration_history]
        Test.@test completed_metadata[:indel_rung_telemetry] == resumed_telemetry
        Test.@test completed_metadata[:final_fastq_file] ==
                   resumed_metadata[:final_fastq_file]
        Test.@test completed_result[:final_assembly] ==
                   resumed_result[:final_assembly]
    end
end

Test.@testset "Iterative checkpoint resumes the next pass exactly" begin
    Base.mktempdir() do temporary_directory
        input_fastq = Base.joinpath(temporary_directory, "input.fastq")
        output_directory = Base.joinpath(temporary_directory, "run")
        Mycelia.write_fastq(
            records = checkpoint_resume_reads(),
            filename = input_fastq,
        )

        first_result = Mycelia.mycelia_iterative_assemble(
            input_fastq;
            max_k = 3,
            max_iterations_per_k = 1,
            improvement_threshold = 0.0,
            stop_on_no_change = false,
            graph_mode = :singlestrand,
            verbose = false,
            enable_checkpointing = true,
            checkpoint_interval = 1,
            output_dir = output_directory,
        )
        first_metadata = first_result[:metadata]
        checkpoint_file = Base.joinpath(
            output_directory, "checkpoints", "latest_checkpoint.json")
        checkpoint = JSON.parsefile(checkpoint_file)

        Test.@test checkpoint["run_complete"]
        Test.@test checkpoint["current_iteration"] == 1
        Test.@test checkpoint["resume_configuration"]["version"] == 2
        Test.@test checkpoint["resume_configuration"][
            "substitution_error_rate"] === nothing
        Test.@test checkpoint["resume_configuration"]["input_fastq_path"] ==
                   Base.realpath(input_fastq)
        Test.@test checkpoint["resume_configuration"]["input_fastq_sha256"] ==
                   Mycelia._stream_sha256_file(input_fastq)
        Test.@test checkpoint["corrector_diagnostics"][
            "window_anchor_rejections"] == 0
        Test.@test checkpoint["corrector_diagnostics"][
            "substitution_length_divergences"] == 0
        Test.@test first(checkpoint["iteration_history"]["3"]) isa
                   Dict{String, Any}
        Test.@test first_metadata[:k_progression] == [3]
        Test.@test length(first_metadata[:indel_rung_telemetry]) == 1

        # Simulate an interruption after pass 1 of a two-pass run. The checkpoint
        # retains pass-1 history and telemetry; only its explicit continuation
        # cursor changes.
        checkpoint["next_k"] = checkpoint["current_k"]
        checkpoint["next_iteration"] = 2
        checkpoint["run_complete"] = false
        checkpoint["resume_configuration"]["max_iterations_per_k"] = 2
        checkpoint["corrector_diagnostics"][
            "substitution_length_divergences"] = 2
        open(checkpoint_file, "w") do io
            JSON.print(io, checkpoint, 2)
        end

        resumed_result = Mycelia.mycelia_iterative_assemble(
            input_fastq;
            max_k = 3,
            max_iterations_per_k = 2,
            improvement_threshold = 0.0,
            stop_on_no_change = false,
            graph_mode = :singlestrand,
            verbose = false,
            enable_checkpointing = true,
            checkpoint_interval = 1,
            output_dir = output_directory,
        )
        resumed_metadata = resumed_result[:metadata]
        resumed_history = resumed_metadata[:iteration_history][3]
        resumed_telemetry = resumed_metadata[:indel_rung_telemetry]

        Test.@test resumed_metadata[:k_progression] == [3]
        Test.@test [Int(row[:iteration]) for row in resumed_history] == [1, 2]
        Test.@test length(resumed_telemetry) == 2
        Test.@test [
            (
                Int(row[:ladder_index]),
                Int(row[:k]),
                Int(row[:iteration]),
            ) for row in resumed_telemetry
        ] == [(1, 3, 1), (1, 3, 2)]
        Test.@test resumed_metadata[:indel_requested] ==
                   sum(Int(row[:requested]) for row in resumed_telemetry)
        Test.@test resumed_metadata[:corrector_errors][
            :substitution_length_divergences] == 2
        Test.@test Base.isfile(resumed_metadata[:final_fastq_file])
        Test.@test Base.occursin(
            "_iter2_", Base.basename(resumed_metadata[:final_fastq_file]))
        Test.@test resumed_result isa Dict{Symbol, Any}

        # A completed checkpoint is idempotent: reopening it does not add a third
        # pass or shift the telemetry cursor.
        completed_result = Mycelia.mycelia_iterative_assemble(
            input_fastq;
            max_k = 3,
            max_iterations_per_k = 2,
            improvement_threshold = 0.0,
            stop_on_no_change = false,
            graph_mode = :singlestrand,
            verbose = false,
            enable_checkpointing = true,
            checkpoint_interval = 1,
            output_dir = output_directory,
        )
        completed_metadata = completed_result[:metadata]
        Test.@test completed_metadata[:k_progression] == [3]
        Test.@test length(completed_metadata[:iteration_history][3]) == 2
        Test.@test completed_metadata[:indel_rung_telemetry] == resumed_telemetry
        Test.@test completed_metadata[:total_runtime] ==
                   resumed_metadata[:total_runtime]
        Test.@test completed_metadata[:final_fastq_file] ==
                   resumed_metadata[:final_fastq_file]
        Test.@test completed_result[:final_assembly] ==
                   resumed_result[:final_assembly]

        valid_checkpoint = JSON.parsefile(checkpoint_file)

        # Schema-v2 checkpoints are complete snapshots, not best-effort legacy
        # records. Missing cursor/root or mandatory per-pass fields fail closed.
        partial_root = Base.deepcopy(valid_checkpoint)
        delete!(partial_root, "next_iteration")
        open(checkpoint_file, "w") do io
            JSON.print(io, partial_root, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint root does not match schema v2",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        missing_profile_disabled_telemetry = Base.deepcopy(valid_checkpoint)
        delete!(
            missing_profile_disabled_telemetry,
            "indel_rung_telemetry",
        )
        open(checkpoint_file, "w") do io
            JSON.print(io, missing_profile_disabled_telemetry, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint root does not match schema v2",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        partial_row = Base.deepcopy(valid_checkpoint)
        delete!(first(partial_row["indel_rung_telemetry"]), "requested")
        open(checkpoint_file, "w") do io
            JSON.print(io, partial_row, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint indel telemetry row lacks mandatory fields: requested",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        missing_full_row_field = Base.deepcopy(valid_checkpoint)
        delete!(
            first(missing_full_row_field["indel_rung_telemetry"]),
            "graph_cleanup",
        )
        open(checkpoint_file, "w") do io
            JSON.print(io, missing_full_row_field, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint indel telemetry row lacks mandatory fields: graph_cleanup",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        invalid_root_timestamp = Base.deepcopy(valid_checkpoint)
        invalid_root_timestamp["timestamp"] = 0
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_root_timestamp, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint timestamp must be a string",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        missing_diagnostic = Base.deepcopy(valid_checkpoint)
        delete!(missing_diagnostic["corrector_diagnostics"], "structural")
        open(checkpoint_file, "w") do io
            JSON.print(io, missing_diagnostic, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint corrector_diagnostics must contain exactly",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        invalid_diagnostic_type = Base.deepcopy(valid_checkpoint)
        invalid_diagnostic_type["corrector_diagnostics"]["structural"] = 0.0
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_diagnostic_type, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint corrector_diagnostics structural must be an integer",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        mismatched_profile_request = Base.deepcopy(valid_checkpoint)
        first(mismatched_profile_request["indel_rung_telemetry"])[
            "profile_requested"] = true
        open(checkpoint_file, "w") do io
            JSON.print(io, mismatched_profile_request, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint indel telemetry profile_requested does not match invocation indel_params",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        # Every aggregate uses checked arithmetic; individually valid counters
        # cannot wrap while their checkpoint total is reconstructed.
        improvement_overflow = Base.deepcopy(valid_checkpoint)
        improvement_rows = improvement_overflow["iteration_history"]["3"]
        improvement_rows[1]["improvements_made"] = typemax(Int)
        improvement_rows[2]["improvements_made"] = 1
        improvement_overflow["total_improvements"] = typemax(Int)
        open(checkpoint_file, "w") do io
            JSON.print(io, improvement_overflow, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint total_improvements exceeds Int range",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        window_counter_overflow = Base.deepcopy(valid_checkpoint)
        overflow_row = first(window_counter_overflow["indel_rung_telemetry"])
        overflow_row["requested"] = typemax(Int)
        overflow_row["admitted_windows"] = typemax(Int)
        overflow_row["rejected_windows"] = 1
        open(checkpoint_file, "w") do io
            JSON.print(io, window_counter_overflow, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint indel telemetry admitted + rejected exceeds Int range",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        # A valid, hash-matching FASTQ from another completed pass in the same
        # output directory cannot be substituted for the checkpoint cursor.
        first_pass_timestamp = first(
            valid_checkpoint["iteration_history"]["3"])["timestamp"]
        wrong_pass_fastq = Base.joinpath(
            output_directory,
            "reads_k3_iter1_$(first_pass_timestamp).fastq",
        )
        Test.@test Base.isfile(wrong_pass_fastq)
        wrong_pass_checkpoint = Base.deepcopy(valid_checkpoint)
        wrong_pass_checkpoint["current_fastq_file"] = wrong_pass_fastq
        wrong_pass_checkpoint["current_fastq_sha256"] =
            Mycelia._stream_sha256_file(wrong_pass_fastq)
        open(checkpoint_file, "w") do io
            JSON.print(io, wrong_pass_checkpoint, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint current_fastq_file does not match the completed pass",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        # The output directory is bound to both the canonical path and content of
        # the original input. Even byte-identical reads at another path cannot
        # inherit a stale run.
        wrong_input_fastq = Base.joinpath(
            temporary_directory, "wrong_input.fastq")
        Mycelia.write_fastq(
            records = checkpoint_resume_reads(),
            filename = wrong_input_fastq,
        )
        checkpoint_test_error(
            ArgumentError,
            "checkpoint resume_configuration does not match this invocation: input_fastq_path",
        ) do
            checkpoint_resume_invocation(
                wrong_input_fastq, output_directory)
        end

        original_input_bytes = Base.read(input_fastq)
        Mycelia.write_fastq(
            records = FASTX.FASTQ.Record[
                FASTX.FASTQ.Record("changed", "ACGT", "IIII"),
            ],
            filename = input_fastq,
        )
        checkpoint_test_error(
            ArgumentError,
            "checkpoint resume_configuration does not match this invocation: input_fastq_sha256",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end
        Base.write(input_fastq, original_input_bytes)

        excessive_metric_samples = Base.deepcopy(valid_checkpoint)
        first(excessive_metric_samples["indel_rung_telemetry"])[
            "raw_frontier_metrics"] = Any[
            Dict{String, Any}(
                "anchored" => false,
                "reason" => "unanchored_start",
            ) for _ in 1:65
        ]
        open(checkpoint_file, "w") do io
            JSON.print(io, excessive_metric_samples, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint indel telemetry raw_frontier_metrics exceeds the sample limit",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        # A checkpoint-selected FASTQ must resolve inside the canonical output
        # directory, even if its forged digest is otherwise correct.
        outside_fastq = Base.joinpath(temporary_directory, "outside.fastq")
        Mycelia.write_fastq(
            records = checkpoint_resume_reads(),
            filename = outside_fastq,
        )
        outside_checkpoint = Base.deepcopy(valid_checkpoint)
        outside_checkpoint["current_fastq_file"] = outside_fastq
        outside_checkpoint["current_fastq_sha256"] =
            Mycelia._stream_sha256_file(outside_fastq)
        open(checkpoint_file, "w") do io
            JSON.print(io, outside_checkpoint, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint current_fastq_file must resolve as a direct child of output_dir",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        # Hash success is insufficient: the resumed FASTQ is parsed through the
        # same descriptor and malformed records fail closed.
        malformed_fastq = String(valid_checkpoint["current_fastq_file"])
        valid_final_fastq_bytes = Base.read(malformed_fastq)
        Base.write(malformed_fastq, "@bad\nACGT\n+\nIII\n")
        malformed_fastq_checkpoint = Base.deepcopy(valid_checkpoint)
        malformed_fastq_checkpoint["current_fastq_sha256"] =
            Mycelia._stream_sha256_file(malformed_fastq)
        open(checkpoint_file, "w") do io
            JSON.print(io, malformed_fastq_checkpoint, 2)
        end
        checkpoint_test_error(
            ErrorException,
            "Length of quality must be identical to length of sequence",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end
        Base.write(malformed_fastq, valid_final_fastq_bytes)

        # All serialized names are whitelisted before conversion to Symbol, so an
        # arbitrary JSON key can never be interned into the process.
        unsupported_history_field = Base.deepcopy(valid_checkpoint)
        first(unsupported_history_field["iteration_history"]["3"])[
            "attacker_controlled_symbol"] = 1
        open(checkpoint_file, "w") do io
            JSON.print(io, unsupported_history_field, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint iteration history row contains too many fields",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        unsupported_metric_field = Base.deepcopy(valid_checkpoint)
        first(unsupported_metric_field["indel_rung_telemetry"])[
            "raw_frontier_metrics"] = Any[
            Dict{String, Any}(
                "anchored" => true,
                "reason" => "complete",
                "attacker_controlled_symbol" => 1,
            ),
        ]
        open(checkpoint_file, "w") do io
            JSON.print(io, unsupported_metric_field, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint indel frontier metric contains unsupported field " *
            "attacker_controlled_symbol",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        # The fixed-toy anchor-loss decision is a supported, round-trippable
        # production enum rather than an arbitrary Symbol conversion.
        anchor_loss_checkpoint = Base.deepcopy(valid_checkpoint)
        first(anchor_loss_checkpoint["indel_rung_telemetry"])[
            "decision_reason"] = "cleaning_lost_window_anchor"
        open(checkpoint_file, "w") do io
            JSON.print(io, anchor_loss_checkpoint, 2)
        end
        anchor_loss_result = checkpoint_resume_invocation(
            input_fastq, output_directory)
        Test.@test first(anchor_loss_result[:metadata][:indel_rung_telemetry])[
            :decision_reason] == :cleaning_lost_window_anchor

        # Exact integer restoration rejects Bool and negative auxiliary counters.
        boolean_diagnostic = Base.deepcopy(valid_checkpoint)
        boolean_diagnostic["corrector_diagnostics"]["structural"] = true
        open(checkpoint_file, "w") do io
            JSON.print(io, boolean_diagnostic, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint corrector_diagnostics structural must be an integer",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        negative_auxiliary = Base.deepcopy(valid_checkpoint)
        negative_auxiliary["cheap_correction_counts"][1] = -1
        open(checkpoint_file, "w") do io
            JSON.print(io, negative_auxiliary, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint cheap_correction_counts[1] must be nonnegative",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        # Size is checked from the open checkpoint descriptor before JSON parsing.
        open(checkpoint_file, "w") do io
            Base.truncate(io, Mycelia._ITERATIVE_CHECKPOINT_MAX_BYTES + 1)
        end
        checkpoint_test_error(
            ArgumentError,
            "maximum is $(Mycelia._ITERATIVE_CHECKPOINT_MAX_BYTES) bytes",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end
        open(checkpoint_file, "w") do io
            JSON.print(io, valid_checkpoint, 2)
        end

        # Resume-critical correction settings are provenance, not caller-local
        # display metadata. A mismatch must fail before any pass is replayed.
        mismatched_configuration = Base.deepcopy(valid_checkpoint)
        mismatched_configuration["resume_configuration"]["indel_schedule"] =
            "frontier_budgeted"
        open(checkpoint_file, "w") do io
            JSON.print(io, mismatched_configuration, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint resume_configuration does not match this invocation: indel_schedule",
        ) do
            Mycelia.mycelia_iterative_assemble(
                input_fastq;
                max_k = 3,
                max_iterations_per_k = 2,
                improvement_threshold = 0.0,
                stop_on_no_change = false,
                graph_mode = :singlestrand,
                verbose = false,
                enable_checkpointing = true,
                checkpoint_interval = 1,
                output_dir = output_directory,
            )
        end

        mismatched_substitution_rate = Base.deepcopy(valid_checkpoint)
        mismatched_substitution_rate["resume_configuration"][
            "substitution_error_rate"] = 0.001
        open(checkpoint_file, "w") do io
            JSON.print(io, mismatched_substitution_rate, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint resume_configuration does not match this invocation: " *
            "substitution_error_rate",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        boolean_as_integer_configuration = Base.deepcopy(valid_checkpoint)
        boolean_as_integer_configuration["resume_configuration"][
            "stop_on_no_change"] = 0
        open(checkpoint_file, "w") do io
            JSON.print(io, boolean_as_integer_configuration, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint resume_configuration does not match this invocation: stop_on_no_change",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        integer_as_float_configuration = Base.deepcopy(valid_checkpoint)
        integer_as_float_configuration["resume_configuration"]["max_k"] = 3.0
        open(checkpoint_file, "w") do io
            JSON.print(io, integer_as_float_configuration, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint resume_configuration does not match this invocation: max_k",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory)
        end

        # JSON arrays deserialize to a different concrete container type than the
        # invocation's Vector{Int}; equal scalar leaves remain compatible.
        valid_ladder_configuration = Base.deepcopy(valid_checkpoint)
        valid_ladder_configuration["resume_configuration"]["k_ladder"] = Any[3]
        valid_ladder_configuration["resume_configuration"][
            "resolved_k_schedule"] = Any[3]
        open(checkpoint_file, "w") do io
            JSON.print(io, valid_ladder_configuration, 2)
        end
        ladder_result = checkpoint_resume_invocation(
            input_fastq, output_directory; k_ladder = [3])
        Test.@test ladder_result[:metadata][:k_progression] == [3]

        invalid_nested_ladder_type = Base.deepcopy(valid_ladder_configuration)
        invalid_nested_ladder_type["resume_configuration"]["k_ladder"] =
            Any[3.0]
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_nested_ladder_type, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint resume_configuration does not match this invocation: k_ladder",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory; k_ladder = [3])
        end

        indel_params = Mycelia.IndelDecodeParams(
            0.05, 0.2, 0.2, 0.1, 0.1, 3, 3, 16)
        invalid_nested_indel_type = Base.deepcopy(valid_checkpoint)
        invalid_nested_indel_type["resume_configuration"]["indel_params"] =
            Dict{String, Any}(
                "base_error_rate" => indel_params.base_error_rate,
                "insertion_fraction" => indel_params.insertion_fraction,
                "deletion_fraction" => indel_params.deletion_fraction,
                "insertion_extend_probability" =>
                    indel_params.insertion_extend_probability,
                "deletion_extend_probability" =>
                    indel_params.deletion_extend_probability,
                "deletion_max_run" => 3.0,
                "max_insertion_run" => indel_params.max_insertion_run,
                "band_width" => indel_params.band_width,
            )
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_nested_indel_type, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint resume_configuration does not match this invocation: indel_params",
        ) do
            checkpoint_resume_invocation(
                input_fastq, output_directory; indel_params = indel_params)
        end

        # Existing-but-corrupt JSON is never treated as permission to restart
        # fresh and silently replay completed work.
        Base.write(checkpoint_file, "{")
        checkpoint_test_error(
            ArgumentError,
            "is invalid; refusing an implicit fresh restart",
        ) do
            Mycelia.mycelia_iterative_assemble(
                input_fastq;
                max_k = 3,
                max_iterations_per_k = 2,
                improvement_threshold = 0.0,
                stop_on_no_change = false,
                graph_mode = :singlestrand,
                verbose = false,
                enable_checkpointing = true,
                checkpoint_interval = 1,
                output_dir = output_directory,
            )
        end
        open(checkpoint_file, "w") do io
            JSON.print(io, valid_checkpoint, 2)
        end

        invalid_fastq_hash = Base.deepcopy(valid_checkpoint)
        invalid_fastq_hash["current_fastq_sha256"] = repeat("0", 64)
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_fastq_hash, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint current_fastq_file SHA-256 does not match saved provenance",
        ) do
            Mycelia.mycelia_iterative_assemble(
                input_fastq;
                max_k = 3,
                max_iterations_per_k = 2,
                improvement_threshold = 0.0,
                stop_on_no_change = false,
                graph_mode = :singlestrand,
                verbose = false,
                enable_checkpointing = true,
                checkpoint_interval = 1,
                output_dir = output_directory,
            )
        end

        # A corrupt explicit cursor must not replay pass 2.
        invalid_cursor = Base.deepcopy(valid_checkpoint)
        invalid_cursor["run_complete"] = false
        invalid_cursor["next_k"] = invalid_cursor["current_k"]
        invalid_cursor["next_iteration"] = invalid_cursor["current_iteration"]
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_cursor, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "same-rung checkpoint next_iteration must equal",
        ) do
            Mycelia.mycelia_iterative_assemble(
                input_fastq;
                max_k = 3,
                max_iterations_per_k = 2,
                improvement_threshold = 0.0,
                stop_on_no_change = false,
                graph_mode = :singlestrand,
                verbose = false,
                enable_checkpointing = true,
                checkpoint_interval = 1,
                output_dir = output_directory,
            )
        end

        # Impossible telemetry cannot seed resumed aggregate diagnostics.
        invalid_telemetry = Base.deepcopy(valid_checkpoint)
        invalid_row = first(invalid_telemetry["indel_rung_telemetry"])
        invalid_row["requested"] = 1
        invalid_row["attempted"] = 3
        invalid_row["completed"] = 2
        invalid_row["truncated"] = 2
        invalid_row["engaged"] = 5
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_telemetry, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint indel telemetry attempted=3 exceeds requested=1",
        ) do
            Mycelia.mycelia_iterative_assemble(
                input_fastq;
                max_k = 3,
                max_iterations_per_k = 2,
                improvement_threshold = 0.0,
                stop_on_no_change = false,
                graph_mode = :singlestrand,
                verbose = false,
                enable_checkpointing = true,
                checkpoint_interval = 1,
                output_dir = output_directory,
            )
        end

        invalid_metric = Base.deepcopy(valid_checkpoint)
        first(invalid_metric["indel_rung_telemetry"])["raw_frontier_metrics"] =
            Any[42]
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_metric, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint indel telemetry raw_frontier_metrics contains a non-object metric",
        ) do
            Mycelia.mycelia_iterative_assemble(
                input_fastq;
                max_k = 3,
                max_iterations_per_k = 2,
                improvement_threshold = 0.0,
                stop_on_no_change = false,
                graph_mode = :singlestrand,
                verbose = false,
                enable_checkpointing = true,
                checkpoint_interval = 1,
                output_dir = output_directory,
            )
        end

        invalid_anchor_rejections = Base.deepcopy(valid_checkpoint)
        invalid_anchor_rejections["corrector_diagnostics"][
            "window_anchor_rejections"] = -1
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_anchor_rejections, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint corrector_diagnostics window_anchor_rejections must be nonnegative",
        ) do
            Mycelia.mycelia_iterative_assemble(
                input_fastq;
                max_k = 3,
                max_iterations_per_k = 2,
                improvement_threshold = 0.0,
                stop_on_no_change = false,
                graph_mode = :singlestrand,
                verbose = false,
                enable_checkpointing = true,
                checkpoint_interval = 1,
                output_dir = output_directory,
            )
        end

        invalid_substitution_divergences = Base.deepcopy(valid_checkpoint)
        invalid_substitution_divergences["corrector_diagnostics"][
            "substitution_length_divergences"] = -1
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_substitution_divergences, 2)
        end
        checkpoint_test_error(
            ArgumentError,
            "checkpoint corrector_diagnostics " *
            "substitution_length_divergences must be nonnegative",
        ) do
            Mycelia.mycelia_iterative_assemble(
                input_fastq;
                max_k = 3,
                max_iterations_per_k = 2,
                improvement_threshold = 0.0,
                stop_on_no_change = false,
                graph_mode = :singlestrand,
                verbose = false,
                enable_checkpointing = true,
                checkpoint_interval = 1,
                output_dir = output_directory,
            )
        end

    end
end

Test.@testset "Checkpoint subdirectories cannot escape output_dir" begin
    Base.mktempdir() do temporary_directory
        input_fastq = Base.joinpath(temporary_directory, "input.fastq")
        Mycelia.write_fastq(
            records = checkpoint_resume_reads(),
            filename = input_fastq,
        )
        for subdirectory in ("checkpoints", "graphs", "progress")
            output_directory = Base.joinpath(
                temporary_directory, "run_$subdirectory")
            outside_directory = Base.joinpath(
                temporary_directory, "outside_$subdirectory")
            Base.mkpath(output_directory)
            Base.mkpath(outside_directory)
            Base.symlink(
                outside_directory,
                Base.joinpath(output_directory, subdirectory),
            )
            checkpoint_test_error(
                ArgumentError,
                "output $(subdirectory) directory must resolve inside output_dir",
            ) do
                Mycelia.mycelia_iterative_assemble(
                    input_fastq;
                    max_k = 3,
                    max_iterations_per_k = 1,
                    improvement_threshold = 0.0,
                    stop_on_no_change = false,
                    graph_mode = :singlestrand,
                    verbose = false,
                    enable_checkpointing = true,
                    checkpoint_interval = 1,
                    output_dir = output_directory,
                )
            end
        end
    end
end

Test.@testset "Iterative checkpoint advances to the next rung once" begin
    Base.mktempdir() do temporary_directory
        input_fastq = Base.joinpath(temporary_directory, "input.fastq")
        output_directory = Base.joinpath(temporary_directory, "run")
        Mycelia.write_fastq(
            records = checkpoint_resume_reads(),
            filename = input_fastq,
        )
        Mycelia.mycelia_iterative_assemble(
            input_fastq;
            max_k = 3,
            max_iterations_per_k = 1,
            improvement_threshold = 0.0,
            stop_on_no_change = false,
            graph_mode = :singlestrand,
            verbose = false,
            enable_checkpointing = true,
            checkpoint_interval = 1,
            output_dir = output_directory,
        )
        checkpoint_file = Base.joinpath(
            output_directory, "checkpoints", "latest_checkpoint.json")
        checkpoint = JSON.parsefile(checkpoint_file)
        checkpoint["next_k"] = 5
        checkpoint["next_iteration"] = 1
        checkpoint["run_complete"] = false
        checkpoint["resume_configuration"]["max_k"] = 5
        open(checkpoint_file, "w") do io
            JSON.print(io, checkpoint, 2)
        end

        result = Mycelia.mycelia_iterative_assemble(
            input_fastq;
            max_k = 5,
            max_iterations_per_k = 1,
            improvement_threshold = 0.0,
            stop_on_no_change = false,
            graph_mode = :singlestrand,
            verbose = false,
            enable_checkpointing = true,
            checkpoint_interval = 1,
            output_dir = output_directory,
        )
        metadata = result[:metadata]
        Test.@test metadata[:k_progression] == [3, 5]
        Test.@test length(metadata[:iteration_history][3]) == 1
        Test.@test length(metadata[:iteration_history][5]) == 1
        Test.@test [
            (
                Int(row[:ladder_index]),
                Int(row[:k]),
                Int(row[:iteration]),
            ) for row in metadata[:indel_rung_telemetry]
        ] == [(1, 3, 1), (2, 5, 1)]
    end
end
