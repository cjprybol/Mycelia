# Checkpoint JSON must round-trip typed history and resume at the saved NEXT-pass
# cursor without replaying a completed pass or duplicating indel telemetry.
#
# Run directly:
#   julia --project=. -e 'include("test/4_assembly/iterative_checkpoint_resume_test.jl")'

import FASTX
import JSON
import Mycelia
import Test

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
        Test.@test_throws ArgumentError Mycelia._load_hashed_input_fastq(
            input_fastq;
            after_initial_hash = () -> Base.open(input_fastq, "w") do io
                Base.write(io, replacement_bytes)
            end,
        )
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

        Test.@test_throws ArgumentError Mycelia._load_validated_checkpoint_fastq(
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
        Test.@test checkpoint["resume_configuration"]["input_fastq_path"] ==
                   Base.realpath(input_fastq)
        Test.@test checkpoint["resume_configuration"]["input_fastq_sha256"] ==
                   Mycelia._stream_sha256_file(input_fastq)
        Test.@test checkpoint["corrector_diagnostics"][
            "window_anchor_rejections"] == 0
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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        missing_profile_disabled_telemetry = Base.deepcopy(valid_checkpoint)
        delete!(
            missing_profile_disabled_telemetry,
            "indel_rung_telemetry",
        )
        open(checkpoint_file, "w") do io
            JSON.print(io, missing_profile_disabled_telemetry, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        partial_row = Base.deepcopy(valid_checkpoint)
        delete!(first(partial_row["indel_rung_telemetry"]), "requested")
        open(checkpoint_file, "w") do io
            JSON.print(io, partial_row, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        missing_full_row_field = Base.deepcopy(valid_checkpoint)
        delete!(
            first(missing_full_row_field["indel_rung_telemetry"]),
            "graph_cleanup",
        )
        open(checkpoint_file, "w") do io
            JSON.print(io, missing_full_row_field, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        invalid_root_timestamp = Base.deepcopy(valid_checkpoint)
        invalid_root_timestamp["timestamp"] = 0
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_root_timestamp, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        missing_diagnostic = Base.deepcopy(valid_checkpoint)
        delete!(missing_diagnostic["corrector_diagnostics"], "structural")
        open(checkpoint_file, "w") do io
            JSON.print(io, missing_diagnostic, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        invalid_diagnostic_type = Base.deepcopy(valid_checkpoint)
        invalid_diagnostic_type["corrector_diagnostics"]["structural"] = 0.0
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_diagnostic_type, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        mismatched_profile_request = Base.deepcopy(valid_checkpoint)
        first(mismatched_profile_request["indel_rung_telemetry"])[
            "profile_requested"] = true
        open(checkpoint_file, "w") do io
            JSON.print(io, mismatched_profile_request, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        window_counter_overflow = Base.deepcopy(valid_checkpoint)
        overflow_row = first(window_counter_overflow["indel_rung_telemetry"])
        overflow_row["requested"] = typemax(Int)
        overflow_row["admitted_windows"] = typemax(Int)
        overflow_row["rejected_windows"] = 1
        open(checkpoint_file, "w") do io
            JSON.print(io, window_counter_overflow, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        # The output directory is bound to both the canonical path and content of
        # the original input. Even byte-identical reads at another path cannot
        # inherit a stale run.
        wrong_input_fastq = Base.joinpath(
            temporary_directory, "wrong_input.fastq")
        Mycelia.write_fastq(
            records = checkpoint_resume_reads(),
            filename = wrong_input_fastq,
        )
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            wrong_input_fastq, output_directory)

        original_input_bytes = Base.read(input_fastq)
        Mycelia.write_fastq(
            records = FASTX.FASTQ.Record[
                FASTX.FASTQ.Record("changed", "ACGT", "IIII"),
            ],
            filename = input_fastq,
        )
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)
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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

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
        Test.@test_throws Exception checkpoint_resume_invocation(
            input_fastq, output_directory)
        Base.write(malformed_fastq, valid_final_fastq_bytes)

        # All serialized names are whitelisted before conversion to Symbol, so an
        # arbitrary JSON key can never be interned into the process.
        unsupported_history_field = Base.deepcopy(valid_checkpoint)
        first(unsupported_history_field["iteration_history"]["3"])[
            "attacker_controlled_symbol"] = 1
        open(checkpoint_file, "w") do io
            JSON.print(io, unsupported_history_field, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        negative_auxiliary = Base.deepcopy(valid_checkpoint)
        negative_auxiliary["cheap_correction_counts"][1] = -1
        open(checkpoint_file, "w") do io
            JSON.print(io, negative_auxiliary, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        # Size is checked from the open checkpoint descriptor before JSON parsing.
        open(checkpoint_file, "w") do io
            Base.truncate(io, Mycelia._ITERATIVE_CHECKPOINT_MAX_BYTES + 1)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)
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
        Test.@test_throws ArgumentError Mycelia.mycelia_iterative_assemble(
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

        boolean_as_integer_configuration = Base.deepcopy(valid_checkpoint)
        boolean_as_integer_configuration["resume_configuration"][
            "stop_on_no_change"] = 0
        open(checkpoint_file, "w") do io
            JSON.print(io, boolean_as_integer_configuration, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

        integer_as_float_configuration = Base.deepcopy(valid_checkpoint)
        integer_as_float_configuration["resume_configuration"]["max_k"] = 3.0
        open(checkpoint_file, "w") do io
            JSON.print(io, integer_as_float_configuration, 2)
        end
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory)

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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory; k_ladder = [3])

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
        Test.@test_throws ArgumentError checkpoint_resume_invocation(
            input_fastq, output_directory; indel_params = indel_params)

        # Existing-but-corrupt JSON is never treated as permission to restart
        # fresh and silently replay completed work.
        Base.write(checkpoint_file, "{")
        Test.@test_throws ArgumentError Mycelia.mycelia_iterative_assemble(
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
        open(checkpoint_file, "w") do io
            JSON.print(io, valid_checkpoint, 2)
        end

        invalid_fastq_hash = Base.deepcopy(valid_checkpoint)
        invalid_fastq_hash["current_fastq_sha256"] = repeat("0", 64)
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_fastq_hash, 2)
        end
        Test.@test_throws ArgumentError Mycelia.mycelia_iterative_assemble(
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

        # A corrupt explicit cursor must not replay pass 2.
        invalid_cursor = Base.deepcopy(valid_checkpoint)
        invalid_cursor["run_complete"] = false
        invalid_cursor["next_k"] = invalid_cursor["current_k"]
        invalid_cursor["next_iteration"] = invalid_cursor["current_iteration"]
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_cursor, 2)
        end
        Test.@test_throws ArgumentError Mycelia.mycelia_iterative_assemble(
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
        Test.@test_throws ArgumentError Mycelia.mycelia_iterative_assemble(
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

        invalid_metric = Base.deepcopy(valid_checkpoint)
        first(invalid_metric["indel_rung_telemetry"])["raw_frontier_metrics"] =
            Any[42]
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_metric, 2)
        end
        Test.@test_throws ArgumentError Mycelia.mycelia_iterative_assemble(
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

        invalid_anchor_rejections = Base.deepcopy(valid_checkpoint)
        invalid_anchor_rejections["corrector_diagnostics"][
            "window_anchor_rejections"] = -1
        open(checkpoint_file, "w") do io
            JSON.print(io, invalid_anchor_rejections, 2)
        end
        Test.@test_throws ArgumentError Mycelia.mycelia_iterative_assemble(
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
            Test.@test_throws ArgumentError Mycelia.mycelia_iterative_assemble(
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
