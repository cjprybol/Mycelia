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

        valid_checkpoint = JSON.parsefile(checkpoint_file)

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

        # A legacy profile-disabled checkpoint can reconstruct honest zero-valued
        # rows, but an indel-enabled resume fails closed because its prefix counters
        # cannot be recovered.
        legacy_checkpoint = Base.deepcopy(valid_checkpoint)
        for key in (
                "next_k",
                "next_iteration",
                "run_complete",
                "indel_rung_telemetry",
        )
            delete!(legacy_checkpoint, key)
        end
        open(checkpoint_file, "w") do io
            JSON.print(io, legacy_checkpoint, 2)
        end
        legacy_result = Mycelia.mycelia_iterative_assemble(
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
        legacy_telemetry = legacy_result[:metadata][:indel_rung_telemetry]
        Test.@test length(legacy_telemetry) == 2
        Test.@test all(!row[:profile_requested] for row in legacy_telemetry)
        Test.@test all(row[:requested] == 0 for row in legacy_telemetry)

        indel_params = Mycelia.IndelDecodeParams(
            0.05, 0.2, 0.2, 0.1, 0.1, 3, 3, 16)
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
            indel_params = indel_params,
        )
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
