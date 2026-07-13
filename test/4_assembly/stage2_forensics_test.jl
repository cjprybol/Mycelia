import CodecZlib
import Test

include(joinpath(
    @__DIR__,
    "..",
    "..",
    "benchmarking",
    "rhizomorph_stage2_forensics.jl",
))

Test.@testset "Stage-2 frozen seed partitions" begin
    Test.@test STAGE2_DEVELOPMENT_SEEDS == (43_003, 43_006, 43_007, 43_008, 43_009)
    Test.@test STAGE2_CONFIRMATION_SEEDS == (44_001, 44_002, 44_003)
    Test.@test isempty(
        intersect(Set(STAGE2_DEVELOPMENT_SEEDS), Set(STAGE2_CONFIRMATION_SEEDS)),
    )
end

Test.@testset "Stage-2 FASTA and FASTQ length summaries" begin
    mktempdir() do directory
        fasta_path = joinpath(directory, "records.fasta")
        fastq_path = joinpath(directory, "records.fastq.gz")
        empty_path = joinpath(directory, "empty.fasta")
        write(fasta_path, ">one\nACGT\n>two\nAACCGG\n")
        open(fastq_path, "w") do raw_stream
            gzip_stream = CodecZlib.GzipCompressorStream(raw_stream)
            write(gzip_stream, "@one\nACG\n+\nIII\n@two\nAACCG\n+\nIIIII\n")
            close(gzip_stream)
        end
        write(empty_path, "")

        fasta = summarize_fasta_record_lengths(fasta_path)
        Test.@test fasta.record_count == 2
        Test.@test fasta.total_bases == 10
        Test.@test fasta.minimum_bases == 4
        Test.@test fasta.maximum_bases == 6
        Test.@test fasta.mean_bases == 5.0

        fastq = summarize_fastq_record_lengths(fastq_path)
        Test.@test fastq.record_count == 2
        Test.@test fastq.total_bases == 8
        Test.@test fastq.minimum_bases == 3
        Test.@test fastq.maximum_bases == 5
        Test.@test fastq.mean_bases == 4.0

        empty = summarize_fasta_record_lengths(empty_path)
        Test.@test empty.record_count == 0
        Test.@test empty.total_bases == 0
        Test.@test ismissing(empty.minimum_bases)
        Test.@test ismissing(empty.maximum_bases)
        Test.@test ismissing(empty.mean_bases)
    end
end

Test.@testset "Stage-2 assembly discrepancy components" begin
    components = assembly_discrepancy_components(8_000, 8_030, 7_900, 7_920, 40)
    Test.@test components.unaligned_reference_bases == 100
    Test.@test components.unaligned_query_bases == 110
    Test.@test components.aligned_error_bases == 40
    Test.@test components.total_discrepancy_bases == 250
    Test.@test components.total_errors_per_100kbp == 3_125.0

    message = try
        assembly_discrepancy_components(100, 100, 101, 100, 0)
        ""
    catch error
        sprint(showerror, error)
    end
    Test.@test occursin("aligned_reference_bases exceeds", message)
end

Test.@testset "Stage-2 candidate assignment summary" begin
    expected = ["primary", "secondary", "tertiary"]
    exact = summarize_candidate_assignments(expected, expected)
    Test.@test exact.exact_recovery
    Test.@test isempty(exact.duplicate_truth_ids)
    Test.@test isempty(exact.missing_truth_ids)
    Test.@test isempty(exact.unexpected_truth_ids)

    collision = summarize_candidate_assignments(
        expected,
        ["primary", "primary", "unexpected"],
    )
    Test.@test !collision.exact_recovery
    Test.@test collision.unique_assigned_truth_ids == ["primary", "unexpected"]
    Test.@test collision.duplicate_truth_ids == ["primary"]
    Test.@test collision.missing_truth_ids == ["secondary", "tertiary"]
    Test.@test collision.unexpected_truth_ids == ["unexpected"]
end

function _sha256(label::AbstractString)::String
    return bytes2hex(SHA.sha256(label))
end

function _replace_observation(
        observation::Stage2ReleaseObservation;
        replacements...
)::Stage2ReleaseObservation
    values = [
        haskey(replacements, field) ? replacements[field] :
        getfield(observation, field)
        for field in fieldnames(Stage2ReleaseObservation)
    ]
    return Stage2ReleaseObservation(values...)
end

function _replace_ledger_row(
        row::Stage2AttemptLedgerRow;
        development_freeze::Union{Nothing, Stage2DevelopmentFreeze} = nothing,
        replacements...
)::Stage2AttemptLedgerRow
    values = [
        haskey(replacements, field) ? replacements[field] : getfield(row, field)
        for field in fieldnames(Stage2AttemptLedgerRow)
    ]
    return Stage2AttemptLedgerRow(values..., development_freeze)
end

function _write_confirmation_fixture!(
        directory::AbstractString,
        manifest::AbstractDict;
        bad_truth_seed::Union{Nothing, Int} = nothing
)::NamedTuple
    truth_reference_paths = Dict{Int, String}()
    truth_reference_sha256 = Dict{Int, String}()
    fixture_truth_paths = Dict{String, String}()
    fixture_truth_sha256 = Dict{String, String}()
    truth_markers = Dict(
        44_001 => "ACG",
        44_002 => "GTA",
        44_003 => "TGC",
    )
    for seed in STAGE2_CONFIRMATION_SEEDS
        relative_path = joinpath("truth", "truth-$(seed).fasta")
        fixture_relative_path = joinpath("results", relative_path)
        absolute_path = joinpath(directory, fixture_relative_path)
        mkpath(dirname(absolute_path))
        truth_ids = seed == bad_truth_seed ?
                    ("primary", "secondary", "unexpected") :
                    ("primary", "secondary", "tertiary")
        marker = truth_markers[seed]
        open(absolute_path, "w") do stream
            for (truth_id, base) in zip(truth_ids, ('A', 'C', 'G'))
                println(stream, ">$(truth_id)")
                println(stream, repeat(string(base), 7_997) * marker)
            end
        end
        digest = bytes2hex(SHA.sha256(read(absolute_path)))
        truth_reference_paths[seed] = relative_path
        truth_reference_sha256[seed] = digest
        fixture_truth_paths[string(seed)] = fixture_relative_path
        fixture_truth_sha256[string(seed)] = digest
    end

    confirmation_cells = Dict{Tuple{Symbol, Int}, NamedTuple}()
    fixture_cells = Dict{String, Any}()
    required_modalities = sort!([
        String(modality_key)
        for (modality_key, modality) in manifest["modalities"]
        if modality["required"]
    ])
    for modality_key in required_modalities
        modality = manifest["modalities"][modality_key]
        components = sort!(String.(modality["input_components"]))
        for seed in STAGE2_CONFIRMATION_SEEDS
            input_component_paths = Dict{String, String}()
            fixture_component_paths = Dict{String, String}()
            input_component_sha256 = Dict{String, String}()
            for component in components
                relative_path = joinpath(
                    "inputs", "$(modality_key)-$(seed)-$(component).txt")
                fixture_relative_path = joinpath("results", relative_path)
                absolute_path = joinpath(directory, fixture_relative_path)
                mkpath(dirname(absolute_path))
                write(absolute_path, "$(modality_key)-$(seed)-$(component)\n")
                input_component_paths[component] = relative_path
                fixture_component_paths[component] = fixture_relative_path
                input_component_sha256[component] =
                    bytes2hex(SHA.sha256(read(absolute_path)))
            end
            input_bundle_sha256 = stage2_input_component_digest(
                input_component_sha256, components)
            error_profile_path = joinpath(
                "errors", "$(modality_key)-$(seed).toml")
            fixture_error_profile_path = joinpath(
                "results", error_profile_path)
            absolute_error_profile_path = joinpath(
                directory, fixture_error_profile_path)
            mkpath(dirname(absolute_error_profile_path))
            open(absolute_error_profile_path, "w") do stream
                TOML.print(stream, Dict(
                    "schema" => "rhizomorph-stage2-error-profile-v1",
                    "modality" => modality_key,
                    "seed" => seed,
                ))
            end
            error_profile_sha256 = bytes2hex(
                SHA.sha256(read(absolute_error_profile_path)))
            modality_symbol = Symbol(modality_key)
            confirmation_cells[(modality_symbol, seed)] = (;
                input_component_paths,
                input_component_sha256,
                input_bundle_sha256,
                error_profile_path,
                error_profile_sha256,
            )
            fixture_cells["$(modality_key):$(seed)"] = Dict(
                "input_component_paths" => fixture_component_paths,
                "input_component_sha256" => input_component_sha256,
                "input_bundle_sha256" => input_bundle_sha256,
                "error_profile_path" => fixture_error_profile_path,
                "error_profile_sha256" => error_profile_sha256,
            )
        end
    end

    fixture_path = joinpath(directory, "fixture-index.toml")
    open(fixture_path, "w") do stream
        TOML.print(stream, Dict(
            "schema" => "rhizomorph-stage2-fixture-index-v1",
            "confirmation_truth_path" => fixture_truth_paths,
            "confirmation_truth_sha256" => fixture_truth_sha256,
            "confirmation_cells" => fixture_cells,
        ))
    end
    return (;
        fixture_path,
        fixture_sha256 = bytes2hex(SHA.sha256(read(fixture_path))),
        truth_reference_paths,
        truth_reference_sha256,
        confirmation_cells,
    )
end

function _write_strict_development_gate!(
        directory::AbstractString,
        result_path::AbstractString,
        manifest::AbstractDict;
        git_sha::AbstractString,
        manifest_sha256::AbstractString,
        configuration_sha256::AbstractString,
        fixture_sha256::AbstractString,
        grid_mutation::Symbol = :none,
        stale_source_provenance::Bool = false
)::String
    absolute_result_path = joinpath(directory, result_path)
    mkpath(dirname(absolute_result_path))
    decisions = Dict{String, Any}[]
    required_modalities = sort!([
        String(modality_key)
        for (modality_key, modality) in manifest["modalities"]
        if modality["required"]
    ])
    for modality_key in required_modalities
        for seed in sort!(Int.(manifest["development_seeds"]))
            cell_id = "$(modality_key):$(seed):rhizomorph-native"
            source_name = "dev-source-$(modality_key)-$(seed).toml"
            source_path = joinpath(dirname(absolute_result_path), source_name)
            source_manifest_sha256 =
                stale_source_provenance && isempty(decisions) ?
                repeat("0", 64) : manifest_sha256
            open(source_path, "w") do stream
                TOML.print(stream, Dict(
                    "schema" => "rhizomorph-stage2-development-cell-v1",
                    "cell_id" => cell_id,
                    "passed" => true,
                    "git_sha" => git_sha,
                    "manifest_sha256" => source_manifest_sha256,
                    "configuration_sha256" => configuration_sha256,
                    "fixture_sha256" => fixture_sha256,
                ))
            end
            push!(decisions, Dict(
                "cell_id" => cell_id,
                "passed" => true,
                "source_report_path" => source_name,
                "source_report_sha256" =>
                    bytes2hex(SHA.sha256(read(source_path))),
            ))
        end
    end
    if grid_mutation == :missing
        pop!(decisions)
    elseif grid_mutation == :duplicate
        push!(decisions, copy(first(decisions)))
    elseif grid_mutation == :extra
        cell_id = "unexpected:43003:rhizomorph-native"
        source_name = "dev-source-unexpected-43003.toml"
        source_path = joinpath(dirname(absolute_result_path), source_name)
        open(source_path, "w") do stream
            TOML.print(stream, Dict(
                "schema" => "rhizomorph-stage2-development-cell-v1",
                "cell_id" => cell_id,
                "passed" => true,
                "git_sha" => git_sha,
                "manifest_sha256" => manifest_sha256,
                "configuration_sha256" => configuration_sha256,
                "fixture_sha256" => fixture_sha256,
            ))
        end
        push!(decisions, Dict(
            "cell_id" => cell_id,
            "passed" => true,
            "source_report_path" => source_name,
            "source_report_sha256" =>
                bytes2hex(SHA.sha256(read(source_path))),
        ))
    elseif grid_mutation != :none
        error("unsupported development-grid mutation: $(grid_mutation)")
    end
    open(absolute_result_path, "w") do stream
        TOML.print(stream, Dict(
            "schema" => "rhizomorph-stage2-development-gate-v1",
            "gate_id" => "strict-development-grid",
            "strict_gate_passed" => true,
            "git_sha" => git_sha,
            "manifest_sha256" => manifest_sha256,
            "configuration_sha256" => configuration_sha256,
            "fixture_sha256" => fixture_sha256,
            "cell_decisions" => decisions,
        ))
    end
    return bytes2hex(SHA.sha256(read(absolute_result_path)))
end

function _build_release_case(
        directory::AbstractString;
        resolve_comparators::Bool = true,
        failing_cell::Union{Nothing, Tuple{Symbol, Int}} = nothing,
        malformed_cell::Union{Nothing, Tuple{Symbol, Int}} = nothing,
        conflicting_native::Bool = false,
        component_mismatch::Bool = false,
        full_baseline_failure::Bool = false,
        failed_attempt::Bool = false,
        forged_assignment::Bool = false,
        failed_subset_attempts::Bool = false,
        incomplete_attempt::Bool = false,
        metric_report_mismatch::Bool = false,
        metric_nonce_mismatch::Bool = false,
        populate_confirmation::Bool = true,
        bad_truth_seed::Union{Nothing, Int} = nothing,
        development_grid_mutation::Symbol = :none,
        stale_development_source_provenance::Bool = false
)::NamedTuple
    source_manifest = joinpath(
        @__DIR__, "..", "..", "benchmarking",
        "rhizomorph_stage2_multimodal_benchmark.toml")
    manifest = TOML.parsefile(source_manifest)
    if resolve_comparators
        for (comparator_id, comparator) in manifest["comparators"]
            comparator["version"] = "locked-$(comparator_id)"
            comparator["command"] = "run-$(comparator_id) --frozen-input"
        end
        manifest["evaluation"]["version"] = "QUAST-5.3.0"
    end
    manifest_path = joinpath(directory, "benchmark.toml")
    open(manifest_path, "w") do stream
        TOML.print(stream, manifest)
    end
    manifest_sha256 = bytes2hex(SHA.sha256(read(manifest_path)))
    configuration_path = joinpath(directory, "configuration.toml")
    open(configuration_path, "w") do stream
        TOML.print(stream, Dict(
            "schema" => "rhizomorph-stage2-configuration-v1",
            "configuration_id" => "frozen-configuration",
        ))
    end
    configuration_sha256 = bytes2hex(SHA.sha256(read(configuration_path)))
    fixture = _write_confirmation_fixture!(
        directory, manifest; bad_truth_seed)
    fixture_path = fixture.fixture_path
    fixture_sha256 = fixture.fixture_sha256
    truth_reference_paths = fixture.truth_reference_paths
    truth_reference_sha256 = fixture.truth_reference_sha256
    confirmation_cells = fixture.confirmation_cells
    git_sha = repeat("d", 40)
    development_ledger_path = joinpath(directory, "development.tsv")
    development_result_path = joinpath("results", "dev-final.toml")
    development_result_sha256 = _write_strict_development_gate!(
        directory,
        development_result_path,
        manifest;
        git_sha,
        manifest_sha256,
        configuration_sha256,
        fixture_sha256,
        grid_mutation = development_grid_mutation,
        stale_source_provenance = stale_development_source_provenance,
    )
    development_row = Stage2AttemptLedgerRow(;
        attempt_id = "dev-final",
        recorded_at_utc = "2026-07-12T16:00:00Z",
        seed_partition = :development,
        seed = 43_003,
        modality = :long_noisy,
        configuration_id = "frozen-configuration",
        git_sha,
        git_dirty = false,
        manifest_sha256,
        configuration_sha256,
        fixture_sha256,
        input_bundle_sha256 = _sha256("development-bundle"),
        consumed_input_sha256 = _sha256("development-bundle"),
        error_profile_sha256 = _sha256("development-errors"),
        comparator_id = "rhizomorph-native",
        comparator_version = git_sha,
        command_sha256 = _sha256("native-command"),
        status = :passed,
        result_path = development_result_path,
        result_sha256 = development_result_sha256,
        notes = "final development attempt",
    )
    append_stage2_attempt!(development_ledger_path, development_row)
    seal = seal_stage2_development_ledger!(
        development_ledger_path;
        frozen_at_utc = "2026-07-12T16:30:00Z",
        manifest_path,
        configuration_path,
        fixture_path,
        git_sha,
        manifest_sha256,
        configuration_sha256,
        fixture_sha256,
        git_dirty = false,
    )
    confirmation_ledger_path = stage2_confirmation_ledger_path(
        development_ledger_path)
    result_paths = String[]
    assignment = summarize_candidate_assignments(
        ["primary", "secondary", "tertiary"],
        ["primary", "secondary", "tertiary"],
    )
    row_index = 0
    confirmation_reservations = Stage2AttemptLedgerRow[]
    for (modality_key, modality) in manifest["modalities"]
        modality["required"] || continue
        modality_symbol = Symbol(modality_key)
        modality_components = Set(String.(modality["input_components"]))
        for seed in STAGE2_CONFIRMATION_SEEDS
            cell = (modality_symbol, seed)
            failing = cell == failing_cell
            malformed = cell == malformed_cell
            confirmation_cell = confirmation_cells[cell]
            input_component_paths = confirmation_cell.input_component_paths
            input_component_sha256 = confirmation_cell.input_component_sha256
            input_bundle_sha256 = confirmation_cell.input_bundle_sha256
            native_result_path = joinpath(
                "raw", "native-$(modality_key)-$(seed).fasta")
            truth_reference_path = truth_reference_paths[seed]
            cell_truth_reference_sha256 = truth_reference_sha256[seed]
            for comparator_id in modality["comparators"]
                row_index += 1
                comparator = manifest["comparators"][comparator_id]
                components = String.(
                    comparator["consumed_input_components"][modality_key])
                is_full = Set(components) == modality_components
                consumed_sha256 = stage2_input_component_digest(
                    input_component_sha256, components)
                if component_mismatch && row_index == 1
                    components = ["not-in-bundle"]
                end
                native_error = failing ? 1_200.0 :
                               (seed == 44_003 ? 950.0 : 800.0)
                comparator_error = full_baseline_failure &&
                                   modality_symbol in (
                    :hybrid_ont_short, :hybrid_hifi_ont) &&
                                   is_full ? 700.0 : 1_000.0
                cell_assignment = failing ? summarize_candidate_assignments(
                    ["primary", "secondary", "tertiary"],
                    ["primary", "primary", "tertiary"],
                ) : assignment
                if forged_assignment && row_index == 1
                    cell_assignment = CandidateAssignmentSummary(
                        ["primary", "secondary", "tertiary"],
                        ["primary", "primary", "tertiary"],
                        ["primary", "tertiary"],
                        String[],
                        String[],
                        String[],
                        true,
                    )
                end
                attempt_id = "confirmation-$(row_index)"
                comparator_result_path = joinpath(
                    "raw", "$(attempt_id)-comparator.fasta")
                native_nga50 = malformed ? -1.0 :
                               (conflicting_native && row_index == 2 ?
                                7_949.0 : 7_950.0)
                truth_coverages = failing ? [0.98, 0.996, 0.994] :
                                  [0.995, 0.996, 0.994]
                abundance_mae = failing ? 0.06 : 0.03
                strain_count_mode = failing ? 2 : 3
                evaluator = manifest["evaluation"]
                evaluator_command_sha256 = _sha256(evaluator["command"])
                metric_report_path = joinpath(
                    "metrics", "$(attempt_id).toml")
                result_path = joinpath("results", "$(attempt_id).toml")
                terminal_status = failed_attempt && row_index == 1 ? :failed :
                                  (failed_subset_attempts && !is_full ?
                                   :failed : :passed)
                metric_mismatch = metric_report_mismatch && row_index == 1
                metric_nonce_mismatch_for_row =
                    metric_nonce_mismatch && row_index == 1
                reservation = Stage2AttemptLedgerRow(;
                    attempt_id,
                    recorded_at_utc = "2026-07-12T16:59:00Z",
                    seed_partition = :confirmation,
                    seed,
                    modality = modality_symbol,
                    configuration_id = "frozen-configuration",
                    git_sha,
                    git_dirty = false,
                    manifest_sha256,
                    configuration_sha256,
                    fixture_sha256,
                    truth_reference_sha256 = cell_truth_reference_sha256,
                    input_bundle_sha256,
                    consumed_input_sha256 = consumed_sha256,
                    error_profile_sha256 =
                        confirmation_cell.error_profile_sha256,
                    comparator_id,
                    comparator_version = comparator["version"],
                    command_sha256 = _sha256(comparator["command"]),
                    comparator_configuration_sha256 =
                        _sha256(comparator["configuration"]),
                    status = :reserved,
                    result_path = "",
                    result_sha256 = "",
                    development_freeze_sha256 = seal.sha256,
                    notes = "reserved before holdout execution",
                    development_freeze = seal.freeze,
                )
                push!(confirmation_reservations, reservation)
                if populate_confirmation
                    if incomplete_attempt && row_index == 1
                        append_stage2_confirmation_attempt!(
                            reservation;
                            development_ledger_path,
                        )
                        continue
                    end
                    run_stage2_confirmation_attempt!(
                        reservation,
                        (_, reservation_nonce_sha256) -> begin
                            absolute_native_result = normpath(joinpath(
                                directory, "results", native_result_path))
                            mkpath(dirname(absolute_native_result))
                            write(
                                absolute_native_result,
                                ">native\n$(repeat("A", 8_000))\n",
                            )
                            native_result_sha256 = bytes2hex(
                                SHA.sha256(read(absolute_native_result)))
                            absolute_comparator_result = normpath(joinpath(
                                directory, "results", comparator_result_path))
                            write(
                                absolute_comparator_result,
                                ">comparator\n$(repeat("A", 8_000))\n",
                            )
                            comparator_result_sha256 = bytes2hex(
                                SHA.sha256(read(absolute_comparator_result)))
                            absolute_metric_report = normpath(joinpath(
                                directory, "results", metric_report_path))
                            mkpath(dirname(absolute_metric_report))
                            open(absolute_metric_report, "w") do stream
                                TOML.print(stream, Dict{String, Any}(
                                    "schema" => evaluator["report_schema"],
                                    "evaluator_id" => evaluator["tool"],
                                    "evaluator_version" => evaluator["version"],
                                    "evaluator_command_sha256" =>
                                        evaluator_command_sha256,
                                    "reservation_nonce_sha256" =>
                                        metric_nonce_mismatch_for_row ?
                                        repeat("0", 64) :
                                        reservation_nonce_sha256,
                                    "truth_reference_sha256" =>
                                        cell_truth_reference_sha256,
                                    "native_result_sha256" =>
                                        native_result_sha256,
                                    "comparator_result_sha256" =>
                                        comparator_result_sha256,
                                    "native_nga50" => native_nga50,
                                    "comparator_nga50" => 8_000.0,
                                    "native_genome_fraction" => 99.7,
                                    "comparator_genome_fraction" => 100.0,
                                    "native_misassemblies" => 0,
                                    "comparator_misassemblies" => 0,
                                    "native_error_per_100kbp" => native_error,
                                    "comparator_error_per_100kbp" =>
                                        comparator_error,
                                    "truth_reference_coverages" =>
                                        truth_coverages,
                                    "abundance_mae" => abundance_mae,
                                    "comparator_abundance_mae" => 0.04,
                                    "strain_count_mode" => strain_count_mode,
                                    "strain_count_set_contains_truth" => !failing,
                                    "primary_rank_agreement" => !failing,
                                    "ensemble_noninferior_or_abstained" => !failing,
                                    "assignment" => Dict{String, Any}(
                                        "expected_truth_ids" =>
                                            cell_assignment.expected_truth_ids,
                                        "assigned_truth_ids" =>
                                            cell_assignment.assigned_truth_ids,
                                        "unique_assigned_truth_ids" =>
                                            cell_assignment.unique_assigned_truth_ids,
                                        "duplicate_truth_ids" =>
                                            cell_assignment.duplicate_truth_ids,
                                        "missing_truth_ids" =>
                                            cell_assignment.missing_truth_ids,
                                        "unexpected_truth_ids" =>
                                            cell_assignment.unexpected_truth_ids,
                                        "exact_recovery" =>
                                            cell_assignment.exact_recovery,
                                    ),
                                ))
                            end
                            metric_report_sha256 = bytes2hex(
                                SHA.sha256(read(absolute_metric_report)))
                            observation = Stage2ReleaseObservation(
                                modality_symbol,
                                seed,
                                attempt_id,
                                reservation_nonce_sha256,
                                native_result_path,
                                native_result_sha256,
                                metric_mismatch ?
                                native_nga50 + 1.0 : native_nga50,
                                8_000.0,
                                99.7,
                                100.0,
                                0,
                                0,
                                native_error,
                                comparator_error,
                                comparator_id,
                                comparator["version"],
                                comparator_result_path,
                                comparator_result_sha256,
                                _sha256(comparator["command"]),
                                _sha256(comparator["configuration"]),
                                evaluator["tool"],
                                evaluator["version"],
                                evaluator_command_sha256,
                                seal.sha256,
                                manifest_sha256,
                                configuration_sha256,
                                fixture_sha256,
                                input_component_paths,
                                input_component_sha256,
                                truth_reference_path,
                                cell_truth_reference_sha256,
                                metric_report_path,
                                metric_report_sha256,
                                input_bundle_sha256,
                                input_bundle_sha256,
                                consumed_sha256,
                                components,
                                cell_assignment,
                                truth_coverages,
                                abundance_mae,
                                0.04,
                                strain_count_mode,
                                !failing,
                                !failing,
                                !failing,
                            )
                            absolute_result_path = joinpath(directory, result_path)
                            write_stage2_release_observation(
                                absolute_result_path, observation)
                            push!(result_paths, absolute_result_path)
                            return _replace_ledger_row(
                                reservation;
                                development_freeze = seal.freeze,
                                recorded_at_utc = "2026-07-12T17:00:00Z",
                                status = terminal_status,
                                result_path,
                                result_sha256 = bytes2hex(
                                    SHA.sha256(read(absolute_result_path))),
                                notes = "frozen confirmation attempt",
                            )
                        end;
                        development_ledger_path,
                    )
                end
            end
        end
    end
    return (;
        manifest_path,
        development_ledger_path,
        confirmation_ledger_path,
        seal,
        result_paths,
        truth_reference_sha256,
        confirmation_reservations,
    )
end

Test.@testset "Stage-2 ledger-owned release gate" begin
    mktempdir() do directory
        case = _build_release_case(directory)
        passing = evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
        Test.@test passing.passed
        Test.@test passing.development_freeze_sha256 == case.seal.sha256
        Test.@test passing.confirmation_ledger_sha256 ==
                   bytes2hex(SHA.sha256(read(case.confirmation_ledger_path)))
        Test.@test length(passing.attempt_result_sha256) == 48
        Test.@test length(passing.diagnostic_attempt_status) == 9
        Test.@test length(passing.diagnostic_result_sha256) == 9
        Test.@test all(==(:passed), values(passing.diagnostic_attempt_status))
        Test.@test all(gates["accuracy_wins"]
            for gates in values(passing.modality_gates))
        Test.@test length(unique(values(case.truth_reference_sha256))) ==
                   length(STAGE2_CONFIRMATION_SEEDS)
        observation = _parse_stage2_release_observation(
            read(first(case.result_paths)))
        different_assignment = summarize_candidate_assignments(
            ["primary", "secondary", "tertiary"],
            ["secondary", "primary", "tertiary"],
        )
        different_component_paths = copy(observation.input_component_paths)
        first_component = first(keys(different_component_paths))
        different_component_paths[first_component] = "different-input.fastq"
        different_component_sha256 = copy(observation.input_component_sha256)
        different_component_sha256[first_component] = repeat("e", 64)
        mutations = (
            _replace_observation(
                observation; native_result_path = "different-native.fasta"),
            _replace_observation(observation; native_result_sha256 = repeat("e", 64)),
            _replace_observation(observation; native_nga50 = 7_949.0),
            _replace_observation(observation; native_genome_fraction = 99.6),
            _replace_observation(observation; native_misassemblies = 1),
            _replace_observation(observation; native_error_per_100kbp = 801.0),
            _replace_observation(observation; assignment = different_assignment),
            _replace_observation(
                observation; truth_reference_coverages = [0.994, 0.996, 0.994]),
            _replace_observation(observation; abundance_mae = 0.031),
            _replace_observation(observation; strain_count_mode = 2),
            _replace_observation(
                observation; strain_count_set_contains_truth = false),
            _replace_observation(observation; primary_rank_agreement = false),
            _replace_observation(
                observation; ensemble_noninferior_or_abstained = false),
            _replace_observation(
                observation; input_bundle_sha256 = repeat("f", 64)),
            _replace_observation(
                observation; input_component_paths = different_component_paths),
            _replace_observation(
                observation; input_component_sha256 = different_component_sha256),
            _replace_observation(
                observation; truth_reference_path = "different-truth.fasta"),
            _replace_observation(
                observation; truth_reference_sha256 = repeat("f", 64)),
            _replace_observation(
                observation; native_consumed_input_sha256 = repeat("f", 64)),
        )
        Test.@test all(
            !isequal(
                _native_scientific_signature(observation),
                _native_scientific_signature(mutation),
            )
            for mutation in mutations
        )
    end
    mktempdir() do directory
        case = _build_release_case(
            directory; failing_cell = (:long_noisy, 44_001))
        failing = evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
        Test.@test !failing.passed
        Test.@test !failing.observation_gates[
            (:long_noisy, 44_001)]["three_truths"]
    end
    for options in (
            (; malformed_cell = (:long_noisy, 44_001)),
            (; conflicting_native = true),
            (; component_mismatch = true),
            (; failed_attempt = true),
            (; forged_assignment = true),
            (; metric_report_mismatch = true),
            (; metric_nonce_mismatch = true),
        )
        mktempdir() do directory
            case = _build_release_case(directory; options...)
            confirmation_before = read(case.confirmation_ledger_path)
            Test.@test_throws ErrorException evaluate_stage2_release_gate(;
                case.manifest_path,
                case.development_ledger_path,
            )
            Test.@test read(case.confirmation_ledger_path) ==
                       confirmation_before
        end
    end
    mktempdir() do directory
        case = _build_release_case(directory; full_baseline_failure = true)
        result = evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
        Test.@test !result.passed
        Test.@test !result.observation_gates[
            (:hybrid_ont_short, 44_001)]["no_large_error_regression"]
        Test.@test !result.observation_gates[
            (:hybrid_hifi_ont, 44_001)]["no_large_error_regression"]
    end
    mktempdir() do directory
        case = _build_release_case(directory; failed_subset_attempts = true)
        result = evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
        Test.@test result.passed
        Test.@test all(==(:failed), values(result.diagnostic_attempt_status))
        Test.@test length(result.diagnostic_result_sha256) == 9
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        manifest = TOML.parsefile(case.manifest_path)
        subset_result = first([
            path for path in case.result_paths if begin
                observation = _parse_stage2_release_observation(read(path))
                Set(observation.comparator_consumed_input_components) != Set(
                    manifest["modalities"][String(observation.modality)][
                        "input_components"])
            end
        ])
        open(subset_result, "a") do stream
            write(stream, "\n# tampered subset report\n")
        end
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory; incomplete_attempt = true)
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        open(first(case.result_paths), "a") do stream
            write(stream, "\n# tampered\n")
        end
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        observation_path = first(case.result_paths)
        Test.@test_throws ErrorException _observation_relative_path(
            observation_path, "../outside")
        Test.@test_throws ErrorException _observation_relative_path(
            observation_path, abspath("outside"))
        observation = _parse_stage2_release_observation(read(observation_path))
        native_result_path = _observation_relative_path(
            observation_path, observation.native_result_path)
        open(native_result_path, "a") do stream
            write(stream, "N")
        end
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        observation_path = first(case.result_paths)
        observation = _parse_stage2_release_observation(read(observation_path))
        comparator_result_path = _observation_relative_path(
            observation_path, observation.comparator_result_path)
        open(comparator_result_path, "a") do stream
            write(stream, "N")
        end
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        target = joinpath(directory, "target.txt")
        link = joinpath(directory, "link.txt")
        write(target, "target")
        symlink(target, link)
        Test.@test_throws ErrorException _validate_bundle_file(
            link, directory, "symlink test")
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        observation_path = first(case.result_paths)
        observation = _parse_stage2_release_observation(read(observation_path))
        component = first(keys(observation.input_component_paths))
        component_path = _observation_relative_path(
            observation_path, observation.input_component_paths[component])
        open(component_path, "a") do stream
            write(stream, "tampered")
        end
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        rows = load_stage2_attempt_ledger(
            case.confirmation_ledger_path;
            development_freeze = case.seal.freeze,
        )
        source_row = first(rows)
        source_observation = _parse_stage2_release_observation(
            read(first(case.result_paths)))
        duplicate_attempt_id = "confirmation-repeat"
        duplicate_observation = _replace_observation(
            source_observation; attempt_id = duplicate_attempt_id)
        duplicate_result_path = joinpath("results", "repeat.toml")
        absolute_result_path = joinpath(directory, duplicate_result_path)
        write_stage2_release_observation(
            absolute_result_path, duplicate_observation)
        duplicate_row = Stage2AttemptLedgerRow(;
            attempt_id = duplicate_attempt_id,
            recorded_at_utc = "2026-07-12T17:30:00Z",
            seed_partition = :confirmation,
            seed = source_row.seed,
            modality = source_row.modality,
            configuration_id = source_row.configuration_id,
            git_sha = source_row.git_sha,
            git_dirty = false,
            manifest_sha256 = source_row.manifest_sha256,
            configuration_sha256 = source_row.configuration_sha256,
            fixture_sha256 = source_row.fixture_sha256,
            truth_reference_sha256 = source_row.truth_reference_sha256,
            input_bundle_sha256 = source_row.input_bundle_sha256,
            consumed_input_sha256 = source_row.consumed_input_sha256,
            error_profile_sha256 = source_row.error_profile_sha256,
            comparator_id = source_row.comparator_id,
            comparator_version = source_row.comparator_version,
            command_sha256 = source_row.command_sha256,
            comparator_configuration_sha256 =
                source_row.comparator_configuration_sha256,
            status = :passed,
            result_path = duplicate_result_path,
            result_sha256 = bytes2hex(SHA.sha256(read(absolute_result_path))),
            development_freeze_sha256 = case.seal.sha256,
            notes = "forbidden holdout rerun",
            development_freeze = case.seal.freeze,
        )
        duplicate_reservation = _replace_ledger_row(
            duplicate_row;
            development_freeze = case.seal.freeze,
            status = :reserved,
            result_path = "",
            result_sha256 = "",
            notes = "forbidden holdout rerun reservation",
        )
        retained_lines = filter(readlines(case.confirmation_ledger_path)) do line
            !startswith(line, source_row.attempt_id * "\t")
        end
        write(case.confirmation_ledger_path, join(retained_lines, '\n') * "\n")
        Test.@test_throws ErrorException append_stage2_confirmation_attempt!(
            duplicate_reservation;
            development_ledger_path = case.development_ledger_path,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        lines = readlines(case.confirmation_ledger_path; keep = true)
        write(case.confirmation_ledger_path, join(lines[1:(end - 1)]))
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        open(case.development_ledger_path, "a") do stream
            write(stream, "tampered-after-freeze\n")
        end
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        open(case.seal.path, "a") do stream
            write(stream, "\n# freeze tamper\n")
        end
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        original_directory = joinpath(directory, "original")
        clone_directory = joinpath(directory, "clone")
        mkpath(original_directory)
        mkpath(clone_directory)
        case = _build_release_case(original_directory)
        cloned_ledger = joinpath(clone_directory, "development.tsv")
        cp(case.development_ledger_path, cloned_ledger)
        cp(case.seal.path, stage2_development_freeze_path(cloned_ledger))
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            development_ledger_path = cloned_ledger,
        )
    end
    mktempdir() do directory
        case = _build_release_case(directory)
        open(case.manifest_path, "a") do stream
            write(stream, "\n# manifest tamper\n")
        end
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
    mktempdir() do directory
        message = try
            _build_release_case(directory; resolve_comparators = false)
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin("resolved evaluator version", message)
    end
end

Test.@testset "Stage-2 strict development grid identity" begin
    for mutation in (:missing, :duplicate, :extra)
        mktempdir() do directory
            Test.@test_throws ErrorException _build_release_case(
                directory;
                populate_confirmation = false,
                development_grid_mutation = mutation,
            )
            Test.@test !isfile(stage2_development_freeze_path(
                joinpath(directory, "development.tsv")))
        end
    end
    mktempdir() do directory
        Test.@test_throws ErrorException _build_release_case(
            directory;
            populate_confirmation = false,
            stale_development_source_provenance = true,
        )
        Test.@test !isfile(stage2_development_freeze_path(
            joinpath(directory, "development.tsv")))
    end
end

Test.@testset "Stage-2 confirmation preflight is mutation-free" begin
    mktempdir() do directory
        case = _build_release_case(directory; populate_confirmation = false)
        valid_reservation = only(filter(case.confirmation_reservations) do row
            row.modality == :long_noisy && row.seed == 44_001 &&
                row.comparator_id == "metaflye-strainy"
        end)
        target_reservation = only(filter(case.confirmation_reservations) do row
            row.modality == :long_noisy && row.seed == 44_002 &&
                row.comparator_id == "metaflye-strainy"
        end)
        subset_target_reservation = only(filter(
            case.confirmation_reservations) do row
            row.modality == :hybrid_ont_short && row.seed == 44_002 &&
                row.comparator_id == "metaflye-strainy"
        end)
        append_stage2_confirmation_attempt!(
            valid_reservation;
            development_ledger_path = case.development_ledger_path,
        )
        ledger_before = read(case.confirmation_ledger_path)
        invalid_reservations = (
            _replace_ledger_row(
                target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-cross-seed-truth",
                truth_reference_sha256 =
                    case.truth_reference_sha256[44_001],
            ),
            _replace_ledger_row(
                target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-configuration-id",
                configuration_id = "wrong-configuration",
            ),
            _replace_ledger_row(
                target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-comparator",
                comparator_id = "unregistered-comparator",
            ),
            _replace_ledger_row(
                target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-comparator-version",
                comparator_version = "wrong-version",
            ),
            _replace_ledger_row(
                target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-command",
                command_sha256 = repeat("a", 64),
            ),
            _replace_ledger_row(
                target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-comparator-configuration",
                comparator_configuration_sha256 = repeat("b", 64),
            ),
            _replace_ledger_row(
                target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-input-bundle",
                input_bundle_sha256 = repeat("c", 64),
            ),
            _replace_ledger_row(
                subset_target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-consumed-input",
                consumed_input_sha256 = repeat("d", 64),
            ),
            _replace_ledger_row(
                target_reservation;
                development_freeze = case.seal.freeze,
                attempt_id = "wrong-error-profile",
                error_profile_sha256 = repeat("e", 64),
            ),
        )
        for invalid_reservation in invalid_reservations
            receipt_directory = _confirmation_reservation_directory(
                case.confirmation_ledger_path, invalid_reservation)
            Test.@test !ispath(receipt_directory)
            Test.@test_throws ErrorException append_stage2_confirmation_attempt!(
                invalid_reservation;
                development_ledger_path = case.development_ledger_path,
            )
            Test.@test read(case.confirmation_ledger_path) == ledger_before
            Test.@test !ispath(receipt_directory)
        end
    end
end

Test.@testset "Stage-2 evaluation rejects wrong truth FASTA IDs" begin
    mktempdir() do directory
        case = _build_release_case(directory; bad_truth_seed = 44_001)
        Test.@test_throws ErrorException evaluate_stage2_release_gate(;
            case.manifest_path,
            case.development_ledger_path,
        )
    end
end

Test.@testset "Stage-2 append-only attempt ledger" begin
    mktempdir() do directory
        source_manifest = joinpath(
            @__DIR__, "..", "..", "benchmarking",
            "rhizomorph_stage2_multimodal_benchmark.toml")
        manifest = TOML.parsefile(source_manifest)
        for (comparator_id, comparator) in manifest["comparators"]
            comparator["version"] = "locked-$(comparator_id)"
            comparator["command"] = "run-$(comparator_id) --frozen-input"
        end
        manifest["evaluation"]["version"] = "QUAST-5.3.0"
        manifest_path = joinpath(directory, "benchmark.toml")
        open(manifest_path, "w") do stream
            TOML.print(stream, manifest)
        end
        configuration_path = joinpath(directory, "configuration.toml")
        open(configuration_path, "w") do stream
            TOML.print(stream, Dict(
                "schema" => "rhizomorph-stage2-configuration-v1",
                "configuration_id" => "frozen-configuration",
            ))
        end
        fixture = _write_confirmation_fixture!(directory, manifest)
        fixture_path = fixture.fixture_path
        confirmation_truth_sha256 = fixture.truth_reference_sha256
        confirmation_cells = fixture.confirmation_cells
        cell_44_001 = confirmation_cells[(:long_noisy, 44_001)]
        cell_44_002 = confirmation_cells[(:long_noisy, 44_002)]
        cell_44_003 = confirmation_cells[(:long_noisy, 44_003)]
        sha_a = bytes2hex(SHA.sha256(read(manifest_path)))
        sha_b = bytes2hex(SHA.sha256(read(configuration_path)))
        sha_c = fixture.fixture_sha256
        git_sha = repeat("d", 40)
        ledger_path = joinpath(directory, "attempts.tsv")
        gate_result_path = joinpath("results", "development-gate.toml")
        gate_result_sha256 = _write_strict_development_gate!(
            directory,
            gate_result_path,
            manifest;
            git_sha,
            manifest_sha256 = sha_a,
            configuration_sha256 = sha_b,
            fixture_sha256 = sha_c,
        )
        first_row = Stage2AttemptLedgerRow(;
            attempt_id = "dev-001",
            recorded_at_utc = "2026-07-12T16:00:00Z",
            seed_partition = :development,
            seed = 43_003,
            modality = :long_noisy,
            configuration_id = "configuration-a",
            git_sha,
            git_dirty = true,
            manifest_sha256 = sha_a,
            configuration_sha256 = sha_b,
            fixture_sha256 = sha_c,
            input_bundle_sha256 = sha_a,
            consumed_input_sha256 = sha_a,
            error_profile_sha256 = sha_b,
            comparator_id = "rhizomorph-native",
            comparator_version = "dev",
            command_sha256 = sha_c,
            status = :failed,
            result_path = "results/dev-001",
            result_sha256 = sha_a,
            notes = "preserve failed development attempt",
        )
        second_row = Stage2AttemptLedgerRow(;
            attempt_id = "dev-002",
            recorded_at_utc = "2026-07-12T16:05:00Z",
            seed_partition = :development,
            seed = 43_006,
            modality = :long_noisy,
            configuration_id = "frozen-configuration",
            git_sha,
            git_dirty = false,
            manifest_sha256 = sha_a,
            configuration_sha256 = sha_b,
            fixture_sha256 = sha_c,
            input_bundle_sha256 = sha_b,
            consumed_input_sha256 = sha_b,
            error_profile_sha256 = sha_b,
            comparator_id = "rhizomorph-native",
            comparator_version = git_sha,
            command_sha256 = sha_c,
            status = :passed,
            result_path = gate_result_path,
            result_sha256 = gate_result_sha256,
            notes = "final strict-gate pass",
        )
        append_stage2_attempt!(ledger_path, first_row)
        original = read(ledger_path, String)
        append_stage2_attempt!(ledger_path, second_row)
        appended = read(ledger_path, String)
        Test.@test startswith(appended, original)
        Test.@test countlines(ledger_path) == 3

        duplicate_message = try
            append_stage2_attempt!(ledger_path, first_row)
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin("duplicate attempt_id: dev-001", duplicate_message)
        Test.@test read(ledger_path, String) == appended

        seed_message = try
            Stage2AttemptLedgerRow(;
                attempt_id = "bad-seed",
                recorded_at_utc = "2026-07-12T16:10:00Z",
                seed_partition = :development,
                seed = 44_001,
                modality = :long_noisy,
                configuration_id = "configuration-c",
                git_sha,
                git_dirty = true,
                manifest_sha256 = sha_a,
                configuration_sha256 = sha_b,
                fixture_sha256 = sha_c,
                input_bundle_sha256 = sha_b,
                consumed_input_sha256 = sha_b,
                error_profile_sha256 = sha_b,
                comparator_id = "rhizomorph-native",
                comparator_version = "dev",
                command_sha256 = sha_c,
                status = :passed,
                result_path = "results/bad-seed",
                result_sha256 = sha_b,
                notes = "confirmation seed must not enter development",
            )
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin("not frozen in the development partition", seed_message)

        missing_freeze_message = try
            Stage2AttemptLedgerRow(;
                attempt_id = "confirmation-001",
                recorded_at_utc = "2026-07-12T17:00:00Z",
                seed_partition = :confirmation,
                seed = 44_001,
                modality = :long_noisy,
                configuration_id = "frozen-configuration",
                git_sha,
                git_dirty = false,
                manifest_sha256 = sha_a,
                configuration_sha256 = sha_b,
                fixture_sha256 = sha_c,
                truth_reference_sha256 = confirmation_truth_sha256[44_001],
                input_bundle_sha256 = cell_44_001.input_bundle_sha256,
                consumed_input_sha256 = cell_44_001.input_bundle_sha256,
                error_profile_sha256 = cell_44_001.error_profile_sha256,
                comparator_id = "rhizomorph-native",
                comparator_version = git_sha,
                command_sha256 = sha_c,
                comparator_configuration_sha256 = sha_c,
                status = :passed,
                result_path = "results/confirmation-001",
                result_sha256 = sha_a,
                notes = "first frozen confirmation attempt",
            )
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin(
            "require a development freeze", missing_freeze_message)

        fingerprint = stage2_attempt_ledger_fingerprint(ledger_path)
        seal = seal_stage2_development_ledger!(
            ledger_path;
            frozen_at_utc = "2026-07-12T16:30:00Z",
            manifest_path,
            configuration_path,
            fixture_path,
            git_sha,
            manifest_sha256 = sha_a,
            configuration_sha256 = sha_b,
            fixture_sha256 = sha_c,
            git_dirty = false,
        )
        freeze = seal.freeze
        freeze_path = seal.path
        Test.@test freeze.development_ledger_sha256 == fingerprint.sha256
        Test.@test freeze_path == stage2_development_freeze_path(ledger_path)
        Test.@test_throws ErrorException seal_stage2_development_ledger!(
            ledger_path;
            frozen_at_utc = "2026-07-12T16:31:00Z",
            manifest_path,
            configuration_path,
            fixture_path,
            git_sha,
            manifest_sha256 = sha_a,
            configuration_sha256 = sha_b,
            fixture_sha256 = sha_c,
            git_dirty = false,
        )
        Test.@test_throws ErrorException append_stage2_attempt!(
            ledger_path, second_row)
        confirmation = Stage2AttemptLedgerRow(;
            attempt_id = "confirmation-001",
            recorded_at_utc = "2026-07-12T17:00:00Z",
            seed_partition = :confirmation,
            seed = 44_001,
            modality = :long_noisy,
            configuration_id = "frozen-configuration",
            git_sha,
            git_dirty = false,
            manifest_sha256 = sha_a,
            configuration_sha256 = sha_b,
            fixture_sha256 = sha_c,
            truth_reference_sha256 = confirmation_truth_sha256[44_001],
            input_bundle_sha256 = cell_44_001.input_bundle_sha256,
            consumed_input_sha256 = cell_44_001.input_bundle_sha256,
            error_profile_sha256 = cell_44_001.error_profile_sha256,
            comparator_id = "metaflye-strainy",
            comparator_version =
                manifest["comparators"]["metaflye-strainy"]["version"],
            command_sha256 = bytes2hex(SHA.sha256(
                manifest["comparators"]["metaflye-strainy"]["command"])),
            comparator_configuration_sha256 = bytes2hex(SHA.sha256(
                manifest["comparators"]["metaflye-strainy"]["configuration"])),
            status = :passed,
            result_path = "results/confirmation-001",
            result_sha256 = sha_a,
            development_freeze_sha256 = seal.sha256,
            notes = "first frozen confirmation attempt",
            development_freeze = freeze,
        )
        confirmation_ledger_path = stage2_confirmation_ledger_path(ledger_path)
        reservation = _replace_ledger_row(
            confirmation;
            development_freeze = freeze,
            recorded_at_utc = "2026-07-12T16:59:00Z",
            status = :reserved,
            result_path = "",
            result_sha256 = "",
            notes = "reserve confirmation cell",
        )
        append_stage2_confirmation_attempt!(
            reservation;
            development_ledger_path = ledger_path,
        )
        reservation_nonce_sha256 =
            _validate_confirmation_reservation_receipt(
                confirmation_ledger_path, reservation)
        absolute_confirmation_result = _stage2_result_artifact_path(
            confirmation_ledger_path, confirmation.result_path)
        mkpath(dirname(absolute_confirmation_result))
        open(absolute_confirmation_result, "w") do stream
            TOML.print(stream, Dict(
                "schema" => "rhizomorph-stage2-test-completion-v1",
                "reservation_nonce_sha256" => reservation_nonce_sha256,
            ))
        end
        confirmation = _replace_ledger_row(
            confirmation;
            development_freeze = freeze,
            result_sha256 = bytes2hex(
                SHA.sha256(read(absolute_confirmation_result))),
        )
        append_stage2_confirmation_attempt!(
            confirmation;
            development_ledger_path = ledger_path,
        )
        Test.@test countlines(ledger_path) == 3
        Test.@test countlines(confirmation_ledger_path) == 3

        wrong_nonce_result_path = "results/precreated-wrong-nonce.toml"
        absolute_wrong_nonce_result = _stage2_result_artifact_path(
            confirmation_ledger_path, wrong_nonce_result_path)
        open(absolute_wrong_nonce_result, "w") do stream
            TOML.print(stream, Dict(
                "schema" => "rhizomorph-stage2-test-completion-v1",
                "reservation_nonce_sha256" => repeat("f", 64),
            ))
        end
        wrong_nonce_comparator = manifest["comparators"]["hairsplitter"]
        wrong_nonce_reservation = _replace_ledger_row(
            reservation;
            development_freeze = freeze,
            attempt_id = "precreated-wrong-nonce",
            comparator_id = "hairsplitter",
            comparator_version = wrong_nonce_comparator["version"],
            command_sha256 = _sha256(wrong_nonce_comparator["command"]),
            comparator_configuration_sha256 =
                _sha256(wrong_nonce_comparator["configuration"]),
        )
        wrong_nonce_completion = _replace_ledger_row(
            wrong_nonce_reservation;
            development_freeze = freeze,
            recorded_at_utc = "2026-07-12T17:01:00Z",
            status = :passed,
            result_path = wrong_nonce_result_path,
            result_sha256 = bytes2hex(
                SHA.sha256(read(absolute_wrong_nonce_result))),
            notes = "pre-created result must not satisfy reservation",
        )
        append_stage2_confirmation_attempt!(
            wrong_nonce_reservation;
            development_ledger_path = ledger_path,
        )
        confirmation_before_wrong_nonce = read(confirmation_ledger_path)
        Test.@test_throws ErrorException append_stage2_confirmation_attempt!(
            wrong_nonce_completion;
            development_ledger_path = ledger_path,
        )
        Test.@test read(confirmation_ledger_path) ==
                   confirmation_before_wrong_nonce

        callback_completion = _replace_ledger_row(
            confirmation;
            development_freeze = freeze,
            attempt_id = "confirmation-callback",
            seed = 44_002,
            truth_reference_sha256 = confirmation_truth_sha256[44_002],
            input_bundle_sha256 = cell_44_002.input_bundle_sha256,
            consumed_input_sha256 = cell_44_002.input_bundle_sha256,
            error_profile_sha256 = cell_44_002.error_profile_sha256,
            result_path = "results/confirmation-callback",
            notes = "callback completion",
        )
        callback_reservation = _replace_ledger_row(
            callback_completion;
            development_freeze = freeze,
            recorded_at_utc = "2026-07-12T16:59:00Z",
            status = :reserved,
            result_path = "",
            result_sha256 = "",
            notes = "callback reservation",
        )
        reservation_was_visible = Ref(false)
        completed = run_stage2_confirmation_attempt!(
            callback_reservation,
            (canonical_ledger, reservation_nonce_sha256) -> begin
            events = load_stage2_attempt_ledger(
                canonical_ledger; development_freeze = freeze)
            reservation_was_visible[] = any(event ->
                event.attempt_id == callback_reservation.attempt_id &&
                event.status == :reserved,
                events,
            )
            absolute_callback_result = _stage2_result_artifact_path(
                canonical_ledger, callback_completion.result_path)
            mkpath(dirname(absolute_callback_result))
            open(absolute_callback_result, "w") do stream
                TOML.print(stream, Dict(
                    "schema" => "rhizomorph-stage2-test-completion-v1",
                    "reservation_nonce_sha256" => reservation_nonce_sha256,
                ))
            end
            return _replace_ledger_row(
                callback_completion;
                development_freeze = freeze,
                result_sha256 = bytes2hex(
                    SHA.sha256(read(absolute_callback_result))),
            )
            end;
            development_ledger_path = ledger_path,
        )
        Test.@test reservation_was_visible[]
        Test.@test completed.status == :passed

        wrong_freeze_reservation = _replace_ledger_row(
            callback_reservation;
            development_freeze = freeze,
            attempt_id = "wrong-freeze",
            seed = 44_003,
            truth_reference_sha256 = confirmation_truth_sha256[44_003],
            input_bundle_sha256 = cell_44_003.input_bundle_sha256,
            consumed_input_sha256 = cell_44_003.input_bundle_sha256,
            error_profile_sha256 = cell_44_003.error_profile_sha256,
            development_freeze_sha256 = repeat("e", 64),
        )
        Test.@test_throws ErrorException append_stage2_confirmation_attempt!(
            wrong_freeze_reservation;
            development_ledger_path = ledger_path,
        )

        mkdir(ledger_path * ".lock")
        locked_reservation = _replace_ledger_row(
            callback_reservation;
            development_freeze = freeze,
            attempt_id = "locked-attempt",
            seed = 44_003,
            truth_reference_sha256 = confirmation_truth_sha256[44_003],
            input_bundle_sha256 = cell_44_003.input_bundle_sha256,
            consumed_input_sha256 = cell_44_003.input_bundle_sha256,
            error_profile_sha256 = cell_44_003.error_profile_sha256,
        )
        Test.@test_throws ErrorException append_stage2_confirmation_attempt!(
            locked_reservation;
            development_ledger_path = ledger_path,
        )
        rm(ledger_path * ".lock"; recursive = true)

        throwing_reservation = _replace_ledger_row(
            callback_reservation;
            development_freeze = freeze,
            attempt_id = "throwing-attempt",
            seed = 44_003,
            truth_reference_sha256 = confirmation_truth_sha256[44_003],
            input_bundle_sha256 = cell_44_003.input_bundle_sha256,
            consumed_input_sha256 = cell_44_003.input_bundle_sha256,
            error_profile_sha256 = cell_44_003.error_profile_sha256,
        )
        Test.@test_throws ErrorException run_stage2_confirmation_attempt!(
            throwing_reservation,
            (_, _) -> error("synthetic runner failure");
            development_ledger_path = ledger_path,
        )
        events_after_throw = load_stage2_attempt_ledger(
            confirmation_ledger_path; development_freeze = freeze)
        Test.@test any(event ->
            event.attempt_id == throwing_reservation.attempt_id &&
            event.status == :error,
            events_after_throw,
        )
        replacement_after_error = _replace_ledger_row(
            throwing_reservation;
            development_freeze = freeze,
            attempt_id = "replacement-after-error",
        )
        Test.@test_throws ErrorException append_stage2_confirmation_attempt!(
            replacement_after_error;
            development_ledger_path = ledger_path,
        )

        malformed_time_message = try
            Stage2DevelopmentFreeze(;
                frozen_at_utc = "not-a-time",
                development_ledger_path = ledger_path,
                confirmation_ledger_path,
                manifest_path,
                configuration_path,
                fixture_path,
                git_sha,
                manifest_sha256 = sha_a,
                configuration_sha256 = sha_b,
                fixture_sha256 = sha_c,
                development_ledger_sha256 = fingerprint.sha256,
                development_ledger_rows = fingerprint.rows,
                final_attempt_id = fingerprint.final_attempt_id,
                git_dirty = false,
            )
            ""
        catch error
            sprint(showerror, error)
        end
        Test.@test occursin("UTC format", malformed_time_message)
        Test.@test_throws ErrorException Stage2DevelopmentFreeze(;
            frozen_at_utc = "2026-07-12T16:30:00Z",
            development_ledger_path = ledger_path,
            confirmation_ledger_path,
            manifest_path,
            configuration_path,
            fixture_path,
            git_sha,
            manifest_sha256 = sha_a,
            configuration_sha256 = sha_b,
            fixture_sha256 = sha_c,
            development_ledger_sha256 = fingerprint.sha256,
            development_ledger_rows = fingerprint.rows,
            final_attempt_id = fingerprint.final_attempt_id,
            git_dirty = true,
        )
        Test.@test_throws MethodError Stage2AttemptLedgerRow((
            getfield(first_row, index)
            for index in 1:fieldcount(Stage2AttemptLedgerRow)
        )...)
    end
end

Test.@testset "Stage-2 development seal requires strict PASS evidence" begin
    mktempdir() do directory
        digest = repeat("a", 64)
        git_sha = repeat("b", 40)
        manifest_path = joinpath(directory, "manifest.toml")
        configuration_path = joinpath(directory, "configuration.toml")
        fixture_path = joinpath(directory, "fixture.toml")
        write(manifest_path, "schema_version = 1\n")
        write(configuration_path, "configuration_id = \"failed\"\n")
        write(fixture_path,
            "schema = \"rhizomorph-stage2-fixture-index-v1\"\n")
        manifest_sha256 = bytes2hex(SHA.sha256(read(manifest_path)))
        configuration_sha256 = bytes2hex(SHA.sha256(read(configuration_path)))
        fixture_sha256 = bytes2hex(SHA.sha256(read(fixture_path)))
        ledger_path = joinpath(directory, "failed-development.tsv")
        failed_row = Stage2AttemptLedgerRow(;
            attempt_id = "failed-final",
            recorded_at_utc = "2026-07-12T15:00:00Z",
            seed_partition = :development,
            seed = 43_003,
            modality = :long_noisy,
            configuration_id = "failed-configuration",
            git_sha,
            git_dirty = false,
            manifest_sha256,
            configuration_sha256,
            fixture_sha256,
            input_bundle_sha256 = digest,
            consumed_input_sha256 = digest,
            error_profile_sha256 = digest,
            comparator_id = "rhizomorph-native",
            comparator_version = git_sha,
            command_sha256 = digest,
            status = :failed,
            result_path = "results/failed.toml",
            result_sha256 = digest,
            notes = "failed strict development gate",
        )
        append_stage2_attempt!(ledger_path, failed_row)
        Test.@test_throws ErrorException seal_stage2_development_ledger!(
            ledger_path;
            frozen_at_utc = "2026-07-12T15:30:00Z",
            manifest_path,
            configuration_path,
            fixture_path,
            git_sha,
            manifest_sha256,
            configuration_sha256,
            fixture_sha256,
            git_dirty = false,
        )
        Test.@test !isfile(stage2_development_freeze_path(ledger_path))

        attestation_ledger_path = joinpath(directory, "self-attested.tsv")
        attestation_path = joinpath(directory, "self-attested.toml")
        write(attestation_path,
            "strict_gate_passed = true\ngit_sha = \"$(git_sha)\"\n")
        self_attested_row = _replace_ledger_row(
            failed_row;
            attempt_id = "self-attested-final",
            status = :passed,
            result_path = basename(attestation_path),
            result_sha256 = bytes2hex(SHA.sha256(read(attestation_path))),
        )
        append_stage2_attempt!(attestation_ledger_path, self_attested_row)
        Test.@test_throws ErrorException seal_stage2_development_ledger!(
            attestation_ledger_path;
            frozen_at_utc = "2026-07-12T15:31:00Z",
            manifest_path,
            configuration_path,
            fixture_path,
            git_sha,
            manifest_sha256,
            configuration_sha256,
            fixture_sha256,
            git_dirty = false,
        )
        Test.@test !isfile(
            stage2_development_freeze_path(attestation_ledger_path))
    end
end
