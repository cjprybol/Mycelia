import CodecZlib
import Dates
import FASTX
import Random
import SHA
import Statistics
import TOML

"""
Development seeds are the only seeds available while tuning Stage-2.

The confirmation seeds are a frozen holdout. They must remain unrun until the
development configuration and decision rules have been locked.
"""
const STAGE2_DEVELOPMENT_SEEDS = (43_003, 43_006, 43_007, 43_008, 43_009)
const STAGE2_CONFIRMATION_SEEDS = (44_001, 44_002, 44_003)

const STAGE2_ATTEMPT_LEDGER_COLUMNS = (
    "attempt_id",
    "recorded_at_utc",
    "seed_partition",
    "seed",
    "modality",
    "configuration_id",
    "git_sha",
    "git_dirty",
    "manifest_sha256",
    "configuration_sha256",
    "fixture_sha256",
    "truth_reference_sha256",
    "input_bundle_sha256",
    "consumed_input_sha256",
    "error_profile_sha256",
    "comparator_id",
    "comparator_version",
    "command_sha256",
    "comparator_configuration_sha256",
    "status",
    "result_path",
    "result_sha256",
    "development_freeze_sha256",
    "notes",
)

const STAGE2_TERMINAL_ATTEMPT_STATUSES = (:passed, :failed, :error)
const STAGE2_ATTEMPT_STATUSES = (:reserved, STAGE2_TERMINAL_ATTEMPT_STATUSES...)
const STAGE2_MODALITIES = (
    :long_noisy,
    :long_accurate,
    :short_paired,
    :short_single,
    :hybrid_ont_short,
    :hybrid_hifi_ont,
)

"""Length-only summary of the records in one FASTA or FASTQ file."""
struct RecordLengthSummary
    record_count::Int
    total_bases::Int
    minimum_bases::Union{Missing, Int}
    maximum_bases::Union{Missing, Int}
    mean_bases::Union{Missing, Float64}
end

"""Auditable components of an assembly-to-truth discrepancy count."""
struct AssemblyDiscrepancyComponents
    total_reference_bases::Int
    total_query_bases::Int
    aligned_reference_bases::Int
    aligned_query_bases::Int
    aligned_error_bases::Int
    unaligned_reference_bases::Int
    unaligned_query_bases::Int
    total_discrepancy_bases::Int
    total_errors_per_100kbp::Float64
end

"""Coverage summary for candidate-to-truth assignments."""
struct CandidateAssignmentSummary
    expected_truth_ids::Vector{String}
    assigned_truth_ids::Vector{String}
    unique_assigned_truth_ids::Vector{String}
    duplicate_truth_ids::Vector{String}
    missing_truth_ids::Vector{String}
    unexpected_truth_ids::Vector{String}
    exact_recovery::Bool
end

"""One modality/seed comparison consumed by the release-gate evaluator."""
struct Stage2ReleaseObservation
    modality::Symbol
    seed::Int
    attempt_id::String
    reservation_nonce_sha256::String
    native_result_path::String
    native_result_sha256::String
    native_nga50::Float64
    comparator_nga50::Float64
    native_genome_fraction::Float64
    comparator_genome_fraction::Float64
    native_misassemblies::Int
    comparator_misassemblies::Int
    native_error_per_100kbp::Float64
    comparator_error_per_100kbp::Float64
    comparator_id::String
    comparator_version::String
    comparator_result_path::String
    comparator_result_sha256::String
    comparator_command_sha256::String
    comparator_configuration_sha256::String
    evaluator_id::String
    evaluator_version::String
    evaluator_command_sha256::String
    development_freeze_sha256::String
    manifest_sha256::String
    configuration_sha256::String
    fixture_sha256::String
    input_component_paths::Dict{String, String}
    input_component_sha256::Dict{String, String}
    truth_reference_path::String
    truth_reference_sha256::String
    metric_report_path::String
    metric_report_sha256::String
    input_bundle_sha256::String
    native_consumed_input_sha256::String
    comparator_consumed_input_sha256::String
    comparator_consumed_input_components::Vector{String}
    assignment::CandidateAssignmentSummary
    truth_reference_coverages::Vector{Float64}
    abundance_mae::Float64
    comparator_abundance_mae::Float64
    strain_count_mode::Int
    strain_count_set_contains_truth::Bool
    primary_rank_agreement::Bool
    ensemble_noninferior_or_abstained::Bool
end

"""Fail-closed release-gate result with auditable component decisions."""
struct Stage2ReleaseGateResult
    passed::Bool
    observation_gates::Dict{Tuple{Symbol, Int}, Dict{String, Bool}}
    modality_gates::Dict{Symbol, Dict{String, Bool}}
    development_freeze_sha256::String
    confirmation_ledger_sha256::String
    attempt_result_sha256::Dict{String, String}
    diagnostic_attempt_status::Dict{Tuple{Symbol, Int, String}, Symbol}
    diagnostic_result_sha256::Dict{String, String}
end

"""Immutable hashes authorizing the first confirmation run after development."""
struct Stage2DevelopmentFreeze
    frozen_at_utc::String
    development_ledger_path::String
    confirmation_ledger_path::String
    manifest_path::String
    configuration_path::String
    fixture_path::String
    git_sha::String
    manifest_sha256::String
    configuration_sha256::String
    fixture_sha256::String
    development_ledger_sha256::String
    development_ledger_rows::Int
    final_attempt_id::String
    git_dirty::Bool

    function Stage2DevelopmentFreeze(
            frozen_at_utc::String,
            development_ledger_path::String,
            confirmation_ledger_path::String,
            manifest_path::String,
            configuration_path::String,
            fixture_path::String,
            git_sha::String,
            manifest_sha256::String,
            configuration_sha256::String,
            fixture_sha256::String,
            development_ledger_sha256::String,
            development_ledger_rows::Int,
            final_attempt_id::String,
            git_dirty::Bool
    )
        git_dirty && error("development freeze requires git_dirty=false")
        development_ledger_rows > 0 ||
            error("development_ledger_rows must be positive")
        canonical_development_path = abspath(normpath(development_ledger_path))
        canonical_confirmation_path = abspath(normpath(confirmation_ledger_path))
        canonical_development_path != canonical_confirmation_path ||
            error("development and confirmation ledgers must be separate")
        return new(
            _validate_utc_timestamp(frozen_at_utc, "frozen_at_utc"),
            canonical_development_path,
            canonical_confirmation_path,
            abspath(normpath(manifest_path)),
            abspath(normpath(configuration_path)),
            abspath(normpath(fixture_path)),
            _validate_hex_digest(git_sha, "git_sha", (40, 64)),
            _validate_hex_digest(manifest_sha256, "manifest_sha256", (64,)),
            _validate_hex_digest(
                configuration_sha256, "configuration_sha256", (64,)),
            _validate_hex_digest(fixture_sha256, "fixture_sha256", (64,)),
            _validate_hex_digest(
                development_ledger_sha256, "development_ledger_sha256", (64,)),
            development_ledger_rows,
            _validate_ledger_text(final_attempt_id, "final_attempt_id"),
            git_dirty,
        )
    end
end

"""Canonical on-disk seal and digest for one frozen development ledger."""
struct Stage2DevelopmentSeal
    freeze::Stage2DevelopmentFreeze
    path::String
    sha256::String
end

"""One immutable row in the Stage-2 tuning-attempt ledger."""
struct Stage2AttemptLedgerRow
    attempt_id::String
    recorded_at_utc::String
    seed_partition::Symbol
    seed::Int
    modality::Symbol
    configuration_id::String
    git_sha::String
    git_dirty::Bool
    manifest_sha256::String
    configuration_sha256::String
    fixture_sha256::String
    truth_reference_sha256::String
    input_bundle_sha256::String
    consumed_input_sha256::String
    error_profile_sha256::String
    comparator_id::String
    comparator_version::String
    command_sha256::String
    comparator_configuration_sha256::String
    status::Symbol
    result_path::String
    result_sha256::String
    development_freeze_sha256::String
    notes::String

    function Stage2AttemptLedgerRow(
            attempt_id::String,
            recorded_at_utc::String,
            seed_partition::Symbol,
            seed::Int,
            modality::Symbol,
            configuration_id::String,
            git_sha::String,
            git_dirty::Bool,
            manifest_sha256::String,
            configuration_sha256::String,
            fixture_sha256::String,
            truth_reference_sha256::String,
            input_bundle_sha256::String,
            consumed_input_sha256::String,
            error_profile_sha256::String,
            comparator_id::String,
            comparator_version::String,
            command_sha256::String,
            comparator_configuration_sha256::String,
            status::Symbol,
            result_path::String,
            result_sha256::String,
            development_freeze_sha256::String,
            notes::String,
            development_freeze::Union{Nothing, Stage2DevelopmentFreeze}
    )
        seed_partition in (:development, :confirmation) ||
            error("seed_partition must be :development or :confirmation")
        allowed_seeds = seed_partition == :development ?
                        STAGE2_DEVELOPMENT_SEEDS : STAGE2_CONFIRMATION_SEEDS
        seed in allowed_seeds ||
            error("seed $(seed) is not frozen in the $(seed_partition) partition")
        modality in STAGE2_MODALITIES ||
            error("unknown Stage-2 modality: $(modality)")
        status in STAGE2_ATTEMPT_STATUSES ||
            error("unknown attempt status: $(status)")
        seed_partition == :development && status == :reserved &&
            error("development attempts cannot use reserved status")
        validated_git_sha = _validate_hex_digest(git_sha, "git_sha", (40, 64))
        validated_manifest = _validate_hex_digest(
            manifest_sha256, "manifest_sha256", (64,))
        validated_configuration = _validate_hex_digest(
            configuration_sha256, "configuration_sha256", (64,))
        validated_fixture = _validate_hex_digest(
            fixture_sha256, "fixture_sha256", (64,))
        if seed_partition == :confirmation
            development_freeze === nothing &&
                error("confirmation rows require a development freeze")
            git_dirty && error("confirmation rows require git_dirty=false")
            validated_git_sha == development_freeze.git_sha ||
                error("confirmation git SHA differs from freeze")
            validated_manifest == development_freeze.manifest_sha256 ||
                error("confirmation manifest hash differs from freeze")
            validated_configuration == development_freeze.configuration_sha256 ||
                error("confirmation configuration hash differs from freeze")
            validated_fixture == development_freeze.fixture_sha256 ||
                error("confirmation fixture hash differs from freeze")
        end
        validated_result_path = status == :reserved ?
                                _validate_empty_text(result_path, "result_path") :
                                _validate_ledger_text(result_path, "result_path")
        validated_result_sha256 = status == :reserved ?
                                  _validate_empty_text(
            result_sha256, "result_sha256") :
                                  _validate_hex_digest(
            result_sha256, "result_sha256", (64,))
        return new(
            _validate_ledger_text(attempt_id, "attempt_id"),
            _validate_utc_timestamp(recorded_at_utc, "recorded_at_utc"),
            seed_partition,
            seed,
            modality,
            _validate_ledger_text(configuration_id, "configuration_id"),
            validated_git_sha,
            git_dirty,
            validated_manifest,
            validated_configuration,
            validated_fixture,
            seed_partition == :confirmation ?
            _validate_hex_digest(
                truth_reference_sha256, "truth_reference_sha256", (64,)) :
            _validate_empty_text(
                truth_reference_sha256, "truth_reference_sha256"),
            _validate_hex_digest(
                input_bundle_sha256, "input_bundle_sha256", (64,)),
            _validate_hex_digest(
                consumed_input_sha256, "consumed_input_sha256", (64,)),
            _validate_hex_digest(
                error_profile_sha256, "error_profile_sha256", (64,)),
            _validate_ledger_text(comparator_id, "comparator_id"),
            _validate_ledger_text(comparator_version, "comparator_version"),
            _validate_hex_digest(command_sha256, "command_sha256", (64,)),
            seed_partition == :confirmation ?
            _validate_hex_digest(
                comparator_configuration_sha256,
                "comparator_configuration_sha256",
                (64,),
            ) : _validate_empty_text(
                comparator_configuration_sha256,
                "comparator_configuration_sha256",
            ),
            status,
            validated_result_path,
            validated_result_sha256,
            seed_partition == :confirmation ?
            _validate_hex_digest(
                development_freeze_sha256,
                "development_freeze_sha256",
                (64,),
            ) : _validate_empty_text(
                development_freeze_sha256, "development_freeze_sha256"),
            _validate_ledger_text(notes, "notes"),
        )
    end
end

function _record_length_summary(lengths::Vector{Int})::RecordLengthSummary
    isempty(lengths) && return RecordLengthSummary(0, 0, missing, missing, missing)
    return RecordLengthSummary(
        length(lengths),
        sum(lengths),
        minimum(lengths),
        maximum(lengths),
        sum(lengths) / length(lengths),
    )
end

function _open_sequence_stream(path::AbstractString)::IO
    raw_stream = open(path, "r")
    if endswith(lowercase(String(path)), ".gz")
        return CodecZlib.GzipDecompressorStream(raw_stream)
    end
    return raw_stream
end

"""Summarize FASTA record lengths without invoking an external program."""
function summarize_fasta_record_lengths(path::AbstractString)::RecordLengthSummary
    stream = _open_sequence_stream(path)
    reader = FASTX.FASTA.Reader(stream)
    try
        lengths = Int[
            length(FASTX.FASTA.sequence(String, record)) for record in reader
        ]
        return _record_length_summary(lengths)
    finally
        close(reader)
    end
end

"""Summarize FASTQ record lengths without invoking an external program."""
function summarize_fastq_record_lengths(path::AbstractString)::RecordLengthSummary
    stream = _open_sequence_stream(path)
    reader = FASTX.FASTQ.Reader(stream)
    try
        lengths = Int[
            length(FASTX.FASTQ.sequence(String, record)) for record in reader
        ]
        return _record_length_summary(lengths)
    finally
        close(reader)
    end
end

"""
Split a full assembly discrepancy into aligned errors and unaligned tails.

The total discrepancy is the sum of aligned error bases, unaligned reference
bases, and unaligned query bases. The normalized rate uses the full reference
length as its denominator.
"""
function assembly_discrepancy_components(
    total_reference_bases::Integer,
    total_query_bases::Integer,
    aligned_reference_bases::Integer,
    aligned_query_bases::Integer,
    aligned_error_bases::Integer,
)::AssemblyDiscrepancyComponents
    counts = (
        total_reference_bases,
        total_query_bases,
        aligned_reference_bases,
        aligned_query_bases,
        aligned_error_bases,
    )
    any(<(0), counts) && error("assembly discrepancy counts must be nonnegative")
    total_reference_bases > 0 ||
        error("total_reference_bases must be positive")
    aligned_reference_bases <= total_reference_bases ||
        error("aligned_reference_bases exceeds total_reference_bases")
    aligned_query_bases <= total_query_bases ||
        error("aligned_query_bases exceeds total_query_bases")

    total_reference = Int(total_reference_bases)
    total_query = Int(total_query_bases)
    aligned_reference = Int(aligned_reference_bases)
    aligned_query = Int(aligned_query_bases)
    aligned_errors = Int(aligned_error_bases)
    unaligned_reference = total_reference - aligned_reference
    unaligned_query = total_query - aligned_query
    total_discrepancy = aligned_errors + unaligned_reference + unaligned_query
    return AssemblyDiscrepancyComponents(
        total_reference,
        total_query,
        aligned_reference,
        aligned_query,
        aligned_errors,
        unaligned_reference,
        unaligned_query,
        total_discrepancy,
        total_discrepancy / total_reference * 100_000.0,
    )
end

function _ordered_unique(values::Vector{String})::Vector{String}
    seen = Set{String}()
    unique_values = String[]
    for value in values
        if value ∉ seen
            push!(seen, value)
            push!(unique_values, value)
        end
    end
    return unique_values
end

"""
Summarize candidate assignments without hiding truth collisions or dropouts.

`exact_recovery` is true only for a one-to-one assignment covering every
expected truth and no unexpected truth.
"""
function summarize_candidate_assignments(
    expected_truth_ids::AbstractVector{<:AbstractString},
    assigned_truth_ids::AbstractVector{<:AbstractString},
)::CandidateAssignmentSummary
    expected = String.(expected_truth_ids)
    assigned = String.(assigned_truth_ids)
    any(isempty, expected) && error("expected truth IDs must be nonempty")
    any(isempty, assigned) && error("assigned truth IDs must be nonempty")
    length(unique(expected)) == length(expected) ||
        error("expected truth IDs must be unique")

    counts = Dict{String, Int}()
    for truth_id in assigned
        counts[truth_id] = get(counts, truth_id, 0) + 1
    end
    unique_assigned = _ordered_unique(assigned)
    duplicates = [truth_id for truth_id in unique_assigned if counts[truth_id] > 1]
    assigned_set = Set(assigned)
    expected_set = Set(expected)
    missing_truths = [truth_id for truth_id in expected if truth_id ∉ assigned_set]
    unexpected_truths = [
        truth_id for truth_id in unique_assigned if truth_id ∉ expected_set
    ]
    exact_recovery = isempty(duplicates) && isempty(missing_truths) &&
                     isempty(unexpected_truths) && length(assigned) == length(expected)
    return CandidateAssignmentSummary(
        expected,
        assigned,
        unique_assigned,
        duplicates,
        missing_truths,
        unexpected_truths,
        exact_recovery,
    )
end

function _assignment_signature(
        summary::CandidateAssignmentSummary
)::NamedTuple
    return (
        expected_truth_ids = Tuple(summary.expected_truth_ids),
        assigned_truth_ids = Tuple(summary.assigned_truth_ids),
        unique_assigned_truth_ids = Tuple(summary.unique_assigned_truth_ids),
        duplicate_truth_ids = Tuple(summary.duplicate_truth_ids),
        missing_truth_ids = Tuple(summary.missing_truth_ids),
        unexpected_truth_ids = Tuple(summary.unexpected_truth_ids),
        exact_recovery = summary.exact_recovery,
    )
end

function _validate_assignment_summary(
        summary::CandidateAssignmentSummary
)::Nothing
    canonical = summarize_candidate_assignments(
        summary.expected_truth_ids, summary.assigned_truth_ids)
    isequal(_assignment_signature(summary), _assignment_signature(canonical)) ||
        error("candidate assignment summary is internally inconsistent")
    return nothing
end

function _native_scientific_signature(
        observation::Stage2ReleaseObservation
)::NamedTuple
    return (
        native_result_path = observation.native_result_path,
        native_result_sha256 = observation.native_result_sha256,
        native_nga50 = observation.native_nga50,
        native_genome_fraction = observation.native_genome_fraction,
        native_misassemblies = observation.native_misassemblies,
        native_error_per_100kbp = observation.native_error_per_100kbp,
        assignment = _assignment_signature(observation.assignment),
        truth_reference_coverages = Tuple(observation.truth_reference_coverages),
        abundance_mae = observation.abundance_mae,
        strain_count_mode = observation.strain_count_mode,
        strain_count_set_contains_truth =
            observation.strain_count_set_contains_truth,
        primary_rank_agreement = observation.primary_rank_agreement,
        ensemble_noninferior_or_abstained =
            observation.ensemble_noninferior_or_abstained,
        input_bundle_sha256 = observation.input_bundle_sha256,
        input_component_paths = Tuple(sort!(collect(
            observation.input_component_paths))),
        input_component_sha256 = Tuple(sort!(collect(
            observation.input_component_sha256))),
        truth_reference_path = observation.truth_reference_path,
        truth_reference_sha256 = observation.truth_reference_sha256,
        native_consumed_input_sha256 =
            observation.native_consumed_input_sha256,
    )
end

function _stage2_release_observation_table(
        observation::Stage2ReleaseObservation
)::Dict{String, Any}
    return Dict{String, Any}(
        "modality" => String(observation.modality),
        "seed" => observation.seed,
        "attempt_id" => observation.attempt_id,
        "reservation_nonce_sha256" => observation.reservation_nonce_sha256,
        "native_result_path" => observation.native_result_path,
        "native_result_sha256" => observation.native_result_sha256,
        "native_nga50" => observation.native_nga50,
        "comparator_nga50" => observation.comparator_nga50,
        "native_genome_fraction" => observation.native_genome_fraction,
        "comparator_genome_fraction" => observation.comparator_genome_fraction,
        "native_misassemblies" => observation.native_misassemblies,
        "comparator_misassemblies" => observation.comparator_misassemblies,
        "native_error_per_100kbp" => observation.native_error_per_100kbp,
        "comparator_error_per_100kbp" =>
            observation.comparator_error_per_100kbp,
        "comparator_id" => observation.comparator_id,
        "comparator_version" => observation.comparator_version,
        "comparator_result_path" => observation.comparator_result_path,
        "comparator_result_sha256" => observation.comparator_result_sha256,
        "comparator_command_sha256" => observation.comparator_command_sha256,
        "comparator_configuration_sha256" =>
            observation.comparator_configuration_sha256,
        "evaluator_id" => observation.evaluator_id,
        "evaluator_version" => observation.evaluator_version,
        "evaluator_command_sha256" => observation.evaluator_command_sha256,
        "development_freeze_sha256" => observation.development_freeze_sha256,
        "manifest_sha256" => observation.manifest_sha256,
        "configuration_sha256" => observation.configuration_sha256,
        "fixture_sha256" => observation.fixture_sha256,
        "input_component_paths" => observation.input_component_paths,
        "input_component_sha256" => observation.input_component_sha256,
        "truth_reference_path" => observation.truth_reference_path,
        "truth_reference_sha256" => observation.truth_reference_sha256,
        "metric_report_path" => observation.metric_report_path,
        "metric_report_sha256" => observation.metric_report_sha256,
        "input_bundle_sha256" => observation.input_bundle_sha256,
        "native_consumed_input_sha256" =>
            observation.native_consumed_input_sha256,
        "comparator_consumed_input_sha256" =>
            observation.comparator_consumed_input_sha256,
        "comparator_consumed_input_components" =>
            observation.comparator_consumed_input_components,
        "truth_reference_coverages" => observation.truth_reference_coverages,
        "abundance_mae" => observation.abundance_mae,
        "comparator_abundance_mae" => observation.comparator_abundance_mae,
        "strain_count_mode" => observation.strain_count_mode,
        "strain_count_set_contains_truth" =>
            observation.strain_count_set_contains_truth,
        "primary_rank_agreement" => observation.primary_rank_agreement,
        "ensemble_noninferior_or_abstained" =>
            observation.ensemble_noninferior_or_abstained,
        "assignment" => Dict{String, Any}(
            "expected_truth_ids" => observation.assignment.expected_truth_ids,
            "assigned_truth_ids" => observation.assignment.assigned_truth_ids,
            "unique_assigned_truth_ids" =>
                observation.assignment.unique_assigned_truth_ids,
            "duplicate_truth_ids" => observation.assignment.duplicate_truth_ids,
            "missing_truth_ids" => observation.assignment.missing_truth_ids,
            "unexpected_truth_ids" => observation.assignment.unexpected_truth_ids,
            "exact_recovery" => observation.assignment.exact_recovery,
        ),
    )
end

function write_stage2_release_observation(
        path::AbstractString,
        observation::Stage2ReleaseObservation
)::String
    output = abspath(normpath(String(path)))
    mkpath(dirname(output))
    open(output, "w") do stream
        TOML.print(stream, _stage2_release_observation_table(observation))
    end
    return output
end

function _parse_stage2_release_observation(
        bytes::Vector{UInt8}
)::Stage2ReleaseObservation
    table = TOML.parse(String(bytes))
    assignment = table["assignment"]
    summary = CandidateAssignmentSummary(
        String.(assignment["expected_truth_ids"]),
        String.(assignment["assigned_truth_ids"]),
        String.(assignment["unique_assigned_truth_ids"]),
        String.(assignment["duplicate_truth_ids"]),
        String.(assignment["missing_truth_ids"]),
        String.(assignment["unexpected_truth_ids"]),
        assignment["exact_recovery"],
    )
    _validate_assignment_summary(summary)
    return Stage2ReleaseObservation(
        Symbol(table["modality"]),
        Int(table["seed"]),
        String(table["attempt_id"]),
        _validate_hex_digest(
            String(table["reservation_nonce_sha256"]),
            "reservation_nonce_sha256",
            (64,),
        ),
        String(table["native_result_path"]),
        String(table["native_result_sha256"]),
        Float64(table["native_nga50"]),
        Float64(table["comparator_nga50"]),
        Float64(table["native_genome_fraction"]),
        Float64(table["comparator_genome_fraction"]),
        Int(table["native_misassemblies"]),
        Int(table["comparator_misassemblies"]),
        Float64(table["native_error_per_100kbp"]),
        Float64(table["comparator_error_per_100kbp"]),
        String(table["comparator_id"]),
        String(table["comparator_version"]),
        String(table["comparator_result_path"]),
        String(table["comparator_result_sha256"]),
        String(table["comparator_command_sha256"]),
        String(table["comparator_configuration_sha256"]),
        String(table["evaluator_id"]),
        String(table["evaluator_version"]),
        String(table["evaluator_command_sha256"]),
        String(table["development_freeze_sha256"]),
        String(table["manifest_sha256"]),
        String(table["configuration_sha256"]),
        String(table["fixture_sha256"]),
        Dict(String(key) => String(value)
            for (key, value) in table["input_component_paths"]),
        Dict(String(key) => String(value)
            for (key, value) in table["input_component_sha256"]),
        String(table["truth_reference_path"]),
        String(table["truth_reference_sha256"]),
        String(table["metric_report_path"]),
        String(table["metric_report_sha256"]),
        String(table["input_bundle_sha256"]),
        String(table["native_consumed_input_sha256"]),
        String(table["comparator_consumed_input_sha256"]),
        String.(table["comparator_consumed_input_components"]),
        summary,
        Float64.(table["truth_reference_coverages"]),
        Float64(table["abundance_mae"]),
        Float64(table["comparator_abundance_mae"]),
        Int(table["strain_count_mode"]),
        table["strain_count_set_contains_truth"],
        table["primary_rank_agreement"],
        table["ensemble_noninferior_or_abstained"],
    )
end

function stage2_input_component_digest(
        component_sha256::AbstractDict{<:AbstractString, <:AbstractString},
        components::AbstractVector{<:AbstractString}
)::String
    component_names = sort!(unique(String.(components)))
    isempty(component_names) && error("input component set must be nonempty")
    records = String[]
    for component in component_names
        haskey(component_sha256, component) ||
            error("missing digest for input component $(component)")
        digest = _validate_hex_digest(
            component_sha256[component], "input component digest", (64,))
        push!(records, "$(component)=$(digest)")
    end
    return bytes2hex(SHA.sha256(join(records, '\n')))
end

function _observation_relative_path(
        observation_path::AbstractString,
        artifact_path::AbstractString
)::String
    isabspath(artifact_path) &&
        error("artifact paths must be relative to the run bundle")
    relative_path = normpath(String(artifact_path))
    escapes_bundle = relative_path == ".." || startswith(
        relative_path, "..$(Base.Filesystem.path_separator)")
    escapes_bundle &&
        error("artifact path escapes the run bundle")
    return normpath(joinpath(dirname(observation_path), relative_path))
end

function _validate_bundle_file(
        path::AbstractString,
        bundle_directory::AbstractString,
        label::String
)::Nothing
    isfile(path) || error("$(label) artifact does not exist")
    islink(path) && error("$(label) artifact must not be a symbolic link")
    bundle_root = realpath(bundle_directory)
    artifact_realpath = realpath(path)
    relpath(artifact_realpath, bundle_root) == ".." &&
        error("$(label) artifact resolves outside the run bundle")
    startswith(
        relpath(artifact_realpath, bundle_root),
        "..$(Base.Filesystem.path_separator)",
    ) && error("$(label) artifact resolves outside the run bundle")
    return nothing
end

function _validate_stage2_configuration_artifact(
        path::AbstractString,
        bundle_directory::AbstractString,
        expected_sha256::AbstractString
)::String
    _validate_bundle_file(path, bundle_directory, "configuration")
    configuration_bytes = read(path)
    expected_digest = _validate_hex_digest(
        expected_sha256, "configuration_sha256", (64,))
    bytes2hex(SHA.sha256(configuration_bytes)) == expected_digest ||
        error("configuration artifact hash differs from frozen campaign")
    configuration = TOML.parse(String(configuration_bytes))
    get(configuration, "schema", "") == "rhizomorph-stage2-configuration-v1" ||
        error("configuration artifact has an unsupported schema")
    configuration_id = get(configuration, "configuration_id", nothing)
    configuration_id isa AbstractString ||
        error("configuration artifact must define a string configuration_id")
    return _validate_ledger_text(
        configuration_id, "configuration artifact configuration_id")
end

function _required_stage2_modalities(
        manifest::AbstractDict
)::Vector{Symbol}
    modality_contracts = get(manifest, "modalities", nothing)
    modality_contracts isa AbstractDict ||
        error("benchmark manifest must define modality contracts")
    modalities = Symbol[]
    for (modality, contract) in modality_contracts
        contract isa AbstractDict ||
            error("modality $(modality) contract must be a table")
        required = get(contract, "required", nothing)
        required isa Bool ||
            error("modality $(modality) required flag must be Boolean")
        required && push!(modalities, Symbol(modality))
    end
    sort!(modalities)
    isempty(modalities) && error("benchmark manifest has no required modalities")
    unknown = setdiff(Set(modalities), Set(STAGE2_MODALITIES))
    isempty(unknown) ||
        error("benchmark manifest has unknown required modalities: " *
              join(sort!(String.(collect(unknown))), ", "))
    return modalities
end

function _stage2_manifest_seeds(
        manifest::AbstractDict,
        field::String,
        frozen_seeds::Tuple{Vararg{Int}}
)::Vector{Int}
    raw_seeds = get(manifest, field, nothing)
    raw_seeds isa AbstractVector ||
        error("benchmark manifest must define $(field)")
    seeds = Int.(raw_seeds)
    length(unique(seeds)) == length(seeds) ||
        error("$(field) must be unique")
    Tuple(seeds) == frozen_seeds ||
        error("benchmark manifest $(field) differ from the frozen seed partition")
    return seeds
end

function _expected_stage2_development_cell_ids(
        manifest::AbstractDict
)::Set{String}
    gates = get(manifest, "gates", nothing)
    gates isa AbstractDict || error("benchmark manifest must define gates")
    development = get(gates, "development", nothing)
    development isa AbstractDict ||
        error("benchmark manifest must define the development gate")
    for (field, expected) in (
            (
                "evaluation_unit",
                "each-required-modality-each-development-seed",
            ),
            ("comparator_id", "rhizomorph-native"),
            ("cell_id_format", "<modality>:<seed>:rhizomorph-native"),
        )
        get(development, field, nothing) == expected ||
            error("development gate $(field) differs from the frozen contract")
    end
    get(development, "require_exact_grid", false) === true ||
        error("development gate must require the exact frozen grid")
    modalities = _required_stage2_modalities(manifest)
    seeds = _stage2_manifest_seeds(
        manifest, "development_seeds", STAGE2_DEVELOPMENT_SEEDS)
    return Set(
        "$(modality):$(seed):rhizomorph-native"
        for modality in modalities for seed in seeds
    )
end

function _validate_fixture_no_symlink_components(
        path::AbstractString,
        bundle_directory::AbstractString,
        label::String
)::Nothing
    bundle_root = abspath(normpath(bundle_directory))
    current = abspath(normpath(path))
    while current != bundle_root
        islink(current) &&
            error("$(label) artifact path must not contain symbolic links")
        parent = dirname(current)
        parent == current &&
            error("$(label) artifact path is outside the run bundle")
        current = parent
    end
    return nothing
end

function _validate_fixture_artifact(
        fixture_path::AbstractString,
        bundle_directory::AbstractString,
        recorded_path::Any,
        recorded_sha256::Any,
        label::String;
        verify_bytes::Bool = true
)::NamedTuple{(:path, :sha256), Tuple{String, String}}
    recorded_path isa AbstractString ||
        error("$(label) path must be a string")
    recorded_sha256 isa AbstractString ||
        error("$(label) SHA256 must be a string")
    digest = _validate_hex_digest(
        recorded_sha256, "$(label) SHA256", (64,))
    artifact_path = _observation_relative_path(fixture_path, recorded_path)
    bundle_root = abspath(normpath(bundle_directory))
    relative_to_bundle = relpath(abspath(artifact_path), bundle_root)
    escapes_bundle = relative_to_bundle == ".." || startswith(
        relative_to_bundle, "..$(Base.Filesystem.path_separator)")
    escapes_bundle && error("$(label) artifact path is outside the run bundle")
    verify_bytes || return (path = artifact_path, sha256 = digest)
    _validate_bundle_file(artifact_path, bundle_directory, label)
    _validate_fixture_no_symlink_components(
        artifact_path, bundle_directory, label)
    bytes2hex(SHA.sha256(read(artifact_path))) == digest ||
        error("$(label) artifact hash differs from fixture index")
    return (path = artifact_path, sha256 = digest)
end

function _validate_stage2_fixture_index(
        manifest::AbstractDict,
        fixture_path::AbstractString,
        bundle_directory::AbstractString;
        target_cell::Union{Nothing, Tuple{Symbol, Int}} = nothing
)::NamedTuple
    fixture = TOML.parsefile(fixture_path)
    get(fixture, "schema", "") == "rhizomorph-stage2-fixture-index-v1" ||
        error("confirmation fixture index has an unsupported schema")
    modalities = _required_stage2_modalities(manifest)
    seeds = _stage2_manifest_seeds(
        manifest, "confirmation_seeds", STAGE2_CONFIRMATION_SEEDS)
    expected_cell_ids = Set(
        "$(modality):$(seed)" for modality in modalities for seed in seeds)
    if target_cell !== nothing
        "$(target_cell[1]):$(target_cell[2])" in expected_cell_ids ||
            error("confirmation target cell is outside the frozen grid")
    end

    truth_paths = get(fixture, "confirmation_truth_path", nothing)
    truth_sha256 = get(fixture, "confirmation_truth_sha256", nothing)
    truth_paths isa AbstractDict ||
        error("confirmation fixture index must define confirmation_truth_path")
    truth_sha256 isa AbstractDict ||
        error("confirmation fixture index must define confirmation_truth_sha256")
    expected_seed_keys = Set(string.(seeds))
    Set(String.(keys(truth_paths))) == expected_seed_keys ||
        error("confirmation truth paths do not exactly cover frozen seeds")
    Set(String.(keys(truth_sha256))) == expected_seed_keys ||
        error("confirmation truth digests do not exactly cover frozen seeds")
    validated_truth_paths = Dict{Int, String}()
    validated_truth_sha256 = Dict{Int, String}()
    for seed in seeds
        seed_key = string(seed)
        truth = _validate_fixture_artifact(
            fixture_path,
            bundle_directory,
            truth_paths[seed_key],
            truth_sha256[seed_key],
            "confirmation truth $(seed)",
            verify_bytes = target_cell === nothing || seed == target_cell[2],
        )
        validated_truth_paths[seed] = truth.path
        validated_truth_sha256[seed] = truth.sha256
    end

    confirmation_cells = get(fixture, "confirmation_cells", nothing)
    confirmation_cells isa AbstractDict ||
        error("confirmation fixture index must define confirmation_cells")
    Set(String.(keys(confirmation_cells))) == expected_cell_ids ||
        error("confirmation fixture cells do not exactly cover the frozen grid")
    validated_cells = Dict{Tuple{Symbol, Int}, NamedTuple}()
    modality_contracts = manifest["modalities"]
    for modality in modalities, seed in seeds
        modality_key = String(modality)
        cell_id = "$(modality):$(seed)"
        cell = confirmation_cells[cell_id]
        cell isa AbstractDict ||
            error("confirmation fixture cell $(cell_id) must be a table")
        component_paths = get(cell, "input_component_paths", nothing)
        component_sha256 = get(cell, "input_component_sha256", nothing)
        component_paths isa AbstractDict ||
            error("confirmation fixture cell $(cell_id) lacks component paths")
        component_sha256 isa AbstractDict ||
            error("confirmation fixture cell $(cell_id) lacks component digests")
        components = String.(
            modality_contracts[modality_key]["input_components"])
        isempty(components) &&
            error("required modality $(modality) has no input components")
        length(unique(components)) == length(components) ||
            error("required modality $(modality) repeats an input component")
        expected_components = Set(components)
        Set(String.(keys(component_paths))) == expected_components ||
            error("confirmation fixture cell $(cell_id) component paths differ " *
                  "from the manifest")
        Set(String.(keys(component_sha256))) == expected_components ||
            error("confirmation fixture cell $(cell_id) component digests differ " *
                  "from the manifest")
        validated_component_paths = Dict{String, String}()
        validated_component_sha256 = Dict{String, String}()
        for component in sort(components)
            artifact = _validate_fixture_artifact(
                fixture_path,
                bundle_directory,
                component_paths[component],
                component_sha256[component],
                "confirmation cell $(cell_id) component $(component)",
                verify_bytes = target_cell === nothing ||
                               (modality, seed) == target_cell,
            )
            validated_component_paths[component] = artifact.path
            validated_component_sha256[component] = artifact.sha256
        end
        input_bundle_sha256 = _validate_hex_digest(
            get(cell, "input_bundle_sha256", ""),
            "confirmation cell $(cell_id) input bundle SHA256",
            (64,),
        )
        stage2_input_component_digest(
            validated_component_sha256, components) == input_bundle_sha256 ||
            error("confirmation fixture cell $(cell_id) bundle digest differs " *
                  "from its components")
        error_profile = _validate_fixture_artifact(
            fixture_path,
            bundle_directory,
            get(cell, "error_profile_path", nothing),
            get(cell, "error_profile_sha256", nothing),
            "confirmation cell $(cell_id) error profile",
            verify_bytes = target_cell === nothing ||
                           (modality, seed) == target_cell,
        )
        validated_cells[(modality, seed)] = (;
            input_component_paths = validated_component_paths,
            input_component_sha256 = validated_component_sha256,
            input_bundle_sha256,
            error_profile_path = error_profile.path,
            error_profile_sha256 = error_profile.sha256,
        )
    end
    return (;
        truth_paths = validated_truth_paths,
        truth_sha256 = validated_truth_sha256,
        confirmation_cells = validated_cells,
    )
end

function _truth_reference_ids(path::AbstractString)::Vector{String}
    reader = FASTX.FASTA.Reader(open(path, "r"))
    try
        return String[String(FASTX.identifier(record)) for record in reader]
    finally
        close(reader)
    end
end

function _validate_observation_raw_artifacts(
        observation::Stage2ReleaseObservation,
        observation_path::AbstractString
)::Nothing
    for (recorded_path, expected_sha256, label) in (
            (
                observation.native_result_path,
                observation.native_result_sha256,
                "native result",
            ),
            (
                observation.comparator_result_path,
                observation.comparator_result_sha256,
                "comparator result",
            ),
            (
                observation.truth_reference_path,
                observation.truth_reference_sha256,
                "truth reference",
            ),
        )
        artifact_path = _observation_relative_path(
            observation_path, recorded_path)
        _validate_bundle_file(
            artifact_path, dirname(observation_path), label)
        bytes2hex(SHA.sha256(read(artifact_path))) == expected_sha256 ||
            error("$(label) artifact hash differs from observation")
    end
    truth_path = _observation_relative_path(
        observation_path, observation.truth_reference_path)
    truth_ids = _truth_reference_ids(truth_path)
    truth_ids == observation.assignment.expected_truth_ids ||
        error("truth FASTA identifiers differ from expected truth IDs")
    length(unique(truth_ids)) == length(truth_ids) ||
        error("truth FASTA identifiers must be unique")
    Set(keys(observation.input_component_paths)) ==
    Set(keys(observation.input_component_sha256)) ||
        error("input component paths and digests have different keys")
    for component in keys(observation.input_component_paths)
        component_path = _observation_relative_path(
            observation_path, observation.input_component_paths[component])
        _validate_bundle_file(
            component_path,
            dirname(observation_path),
            "input component $(component)",
        )
        bytes2hex(SHA.sha256(read(component_path))) ==
        observation.input_component_sha256[component] ||
            error("input component artifact hash differs: $(component)")
    end
    return nothing
end

function _validate_observation_metric_report(
        observation::Stage2ReleaseObservation,
        observation_path::AbstractString,
        evaluation_contract::AbstractDict
)::Nothing
    evaluation_contract["version"] != "UNRESOLVED" ||
        error("evaluation tool version is unresolved")
    observation.evaluator_id == evaluation_contract["tool"] ||
        error("metric report evaluator differs from manifest")
    observation.evaluator_version == evaluation_contract["version"] ||
        error("metric report evaluator version differs from manifest")
    expected_command_sha256 = bytes2hex(SHA.sha256(
        evaluation_contract["command"]))
    observation.evaluator_command_sha256 == expected_command_sha256 ||
        error("metric evaluator command hash differs from manifest")
    report_path = _observation_relative_path(
        observation_path, observation.metric_report_path)
    _validate_bundle_file(
        report_path, dirname(observation_path), "metric report")
    report_bytes = read(report_path)
    bytes2hex(SHA.sha256(report_bytes)) == observation.metric_report_sha256 ||
        error("metric report artifact hash differs from observation")
    report = TOML.parse(String(report_bytes))
    report["schema"] == evaluation_contract["report_schema"] ||
        error("metric report schema differs from manifest")
    for (field, expected) in (
            ("evaluator_id", observation.evaluator_id),
            ("evaluator_version", observation.evaluator_version),
            ("evaluator_command_sha256", observation.evaluator_command_sha256),
            ("reservation_nonce_sha256",
                observation.reservation_nonce_sha256),
            ("truth_reference_sha256", observation.truth_reference_sha256),
            ("native_result_sha256", observation.native_result_sha256),
            ("comparator_result_sha256", observation.comparator_result_sha256),
            ("native_nga50", observation.native_nga50),
            ("comparator_nga50", observation.comparator_nga50),
            ("native_genome_fraction", observation.native_genome_fraction),
            ("comparator_genome_fraction", observation.comparator_genome_fraction),
            ("native_misassemblies", observation.native_misassemblies),
            ("comparator_misassemblies", observation.comparator_misassemblies),
            ("native_error_per_100kbp", observation.native_error_per_100kbp),
            ("comparator_error_per_100kbp",
                observation.comparator_error_per_100kbp),
            ("abundance_mae", observation.abundance_mae),
            ("comparator_abundance_mae", observation.comparator_abundance_mae),
            ("strain_count_mode", observation.strain_count_mode),
            ("strain_count_set_contains_truth",
                observation.strain_count_set_contains_truth),
            ("primary_rank_agreement", observation.primary_rank_agreement),
            ("ensemble_noninferior_or_abstained",
                observation.ensemble_noninferior_or_abstained),
        )
        isequal(report[field], expected) ||
            error("metric report field $(field) differs from observation")
    end
    report_assignment = report["assignment"]
    reported_summary = CandidateAssignmentSummary(
        String.(report_assignment["expected_truth_ids"]),
        String.(report_assignment["assigned_truth_ids"]),
        String.(report_assignment["unique_assigned_truth_ids"]),
        String.(report_assignment["duplicate_truth_ids"]),
        String.(report_assignment["missing_truth_ids"]),
        String.(report_assignment["unexpected_truth_ids"]),
        report_assignment["exact_recovery"],
    )
    isequal(
        _assignment_signature(reported_summary),
        _assignment_signature(observation.assignment),
    ) || error("metric report assignment differs from observation")
    isequal(
        Float64.(report["truth_reference_coverages"]),
        observation.truth_reference_coverages,
    ) || error("metric report truth coverages differ from observation")
    return nothing
end

"""
Evaluate the tight practical release gate without pooling modalities or seeds.

Every required modality/seed cell must be present exactly once. Structural,
haplotype, abundance, count, primary, and ensemble gates apply to every cell;
accuracy wins and the median reduction are evaluated independently within each
modality across its three confirmation seeds.
"""
function _evaluate_stage2_release_observations(
        observations::AbstractVector{Stage2ReleaseObservation};
        manifest::AbstractDict,
        development_freeze_sha256::String,
        confirmation_ledger_sha256::String,
        attempt_result_sha256::Dict{String, String},
        diagnostic_attempt_status::Dict{Tuple{Symbol, Int, String}, Symbol},
        diagnostic_result_sha256::Dict{String, String}
)::Stage2ReleaseGateResult
    release_gate = manifest["gates"]["release"]
    release_gate["evaluation_unit"] == "each-modality-each-confirmation-seed" ||
        error("release manifest has an unsupported evaluation unit")
    release_gate["baseline"] ==
    "best-compatible-locked-comparator-within-modality-and-input-stratum" ||
        error("release manifest has an unsupported comparator baseline")
    release_gate["required_input_stratum"] == "full" ||
        error("release manifest must gate the full-input stratum")
    get(release_gate, "subset_input_strata_are_gating", true) == false ||
        error("subset input strata must be non-gating")
    for policy in (
            "require_three_distinct_truths",
            "require_exactly_three_supported_haplotypes",
            "require_no_spurious_or_duplicate_haplotypes",
            "require_90pct_strain_count_set_contains_truth",
            "require_primary_rank_agreement",
            "require_every_modality",
        )
        get(release_gate, policy, false) ||
            error("release manifest must enable $(policy)")
    end
    release_gate["required_strain_count_mode"] == 3 ||
        error("release manifest must require strain-count mode three")
    release_gate["abundance_baseline_relation"] ==
    "mae-no-worse-than-best-compatible-comparator" ||
        error("release manifest has an unsupported abundance relation")
    release_gate["ensemble_policy"] == "noninferior-to-native-or-abstain" ||
        error("release manifest has an unsupported ensemble policy")
    modality_contracts = manifest["modalities"]
    comparator_contracts = manifest["comparators"]
    required_modalities = Tuple(sort!(Symbol[
        Symbol(modality) for (modality, contract) in modality_contracts
        if get(contract, "required", false)
    ]))
    required_seeds = Tuple(Int.(manifest["confirmation_seeds"]))
    length(unique(required_seeds)) == length(required_seeds) ||
        error("confirmation seeds must be unique")
    required_seeds == STAGE2_CONFIRMATION_SEEDS ||
        error("release manifest confirmation seeds differ from the frozen holdout")
    if get(release_gate, "require_every_modality", false)
        Set(required_modalities) == Set(STAGE2_MODALITIES) ||
            error("release manifest does not require every frozen modality")
    end
    expected_keys = Set((modality, seed)
        for modality in required_modalities for seed in required_seeds)
    observed_keys = Set((observation.modality, observation.seed)
        for observation in observations)
    observed_keys == expected_keys ||
        error("release observations do not exactly cover required modality/seed cells")

    observation_gates = Dict{Tuple{Symbol, Int}, Dict{String, Bool}}()
    reductions_by_key = Dict{Tuple{Symbol, Int}, Float64}()
    for key in sort!(collect(expected_keys))
        modality, _ = key
        modality_key = String(modality)
        modality_contract = modality_contracts[modality_key]
        modality_components = Set(String.(modality_contract["input_components"]))
        required_comparators = Set(String(comparator) for comparator in
            modality_contract["comparators"] if Set(String.(
                comparator_contracts[String(comparator)][
                    "consumed_input_components"][modality_key])) ==
            modality_components)
        rows = filter(observation ->
            (observation.modality, observation.seed) == key, observations)
        length(rows) == length(required_comparators) ||
            error("release cell $(key) lacks a complete comparator panel")
        length(unique(row.comparator_id for row in rows)) == length(rows) ||
            error("release cell $(key) contains duplicate comparator rows")
        Set(row.comparator_id for row in rows) == required_comparators ||
            error("release cell $(key) comparator panel differs from manifest")

        first_row = first(rows)
        _validate_assignment_summary(first_row.assignment)
        native_signature = _native_scientific_signature(first_row)
        isempty(modality_components) &&
            error("release modality $(modality) has no input components")
        full_input_rows = Stage2ReleaseObservation[]
        for row in rows
            values_to_check = (
                row.native_nga50,
                row.comparator_nga50,
                row.native_genome_fraction,
                row.comparator_genome_fraction,
                row.native_error_per_100kbp,
                row.comparator_error_per_100kbp,
                row.abundance_mae,
                row.comparator_abundance_mae,
                row.truth_reference_coverages...,
            )
            all(isfinite, values_to_check) || error("release metrics must be finite")
            row.native_nga50 >= 0.0 && row.comparator_nga50 >= 0.0 ||
                error("NGA50 values must be nonnegative")
            0.0 <= row.native_genome_fraction <= 100.0 &&
                0.0 <= row.comparator_genome_fraction <= 100.0 ||
                error("genome fractions must be between zero and 100")
            row.native_misassemblies >= 0 && row.comparator_misassemblies >= 0 ||
                error("misassembly counts must be nonnegative")
            row.native_error_per_100kbp >= 0.0 &&
                row.comparator_error_per_100kbp > 0.0 ||
                error("error rates must be nonnegative and comparator error positive")
            0.0 <= row.abundance_mae <= 1.0 &&
                0.0 <= row.comparator_abundance_mae <= 1.0 ||
                error("abundance MAE must be between zero and one")
            all(coverage -> 0.0 <= coverage <= 1.0,
                row.truth_reference_coverages) ||
                error("truth reference coverage must be between zero and one")
            row.strain_count_mode >= 0 || error("strain count mode must be nonnegative")
            length(row.truth_reference_coverages) ==
            length(row.assignment.expected_truth_ids) ||
                error("truth coverage count differs from expected truth count")
            _validate_assignment_summary(row.assignment)
            isequal(_native_scientific_signature(row), native_signature) ||
                error("native evidence disagrees across comparator rows")
            Set(keys(row.input_component_sha256)) == modality_components ||
                error("observation does not bind every modality input component")
            stage2_input_component_digest(
                row.input_component_sha256,
                collect(modality_components),
            ) == row.input_bundle_sha256 ||
                error("modality bundle hash differs from component digests")
            row.native_consumed_input_sha256 == row.input_bundle_sha256 ||
                error("native Stage-2 did not consume the full modality bundle")

            comparator_contract = comparator_contracts[row.comparator_id]
            comparator_contract["version"] != "UNRESOLVED" ||
                error("comparator $(row.comparator_id) version is unresolved")
            comparator_contract["version"] == row.comparator_version ||
                error("comparator version differs from locked manifest")
            occursin("documented", lowercase(comparator_contract["command"])) &&
                error("comparator command is not executable")
            expected_command_hash = bytes2hex(SHA.sha256(
                comparator_contract["command"]))
            expected_configuration_hash = bytes2hex(SHA.sha256(
                comparator_contract["configuration"]))
            row.comparator_command_sha256 == expected_command_hash ||
                error("comparator command hash differs from manifest")
            row.comparator_configuration_sha256 == expected_configuration_hash ||
                error("comparator configuration hash differs from manifest")
            _validate_hex_digest(
                row.comparator_result_sha256, "comparator_result_sha256", (64,))
            _validate_hex_digest(
                row.input_bundle_sha256, "input_bundle_sha256", (64,))
            _validate_hex_digest(
                row.comparator_consumed_input_sha256,
                "comparator_consumed_input_sha256",
                (64,),
            )
            expected_components = Set(String.(
                comparator_contract["consumed_input_components"][modality_key]))
            observed_components = Set(row.comparator_consumed_input_components)
            observed_components == expected_components ||
                error("comparator consumed-input components differ from manifest")
            isempty(observed_components) &&
                error("comparator consumed-input components must be nonempty")
            issubset(observed_components, modality_components) ||
                error("comparator consumes input outside the frozen modality bundle")
            stage2_input_component_digest(
                row.input_component_sha256,
                row.comparator_consumed_input_components,
            ) == row.comparator_consumed_input_sha256 ||
                error("comparator consumed-input hash differs from component digests")
            if observed_components == modality_components
                row.comparator_consumed_input_sha256 ==
                row.input_bundle_sha256 ||
                    error("full-input comparator hash differs from modality bundle")
                push!(full_input_rows, row)
            else
                row.comparator_consumed_input_sha256 !=
                row.input_bundle_sha256 ||
                    error("subset comparator is mislabeled as a full-input consumer")
            end
        end
        length(unique(row.input_bundle_sha256 for row in rows)) == 1 ||
            error("release cell comparator rows use different modality bundles")
        isempty(full_input_rows) &&
            error("release cell $(key) has no full-input comparator")

        comparator_nga50 = maximum(row.comparator_nga50 for row in full_input_rows)
        comparator_genome_fraction = maximum(
            row.comparator_genome_fraction for row in full_input_rows)
        comparator_misassemblies = minimum(
            row.comparator_misassemblies for row in full_input_rows)
        comparator_error = minimum(
            row.comparator_error_per_100kbp for row in full_input_rows)
        comparator_abundance_mae = minimum(
            row.comparator_abundance_mae for row in full_input_rows)
        reductions_by_key[key] =
            (comparator_error - first_row.native_error_per_100kbp) / comparator_error
        observation_gates[key] = Dict(
            "misassemblies" => first_row.native_misassemblies <=
                               comparator_misassemblies +
                               release_gate["maximum_extra_misassemblies"],
            "nga50" => first_row.native_nga50 >=
                       release_gate["minimum_nga50_ratio"] *
                       comparator_nga50,
            "genome_fraction" => first_row.native_genome_fraction >=
                                 comparator_genome_fraction -
                                 release_gate[
                "maximum_genome_fraction_deficit_points"],
            "three_truths" => first_row.assignment.exact_recovery &&
                              length(first_row.assignment.assigned_truth_ids) == 3,
            "truth_coverage" => all(coverage -> coverage >=
                                    release_gate["minimum_truth_reference_coverage"],
                first_row.truth_reference_coverages),
            "abundance_absolute" => first_row.abundance_mae <=
                                    release_gate["maximum_abundance_mae"],
            "abundance_comparator" => first_row.abundance_mae <=
                                      comparator_abundance_mae,
            "strain_count" => first_row.strain_count_mode == 3 &&
                              first_row.strain_count_set_contains_truth,
            "primary" => first_row.primary_rank_agreement,
            "ensemble" => first_row.ensemble_noninferior_or_abstained,
            "no_large_error_regression" => first_row.native_error_per_100kbp <=
                                           (1.0 + release_gate[
                "maximum_single_seed_error_regression_fraction"]) * comparator_error,
        )
    end

    modality_gates = Dict{Symbol, Dict{String, Bool}}()
    for modality in required_modalities
        reductions = [reductions_by_key[(modality, seed)] for seed in required_seeds]
        modality_gates[modality] = Dict(
            "accuracy_wins" => count(>(0.0), reductions) >=
                               release_gate[
                "minimum_accuracy_win_seeds_per_modality"],
            "median_error_reduction" => Statistics.median(reductions) >=
                                        release_gate[
                "minimum_median_error_reduction_fraction"],
        )
    end
    passed = all(all(values(gates)) for gates in values(observation_gates)) &&
             all(all(values(gates)) for gates in values(modality_gates))
    return Stage2ReleaseGateResult(
        passed,
        observation_gates,
        modality_gates,
        development_freeze_sha256,
        confirmation_ledger_sha256,
        attempt_result_sha256,
        diagnostic_attempt_status,
        diagnostic_result_sha256,
    )
end

function _validate_ledger_text(value::AbstractString, field::String)::String
    normalized = String(value)
    isempty(normalized) && error("$(field) must be nonempty")
    any(character -> character == '\t' || character == '\n' || character == '\r',
        normalized) && error("$(field) must not contain tabs or newlines")
    return normalized
end

function _validate_empty_text(value::AbstractString, field::String)::String
    isempty(value) || error("$(field) must be empty for development attempts")
    return String(value)
end

function _validate_utc_timestamp(value::AbstractString, field::String)::String
    timestamp = _validate_ledger_text(value, field)
    endswith(timestamp, "Z") ||
        error("$(field) must use UTC format yyyy-mm-ddTHH:MM:SSZ")
    try
        Dates.DateTime(
            chop(timestamp; tail = 1), Dates.dateformat"yyyy-mm-ddTHH:MM:SS")
    catch
        error("$(field) must use UTC format yyyy-mm-ddTHH:MM:SSZ")
    end
    return timestamp
end

function _validate_hex_digest(
        value::AbstractString,
        field::String,
        lengths::Tuple{Vararg{Int}}
)::String
    digest = lowercase(_validate_ledger_text(value, field))
    length(digest) in lengths ||
        error("$(field) must have length $(join(lengths, " or "))")
    all(isxdigit, digest) || error("$(field) must be hexadecimal")
    return digest
end

function Stage2DevelopmentFreeze(;
        frozen_at_utc::AbstractString,
        development_ledger_path::AbstractString,
        confirmation_ledger_path::AbstractString,
        manifest_path::AbstractString,
        configuration_path::AbstractString,
        fixture_path::AbstractString,
        git_sha::AbstractString,
        manifest_sha256::AbstractString,
        configuration_sha256::AbstractString,
        fixture_sha256::AbstractString,
        development_ledger_sha256::AbstractString,
        development_ledger_rows::Integer,
        final_attempt_id::AbstractString,
        git_dirty::Bool
)::Stage2DevelopmentFreeze
    return Stage2DevelopmentFreeze(
        String(frozen_at_utc),
        String(development_ledger_path),
        String(confirmation_ledger_path),
        String(manifest_path),
        String(configuration_path),
        String(fixture_path),
        String(git_sha),
        String(manifest_sha256),
        String(configuration_sha256),
        String(fixture_sha256),
        String(development_ledger_sha256),
        Int(development_ledger_rows),
        String(final_attempt_id),
        git_dirty,
    )
end

function stage2_development_freeze_path(
        development_ledger_path::AbstractString
)::String
    return abspath(normpath(String(development_ledger_path) * ".freeze.toml"))
end

function stage2_confirmation_ledger_path(
        development_ledger_path::AbstractString
)::String
    return abspath(normpath(String(development_ledger_path) * ".confirmation.tsv"))
end

function _validate_development_gate_artifact(
        development_ledger_path::AbstractString,
        final_row::Stage2AttemptLedgerRow,
        manifest_path::AbstractString
)::Nothing
    final_row.status == :passed ||
        error("development ledger final attempt must pass the strict gate")
    result_path = _stage2_result_artifact_path(
        development_ledger_path, final_row.result_path)
    _validate_bundle_file(
        result_path, dirname(development_ledger_path), "development gate")
    result_bytes = read(result_path)
    bytes2hex(SHA.sha256(result_bytes)) == final_row.result_sha256 ||
        error("development gate artifact hash differs from final ledger row")
    gate = TOML.parse(String(result_bytes))
    get(gate, "schema", "") == "rhizomorph-stage2-development-gate-v1" ||
        error("development gate artifact has an unsupported schema")
    get(gate, "gate_id", "") == "strict-development-grid" ||
        error("development gate artifact has an unsupported gate ID")
    decisions = get(gate, "cell_decisions", Any[])
    decisions isa AbstractVector ||
        error("development gate cell decisions must be an array")
    isempty(decisions) &&
        error("development gate artifact has no cell decisions")
    manifest = TOML.parsefile(manifest_path)
    expected_cell_ids = _expected_stage2_development_cell_ids(manifest)
    observed_cell_ids = String[]
    for decision in decisions
        decision isa AbstractDict ||
            error("development gate cell decisions must be tables")
        cell_id = get(decision, "cell_id", nothing)
        cell_id isa AbstractString ||
            error("development cell decision lacks a cell ID")
        isempty(cell_id) &&
            error("development cell decision lacks a cell ID")
        decision_passed = get(decision, "passed", nothing)
        decision_passed isa Bool ||
            error("development cell decision passed value must be Boolean")
        push!(observed_cell_ids, String(cell_id))
    end
    length(unique(observed_cell_ids)) == length(observed_cell_ids) ||
        error("development gate artifact contains duplicate cell IDs")
    Set(observed_cell_ids) == expected_cell_ids ||
        error("development gate artifact does not exactly cover the frozen grid")
    recomputed_pass = all(decision["passed"] for decision in decisions)
    strict_gate_passed = get(gate, "strict_gate_passed", nothing)
    strict_gate_passed isa Bool ||
        error("development gate aggregate must be Boolean")
    strict_gate_passed === recomputed_pass ||
        error("development gate aggregate differs from cell decisions")
    recomputed_pass ||
        error("development artifact does not attest a strict-gate PASS")
    for (field, expected) in (
            ("git_sha", final_row.git_sha),
            ("manifest_sha256", final_row.manifest_sha256),
            ("configuration_sha256", final_row.configuration_sha256),
            ("fixture_sha256", final_row.fixture_sha256),
        )
        get(gate, field, "") == expected ||
            error("development gate artifact $(field) differs from final ledger row")
    end
    for decision in decisions
        cell_id = String(decision["cell_id"])
        decision_passed = get(decision, "passed", false)
        decision_passed isa Bool ||
            error("development cell decision passed value must be Boolean")
        source_path = get(decision, "source_report_path", "")
        isempty(source_path) &&
            error("development cell decision lacks a source report")
        source_sha256 = _validate_hex_digest(
            get(decision, "source_report_sha256", ""),
            "development source report SHA256",
            (64,),
        )
        absolute_source_path = _observation_relative_path(result_path, source_path)
        _validate_bundle_file(
            absolute_source_path,
            dirname(result_path),
            "development source report",
        )
        source_bytes = read(absolute_source_path)
        bytes2hex(SHA.sha256(source_bytes)) == source_sha256 ||
            error("development source report hash differs from gate artifact")
        source_report = TOML.parse(String(source_bytes))
        get(source_report, "schema", "") ==
        "rhizomorph-stage2-development-cell-v1" ||
            error("development source report has an unsupported schema")
        get(source_report, "cell_id", "") == cell_id ||
            error("development source report cell differs from gate decision")
        source_passed = get(source_report, "passed", nothing)
        source_passed isa Bool ||
            error("development source report passed value must be Boolean")
        source_passed === decision_passed ||
            error("development source report decision differs from gate decision")
        for (field, expected) in (
                ("git_sha", final_row.git_sha),
                ("manifest_sha256", final_row.manifest_sha256),
                ("configuration_sha256", final_row.configuration_sha256),
                ("fixture_sha256", final_row.fixture_sha256),
            )
            get(source_report, field, nothing) == expected ||
                error("development source report $(field) differs from campaign")
        end
    end
    return nothing
end

function seal_stage2_development_ledger!(
        development_ledger_path::AbstractString;
        frozen_at_utc::AbstractString,
        manifest_path::AbstractString,
        configuration_path::AbstractString,
        fixture_path::AbstractString,
        git_sha::AbstractString,
        manifest_sha256::AbstractString,
        configuration_sha256::AbstractString,
        fixture_sha256::AbstractString,
        git_dirty::Bool
)::Stage2DevelopmentSeal
    ledger_path = abspath(normpath(String(development_ledger_path)))
    canonical_manifest_path = abspath(normpath(String(manifest_path)))
    canonical_configuration_path = abspath(normpath(String(configuration_path)))
    canonical_fixture_path = abspath(normpath(String(fixture_path)))
    output = stage2_development_freeze_path(ledger_path)
    lock_path = ledger_path * ".lock"
    try
        mkdir(lock_path)
    catch
        error("development ledger is locked: $(ledger_path)")
    end
    try
        for (artifact_path, expected_sha256, label) in (
                (canonical_manifest_path, manifest_sha256, "manifest"),
                (canonical_configuration_path, configuration_sha256, "configuration"),
                (canonical_fixture_path, fixture_sha256, "fixture"),
            )
            _validate_bundle_file(
                artifact_path, dirname(ledger_path), label)
            bytes2hex(SHA.sha256(read(artifact_path))) == expected_sha256 ||
                error("$(label) artifact hash differs from seal request")
        end
        frozen_configuration_id = _validate_stage2_configuration_artifact(
            canonical_configuration_path,
            dirname(ledger_path),
            configuration_sha256,
        )
        manifest = TOML.parsefile(canonical_manifest_path)
        _validate_stage2_fixture_index(
            manifest,
            canonical_fixture_path,
            dirname(ledger_path),
        )
        fingerprint = stage2_attempt_ledger_fingerprint(ledger_path)
        ispath(output) && error("development freeze already exists: $(output)")
        freeze = Stage2DevelopmentFreeze(;
            frozen_at_utc,
            development_ledger_path = ledger_path,
            confirmation_ledger_path = stage2_confirmation_ledger_path(
                ledger_path),
            manifest_path = canonical_manifest_path,
            configuration_path = canonical_configuration_path,
            fixture_path = canonical_fixture_path,
            git_sha,
            manifest_sha256,
            configuration_sha256,
            fixture_sha256,
            development_ledger_sha256 = fingerprint.sha256,
            development_ledger_rows = fingerprint.rows,
            final_attempt_id = fingerprint.final_attempt_id,
            git_dirty,
        )
        development_rows = load_stage2_attempt_ledger(ledger_path)
        final_row = last(development_rows)
        final_row.attempt_id == freeze.final_attempt_id ||
            error("development ledger final row differs from fingerprint")
        final_row.git_dirty &&
            error("development ledger final attempt must record git_dirty=false")
        final_row.git_sha == freeze.git_sha ||
            error("development ledger final git SHA differs from freeze")
        final_row.manifest_sha256 == freeze.manifest_sha256 ||
            error("development ledger final manifest hash differs from freeze")
        final_row.configuration_sha256 == freeze.configuration_sha256 ||
            error("development ledger final configuration hash differs from freeze")
        final_row.configuration_id == frozen_configuration_id ||
            error("development ledger final configuration ID differs from artifact")
        final_row.fixture_sha256 == freeze.fixture_sha256 ||
            error("development ledger final fixture hash differs from freeze")
        _validate_development_gate_artifact(
            ledger_path, final_row, canonical_manifest_path)
        mktemp(dirname(output)) do temporary_path, stream
            TOML.print(stream, Dict(
                "frozen_at_utc" => freeze.frozen_at_utc,
                "development_ledger_path" => freeze.development_ledger_path,
                "confirmation_ledger_path" => freeze.confirmation_ledger_path,
                "manifest_path" => freeze.manifest_path,
                "configuration_path" => freeze.configuration_path,
                "fixture_path" => freeze.fixture_path,
                "git_sha" => freeze.git_sha,
                "git_dirty" => freeze.git_dirty,
                "manifest_sha256" => freeze.manifest_sha256,
                "configuration_sha256" => freeze.configuration_sha256,
                "fixture_sha256" => freeze.fixture_sha256,
                "development_ledger_sha256" => freeze.development_ledger_sha256,
                "development_ledger_rows" => freeze.development_ledger_rows,
                "final_attempt_id" => freeze.final_attempt_id,
            ))
            close(stream)
            mv(temporary_path, output)
        end
        return Stage2DevelopmentSeal(
            freeze,
            output,
            bytes2hex(SHA.sha256(read(output))),
        )
    finally
        rm(lock_path; force = true, recursive = true)
    end
end

function load_stage2_development_freeze(
        path::AbstractString
)::Stage2DevelopmentFreeze
    isfile(path) || error("development freeze does not exist: $(path)")
    table = TOML.parsefile(path)
    get(table, "git_dirty", true) == false ||
        error("development freeze must record git_dirty=false")
    return Stage2DevelopmentFreeze(;
        frozen_at_utc = table["frozen_at_utc"],
        development_ledger_path = table["development_ledger_path"],
        confirmation_ledger_path = table["confirmation_ledger_path"],
        manifest_path = table["manifest_path"],
        configuration_path = table["configuration_path"],
        fixture_path = table["fixture_path"],
        git_sha = table["git_sha"],
        manifest_sha256 = table["manifest_sha256"],
        configuration_sha256 = table["configuration_sha256"],
        fixture_sha256 = table["fixture_sha256"],
        development_ledger_sha256 = table["development_ledger_sha256"],
        development_ledger_rows = table["development_ledger_rows"],
        final_attempt_id = table["final_attempt_id"],
        git_dirty = table["git_dirty"],
    )
end

function Stage2AttemptLedgerRow(;
        attempt_id::AbstractString,
        recorded_at_utc::AbstractString,
        seed_partition::Symbol,
        seed::Integer,
        modality::Symbol,
        configuration_id::AbstractString,
        git_sha::AbstractString,
        git_dirty::Bool,
        manifest_sha256::AbstractString,
        configuration_sha256::AbstractString,
        fixture_sha256::AbstractString,
        truth_reference_sha256::AbstractString = "",
        input_bundle_sha256::AbstractString,
        consumed_input_sha256::AbstractString,
        error_profile_sha256::AbstractString,
        comparator_id::AbstractString,
        comparator_version::AbstractString,
        command_sha256::AbstractString,
        comparator_configuration_sha256::AbstractString = "",
        status::Symbol,
        result_path::AbstractString,
        result_sha256::AbstractString,
        development_freeze_sha256::AbstractString = "",
        notes::AbstractString,
        development_freeze::Union{Nothing, Stage2DevelopmentFreeze} = nothing
)::Stage2AttemptLedgerRow
    return Stage2AttemptLedgerRow(
        String(attempt_id),
        String(recorded_at_utc),
        seed_partition,
        Int(seed),
        modality,
        String(configuration_id),
        String(git_sha),
        git_dirty,
        String(manifest_sha256),
        String(configuration_sha256),
        String(fixture_sha256),
        String(truth_reference_sha256),
        String(input_bundle_sha256),
        String(consumed_input_sha256),
        String(error_profile_sha256),
        String(comparator_id),
        String(comparator_version),
        String(command_sha256),
        String(comparator_configuration_sha256),
        status,
        String(result_path),
        String(result_sha256),
        String(development_freeze_sha256),
        String(notes),
        development_freeze,
    )
end

function _attempt_ledger_values(
        row::Stage2AttemptLedgerRow
)::NTuple{24, String}
    return (
        row.attempt_id,
        row.recorded_at_utc,
        String(row.seed_partition),
        string(row.seed),
        String(row.modality),
        row.configuration_id,
        row.git_sha,
        string(row.git_dirty),
        row.manifest_sha256,
        row.configuration_sha256,
        row.fixture_sha256,
        row.truth_reference_sha256,
        row.input_bundle_sha256,
        row.consumed_input_sha256,
        row.error_profile_sha256,
        row.comparator_id,
        row.comparator_version,
        row.command_sha256,
        row.comparator_configuration_sha256,
        String(row.status),
        row.result_path,
        row.result_sha256,
        row.development_freeze_sha256,
        row.notes,
    )
end

function _validate_confirmation_freeze_locked(
        row::Stage2AttemptLedgerRow,
        development_ledger_path::AbstractString
)::Stage2DevelopmentFreeze
    row.seed_partition == :confirmation ||
        error("confirmation authorization requires a confirmation row")
    freeze_path = stage2_development_freeze_path(development_ledger_path)
    freeze = load_stage2_development_freeze(freeze_path)
    abspath(normpath(String(development_ledger_path))) ==
    freeze.development_ledger_path ||
        error("development ledger path differs from frozen campaign identity")
    stage2_confirmation_ledger_path(development_ledger_path) ==
    freeze.confirmation_ledger_path ||
        error("confirmation ledger path differs from frozen campaign identity")
    fingerprint = stage2_attempt_ledger_fingerprint(development_ledger_path)
    fingerprint.sha256 == freeze.development_ledger_sha256 ||
        error("development ledger hash differs from freeze")
    fingerprint.rows == freeze.development_ledger_rows ||
        error("development ledger row count differs from freeze")
    fingerprint.final_attempt_id == freeze.final_attempt_id ||
        error("development ledger final attempt differs from freeze")
    row.git_dirty && error("confirmation attempts require git_dirty=false")
    row.git_sha == freeze.git_sha || error("confirmation git SHA differs from freeze")
    row.manifest_sha256 == freeze.manifest_sha256 ||
        error("confirmation manifest hash differs from freeze")
    row.configuration_sha256 == freeze.configuration_sha256 ||
        error("confirmation configuration hash differs from freeze")
    row.fixture_sha256 == freeze.fixture_sha256 ||
        error("confirmation fixture hash differs from freeze")
    freeze_sha256 = bytes2hex(SHA.sha256(read(freeze_path)))
    row.development_freeze_sha256 == freeze_sha256 ||
        error("confirmation development freeze hash differs from canonical freeze")
    return freeze
end

function _validate_confirmation_reservation_preflight(
        row::Stage2AttemptLedgerRow,
        freeze::Stage2DevelopmentFreeze
)::Nothing
    row.status == :reserved || return nothing
    for (artifact_path, expected_sha256, label) in (
            (freeze.manifest_path, freeze.manifest_sha256, "manifest"),
            (freeze.configuration_path,
                freeze.configuration_sha256, "configuration"),
            (freeze.fixture_path, freeze.fixture_sha256, "fixture"),
        )
        _validate_bundle_file(
            artifact_path, dirname(freeze.development_ledger_path), label)
        bytes2hex(SHA.sha256(read(artifact_path))) == expected_sha256 ||
            error("confirmation preflight found changed $(label) bytes")
    end
    frozen_configuration_id = _validate_stage2_configuration_artifact(
        freeze.configuration_path,
        dirname(freeze.development_ledger_path),
        freeze.configuration_sha256,
    )
    row.configuration_id == frozen_configuration_id ||
        error("confirmation configuration ID differs from frozen artifact")
    manifest = TOML.parsefile(freeze.manifest_path)
    manifest["evaluation"]["version"] != "UNRESOLVED" ||
        error("confirmation preflight requires a resolved evaluator version")
    modality_key = String(row.modality)
    haskey(manifest["modalities"], modality_key) ||
        error("confirmation preflight found an unregistered modality")
    modality = manifest["modalities"][modality_key]
    row.comparator_id in String.(modality["comparators"]) ||
        error("confirmation preflight found an unregistered comparator cell")
    comparator = manifest["comparators"][row.comparator_id]
    comparator["version"] != "UNRESOLVED" ||
        error("confirmation preflight requires a resolved comparator version")
    row.comparator_version == comparator["version"] ||
        error("confirmation comparator version differs from frozen manifest")
    comparator_command = get(comparator, "command", nothing)
    comparator_command isa AbstractString ||
        error("confirmation comparator command must be a string")
    isempty(strip(comparator_command)) &&
        error("confirmation comparator command must be nonempty")
    occursin("documented", lowercase(comparator_command)) &&
        error("confirmation comparator command is not executable")
    row.command_sha256 == bytes2hex(SHA.sha256(comparator_command)) ||
        error("confirmation command hash differs from frozen manifest")
    row.comparator_configuration_sha256 == bytes2hex(SHA.sha256(
        comparator["configuration"])) ||
        error("confirmation comparator configuration differs from manifest")
    cell_key = (row.modality, row.seed)
    fixture_contract = _validate_stage2_fixture_index(
        manifest,
        freeze.fixture_path,
        dirname(freeze.development_ledger_path),
        target_cell = cell_key,
    )
    haskey(fixture_contract.confirmation_cells, cell_key) ||
        error("confirmation preflight cell is absent from the frozen fixture grid")
    cell = fixture_contract.confirmation_cells[cell_key]
    row.truth_reference_sha256 == fixture_contract.truth_sha256[row.seed] ||
        error("confirmation truth hash differs from frozen fixture index")
    row.input_bundle_sha256 == cell.input_bundle_sha256 ||
        error("confirmation full input bundle differs from frozen fixture index")
    row.error_profile_sha256 == cell.error_profile_sha256 ||
        error("confirmation error profile differs from frozen fixture index")
    consumed_contracts = get(
        comparator, "consumed_input_components", nothing)
    consumed_contracts isa AbstractDict ||
        error("confirmation comparator lacks consumed-input contracts")
    haskey(consumed_contracts, modality_key) ||
        error("confirmation comparator lacks a modality input contract")
    consumed_components = String.(consumed_contracts[modality_key])
    isempty(consumed_components) &&
        error("confirmation comparator consumes no input components")
    length(unique(consumed_components)) == length(consumed_components) ||
        error("confirmation comparator repeats an input component")
    issubset(
        Set(consumed_components), Set(keys(cell.input_component_sha256))) ||
        error("confirmation comparator consumes an unfrozen input component")
    expected_consumed_sha256 = stage2_input_component_digest(
        cell.input_component_sha256, consumed_components)
    row.consumed_input_sha256 == expected_consumed_sha256 ||
        error("confirmation consumed-input digest differs from frozen fixture index")
    return nothing
end

function stage2_attempt_ledger_fingerprint(
        path::AbstractString
)::NamedTuple{(:sha256, :rows, :final_attempt_id), Tuple{String, Int, String}}
    isfile(path) || error("attempt ledger does not exist: $(path)")
    lines = readlines(path)
    length(lines) >= 2 || error("attempt ledger has no attempt rows")
    header = split(first(lines), '\t')
    Tuple(header) == STAGE2_ATTEMPT_LEDGER_COLUMNS ||
        error("attempt ledger schema mismatch")
    final_fields = split(last(lines), '\t'; keepempty = true)
    length(final_fields) == length(STAGE2_ATTEMPT_LEDGER_COLUMNS) ||
        error("attempt ledger final row is malformed")
    return (
        sha256 = bytes2hex(SHA.sha256(read(path))),
        rows = length(lines) - 1,
        final_attempt_id = String(final_fields[1]),
    )
end

function _parse_ledger_bool(value::AbstractString, field::String)::Bool
    value == "true" && return true
    value == "false" && return false
    error("$(field) must be true or false")
end

function load_stage2_attempt_ledger(
        path::AbstractString;
        development_freeze::Union{Nothing, Stage2DevelopmentFreeze} = nothing
)::Vector{Stage2AttemptLedgerRow}
    isfile(path) || error("attempt ledger does not exist: $(path)")
    lines = readlines(path)
    isempty(lines) && error("attempt ledger is empty: $(path)")
    Tuple(split(first(lines), '\t')) == STAGE2_ATTEMPT_LEDGER_COLUMNS ||
        error("attempt ledger schema mismatch")
    rows = Stage2AttemptLedgerRow[]
    for line in Iterators.drop(lines, 1)
        fields = split(line, '\t'; keepempty = true)
        length(fields) == length(STAGE2_ATTEMPT_LEDGER_COLUMNS) ||
            error("attempt ledger contains a malformed row")
        row = Dict(zip(STAGE2_ATTEMPT_LEDGER_COLUMNS, fields))
        push!(rows, Stage2AttemptLedgerRow(;
            attempt_id = row["attempt_id"],
            recorded_at_utc = row["recorded_at_utc"],
            seed_partition = Symbol(row["seed_partition"]),
            seed = parse(Int, row["seed"]),
            modality = Symbol(row["modality"]),
            configuration_id = row["configuration_id"],
            git_sha = row["git_sha"],
            git_dirty = _parse_ledger_bool(row["git_dirty"], "git_dirty"),
            manifest_sha256 = row["manifest_sha256"],
            configuration_sha256 = row["configuration_sha256"],
            fixture_sha256 = row["fixture_sha256"],
            truth_reference_sha256 = row["truth_reference_sha256"],
            input_bundle_sha256 = row["input_bundle_sha256"],
            consumed_input_sha256 = row["consumed_input_sha256"],
            error_profile_sha256 = row["error_profile_sha256"],
            comparator_id = row["comparator_id"],
            comparator_version = row["comparator_version"],
            command_sha256 = row["command_sha256"],
            comparator_configuration_sha256 =
                row["comparator_configuration_sha256"],
            status = Symbol(row["status"]),
            result_path = row["result_path"],
            result_sha256 = row["result_sha256"],
            development_freeze_sha256 = row["development_freeze_sha256"],
            notes = row["notes"],
            development_freeze,
        ))
    end
    return rows
end

function _stage2_result_artifact_path(
        confirmation_ledger_path::AbstractString,
        result_path::AbstractString
)::String
    return _observation_relative_path(confirmation_ledger_path, result_path)
end

function _validate_release_observation_binding(
        observation::Stage2ReleaseObservation,
        row::Stage2AttemptLedgerRow,
        freeze::Stage2DevelopmentFreeze,
        freeze_sha256::String
)::Nothing
    observation.attempt_id == row.attempt_id ||
        error("release observation attempt ID differs from confirmation ledger")
    observation.modality == row.modality && observation.seed == row.seed ||
        error("release observation cell differs from confirmation ledger")
    observation.comparator_id == row.comparator_id ||
        error("release observation comparator differs from confirmation ledger")
    observation.comparator_version == row.comparator_version ||
        error("release observation comparator version differs from confirmation ledger")
    observation.comparator_command_sha256 == row.command_sha256 ||
        error("release observation command hash differs from confirmation ledger")
    observation.comparator_configuration_sha256 ==
    row.comparator_configuration_sha256 ||
        error("release observation comparator configuration differs from ledger")
    observation.input_bundle_sha256 == row.input_bundle_sha256 ||
        error("release observation input bundle differs from confirmation ledger")
    observation.comparator_consumed_input_sha256 == row.consumed_input_sha256 ||
        error("release observation consumed input differs from confirmation ledger")
    observation.truth_reference_sha256 == row.truth_reference_sha256 ||
        error("release observation truth reference differs from confirmation ledger")
    observation.development_freeze_sha256 == freeze_sha256 ||
        error("release observation differs from the canonical development freeze")
    observation.manifest_sha256 == freeze.manifest_sha256 ||
        error("release observation manifest hash differs from development freeze")
    observation.configuration_sha256 == freeze.configuration_sha256 ||
        error("release observation configuration hash differs from development freeze")
    observation.fixture_sha256 == freeze.fixture_sha256 ||
        error("release observation fixture hash differs from development freeze")
    for (value, field) in (
            (observation.native_result_sha256, "native_result_sha256"),
            (observation.comparator_result_sha256, "comparator_result_sha256"),
            (observation.native_consumed_input_sha256,
                "native_consumed_input_sha256"),
        )
        _validate_hex_digest(value, field, (64,))
    end
    return nothing
end

function _validate_subset_observation_contract(
        observation::Stage2ReleaseObservation,
        manifest::AbstractDict
)::Nothing
    modality_key = String(observation.modality)
    modality = manifest["modalities"][modality_key]
    observation.comparator_id in String.(modality["comparators"]) ||
        error("diagnostic subset comparator is not registered for its modality")
    comparator = manifest["comparators"][observation.comparator_id]
    comparator["version"] != "UNRESOLVED" ||
        error("diagnostic subset comparator version is unresolved")
    observation.comparator_version == comparator["version"] ||
        error("diagnostic subset comparator version differs from manifest")
    observation.comparator_command_sha256 == bytes2hex(SHA.sha256(
        comparator["command"])) ||
        error("diagnostic subset command hash differs from manifest")
    observation.comparator_configuration_sha256 == bytes2hex(SHA.sha256(
        comparator["configuration"])) ||
        error("diagnostic subset configuration hash differs from manifest")

    modality_components = Set(String.(modality["input_components"]))
    comparator_components = Set(String.(
        comparator["consumed_input_components"][modality_key]))
    comparator_components != modality_components ||
        error("diagnostic subset observation unexpectedly consumes the full bundle")
    Set(observation.comparator_consumed_input_components) ==
    comparator_components ||
        error("diagnostic subset components differ from manifest")
    Set(keys(observation.input_component_sha256)) == modality_components ||
        error("diagnostic subset does not bind the full modality bundle")
    stage2_input_component_digest(
        observation.input_component_sha256,
        collect(modality_components),
    ) == observation.input_bundle_sha256 ||
        error("diagnostic subset modality bundle differs from component digests")
    observation.native_consumed_input_sha256 == observation.input_bundle_sha256 ||
        error("diagnostic subset native run did not consume the full modality bundle")
    stage2_input_component_digest(
        observation.input_component_sha256,
        observation.comparator_consumed_input_components,
    ) == observation.comparator_consumed_input_sha256 ||
        error("diagnostic subset hash differs from component digests")
    observation.comparator_consumed_input_sha256 != observation.input_bundle_sha256 ||
        error("diagnostic subset is mislabeled as a full-input consumer")
    return nothing
end

"""
Evaluate the Stage-2 release gate from the complete sealed confirmation ledger.

The API intentionally accepts no caller-selected observation vector. Every
preregistered cell is loaded from the append-only ledger, every result artifact
is hash-verified before parsing, and any failed, repeated, missing, or extra
attempt blocks evaluation.
"""
function evaluate_stage2_release_gate(;
        manifest_path::AbstractString,
        development_ledger_path::AbstractString
)::Stage2ReleaseGateResult
    manifest_file = abspath(normpath(String(manifest_path)))
    development_path = abspath(normpath(String(development_ledger_path)))
    confirmation_path = stage2_confirmation_ledger_path(development_path)
    development_lock = development_path * ".lock"
    try
        mkdir(development_lock)
    catch
        error("development ledger is locked: $(development_path)")
    end
    try
        freeze_path = stage2_development_freeze_path(development_path)
        freeze = load_stage2_development_freeze(freeze_path)
        development_path == freeze.development_ledger_path ||
            error("development ledger path differs from frozen campaign identity")
        confirmation_path == freeze.confirmation_ledger_path ||
            error("confirmation ledger path differs from frozen campaign identity")
        manifest_file == freeze.manifest_path ||
            error("release manifest path differs from frozen campaign identity")
        for (artifact_path, expected_sha256, label) in (
                (freeze.configuration_path,
                    freeze.configuration_sha256, "configuration"),
                (freeze.fixture_path, freeze.fixture_sha256, "fixture"),
            )
            _validate_bundle_file(
                artifact_path, dirname(development_path), label)
            bytes2hex(SHA.sha256(read(artifact_path))) == expected_sha256 ||
                error("frozen $(label) artifact hash has changed")
        end
        frozen_configuration_id = _validate_stage2_configuration_artifact(
            freeze.configuration_path,
            dirname(development_path),
            freeze.configuration_sha256,
        )
        freeze_sha256 = bytes2hex(SHA.sha256(read(freeze_path)))
        fingerprint = stage2_attempt_ledger_fingerprint(development_path)
        fingerprint.sha256 == freeze.development_ledger_sha256 ||
            error("development ledger hash differs from freeze")
        fingerprint.rows == freeze.development_ledger_rows ||
            error("development ledger row count differs from freeze")
        fingerprint.final_attempt_id == freeze.final_attempt_id ||
            error("development ledger final attempt differs from freeze")
        manifest_bytes = read(manifest_file)
        bytes2hex(SHA.sha256(manifest_bytes)) == freeze.manifest_sha256 ||
            error("release manifest bytes differ from development freeze")
        manifest = TOML.parse(String(manifest_bytes))
        development_rows = load_stage2_attempt_ledger(development_path)
        last(development_rows).configuration_id == frozen_configuration_id ||
            error("development configuration ID differs from frozen artifact")
        _validate_development_gate_artifact(
            development_path, last(development_rows), manifest_file)

        confirmation_lock = confirmation_path * ".lock"
        try
            mkdir(confirmation_lock)
        catch
            error("confirmation ledger is locked: $(confirmation_path)")
        end
        try
            confirmation_bytes = read(confirmation_path)
            confirmation_sha256 = bytes2hex(SHA.sha256(confirmation_bytes))
            rows = load_stage2_attempt_ledger(
                confirmation_path; development_freeze = freeze)
            isempty(rows) && error("confirmation ledger has no attempts")
            all(row -> row.seed_partition == :confirmation, rows) ||
                error("confirmation ledger contains development attempts")
            all(row -> row.development_freeze_sha256 == freeze_sha256, rows) ||
                error("confirmation ledger contains a noncanonical freeze hash")
            all(row -> row.configuration_id == frozen_configuration_id, rows) ||
                error("confirmation ledger configuration IDs differ from artifact")

            events_by_attempt = Dict{String, Vector{Stage2AttemptLedgerRow}}()
            for row in rows
                push!(get!(events_by_attempt, row.attempt_id) do
                    Stage2AttemptLedgerRow[]
                end, row)
            end
            terminal_rows = Stage2AttemptLedgerRow[]
            for (attempt_id, events) in events_by_attempt
                length(events) == 2 ||
                    error("confirmation attempt $(attempt_id) is incomplete")
                first(events).status == :reserved &&
                    last(events).status in STAGE2_TERMINAL_ATTEMPT_STATUSES ||
                    error("confirmation attempt $(attempt_id) has invalid lifecycle")
                isequal(
                    _attempt_reservation_signature(first(events)),
                    _attempt_reservation_signature(last(events)),
                ) || error("confirmation attempt differs from its reservation")
                reservation_nonce_sha256 =
                    _validate_confirmation_reservation_receipt(
                        confirmation_path, first(events))
                _validate_confirmation_result_reservation_nonce(
                    confirmation_path,
                    last(events),
                    reservation_nonce_sha256,
                )
                first(events).recorded_at_utc <= last(events).recorded_at_utc ||
                    error("confirmation completion predates its reservation")
                push!(terminal_rows, last(events))
            end

            modality_contracts = manifest["modalities"]
            comparator_contracts = manifest["comparators"]
            required_modalities = [
                Symbol(modality) for (modality, contract) in modality_contracts
                if get(contract, "required", false)
            ]
            required_seeds = Int.(manifest["confirmation_seeds"])
            full_comparators = Dict{Symbol, Vector{String}}()
            subset_comparators = Dict{Symbol, Vector{String}}()
            for modality in required_modalities
                modality_key = String(modality)
                modality_components = Set(String.(
                    modality_contracts[modality_key]["input_components"]))
                full_comparators[modality] = String[]
                subset_comparators[modality] = String[]
                for comparator in
                    modality_contracts[modality_key]["comparators"]
                    comparator_id = String(comparator)
                    consumed_components = Set(String.(
                        comparator_contracts[comparator_id][
                            "consumed_input_components"][modality_key]))
                    target = consumed_components == modality_components ?
                             full_comparators : subset_comparators
                    push!(target[modality], comparator_id)
                end
                isempty(full_comparators[modality]) &&
                    error("required modality $(modality) has no full-input comparator")
            end
            expected_keys = Set(
                (modality, seed, String(comparator))
                for modality in required_modalities
                for seed in required_seeds
                for comparator in full_comparators[modality]
            )
            allowed_keys = union(expected_keys, Set(
                (modality, seed, String(comparator))
                for modality in required_modalities
                for seed in required_seeds
                for comparator in subset_comparators[modality]
            ))
            observed_keys = [
                (row.modality, row.seed, row.comparator_id)
                for row in terminal_rows
            ]
            length(unique(observed_keys)) == length(observed_keys) ||
                error("confirmation ledger repeats a modality/seed/comparator cell")
            issubset(Set(observed_keys), allowed_keys) ||
                error("confirmation ledger contains an unregistered cell")
            issubset(expected_keys, Set(observed_keys)) ||
                error("confirmation ledger does not exactly cover the frozen grid")
            gating_rows = filter(row ->
                (row.modality, row.seed, row.comparator_id) in expected_keys,
                terminal_rows,
            )
            all(row -> row.status == :passed, gating_rows) ||
                error("a required full-input confirmation attempt did not pass")
            length(unique(row.result_path for row in gating_rows)) ==
            length(gating_rows) ||
                error("confirmation attempts must use distinct result artifacts")
            for modality in required_modalities, seed in required_seeds
                cell_rows = filter(row ->
                    row.modality == modality && row.seed == seed, gating_rows)
                length(unique(row.configuration_id for row in cell_rows)) == 1 ||
                    error("confirmation cell uses different configuration IDs")
                length(unique(row.error_profile_sha256 for row in cell_rows)) == 1 ||
                    error("confirmation cell uses different error profiles")
                length(unique(row.input_bundle_sha256 for row in cell_rows)) == 1 ||
                    error("confirmation cell uses different input bundles")
            end

            diagnostic_rows = filter(row ->
                (row.modality, row.seed, row.comparator_id) ∉ expected_keys,
                terminal_rows,
            )
            diagnostic_status = Dict(
                (row.modality, row.seed, row.comparator_id) => row.status
                for row in diagnostic_rows
            )
            diagnostic_hashes = Dict{String, String}()
            for row in diagnostic_rows
                result_path = _stage2_result_artifact_path(
                    confirmation_path, row.result_path)
                _validate_bundle_file(
                    result_path,
                    dirname(confirmation_path),
                    "diagnostic result",
                )
                result_bytes = read(result_path)
                result_sha256 = bytes2hex(SHA.sha256(result_bytes))
                result_sha256 == row.result_sha256 ||
                    error("diagnostic result artifact hash differs from ledger")
                diagnostic_hashes[row.attempt_id] = result_sha256
                row.status == :passed || continue
                observation = _parse_stage2_release_observation(result_bytes)
                _validate_observation_raw_artifacts(observation, result_path)
                _validate_observation_metric_report(
                    observation, result_path, manifest["evaluation"])
                _validate_release_observation_binding(
                    observation, row, freeze, freeze_sha256)
                _validate_subset_observation_contract(observation, manifest)
            end

            observations = Stage2ReleaseObservation[]
            result_hashes = Dict{String, String}()
            for row in gating_rows
                result_path = _stage2_result_artifact_path(
                    confirmation_path, row.result_path)
                _validate_bundle_file(
                    result_path,
                    dirname(confirmation_path),
                    "confirmation result",
                )
                result_bytes = read(result_path)
                result_sha256 = bytes2hex(SHA.sha256(result_bytes))
                result_sha256 == row.result_sha256 ||
                    error("confirmation result artifact hash differs from ledger")
                observation = _parse_stage2_release_observation(result_bytes)
                _validate_observation_raw_artifacts(
                    observation, result_path)
                _validate_observation_metric_report(
                    observation, result_path, manifest["evaluation"])
                _validate_release_observation_binding(
                    observation, row, freeze, freeze_sha256)
                push!(observations, observation)
                result_hashes[row.attempt_id] = result_sha256
            end
            return _evaluate_stage2_release_observations(
                observations;
                manifest,
                development_freeze_sha256 = freeze_sha256,
                confirmation_ledger_sha256 = confirmation_sha256,
                attempt_result_sha256 = result_hashes,
                diagnostic_attempt_status = diagnostic_status,
                diagnostic_result_sha256 = diagnostic_hashes,
            )
        finally
            rm(confirmation_lock; force = true, recursive = true)
        end
    finally
        rm(development_lock; force = true, recursive = true)
    end
end

function _append_stage2_attempt_locked!(
        ledger_path::String,
        row::Stage2AttemptLedgerRow;
        allow_existing_attempt_id::Bool = false
)::String
    expected_header = join(STAGE2_ATTEMPT_LEDGER_COLUMNS, '\t')
    if isfile(ledger_path)
        open(ledger_path, "r") do stream
            eof(stream) && error("existing attempt ledger is empty")
            header = readline(stream)
            header == expected_header ||
                error("existing attempt ledger schema mismatch")
            for line in eachline(stream)
                fields = split(line, '\t'; keepempty = true)
                length(fields) == length(STAGE2_ATTEMPT_LEDGER_COLUMNS) ||
                    error("existing attempt ledger contains a malformed row")
                fields[1] == row.attempt_id && !allow_existing_attempt_id &&
                    error("duplicate attempt_id: $(row.attempt_id)")
            end
        end
    end
    new_ledger = !isfile(ledger_path)
    open(ledger_path, "a") do stream
        new_ledger && println(stream, expected_header)
        println(stream, join(_attempt_ledger_values(row), '\t'))
    end
    return ledger_path
end

function _attempt_reservation_signature(
        row::Stage2AttemptLedgerRow
)::NamedTuple
    return (
        attempt_id = row.attempt_id,
        seed_partition = row.seed_partition,
        seed = row.seed,
        modality = row.modality,
        configuration_id = row.configuration_id,
        git_sha = row.git_sha,
        git_dirty = row.git_dirty,
        manifest_sha256 = row.manifest_sha256,
        configuration_sha256 = row.configuration_sha256,
        fixture_sha256 = row.fixture_sha256,
        truth_reference_sha256 = row.truth_reference_sha256,
        input_bundle_sha256 = row.input_bundle_sha256,
        consumed_input_sha256 = row.consumed_input_sha256,
        error_profile_sha256 = row.error_profile_sha256,
        comparator_id = row.comparator_id,
        comparator_version = row.comparator_version,
        command_sha256 = row.command_sha256,
        comparator_configuration_sha256 =
            row.comparator_configuration_sha256,
        development_freeze_sha256 = row.development_freeze_sha256,
    )
end

function _confirmation_reservation_directory(
        confirmation_ledger_path::AbstractString,
        row::Stage2AttemptLedgerRow
)::String
    cell_key = "$(row.modality)\t$(row.seed)\t$(row.comparator_id)"
    cell_sha256 = bytes2hex(SHA.sha256(cell_key))
    return String(confirmation_ledger_path) * ".reservations/$(cell_sha256)"
end

function _reservation_sha256(row::Stage2AttemptLedgerRow)::String
    row.status == :reserved || error("reservation digest requires reserved status")
    return bytes2hex(SHA.sha256(join(_attempt_ledger_values(row), '\t')))
end

function _write_confirmation_reservation_receipt!(
        confirmation_ledger_path::AbstractString,
        row::Stage2AttemptLedgerRow
)::String
    reservation_nonce_sha256 = bytes2hex(SHA.sha256(
        Random.rand(Random.RandomDevice(), UInt8, 32)))
    receipt_directory = _confirmation_reservation_directory(
        confirmation_ledger_path, row)
    mkpath(dirname(receipt_directory))
    try
        mkdir(receipt_directory)
    catch
        error("confirmation cell already has an immutable reservation receipt")
    end
    receipt_path = joinpath(receipt_directory, "reservation.toml")
    open(receipt_path, "w") do stream
        TOML.print(stream, Dict(
            "attempt_id" => row.attempt_id,
            "modality" => String(row.modality),
            "seed" => row.seed,
            "comparator_id" => row.comparator_id,
            "development_freeze_sha256" => row.development_freeze_sha256,
            "reservation_nonce_sha256" => reservation_nonce_sha256,
            "reservation_sha256" => _reservation_sha256(row),
        ))
    end
    return receipt_path
end

function _validate_confirmation_reservation_receipt(
        confirmation_ledger_path::AbstractString,
        reservation::Stage2AttemptLedgerRow
)::String
    receipt_path = joinpath(
        _confirmation_reservation_directory(
            confirmation_ledger_path, reservation),
        "reservation.toml",
    )
    isfile(receipt_path) ||
        error("confirmation reservation receipt is missing")
    receipt = TOML.parsefile(receipt_path)
    receipt["attempt_id"] == reservation.attempt_id ||
        error("confirmation reservation receipt attempt differs")
    receipt["development_freeze_sha256"] ==
    reservation.development_freeze_sha256 ||
        error("confirmation reservation receipt freeze differs")
    receipt["reservation_sha256"] == _reservation_sha256(reservation) ||
        error("confirmation reservation receipt digest differs")
    reservation_nonce = get(receipt, "reservation_nonce_sha256", nothing)
    reservation_nonce isa AbstractString ||
        error("confirmation reservation receipt nonce must be a string")
    return _validate_hex_digest(
        reservation_nonce, "reservation_nonce_sha256", (64,))
end

function _validate_confirmation_result_reservation_nonce(
        confirmation_ledger_path::AbstractString,
        completion::Stage2AttemptLedgerRow,
        reservation_nonce_sha256::AbstractString
)::Nothing
    result_path = _stage2_result_artifact_path(
        confirmation_ledger_path, completion.result_path)
    _validate_bundle_file(
        result_path, dirname(confirmation_ledger_path), "confirmation result")
    result_bytes = read(result_path)
    bytes2hex(SHA.sha256(result_bytes)) == completion.result_sha256 ||
        error("confirmation result artifact hash differs from completion row")
    result = TOML.parse(String(result_bytes))
    result_nonce = get(result, "reservation_nonce_sha256", nothing)
    result_nonce isa AbstractString ||
        error("confirmation terminal result lacks a reservation nonce")
    validated_result_nonce = _validate_hex_digest(
        result_nonce, "reservation_nonce_sha256", (64,))
    validated_expected_nonce = _validate_hex_digest(
        reservation_nonce_sha256, "reservation_nonce_sha256", (64,))
    validated_result_nonce == validated_expected_nonce ||
        error("confirmation terminal result reservation nonce differs from receipt")
    return nothing
end

function _validate_confirmation_event_append(
        existing_rows::Vector{Stage2AttemptLedgerRow},
        row::Stage2AttemptLedgerRow,
        confirmation_ledger_path::AbstractString
)::Bool
    cell = (row.modality, row.seed, row.comparator_id)
    attempt_rows = filter(
        existing -> existing.attempt_id == row.attempt_id, existing_rows)
    cell_rows = filter(existing ->
        (existing.modality, existing.seed, existing.comparator_id) == cell,
        existing_rows,
    )
    if row.status == :reserved
        isempty(attempt_rows) ||
            error("confirmation attempt is already reserved: $(row.attempt_id)")
        isempty(cell_rows) ||
            error("confirmation cell has already been consumed: $(cell)")
        _write_confirmation_reservation_receipt!(
            confirmation_ledger_path, row)
        return false
    end
    row.status in STAGE2_TERMINAL_ATTEMPT_STATUSES ||
        error("confirmation completion must use a terminal status")
    length(attempt_rows) == 1 && first(attempt_rows).status == :reserved ||
        error("confirmation completion requires one prior reservation")
    isequal(
        _attempt_reservation_signature(first(attempt_rows)),
        _attempt_reservation_signature(row),
    ) || error("confirmation completion differs from its reservation")
    length(cell_rows) == 1 ||
        error("confirmation cell contains an invalid event sequence")
    reservation_nonce_sha256 = _validate_confirmation_reservation_receipt(
        confirmation_ledger_path, first(attempt_rows))
    _validate_confirmation_result_reservation_nonce(
        confirmation_ledger_path, row, reservation_nonce_sha256)
    return true
end

"""Append one development row, unless the canonical development seal exists."""
function append_stage2_attempt!(
        path::AbstractString,
        row::Stage2AttemptLedgerRow
)::String
    row.seed_partition == :development ||
        error("confirmation rows require append_stage2_confirmation_attempt!")
    ledger_path = abspath(normpath(String(path)))
    parent = dirname(ledger_path)
    isempty(parent) || mkpath(parent)
    lock_path = ledger_path * ".lock"
    try
        mkdir(lock_path)
    catch
        error("attempt ledger is locked: $(ledger_path)")
    end
    try
        isfile(stage2_development_freeze_path(ledger_path)) &&
            error("development ledger is sealed by its canonical freeze")
        return _append_stage2_attempt_locked!(ledger_path, row)
    finally
        rm(lock_path; force = true, recursive = true)
    end
end

"""
Append one confirmation row while holding the sealed development-ledger lock.

The fixed lock order is development then confirmation. This makes the
authorization check and append one critical section and prevents post-freeze
tuning or a TOCTOU change to the sealed development ledger.
"""
function append_stage2_confirmation_attempt!(
        row::Stage2AttemptLedgerRow;
        development_ledger_path::AbstractString
)::String
    row.seed_partition == :confirmation ||
        error("development rows require append_stage2_attempt!")
    development_path = abspath(normpath(String(development_ledger_path)))
    confirmation_path = stage2_confirmation_ledger_path(development_path)
    mkpath(dirname(confirmation_path))
    development_lock = development_path * ".lock"
    try
        mkdir(development_lock)
    catch
        error("development ledger is locked: $(development_path)")
    end
    try
        freeze = _validate_confirmation_freeze_locked(row, development_path)
        _validate_confirmation_reservation_preflight(row, freeze)
        confirmation_lock = confirmation_path * ".lock"
        try
            mkdir(confirmation_lock)
        catch
            error("attempt ledger is locked: $(confirmation_path)")
        end
        try
            existing_rows = isfile(confirmation_path) ?
                            load_stage2_attempt_ledger(
                confirmation_path; development_freeze = freeze) :
                            Stage2AttemptLedgerRow[]
            allow_existing_attempt_id = _validate_confirmation_event_append(
                existing_rows, row, confirmation_path)
            return _append_stage2_attempt_locked!(
                confirmation_path,
                row;
                allow_existing_attempt_id,
            )
        finally
            rm(confirmation_lock; force = true, recursive = true)
        end
    finally
        rm(development_lock; force = true, recursive = true)
    end
end

"""
Reserve a confirmation cell before invoking its runner, then append completion.

If the runner throws, an error completion consumes the reservation and keeps the
release gate blocked. The callback receives the canonical confirmation ledger
path and its unpredictable reservation nonce digest. Every terminal result must
bind that digest, proving it was created after the reservation receipt.
"""
function run_stage2_confirmation_attempt!(
        reservation::Stage2AttemptLedgerRow,
        runner::Function;
        development_ledger_path::AbstractString
)::Stage2AttemptLedgerRow
    reservation.status == :reserved ||
        error("confirmation runner requires a reserved row")
    ledger_path = append_stage2_confirmation_attempt!(
        reservation; development_ledger_path)
    reservation_nonce_sha256 = _validate_confirmation_reservation_receipt(
        ledger_path, reservation)
    completion = try
        runner(ledger_path, reservation_nonce_sha256)
    catch caught
        development_path = abspath(normpath(String(development_ledger_path)))
        freeze = load_stage2_development_freeze(
            stage2_development_freeze_path(development_path))
        result_path = joinpath("errors", "$(reservation.attempt_id).toml")
        absolute_result_path = joinpath(dirname(ledger_path), result_path)
        mkpath(dirname(absolute_result_path))
        open(absolute_result_path, "w") do stream
            TOML.print(stream, Dict(
                "schema" => "rhizomorph-stage2-attempt-error-v1",
                "attempt_id" => reservation.attempt_id,
                "error" => sprint(showerror, caught),
                "reservation_nonce_sha256" => reservation_nonce_sha256,
            ))
        end
        timestamp = Dates.format(
            Dates.now(Dates.UTC), Dates.dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"
        error_completion = Stage2AttemptLedgerRow(;
            attempt_id = reservation.attempt_id,
            recorded_at_utc = timestamp,
            seed_partition = reservation.seed_partition,
            seed = reservation.seed,
            modality = reservation.modality,
            configuration_id = reservation.configuration_id,
            git_sha = reservation.git_sha,
            git_dirty = reservation.git_dirty,
            manifest_sha256 = reservation.manifest_sha256,
            configuration_sha256 = reservation.configuration_sha256,
            fixture_sha256 = reservation.fixture_sha256,
            truth_reference_sha256 = reservation.truth_reference_sha256,
            input_bundle_sha256 = reservation.input_bundle_sha256,
            consumed_input_sha256 = reservation.consumed_input_sha256,
            error_profile_sha256 = reservation.error_profile_sha256,
            comparator_id = reservation.comparator_id,
            comparator_version = reservation.comparator_version,
            command_sha256 = reservation.command_sha256,
            comparator_configuration_sha256 =
                reservation.comparator_configuration_sha256,
            status = :error,
            result_path,
            result_sha256 = bytes2hex(SHA.sha256(read(absolute_result_path))),
            development_freeze_sha256 =
                reservation.development_freeze_sha256,
            notes = "runner threw; reservation permanently consumed",
            development_freeze = freeze,
        )
        append_stage2_confirmation_attempt!(
            error_completion; development_ledger_path)
        rethrow()
    end
    completion isa Stage2AttemptLedgerRow ||
        error("confirmation runner must return Stage2AttemptLedgerRow")
    append_stage2_confirmation_attempt!(completion; development_ledger_path)
    return completion
end
