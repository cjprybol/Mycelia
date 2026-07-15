import CSV
import DataFrames
import Logging
import SHA
import Statistics

const ORACLE_MATRIX_SUMMARY_NAME = "oracle_matrix_summary.csv"
const ORACLE_MATRIX_MANIFEST_NAME = "oracle_matrix_manifest.csv"
const ORACLE_MATRIX_EXPECTED_SEEDS = collect(1:5)

mutable struct OracleMatrixCaptureLogger <: Logging.AbstractLogger
    delegate::Logging.ConsoleLogger
    captures::Vector{NamedTuple}
end

Logging.min_enabled_level(::OracleMatrixCaptureLogger)::Logging.LogLevel =
    Logging.Info

function Logging.shouldlog(
        ::OracleMatrixCaptureLogger,
        level::Logging.LogLevel,
        ::Module,
        ::Any,
        ::Any,
)::Bool
    return level >= Logging.Info
end

Logging.catch_exceptions(::OracleMatrixCaptureLogger)::Bool = false

function Logging.handle_message(
        logger::OracleMatrixCaptureLogger,
        level::Logging.LogLevel,
        message::Any,
        module_::Module,
        group::Any,
        id::Any,
        file::Any,
        line::Any;
        kwargs...,
)::Nothing
    if string(message) == "E2E nanopore-vs-illumina identity"
        values = (; kwargs...)
        push!(
            logger.captures,
            (
                seeds = Int.(collect(values.seeds)),
                nanopore_identity = Float64.(collect(values.id_np)),
                illumina_identity = Float64.(collect(values.id_il)),
                mean_nanopore_identity = Float64(values.mean_np),
                mean_illumina_identity = Float64(values.mean_il),
                wins_or_ties = Int(values.wins),
            ),
        )
    end
    Logging.handle_message(
        logger.delegate,
        level,
        message,
        module_,
        group,
        id,
        file,
        line;
        kwargs...,
    )
    return nothing
end

function oracle_matrix_file_sha256(path::String)::String
    return Base.bytes2hex(SHA.sha256(Base.read(path)))
end

function oracle_matrix_optional_file_sha256(path::String)::String
    return isfile(path) ? oracle_matrix_file_sha256(path) : "MISSING"
end

function oracle_matrix_provenance(test_path::String)::NamedTuple
    repository_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    git_head_sha = strip(Base.read(
        `git -C $repository_root rev-parse HEAD`,
        String,
    ))
    tracked_diff = Base.read(
        `git -C $repository_root diff --binary --no-ext-diff HEAD --`
    )
    tracked_diff_sha256 = Base.bytes2hex(SHA.sha256(tracked_diff))
    active_project = Base.active_project()
    active_project isa String || error(
        "an active Project.toml is required for oracle-matrix provenance"
    )
    realpath(dirname(active_project)) == realpath(repository_root) || error(
        "active project must be the audited Mycelia repository: " *
        "$(active_project)"
    )
    manifest_path = joinpath(dirname(active_project), "Manifest.toml")
    test_source_sha256 = oracle_matrix_file_sha256(test_path)
    capture_runner_sha256 = oracle_matrix_file_sha256(@__FILE__)
    project_toml_sha256 = oracle_matrix_file_sha256(active_project)
    manifest_toml_sha256 = oracle_matrix_optional_file_sha256(manifest_path)
    code_environment_components = (
        git_head_sha,
        tracked_diff_sha256,
        test_source_sha256,
        capture_runner_sha256,
        project_toml_sha256,
        manifest_toml_sha256,
        string(VERSION),
        string(Threads.nthreads()),
        string(Sys.CPU_NAME),
        string(Sys.ARCH),
        string(Sys.KERNEL),
        string(Sys.CPU_THREADS),
    )
    code_environment_fingerprint = Base.bytes2hex(SHA.sha256(codeunits(join(
        code_environment_components,
        ":",
    ))))
    generation_components = (
        code_environment_fingerprint,
        join(ORACLE_MATRIX_EXPECTED_SEEDS, ";"),
        "60",
        "40",
        "30",
        "0.05",
        "17",
    )
    generation_id = Base.bytes2hex(SHA.sha256(codeunits(join(
        generation_components,
        ":",
    ))))
    return (
        manifest_schema_version = 1,
        generation_id = generation_id,
        code_environment_fingerprint = code_environment_fingerprint,
        code_sha = git_head_sha,
        git_tracked_worktree_dirty = !isempty(tracked_diff),
        git_tracked_diff_sha256 = tracked_diff_sha256,
        test_source_path = relpath(test_path, repository_root),
        test_source_sha256 = test_source_sha256,
        capture_runner_sha256 = capture_runner_sha256,
        active_project_path = normpath(active_project),
        project_toml_sha256 = project_toml_sha256,
        manifest_toml_present = isfile(manifest_path),
        manifest_toml_sha256 = manifest_toml_sha256,
        julia_version = string(VERSION),
        julia_threads = Threads.nthreads(),
        cpu_name = string(Sys.CPU_NAME),
        architecture = string(Sys.ARCH),
        kernel = string(Sys.KERNEL),
        cpu_threads = Sys.CPU_THREADS,
        seeds = join(ORACLE_MATRIX_EXPECTED_SEEDS, ";"),
        genome_length = 60,
        source_read_length = 40,
        coverage = 30,
        error_rate = 0.05,
        k = 17,
        accuracy_used_for_classifier = false,
    )
end

function oracle_matrix_validate_capture(capture::NamedTuple)::NamedTuple
    capture.seeds == ORACLE_MATRIX_EXPECTED_SEEDS || error(
        "oracle matrix seeds changed: $(capture.seeds)"
    )
    length(capture.nanopore_identity) == length(capture.seeds) || error(
        "nanopore identity vector length does not match seeds"
    )
    length(capture.illumina_identity) == length(capture.seeds) || error(
        "Illumina identity vector length does not match seeds"
    )
    all(isfinite, capture.nanopore_identity) || error(
        "nanopore oracle matrix contains a non-finite identity"
    )
    all(isfinite, capture.illumina_identity) || error(
        "Illumina oracle matrix contains a non-finite identity"
    )
    mean_nanopore = Statistics.mean(capture.nanopore_identity)
    mean_illumina = Statistics.mean(capture.illumina_identity)
    wins_or_ties = count(
        capture.nanopore_identity .>= capture.illumina_identity
    )
    mean_nanopore == capture.mean_nanopore_identity || error(
        "logged nanopore mean does not match the captured per-seed values"
    )
    mean_illumina == capture.mean_illumina_identity || error(
        "logged Illumina mean does not match the captured per-seed values"
    )
    wins_or_ties == capture.wins_or_ties || error(
        "logged win count does not match the captured per-seed values"
    )
    mean_nanopore > mean_illumina || error(
        "nanopore mean identity did not beat Illumina"
    )
    wins_or_ties >= 3 || error(
        "nanopore did not win or tie a majority of seeds"
    )
    mean_nanopore > 0.80 || error(
        "nanopore mean identity did not exceed the recovery floor"
    )
    return (
        mean_nanopore_identity = mean_nanopore,
        mean_illumina_identity = mean_illumina,
        wins_or_ties = wins_or_ties,
        mean_nanopore_beats_illumina = true,
        majority_wins_or_ties = true,
        mean_nanopore_above_0_80 = true,
    )
end

function oracle_matrix_summary(
        capture::NamedTuple,
        validation::NamedTuple,
)::DataFrames.DataFrame
    rows = NamedTuple[]
    for index in eachindex(capture.seeds)
        push!(
            rows,
            (
                seed = capture.seeds[index],
                nanopore_identity = capture.nanopore_identity[index],
                illumina_identity = capture.illumina_identity[index],
                nanopore_wins_or_ties = capture.nanopore_identity[index] >=
                                         capture.illumina_identity[index],
                identity_delta = capture.nanopore_identity[index] -
                                 capture.illumina_identity[index],
                mean_nanopore_identity = validation.mean_nanopore_identity,
                mean_illumina_identity = validation.mean_illumina_identity,
                wins_or_ties = validation.wins_or_ties,
                mean_nanopore_beats_illumina =
                    validation.mean_nanopore_beats_illumina,
                majority_wins_or_ties = validation.majority_wins_or_ties,
                mean_nanopore_above_0_80 =
                    validation.mean_nanopore_above_0_80,
                accuracy_used_for_classifier = false,
            ),
        )
    end
    return DataFrames.DataFrame(rows)
end

function main()::Nothing
    summary_path = joinpath(@__DIR__, ORACLE_MATRIX_SUMMARY_NAME)
    manifest_path = joinpath(@__DIR__, ORACLE_MATRIX_MANIFEST_NAME)
    Base.rm(summary_path; force = true)
    Base.rm(manifest_path; force = true)

    test_path = normpath(joinpath(
        @__DIR__,
        "..",
        "..",
        "..",
        "test",
        "4_assembly",
        "indel_profile_assembly_test.jl",
    ))
    provenance = oracle_matrix_provenance(test_path)
    provenance.git_tracked_worktree_dirty && error(
        "tracked worktree must be clean before the oracle-matrix run"
    )

    logger = OracleMatrixCaptureLogger(
        Logging.ConsoleLogger(stderr, Logging.Info),
        NamedTuple[],
    )
    Logging.with_logger(logger) do
        Base.include(Main, test_path)
    end
    length(logger.captures) == 1 || error(
        "expected exactly one oracle-matrix log capture, found " *
        "$(length(logger.captures))"
    )
    capture = only(logger.captures)
    validation = oracle_matrix_validate_capture(capture)
    summary = oracle_matrix_summary(capture, validation)

    staging_dir = Base.Filesystem.mktempdir(
        @__DIR__;
        prefix = ".oracle-matrix-staging-",
    )
    try
        staged_summary = joinpath(staging_dir, ORACLE_MATRIX_SUMMARY_NAME)
        staged_manifest = joinpath(staging_dir, ORACLE_MATRIX_MANIFEST_NAME)
        CSV.write(staged_summary, summary)
        current = oracle_matrix_provenance(test_path)
        current.code_environment_fingerprint ==
        provenance.code_environment_fingerprint || error(
            "code/worktree/environment fingerprint changed during the " *
            "oracle-matrix run"
        )
        manifest = DataFrames.DataFrame([
            merge(
                provenance,
                validation,
                (
                    capture_count = length(logger.captures),
                    test_suite_passed = true,
                    summary_sha256 = oracle_matrix_file_sha256(staged_summary),
                ),
            ),
        ])
        CSV.write(staged_manifest, manifest)
        Base.mv(staged_summary, summary_path; force = true)
        Base.mv(staged_manifest, manifest_path; force = true)
    finally
        Base.rm(staging_dir; recursive = true, force = true)
    end

    println("td-jt7r.2 five-seed oracle matrix: PASS")
    println("  summary:  $(summary_path)")
    println("  manifest: $(manifest_path)")
    return nothing
end

main()
