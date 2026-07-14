# Fixed 2 kb / 1.2 kb-read / 8x / 5%-error proof for td-jt7r.2.
#
# Full proof (intentionally expensive; coordinate before running):
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/indel_frontier_fixed_toy_proof.jl
#
# Fixture-only smoke (does not assemble or decode):
#   LD_LIBRARY_PATH='' julia --project=. \
#     benchmarking/indel_frontier_fixed_toy_proof.jl --fixture-only
#
# Both correction arms receive byte-identical reads. The nanopore arm runs first,
# so its <120 s acceptance check conservatively includes pair-HMM compilation.
# Explicit `:illumina` is compared byte-for-byte with the default substitution-
# only oracle. No classifier threshold is selected from the accuracy result.

import BioSequences
import CSV
import DataFrames
import FASTX
import Mycelia
import Random
import SHA

const INDEL_TOY_GENOME_LENGTH = 2_000
const INDEL_TOY_SOURCE_READ_LENGTH = 1_200
const INDEL_TOY_MIN_OBSERVED_READ_LENGTH = 1_000
const INDEL_TOY_COVERAGE = 8
const INDEL_TOY_ERROR_RATE = 0.05
const INDEL_TOY_FIXTURE_SEED = 42
const INDEL_TOY_CORRECTOR_SEED = 1_042
const INDEL_TOY_MAX_K = 31
const INDEL_TOY_MAX_NANOPORE_WALL_SECONDS = 120.0
const INDEL_TOY_MISSING_DEPENDENCY_SENTINEL = "MISSING"
# Detached origin/master at the implementation base (548dc984) produced this
# deterministic explicit-Illumina assembly byte stream. Keeping the golden hash
# separate from the current default-profile comparison prevents both current arms
# from drifting together while still claiming pre-wiring byte identity.
const INDEL_TOY_PREWIRING_ILLUMINA_SHA256 =
    "d36e3b6a10685346aa7b0238b48b4ab7fcefbed88f82cad7d959b0a831cdd311"
const INDEL_TOY_DEFAULT_OUTPUT_DIR = joinpath(
    @__DIR__, "results", "td-jt7r-2-fixed-toy"
)
const INDEL_TOY_ARTIFACT_NAMES = (
    "fixed_toy_arm_summary.csv",
    "fixed_toy_rung_telemetry.csv",
    "fixed_toy_acceptance_checks.csv",
    "fixed_toy_run_manifest.csv",
)

function main(args::Vector{String} = ARGS)::Nothing
    options = _indel_toy_parse_args(args)
    if !options.fixture_only
        _indel_toy_remove_prior_artifacts(
            options.output_dir, INDEL_TOY_ARTIFACT_NAMES
        )
    end
    provenance = options.fixture_only ? nothing : _indel_toy_run_provenance()
    reads, reference = _indel_toy_make_fixture()
    observed_lengths = length.(FASTX.sequence.(String, reads))
    _indel_toy_print_fixture(reads, observed_lengths)
    if options.fixture_only
        minimum(observed_lengths) >= INDEL_TOY_MIN_OBSERVED_READ_LENGTH || error(
            "fixed fixture produced a read below " *
            "$(INDEL_TOY_MIN_OBSERVED_READ_LENGTH) bp"
        )
        println("Fixture-only smoke: PASS (no assembly or pair-HMM decode run)")
        return nothing
    end

    nanopore = _indel_toy_run_arm(reads, reference, :nanopore, "nanopore")
    illumina = _indel_toy_run_arm(reads, reference, :illumina, "illumina")
    oracle = _indel_toy_run_arm(reads, reference, nothing, "default_illumina_oracle")
    arms = [nanopore, illumina, oracle]
    for arm in arms
        _indel_toy_print_arm(arm)
    end

    summary = DataFrames.DataFrame(_indel_toy_summary_row.(arms))
    telemetry = _indel_toy_telemetry_table(arms)
    checks = _indel_toy_checks(reads, nanopore, illumina, oracle)
    println("\nFixed-toy acceptance checks")
    for row in DataFrames.eachrow(checks)
        status = row.passed ? "PASS" : "FAIL"
        println("  $(row.check): $(status) — $(row.detail)")
    end
    failed = String[
        row.check for row in DataFrames.eachrow(checks) if !row.passed
    ]
    if !isempty(failed)
        error("td-jt7r.2 fixed-toy proof failed: $(join(failed, ", "))")
    end

    staging_dir = Base.Filesystem.mktempdir(
        options.output_dir; prefix = ".fixed-toy-staging-"
    )
    try
        summary_staging_path = joinpath(
            staging_dir, INDEL_TOY_ARTIFACT_NAMES[1]
        )
        telemetry_staging_path = joinpath(
            staging_dir, INDEL_TOY_ARTIFACT_NAMES[2]
        )
        checks_staging_path = joinpath(
            staging_dir, INDEL_TOY_ARTIFACT_NAMES[3]
        )
        CSV.write(summary_staging_path, summary)
        CSV.write(telemetry_staging_path, telemetry)
        CSV.write(checks_staging_path, checks)
        manifest = DataFrames.DataFrame([
            merge(
                provenance,
                (
                    summary_sha256 = _indel_toy_file_sha256(
                        summary_staging_path
                    ),
                    telemetry_sha256 = _indel_toy_file_sha256(
                        telemetry_staging_path
                    ),
                    checks_sha256 = _indel_toy_file_sha256(
                        checks_staging_path
                    ),
                ),
            ),
        ])
        CSV.write(joinpath(staging_dir, INDEL_TOY_ARTIFACT_NAMES[4]), manifest)
        _indel_toy_assert_provenance_unchanged(something(provenance))
        _indel_toy_publish_artifacts(
            staging_dir, options.output_dir, INDEL_TOY_ARTIFACT_NAMES
        )
    finally
        Base.rm(staging_dir; recursive = true, force = true)
    end

    summary_path = joinpath(options.output_dir, INDEL_TOY_ARTIFACT_NAMES[1])
    telemetry_path = joinpath(options.output_dir, INDEL_TOY_ARTIFACT_NAMES[2])
    checks_path = joinpath(options.output_dir, INDEL_TOY_ARTIFACT_NAMES[3])
    manifest_path = joinpath(options.output_dir, INDEL_TOY_ARTIFACT_NAMES[4])

    println("\ntd-jt7r.2 fixed-toy/oracle proof: PASS")
    println("  summary:   $(summary_path)")
    println("  telemetry: $(telemetry_path)")
    println("  checks:    $(checks_path)")
    println("  manifest:  $(manifest_path)")
    return nothing
end

function _indel_toy_make_fixture()::Tuple{
        Vector{FASTX.FASTQ.Record}, BioSequences.LongDNA{4}}
    reference_record = Mycelia.random_fasta_record(
        moltype = :DNA,
        seed = INDEL_TOY_FIXTURE_SEED,
        L = INDEL_TOY_GENOME_LENGTH,
    )
    reference = FASTX.sequence(BioSequences.LongDNA{4}, reference_record)
    rng = Random.MersenneTwister(INDEL_TOY_FIXTURE_SEED)
    Random.seed!(INDEL_TOY_FIXTURE_SEED)

    n_reads = ceil(
        Int,
        INDEL_TOY_COVERAGE * INDEL_TOY_GENOME_LENGTH /
        INDEL_TOY_SOURCE_READ_LENGTH,
    )
    reads = FASTX.FASTQ.Record[]
    for read_index in 1:n_reads
        start_position = rand(
            rng,
            1:(INDEL_TOY_GENOME_LENGTH - INDEL_TOY_SOURCE_READ_LENGTH + 1),
        )
        fragment = reference[
            start_position:(start_position + INDEL_TOY_SOURCE_READ_LENGTH - 1)
        ]
        if rand(rng, Bool)
            fragment = BioSequences.reverse_complement(fragment)
        end
        observed, qualities = Mycelia.observe(
            fragment; error_rate = INDEL_TOY_ERROR_RATE, tech = :nanopore
        )
        isempty(observed) && continue
        quality_string = String([Char(quality + 33) for quality in qualities])
        push!(
            reads,
            FASTX.FASTQ.Record(
                "nanopore_read_$(read_index)", string(observed), quality_string
            ),
        )
    end
    return reads, reference
end

function _indel_toy_run_arm(
        reads::Vector{FASTX.FASTQ.Record},
        reference::BioSequences.LongDNA{4},
        sequencing_tech::Union{Symbol, Nothing},
        label::String,
)::NamedTuple
    Random.seed!(INDEL_TOY_CORRECTOR_SEED)
    local assembly::Mycelia.Rhizomorph.AssemblyResult
    wall_seconds = @elapsed begin
        if sequencing_tech === nothing
            assembly = Mycelia.Rhizomorph.assemble_genome(
                deepcopy(reads);
                k = INDEL_TOY_MAX_K,
                corrector = :iterative,
                strategy = :scalable,
            )
        else
            assembly = Mycelia.Rhizomorph.assemble_genome(
                deepcopy(reads);
                k = INDEL_TOY_MAX_K,
                corrector = :iterative,
                strategy = :scalable,
                sequencing_tech = sequencing_tech,
            )
        end
    end

    alignment = _indel_toy_best_reference_alignment(assembly.contigs, reference)
    stats = assembly.assembly_stats
    telemetry = get(stats, "indel_rung_telemetry", Dict{Symbol, Any}[])
    return (
        label = label,
        sequencing_tech = sequencing_tech === nothing ? :default : sequencing_tech,
        wall_seconds = wall_seconds,
        identity = alignment.identity,
        edit_distance = alignment.edit_distance,
        matches = alignment.matches,
        aligned_bases = alignment.aligned_bases,
        best_contig_length = alignment.contig_length,
        best_orientation = alignment.orientation,
        n_contigs = length(assembly.contigs),
        total_assembled_bases = sum(length, assembly.contigs; init = 0),
        largest_contig = isempty(assembly.contigs) ?
                         0 : maximum(length, assembly.contigs),
        n50 = _indel_toy_n50(assembly.contigs),
        k_progression = get(stats, "k_progression", Int[]),
        rung_vertex_counts = get(
            stats, "rung_vertex_counts", Dict{Int, Vector{Int}}()
        ),
        telemetry = telemetry,
        indel_requested = get(stats, "indel_requested", 0),
        indel_attempted = get(stats, "indel_attempted", 0),
        indel_completed = get(stats, "indel_completed", 0),
        indel_truncated = get(stats, "indel_truncated", 0),
        indel_engaged = get(stats, "indel_engaged", 0),
        trace_contract_errors = get(stats, "trace_contract_errors", 0),
        window_anchor_rejections = get(
            stats, "window_anchor_rejections", 0),
        window_divergences = get(stats, "window_divergences", 0),
        assembly_bytes = _indel_toy_assembly_bytes(assembly),
    )
end

function _indel_toy_best_reference_alignment(
        contigs::Vector{String},
        reference::BioSequences.LongDNA{4},
)::NamedTuple
    reference_forward = string(reference)
    reference_reverse = string(BioSequences.reverse_complement(reference))
    best = (
        identity = 0.0,
        edit_distance = typemax(Int),
        matches = 0,
        aligned_bases = 0,
        contig_length = 0,
        orientation = :none,
    )

    for contig in contigs
        for (orientation, target) in (
                (:forward, reference_forward),
                (:reverse, reference_reverse),
        )
            alignment = Mycelia.assess_alignment(contig, target)
            aligned_bases = alignment.total_matches + alignment.total_edits
            identity = aligned_bases == 0 ?
                       0.0 : alignment.total_matches / aligned_bases
            candidate = (
                identity = identity,
                edit_distance = alignment.total_edits,
                matches = alignment.total_matches,
                aligned_bases = aligned_bases,
                contig_length = length(contig),
                orientation = orientation,
            )
            if candidate.identity > best.identity ||
               (candidate.identity == best.identity &&
                candidate.matches > best.matches) ||
               (candidate.identity == best.identity &&
                candidate.matches == best.matches &&
                candidate.edit_distance < best.edit_distance)
                best = candidate
            end
        end
    end
    return best
end

function _indel_toy_n50(contigs::Vector{String})::Int
    isempty(contigs) && return 0
    lengths = sort(length.(contigs); rev = true)
    target = cld(sum(lengths), 2)
    cumulative = 0
    for contig_length in lengths
        cumulative += contig_length
        cumulative >= target && return contig_length
    end
    return 0
end

function _indel_toy_assembly_bytes(
        assembly::Mycelia.Rhizomorph.AssemblyResult
)::Vector{UInt8}
    buffer = IOBuffer()
    for (name, contig) in zip(assembly.contig_names, assembly.contigs)
        println(buffer, ">$(name)")
        println(buffer, contig)
    end
    return Base.take!(buffer)
end

function _indel_toy_summary_row(arm::NamedTuple)::NamedTuple
    return (
        arm = arm.label,
        sequencing_tech = string(arm.sequencing_tech),
        wall_seconds = arm.wall_seconds,
        identity = arm.identity,
        edit_distance = arm.edit_distance,
        matches = arm.matches,
        aligned_bases = arm.aligned_bases,
        best_contig_length = arm.best_contig_length,
        best_orientation = string(arm.best_orientation),
        n_contigs = arm.n_contigs,
        total_assembled_bases = arm.total_assembled_bases,
        largest_contig = arm.largest_contig,
        n50 = arm.n50,
        k_progression = join(arm.k_progression, ";"),
        indel_requested = arm.indel_requested,
        indel_attempted = arm.indel_attempted,
        indel_completed = arm.indel_completed,
        indel_truncated = arm.indel_truncated,
        indel_engaged = arm.indel_engaged,
        trace_contract_errors = arm.trace_contract_errors,
        window_anchor_rejections = arm.window_anchor_rejections,
        window_divergences = arm.window_divergences,
        assembly_sha256 = Base.bytes2hex(SHA.sha256(arm.assembly_bytes)),
    )
end

function _indel_toy_telemetry_table(
        arms::AbstractVector
)::DataFrames.DataFrame
    rows = NamedTuple[]
    for arm in arms
        for rung in arm.telemetry
            push!(
                rows,
                (
                    arm = arm.label,
                    ladder_index = _indel_toy_rung_value(rung, :ladder_index, missing),
                    k = _indel_toy_rung_value(rung, :k, missing),
                    iteration = _indel_toy_rung_value(rung, :iteration, missing),
                    profile_requested = _indel_toy_rung_value(
                        rung, :profile_requested, false
                    ),
                    requested = _indel_toy_rung_value(rung, :requested, 0),
                    attempted = _indel_toy_rung_value(rung, :attempted, 0),
                    completed = _indel_toy_rung_value(rung, :completed, 0),
                    truncated = _indel_toy_rung_value(rung, :truncated, 0),
                    engaged = _indel_toy_rung_value(rung, :engaged, 0),
                    admitted = _indel_toy_rung_value(rung, :admitted, false),
                    graph_source = string(
                        _indel_toy_rung_value(rung, :graph_source, :missing)
                    ),
                    decision_reason = string(
                        _indel_toy_rung_value(rung, :decision_reason, :missing)
                    ),
                    frontier_work_limit = _indel_toy_rung_value(
                        rung, :frontier_work_limit, missing
                    ),
                ),
            )
        end
    end
    return DataFrames.DataFrame(rows)
end

function _indel_toy_rung_value(
        rung::Any,
        key::Symbol,
        default::Any,
)::Any
    if rung isa AbstractDict
        return get(rung, key, get(rung, string(key), default))
    elseif rung isa NamedTuple
        return get(rung, key, default)
    end
    return default
end

function _indel_toy_rung_has_key(rung::Any, key::Symbol)::Bool
    if rung isa AbstractDict
        return haskey(rung, key) || haskey(rung, string(key))
    elseif rung isa NamedTuple
        return haskey(rung, key)
    end
    return false
end

function _indel_toy_telemetry_total(
        arm::NamedTuple,
        key::Symbol,
)::Union{Int, Nothing}
    values = Union{Int, Nothing}[
        _indel_toy_rung_counter(rung, key) for rung in arm.telemetry
    ]
    any(isnothing, values) && return nothing
    return sum(something(value) for value in values; init = 0)
end

function _indel_toy_rung_counter(
        rung::Any,
        key::Symbol,
)::Union{Int, Nothing}
    value = _indel_toy_rung_value(rung, key, nothing)
    return value isa Int && value >= 0 ? value : nothing
end

function _indel_toy_validate_rung_telemetry(
        arms::Tuple{Vararg{NamedTuple}}
)::NamedTuple
    required_keys = (:requested, :attempted, :completed, :truncated, :engaged)
    for arm in arms
        for (row_index, rung) in enumerate(arm.telemetry)
            counters = Dict(
                key => _indel_toy_rung_counter(rung, key)
                for key in required_keys
            )
            invalid_keys = Symbol[
                key for key in required_keys if isnothing(counters[key])
            ]
            if !isempty(invalid_keys)
                return (
                    passed = false,
                    detail = "arm=$(arm.label), row=$(row_index): exact " *
                             "nonnegative Int required for " *
                             "$(join(invalid_keys, ","))",
                )
            end
            requested = something(counters[:requested])
            attempted = something(counters[:attempted])
            completed = something(counters[:completed])
            truncated = something(counters[:truncated])
            engaged = something(counters[:engaged])
            attempted <= requested || return (
                passed = false,
                detail = "arm=$(arm.label), row=$(row_index): " *
                         "attempted=$(attempted) > requested=$(requested)",
            )
            completed + truncated <= attempted || return (
                passed = false,
                detail = "arm=$(arm.label), row=$(row_index): completed+" *
                         "truncated=$(completed + truncated) > " *
                         "attempted=$(attempted)",
            )
            engaged <= completed || return (
                passed = false,
                detail = "arm=$(arm.label), row=$(row_index): " *
                         "engaged=$(engaged) > completed=$(completed)",
            )
        end
    end
    return (
        passed = true,
        detail = "all per-rung counters are exact nonnegative Int values " *
                 "with requested/attempted/completed/truncated/engaged " *
                 "inequalities satisfied",
    )
end

function _indel_toy_totals_consistent(arm::NamedTuple)::Bool
    return all(
        _indel_toy_telemetry_total(arm, key) ==
        getproperty(arm, Symbol("indel_$(key)"))
        for key in (:requested, :attempted, :completed, :truncated, :engaged)
    )
end

function _indel_toy_totals_zero(arm::NamedTuple)::Bool
    return all(
        getproperty(arm, Symbol("indel_$(key)")) == 0 &&
        _indel_toy_telemetry_total(arm, key) == 0
        for key in (:requested, :attempted, :completed, :truncated, :engaged)
    )
end

function _indel_toy_totals_detail(arm::NamedTuple)::String
    labels = "requested/attempted/completed/truncated/engaged"
    stats = "$(arm.indel_requested)/$(arm.indel_attempted)/" *
            "$(arm.indel_completed)/$(arm.indel_truncated)/" *
            "$(arm.indel_engaged)"
    telemetry = "$(_indel_toy_telemetry_total(arm, :requested))/" *
                "$(_indel_toy_telemetry_total(arm, :attempted))/" *
                "$(_indel_toy_telemetry_total(arm, :completed))/" *
                "$(_indel_toy_telemetry_total(arm, :truncated))/" *
                "$(_indel_toy_telemetry_total(arm, :engaged))"
    return "$(labels): stats=$(stats), telemetry=$(telemetry)"
end

function _indel_toy_checks(
        reads::Vector{FASTX.FASTQ.Record},
        nanopore::NamedTuple,
        illumina::NamedTuple,
        oracle::NamedTuple,
)::DataFrames.DataFrame
    observed_lengths = length.(FASTX.sequence.(String, reads))
    telemetry_validation = _indel_toy_validate_rung_telemetry(
        (nanopore, illumina, oracle)
    )
    initial_k = isempty(nanopore.k_progression) ?
                nothing : first(nanopore.k_progression)
    has_noninitial_completion = telemetry_validation.passed &&
                                initial_k !== nothing && any(
        Int(_indel_toy_rung_value(rung, :k, initial_k)) > initial_k &&
        Int(_indel_toy_rung_value(rung, :attempted, 0)) > 0 &&
        Int(_indel_toy_rung_value(rung, :completed, 0)) > 0
        for rung in nanopore.telemetry
    )
    required_telemetry_keys = (
        :requested, :attempted, :completed, :truncated, :engaged
    )
    telemetry_schema_complete = all(
        !isempty(arm.telemetry) && all(
            all(_indel_toy_rung_has_key(rung, key) for key in required_telemetry_keys)
            for rung in arm.telemetry
        )
        for arm in (nanopore, illumina, oracle)
    )
    oracle_byte_identical = illumina.assembly_bytes == oracle.assembly_bytes
    illumina_sha256 = Base.bytes2hex(SHA.sha256(illumina.assembly_bytes))

    rows = [
        (
            check = "reads_are_1000bp_plus",
            passed = minimum(observed_lengths) >= INDEL_TOY_MIN_OBSERVED_READ_LENGTH,
            detail = "observed range=$(minimum(observed_lengths))-" *
                     "$(maximum(observed_lengths)) bp",
        ),
        (
            check = "nanopore_requested_indel_decode",
            passed = nanopore.indel_requested > 0,
            detail = "requested=$(nanopore.indel_requested)",
        ),
        (
            check = "noninitial_rung_attempted_and_completed",
            passed = has_noninitial_completion,
            detail = "initial_k=$(initial_k), attempted=$(nanopore.indel_attempted), " *
                     "completed=$(nanopore.indel_completed)",
        ),
        (
            check = "nanopore_beats_identical_read_illumina",
            passed = nanopore.identity > illumina.identity,
            detail = "nanopore=$(nanopore.identity), illumina=$(illumina.identity)",
        ),
        (
            check = "nanopore_under_120_seconds",
            passed = nanopore.wall_seconds < INDEL_TOY_MAX_NANOPORE_WALL_SECONDS,
            detail = "wall=$(nanopore.wall_seconds) s",
        ),
        (
            check = "nanopore_decode_not_truncated",
            passed = nanopore.indel_truncated == 0,
            detail = "truncated=$(nanopore.indel_truncated)",
        ),
        (
            check = "nanopore_substitution_window_contract_clean",
            passed = nanopore.window_divergences == 0,
            detail = "window_divergences=$(nanopore.window_divergences)",
        ),
        (
            check = "illumina_all_indel_counters_zero",
            passed = _indel_toy_totals_zero(illumina),
            detail = _indel_toy_totals_detail(illumina),
        ),
        (
            check = "default_oracle_all_indel_counters_zero",
            passed = _indel_toy_totals_zero(oracle),
            detail = _indel_toy_totals_detail(oracle),
        ),
        (
            check = "illumina_byte_identical_to_default_oracle",
            passed = oracle_byte_identical,
            detail = "illumina_bytes=$(length(illumina.assembly_bytes)), " *
                     "oracle_bytes=$(length(oracle.assembly_bytes))",
        ),
        (
            check = "illumina_byte_identical_to_prewiring_oracle",
            passed = illumina_sha256 == INDEL_TOY_PREWIRING_ILLUMINA_SHA256,
            detail = "sha256=$(illumina_sha256), " *
                     "origin_master=$(INDEL_TOY_PREWIRING_ILLUMINA_SHA256)",
        ),
        (
            check = "per_rung_telemetry_schema_complete",
            passed = telemetry_schema_complete,
            detail = "required keys=requested/attempted/completed/truncated/engaged",
        ),
        (
            check = "per_rung_telemetry_values_valid",
            passed = telemetry_validation.passed,
            detail = telemetry_validation.detail,
        ),
        (
            check = "nanopore_per_rung_totals_consistent",
            passed = _indel_toy_totals_consistent(nanopore),
            detail = _indel_toy_totals_detail(nanopore),
        ),
        (
            check = "illumina_per_rung_totals_consistent",
            passed = _indel_toy_totals_consistent(illumina),
            detail = _indel_toy_totals_detail(illumina),
        ),
        (
            check = "default_oracle_per_rung_totals_consistent",
            passed = _indel_toy_totals_consistent(oracle),
            detail = _indel_toy_totals_detail(oracle),
        ),
        (
            check = "nanopore_attempts_all_classified",
            passed = nanopore.indel_attempted ==
                     nanopore.indel_completed + nanopore.indel_truncated,
            detail = "attempted=$(nanopore.indel_attempted), " *
                     "completed=$(nanopore.indel_completed), " *
                     "truncated=$(nanopore.indel_truncated)",
        ),
        (
            check = "nanopore_trace_contract_clean",
            passed = nanopore.trace_contract_errors == 0,
            detail = "trace_contract_errors=$(nanopore.trace_contract_errors)",
        ),
        (
            check = "both_assemblies_nonempty",
            passed = nanopore.n_contigs > 0 && illumina.n_contigs > 0,
            detail = "nanopore=$(nanopore.n_contigs), illumina=$(illumina.n_contigs)",
        ),
    ]
    return DataFrames.DataFrame(rows)
end

function _indel_toy_print_fixture(
        reads::Vector{FASTX.FASTQ.Record},
        observed_lengths::Vector{Int},
)::Nothing
    println("td-jt7r.2 fixed toy fixture")
    println("  fixture_seed:          $(INDEL_TOY_FIXTURE_SEED)")
    println("  corrector_seed:        $(INDEL_TOY_CORRECTOR_SEED)")
    println(
        "  reference/source read: $(INDEL_TOY_GENOME_LENGTH)/" *
        "$(INDEL_TOY_SOURCE_READ_LENGTH) bp"
    )
    println(
        "  reads/coverage/error:  $(length(reads))/$(INDEL_TOY_COVERAGE)x/" *
        "$(INDEL_TOY_ERROR_RATE)"
    )
    println(
        "  observed read lengths: min=$(minimum(observed_lengths)), " *
        "max=$(maximum(observed_lengths))"
    )
    println("  max_k:                 $(INDEL_TOY_MAX_K)")
    println("  nanopore budget:       $(INDEL_TOY_MAX_NANOPORE_WALL_SECONDS) s")
    return nothing
end

function _indel_toy_print_arm(arm::NamedTuple)::Nothing
    println("\n$(arm.label) correction arm")
    println("  wall_seconds:          $(round(arm.wall_seconds; digits = 3))")
    println("  identity:              $(round(arm.identity; digits = 6))")
    println("  edit_distance:         $(arm.edit_distance)")
    println("  matches/aligned:       $(arm.matches)/$(arm.aligned_bases)")
    println("  contigs/total/largest: $(arm.n_contigs)/" *
            "$(arm.total_assembled_bases)/$(arm.largest_contig)")
    println("  N50:                   $(arm.n50)")
    println("  k_progression:         $(arm.k_progression)")
    println("  graph vertices/pass:   $(arm.rung_vertex_counts)")
    println("  indel requested:       $(arm.indel_requested)")
    println("  indel attempted:       $(arm.indel_attempted)")
    println("  indel completed:       $(arm.indel_completed)")
    println("  indel truncated:       $(arm.indel_truncated)")
    println("  indel engaged:         $(arm.indel_engaged)")
    println("  trace_contract_errors: $(arm.trace_contract_errors)")
    println("  window anchor rejects: $(arm.window_anchor_rejections)")
    println("  window_divergences:    $(arm.window_divergences)")
    println("  per-rung telemetry:")
    if isempty(arm.telemetry)
        println("    none")
    else
        for rung in arm.telemetry
            println(
                "    k=$(_indel_toy_rung_value(rung, :k, missing)) " *
                "iter=$(_indel_toy_rung_value(rung, :iteration, missing)) " *
                "requested=$(_indel_toy_rung_value(rung, :requested, 0)) " *
                "attempted=$(_indel_toy_rung_value(rung, :attempted, 0)) " *
                "completed=$(_indel_toy_rung_value(rung, :completed, 0)) " *
                "truncated=$(_indel_toy_rung_value(rung, :truncated, 0)) " *
                "engaged=$(_indel_toy_rung_value(rung, :engaged, 0)) " *
                "admitted=$(_indel_toy_rung_value(rung, :admitted, false)) " *
                "reason=$(_indel_toy_rung_value(rung, :decision_reason, :missing))"
            )
        end
    end
    return nothing
end

function _indel_toy_parse_args(args::Vector{String})::NamedTuple
    fixture_only = false
    output_dir = INDEL_TOY_DEFAULT_OUTPUT_DIR
    seen = Set{String}()
    index = 1
    while index <= length(args)
        flag = args[index]
        flag in seen && throw(ArgumentError("duplicate argument: $(flag)"))
        if flag == "--fixture-only"
            fixture_only = true
            push!(seen, flag)
            index += 1
        elseif flag == "--output-dir"
            push!(seen, flag)
            index == length(args) && throw(
                ArgumentError("--output-dir requires a value")
            )
            value = args[index + 1]
            (isempty(value) || startswith(value, "--")) && throw(
                ArgumentError("--output-dir requires a nonempty value")
            )
            output_dir = value
            index += 2
        else
            throw(ArgumentError("unknown argument: $(flag)"))
        end
    end
    return (fixture_only = fixture_only, output_dir = output_dir)
end

function _indel_toy_remove_prior_artifacts(
        output_dir::String,
        artifact_names::Tuple{Vararg{String}},
)::Nothing
    Base.Filesystem.mkpath(output_dir)
    # Invalidate the completion manifest and PASS-bearing checks before data.
    for artifact_index in length(artifact_names):-1:1
        artifact_name = artifact_names[artifact_index]
        artifact_path = joinpath(output_dir, artifact_name)
        isdir(artifact_path) && error(
            "refusing to replace artifact directory: $(artifact_path)"
        )
        Base.rm(artifact_path; force = true)
    end
    return nothing
end

function _indel_toy_publish_artifacts(
        staging_dir::String,
        output_dir::String,
        artifact_names::Tuple{Vararg{String}},
)::Nothing
    for artifact_name in artifact_names
        staging_path = joinpath(staging_dir, artifact_name)
        isfile(staging_path) || error(
            "staged artifact is missing: $(staging_path)"
        )
    end
    for artifact_name in artifact_names
        Base.Filesystem.rename(
            joinpath(staging_dir, artifact_name),
            joinpath(output_dir, artifact_name),
        )
    end
    return nothing
end

function _indel_toy_run_provenance()::NamedTuple
    repository_root = normpath(joinpath(@__DIR__, ".."))
    git_head_sha = strip(
        Base.read(
            `git -C $repository_root rev-parse HEAD`,
            String,
        ),
    )
    tracked_diff = Base.read(
        `git -C $repository_root diff --binary --no-ext-diff HEAD --`
    )
    tracked_diff_sha256 = Base.bytes2hex(SHA.sha256(tracked_diff))
    benchmark_source_sha256 = _indel_toy_file_sha256(@__FILE__)
    dependency = _indel_toy_dependency_provenance()
    code_environment_components = (
        git_head_sha,
        tracked_diff_sha256,
        benchmark_source_sha256,
        dependency.project_toml_sha256,
        dependency.manifest_toml_sha256,
        string(VERSION),
        string(Threads.nthreads()),
        string(Sys.CPU_NAME),
        string(Sys.ARCH),
        string(Sys.KERNEL),
        string(Sys.CPU_THREADS),
    )
    code_environment_fingerprint = Base.bytes2hex(SHA.sha256(codeunits(join(
        code_environment_components, ":"
    ))))
    run_fingerprint = join(
        (
            code_environment_fingerprint,
            string(INDEL_TOY_FIXTURE_SEED),
            string(INDEL_TOY_CORRECTOR_SEED),
        ),
        ":",
    )
    generation_id = Base.bytes2hex(SHA.sha256(codeunits(run_fingerprint)))
    return (
        manifest_schema_version = 1,
        generation_id = generation_id,
        code_environment_fingerprint = code_environment_fingerprint,
        code_sha = git_head_sha,
        git_tracked_worktree_dirty = !isempty(tracked_diff),
        git_tracked_diff_sha256 = tracked_diff_sha256,
        benchmark_source_sha256 = benchmark_source_sha256,
        active_project_path = dependency.active_project_path,
        project_toml_sha256 = dependency.project_toml_sha256,
        manifest_toml_present = dependency.manifest_toml_present,
        manifest_toml_sha256 = dependency.manifest_toml_sha256,
        julia_version = string(VERSION),
        julia_threads = Threads.nthreads(),
        cpu_name = string(Sys.CPU_NAME),
        architecture = string(Sys.ARCH),
        kernel = string(Sys.KERNEL),
        cpu_threads = Sys.CPU_THREADS,
        genome_length = INDEL_TOY_GENOME_LENGTH,
        source_read_length = INDEL_TOY_SOURCE_READ_LENGTH,
        minimum_observed_read_length = INDEL_TOY_MIN_OBSERVED_READ_LENGTH,
        coverage = INDEL_TOY_COVERAGE,
        error_rate = INDEL_TOY_ERROR_RATE,
        fixture_seed = INDEL_TOY_FIXTURE_SEED,
        corrector_seed = INDEL_TOY_CORRECTOR_SEED,
        max_k = INDEL_TOY_MAX_K,
        nanopore_wall_ceiling_seconds = INDEL_TOY_MAX_NANOPORE_WALL_SECONDS,
        prewiring_illumina_sha256 = INDEL_TOY_PREWIRING_ILLUMINA_SHA256,
    )
end

function _indel_toy_assert_provenance_unchanged(
        initial::NamedTuple
)::Nothing
    current = _indel_toy_run_provenance()
    current.project_toml_sha256 == initial.project_toml_sha256 || error(
        "active Project.toml changed during the fixed-toy run; refusing to " *
        "publish artifacts"
    )
    current.manifest_toml_sha256 == initial.manifest_toml_sha256 || error(
        "active Manifest.toml changed during the fixed-toy run; refusing to " *
        "publish artifacts"
    )
    current.code_environment_fingerprint ==
    initial.code_environment_fingerprint || error(
        "code/worktree/environment fingerprint changed during the fixed-toy " *
        "run; refusing to publish artifacts"
    )
    return nothing
end

function _indel_toy_dependency_provenance()::NamedTuple
    active_project = Base.active_project()
    active_project isa String || error(
        "an active Project.toml is required for benchmark provenance"
    )
    isfile(active_project) || error(
        "active Project.toml does not exist: $(active_project)"
    )
    manifest_path = joinpath(dirname(active_project), "Manifest.toml")
    manifest_present = isfile(manifest_path)
    return (
        active_project_path = normpath(active_project),
        project_toml_sha256 = _indel_toy_file_sha256(active_project),
        manifest_toml_present = manifest_present,
        manifest_toml_sha256 = _indel_toy_optional_dependency_sha256(
            manifest_path
        ),
    )
end

function _indel_toy_optional_dependency_sha256(path::String)::String
    return isfile(path) ?
           _indel_toy_file_sha256(path) :
           INDEL_TOY_MISSING_DEPENDENCY_SENTINEL
end

function _indel_toy_file_sha256(path::String)::String
    return Base.bytes2hex(SHA.sha256(Base.read(path)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
