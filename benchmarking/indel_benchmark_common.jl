# Shared helpers for the td-jt7r indel-aware pair-HMM benchmark family.
# =====================================================================
#
# `benchmarking/indel_frontier_fixed_toy_proof.jl` and
# `benchmarking/indel_frontier_runtime.jl` grew byte-for-byte duplicated
# provenance, artifact-publication, and fixture-simulation helpers that differed
# only by identifier prefix (`_indel_toy_*` vs `_indel_frontier_*`). A third
# consumer is queued, so the duplication is consolidated here.
#
# Convention matches the rest of `benchmarking/` (`benchmark_artifacts.jl`,
# `calibration_metrics.jl`, `quast_report_parsing.jl`): plain top-level
# definitions with no `module` and no include guard, consumed via
#
#     include(joinpath(@__DIR__, "indel_benchmark_common.jl"))
#
# with a `*_test.jl` companion.
#
# SCOPED DEFERRAL — not merged into `benchmark_artifacts.jl`.
# `benchmark_artifacts.jl` retains a SEPARATE lineage of provenance machinery
# (`collect_benchmark_provenance`, `write_benchmark_artifacts`, a
# tables/plots/logs/provenance subdirectory layout, JSON sidecars). The two
# lineages record overlapping facts but are not interchangeable: this one emits a
# single-row CSV manifest with a code-environment fingerprint that is re-asserted
# immediately before publication, and it publishes through an atomic staging
# rename. Unifying them would change artifact layout and manifest schema for
# `track_a_baseline_benchmark.jl`, `rhizomorph_prime_k_ablation_benchmark.jl`,
# and `rhizomorph_benchmark_harness.jl`, none of which are in scope for td-jt7r.
# The merge is deliberately deferred rather than overlooked; keeping the blast
# radius inside the indel benchmark family is the point.
#
# BYTE-IDENTITY WARNING — `indel_bench_simulate_reads`.
# `Mycelia.observe` draws its quality model from the GLOBAL RNG, so the simulated
# fixture's byte identity depends on the exact interleaving of local
# `MersenneTwister` draws and global draws. The draw order in
# `indel_bench_simulate_reads` is load-bearing and pinned by a golden hash in
# `indel_benchmark_common_test.jl`. Do not reorder, add, or remove RNG draws.

import BioAlignments
import BioSequences
import Dates
import FASTX
import Mycelia
import Random
import SHA

const INDEL_BENCH_MISSING_DEPENDENCY_SENTINEL = "MISSING"

# This file's own path, captured at include time. `indel_bench_code_environment`
# digests it alongside the calling script: a benchmark's behaviour now depends on
# BOTH files, so a manifest that named only the script would be a self-describing
# record with a hole in it.
const INDEL_BENCH_COMMON_SOURCE_PATH = @__FILE__

# Sidecar written into a benchmark's output directory. Publication is a
# whole-generation replace, so without this marker an aborted run and a
# successful one leave byte-indistinguishable output directories.
const INDEL_BENCH_RUN_STATUS_FILENAME = ".last_run_status"

# ---------------------------------------------------------------------------
# Hashing and dependency provenance
# ---------------------------------------------------------------------------

"""
    indel_bench_file_sha256(path) -> String

Hex-encoded SHA256 of the file at `path`.
"""
function indel_bench_file_sha256(path::String)::String
    return Base.bytes2hex(SHA.sha256(Base.read(path)))
end

"""
    indel_bench_optional_dependency_sha256(path) -> String

Hex-encoded SHA256 of `path`, or `INDEL_BENCH_MISSING_DEPENDENCY_SENTINEL` when
the file does not exist. Used for `Manifest.toml`, which is legitimately absent
in some environments.
"""
function indel_bench_optional_dependency_sha256(path::String)::String
    return isfile(path) ?
           indel_bench_file_sha256(path) :
           INDEL_BENCH_MISSING_DEPENDENCY_SENTINEL
end

"""
    indel_bench_dependency_provenance() -> NamedTuple

Active `Project.toml` path plus `Project.toml`/`Manifest.toml` digests. Errors
when there is no active project, because an unpinned dependency set makes the
recorded provenance unverifiable.
"""
function indel_bench_dependency_provenance()::NamedTuple
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
        project_toml_sha256 = indel_bench_file_sha256(active_project),
        manifest_toml_present = manifest_present,
        manifest_toml_sha256 = indel_bench_optional_dependency_sha256(
            manifest_path
        )
    )
end

"""
    indel_bench_code_environment(source_path) -> NamedTuple

Shared code/worktree/host environment core of a benchmark run manifest, for the
benchmark script at `source_path`.

This is deliberately NOT a whole run-provenance record: each benchmark keeps a
thin `_*_run_provenance` that merges this core with its own experiment constants
and derives its own `generation_id`. Only the environment is shared, because only
the environment is genuinely common.

`code_environment_fingerprint` digests the git HEAD, the tracked-worktree diff,
the benchmark source, THIS shared helper file, the dependency digests, and the
host/Julia identity. It is the value re-asserted by
[`indel_bench_assert_environment_unchanged`](@ref) immediately before artifacts
are published.

`benchmark_source_sha256` and `common_source_sha256` are recorded separately
because the executed logic is split across the two files; digesting only the
script would leave a self-describing manifest that omits half the code it
describes. The fingerprint already covered the shared file transitively via
`git_head_sha` + `git_tracked_diff_sha256`, so this makes an existing guarantee
explicit rather than adding one.
"""
function indel_bench_code_environment(source_path::String)::NamedTuple
    repository_root = normpath(joinpath(@__DIR__, ".."))
    git_head_sha = strip(
        Base.read(
        `git -C $repository_root rev-parse HEAD`,
        String
    ),
    )
    tracked_diff = Base.read(
        `git -C $repository_root diff --binary --no-ext-diff HEAD --`
    )
    tracked_diff_sha256 = Base.bytes2hex(SHA.sha256(tracked_diff))
    benchmark_source_sha256 = indel_bench_file_sha256(source_path)
    common_source_sha256 = indel_bench_file_sha256(
        INDEL_BENCH_COMMON_SOURCE_PATH
    )
    dependency = indel_bench_dependency_provenance()
    code_environment_components = (
        git_head_sha,
        tracked_diff_sha256,
        benchmark_source_sha256,
        common_source_sha256,
        dependency.project_toml_sha256,
        dependency.manifest_toml_sha256,
        string(VERSION),
        string(Threads.nthreads()),
        string(Sys.CPU_NAME),
        string(Sys.ARCH),
        string(Sys.KERNEL),
        string(Sys.CPU_THREADS)
    )
    code_environment_fingerprint = Base.bytes2hex(SHA.sha256(codeunits(join(
        code_environment_components, ":"
    ))))
    return (
        code_environment_fingerprint = code_environment_fingerprint,
        code_sha = git_head_sha,
        git_tracked_worktree_dirty = !isempty(tracked_diff),
        git_tracked_diff_sha256 = tracked_diff_sha256,
        benchmark_source_sha256 = benchmark_source_sha256,
        common_source_sha256 = common_source_sha256,
        active_project_path = dependency.active_project_path,
        project_toml_sha256 = dependency.project_toml_sha256,
        manifest_toml_present = dependency.manifest_toml_present,
        manifest_toml_sha256 = dependency.manifest_toml_sha256,
        julia_version = string(VERSION),
        julia_threads = Threads.nthreads(),
        cpu_name = string(Sys.CPU_NAME),
        architecture = string(Sys.ARCH),
        kernel = string(Sys.KERNEL),
        cpu_threads = Sys.CPU_THREADS
    )
end

"""
    indel_bench_assert_environment_unchanged(initial, source_path; context)

Error unless the dependency digests and the code-environment fingerprint still
match `initial`. `context` names the run in the error message (for example
`"fixed-toy run"`). Call this immediately before publishing so a mid-run edit to
the source, the worktree, or the environment cannot be silently attributed to the
manifest captured at start-of-run.
"""
function indel_bench_assert_environment_unchanged(
        initial::NamedTuple,
        source_path::String;
        context::String
)::Nothing
    current = indel_bench_code_environment(source_path)
    current.project_toml_sha256 == initial.project_toml_sha256 || error(
        "active Project.toml changed during the $(context); refusing to " *
        "publish artifacts"
    )
    current.manifest_toml_sha256 == initial.manifest_toml_sha256 || error(
        "active Manifest.toml changed during the $(context); refusing to " *
        "publish artifacts"
    )
    current.common_source_sha256 == initial.common_source_sha256 || error(
        "indel_benchmark_common.jl changed during the $(context); refusing " *
        "to publish artifacts"
    )
    current.code_environment_fingerprint ==
    initial.code_environment_fingerprint || error(
        "code/worktree/environment fingerprint changed during the " *
        "$(context); refusing to publish artifacts"
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Run status sidecar and start-of-run prechecks
# ---------------------------------------------------------------------------

"""
    indel_bench_run_status_path(output_dir) -> String

Path of the `.last_run_status` sidecar for `output_dir`.
"""
function indel_bench_run_status_path(output_dir::String)::String
    return joinpath(output_dir, INDEL_BENCH_RUN_STATUS_FILENAME)
end

"""
    indel_bench_read_run_status(output_dir) -> Union{Dict{String,String}, Nothing}

Parse the `.last_run_status` sidecar, or `nothing` when absent. The format is
`key=value` lines; unparseable lines are skipped, because a corrupt marker must
not be able to abort a run.
"""
function indel_bench_read_run_status(
        output_dir::String
)::Union{Dict{String, String}, Nothing}
    status_path = indel_bench_run_status_path(output_dir)
    isfile(status_path) || return nothing
    status = Dict{String, String}()
    for line in eachline(status_path)
        isempty(strip(line)) && continue
        separator = findfirst(isequal('='), line)
        separator === nothing && continue
        status[String(line[1:(separator - 1)])] =
            String(line[(separator + 1):end])
    end
    return status
end

"""
    indel_bench_write_run_status(output_dir, fields)

Write `fields` as `key=value` lines to the `.last_run_status` sidecar,
replacing any previous marker.
"""
function indel_bench_write_run_status(
        output_dir::String,
        fields::Vector{Pair{String, String}}
)::Nothing
    Base.Filesystem.mkpath(output_dir)
    open(indel_bench_run_status_path(output_dir), "w") do io
        for (key, value) in fields
            println(io, "$(key)=$(value)")
        end
    end
    return nothing
end

"""
    indel_bench_prior_generation_id(status) -> String

Generation identifier described by a `.last_run_status` reading, or `"none"`
when there is none.

A run in flight has already replaced the marker with `generation_id=pending`,
so the generation actually occupying the directory is carried in
`previous_generation_id`. Without that carry-through, `begin_run` would erase
the very identifier the publish-time report exists to surface.
"""
function indel_bench_prior_generation_id(
        status::Union{Dict{String, String}, Nothing}
)::String
    status === nothing && return "none"
    generation_id = get(status, "generation_id", "unknown")
    return generation_id == "pending" ?
           get(status, "previous_generation_id", "unknown") : generation_id
end

"""
    indel_bench_assert_publishable(output_dir, artifact_names)

Non-destructive precheck that publication into `output_dir` can succeed.

[`indel_bench_remove_prior_artifacts`](@ref) refuses to delete a DIRECTORY
squatting on an artifact name, which is correct but fires at publish time — at
the end of a run that may have taken hours. This runs the same read-only check
at start-of-run so the failure is immediate. It deletes and creates nothing, so
calling it cannot itself destroy a prior generation.
"""
function indel_bench_assert_publishable(
        output_dir::String,
        artifact_names::Tuple{Vararg{String}}
)::Nothing
    isdir(output_dir) || return nothing
    for artifact_name in artifact_names
        artifact_path = joinpath(output_dir, artifact_name)
        isdir(artifact_path) && error(
            "refusing to start: a directory squats on the artifact name " *
            "$(artifact_path); publication would fail at the end of the run"
        )
    end
    return nothing
end

"""
    indel_bench_begin_run(output_dir, artifact_names; context)

Start-of-run entry point: precheck publishability, report the generation
currently occupying `output_dir`, and mark the directory as having a run in
flight.

The marker exists because publication replaces a whole generation in place. An
aborted run leaves the PREVIOUS generation intact — which is the desired
durability property, but it also means the output directory of a failed run is
byte-indistinguishable from that of a successful one. A `status=running` marker
that is never advanced to `status=complete` is the difference.
"""
function indel_bench_begin_run(
        output_dir::String,
        artifact_names::Tuple{Vararg{String}};
        context::String
)::Nothing
    indel_bench_assert_publishable(output_dir, artifact_names)
    prior = indel_bench_read_run_status(output_dir)
    if prior === nothing
        println(
            "run status: no prior $(INDEL_BENCH_RUN_STATUS_FILENAME) in " *
            "$(output_dir)"
        )
    else
        println(
            "run status: prior status=$(get(prior, "status", "unknown")), " *
            "generation_id=$(indel_bench_prior_generation_id(prior)), " *
            "context=$(get(prior, "context", "unknown"))"
        )
        if get(prior, "status", "") == "running"
            @warn(
                "the previous run against this output directory never " *
                "published; the artifacts present are from an older " *
                "generation than that run",
                output_dir = output_dir
            )
        end
    end
    indel_bench_write_run_status(
        output_dir,
        [
            "status" => "running",
            "context" => context,
            "started_at" => string(Dates.now()),
            "generation_id" => "pending",
            "previous_generation_id" =>
                indel_bench_prior_generation_id(prior)
        ]
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Artifact publication
# ---------------------------------------------------------------------------

"""
    indel_bench_remove_prior_artifacts(output_dir, artifact_names)

Remove a prior generation's artifacts from `output_dir`, LAST name first.

`artifact_names` is ordered so the completion manifest is last; removing in
reverse invalidates the manifest before any data it vouches for disappears, so a
crash mid-removal can never leave a manifest that claims to describe files that
are gone.

This is called by [`indel_bench_publish_artifacts`](@ref) at publish time, NOT
before the run. Removing eagerly at start-of-run means an abort mid-run destroys
the last complete generation and leaves nothing in its place.
"""
function indel_bench_remove_prior_artifacts(
        output_dir::String,
        artifact_names::Tuple{Vararg{String}}
)::Nothing
    Base.Filesystem.mkpath(output_dir)
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

"""
    indel_bench_publish_artifacts(staging_dir, output_dir, artifact_names;
                                  context, generation_id)

Atomically promote a fully staged generation from `staging_dir` into
`output_dir`.

Every name in `artifact_names` must already exist in `staging_dir`; the
completeness check runs BEFORE anything in `output_dir` is touched. Only once the
new generation is known complete is the prior generation invalidated
(manifest-first) and the staged files renamed into place. A failure at any point
BEFORE this call therefore leaves the previous complete generation intact.

The window this does NOT cover is a crash DURING publication: the prior
generation is removed manifest-first and only then are the staged files renamed
in, so a crash between those two loops leaves a partial generation with no
manifest. That is deliberate — manifest-first removal makes the partial state
DETECTABLE rather than plausible — but it is not the same guarantee as the
pre-publish window, and the `.last_run_status` marker left at `status=running`
is what distinguishes it after the fact.

This is a BLIND overwrite of whatever occupies `output_dir`: nothing here reads
the prior manifest, so it cannot tell that the target is a newer or more
expensive generation than the one being published. `generation_id` is printed
before the replace so the overwrite is at least visible in the run log; callers
that must not clobber a different KIND of run (a full calibration, say) are
responsible for not pointing two kinds of run at one directory.

On success the `.last_run_status` marker is ALWAYS advanced to
`status=complete`, whether or not `context` and `generation_id` were supplied;
the descriptive fields fall back to `"unknown"` when they were not. Advancing
conditionally would mean a keyword-less publish left the `status=running` marker
written by [`indel_bench_begin_run`](@ref) in place, so a directory holding a
complete generation would report a run that never published — which is exactly
the aborted-run state the marker exists to distinguish. The keywords describe
the generation; they do not decide whether it finished.
"""
function indel_bench_publish_artifacts(
        staging_dir::String,
        output_dir::String,
        artifact_names::Tuple{Vararg{String}};
        context::Union{String, Nothing} = nothing,
        generation_id::Union{String, Nothing} = nothing
)::Nothing
    for artifact_name in artifact_names
        staging_path = joinpath(staging_dir, artifact_name)
        isfile(staging_path) || error(
            "staged artifact is missing: $(staging_path)"
        )
    end
    prior = indel_bench_read_run_status(output_dir)
    if prior !== nothing
        println(
            "publishing over prior generation " *
            "$(indel_bench_prior_generation_id(prior)) " *
            "(context=$(get(prior, "context", "unknown")))"
        )
    end
    indel_bench_remove_prior_artifacts(output_dir, artifact_names)
    for artifact_name in artifact_names
        Base.Filesystem.rename(
            joinpath(staging_dir, artifact_name),
            joinpath(output_dir, artifact_name)
        )
    end
    indel_bench_write_run_status(
        output_dir,
        [
            "status" => "complete",
            "context" => something(context, "unknown"),
            "completed_at" => string(Dates.now()),
            "generation_id" => something(generation_id, "unknown")
        ]
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Deterministic read simulation
# ---------------------------------------------------------------------------

"""
    indel_bench_simulate_reads(; genome_length, source_read_length, coverage,
                                 error_rate, seed, tech = :nanopore,
                                 identifier_prefix = nothing)

Simulate the shared deterministic long-read fixture and return
`(reads, reference)`.

`identifier_prefix` defaults to `"\$(tech)_read"` rather than to a hardcoded
`"nanopore_read"`. Read identifiers are part of the digested fixture bytes, so
two independent defaults meant a caller passing `tech = :illumina` silently got
`nanopore_read_*` identifiers baked into the golden digest. Deriving one from
the other makes that combination unrepresentable. The `:nanopore` default is
unchanged, so the pinned golden digest is unaffected.

BYTE IDENTITY IS LOAD-BEARING. `Mycelia.observe` samples its quality model from
the GLOBAL RNG, so the output depends on the exact interleaving of local
`MersenneTwister` draws with global ones. The following order is pinned and must
not change:

1. `Mycelia.random_fasta_record` BEFORE any seeding;
2. construct the local `MersenneTwister(seed)`, THEN `Random.seed!(seed)`;
3. `n_reads = ceil(Int, coverage * genome_length / source_read_length)`;
4. per read: local draw for the start position, then a local `Bool` draw for the
   reverse-complement flip, then `Mycelia.observe`.

Any reordering or extra draw silently changes the fixture bytes and invalidates
the downstream assembly oracles. `indel_benchmark_common_test.jl` pins the
resulting digest.
"""
function indel_bench_simulate_reads(;
        genome_length::Int,
        source_read_length::Int,
        coverage::Real,
        error_rate::Real,
        seed::Int,
        tech::Symbol = :nanopore,
        identifier_prefix::Union{String, Nothing} = nothing
)::Tuple{Vector{FASTX.FASTQ.Record}, BioSequences.LongDNA{4}}
    read_identifier_prefix = something(identifier_prefix, "$(tech)_read")
    reference_record = Mycelia.random_fasta_record(
        moltype = :DNA,
        seed = seed,
        L = genome_length
    )
    reference = FASTX.sequence(BioSequences.LongDNA{4}, reference_record)
    rng = Random.MersenneTwister(seed)
    Random.seed!(seed)

    n_reads = ceil(Int, coverage * genome_length / source_read_length)
    reads = FASTX.FASTQ.Record[]
    for read_index in 1:n_reads
        start_position = rand(
            rng,
            1:(genome_length - source_read_length + 1)
        )
        fragment = reference[
        start_position:(start_position + source_read_length - 1)
]
        if rand(rng, Bool)
            fragment = BioSequences.reverse_complement(fragment)
        end
        observed,
        qualities = Mycelia.observe(
            fragment; error_rate = error_rate, tech = tech
        )
        isempty(observed) && continue
        quality_string = String([Char(quality + 33) for quality in qualities])
        push!(
            reads,
            FASTX.FASTQ.Record(
                "$(read_identifier_prefix)_$(read_index)",
                string(observed),
                quality_string
            )
        )
    end
    return reads, reference
end

# ---------------------------------------------------------------------------
# Assembly summarisation
# ---------------------------------------------------------------------------

"""
    indel_bench_n50(contigs) -> Int

N50 of `contigs`: the length of the shortest contig in the set of longest contigs
that together cover at least half the assembled bases. Returns 0 for an empty
assembly.
"""
function indel_bench_n50(contigs::Vector{String})::Int
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

"""
    indel_bench_assembly_bytes(assembly) -> Vector{UInt8}

Canonical FASTA byte stream of an assembly, in emitted contig order. This is the
byte stream the fixed-toy oracle digests, so its formatting is load-bearing:
`>name\\n` then `sequence\\n`, no wrapping.
"""
function indel_bench_assembly_bytes(
        assembly::Mycelia.Rhizomorph.AssemblyResult
)::Vector{UInt8}
    buffer = IOBuffer()
    for (name, contig) in zip(assembly.contig_names, assembly.contigs)
        println(buffer, ">$(name)")
        println(buffer, contig)
    end
    return Base.take!(buffer)
end

"""
    indel_bench_best_contig_fit_alignment(contig, target) -> NamedTuple

Unit-cost FIT (semi-global) alignment of `contig` into `target`, returned as
`(matches, edit_distance, identity)`. Reference bases before and after the
contig's placement are free; everything inside the placement — mismatches,
insertions, and internal deletions — costs exactly 1, so `edit_distance` is the
Levenshtein distance between the contig and the reference window it covers, and

    identity = matches / (matches + edit_distance)

is per-base sequence accuracy over that window. Returns the all-zero record for
an empty contig.

`match = 0, mismatch = -1, gap_open = 0, gap_extend = -1` is not a tuned
scoring choice: it makes the alignment score the negated unit-cost edit
distance, so this is the edit-distance formulation, expressed in the score
model `BioAlignments.pairalign` requires.

This function exists because the global alignment used for contig SELECTION
cannot answer an accuracy question at all — see
`indel_bench_best_reference_alignment` below.
"""
function indel_bench_best_contig_fit_alignment(
        contig::AbstractString,
        target::AbstractString
)::NamedTuple
    isempty(contig) &&
        return (matches = 0, edit_distance = 0, identity = 0.0)
    model = BioAlignments.AffineGapScoreModel(
        match = 0, mismatch = -1, gap_open = 0, gap_extend = -1
    )
    result = BioAlignments.pairalign(
        BioAlignments.SemiGlobalAlignment(), contig, target, model
    )
    matches = Int(BioAlignments.count_matches(BioAlignments.alignment(result)))
    edit_distance = -Int(BioAlignments.score(result))
    denominator = matches + edit_distance
    identity = denominator == 0 ? 0.0 : matches / denominator
    return (
        matches = matches,
        edit_distance = edit_distance,
        identity = identity
    )
end

"""
    indel_bench_best_reference_alignment(contigs, reference) -> NamedTuple

Best contig-to-reference alignment over both orientations, returned as
`(best_contig_reference_coverage, best_contig_fit_identity,
best_contig_fit_matches, best_contig_fit_edit_distance, edit_distance, matches,
aligned_bases, contig_length, orientation)`. Selection maximises
`best_contig_reference_coverage`, then `matches`, then lower `edit_distance` —
unchanged behaviour. Returns the all-zero `:none` sentinel for an empty
assembly.

METRIC SEMANTICS — read this before quoting any field.

`Mycelia.assess_alignment` runs a GLOBAL (Levenshtein) alignment of one contig
against the WHOLE reference. A contig shorter than the reference is padded out
with reference-only deletion columns, so `aligned_bases` is the reference
length, and

    best_contig_reference_coverage = matches / aligned_bases

is NORMALISED BEST-CONTIG REFERENCE COVERAGE — a contiguity measure — NOT
sequence accuracy. Through commit 527c2d67 this field was named `identity`, and
a headline "nanopore 0.1 vs illumina 0.0655" was read as a 10%-vs-6.6%
accuracy claim it never supported. Those two numbers are 200/2000 and 131/2000:
best contigs of 200 bp and 131 bp on a 2 kb reference.

The degeneracy is exact, not incidental. Minimising global Levenshtein cost
maximises `matches + paired_columns`, and a short contig on a long reference
can always reach the maximum by scattering exactly-matching blocks across the
reference. So `matches` SATURATES at the contig length no matter how wrong the
contig is — a 200 bp contig carrying a substitution still scores 200 matches —
and therefore

    best_contig_reference_coverage == contig_length / reference_length

identically, in this regime. `matches`, `edit_distance`, and `aligned_bases`
are retained for continuity with the published artifacts, but they carry no
accuracy information here; they are functions of the contig and reference
lengths alone. Saturation also makes `orientation` arbitrary — a contig and
its reverse complement both reach the maximum — so that field records which
orientation won an effectively tied comparison, not which strand the contig
came from.

The accuracy statement comes instead from `indel_bench_best_contig_fit_alignment`
above, run once on the selected contig:

    best_contig_fit_identity = best_contig_fit_matches /
        (best_contig_fit_matches + best_contig_fit_edit_distance)

A fit alignment pins the contig to one reference window — internal indels and
mismatches cost — so it cannot be gamed by scattering. Report the two ratios
together: coverage answers "how much of the reference did the best contig
reach", fit identity answers "was what it reached correct".
"""
function indel_bench_best_reference_alignment(
        contigs::Vector{String},
        reference::BioSequences.LongDNA{4}
)::NamedTuple
    reference_forward = string(reference)
    reference_reverse = string(BioSequences.reverse_complement(reference))
    best = (
        best_contig_reference_coverage = 0.0,
        best_contig_fit_identity = 0.0,
        best_contig_fit_matches = 0,
        best_contig_fit_edit_distance = 0,
        edit_distance = typemax(Int),
        matches = 0,
        aligned_bases = 0,
        contig_length = 0,
        orientation = :none
    )
    best_contig = ""

    for contig in contigs
        for (orientation, target) in (
            (:forward, reference_forward),
            (:reverse, reference_reverse)
        )
            alignment = Mycelia.assess_alignment(contig, target)
            aligned_bases = alignment.total_matches + alignment.total_edits
            coverage = aligned_bases == 0 ?
                       0.0 : alignment.total_matches / aligned_bases
            candidate = (
                best_contig_reference_coverage = coverage,
                best_contig_fit_identity = 0.0,
                best_contig_fit_matches = 0,
                best_contig_fit_edit_distance = 0,
                edit_distance = alignment.total_edits,
                matches = alignment.total_matches,
                aligned_bases = aligned_bases,
                contig_length = length(contig),
                orientation = orientation
            )
            best_coverage = best.best_contig_reference_coverage
            if coverage > best_coverage ||
               (coverage == best_coverage &&
                candidate.matches > best.matches) ||
               (coverage == best_coverage &&
                candidate.matches == best.matches &&
                candidate.edit_distance < best.edit_distance)
                best = candidate
                best_contig = contig
            end
        end
    end

    best.orientation === :none && return best
    # Two extra alignments, for the winning contig only. The fit is
    # recomputed over BOTH orientations rather than reusing the selected
    # `best_target`, because saturation leaves `orientation` arbitrary: a
    # reverse-strand contig ties its own reverse complement at the saturated
    # coverage maximum, so the global alignment cannot resolve strand and the
    # recorded orientation may be the wrong one. Assemblies are
    # strand-agnostic, so the better fit is the accuracy answer.
    forward_fit = indel_bench_best_contig_fit_alignment(
        best_contig, reference_forward
    )
    reverse_fit = indel_bench_best_contig_fit_alignment(
        best_contig, reference_reverse
    )
    fit = forward_fit.edit_distance <= reverse_fit.edit_distance ?
          forward_fit : reverse_fit
    return merge(
        best,
        (
            best_contig_fit_identity = fit.identity,
            best_contig_fit_matches = fit.matches,
            best_contig_fit_edit_distance = fit.edit_distance
        )
    )
end

# ---------------------------------------------------------------------------
# Per-rung telemetry accessors
# ---------------------------------------------------------------------------

"""
    indel_bench_rung_value(rung, key, default) -> Any

Read `key` from a per-rung telemetry record, tolerating both `Dict` (symbol- or
string-keyed) and `NamedTuple` payloads. Returns `default` for anything else.
"""
function indel_bench_rung_value(
        rung::Any,
        key::Symbol,
        default::Any
)::Any
    if rung isa AbstractDict
        return get(rung, key, get(rung, string(key), default))
    elseif rung isa NamedTuple
        return get(rung, key, default)
    end
    return default
end

"""
    indel_bench_rung_has_key(rung, key) -> Bool

Whether a per-rung telemetry record carries `key` at all. Distinct from
[`indel_bench_rung_value`](@ref) returning a default, which cannot tell a missing
key from a key whose value happens to equal the default.
"""
function indel_bench_rung_has_key(rung::Any, key::Symbol)::Bool
    if rung isa AbstractDict
        return haskey(rung, key) || haskey(rung, string(key))
    elseif rung isa NamedTuple
        return haskey(rung, key)
    end
    return false
end

"""
    indel_bench_rung_counter(rung, key) -> Union{Int, Nothing}

Strict counter accessor: returns the value only when it is an exact nonnegative
`Int`, and `nothing` otherwise. Missing, negative, and non-`Int` (including
`Float64`) values all read as `nothing` so a telemetry schema regression cannot
be laundered into a plausible-looking total.
"""
function indel_bench_rung_counter(
        rung::Any,
        key::Symbol
)::Union{Int, Nothing}
    value = indel_bench_rung_value(rung, key, nothing)
    return value isa Int && value >= 0 ? value : nothing
end
