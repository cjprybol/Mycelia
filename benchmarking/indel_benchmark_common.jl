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

import BioSequences
import Distributions
import FASTX
import Mycelia
import Random
import SHA

const INDEL_BENCH_MISSING_DEPENDENCY_SENTINEL = "MISSING"

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
the benchmark source, the dependency digests, and the host/Julia identity. It is
the value re-asserted by [`indel_bench_assert_environment_unchanged`](@ref)
immediately before artifacts are published.
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
    dependency = indel_bench_dependency_provenance()
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
    current.code_environment_fingerprint ==
    initial.code_environment_fingerprint || error(
        "code/worktree/environment fingerprint changed during the " *
        "$(context); refusing to publish artifacts"
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
    indel_bench_publish_artifacts(staging_dir, output_dir, artifact_names)

Atomically promote a fully staged generation from `staging_dir` into
`output_dir`.

Every name in `artifact_names` must already exist in `staging_dir`; the
completeness check runs BEFORE anything in `output_dir` is touched. Only once the
new generation is known complete is the prior generation invalidated
(manifest-first) and the staged files renamed into place. A failure at any point
before this call therefore leaves the previous complete generation intact.
"""
function indel_bench_publish_artifacts(
        staging_dir::String,
        output_dir::String,
        artifact_names::Tuple{Vararg{String}}
)::Nothing
    for artifact_name in artifact_names
        staging_path = joinpath(staging_dir, artifact_name)
        isfile(staging_path) || error(
            "staged artifact is missing: $(staging_path)"
        )
    end
    indel_bench_remove_prior_artifacts(output_dir, artifact_names)
    for artifact_name in artifact_names
        Base.Filesystem.rename(
            joinpath(staging_dir, artifact_name),
            joinpath(output_dir, artifact_name)
        )
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Deterministic read simulation
# ---------------------------------------------------------------------------

"""
    indel_bench_simulate_reads(; genome_length, source_read_length, coverage,
                                 error_rate, seed, tech = :nanopore,
                                 identifier_prefix = "nanopore_read")

Simulate the shared deterministic long-read fixture and return
`(reads, reference)`.

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
        identifier_prefix::String = "nanopore_read"
)::Tuple{Vector{FASTX.FASTQ.Record}, BioSequences.LongDNA{4}}
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
                "$(identifier_prefix)_$(read_index)",
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
    indel_bench_best_reference_alignment(contigs, reference) -> NamedTuple

Best contig-to-reference alignment over both orientations, returned as
`(identity, edit_distance, matches, aligned_bases, contig_length, orientation)`.
Ties break on `matches`, then on lower `edit_distance`. Returns the zero-identity
`:none` sentinel for an empty assembly.
"""
function indel_bench_best_reference_alignment(
        contigs::Vector{String},
        reference::BioSequences.LongDNA{4}
)::NamedTuple
    reference_forward = string(reference)
    reference_reverse = string(BioSequences.reverse_complement(reference))
    best = (
        identity = 0.0,
        edit_distance = typemax(Int),
        matches = 0,
        aligned_bases = 0,
        contig_length = 0,
        orientation = :none
    )

    for contig in contigs
        for (orientation, target) in (
            (:forward, reference_forward),
            (:reverse, reference_reverse)
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
                orientation = orientation
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

# ---------------------------------------------------------------------------
# Paired Wilcoxon signed-rank test
# ---------------------------------------------------------------------------

"""
    indel_bench_signed_rank(deltas) -> NamedTuple

Two-sided paired Wilcoxon signed-rank test on the paired differences `deltas`,
returning `(n_nonzero, statistic, p_two_sided, method)`.

- `n_nonzero` — count of nonzero differences. Zeros are dropped (Wilcoxon's
  original handling), so `n_nonzero` is the effective sample size.
- `statistic` — `V`, the sum of the ranks of the POSITIVE differences, ranking
  `abs.(deltas)` ascending with midranks for ties. `Float64` because midranks are
  half-integers.
- `p_two_sided` — two-sided p-value; `NaN` when `n_nonzero == 0`.
- `method` — `:exact`, `:normal_approximation`, or `:undefined`.

For `n_nonzero <= 20` the exact conditional distribution of `V` is enumerated by
subset-sum dynamic programming over the observed (possibly tied) ranks, which is
the exact randomization distribution given `abs.(deltas)`; the two-sided p-value
is the symmetric tail `P(|V - T/2| >= |v_obs - T/2|)` where `T` is the total rank
sum. Above 20 a normal approximation is used, with the standard tie correction
`-sum(t^3 - t)/48` on the variance and a continuity correction of 0.5 on the
numerator.

Implemented directly rather than via `HypothesisTests` so the benchmark family
does not acquire a dependency for one test.

Non-finite `deltas` are rejected — a silent `NaN` would otherwise be dropped as
"not positive" and quietly shrink the sample.
"""
function indel_bench_signed_rank(deltas::AbstractVector{<:Real})::NamedTuple
    all(isfinite, deltas) || throw(
        ArgumentError("deltas must all be finite")
    )
    nonzero = Float64[Float64(delta) for delta in deltas if delta != 0]
    n_nonzero = length(nonzero)
    if n_nonzero == 0
        return (
            n_nonzero = 0,
            statistic = 0.0,
            p_two_sided = NaN,
            method = :undefined
        )
    end

    magnitudes = abs.(nonzero)
    order = sortperm(magnitudes)
    ranks = Vector{Float64}(undef, n_nonzero)
    tie_sizes = Int[]
    index = 1
    while index <= n_nonzero
        stop = index
        while stop < n_nonzero &&
            magnitudes[order[stop + 1]] == magnitudes[order[index]]
            stop += 1
        end
        midrank = (index + stop) / 2
        for tied_index in index:stop
            ranks[order[tied_index]] = midrank
        end
        push!(tie_sizes, stop - index + 1)
        index = stop + 1
    end

    statistic = sum(
        ranks[position] for position in 1:n_nonzero if nonzero[position] > 0;
        init = 0.0
    )

    if n_nonzero <= 20
        # Midranks are half-integers; double them so the subset-sum DP is exact
        # in integer arithmetic.
        scaled_ranks = Int[round(Int, 2 * rank) for rank in ranks]
        total = sum(scaled_ranks)
        counts = zeros(Float64, total + 1)
        counts[1] = 1.0
        for scaled_rank in scaled_ranks
            for target in total:-1:scaled_rank
                counts[target + 1] += counts[target - scaled_rank + 1]
            end
        end
        observed = round(Int, 2 * statistic)
        center = total / 2
        deviation = abs(observed - center)
        tail = sum(
            counts[value + 1]
            for value in 0:total if abs(value - center) >= deviation - 1e-9;
            init = 0.0
        )
        p_two_sided = min(1.0, tail / 2.0^n_nonzero)
        return (
            n_nonzero = n_nonzero,
            statistic = statistic,
            p_two_sided = p_two_sided,
            method = :exact
        )
    end

    mean_statistic = n_nonzero * (n_nonzero + 1) / 4
    tie_adjustment = sum(
        Float64(size)^3 - Float64(size) for size in tie_sizes; init = 0.0
    )
    variance = n_nonzero * (n_nonzero + 1) * (2 * n_nonzero + 1) / 24 -
               tie_adjustment / 48
    variance > 0 || return (
        n_nonzero = n_nonzero,
        statistic = statistic,
        p_two_sided = NaN,
        method = :normal_approximation
    )
    deviation = abs(statistic - mean_statistic)
    corrected = max(deviation - 0.5, 0.0)
    z = corrected / sqrt(variance)
    p_two_sided = min(
        1.0, 2 * Distributions.ccdf(Distributions.Normal(), z)
    )
    return (
        n_nonzero = n_nonzero,
        statistic = statistic,
        p_two_sided = p_two_sided,
        method = :normal_approximation
    )
end
