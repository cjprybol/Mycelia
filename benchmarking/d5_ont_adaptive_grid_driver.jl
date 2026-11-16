# D5 ONT adaptive-grid DRIVER (td-yddpr.5, pre-registration lock step 2).
#
# Threads the axis registry (d5_ont_adaptive_grid_axes.jl) through
# simulate -> survival -> assemble -> QUAST. The CLI lives in
# d5_ont_adaptive_grid_run.jl; this file is include-safe (no ARGS, no side
# effects), so tests can load every function here.
#
# Design (bead td-yddpr.5, 2026-11-16 notes, all binding):
#   * Survival everywhere. Every read set records S(k) and S_hat(k) over ALL 23
#     D5_CANDIDATE_LADDER rungs. Survival is computed on the RAW reads once per
#     read set (technology, identity, coverage, organism, seed) and shared by both
#     correctors and every k-policy; it is the k* input. The corrector therefore
#     does not move k*. Assembled :iterative cells additionally record survival on
#     their own corrected reads, labelled diagnostic.
#   * Phase A: the fixed-rung ONT cells, the validation stratum, the Illumina
#     control's fixed-rung cells, and every survival vector.
#   * Phase B (after a human sets tau / tau_hat): the :adaptive cells, assembled
#     at k* = the largest rung with S(k) >= tau. NoFeasibleK, never a fallback.
#   * estimate_residual_error is recorded as a DIAGNOSTIC and never reaches k.
#   * k > 64 is rejected by assemble_genome; such cells record
#     STATUS_INFEASIBLE_K with null metrics (user decision 2026-11-16).
#
# Definitions are from rhizomorph-paper DEC-2026-08-09:
#   S(k)     = fraction of the G-k+1 reference k-mer positions observed
#              error-free in at least one read (canonical: reads come from both
#              strands).
#   S_hat(k) = criterion 5's reference-free surrogate; `s_hat` is an exact port
#              of decisions/surrogate_spectrum.py (rhizomorph-paper f3a5ed5).
#   k*       = max{k in ladder : S(k) >= tau}, subject to the collision floor
#              k >= ceil(log_|A|(G / p)), p = 0.01.

module D5OntAdaptiveGridDriver

import BioSequences
import FASTX
import JSON
import Kmers
import Mycelia

include(joinpath(@__DIR__, "d5_ont_adaptive_grid_axes.jl"))
import .D5OntAdaptiveGridAxes

include(joinpath(@__DIR__, "quast_censoring.jl"))
include(joinpath(@__DIR__, "quast_min_contig.jl"))
include(joinpath(@__DIR__, "peak_rss_sampler.jl"))

# identities_from_sam lives in ont_read_identity.jl, whose top level defines
# Lambda-specific consts; isolate it so none of them leak into this namespace.
module ReadIdentity
include(joinpath(@__DIR__, "ont_read_identity.jl"))
end

# === Pins (user decision 2026-11-16) ===
#
# Identity is the only factor varying along the identity axis: one error model
# and one qscore model for every level. The strings are PINNED here and recorded
# per read set from the argv actually built; they are never inferred from the
# target percentage.
const D5_BADREAD_ERROR_MODEL = "nanopore2023"
const D5_BADREAD_QSCORE_MODEL = "nanopore2023"
const D5_BADREAD_IDENTITY_PARAMS = Dict(
    90 => "90,98,5",
    95 => "95,99,2.5",  # the existing default; measured mean BLAST identity 94.4%
    99 => "99,100,1"
)
const D5_BADREAD_LENGTH_DIST = Mycelia.BADREAD_NANOPORE_DEFAULTS.length_dist

# ART parameters passed verbatim to simulate_illumina_reads; recorded as given.
const D5_ILLUMINA_PARAMS = (; read_length = 150, mflen = 500, sdev = 10,
    seqSys = "HS25", paired = true)

# Accession and expected size per registered organism (same as Track A).
const D5_ORGANISM_REFERENCES = Dict(
    :Lambda => (accession = "NC_001416", expected_bp = 48_502),
    :T4 => (accession = "NC_000866", expected_bp = 168_903)
)

# === Illumina control (user decision 2026-11-16) ===
#
# Derived from the registry; no literals. No identity axis. Crosses the corrector
# axis; tau / tau_hat are calibrated on its :none stratum by a human.
const D5_ILLUMINA_COVERAGES = D5OntAdaptiveGridAxes.D5_COVERAGES
const D5_ILLUMINA_CORRECTORS = D5OntAdaptiveGridAxes.D5_CORRECTORS
const D5_ILLUMINA_RUNGS = D5OntAdaptiveGridAxes.D5_COMPARATOR_RUNGS
const D5_ILLUMINA_K_POLICIES = D5OntAdaptiveGridAxes.D5_K_POLICIES
const D5_ILLUMINA_ORGANISMS = D5OntAdaptiveGridAxes.D5_ORGANISMS
const D5_ILLUMINA_SEEDS = D5OntAdaptiveGridAxes.D5_SEEDS
const D5_ILLUMINA_CALIBRATION_CORRECTOR = first(D5OntAdaptiveGridAxes.D5_CORRECTORS)  # :none

function n_illumina_phase_a_cells()
    length(D5_ILLUMINA_COVERAGES) * length(D5_ILLUMINA_CORRECTORS) *
    length(D5_ILLUMINA_RUNGS) * length(D5_ILLUMINA_ORGANISMS) * length(D5_ILLUMINA_SEEDS)
end
function n_illumina_phase_b_cells()
    length(D5_ILLUMINA_COVERAGES) * length(D5_ILLUMINA_CORRECTORS) *
    length(D5_ILLUMINA_ORGANISMS) * length(D5_ILLUMINA_SEEDS)
end

# === Statuses ===

const STATUS_OK = "ok"
const STATUS_ERROR = "error"
const STATUS_QUAST_FAILED = "quast_failed"
const STATUS_INFEASIBLE_K = "infeasible_k_gt_assembler_max"
const STATUS_NO_FEASIBLE_K = "no_feasible_k"
const NO_FEASIBLE_K = :NoFeasibleK
const RETRYABLE_STATUSES = (STATUS_ERROR, STATUS_QUAST_FAILED)

const D5_STRATEGY = :scalable
const D5_COLLISION_P = 0.01

# === Enumeration ===

"""
    ont_cells() -> Vector{NamedTuple}

Every registered ONT cell (grid over all 8 k-policies, then validation), built
with the registry's own constructors. `length(ont_cells()) == n_d5_cells()`.
"""
function ont_cells()
    cells = NamedTuple[]
    for identity in D5OntAdaptiveGridAxes.D5_IDENTITIES, coverage in D5OntAdaptiveGridAxes.D5_COVERAGES,
        corrector in D5OntAdaptiveGridAxes.D5_CORRECTORS, k_policy in D5OntAdaptiveGridAxes.D5_K_POLICIES,
        organism in D5OntAdaptiveGridAxes.D5_ORGANISMS, seed in D5OntAdaptiveGridAxes.D5_SEEDS
        push!(cells,
            merge((; technology = :ont),
                D5OntAdaptiveGridAxes.ont_grid_cell(identity, coverage, corrector, k_policy, organism, seed)))
    end
    for rung in D5OntAdaptiveGridAxes.D5_DROPPED_RUNGS, organism in D5OntAdaptiveGridAxes.D5_ORGANISMS, seed in D5OntAdaptiveGridAxes.D5_SEEDS
        push!(cells, merge((; technology = :ont), D5OntAdaptiveGridAxes.validation_cell(rung, organism, seed)))
    end
    return cells
end

"""
    illumina_cells() -> Vector{NamedTuple}

The Illumina control: coverage x corrector x (adaptive + 7 comparator rungs) x
organism x seed. 336 fixed-rung (Phase A) + 48 adaptive (Phase B).
"""
function illumina_cells()
    cells = NamedTuple[]
    for coverage in D5_ILLUMINA_COVERAGES, corrector in D5_ILLUMINA_CORRECTORS,
        k_policy in D5_ILLUMINA_K_POLICIES, organism in D5_ILLUMINA_ORGANISMS,
        seed in D5_ILLUMINA_SEEDS
        push!(cells,
            (; technology = :illumina, identity = nothing, coverage, corrector,
                k_policy, organism, seed, stratum = :illumina_control))
    end
    return cells
end

all_cells() = vcat(ont_cells(), illumina_cells())

"""Phase of a cell: `:B` for `:adaptive` (needs tau), `:A` for a fixed rung."""
cell_phase(cell) = cell.k_policy === :adaptive ? :B : :A

phase_a_cells() = filter(c -> cell_phase(c) === :A, all_cells())
phase_b_cells() = filter(c -> cell_phase(c) === :B, all_cells())

"""The read set a cell draws from. Survival is keyed on exactly these fields."""
function read_set_of(cell)
    (; technology = cell.technology, identity = cell.identity,
        coverage = cell.coverage, organism = cell.organism, seed = cell.seed)
end

"""
    read_sets() -> Vector{NamedTuple}

Every distinct read set, in a fixed order (ONT then Illumina). The SLURM array
index i (0-based) maps to `read_sets()[i + 1]`.
"""
read_sets() = unique(read_set_of.(all_cells()))

function cells_for_read_set(read_set; phase::Symbol)
    phase in (:A, :B) || error("phase must be :A or :B, got :$(phase)")
    return filter(c -> read_set_of(c) == read_set && cell_phase(c) === phase, all_cells())
end

function read_set_id(rs)
    base = rs.technology === :ont ? "ont__id$(rs.identity)" : "illumina"
    return "$(base)__cov$(rs.coverage)x__$(rs.organism)__seed$(rs.seed)"
end

function cell_id(cell; strategy::Symbol = D5_STRATEGY)
    k = cell.k_policy === :adaptive ? "adaptive" : "k$(cell.k_policy)"
    return "$(read_set_id(read_set_of(cell)))__$(cell.corrector)__$(k)__" *
           "$(cell.stratum)__$(strategy)"
end

"""
    toy_read_sets() / toy_cells(read_set)

The pre-declared toy subset (bead td-yddpr.5 design note): Lambda, first seed;
ONT 95%/10x, ONT 99%/100x (validation stratum), Illumina 30x. Cells are SELECTED
from the registered enumeration, so the toy can never contain an unregistered
cell; the toy tree is written to its own output directory.
"""
function toy_read_sets()
    seed = first(D5OntAdaptiveGridAxes.D5_SEEDS)
    keep(rs) = rs.organism === :Lambda && rs.seed == seed &&
               (
                   (rs.technology === :ont && rs.identity == 95 && rs.coverage == 10) ||
                   (rs.technology === :ont && rs.identity == 99 && rs.coverage == 100) ||
                   (rs.technology === :illumina && rs.coverage == 30))
    return filter(keep, read_sets())
end

function toy_cells(read_set)
    keep(c) =
        if c.technology === :ont && c.identity == 95
            c.stratum === :grid && c.k_policy in (21, 31)
        elseif c.technology === :ont && c.identity == 99
            c.stratum === :validation && c.k_policy == 29
        else
            c.k_policy == 21
        end
    return filter(keep, cells_for_read_set(read_set; phase = :A))
end

# === Survival ===

"""
    kmer_survival(reads, reference, k) -> NamedTuple

Measured survival at one rung (DEC-2026-08-09). `reads` and `reference` are
vectors of `BioSequences.LongDNA`. Canonical k-mers; windows containing an
ambiguous base are skipped on both sides.

  * `S` = observed / `n_ref_positions`, where a reference position is observed
    when its canonical k-mer occurs exactly (error-free) in at least one read.
  * `histogram` maps exact multiplicity -> number of DISTINCT read k-mers.
  * `S_hat`, `s_hat_v`, `s_hat_m` from `s_hat(histogram)`; `nothing` when the
    surrogate is undefined at this rung (`s_hat_status == "undefined"`).
"""
function kmer_survival(reads::AbstractVector{<:BioSequences.LongDNA},
        reference::AbstractVector{<:BioSequences.LongDNA}, k::Integer)
    k >= 1 || error("k must be positive, got k=$(k)")
    T = Kmers.derive_type(Kmers.Kmer{BioSequences.DNAAlphabet{2}, Int(k)})
    return _kmer_survival(T, Val(Int(k)), reads, reference)
end

function _kmer_survival(::Type{T}, ::Val{K}, reads, reference) where {T, K}
    counts = Dict{T, Int}()
    for read in reads
        for (kmer, _) in Kmers.UnambiguousDNAMers{K}(read)
            canonical = BioSequences.canonical(kmer)
            counts[canonical] = get(counts, canonical, 0) + 1
        end
    end
    n_ref_positions = 0
    n_observed = 0
    for contig in reference
        for (kmer, _) in Kmers.UnambiguousDNAMers{K}(contig)
            n_ref_positions += 1
            haskey(counts, BioSequences.canonical(kmer)) && (n_observed += 1)
        end
    end
    histogram = Dict{Int, Int}()
    for multiplicity in values(counts)
        histogram[multiplicity] = get(histogram, multiplicity, 0) + 1
    end
    surrogate = s_hat(histogram)
    return (; k = K, n_ref_positions, n_ref_positions_observed = n_observed,
        S = n_ref_positions == 0 ? nothing : n_observed / n_ref_positions,
        n_distinct_read_kmers = length(counts), histogram,
        S_hat = surrogate.s_hat, s_hat_v = surrogate.v, s_hat_m = surrogate.m,
        s_hat_status = surrogate.v === nothing ? "undefined" : "defined")
end

"""
    survival_ladder(reads, reference; ladder = D5_CANDIDATE_LADDER) -> Vector

One `kmer_survival` row per rung, in ladder order. Rungs run sequentially so
only one rung's k-mer table is resident at a time.
"""
function survival_ladder(reads, reference;
        ladder = D5OntAdaptiveGridAxes.D5_CANDIDATE_LADDER)
    return [kmer_survival(reads, reference, k) for k in ladder]
end

"""
    s_hat(counts) -> (; v, m, s_hat)

DEC-2026-08-09 criterion 5, ported exactly from rhizomorph-paper
decisions/surrogate_spectrum.py (f3a5ed5). `counts` maps exact multiplicity ->
number of distinct read k-mers; absent keys are zero.

  1. Valley `v`: the smallest i >= 2 with c(i) < c(i+1), STRICT.
  2. Mode `m`: the smallest maximiser of c over multiplicities strictly above v.
  3. Undefined unless c(m) > c(v). S_hat = fraction of distinct k-mers with
     multiplicity >= v.

Undefined returns all three fields as `nothing`.
"""
function s_hat(counts::AbstractDict{<:Integer, <:Integer})
    undefined = (; v = nothing, m = nothing, s_hat = nothing)
    isempty(counts) && return undefined
    c(i) = get(counts, i, 0)
    hi = maximum(keys(counts))
    v = nothing
    for i in 2:(hi - 1)
        if c(i) < c(i + 1)
            v = i
            break
        end
    end
    v === nothing && return undefined
    above = [(j, n) for (j, n) in counts if j > v]
    isempty(above) && return undefined
    peak = maximum(last, above)
    m = minimum(first(p) for p in above if last(p) == peak)
    c(m) <= c(v) && return undefined
    total = sum(values(counts))
    return (; v, m, s_hat = sum(n for (i, n) in counts if i >= v) / total)
end

# === k* lookup (Phase B) ===

"""Collision floor ceil(log_|A|(G / p)) from DEC-2026-08-09 (|A| = 4, p = 0.01)."""
function collision_floor_k(genome_size::Integer; alphabet_size::Integer = 4,
        p::Real = D5_COLLISION_P)
    ceil(Int, log(alphabet_size, genome_size / p))
end

"""
    select_k_star(values, tau; ladder, genome_size) -> (; k_star, reason, floor_k)

The registered level-set rule. `values[i]` is the survival (S or S_hat) at
`ladder[i]`; `nothing` / NaN entries (an undefined S_hat rung) are dropped.
Returns the largest rung at or above the collision floor whose value is >= tau,
or `NO_FEASIBLE_K` with the reason. Never falls back to the floor, the ceiling,
or an error-rate-derived k. Takes no error-rate input by construction.
"""
function select_k_star(values::AbstractVector, tau::Real;
        ladder = D5OntAdaptiveGridAxes.D5_CANDIDATE_LADDER, genome_size::Integer)
    length(values) == length(ladder) ||
        error("survival vector has $(length(values)) entries for a " *
              "$(length(ladder))-rung ladder")
    floor_k = collision_floor_k(genome_size)
    above_floor = [i for i in eachindex(ladder) if ladder[i] >= floor_k]
    isempty(above_floor) &&
        return (; k_star = NO_FEASIBLE_K, reason = "no_rung_above_collision_floor", floor_k)
    admissible = [ladder[i]
                  for i in above_floor
                  if values[i] !== nothing && !isnan(values[i]) && values[i] >= tau]
    isempty(admissible) &&
        return (; k_star = NO_FEASIBLE_K, reason = "no_rung_meets_threshold", floor_k)
    return (; k_star = maximum(admissible), reason = "admissible", floor_k)
end

"""
    assembler_accepts_k(k) -> Bool

Asks the assembler itself rather than copying its limit: true when
`AssemblyConfig` accepts `k`. Today that is k <= 64 (assembly.jl); when the cap
is lifted (td-yddpr.6) the driver follows without an edit.
"""
function assembler_accepts_k(k::Integer)
    try
        Mycelia.Rhizomorph.AssemblyConfig(; k = Int(k))
        return true
    catch e
        occursin("k-mer size must be between", sprint(showerror, e)) && return false
        rethrow()
    end
end

# === Simulation parameters (recorded literally) ===

"""The Badread argv for an ONT read set, exactly as the simulator will build it."""
function badread_argv(read_set, reference_fasta::AbstractString)
    read_set.technology === :ont || error("badread_argv is for ONT read sets only")
    return Mycelia._badread_nanopore_args(; fasta = reference_fasta,
        quantity = "$(read_set.coverage)x",
        error_model = D5_BADREAD_ERROR_MODEL, qscore_model = D5_BADREAD_QSCORE_MODEL,
        identity = D5_BADREAD_IDENTITY_PARAMS[read_set.identity],
        length_dist = D5_BADREAD_LENGTH_DIST, seed = read_set.seed)
end

sequencing_tech_for(read_set) = read_set.technology === :ont ? :nanopore : :illumina

# === Checkpoint I/O ===

"""Atomic JSON write: temp file in the same directory, then rename."""
function save_json_atomically(path::AbstractString, data)
    mkpath(dirname(path))
    tmp = tempname(dirname(path))
    try
        open(tmp, "w") do io
            JSON.print(io, data, 2)
        end
        mv(tmp, path; force = true)
    catch
        isfile(tmp) && rm(tmp; force = true)
        rethrow()
    end
    return path
end

"""
JSON-ready value: `missing` and `nothing` become `nothing` (JSON null). A
censored metric is null on disk, never 0.0.
"""
jsonable(x) = (x === missing || x === nothing) ? nothing : x
jsonable(x::Symbol) = String(x)

row_to_json(row::NamedTuple) = Dict(String(k) => jsonable(v) for (k, v) in pairs(row))

"""Read a checkpoint back with nulls as `missing`."""
function load_row(path::AbstractString)
    raw = JSON.parsefile(path)
    return Dict(k => (v === nothing ? missing : v) for (k, v) in raw)
end

# === Per-cell row ===

const METRIC_KEYS = (:NGA50, :genome_fraction, :misassemblies, :duplication_ratio,
    :largest_contig, :N50, :total_length, :quast_contigs)

"""
    cell_row(cell, context; kwargs...) -> NamedTuple

The per-cell record. `context` carries read-set and run provenance. Metrics come
from `parse_quast_metrics` and stay `missing` when QUAST could not compute them.
"""
function cell_row(cell, context; phase::Symbol, status::AbstractString,
        k_requested = cell.k_policy === :adaptive ? missing : cell.k_policy,
        k_selected = k_requested, k_star_reason = missing, k_star_hat = missing,
        reassembly_k = missing, reassembly_k_policy = missing,
        reused_from_cell_id = missing, n_contigs = missing, contigs_ge_min = missing,
        quast_min_contig = missing, metrics = empty_metrics(), nga50_status = missing,
        wall_seconds = missing, peak_rss_bytes = missing, peak_rss_method = missing,
        error_message = missing)
    rs = read_set_of(cell)
    metric_values = NamedTuple{METRIC_KEYS}(Tuple(getproperty(metrics, k)
    for k in METRIC_KEYS))
    return (;
        cell_id = cell_id(cell), phase = String(phase), stratum = String(cell.stratum),
        technology = String(cell.technology), organism = String(cell.organism),
        accession = D5_ORGANISM_REFERENCES[cell.organism].accession,
        seed = cell.seed, coverage = cell.coverage,
        identity_target = something(cell.identity, missing),
        identity_params = context.identity_params, simulator_argv = context.simulator_argv,
        error_model = context.error_model, qscore_model = context.qscore_model,
        read_set_id = read_set_id(rs), corrector = String(cell.corrector),
        strategy = String(D5_STRATEGY),
        k_policy = cell.k_policy === :adaptive ? "adaptive" : "fixed",
        k_requested, k_selected, k_star_reason, k_star_hat, reassembly_k,
        reassembly_k_policy, reused_from_cell_id,
        n_reads = context.n_reads, read_bases = context.read_bases,
        n_contigs, contigs_ge_min, quast_min_contig, metric_values...,
        nga50_status, status = String(status),
        error_rate_estimate_diagnostic = context.error_rate_estimate_diagnostic,
        wall_seconds, peak_rss_bytes, peak_rss_method, error_message,
        mycelia_git_sha = context.mycelia_git_sha, git_dirty = context.git_dirty,
        host = context.host, slurm_job_id = context.slurm_job_id,
        julia_version = string(VERSION))
end

"""Row for a cell whose k the assembler cannot take (k > 64 today)."""
function infeasible_row(cell, context; phase::Symbol, k_selected = cell.k_policy, kwargs...)
    cell_row(cell, context; phase, status = STATUS_INFEASIBLE_K, k_selected, kwargs...)
end

# === CLI parsing (pure, so the tau requirement is testable) ===

"""
    parse_cli(args) -> NamedTuple

`--phase A --read-set-index i` | `--phase B --tau x --tau-hat y` |
`--phase merge` | `--list`; plus `--output-dir`, `--toy`, `--allow-sha-mismatch`.
Phase B refuses to run without tau: the grid never chooses tau itself.
"""
function parse_cli(args::AbstractVector{<:AbstractString})
    value(flag) =
        let i = findfirst(==(flag), args)
            i === nothing ? nothing :
            i == length(args) ? error("$(flag) requires a value") : args[i + 1]
        end
    list = "--list" in args
    phase_raw = value("--phase")
    phase = phase_raw === nothing ? (list ? :list : error("--phase is required")) :
            Symbol(phase_raw)
    phase in (:A, :B, :merge, :list) ||
        error("--phase must be A, B or merge, got $(repr(phase_raw))")
    index_raw = value("--read-set-index")
    read_set_index = index_raw === nothing ? nothing : parse(Int, index_raw)
    tau_raw = value("--tau")
    tau_hat_raw = value("--tau-hat")
    if phase === :B && (tau_raw === nothing || tau_hat_raw === nothing)
        error("--phase B requires both --tau and --tau-hat. tau and tau_hat are set " *
              "by a human from the Illumina :none stratum (lock step 2); the driver " *
              "never chooses them.")
    end
    tau = tau_raw === nothing ? nothing : parse(Float64, tau_raw)
    tau_hat = tau_hat_raw === nothing ? nothing : parse(Float64, tau_hat_raw)
    for (name, t) in (("--tau", tau), ("--tau-hat", tau_hat))
        t === nothing || (0.0 <= t <= 1.0) || error("$(name) must be in [0, 1], got $(t)")
    end
    return (; phase, read_set_index, tau, tau_hat, toy = "--toy" in args,
        output_dir = value("--output-dir"),
        allow_sha_mismatch = "--allow-sha-mismatch" in args)
end

end # module
