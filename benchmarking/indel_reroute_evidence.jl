# Indel reroute evidence harness. Written for td-4mbg (does gating complete-span
# windowing on a read-level indel predicate change end-to-end correction
# quality?) and reused for td-2rxh (does :nanopore beat :illumina at 1000bp+?).
#
# WHY THIS FILE EXISTS. The first version of this evidence was produced by an
# uncommitted scratch script, so nothing in the repo reproduced the numbers that
# argued for the change. A reviewer could not tell which read generator had been
# used -- and the two committed candidates differ in whether the reads contain
# indels at all, which is decisive for an indel-decode fix. This harness makes
# the generator, the cells, the seeds and the metrics explicit and re-runnable.
#
# WHAT IT MEASURES. For each (genome, read length, coverage, error rate, seed)
# cell it simulates ONE read set and runs it through `assemble_genome` twice,
# varying ONLY `sequencing_tech`:
#
#   illumina  -- substitution-only profile; builds no window-source map, so it is
#                the control arm and the baseline this change tries to match
#   nanopore  -- indel profile; under `:frontier_budgeted` its scheduler emits a
#                window-source map, which is what the change reroutes
#
# Reads are simulated with `Mycelia.observe(...; tech = :nanopore)`, which per
# `src/simulation.jl` is 40% mismatch / 30% insertion / 30% deletion -- so the
# reads DO carry indels and the indel arm has something real to recover. Both
# arms receive byte-identical reads.
#
# MULTI-SEED BY DESIGN. Largest-contig is a discontinuous order statistic: one
# correction can join or split a contig, so a single seed cannot separate an
# effect from fixture noise. Every cell is replicated across seeds and the
# summary reports the spread, not just a point.
#
# Run against a specific checkout (this is how the before/after comparison is
# made -- point it at an unfixed checkout, then at the branch):
#
#   MYCELIA_PROJECT=/path/to/checkout LD_LIBRARY_PATH='' \
#     julia --startup-file=no --project=/path/to/checkout \
#       benchmarking/indel_reroute_evidence.jl --seeds 42,43,44 --out results.csv
#
# MYCELIA_PROJECT (default: this file's checkout) is what selects the code: the
# script activates it, which supersedes --project. Before running anything it
# checks the resulting state (`check_provenance`): Mycelia must load from
# <checkout>/src/Mycelia.jl and the checkout must be first on the load path, so
# no dependency can resolve from another environment. --startup-file=no stops a
# depot startup.jl from editing the load path; the check refuses if one did.
#
# Options:
#   --seeds a,b,c    replicate seeds            (default 42,43,44)
#   --cells SET      cell set: all (default, the three td-4mbg cells), smoke
#                    (the two 2kb cells) or td2rxh (the td-2rxh 1000bp+ grid)
#   --out PATH       CSV output path            (default indel_reroute_evidence.csv)
#   --label NAME     provenance tag written into every row (e.g. "master", "fix")
#   --warmup BOOL    untimed warm-up of both arms before timing: true (default)
#                    or false; any other value is an error
#
# td-2rxh grid (`--cells td2rxh`): 6 kb genome x read length {1000, 2000} x 20x
# x error {0.01, 0.02, 0.05}. It asks whether the merged sparse-rung staging
# (#420) plus the #450 reroute fix deliver the :nanopore > :illumina win at
# 1000bp+. Its decision rule (and amendments 1-3) was pre-registered on bead
# td-2rxh before any run; benchmarking/td2rxh_paired_analysis.jl applies it.
#
# Timing fairness: both arms of a pair run in one process; `GC.gc()` runs before
# each timed arm so neither inherits the other's garbage, and the arm order
# alternates with seed parity so any order effect is balanced across seeds (8 vs 8).

import Pkg
const _PROJECT = get(ENV, "MYCELIA_PROJECT", dirname(@__DIR__))
# Only re-activate when run as a script: `include`-ing this file from a test must
# not swap the test environment out from under the test runner.
if abspath(PROGRAM_FILE) == @__FILE__
    Pkg.activate(_PROJECT; io = devnull)
end

import Mycelia
import FASTX
import BioSequences
import Random
import Statistics
import Dates

# `indel_bench_best_reference_alignment`: best-contig reference COVERAGE
# (contiguity) and FIT identity (accuracy). Its docstring explains why the two
# must be reported together and why neither alone is an accuracy claim.
include(joinpath(@__DIR__, "indel_benchmark_common.jl"))

"""
One measurement cell. `error_rate` is the per-base error rate handed to
`Mycelia.observe`; `kmax` is the assembler's top k.
"""
struct EvidenceCell
    genome_length::Int
    read_length::Int
    coverage::Int
    error_rate::Float64
    kmax::Int
end

const ALL_CELLS = [
    EvidenceCell(2000, 1200, 8, 0.01, 31),
    EvidenceCell(2000, 1200, 20, 0.01, 31),
    EvidenceCell(6000, 1200, 20, 0.01, 31)
]
const SMOKE_CELLS = ALL_CELLS[1:2]
const TD2RXH_CELLS = [EvidenceCell(6000, read_length, 20, error_rate, 31)
                      for read_length in (1000, 2000)
                      for error_rate in (0.01, 0.02, 0.05)]

"""
Resolve a `--cells` value to its cell set. An unrecognised name throws instead
of falling back to `ALL_CELLS`: a typo in a launcher must not quietly run a
different grid than the one the pre-registration names.
"""
function select_cells(name::AbstractString)::Vector{EvidenceCell}
    name == "all" && return ALL_CELLS
    name == "smoke" && return SMOKE_CELLS
    name == "td2rxh" && return TD2RXH_CELLS
    throw(ArgumentError("unknown --cells \"$(name)\"; expected all, smoke or td2rxh"))
end

cell_id(cell::EvidenceCell)::String = string(
    cell.genome_length, "/", cell.read_length, "/", cell.coverage, "x/e",
    cell.error_rate)

"""
Simulate one read set. Mirrors `indel_bench_simulate_reads` /
`indel_toy_overcorrection_diagnostic.jl` so the fixture matches the committed
harnesses rather than inventing a third convention.

Returns the reads and the reference, so downstream identity metrics have a
ground truth to align against.
"""
function simulate_reads(cell::EvidenceCell, seed::Int)
    reference_record = Mycelia.random_fasta_record(
        moltype = :DNA, seed = seed, L = cell.genome_length)
    reference = FASTX.sequence(BioSequences.LongDNA{4}, reference_record)
    rng = Random.MersenneTwister(seed)
    Random.seed!(seed)
    n_reads = ceil(Int, cell.coverage * cell.genome_length / cell.read_length)
    reads = FASTX.FASTQ.Record[]
    for read_index in 1:n_reads
        start_position = rand(
            rng, 1:(cell.genome_length - cell.read_length + 1))
        fragment = reference[start_position:(start_position + cell.read_length - 1)]
        if rand(rng, Bool)
            fragment = BioSequences.reverse_complement(fragment)
        end
        observed,
        qualities = Mycelia.observe(
            fragment; error_rate = cell.error_rate, tech = :nanopore)
        isempty(observed) && continue
        quality_string = String([Char(quality + 33) for quality in qualities])
        push!(reads, FASTX.FASTQ.Record(
            "read_$(read_index)", string(observed), quality_string))
    end
    return reads, reference
end

"""
Run one arm. `tech` selects the corrector profile; everything else is held fixed.
The RNG is reseeded identically per arm so the arms differ only in `tech`.

`wall_seconds_per_read` is END-TO-END wall time over the read count, graph
construction included. It is not decode time on its own; the indel profile's
incremental cost is the paired nanopore-minus-illumina difference on the same
reads. The alignment metrics run after the timer stops, so they do not inflate it.

An empty assembly reports `best_contig_fit_identity = NaN`: there is no contig
whose accuracy could be measured, and the alignment helper's `0.0` sentinel would
otherwise score as a maximally wrong contig. Its reference coverage stays `0.0`,
which is a real contiguity outcome.
"""
function run_arm(reads::Vector{FASTX.FASTQ.Record}, cell::EvidenceCell,
        tech::Symbol, reference::BioSequences.LongDNA{4})
    Random.seed!(1_042)
    GC.gc()
    wall_seconds = @elapsed assembly = Mycelia.Rhizomorph.assemble_genome(
        deepcopy(reads);
        k = cell.kmax,
        corrector = :iterative,
        strategy = :scalable,
        sequencing_tech = tech)
    stats = assembly.assembly_stats
    lengths = [length(contig) for contig in assembly.contigs]
    largest = isempty(lengths) ? 0 : maximum(lengths)
    alignment = indel_bench_best_reference_alignment(
        String[string(contig) for contig in assembly.contigs], reference)
    return (;
        substitution_length_divergences = get(
            stats, "substitution_length_divergences", missing),
        window_divergences = get(stats, "window_divergences", missing),
        window_anchor_rejections = get(
            stats, "window_anchor_rejections", missing),
        trace_contract_errors = get(stats, "trace_contract_errors", missing),
        indel_engaged = get(stats, "indel_engaged", missing),
        indel_requested = get(stats, "indel_requested", missing),
        reassembly_k = get(stats, "reassembly_k", missing),
        n_contigs = length(assembly.contigs),
        largest_contig = largest,
        largest_ratio = largest / cell.genome_length,
        total_contig_bases = sum(lengths; init = 0),
        wall_seconds = wall_seconds,
        wall_seconds_per_read = per_read_seconds(wall_seconds, length(reads)),
        best_contig_reference_coverage = alignment.best_contig_reference_coverage,
        best_contig_fit_identity = isempty(assembly.contigs) ? NaN :
                                   alignment.best_contig_fit_identity)
end

"""
Arm order for one seed: illumina first on odd seeds, nanopore first on even
seeds, so any order effect on the paired wall-time difference is balanced across seeds.
"""
arm_order(seed::Integer)::Tuple{
    Symbol, Symbol} = iseven(seed) ? (:nanopore, :illumina) : (:illumina, :nanopore)

"""Parse a strict boolean flag value; anything but "true"/"false" is an error."""
function parse_bool_flag(name::AbstractString, value::AbstractString)::Bool
    value == "true" && return true
    value == "false" && return false
    throw(ArgumentError("--$(name) must be true or false, got \"$(value)\""))
end

"""
Wall seconds per read, or `NaN` for an empty read set. `NaN` keeps a degenerate
cell visible in the CSV; a `0` here would look like a very fast arm.
"""
per_read_seconds(wall_seconds::Real,
    n_reads::Integer)::Float64 = n_reads == 0 ? NaN : Float64(wall_seconds) / n_reads

"""
Parse `--key value` pairs. Unknown keys and a dangling flag are errors, so a typo
(`--cell td2rxh`, `--seed 3`) cannot silently run the default grid or seeds.
"""
function parse_args(args::Vector{String})::Dict{String, String}
    allowed = ("seeds", "cells", "out", "label", "warmup")
    isodd(length(args)) && throw(ArgumentError("expected --key value pairs, got $(args)"))
    parsed = Dict{String, String}()
    for index in 1:2:length(args)
        startswith(args[index], "--") ||
            throw(ArgumentError("expected a --flag, got \"$(args[index])\""))
        key = args[index][3:end]
        key in allowed || throw(ArgumentError("unknown flag --$(key)"))
        startswith(args[index + 1], "--") &&
            throw(ArgumentError("--$(key) is missing its value (got $(args[index + 1]))"))
        parsed[key] = args[index + 1]
    end
    return parsed
end

# New columns are APPENDED so CSVs written before td-2rxh keep a prefix-compatible
# header and still parse positionally.
const CSV_HEADER = [
    "label", "cell", "seed", "arm", "n_reads",
    "substitution_length_divergences", "window_divergences",
    "window_anchor_rejections", "trace_contract_errors",
    "indel_engaged", "indel_requested", "reassembly_k",
    "n_contigs", "largest_contig", "largest_ratio",
    "total_contig_bases", "wall_seconds",
    "wall_seconds_per_read", "best_contig_reference_coverage",
    "best_contig_fit_identity"]

"""
Format one arm's result as a CSV row in `CSV_HEADER` order.

Precision: whole-arm `wall_seconds` keeps 3 decimals (CSVs written before td-2rxh
carry 1), `wall_seconds_per_read` 4, and both alignment metrics 6, so the
td-2rxh analysis can round paired deltas to 1e-6 without the CSV having already
quantised them. An absent telemetry counter is written as an empty field (CSV
`missing`), not the literal string "missing".
"""
function csv_row(label::AbstractString, cell::EvidenceCell, seed::Integer,
        arm::Symbol, n_reads::Integer, result::NamedTuple)::String
    return join(
        _csv_field.([
            label, cell_id(cell), seed, arm, n_reads,
            result.substitution_length_divergences,
            result.window_divergences,
            result.window_anchor_rejections,
            result.trace_contract_errors,
            result.indel_engaged, result.indel_requested,
            result.reassembly_k, result.n_contigs,
            result.largest_contig,
            round(result.largest_ratio, digits = 4),
            result.total_contig_bases,
            round(result.wall_seconds, digits = 3),
            round(result.wall_seconds_per_read, digits = 4),
            round(result.best_contig_reference_coverage, digits = 6),
            round(result.best_contig_fit_identity, digits = 6)]),
        ",")
end

_csv_field(value) = value === missing ? "" : string(value)

function write_rows(output_path::AbstractString, rows::Vector{String})::Nothing
    open(output_path, "w") do io
        for row in rows
            println(io, row)
        end
    end
    return nothing
end

"""
    main(args; runner = run_arm) -> Int

Run every (cell, seed, arm) and return the number of arms that threw. A failed
arm is logged with its backtrace and skipped, so the other arms still produce
rows; `script_main` turns any failure into exit 1, so a short CSV can never come
with exit 0. `runner` exists so tests can force failures without running an
assembly. The warm-up always calls `run_arm` directly (an error there is
uncaught and aborts the run), so tests that pass a `runner` use `--warmup false`.
"""
function main(args::Vector{String} = ARGS; runner = run_arm)::Int
    parsed = parse_args(args)
    seeds = [parse(Int, s) for s in split(get(parsed, "seeds", "42,43,44"), ",")]
    cells_name = get(parsed, "cells", "all")
    cells = select_cells(cells_name)
    output_path = get(parsed, "out", "indel_reroute_evidence.csv")
    label = get(parsed, "label", "unlabeled")
    warmup = parse_bool_flag("warmup", get(parsed, "warmup", "true"))

    println("# indel reroute evidence cells=", cells_name)
    println("# label=", label, " project=", _PROJECT)
    println("# seeds=", join(seeds, ","), " cells=", length(cells), " warmup=", warmup)
    println("# started=", string(Dates.now()))
    flush(stdout)

    # Write the header immediately: a re-run into an existing path must never
    # leave the previous run's rows in place if every arm of this run fails.
    rows = String[join(CSV_HEADER, ",")]
    write_rows(output_path, rows)

    # Untimed warm-up of BOTH arms. Without it, whichever (cell, seed) runs first
    # in the process -- under the one-seed-per-task sbatch, always the first
    # td2rxh cell -- pays for JIT compilation, biasing the paired wall-time
    # difference (the H1' quantity) there. The warm-up uses a small smoke cell; it
    # is not proven to compile every path a 6 kb cell reaches. Results discarded.
    if warmup
        warm_cell = SMOKE_CELLS[1]
        warm_reads, warm_reference = simulate_reads(warm_cell, 0)
        for arm in (:illumina, :nanopore)
            warm_seconds = @elapsed run_arm(warm_reads, warm_cell, arm, warm_reference)
            println("# warmup ", arm, " ", round(warm_seconds, digits = 1), "s (discarded)")
        end
        flush(stdout)
    end

    n_failed = 0
    for cell in cells
        for seed in seeds
            reads, reference = simulate_reads(cell, seed)
            for arm in arm_order(seed)
                result = try
                    runner(reads, cell, arm, reference)
                catch exception
                    n_failed += 1
                    @error "arm failed" cell=cell_id(cell) seed arm exception=(
                        exception, catch_backtrace())
                    continue
                end
                push!(rows, csv_row(label, cell, seed, arm, length(reads), result))
                println("  ", cell_id(cell), " seed=", seed, " ", arm,
                    " div=", result.substitution_length_divergences,
                    " anchor=", result.window_anchor_rejections,
                    " contigs=", result.n_contigs,
                    " largest=", result.largest_contig,
                    " (", round(result.largest_ratio, digits = 3), ")",
                    " cov=", round(result.best_contig_reference_coverage, digits = 3),
                    " fit_id=", round(result.best_contig_fit_identity, digits = 4),
                    " s/read=", round(result.wall_seconds_per_read, digits = 3))
                flush(stdout)
                write_rows(output_path, rows)
            end
        end
    end

    println("# finished=", string(Dates.now()))
    println("# wrote ", output_path, " (", length(rows) - 1, " rows)")
    n_failed > 0 &&
        println(stderr, "# FAILED: ", n_failed, " arm(s) threw; see the errors above")
    return n_failed
end

_resolved(path::AbstractString)::String = ispath(path) ? realpath(path) : String(path)

"""
    check_provenance(project; loaded = pathof(Mycelia), load_path = Base.load_path())
        -> String

Check the loaded STATE, not the environment variables that could have produced
it: Mycelia must be exactly `<project>/src/Mycelia.jl` (a sibling or nested
checkout fails), and `<project>/Project.toml` must be FIRST on the load path.
The project's Manifest lists every dependency, so with the project first each
one resolves from it; an earlier entry -- from `JULIA_LOAD_PATH`, a depot
`startup.jl`, or anything else -- could otherwise supply a dependency from
another environment while Mycelia itself still loads from the right place.
Paths are compared after `realpath`, so symlinks cannot fake or hide a match.

Returns the Mycelia path. Throws `ArgumentError` on a mismatch, and
`Base.IOError` (from `realpath`) when `project` does not exist. Not covered: a
precompile cache built from another checkout, and Manifest `path =` dependencies
(hashed by path, not content; the current Manifest has none).
"""
function check_provenance(project::AbstractString;
        loaded::AbstractString = pathof(Mycelia),
        load_path::AbstractVector{<:AbstractString} = Base.load_path())::String
    root = realpath(project)
    expected = joinpath(root, "src", "Mycelia.jl")
    actual = _resolved(loaded)
    actual == expected ||
        throw(ArgumentError("Mycelia was loaded from $(actual), not $(expected); " *
                            "check MYCELIA_PROJECT and the load path"))
    first_env = isempty(load_path) ? "<empty>" : _resolved(first(load_path))
    first_env == joinpath(root, "Project.toml") ||
        throw(ArgumentError("the project is not first on the load path (first: " *
                            "$(first_env)); dependencies could resolve from another " *
                            "environment. Run with --startup-file=no and without " *
                            "JULIA_LOAD_PATH"))
    return expected
end

"""
    script_main(args; runner = run_arm, provenance = () -> check_provenance(_PROJECT))
        -> Int

The script's exit status: 1 when any arm failed, else 0. The provenance check
runs FIRST, before any arm or output file, so a refused run produces nothing;
a refusal (or a warm-up error) throws rather than returning, and the uncaught
exception still makes julia exit non-zero. Tests inject `provenance`, since a
test environment is legitimately not the script's environment.
"""
function script_main(args::Vector{String} = ARGS; runner = run_arm,
        provenance = () -> check_provenance(_PROJECT))::Int
    println("# mycelia=", provenance())
    return main(args; runner = runner) == 0 ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(script_main())
end
