# td-4mbg evidence harness: does gating complete-span windowing on a read-level
# indel predicate change end-to-end correction quality?
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
#     julia --project=/path/to/checkout benchmarking/indel_reroute_evidence.jl \
#       --seeds 42,43,44 --out results.csv
#
# Options:
#   --seeds a,b,c    replicate seeds            (default 42,43,44)
#   --cells smoke    restrict to the 2kb cells  (default all three)
#   --out PATH       CSV output path            (default indel_reroute_evidence.csv)
#   --label NAME     provenance tag written into every row (e.g. "master", "fix")

import Pkg
const _PROJECT = get(ENV, "MYCELIA_PROJECT", dirname(@__DIR__))
Pkg.activate(_PROJECT; io = devnull)

import Mycelia
import FASTX
import BioSequences
import Random
import Statistics
import Dates

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
"""
function run_arm(reads::Vector{FASTX.FASTQ.Record}, cell::EvidenceCell,
        tech::Symbol)
    Random.seed!(1_042)
    wall_seconds = @elapsed assembly = Mycelia.Rhizomorph.assemble_genome(
        deepcopy(reads);
        k = cell.kmax,
        corrector = :iterative,
        strategy = :scalable,
        sequencing_tech = tech)
    stats = assembly.assembly_stats
    lengths = [length(contig) for contig in assembly.contigs]
    largest = isempty(lengths) ? 0 : maximum(lengths)
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
        wall_seconds = wall_seconds)
end

function parse_args(args::Vector{String})
    parsed = Dict{String, String}()
    index = 1
    while index <= length(args)
        argument = args[index]
        if startswith(argument, "--") && index < length(args)
            parsed[argument[3:end]] = args[index + 1]
            index += 2
        else
            index += 1
        end
    end
    return parsed
end

function main(args::Vector{String} = ARGS)::Nothing
    parsed = parse_args(args)
    seeds = [parse(Int, s) for s in split(get(parsed, "seeds", "42,43,44"), ",")]
    cells = get(parsed, "cells", "all") == "smoke" ? SMOKE_CELLS : ALL_CELLS
    output_path = get(parsed, "out", "indel_reroute_evidence.csv")
    label = get(parsed, "label", "unlabeled")

    println("# td-4mbg reroute evidence")
    println("# label=", label, " project=", _PROJECT)
    println("# seeds=", join(seeds, ","), " cells=", length(cells))
    println("# started=", string(Dates.now()))
    flush(stdout)

    rows = String[]
    push!(rows,
        join(
            [
                "label", "cell", "seed", "arm", "n_reads",
                "substitution_length_divergences", "window_divergences",
                "window_anchor_rejections", "trace_contract_errors",
                "indel_engaged", "indel_requested", "reassembly_k",
                "n_contigs", "largest_contig", "largest_ratio",
                "total_contig_bases", "wall_seconds"],
            ","))

    for cell in cells
        for seed in seeds
            reads, _reference = simulate_reads(cell, seed)
            for arm in (:illumina, :nanopore)
                result = try
                    run_arm(reads, cell, arm)
                catch exception
                    @warn "arm failed" cell=cell_id(cell) seed arm exception
                    continue
                end
                push!(rows,
                    join(
                        string.([
                            label, cell_id(cell), seed, arm, length(reads),
                            result.substitution_length_divergences,
                            result.window_divergences,
                            result.window_anchor_rejections,
                            result.trace_contract_errors,
                            result.indel_engaged, result.indel_requested,
                            result.reassembly_k, result.n_contigs,
                            result.largest_contig,
                            round(result.largest_ratio, digits = 4),
                            result.total_contig_bases,
                            round(result.wall_seconds, digits = 1)]),
                        ","))
                println("  ", cell_id(cell), " seed=", seed, " ", arm,
                    " div=", result.substitution_length_divergences,
                    " anchor=", result.window_anchor_rejections,
                    " contigs=", result.n_contigs,
                    " largest=", result.largest_contig,
                    " (", round(result.largest_ratio, digits = 3), ")")
                flush(stdout)
                open(output_path, "w") do io
                    for row in rows
                        println(io, row)
                    end
                end
            end
        end
    end

    println("# finished=", string(Dates.now()))
    println("# wrote ", output_path, " (", length(rows) - 1, " rows)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
