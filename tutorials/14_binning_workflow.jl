# # Tutorial 14: Binning Workflow with VAMB, MetaBAT2, and dRep
#
# This tutorial demonstrates a practical metagenomic binning workflow:
# 1. Prepare contigs and a coverage/depth table
# 2. Run `VAMB` and `MetaBAT2` on the same assembly
# 3. Inspect the resulting contig-to-bin assignments and bin FASTAs
# 4. Dereplicate candidate MAGs with `dRep`
#
# The tutorial is designed to be useful in two modes:
# - **Default mode**: lightweight parser demos and exact setup guidance
# - **External mode** (`MYCELIA_RUN_EXTERNAL=true`): run the full workflow
#
# If you do not already have contigs and a depth table, the external mode can
# generate reproducible synthetic inputs via `Mycelia.get_binning_test_inputs()`.
#
# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/14_binning_workflow.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import DataFrames

const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
const CONTIGS_FASTA_ENV = get(ENV, "MYCELIA_BINNING_CONTIGS", "")
const DEPTH_FILE_ENV = get(ENV, "MYCELIA_BINNING_DEPTH", "")
const GENOMES_ENV = filter(!isempty, strip.(split(get(ENV, "MYCELIA_BINNING_GENOMES", ""), ',')))

function _preview_lines(path::String; n::Int = 5)
    isfile(path) || return String[]
    lines = String[]
    open(path, "r") do io
        while !eof(io) && length(lines) < n
            push!(lines, readline(io))
        end
    end
    return lines
end

function _find_first_matching_file(dir::String, patterns::Vector{Regex})
    isdir(dir) || return nothing
    for (dirpath, _, filenames) in walkdir(dir)
        for name in sort(filenames)
            for pattern in patterns
                if occursin(pattern, name)
                    return joinpath(dirpath, name)
                end
            end
        end
    end
    return nothing
end

function _fasta_paths(dir::String)
    isdir(dir) || return String[]
    fasta_paths = String[]
    for name in sort(readdir(dir))
        path = joinpath(dir, name)
        isfile(path) || continue
        lower_name = lowercase(name)
        if endswith(lower_name, ".fa") ||
           endswith(lower_name, ".fna") ||
           endswith(lower_name, ".fasta")
            push!(fasta_paths, path)
        end
    end
    return fasta_paths
end

function _summarize_fasta_files(paths::Vector{String})
    rows = NamedTuple[]
    for path in paths
        record_count = 0
        total_bases = 0
        for record in Mycelia.open_fastx(path)
            record_count += 1
            total_bases += length(Mycelia.FASTX.sequence(record))
        end
        push!(rows, (
            genome = basename(path),
            records = record_count,
            total_bases = total_bases
        ))
    end

    if isempty(rows)
        return DataFrames.DataFrame(
            genome = String[],
            records = Int[],
            total_bases = Int[]
        )
    end

    return DataFrames.DataFrame(rows)
end

function _summarize_assignments(df::DataFrames.DataFrame)
    if DataFrames.nrow(df) == 0
        return DataFrames.DataFrame(bin = String[], n_contigs = Int[])
    end
    summary = DataFrames.combine(
        DataFrames.groupby(df, :bin),
        :contig => length => :n_contigs
    )
    DataFrames.sort!(summary, :n_contigs, rev = true)
    return summary
end

function _print_section(title::String)
    println()
    println(title)
    println(repeat("=", length(title)))
end

println("=== Binning Workflow Tutorial ===")
println("Workflow: contigs + depth -> VAMB / MetaBAT2 -> candidate MAGs -> dRep")

workdir = mktempdir()
println("Scratch workspace: ", workdir)

# ## Step 1: Learn the parser interfaces first
#
# The workflow wrappers emit files that you typically inspect or pass into later
# stages. These parsers are cheap to run and clarify what the downstream tables
# should look like before invoking external tools.

_print_section("Step 1: Parser utilities")

toy_assignments = joinpath(workdir, "toy_assignments.tsv")
open(toy_assignments, "w") do io
    write(io, "contig\tbin\tlength\n")
    write(io, "contig_001\tbin_A\t120000\n")
    write(io, "contig_002\tbin_A\t95000\n")
    write(io, "contig_003\tbin_B\t87000\n")
end

assignments_df = Mycelia.parse_bin_assignments(toy_assignments)
println("Parsed contig-to-bin assignments:")
println(assignments_df)
println("Assignment summary:")
println(_summarize_assignments(assignments_df))

toy_drep_clusters = joinpath(workdir, "toy_drep_clusters.csv")
open(toy_drep_clusters, "w") do io
    write(io, "genome,secondary_cluster,representative,ani\n")
    write(io, "metabat_bin_01.fna,1,metabat_bin_01.fna,0.999\n")
    write(io, "vamb_bin_02.fna,1,metabat_bin_01.fna,0.998\n")
    write(io, "metabat_bin_03.fna,2,metabat_bin_03.fna,0.995\n")
end

drep_clusters_df = Mycelia.parse_drep_clusters(toy_drep_clusters)
println("Parsed dRep clustering table:")
println(drep_clusters_df)

# ## Step 2: Resolve workflow inputs
#
# External mode can use either:
# - user-provided data from environment variables, or
# - reproducible synthetic inputs generated by Mycelia
#
# User-provided inputs:
# - `MYCELIA_BINNING_CONTIGS`
# - `MYCELIA_BINNING_DEPTH`
# - `MYCELIA_BINNING_GENOMES` (optional comma-separated FASTAs for dRep)

_print_section("Step 2: Resolve workflow inputs")

inputs = nothing
inputs_source = "none"

if RUN_EXTERNAL
    if isfile(CONTIGS_FASTA_ENV) && isfile(DEPTH_FILE_ENV)
        inputs = (
            contigs_fasta = CONTIGS_FASTA_ENV,
            depth_file = DEPTH_FILE_ENV,
            genomes = filter(isfile, copy(GENOMES_ENV))
        )
        inputs_source = "environment variables"
    else
        synthetic_inputs = Mycelia.get_binning_test_inputs()
        inputs = (
            contigs_fasta = synthetic_inputs.contigs_fasta,
            depth_file = synthetic_inputs.depth_file,
            genomes = synthetic_inputs.genomes
        )
        inputs_source = "Mycelia.get_binning_test_inputs()"
    end

    println("External execution enabled.")
    println("Input source: ", inputs_source)
    println("Contigs FASTA: ", inputs.contigs_fasta)
    println("Depth table:   ", inputs.depth_file)
    println("Depth preview:")
    for line in _preview_lines(inputs.depth_file; n = 4)
        println("  ", line)
    end
else
    println("External execution disabled.")
    println("Set MYCELIA_RUN_EXTERNAL=true to run VAMB, MetaBAT2, and dRep.")
    println("If you already have data, also set:")
    println("  MYCELIA_BINNING_CONTIGS=/path/to/contigs.fna")
    println("  MYCELIA_BINNING_DEPTH=/path/to/jgi_depth.tsv")
    println("  MYCELIA_BINNING_GENOMES=/path/to/bin1.fna,/path/to/bin2.fna")
    println("Otherwise this tutorial will auto-generate synthetic inputs in external mode.")
end

# ## Step 3: Run VAMB and inspect its outputs
#
# `run_vamb` accepts a standard contig FASTA plus a JGI depth table. Mycelia
# converts the JGI table to the abundance TSV format VAMB expects.

vamb_result = nothing
vamb_bin_fastas = String[]

if RUN_EXTERNAL && inputs !== nothing
    _print_section("Step 3: Run VAMB")

    vamb_outdir = joinpath(workdir, "vamb_out")
    vamb_result = Mycelia.run_vamb(
        contigs_fasta = inputs.contigs_fasta,
        depth_file = inputs.depth_file,
        outdir = vamb_outdir
    )

    println("VAMB output directory: ", vamb_result.outdir)
    println("VAMB clusters table:   ", vamb_result.clusters_tsv)

    if vamb_result.clusters_tsv !== nothing
        println("VAMB cluster preview:")
        for line in _preview_lines(vamb_result.clusters_tsv; n = 5)
            println("  ", line)
        end
    end

    vamb_bin_fastas = _fasta_paths(joinpath(vamb_result.outdir, "bins"))
    if isempty(vamb_bin_fastas)
        println("VAMB did not emit bin FASTAs in this run; continuing with any clusters table and MetaBAT2 bins.")
    else
        println("VAMB bin FASTA summary:")
        println(_summarize_fasta_files(vamb_bin_fastas))
    end
end

# ## Step 4: Run MetaBAT2 and inspect produced bins
#
# MetaBAT2 emits bin FASTAs using the prefix returned in `bins_prefix`. For
# post-binning work, those FASTAs are usually the direct handoff into QC or
# dereplication tools.

metabat_result = nothing
metabat_bin_fastas = String[]

if RUN_EXTERNAL && inputs !== nothing
    _print_section("Step 4: Run MetaBAT2")

    metabat_outdir = joinpath(workdir, "metabat2_out")
    metabat_result = Mycelia.run_metabat2(
        contigs_fasta = inputs.contigs_fasta,
        depth_file = inputs.depth_file,
        outdir = metabat_outdir
    )

    println("MetaBAT2 output directory: ", metabat_result.outdir)
    println("Bin prefix: ", metabat_result.bins_prefix)

    metabat_bin_fastas = _fasta_paths(metabat_result.outdir)
    if isempty(metabat_bin_fastas)
        println("MetaBAT2 did not emit bin FASTAs in this run.")
    else
        println("MetaBAT2 bin FASTA summary:")
        println(_summarize_fasta_files(metabat_bin_fastas))
    end
end

# ## Step 5: Assemble the dRep candidate set
#
# A common pattern is to dereplicate bins collected from multiple binners. This
# removes near-duplicate MAGs and chooses representatives before downstream
# refinement, QC, or catalog building.

drep_inputs = String[]

if RUN_EXTERNAL && inputs !== nothing
    _print_section("Step 5: Collect candidate MAGs for dRep")

    append!(drep_inputs, vamb_bin_fastas)
    append!(drep_inputs, metabat_bin_fastas)
    unique!(drep_inputs)

    if length(drep_inputs) < 2 && !isempty(inputs.genomes)
        println("Using tutorial genome FASTAs as dRep inputs because binner FASTAs are limited.")
        append!(drep_inputs, inputs.genomes)
        unique!(drep_inputs)
    end

    println("Candidate genomes for dereplication: ", length(drep_inputs))
    if isempty(drep_inputs)
        println("No candidate MAG FASTAs available yet.")
    else
        println(_summarize_fasta_files(drep_inputs))
    end
end

# ## Step 6: Dereplicate MAGs with dRep
#
# The simulated tutorial genomes do not ship with completeness/contamination
# metadata, so the workflow uses `--ignoreGenomeQuality` for the synthetic case.

if RUN_EXTERNAL && !isempty(drep_inputs)
    _print_section("Step 6: Run dRep")

    drep_outdir = joinpath(workdir, "drep_out")
    drep_extra_args = String[]
    if inputs_source == "Mycelia.get_binning_test_inputs()"
        drep_extra_args = ["--ignoreGenomeQuality", "-l", "1000"]
    end

    drep_result = Mycelia.run_drep_dereplicate(
        genomes = drep_inputs,
        outdir = drep_outdir,
        extra_args = drep_extra_args
    )

    println("dRep output directory: ", drep_result.outdir)
    println("Winning genomes table: ", drep_result.winning_genomes)

    if drep_result.winning_genomes !== nothing
        println("Winning genomes preview:")
        for line in _preview_lines(drep_result.winning_genomes; n = 5)
            println("  ", line)
        end
    end

    cluster_file = _find_first_matching_file(
        drep_result.outdir,
        [r"Cdb\.csv$", r"clusters.*\.csv$"]
    )
    if cluster_file !== nothing
        println("dRep cluster table: ", cluster_file)
        parsed_clusters = Mycelia.parse_drep_clusters(cluster_file)
        println(parsed_clusters)
    end
end

# ## Step 7: Interpret the workflow
#
# At this point you can branch in a few directions:
# - compare VAMB and MetaBAT2 bins against QC metrics
# - pass dereplicated MAGs into CheckM/BUSCO/GTDB-Tk style validation
# - keep both raw binning outputs and dRep representatives for provenance

_print_section("Step 7: Next steps")
println("Suggested downstream steps:")
println("  1. Validate bin quality with completeness/contamination tools.")
println("  2. Track provenance: which representative came from which binner.")
println("  3. Feed dereplicated MAGs into annotation, pangenome, or comparative workflows.")
println("Tutorial workspace: ", workdir)
