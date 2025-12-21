# # Tutorial 14: Binning and Post-binning Workflow
#
# This tutorial focuses on metagenomic binning and post-binning utilities.
# External tool execution is opt-in via environment variables.
#
# ## Setup
import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import DataFrames

# ## Configuration (opt-in external tools)
#
# To run external binners, set:
# - `MYCELIA_RUN_EXTERNAL=true`
# - `MYCELIA_BINNING_CONTIGS` (FASTA)
# - `MYCELIA_BINNING_DEPTH` (JGI depth table)
# - `MYCELIA_BINNING_COVERAGE_TABLE` (coverage table)
# - `MYCELIA_BINNING_MARKERS` (marker table for COMEBin)
# - `MYCELIA_BINNING_GRAPH` (assembly graph for MetaCoAG)
# - `MYCELIA_BINNING_MAPPING` (mapping/coverage file for MetaCoAG)
# - `MYCELIA_BINNING_GENOMES` (comma-separated FASTA list for dRep)
# - `MYCELIA_BINNING_BINS_DIRS` (comma-separated bin directories for MAGmax)
#
# These inputs are intentionally user-supplied to keep the tutorial lightweight.

const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
const contigs_fasta = get(ENV, "MYCELIA_BINNING_CONTIGS", "")
const depth_file = get(ENV, "MYCELIA_BINNING_DEPTH", "")
const coverage_table = get(ENV, "MYCELIA_BINNING_COVERAGE_TABLE", "")
const marker_file = get(ENV, "MYCELIA_BINNING_MARKERS", "")
const assembly_graph = get(ENV, "MYCELIA_BINNING_GRAPH", "")
const mapping_file = get(ENV, "MYCELIA_BINNING_MAPPING", "")
const genomes = filter(!isempty, strip.(split(get(ENV, "MYCELIA_BINNING_GENOMES", ""), ',')))
const bins_dirs = filter(!isempty, strip.(split(get(ENV, "MYCELIA_BINNING_BINS_DIRS", ""), ',')))

# ## Part 1: Parser Utilities (no external tools needed)
#
# The binning module includes parsers for common output tables.

println("=== Binning Tutorial ===")
println("Parser utilities (no external tools required)")

tmp_assignments = tempname()
open(tmp_assignments, "w") do io
    write(io, "contig\tbin\tlength\n")
    write(io, "ctg1\tbinA\t1500\n")
    write(io, "ctg2\tbinB\t2400\n")
end

assignments_df = Mycelia.parse_bin_assignments(tmp_assignments)
println("Parsed assignments:")
println(assignments_df)

tmp_drep = tempname()
open(tmp_drep, "w") do io
    write(io, "genome,secondary_cluster,representative,ani\n")
    write(io, "g1,1,g1,0.99\n")
    write(io, "g2,1,g1,0.99\n")
end

drep_df = Mycelia.parse_drep_clusters(tmp_drep)
println("Parsed dRep clusters:")
println(drep_df)

# ## Part 2: Binning (opt-in external tools)
#
# The wrapper calls below require external tools and inputs.

if RUN_EXTERNAL
    println("\n=== External binning runs (opt-in) ===")

    if isfile(contigs_fasta) && isfile(depth_file)
        println("Running VAMB, MetaBAT2, TaxVAMB, Taxometer...")
        outdir_vamb = mktempdir()
        outdir_metabat = mktempdir()
        outdir_taxvamb = mktempdir()
        outdir_taxometer = mktempdir()
        try
            Mycelia.run_vamb(contigs_fasta=contigs_fasta, depth_file=depth_file, outdir=outdir_vamb)
            Mycelia.run_metabat2(contigs_fasta=contigs_fasta, depth_file=depth_file, outdir=outdir_metabat)
            Mycelia.run_taxvamb(contigs_fasta=contigs_fasta, depth_file=depth_file, outdir=outdir_taxvamb)
            Mycelia.run_taxometer(contigs_fasta=contigs_fasta, depth_file=depth_file, outdir=outdir_taxometer)
        finally
            rm(outdir_vamb; recursive=true, force=true)
            rm(outdir_metabat; recursive=true, force=true)
            rm(outdir_taxvamb; recursive=true, force=true)
            rm(outdir_taxometer; recursive=true, force=true)
        end
    else
        println("Skipping VAMB/MetaBAT2/TaxVAMB/Taxometer (missing contigs/depth inputs).")
    end

    if isfile(contigs_fasta) && isfile(assembly_graph) && isfile(mapping_file)
        println("Running MetaCoAG...")
        outdir_metacoag = mktempdir()
        try
            Mycelia.run_metacoag(
                contigs_fasta=contigs_fasta,
                assembly_graph=assembly_graph,
                mapping_file=mapping_file,
                outdir=outdir_metacoag
            )
        finally
            rm(outdir_metacoag; recursive=true, force=true)
        end
    else
        println("Skipping MetaCoAG (missing contigs/graph/mapping inputs).")
    end

    if isfile(contigs_fasta) && isfile(coverage_table)
        println("Running GenomeFace...")
        outdir_genomeface = mktempdir()
        try
            Mycelia.run_genomeface(
                contigs_fasta=contigs_fasta,
                coverage_table=coverage_table,
                outdir=outdir_genomeface
            )
        finally
            rm(outdir_genomeface; recursive=true, force=true)
        end
    else
        println("Skipping GenomeFace (missing contigs/coverage inputs).")
    end

    if isfile(contigs_fasta) && isfile(coverage_table) && isfile(marker_file)
        println("Running COMEBin...")
        outdir_comebin = mktempdir()
        try
            Mycelia.run_comebin(
                contigs_fasta=contigs_fasta,
                coverage_table=coverage_table,
                marker_file=marker_file,
                outdir=outdir_comebin
            )
        finally
            rm(outdir_comebin; recursive=true, force=true)
        end
    else
        println("Skipping COMEBin (missing contigs/coverage/marker inputs).")
    end
else
    println("\nExternal tools are opt-in. Set MYCELIA_RUN_EXTERNAL=true to run binning wrappers.")
end

# ## Part 3: Post-binning (opt-in external tools)
#
# dRep and MAGmax require external tools and user-provided inputs.

if RUN_EXTERNAL
    if !isempty(genomes) && all(isfile, genomes)
        println("Running dRep...")
        outdir_drep = mktempdir()
        try
            Mycelia.run_drep_dereplicate(genomes=genomes, outdir=outdir_drep)
        finally
            rm(outdir_drep; recursive=true, force=true)
        end
    else
        println("Skipping dRep (missing MYCELIA_BINNING_GENOMES inputs).")
    end

    if !isempty(bins_dirs) && all(isdir, bins_dirs)
        println("Running MAGmax...")
        outdir_magmax = mktempdir()
        try
            Mycelia.run_magmax_merge(bins_dirs=bins_dirs, outdir=outdir_magmax)
        finally
            rm(outdir_magmax; recursive=true, force=true)
        end
    else
        println("Skipping MAGmax (missing MYCELIA_BINNING_BINS_DIRS inputs).")
    end
end
