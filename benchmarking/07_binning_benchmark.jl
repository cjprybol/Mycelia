# # Binning Benchmark
#
# This benchmark exercises binning and post-binning wrappers with user-supplied inputs.
# External tool execution is opt-in via `MYCELIA_RUN_EXTERNAL=true`.

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import Dates

println("=== Binning Benchmark ===")
println("Start time: $(Dates.now())")

const RUN_EXTERNAL = get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
const contigs_fasta = get(ENV, "MYCELIA_BINNING_CONTIGS", "")
const depth_file = get(ENV, "MYCELIA_BINNING_DEPTH", "")
const taxonomy_file = get(ENV, "MYCELIA_BINNING_TAXONOMY", "")
const bam_path = get(ENV, "MYCELIA_BINNING_BAM_PATH", "")
const coverage_table = get(ENV, "MYCELIA_BINNING_COVERAGE_TABLE", "")
const assembly_graph = get(ENV, "MYCELIA_BINNING_GRAPH", "")
const mapping_file = get(ENV, "MYCELIA_BINNING_MAPPING", "")
const genomes = filter(!isempty, strip.(split(get(ENV, "MYCELIA_BINNING_GENOMES", ""), ',')))
const bins_dirs = filter(!isempty, strip.(split(get(ENV, "MYCELIA_BINNING_BINS_DIRS", ""), ',')))

if RUN_EXTERNAL
    println("External binning benchmarks enabled.")

    if isfile(contigs_fasta) && isfile(depth_file)
        println("\n--- VAMB / MetaBAT2 ---")
        outdir_vamb = mktempdir()
        outdir_metabat = mktempdir()
        try
            start = time()
            Mycelia.run_vamb(contigs_fasta=contigs_fasta, depth_file=depth_file, outdir=outdir_vamb)
            println("VAMB elapsed: $(round(time() - start, digits=2))s")

            start = time()
            Mycelia.run_metabat2(contigs_fasta=contigs_fasta, depth_file=depth_file, outdir=outdir_metabat)
            println("MetaBAT2 elapsed: $(round(time() - start, digits=2))s")
        finally
            rm(outdir_vamb; recursive=true, force=true)
            rm(outdir_metabat; recursive=true, force=true)
        end
    else
        println("Skipping VAMB/MetaBAT2 (missing contigs/depth inputs).")
    end

    if isfile(contigs_fasta) && isfile(depth_file) && isfile(taxonomy_file)
        println("\n--- TaxVAMB / Taxometer ---")
        outdir_taxvamb = mktempdir()
        outdir_taxometer = mktempdir()
        try
            start = time()
            Mycelia.run_taxvamb(
                contigs_fasta=contigs_fasta,
                depth_file=depth_file,
                taxonomy_file=taxonomy_file,
                outdir=outdir_taxvamb
            )
            println("TaxVAMB elapsed: $(round(time() - start, digits=2))s")

            start = time()
            Mycelia.run_taxometer(
                contigs_fasta=contigs_fasta,
                depth_file=depth_file,
                taxonomy_file=taxonomy_file,
                outdir=outdir_taxometer
            )
            println("Taxometer elapsed: $(round(time() - start, digits=2))s")
        finally
            rm(outdir_taxvamb; recursive=true, force=true)
            rm(outdir_taxometer; recursive=true, force=true)
        end
    else
        println("Skipping TaxVAMB/Taxometer (missing contigs/depth/taxonomy inputs).")
    end

    if isfile(contigs_fasta) && isfile(assembly_graph) && isfile(mapping_file)
        println("\n--- MetaCoAG ---")
        outdir_metacoag = mktempdir()
        try
            start = time()
            Mycelia.run_metacoag(
                contigs_fasta=contigs_fasta,
                assembly_graph=assembly_graph,
                mapping_file=mapping_file,
                outdir=outdir_metacoag
            )
            println("MetaCoAG elapsed: $(round(time() - start, digits=2))s")
        finally
            rm(outdir_metacoag; recursive=true, force=true)
        end
    else
        println("Skipping MetaCoAG (missing contigs/graph/mapping inputs).")
    end

    if isfile(contigs_fasta) && isfile(coverage_table)
        println("\n--- GenomeFace ---")
        println("Skipping GenomeFace; wrapper disabled until genomeface executable is available again.")
    else
        println("Skipping GenomeFace (missing contigs/coverage inputs; wrapper disabled until executable is available).")
    end

    if isfile(contigs_fasta) && (!isempty(bam_path)) && (isfile(bam_path) || isdir(bam_path))
        println("\n--- COMEBin ---")
        outdir_comebin = mktempdir()
        try
            start = time()
            Mycelia.run_comebin(
                contigs_fasta=contigs_fasta,
                bam_path=bam_path,
                outdir=outdir_comebin
            )
            println("COMEBin elapsed: $(round(time() - start, digits=2))s")
        finally
            rm(outdir_comebin; recursive=true, force=true)
        end
    else
        println("Skipping COMEBin (missing contigs/BAM inputs).")
    end

    if !isempty(genomes) && all(isfile, genomes)
        println("\n--- dRep ---")
        outdir_drep = mktempdir()
        try
            start = time()
            Mycelia.run_drep_dereplicate(genomes=genomes, outdir=outdir_drep)
            println("dRep elapsed: $(round(time() - start, digits=2))s")
        finally
            rm(outdir_drep; recursive=true, force=true)
        end
    else
        println("Skipping dRep (missing MYCELIA_BINNING_GENOMES inputs).")
    end

    if !isempty(bins_dirs) && all(isdir, bins_dirs)
        println("\n--- MAGmax ---")
        outdir_magmax = mktempdir()
        try
            start = time()
            Mycelia.run_magmax_merge(bins_dirs=bins_dirs, outdir=outdir_magmax)
            println("MAGmax elapsed: $(round(time() - start, digits=2))s")
        finally
            rm(outdir_magmax; recursive=true, force=true)
        end
    else
        println("Skipping MAGmax (missing MYCELIA_BINNING_BINS_DIRS inputs).")
    end
else
    println("External benchmarks are opt-in; set MYCELIA_RUN_EXTERNAL=true to enable.")
end

println("End time: $(Dates.now())")
