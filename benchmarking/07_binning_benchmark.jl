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

include("benchmark_utils.jl")

println("=== Binning Benchmark ===")
println("Start time: $(Dates.now())")

benchmark_suite = BenchmarkSuite("Binning Benchmark")

function record_profiled_benchmark!(suite::BenchmarkSuite, name::String, func; metadata=Dict{String, Any}())
    _, profile_stats = profile_execution(func)
    add_profiled_result!(suite, name, profile_stats; metadata=metadata)
    println("$(name) elapsed: $(round(profile_stats["wall_time_seconds"], digits=2))s, allocated $(round(profile_stats["allocated_mb"], digits=2)) MB")
end

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
        tmp_root = mktempdir()
        outdir_vamb = joinpath(tmp_root, "vamb_out")
        outdir_metabat = joinpath(tmp_root, "metabat_out")
        try
            record_profiled_benchmark!(
                benchmark_suite,
                "vamb",
                () -> Mycelia.run_vamb(contigs_fasta=contigs_fasta, depth_file=depth_file, outdir=outdir_vamb);
                metadata=Dict("outdir" => outdir_vamb)
            )

            record_profiled_benchmark!(
                benchmark_suite,
                "metabat2",
                () -> Mycelia.run_metabat2(contigs_fasta=contigs_fasta, depth_file=depth_file, outdir=outdir_metabat);
                metadata=Dict("outdir" => outdir_metabat)
            )
        finally
            rm(tmp_root; recursive=true, force=true)
        end
    else
        println("Skipping VAMB/MetaBAT2 (missing contigs/depth inputs).")
    end

    if isfile(contigs_fasta) && isfile(depth_file) && isfile(taxonomy_file)
        println("\n--- TaxVAMB / Taxometer ---")
        tmp_root = mktempdir()
        outdir_taxvamb = joinpath(tmp_root, "taxvamb_out")
        outdir_taxometer = joinpath(tmp_root, "taxometer_out")
        try
            record_profiled_benchmark!(
                benchmark_suite,
                "taxvamb",
                () -> Mycelia.run_taxvamb(
                    contigs_fasta=contigs_fasta,
                    depth_file=depth_file,
                    taxonomy_file=taxonomy_file,
                    outdir=outdir_taxvamb
                );
                metadata=Dict("outdir" => outdir_taxvamb)
            )

            record_profiled_benchmark!(
                benchmark_suite,
                "taxometer",
                () -> Mycelia.run_taxometer(
                    contigs_fasta=contigs_fasta,
                    depth_file=depth_file,
                    taxonomy_file=taxonomy_file,
                    outdir=outdir_taxometer
                );
                metadata=Dict("outdir" => outdir_taxometer)
            )
        finally
            rm(tmp_root; recursive=true, force=true)
        end
    else
        println("Skipping TaxVAMB/Taxometer (missing contigs/depth/taxonomy inputs).")
    end

    if isfile(contigs_fasta) && isfile(assembly_graph) && isfile(mapping_file)
        println("\n--- MetaCoAG ---")
        outdir_metacoag = mktempdir()
        try
            record_profiled_benchmark!(
                benchmark_suite,
                "metacoag",
                () -> Mycelia.run_metacoag(
                    contigs_fasta=contigs_fasta,
                    assembly_graph=assembly_graph,
                    mapping_file=mapping_file,
                    outdir=outdir_metacoag
                );
                metadata=Dict("outdir" => outdir_metacoag)
            )
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
            record_profiled_benchmark!(
                benchmark_suite,
                "comebin",
                () -> Mycelia.run_comebin(
                    contigs_fasta=contigs_fasta,
                    bam_path=bam_path,
                    outdir=outdir_comebin
                );
                metadata=Dict("outdir" => outdir_comebin)
            )
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
            record_profiled_benchmark!(
                benchmark_suite,
                "drep_dereplicate",
                () -> Mycelia.run_drep_dereplicate(genomes=genomes, outdir=outdir_drep);
                metadata=Dict("outdir" => outdir_drep, "n_genomes" => length(genomes))
            )
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
            record_profiled_benchmark!(
                benchmark_suite,
                "magmax_merge",
                () -> Mycelia.run_magmax_merge(bins_dirs=bins_dirs, outdir=outdir_magmax);
                metadata=Dict("outdir" => outdir_magmax, "n_bins_dirs" => length(bins_dirs))
            )
        finally
            rm(outdir_magmax; recursive=true, force=true)
        end
    else
        println("Skipping MAGmax (missing MYCELIA_BINNING_BINS_DIRS inputs).")
    end
else
    println("External benchmarks are opt-in; set MYCELIA_RUN_EXTERNAL=true to enable.")
end

results_dir = mkpath("results")
results_file = joinpath(results_dir, "binning_benchmark_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).json")
save_benchmark_results(benchmark_suite, results_file)
format_benchmark_summary(benchmark_suite)
println("End time: $(Dates.now())")
println("Results saved to: $(results_file)")
