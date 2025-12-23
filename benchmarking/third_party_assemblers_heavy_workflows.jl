# Heavy third-party assembler workflows (single-node benchmarking).
#
# This script is intentionally NOT part of `benchmarking/run_all_benchmarks.jl`.
# It runs workflows that can be slow and/or memory-intensive even on small inputs:
# - Canu (long-read isolate)
# - hifiasm-meta (long-read metagenomic)
# - metaMDBG (long-read metagenomic, HiFi + ONT)
# - STRONG + Strainy (strain-aware post-assembly workflows)
#
# Usage (examples):
#   julia --project=. benchmarking/third_party_assemblers_heavy_workflows.jl
#   MYCELIA_BENCHMARK_THREADS=32 julia --project=. benchmarking/third_party_assemblers_heavy_workflows.jl

import Pkg
Pkg.activate(dirname(@__DIR__))

import Dates
import Mycelia
import StableRNGs
import BioSequences
import FASTX

threads = clamp(
    something(tryparse(Int, get(ENV, "MYCELIA_BENCHMARK_THREADS", string(Mycelia.get_default_threads()))), Mycelia.get_default_threads()),
    1,
    Sys.CPU_THREADS > 0 ? Sys.CPU_THREADS : typemax(Int)
)

println("=== Third-party heavy workflows ===")
println("Start: $(Dates.now())")
println("Threads: $(threads)")

format_genome_size(bp::Integer) = bp < 1_000_000 ? string(cld(bp, 1_000), "k") : string(cld(bp, 1_000_000), "m")

mktempdir() do dir
    println("\nWorkdir: $dir")

    # ---------------------------------------------------------------------
    # Canu (long-read isolate)
    # ---------------------------------------------------------------------
    let
        genome_len = something(tryparse(Int, get(ENV, "MYCELIA_BENCHMARK_CANU_GENOME_LEN", "25000")), 25_000)
        coverage = get(ENV, "MYCELIA_BENCHMARK_CANU_COVERAGE", "15x")

        ref_fasta = joinpath(dir, "canu_reference.fasta")
        rng = StableRNGs.StableRNG(42)
        genome = BioSequences.randdnaseq(rng, genome_len)
        Mycelia.write_fasta(outfile=ref_fasta, records=[FASTX.FASTA.Record("canu_test_genome", genome)])

        reads_gz = Mycelia.simulate_pacbio_reads(fasta=ref_fasta, quantity=coverage, quiet=true)
        reads_fastq = joinpath(dir, "canu_reads.fq")
        run(pipeline(`gunzip -c $(reads_gz)`, reads_fastq))

        outdir = joinpath(dir, "canu_assembly")
        t0 = time()
        result = Mycelia.run_canu(
            fastq=reads_fastq,
            outdir=outdir,
            genome_size=format_genome_size(genome_len),
            stopOnLowCoverage=8,
            threads=threads,
            cor_threads=min(threads, 8)
        )
        elapsed = round(time() - t0, digits=1)
        println("\nCanu finished in $(elapsed)s")
        println("  Assembly: $(result.assembly)")
    end

    # ---------------------------------------------------------------------
    # hifiasm-meta + metaMDBG (long-read metagenomic)
    # ---------------------------------------------------------------------
    let
        genome_lens = (
            something(tryparse(Int, get(ENV, "MYCELIA_BENCHMARK_META_GENOME1_LEN", "35000")), 35_000),
            something(tryparse(Int, get(ENV, "MYCELIA_BENCHMARK_META_GENOME2_LEN", "30000")), 30_000),
            something(tryparse(Int, get(ENV, "MYCELIA_BENCHMARK_META_GENOME3_LEN", "25000")), 25_000),
        )
        coverage = get(ENV, "MYCELIA_BENCHMARK_META_COVERAGE", "20x")

        ref_fasta = joinpath(dir, "metagenomic_long_ref.fasta")
        rngs = (StableRNGs.StableRNG(789), StableRNGs.StableRNG(790), StableRNGs.StableRNG(791))
        genomes = (
            BioSequences.randdnaseq(rngs[1], genome_lens[1]),
            BioSequences.randdnaseq(rngs[2], genome_lens[2]),
            BioSequences.randdnaseq(rngs[3], genome_lens[3]),
        )
        Mycelia.write_fasta(
            outfile=ref_fasta,
            records=[
                FASTX.FASTA.Record("meta_genome_1", genomes[1]),
                FASTX.FASTA.Record("meta_genome_2", genomes[2]),
                FASTX.FASTA.Record("meta_genome_3", genomes[3]),
            ]
        )

        hifi_reads_gz = Mycelia.simulate_pacbio_reads(fasta=ref_fasta, quantity=coverage, quiet=true)
        hifi_fastq = joinpath(dir, "meta_hifi_reads.fq")
        run(pipeline(`gunzip -c $(hifi_reads_gz)`, hifi_fastq))

        let
            outdir = joinpath(dir, "hifiasm_meta_assembly")
            t0 = time()
            result = Mycelia.run_hifiasm_meta(
                fastq=hifi_fastq,
                outdir=outdir,
                bloom_filter=0,
                read_selection=true,
                threads=threads
            )
            elapsed = round(time() - t0, digits=1)
            println("\nhifiasm-meta finished in $(elapsed)s")
            println("  Outprefix: $(result.hifiasm_outprefix)")
        end

        let
            outdir = joinpath(dir, "metamdbg_hifi_assembly")
            t0 = time()
            result = Mycelia.run_metamdbg(hifi_reads=hifi_fastq, outdir=outdir, abundance_min=2, threads=threads)
            elapsed = round(time() - t0, digits=1)
            println("\nmetaMDBG (HiFi) finished in $(elapsed)s")
            println("  Contigs: $(result.contigs)")
            println("  Graph: $(result.graph)")
        end

        let
            ont_cov = get(ENV, "MYCELIA_BENCHMARK_META_ONT_COVERAGE", "10x")
            ont_reads_gz = Mycelia.simulate_nanopore_reads(fasta=ref_fasta, quantity=ont_cov, quiet=true)
            ont_fastq = joinpath(dir, "meta_ont_reads.fq")
            run(pipeline(`gunzip -c $(ont_reads_gz)`, ont_fastq))

            outdir = joinpath(dir, "metamdbg_ont_assembly")
            t0 = time()
            result = Mycelia.run_metamdbg(ont_reads=ont_fastq, outdir=outdir, abundance_min=2, threads=threads)
            elapsed = round(time() - t0, digits=1)
            println("\nmetaMDBG (ONT) finished in $(elapsed)s")
            println("  Contigs: $(result.contigs)")
            println("  Graph: $(result.graph)")
        end
    end

    # ---------------------------------------------------------------------
    # STRONG + Strainy (strain-aware workflows)
    # ---------------------------------------------------------------------
    let
        base_ref_fasta = joinpath(dir, "base_strain.fasta")
        variant_ref_fasta = joinpath(dir, "variant_strain.fasta")
        rng_base = StableRNGs.StableRNG(901)
        rng_variant = StableRNGs.StableRNG(902)
        base_genome = BioSequences.randdnaseq(rng_base, 5_000)
        variant_genome = BioSequences.randdnaseq(rng_variant, 5_000)
        Mycelia.write_fasta(outfile=base_ref_fasta, records=[FASTX.FASTA.Record("base_strain", base_genome)])
        Mycelia.write_fasta(outfile=variant_ref_fasta, records=[FASTX.FASTA.Record("variant_strain", variant_genome)])

        base_reads_gz = Mycelia.simulate_nanopore_reads(fasta=base_ref_fasta, quantity="12x", quiet=true)
        variant_reads_gz = Mycelia.simulate_nanopore_reads(fasta=variant_ref_fasta, quantity="4x", quiet=true)

        base_fastq = joinpath(dir, "base_strain_reads.fq")
        variant_fastq = joinpath(dir, "variant_strain_reads.fq")
        mixed_fastq = joinpath(dir, "mixed_strain_reads.fq")
        run(pipeline(`gunzip -c $(base_reads_gz)`, base_fastq))
        run(pipeline(`gunzip -c $(variant_reads_gz)`, variant_fastq))
        run(pipeline(`cat $(base_fastq) $(variant_fastq)`, mixed_fastq))

        # metaFlye graph generation
        metaflye_outdir = joinpath(dir, "metaflye_for_strain_workflows")
        t0 = time()
        Mycelia.run_metaflye(
            fastq=mixed_fastq,
            outdir=metaflye_outdir,
            genome_size="5k",
            read_type="nano-raw",
            threads=threads
        )
        elapsed = round(time() - t0, digits=1)
        println("\nmetaFlye (strain-workflows) finished in $(elapsed)s")

        assembly_graph_gfa = joinpath(metaflye_outdir, "assembly_graph.gfa")
        assembly_fasta = joinpath(metaflye_outdir, "assembly.fasta")

        if isfile(assembly_graph_gfa)
            strong_outdir = joinpath(dir, "strong_assembly")
            t0 = time()
            result = Mycelia.run_strong(assembly_graph_gfa, mixed_fastq, outdir=strong_outdir, nb_strains=2)
            elapsed = round(time() - t0, digits=1)
            println("\nSTRONG finished in $(elapsed)s")
            println("  Strain unitigs: $(result.strain_unitigs)")
        else
            println("\nSTRONG skipped (no assembly graph at $(assembly_graph_gfa))")
        end

        if isfile(assembly_fasta)
            strainy_outdir = joinpath(dir, "strainy_assembly")
            t0 = time()
            result = Mycelia.run_strainy(assembly_fasta, mixed_fastq, outdir=strainy_outdir, mode="phase", threads=threads)
            elapsed = round(time() - t0, digits=1)
            println("\nStrainy finished in $(elapsed)s")
            println("  Strain assemblies: $(result.strain_assemblies)")
        else
            println("\nStrainy skipped (no assembly at $(assembly_fasta))")
        end
    end
end

println("\nDone: $(Dates.now())")

