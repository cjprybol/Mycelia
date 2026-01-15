# # Tutorial 12: Coverage Profiling with CoverM
#
# CoverM computes coverage and abundance summaries directly from BAMs.
# This tutorial shows how to:
# - Install/activate the CoverM Bioconda environment
# - Run contig-level coverage (`run_coverm_contig`)
# - Run genome/bin-level abundance (`run_coverm_genome`)
# - Reuse cached TSVs for reproducible workflows
#
# The workflow uses tiny synthetic data so it finishes quickly on a laptop or login node.

# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/12_coverm_coverage.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import CSV
import DataFrames
import StableRNGs
import Test

## Set up a scratch workspace
workdir = mktempdir()
println("Workspace: ", workdir)

## Create a small reference FASTA
rng = StableRNGs.StableRNG(42)
ref_record = Mycelia.random_fasta_record(L=5_000, seed=rand(rng, 0:typemax(Int)))
ref_fasta = joinpath(workdir, "reference.fa")
Mycelia.write_fasta(outfile=ref_fasta, records=[ref_record])

## Simulate short reads against the reference
reads = Mycelia.simulate_illumina_reads(
    fasta=ref_fasta,
    read_count=500,
    len=150,
    rndSeed=rand(rng, 0:typemax(Int)),
    quiet=true
)

## Map reads with minimap2 to produce a BAM
map_result = Mycelia.minimap_map(
    fasta=ref_fasta,
    fastq=reads.forward_reads,
    mapping_type="sr",
    threads=4
)
run(map_result.cmd)
bam_path = map_result.outfile

# ## Contig mode: per-contig coverage
contig_df = Mycelia.run_coverm_contig(
    bam_files=[bam_path],
    reference_fasta=ref_fasta,
    methods=["mean", "covered_fraction"],
    threads=4,
    outdir=joinpath(workdir, "coverm_contig")
)
println("Contig coverage (first rows):")
println(first(contig_df, 5))

# ## Genome mode: per-bin abundance
## Make a simple bin directory (one genome for this toy example)
bins_dir = joinpath(workdir, "bins")
mkpath(bins_dir)
bin_fasta = joinpath(bins_dir, "bin1.fa")
cp(ref_fasta, bin_fasta; force=true)

genome_df = Mycelia.run_coverm_genome(
    bam_files=[bam_path],
    genome_directory=bins_dir,
    genome_extension="fa",
    methods=["relative_abundance", "mean_coverage"],
    threads=4,
    outdir=joinpath(workdir, "coverm_genome")
)
println("Genome abundance (first rows):")
println(first(genome_df, 5))

# ## Reusing cached outputs
# Both wrappers skip execution when the output TSV already exists and is non-empty.
# You can point them at a shared path for deterministic reruns:
cached_tsv = joinpath(workdir, "coverm_genome", "coverm_genome.tsv")
genome_df_cached = Mycelia.run_coverm_genome(
    bam_files=[bam_path],
    genome_directory=bins_dir,
    output_tsv=cached_tsv
)
Test.@test genome_df == genome_df_cached

println("CoverM tutorial complete. Results in: ", workdir)
