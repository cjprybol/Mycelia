# # Tutorial 13b: Rhizomorph Corrector (Read Correction with `:scalable`)
#
# Tutorial 13 assembled reads directly. This companion tutorial adds the
# **Rhizomorph read corrector** — an opt-in front-end that decodes each read
# against a graph-structured hidden Markov model (maximum-likelihood correction)
# before assembly. On short-read data it collapses the fragmented, redundant
# contigs of an uncorrected assembly into a single near-full-length contig.
#
# We run the *same reads* through two arms and compare:
#
# - **baseline** — `assemble_genome(reads; k=21, corrector=:none)`
# - **corrected** — `assemble_genome(reads; k=21, corrector=:iterative, strategy=:scalable)`
#
# See the validated real-genome numbers on the
# [Benchmarks](../../benchmarks.md) page and the algorithm write-up in
# `docs/design/2026-07-graph-as-hmm-corrector-methods.md`.

# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/13_rhizomorph_corrector.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import Random

Random.seed!(42)

println("=== Rhizomorph Corrector Tutorial ===")

# ## Step 1: Download a small reference genome
#
# We use phiX174 (5,386 bp) so the tutorial runs in well under a minute.

println("Downloading reference genome (phiX174)...")
reference_gz = Mycelia.download_genome_by_accession(accession = "NC_001422.1")

## The download is a .fna.gz. Decompress to a clean path inside a work dir so
## downstream filenames (ART outputs) never embed ".gz" in the middle.
work_dir = mktempdir()
reference_file = joinpath(work_dir, "phix174.fna")
run(pipeline(`gzip -dc $(reference_gz)`, stdout = reference_file))

ref_sequence = FASTX.sequence(String, first(Mycelia.open_fastx(reference_file)))
ref_length = length(ref_sequence)
println("Reference length: $ref_length bp")

# ## Step 2: Simulate short reads
#
# We simulate paired-end Illumina reads with ART's HS25 (HiSeq 2500) empirical
# error model — the same read model used for the committed real-genome
# validation. Short reads are what make read correction matter: an uncorrected
# short-read assembly of this genome shatters into hundreds of contigs.
#
# For a real project, replace this block with your own FASTQ inputs.

println("Simulating Illumina reads (ART HS25, 50x, 150 bp PE)...")
art = Mycelia.simulate_illumina_reads(
    fasta = reference_file,
    coverage = 50,
    read_length = 150,
    seqSys = "HS25",
    paired = true,
    rndSeed = 42
)

## ART writes gzipped FASTQ; decompress and load both mates.
read_paths = String[]
for gz in (art.forward_reads, art.reverse_reads)
    fq = replace(gz, r"\.gz$" => "")
    isfile(fq) || run(`gunzip -k -f $(gz)`)
    push!(read_paths, fq)
end
reads = reduce(vcat, [collect(Mycelia.open_fastx(p)) for p in read_paths])
total_bases = sum(length(FASTX.sequence(String, r)) for r in reads)
println("Reads loaded: $(length(reads))  (~$(round(total_bases / ref_length, digits=1))x)")

# ## Step 3: Baseline assembly (no correction)
#
# `corrector=:none` is the default; we name it here to make the contrast
# explicit. This is the byte-identical single-k pipeline from Tutorial 13.

println("\n[arm 1/2] Baseline assembly (corrector=:none)...")
t0 = time()
baseline = Mycelia.Rhizomorph.assemble_genome(
    reads;
    k = 21,
    corrector = :none
)
baseline_elapsed = round(time() - t0, digits = 1)

# ## Step 4: Corrected assembly (`corrector=:iterative, strategy=:scalable`)
#
# `corrector=:iterative` diverts to the read corrector, then re-assembles the
# corrected reads. `strategy=:scalable` (the default when a corrector is
# requested) is the real-scale tier: a coarse 3-rung k-ladder with a low
# per-k iteration cap, a Stage-0 linear k-mer-spectrum pre-correction, skip-solid
# volume reduction, and a size-aware Viterbi beam. The alternative,
# `strategy=:exhaustive`, does an exact prime-by-prime k-walk with an unbounded
# beam — higher sensitivity but only appropriate for small inputs (it can run out
# of memory on large read sets).
#
# Note: the corrector's emission model needs per-base Phred quality, so
# quality-less input (FASTA) is assigned a placeholder Q40 with a warning.
# `k` is treated as the ceiling of the ladder and floored internally at 13.

println("[arm 2/2] Corrected assembly (corrector=:iterative, strategy=:scalable)...")
t0 = time()
corrected = Mycelia.Rhizomorph.assemble_genome(
    reads;
    k = 21,
    corrector = :iterative,
    strategy = :scalable
)
corrected_elapsed = round(time() - t0, digits = 1)

# ## Step 5: Compare the two arms
#
# The corrector's provenance is stamped onto `assembly_stats`.

function contig_summary(assembly)
    contigs = assembly.contigs
    n = length(contigs)
    largest = isempty(contigs) ? 0 : maximum(length, contigs)
    total = isempty(contigs) ? 0 : sum(length, contigs)
    return (; n, largest, total)
end

b = contig_summary(baseline)
c = contig_summary(corrected)

println("\n" * "="^64)
println("RESULT — same reads, two arms (phiX174 reference = $ref_length bp)")
println("="^64)
println("  baseline  (:none)      contigs=$(b.n)  largest=$(b.largest) bp  total=$(b.total) bp  t=$(baseline_elapsed)s")
println("  corrected (:scalable)  contigs=$(c.n)  largest=$(c.largest) bp  total=$(c.total) bp  t=$(corrected_elapsed)s")

# The corrector reports what it did:
println("\nCorrector provenance (assembly_stats):")
for key in ("corrector", "strategy", "skip_solid_effective", "k_progression",
    "corrected_read_count", "reassembly_graph_reused")
    if haskey(corrected.assembly_stats, key)
        println("  $key = $(corrected.assembly_stats[key])")
    end
end

# ## Interpretation
#
# On short reads the baseline assembly fragments into many redundant contigs
# (each strand assembled separately, plus error debris), while the corrected
# assembly resolves to a single near-full-length contig. The trade-off is
# runtime: correction adds a per-read maximum-likelihood decode, so the corrected
# arm is slower. On the committed real-genome benchmark (phiX174 + lambda,
# ART HS25 50x), the `:scalable` corrector collapses thousands of naive contigs
# (N50 = 41) into one contig at the full genome length while holding ~99.3–100%
# genome fraction and ~99.99–100% identity. See the
# [Benchmarks](../../benchmarks.md) page for the full table.

# ## Next steps
#
# - Swap the simulated reads for your own short-read FASTQ.
# - Try `strategy=:exhaustive` on a very small input for maximum sensitivity.
# - Tune `k` (the ladder ceiling) and inspect how `k_progression` changes.
# - Inspect `corrected.graph` and export GFA via
#   `Mycelia.Rhizomorph.write_gfa` for downstream graph analysis.
