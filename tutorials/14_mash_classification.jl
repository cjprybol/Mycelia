# # Tutorial 14: Mash Classification and Screening
#
# This tutorial demonstrates how to use Mash sketches for quick containment
# and distance screening against reference sketches.
#
# ## Learning Objectives
#
# By the end of this tutorial, you will understand:
# - How to sketch reference genomes with Mash
# - How to sketch reads with Mash (`-r`, `-m`)
# - How to screen reads against a sketch database
# - How to interpret identity and shared-hash outputs

# ## Setup
# From the Mycelia base directory, convert this tutorial to a notebook:
#
# ```bash
# julia --project=. -e 'import Literate; Literate.notebook("tutorials/14_mash_classification.jl", "tutorials/notebooks", execute=false)'
# ```

import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import Random

Random.seed!(42)

println("=== Mash Classification Tutorial ===")

# ## Part 1: Sketching Reference Genomes

println("\n--- Reference Sketching ---")

ref_records = [
    Mycelia.random_fasta_record(moltype=:DNA, seed=1, L=800),
    Mycelia.random_fasta_record(moltype=:DNA, seed=2, L=800),
    Mycelia.random_fasta_record(moltype=:DNA, seed=3, L=800)
]

ref_files = String[]
for (i, record) in enumerate(ref_records)
    filename = "mash_reference_$(i).fasta"
    Mycelia.write_fasta(outfile=filename, records=[record])
    push!(ref_files, filename)
end

if get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
    reference_sketches = Mycelia.run_mash_sketch(
        input_files=ref_files,
        outdir=joinpath(pwd(), "mash_sketches"),
        k=21,
        s=1000
    ).sketches

    reference_db = Mycelia.run_mash_paste(
        out_file=joinpath(pwd(), "mash_sketches", "reference_db.msh"),
        in_files=reference_sketches
    )

    println("Reference sketch database: $(reference_db)")
else
    println("Skipping Mash reference sketching; set MYCELIA_RUN_EXTERNAL=true to run.")
end

# ## Part 2: Sketching Reads

println("\n--- Read Sketching ---")

reads = Mycelia.create_test_reads(String(FASTX.sequence(ref_records[1])), 25, 0.01)
reads_fastq = "mash_reads.fastq"
Mycelia.write_fastq(records=reads, filename=reads_fastq)

if get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
    read_sketch = Mycelia.run_mash_sketch(
        input_files=[reads_fastq],
        outdir=joinpath(pwd(), "mash_sketches"),
        k=21,
        s=1000,
        r=true,
        min_copies=2
    ).sketches[1]
    println("Read sketch: $(read_sketch)")
else
    println("Skipping Mash execution; set MYCELIA_RUN_EXTERNAL=true to run sketching.")
end

# ## Part 3: Screening Reads Against the Reference Database

println("\n--- Screening ---")

if get(ENV, "MYCELIA_RUN_EXTERNAL", "false") == "true"
    reference_db = joinpath(pwd(), "mash_sketches", "reference_db.msh")
    screen_result = Mycelia.run_mash_screen(
        reference=reference_db,
        query=reads_fastq,
        outdir=joinpath(pwd(), "mash_screen"),
        winner_takes_all=true
    )

    screen_table = Mycelia.parse_mash_screen_output(screen_result.results_tsv)
    println(screen_table)
else
    println("Skipping Mash screening; set MYCELIA_RUN_EXTERNAL=true to run screening.")
end

# ## Part 4: Interpretation Notes
#
# - **Identity** is the estimated ANI-like similarity (1.0 is identical).
# - **Shared hashes** are the number of shared MinHash sketches; more shared
#   hashes generally imply stronger containment.
# - Combine Mash screening with sourmash and sylph to triangulate the most
#   supported pangenome contexts before mapping with minimap2.
