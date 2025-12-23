# # Tutorial 13: Rhizomorph Assembly (Assembly in 5 Minutes)
#
# This tutorial is the updated "assembly in 5 minutes" workflow using Mycelia's
# Rhizomorph assembly pipeline. It walks from data acquisition to assembled
# contigs with a small viral genome example.

# ## Setup
import Pkg
if isinteractive()
    Pkg.activate("..")
end

import Mycelia
import FASTX
import Random
import Statistics
import Kmers

Random.seed!(42)

println("=== Rhizomorph Assembly Tutorial ===")

# ## Step 1: Installation (one-time setup)
#
# Uncomment if Mycelia is not installed yet.
#
# import Pkg
# Pkg.add(url="https://github.com/cjprybol/Mycelia.git")

# ## Step 2: Download a small reference genome
#
# We'll assemble phiX174 (5,386 bp) so the tutorial runs quickly.

println("Downloading reference genome (phiX174)...")
reference_file = Mycelia.download_genome_by_accession(accession="NC_001422.1")
println("Reference FASTA: $reference_file")

ref_record = first(Mycelia.open_fastx(reference_file))
ref_length = length(FASTX.sequence(ref_record))
println("Reference length: $ref_length bp")

# ## Step 3: Simulate sequencing reads
#
# For real projects, replace this with your own FASTQ.

println("Simulating PacBio-style reads...")
reads_file = Mycelia.simulate_pacbio_reads(
    fasta=reference_file,
    quantity="20x"
)
println("Reads FASTQ: $reads_file")

reads = collect(Mycelia.open_fastx(reads_file))
mean_read_length = round(Statistics.mean(length(FASTX.sequence(r)) for r in reads))
total_bases = sum(length(FASTX.sequence(r)) for r in reads)
println("Reads loaded: $(length(reads))")
println("Mean read length: $mean_read_length bp")
println("Total bases: $total_bases bp")

# ## Step 4: Run Rhizomorph assembly
#
# The unified Rhizomorph assembly API auto-detects read type and builds the
# appropriate graph (qualmer graph for FASTQ).

println("Running Rhizomorph assembly...")
start_time = time()
assembly = Mycelia.Rhizomorph.assemble_genome(
    reads;
    k=31,
    error_rate=0.01,
    min_coverage=3
)
elapsed = round(time() - start_time, digits=1)
println("Assembly completed in $elapsed seconds")

# ## Step 5: Inspect results

contigs = assembly.contigs
contig_names = assembly.contig_names
println("Contigs assembled: $(length(contigs))")

if !isempty(contigs)
    contig_lengths = [length(contig) for contig in contigs]
    println("Total contig length: $(sum(contig_lengths)) bp")
    println("Longest contig: $(maximum(contig_lengths)) bp")
else
    println("No contigs were produced; try adjusting parameters.")
end

if !isempty(contigs)
    # Validate against the reference (placeholder metrics for now).
    metrics = Mycelia.Rhizomorph.validate_assembly(assembly; reference=FASTX.sequence(ref_record))
    println("Validation metrics:")
    for (key, value) in sort(collect(metrics))
        println("  $key: $value")
    end
end

# ## Step 6: Save the assembly
#
# Write contigs to FASTA and, when available, FASTQ with quality scores.

if !isempty(contigs)
    fasta_records = FASTX.FASTA.Record[]
    for (name, seq) in zip(contig_names, contigs)
        push!(fasta_records, FASTX.FASTA.Record(name, seq))
    end
    output_fasta = Mycelia.write_fasta(records=fasta_records, outfile="phiX174_rhizomorph_assembly.fasta")
    println("Saved contigs: $output_fasta")

    if Mycelia.Rhizomorph.has_quality_information(assembly)
        output_fastq = "phiX174_rhizomorph_assembly.fastq"
        Mycelia.Rhizomorph.write_fastq_contigs(assembly, output_fastq)
        println("Saved quality-aware contigs: $output_fastq")
    else
        println("No quality-aware contigs available; FASTQ export skipped.")
    end
end

# ## Step 7: Optional k-mer recovery check
#
# Compare 21-mer content between reference and assembly.

if !isempty(contigs)
    output_fasta = "phiX174_rhizomorph_assembly.fasta"
    if isfile(output_fasta)
        ref_kmers = Mycelia.count_canonical_kmers(Kmers.DNAKmer{21}, reference_file)
        asm_kmers = Mycelia.count_canonical_kmers(Kmers.DNAKmer{21}, output_fasta)
        shared_kmers = length(intersect(keys(ref_kmers), keys(asm_kmers)))
        kmer_recovery = round(100 * shared_kmers / length(ref_kmers), digits=1)
        println("K-mer recovery rate (21-mers): $kmer_recovery%")
    else
        println("K-mer recovery skipped; assembly FASTA not found.")
    end
end

# ## Next steps
#
# - Replace the simulated reads with your own FASTQ inputs.
# - Try alternate parameters (k, min_coverage, error_rate).
# - Explore Rhizomorph graph building functions for debugging and visualization.
