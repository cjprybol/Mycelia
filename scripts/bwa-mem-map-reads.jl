# map reads to the assembly and run qualimap QC
bwt_index = "$(assembled_fasta).bwt"
if !isfile(bwt_index)
    run(`bwa index $(assembled_fasta)`)
end

mapped_reads_bam = "$(assembled_fasta).bwa.bam"
if !isfile(mapped_reads_bam)
    run(pipeline(
        `bwa mem -R "@RG\tID:$(config["sample identifier"])\tSM:bar" -t $(Sys.CPU_THREADS) $(assembled_fasta) $(TRIMMED_FORWARD) $(TRIMMED_REVERSE)`,
        `samtools collate -O - -`,
        `samtools fixmate -m - -`,
        `samtools sort`,
        `samtools markdup - -`,
        `samtools view -buh`,
        mapped_reads_bam))
end

if !isfile("$(mapped_reads_bam).bai")
    run(`samtools index $(mapped_reads_bam)`)
end

qualimap_report_pdf = "$(assembly_dir)/qualimap/report.pdf"
qualimap_report_txt = "$(assembly_dir)/qualimap/genome_results.txt"

if !isfile(qualimap_report_pdf) || !isfile(qualimap_report_txt)
    run(`
        qualimap bamqc
        -nt $(Sys.CPU_THREADS)
        -bam $(mapped_reads_bam)
        -outdir $(assembly_dir)/qualimap
        -outformat PDF:HTML
        --output-genome-coverage $(mapped_reads_bam).genome_coverage.txt
        `)
end

qualimap_contig_coverage_table = Mycelia.parse_qualimap_contig_coverage(qualimap_report_txt)
qualimap_contig_coverage_table[!, "% Mapped bases"] = round.(qualimap_contig_coverage_table[!, "Mapped bases"] ./ sum(qualimap_contig_coverage_table[!, "Mapped bases"]) .* 100, digits=3);