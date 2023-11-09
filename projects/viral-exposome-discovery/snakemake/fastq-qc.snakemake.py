import os
    
sample_directory = config["directory"]
sample_id = config["sample_id"]
forward_reads = config["forward"]
reverse_reads = config["reverse"]
fastqc_outdir = sample_directory + '/fastqc'

forward_base = os.path.basename(forward_reads).replace(".fastq.gz", "")
reverse_base = os.path.basename(reverse_reads).replace(".fastq.gz", "")

forward_html_report = fastqc_outdir + "/" + forward_base + "_fastqc.html"
reverse_html_report = fastqc_outdir + "/" + reverse_base + "_fastqc.html"
forward_zip = fastqc_outdir + "/" + forward_base + "_fastqc.zip"
reverse_zip = fastqc_outdir + "/" + reverse_base + "_fastqc.zip"

# snakemake --use-conda  --snakefile snakemake/fastq-qc.snakemake.py --cores 1 fastqc --config directory=data/SRA/SRR6399459 sample_id=SRR6399459 forward=data/SRA/SRR6399459/SRR6399459_1.fastq.gz reverse=data/SRA/SRR6399459/SRR6399459_2.fastq.gz
rule fastqc:
    conda:
        "../environment.yml"
    input: 
        forward_reads,
        reverse_reads
    output: 
        forward_html_report,
        reverse_html_report,
        forward_zip,
        reverse_zip
    shell:
        """
        mkdir -p {fastqc_outdir}
        fastqc --outdir {fastqc_outdir} {forward_reads} {reverse_reads}
        """

trim_galore_outdir = sample_directory + '/trim_galore'
trimmed_forward_reads = trim_galore_outdir + f'/{forward_base}_val_1.fq.gz'
trimmed_reverse_reads = trim_galore_outdir + f'/{reverse_base}_val_2.fq.gz'
forward_trimming_report = trim_galore_outdir + f'/{forward_base}.fastq.gz_trimming_report.txt'
reverse_trimming_report = trim_galore_outdir + f'/{reverse_base}.fastq.gz_trimming_report.txt'

# snakemake --use-conda  --snakefile snakemake/fastq-qc.snakemake.py --cores 1 trim_galore --config directory=data/SRA/SRR6399459 sample_id=SRR6399459 forward=data/SRA/SRR6399459/SRR6399459_1.fastq.gz reverse=data/SRA/SRR6399459/SRR6399459_2.fastq.gz
rule trim_galore:
    conda:
        "../environment.yml"
    input:
        forward_reads,
        reverse_reads
    output: 
        trimmed_forward_reads,
        trimmed_reverse_reads,
        forward_trimming_report,
        reverse_trimming_report
    shell:
        """
        trim_galore --output_dir {trim_galore_outdir} --paired {forward_reads} {reverse_reads}
        """

post_trim_galore_fastqc_outdir = sample_directory + '/post_trim_galore_fastqc'
trimmed_forward_base = os.path.basename(trimmed_forward_reads).replace(".fq.gz", "")
trimmed_reverse_base = os.path.basename(trimmed_reverse_reads).replace(".fq.gz", "")
post_trim_forward_html_report = post_trim_galore_fastqc_outdir + "/" + trimmed_forward_base + "_fastqc.html"
post_trim_reverse_html_report = post_trim_galore_fastqc_outdir + "/" + trimmed_reverse_base + "_fastqc.html"
post_trim_forward_zip = post_trim_galore_fastqc_outdir + "/" + trimmed_forward_base + "_fastqc.zip"
post_trim_reverse_zip = post_trim_galore_fastqc_outdir + "/" + trimmed_reverse_base + "_fastqc.zip"

# snakemake --use-conda  --snakefile snakemake/fastq-qc.snakemake.py --cores 1 post_trim_galore_fastqc --config directory=data/SRA/SRR6399459 sample_id=SRR6399459 forward=data/SRA/SRR6399459/SRR6399459_1.fastq.gz reverse=data/SRA/SRR6399459/SRR6399459_2.fastq.gz
rule post_trim_galore_fastqc:
    conda:
        "../environment.yml"
    input:
        trimmed_forward_reads,
        trimmed_reverse_reads
    output:
        post_trim_forward_html_report,
        post_trim_reverse_html_report,
        post_trim_forward_zip,
        post_trim_reverse_zip
    shell:
        """
        mkdir -p {post_trim_galore_fastqc_outdir}
        fastqc --outdir {post_trim_galore_fastqc_outdir} {trimmed_forward_reads} {trimmed_reverse_reads}
        """

# snakemake --use-conda  --snakefile snakemake/fastq-qc.snakemake.py --cores 1 all --config directory=data/SRA/SRR6399459 sample_id=SRR6399459 forward=data/SRA/SRR6399459/SRR6399459_1.fastq.gz reverse=data/SRA/SRR6399459/SRR6399459_2.fastq.gz
rule all:
    conda:
        "../environment.yml"
    input:
        forward_html_report,
        reverse_html_report,
        forward_zip,
        reverse_zip,
        trimmed_forward_reads,
        trimmed_reverse_reads,
        forward_trimming_report,
        reverse_trimming_report,
        post_trim_forward_html_report,
        post_trim_reverse_html_report,
        post_trim_forward_zip,
        post_trim_reverse_zip
    shell:
        """
        echo "all done!"
        """