# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

import os

# ################################################################################################
# # SUBSAMPLE READS FOR DEVELOPMENT TESTING
# ################################################################################################

# snakemake --use-conda  --snakefile snakemake/fastq-subset.snakemake.py --cores 1 subsample_paired_fastqs --config forward=data/SRA/SRR6399459/SRR6399459_1.fastq.gz reverse=data/SRA/SRR6399459/SRR6399459_2.fastq.gz sampling_rate=0.0001
# snakemake --use-conda  --snakefile snakemake/fastq-subset.snakemake.py --cores 1 subsample_paired_fastqs --config forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.fq.gz sampling_rate=0.001
# # snakemake --use-conda  --snakefile snakemake/fastq-subset.snakemake.py --cores 1 subsample_paired_fastqs --config forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.fq.gz sampling_rate=0.0001
forward = config["forward"]
reverse = config["reverse"]
sampling_rate = config["sampling_rate"]
sampling_rate_suffix = str(sampling_rate).split('.')[1]

if ".fq.gz" in forward:
    forward_subsampled = forward.replace(".fq.gz", "." + sampling_rate_suffix + ".fq.gz")
else:
    forward_subsampled = forward.replace(".fastq.gz", "." + sampling_rate_suffix + ".fastq.gz")

if ".fq.gz" in reverse:
    reverse_subsampled = reverse.replace(".fq.gz", "." + sampling_rate_suffix + ".fq.gz")    
else:
    reverse_subsampled = reverse.replace(".fastq.gz", "." + sampling_rate_suffix + ".fastq.gz")

rule subsample_paired_fastqs:
    conda:
        "../environment.yml"
    input:
        forward,
        reverse
    output:
        forward_subsampled,
        reverse_subsampled
    shell:
        """
        seqtk sample -s0 {forward} {sampling_rate} | pigz -c > {forward_subsampled}
        seqtk sample -s0 {reverse} {sampling_rate} | pigz -c > {reverse_subsampled}
        """