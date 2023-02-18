# document
# https://snakemake.readthedocs.io/en/v6.0.3/executing/cli.html#visualization
# snakemake --snakefile snakemake/assemble.snakemake.py --cores 1 --forceall --dag | dot -Tpdf > snakemake/assemble.snakemake.py.dag.pdf

# lint
# snakemake --lint --snakefile snakemake/assemble.snakemake.py

import os

forward_reads = config["forward"]
reverse_reads = config["reverse"]
out_directory = config["out_directory"]

########################################################################################
# assemble into contigs with megahit
########################################################################################

# snakemake --use-conda  --snakefile snakemake/assemble.snakemake.py --cores all megahit_assemble --config forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.0001.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.0001.fq.gz out_directory=data/SRA/SRR6399459/
rule megahit_assemble:
    conda:
        "../environment.yml"
    input:
        forward_reads,
        reverse_reads
    # output:
    shell:
        """
        mkdir -p {out_directory}
        megahit \
            -1 {forward_reads} \
            -2 {reverse_reads} \
            -o {out_directory}/megahit
        """

########################################################################################
# assemble into contigs with spades
########################################################################################

#   -t <int>, --threads <int>   number of threads. [default: 16]
#   -m <int>, --memory <int>    RAM limit for SPAdes in Gb (terminates if exceeded). [default: 250]
# snakemake --use-conda  --snakefile snakemake/assemble.snakemake.py --cores all spades_assemble --config forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.0001.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.0001.fq.gz out_directory=data/SRA/SRR6399459/
rule spades_assemble:
    input:
        forward_reads,
        reverse_reads
    conda:
        "../environment.yml"
    shell:
        """
        spades.py \
            --phred-offset 33 \
            --meta \
            -1 {forward_reads} \
            -2 {reverse_reads} \
            -o {out_directory}/spades
        """