import os

########################################################################################
# map assembled contigs with minimap
########################################################################################

# map contigs to reference genomes

reference_fasta = config["reference_fasta"]
query_fasta = config["query_fasta"]

# snakemake --snakefile Snakefile.py --cores 1 minimap_megahit_assembled_contigs_SRR6399459_100k
rule minimap_assembled_contigs:
    # log:
    #     f"snakemake/logs/{invoked_timestamp}.log"
    conda:
        "environment.yml"
    shell:
        """
        minimap2 -ax asm5 data/taxon_10239.genbank/joint.fna data/exposome/SRR6399459/megahit_100k/final.contigs.fa > data/exposome/SRR6399459/megahit_100k/final.contigs.fa.joint.fna.sam
        """

# snakemake --snakefile Snakefile.py --cores 1 samtools_filter_minimap_megahit_assembled_contigs_SRR6399459_100k
rule samtools_filter:
    # log:
        # f"snakemake/logs/{invoked_timestamp}.log"
    conda:
        "environment.yml"
    shell:
        """
        samtools \
            view \
            --with-header \
            --excl-flags 3844 \
            data/exposome/SRR6399459/megahit_100k/final.contigs.fa.joint.fna.sam \
            > data/exposome/SRR6399459/megahit_100k/final.contigs.fa.joint.fna.sam.filtered.sam
        """