import os

########################################################################################
# map reads with bwa
########################################################################################

reference_fasta = config["reference_fasta"]
forward = config["forward"]
reverse = config["reverse"]
joint_identifier = os.path.basename(forward).replace(".fq.gz", "") + "__" + os.path.basename(reverse).replace(".fq.gz", "")
out_bam_base = os.path.dirname(forward) + "/" + joint_identifier + "__" + os.path.basename(reference_fasta)

# bowtie_indices = "data/taxon_10239.genbank/10239.deduped.fna.*.bt2"

bwa_mem_out_bam = f"{out_bam_base}.bwa.bam"
bowtie2_out_bam = f"{out_bam_base}.bowtie2.bam"

# print(reference_fasta)
# print(forward)
# print(reverse)
# print(joint_identifier)
# print(out_sam_base)

# snakemake --use-conda  --snakefile snakemake/map-reads.snakemake.py --cores all bwa_index_reference_fasta --config reference_fasta=data/taxon_10239.genbank/10239.fna
rule bwa_index_reference_fasta:
    conda:
        "../environment.yml"
    input:
        reference_fasta
    output:
        reference_fasta + ".amb",
        reference_fasta + ".ann",
        reference_fasta + ".bwt",
        reference_fasta + ".pac"
    shell:
        """
        bwa index {input}
        """

# snakemake --use-conda  --snakefile snakemake/map-reads.snakemake.py --cores all bwa_mem_align --config reference_fasta=data/taxon_10239.genbank/10239.deduped.fna forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.0001.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.0001.fq.gz
rule bwa_mem_align:
    conda:
        "../environment.yml"
    output:
        bwa_mem_out_bam
    input:
        reference_fasta,
        reference_fasta + ".amb",
        reference_fasta + ".ann",
        reference_fasta + ".bwt",
        reference_fasta + ".pac",
        forward,
        reverse
    shell:
        """
        bwa \
            mem \
            -t `nproc` \
            {reference_fasta} \
            {forward} \
            {reverse} \
            | samtools view -bhF 3844 \
            > {output}
        """

########################################################################################
# map reads with bowtie2/hisat
########################################################################################


# snakemake --use-conda --snakefile snakemake/map-reads.snakemake.py --cores all bowtie2_index_reference_fasta --config reference_fasta=data/taxon_10239.genbank/10239.deduped.fna forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.0001.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.0001.fq.gz
rule bowtie2_index_reference_fasta:
    conda:
        "../environment.yml"
    input:
        reference_fasta
    output:
        reference_fasta + ".1.bt2"
        reference_fasta + ".rev.1.bt2"
    shell:
        """
        bowtie2-build --threads `nproc` {input} {input}
        """


# bowtie2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
# --al <path>, --al-gz <path>, --al-bz2 <path>

# snakemake --use-conda  --snakefile snakemake/map-reads.snakemake.py --cores all bwa_index_reference_fasta --config reference_fasta=data/taxon_10239.genbank/10239.fna
# rule bwa_index_reference_fasta:
#     conda:
#         "../environment.yml"
#     input:
#         reference_fasta
#     output:
#         reference_fasta + ".amb",
#         reference_fasta + ".ann",
#         reference_fasta + ".bwt",
#         reference_fasta + ".pac"
#     shell:
#         """
#         bwa index {input}
#         """

# snakemake --use-conda  --snakefile snakemake/map-reads.snakemake.py --cores all bwa_mem_align --config reference_fasta=data/taxon_10239.genbank/10239.deduped.fna forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.0001.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.0001.fq.gz
# rule bwa_mem_align:
#     conda:
#         "../environment.yml"
#     output:
#         bwa_mem_out_bam
#     input:
#         reference_fasta,
#         reference_fasta + ".amb",
#         reference_fasta + ".ann",
#         reference_fasta + ".bwt",
#         reference_fasta + ".pac",
#         forward,
#         reverse
#     shell:
#         """
#         bwa \
#             mem \
#             -t `nproc` \
#             {reference_fasta} \
#             {forward} \
#             {reverse} \
#             | samtools view -bhF 3844 \
#             > {output}
#         """