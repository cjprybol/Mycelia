import os

################################################################################################
# CLASSIFY READS WITH KRAKEN
################################################################################################

kraken_db = config["kraken_db"]
forward = config["forward"]
reverse = config["reverse"]
out_directory = config["out_directory"]
sample_id = config["sample_id"]
out_report = out_directory + "/" + sample_id + "." + os.path.basename(kraken_db) + ".kraken-output.txt"
classified_out = out_directory + "/" + sample_id + "#." + os.path.basename(kraken_db) + ".classified.fastq"

# print(kraken_db)
# print(forward)
# print(reverse)
# print(out_directory)
# print(out_report)
# print(classified_out)


# Loading database information... done.
# 5666 sequences (1.67 Mbp) processed in 0.086s (3952.6 Kseq/m, 1163.88 Mbp/m).
#   790 sequences classified (13.94%)
#   4876 sequences unclassified (86.06%)

# snakemake --use-conda  --snakefile snakemake/classify-reads.snakemake.py --cores all classify_kraken --config forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.0001.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.0001.fq.gz kraken_db=data/kraken-databases/k2_pluspfp_08gb_20220607 out_directory=data/SRA/SRR6399459/kraken sample_id=SRR6399459


# Loading database information... done.
# 5666 sequences (1.67 Mbp) processed in 0.106s (3221.1 Kseq/m, 948.50 Mbp/m).
#   905 sequences classified (15.97%)
#   4761 sequences unclassified (84.03%)

# snakemake --use-conda  --snakefile snakemake/classify-reads.snakemake.py --cores all classify_kraken --config forward=data/SRA/SRR6399459/trim_galore/SRR6399459_1_val_1.0001.fq.gz reverse=data/SRA/SRR6399459/trim_galore/SRR6399459_2_val_2.0001.fq.gz kraken_db=data/kraken-databases/k2_pluspfp_16gb_20220607 out_directory=data/SRA/SRR6399459/kraken sample_id=SRR6399459
rule classify_kraken:
    input:
        forward,
        reverse
    conda:
        "../environment.yml"
    output:
        out_report
    shell:
        """
        kraken2 \
            --report-zero-counts \
            --use-names \
            --threads `nproc` \
            --db {kraken_db} \
            --output {output} \
            --report {out_report} \
            --gzip-compressed \
            --classified-out {classified_out} \
            --paired {input}
        """
        
########################################################################################
# classify reads with centriguge
########################################################################################
# couldn't get to work

# mkdir -p data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge
# centrifuge \
#     --threads `nproc` \
#     -x data/centrifuge-databases/nt_2018_3_3/nt \
#     -1 data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_1.fastq.gz \
#     -2 data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/SRR6476469_2.fastq.gz \
#     --al-conc-gz data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/SRR6476469.aligned.fastq.gz \
#     -S data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/classification.txt \
#     --report-file data/usa_mt-pleasant-research-farm_cornell-university_new-york/SRR6476469/centrifuge/centrifuge_report.tsv