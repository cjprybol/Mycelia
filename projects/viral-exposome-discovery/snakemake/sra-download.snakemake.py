# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

# document
# https://snakemake.readthedocs.io/en/v6.0.3/executing/cli.html#visualization
# snakemake --snakefile snakemake/0.data-acquisition.snakemake.py --cores 1 --forceall --dag | dot -Tpdf > snakemake/0.data-acquisition.snakemake.py.dag.pdf

# lint
# snakemake --lint --snakefile snakemake/0.data-acquisition.snakemake.py

import os
# import datetime

# invoked_datetime = datetime.datetime.now()
# invoked_timestamp = f"{invoked_datetime.year}{invoked_datetime.month}{invoked_datetime.day}{invoked_datetime.hour}{invoked_datetime.second}"
    
################################################################################################
# DOWNLOAD DATA
################################################################################################

# download IMG/VR
# https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=IMG_VR
# need globus to do automatically
# don't do that, just download, upload to google drive, and then pull from that as a step

# converting jsonl to tables
# dataformat tsv genome --inputfile human/ncbi_dataset/data/assembly_data_report.jsonl

srr_identifier = str(config["srr_identifier"])
base_directory = config["directory"]
srr_directory = base_directory + "/" + srr_identifier
sra_lite_file = srr_directory + "/" + srr_identifier + ".sralite"
forward_reads_gz = srr_directory + "/" + srr_identifier + "_1.fastq.gz"
reverse_reads_gz = srr_directory + "/" + srr_identifier + "_2.fastq.gz"
        
# snakemake --use-conda --snakefile snakemake/0.sra-download.snakemake.py --cores all fasterq_dump_direct --config srr_identifier="SRR6399460" directory="data/SRA"
rule fasterq_dump_direct:
    output:
        forward_reads_gz,
        reverse_reads_gz
    conda:
        "../environment.yml"
    # resources:
    #     tmpdir="data/"
    # log:
    #     f"snakemake/logs/{invoked_timestamp}.log"
    shell:
        """
        fasterq-dump \
            --outdir {srr_directory}\
            --mem 1G \
            --split-3 \
            --threads `nproc` \
            --progress \
            --skip-technical \
            {srr_identifier}
        pigz {srr_directory}/*.fastq
        """
        
# alternatively to the above, can run both of these
# I don't like this approach as much because it leaves
# more intermediate files
# snakemake --use-conda --snakefile snakemake/0.sra-download.snakemake.py --cores 1 sra_prefetch --config srr_identifier="SRR6399459" directory="data/SRA"
# rule sra_prefetch:
#     output:
#         sra_lite_file
#     conda:
#         "../environment.yml"
#     resources:
#         tmpdir="data/"
#     log:
#         f"snakemake/logs/{invoked_timestamp}.log"
#     shell:
#         """
#         prefetch \
#             --max-size u \
#             --progress \
#             --output-directory {base_directory} \
#             {srr_identifier}
#         """

# snakemake --snakefile snakemake/0.sra-download.snakemake.py --cleanup-metadata data/SRA/SRR6399459/SRR6399459_2.fastq.gz --config srr_identifier="SRR6399459" directory="data/SRA"
# snakemake --use-conda --snakefile snakemake/0.sra-download.snakemake.py --cores all fasterq_dump_prefetch --config srr_identifier="SRR6399459" directory="data/SRA"
# rule fasterq_dump_prefetch:
#     input:
#         sra_lite_file
#     output:
#         forward_reads_gz,
#         reverse_reads_gz
#     conda:
#         "../environment.yml"
#     resources:
#         tmpdir="data/"
#     log:
#         f"snakemake/logs/{invoked_timestamp}.log"
#     shell:
#         """
#         fasterq-dump \
#             --outdir {srr_directory}\
#             --mem 1G \
#             --split-3 \
#             --threads `nproc` \
#             --progress \
#             --skip-technical \
#             {input}
#         pigz {srr_directory}/*.fastq
#         """