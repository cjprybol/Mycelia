# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

# snakemake --lint

import os

FASTA_IDs=[file.name.replace(".fna", "") for file in os.scandir("data/genomes") if file.is_file() and file.name.endswith(".fna")]

rule all:
    input:
        expand("data/genomes/{fasta_id}.fna.fai", fasta_id=FASTA_IDs)

# snakemake index_genomes --cores all
rule index_genomes:
    input:
        "data/genomes/{fasta_id}.fna"
    output:
        "data/genomes/{fasta_id}.fna.fai"
    params:
        fasta_id="{fasta_id}"
    shell:
        """
        samtools faidx {params.fasta_id}
        """

# # snakemake decompress_genomes --cores 1
# rule recompress_genomes:
#     input:
#         expand("{id}.fna.gz"), fasta=[file.name.replace(".fna.gz", "") for file in os.scandir("data/genomes") if file.is_file() and file.name.endswith(".fna.gz")])
#     output:
#         "data/genomes/{id}.fna.gz"
#     shell:
#         """
#         gzip -d {input}
#         """

# https://papermill.readthedocs.io/en/latest/usage-cli.html
# snakemake download_ncbi_reference_genomes --cores 1
rule download_ncbi_reference_genomes:
    input:
        "/home/jovyan/rclone-mounts/drive.linked"
    output:
        "data/genomes/downloaded.done"
    shell:
        # 10239
        # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10239
        # have the report title include a timestamp and parameters.yaml hash
        """
        papermill --log-output notebooks/scripts/01.download-ncbi-phage-genomes.ipynb reports/01.download-ncbi-phage-genomes.ipynb -p taxon_id 10239 -p data_dir $PWD/data/genomes -p ncbi_database refseq
        touch {output}
        """

# https://snakemake.readthedocs.io/en/v6.0.3/executing/cli.html#visualization
# snakemake document --cores 1
rule document:
    output:
        "dag.pdf"
    shell:
        """
        snakemake --forceall --dag | dot -Tpdf > dag.pdf
        """

# snakemake mount_google_drive --cores 1
rule mount_google_drive:
    output:
        "/home/jovyan/rclone-mounts/drive.mounted"
    shell:
        """
        mkdir -p /home/jovyan/rclone-mounts/drive
        rclone mount drive: /home/jovyan/rclone-mounts/drive &
        touch /home/jovyan/rclone-mounts/drive.mounted
        """
# snakemake unmount_google_drive --cores 1
rule unmount_google_drive:
    shell:
        """
        fusermount -u /home/jovyan/rclone-mounts/drive
        rm -f /home/jovyan/rclone-mounts/drive.mounted
        """

# snakemake link_data_directory --cores 1
rule link_data_directory:
    input:
        "/home/jovyan/rclone-mounts/drive.mounted"
    output:
        "/home/jovyan/rclone-mounts/drive.linked"
    shell:
        """
        mkdir -p /home/jovyan/rclone-mounts/drive/Projects/$RepositoryName
        ln -s /home/jovyan/rclone-mounts/drive/Projects/$RepositoryName /workspaces/$RepositoryName/data
        touch /home/jovyan/rclone-mounts/drive.linked
        """        

# snakemake unlink_data_directory --cores 1
rule unlink_data_directory:
    input:
        "/home/jovyan/rclone-mounts/drive.linked"
    shell:
        """
        rm /workspaces/$RepositoryName/data
        rm -f /home/jovyan/rclone-mounts/drive.linked
        """