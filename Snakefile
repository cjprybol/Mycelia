# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

# snakemake --lint

import os

rule download_hmp_samples:
    input:
        "metadata/hmp_wgs_fastq/hmp_manifest_5555e11088.tsv"
        "metadata/hmp_wgs_fastq/hmp_manifest_metadata_5c2e51f49b.tsv"
        "/home/jovyan/rclone-mounts/drive.linked"
    output:
        "data/fastqs/donwloaded.done"
    shell:
        """

        """

rule download_kraken_indices:
    input:
        "/home/jovyan/rclone-mounts/drive.linked"
    output:
        "data/indices/"
    shell:
        """
        mkdir -p data/indices
        wget --directory-prefix data/kraken-indices https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20220607.tar.gz
        wget --directory-prefix data/kraken-indices https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20220607.tar.gz
        wget --directory-prefix data/kraken-indices https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20220607.tar.gz
        wget --directory-prefix data/kraken-indices https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20220607.tar.gz
        """

# FASTA_IDs=[file.name.replace(".fna", "") for file in os.scandir("data/genomes") if file.is_file() and file.name.endswith(".fna")]

# rule all:
#     input:
#         expand("data/genomes/{fasta_id}.fna.fai", fasta_id=FASTA_IDs)

# # snakemake index_genomes --cores all
# rule index_genomes:
#     input:
#         "data/genomes/{fasta_id}.fna"
#     output:
#         "data/genomes/{fasta_id}.fna.fai"
#     params:
#         fasta_id="{fasta_id}"
#     shell:
#         """
#         samtools faidx {params.fasta_id}
#         """

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
# snakemake download_all_ncbi_genomes --cores 1
root_taxon_id=1
ncbi_db = "refseq"
rule download_all_ncbi_genomes:
    input:
        "/home/jovyan/rclone-mounts/drive.mounted.linked"
    output:
        f"data/genomes/{root_taxon_id}/downloaded.done"
    shell:
        """
        papermill \
            --log-output \
            notebooks/scripts/download-ncbi-genomes.ipynb \
            reports/$(date "+%Y%m%d%H%M%S").download-ncbi-genomes.ipynb \
            -p ncbi_database {root_taxon_id} \
            -p data_dir $PWD/data \
            -p ncbi_database {ncbi_db}
        touch {output}
        """

# https://papermill.readthedocs.io/en/latest/usage-cli.html
# snakemake download_viral_ncbi_genomes --cores 1
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10239
viral_taxon_id=10239
ncbi_db = "refseq"
rule download_viral_ncbi_genomes:
    input:
        "/home/jovyan/rclone-mounts/drive.mounted.linked"
    output:
        f"data/genomes/{viral_taxon_id}/downloaded.done"
    shell:
        """
        papermill \
            --log-output \
            notebooks/scripts/download-ncbi-genomes.ipynb \
            reports/$(date "+%Y%m%d%H%M%S").download-ncbi-genomes.ipynb \
            -p ncbi_database {viral_taxon_id} \
            -p data_dir $PWD/data \
            -p ncbi_database {ncbi_db}
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

# snakemake mount_storage --cores 1
rule mount_storage:
    output:
        "/home/jovyan/rclone-mounts/storage.mounted"
    shell:
        """
        mkdir -p /home/jovyan/rclone-mounts/storage
        rclone mount mycelia: /home/jovyan/rclone-mounts/storage &
        sleep 5
        touch {output}
        """

# snakemake unmount_and_unlink_storage --cores 1 && snakemake link_storage --cores 1

# snakemake link_storage --cores 1
rule link_storage:
    input:
        "/home/jovyan/rclone-mounts/storage.mounted"
    output:
        "/home/jovyan/rclone-mounts/storage.mounted.linked"
    shell:
        """
        [ -d "/workspaces/$RepositoryName/data" ] && rm /workspaces/$RepositoryName/data
        mkdir -p /home/jovyan/rclone-mounts/storage
        ln -s /home/jovyan/rclone-mounts/storage/mycelia-storage /workspaces/$RepositoryName/data
        touch {output}
        """

# snakemake unmount_and_unlink_storage --cores 1
rule unmount_and_unlink_storage:
    input:
        mounted="/home/jovyan/rclone-mounts/storage.mounted",
        linked="/home/jovyan/rclone-mounts/storage.mounted.linked"
    shell:
        """
        fusermount -u /home/jovyan/rclone-mounts/storage || echo "storage not mounted"
        rm -f {input.mounted}
        rm -f {input.linked}
        """