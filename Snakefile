# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

# rule hello_world:
#     output:
#         "data/test.txt"
#     shell:
#         "mkdir -p data && echo 'this is a test' > data/test.txt"

# https://snakemake.readthedocs.io/en/v6.0.3/executing/cli.html#visualization
rule document:
    output:
        "dag.pdf"
    shell:
        """
        snakemake --forceall --dag | dot -Tpdf > dag.pdf
        """

# snakemake initialize --cores 1
rule initialize:
    output:
        "/home/jovyan/rclone-mounts/drive.unmounted"
    shell:
        """
        touch /home/jovyan/rclone-mounts/drive.unmounted
        """

# snakemake mount_google_drive --cores 1
rule mount_google_drive:
    input:
        "/home/jovyan/rclone-mounts/drive.unmounted"
    output:
        "/home/jovyan/rclone-mounts/drive.mounted"
    shell:
        """
        mkdir -p /home/jovyan/rclone-mounts/drive
        rclone mount drive: /home/jovyan/rclone-mounts/drive &
        mv /home/jovyan/rclone-mounts/drive.unmounted /home/jovyan/rclone-mounts/drive.mounted
        """
# snakemake unmount_google_drive --cores 1
rule unmount_google_drive:
    input:
        "/home/jovyan/rclone-mounts/drive.mounted"
    output:
        "/home/jovyan/rclone-mounts/drive.unmounted"
    shell:
        """
        fusermount -u /home/jovyan/rclone-mounts/drive
        mv /home/jovyan/rclone-mounts/drive.mounted /home/jovyan/rclone-mounts/drive.unmounted
        """

# snakemake link_data_directory --cores 1
rule link_data_directory:
    input:
        "/home/jovyan/rclone-mounts/drive.mounted"
    output:
        "/home/jovyan/rclone-mounts/drive.linked"
    shell:
        """
        ln -s /workspace/Mycelia/data /home/jovyan/rclone-mounts/drive
        touch /home/jovyan/rclone-mounts/drive.linked
        """        

# snakemake unlink_data_directory --cores 1
rule unlink_data_directory:
    input:
        "/home/jovyan/rclone-mounts/drive.linked"
    output:
        "/home/jovyan/rclone-mounts/drive.unlinked"
    shell:
        """
        rm /workspace/Mycelia/data
        mv /home/jovyan/rclone-mounts/drive.linked /home/jovyan/rclone-mounts/drive.unlinked
        """

# https://papermill.readthedocs.io/en/latest/usage-cli.html
# snakemake download_ncbi_reference_genomes --cores 1
rule download_ncbi_reference_genomes:
    input:
        "/home/jovyan/rclone-mounts/drive.mounted"
    output:
        "data/genomes/done.txt"
    shell:
        # 10239
        # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10239
        # have the report title include a timestamp and parameters.yaml hash
        "papermill notebooks/scripts/01.download-ncbi-phage-genomes.ipynb reports/01.download-ncbi-phage-genomes.ipynb -p taxon_id 10239 -p data_dir $PWD/data/genomes"
