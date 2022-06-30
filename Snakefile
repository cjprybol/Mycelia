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

# https://papermill.readthedocs.io/en/latest/usage-cli.html
rule download_ncbi_reference_genomes:
    output:
        "data/genomes/done.txt"
    shell:
        # 10239
        # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10239
        # have the report title include a timestamp and parameters.yaml hash
        "papermill notebooks/scripts/01.download-ncbi-phage-genomes.ipynb reports/01.download-ncbi-phage-genomes.ipynb -p taxon_id 10239 -p data_dir $PWD/data/genomes"


# rule bwa_map:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq"
#     output:
#         "mapped_reads/{sample}.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"

# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} "
#         "-O bam {input} > {output}"