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
import datetime

invoked_datetime = datetime.datetime.now()
invoked_timestamp = f"{invoked_datetime.year}{invoked_datetime.month}{invoked_datetime.day}{invoked_datetime.hour}{invoked_datetime.second}"

taxon = str(config["taxon"])
source = config["source"]
dehydrated_taxon_directory = "data/taxon_" + taxon + "." + source
dehydrated_taxon_zip = dehydrated_taxon_directory + '.zip'

################################################################################################
# BY TAXON ID GENBANK
# Collecting 50,633 genome accessions [===========================>--------------------]  59% 30000/50633
################################################################################################

# snakemake --use-conda --snakefile snakemake/0.ncbi-datasets.snakemake.py --cores 1 download_dehydrated_taxon --config taxon="10239" source="genbank"
rule download_dehydrated_taxon:
    output:
        dehydrated_taxon_zip
    resources:
        tmpdir="data/"
    log:
        f"snakemake/logs/{invoked_timestamp}.log"
    conda:
        "../environment.yml"
    shell:
        """
        datasets download genome taxon {taxon} --dehydrated --assembly-source {source} --filename {output}
        """

# snakemake --use-conda --snakefile snakemake/0.ncbi-datasets.snakemake.py --cores 1 unzip_taxon --config taxon="10239" source="genbank"
rule unzip_taxon:
    input:
        dehydrated_taxon_zip
    resources:
        tmpdir="data/"
    log:
        f"snakemake/logs/{invoked_timestamp}.log"
    conda:
        "../environment.yml"
    output:
        directory(dehydrated_taxon_directory)
    shell:
        """
        unzip -d {output} {input}
        """
        
# snakemake --use-conda --snakefile snakemake/0.ncbi-datasets.snakemake.py --cores 1 rehydrate_taxon --config taxon="10239" source="genbank"
rehydrated_taxon_output = dehydrated_taxon_directory + "/rehydrated.done"
rule rehydrate_taxon:
    output:
        rehydrated_taxon_output
    resources:
        tmpdir="data/"
    conda:
        "../environment.yml"
    shell:
        """
        datasets rehydrate --directory {dehydrated_taxon_directory}
        touch {output}
        """

# snakemake --use-conda --snakefile snakemake/0.ncbi-datasets.snakemake.py --cores 1 merge_fastas --config taxon="10239" source="genbank"
joint_fasta_outfile = dehydrated_taxon_directory + "/" + taxon + ".fna"
rule merge_fastas:
    input:
        rehydrated_taxon_output
    output:
        joint_fasta_outfile
    conda:
        "../environment.yml"
    shell:
        """
        find \
            {dehydrated_taxon_directory}/ncbi_dataset/data \
            -maxdepth 3 \
            -type f \
            -name '*_genomic.fna' \
            -print0 \
            | xargs -0 cat > {output}
        """