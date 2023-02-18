# https://papermill.readthedocs.io/en/latest/usage-cli.html
# run `snakemake some_target --delete-all-output` to clean outputs from target
# run `snakemake --delete-all-output --cores all` to clean all outputs
# run `snakemake --cores all` to run whole workflow
# run `snakemake some_target --cores all` to run up to and through target

# snakemake --lint

import os

# conda activate sars-cov2-pangenome-analysis

# to specificy conda env as part of the rule
# conda:
#     "environment.yaml"

# https://snakemake.readthedocs.io/en/v6.0.3/executing/cli.html#visualization
# snakemake document --cores 1
rule document:
    output:
        "dag.pdf"
    shell:
        """
        snakemake --forceall --dag | dot -Tpdf > dag.pdf
        """

################################################################################################
# FULL COVID DATASET
# holding off on this for now until I finish working with the smaller subsets
################################################################################################

# # snakemake --snakefile Snakefile.py --cores 1 download_ncbi_dataset_zip
# # very large
# rule download_ncbi_dataset_zip:
#     output:
#         "data/sars-cov-2.zip"
    # resources:
    #     tmpdir="data/"
#     shell:
#         """
#         datasets download virus genome taxon sars-cov-2 --filename {output}
#         """

# # raw export was too large to work with on a 1Tb machine
# # snakemake --snakefile Snakefile.py --cores 1 unzip_ncbi_dataset_zip
# rule unzip_ncbi_dataset_zip:
#     input:
#         "data/sars-cov-2.zip"
#     output:
#         "data/sars-cov-2"
    # resources:
    #     tmpdir="data/"
#     shell:
#         """
#         unzip -d {output} {input}
#         """

# ################################################################################################
# # ANNOTATED AND COMPLETE ONLY
# couldn't get this to download successfully to save my life
# ~ 40Gb zip file
# ################################################################################################

# # snakemake --snakefile Snakefile.py --cores 1 download_covid_dataset_annotated_complete
# rule download_covid_dataset_annotated_complete:
#     output:
#         "data/sars-cov-2.annotated.complete.zip"
#     resources:
#         tmpdir="data/"
#     shell:
#         """
#         datasets download virus genome taxon sars-cov-2 --filename {output} --annotated --complete-only
#         """
        
# # snakemake --snakefile Snakefile.py --cores 1 unzip_covid_dataset_annotated_complete
# rule unzip_covid_dataset_annotated_complete:
#     input:
#         "data/sars-cov-2.annotated.complete.zip"
#     output:
#         "data/sars-cov-2.annotated.complete"
#     resources:
#         tmpdir="data/"
#     shell:
#         """
#         unzip -d {output} {input}
#         """

################################################################################################
# ANNOTATED AND COMPLETE REFSEQ ONLY
################################################################################################

# don't need to run
# datasets summary virus genome taxon sars-cov-2 --annotated --complete-only --refseq > dataset-summary.json
# because the data_report.jsonl in the download zip is the same content

# snakemake --snakefile Snakefile.py --cores 1 download_covid_dataset_annotated_complete_refseq
rule download_covid_dataset_annotated_complete_refseq:
    output:
        "data/sars-cov-2.annotated.complete.refseq.zip"
    resources:
        tmpdir="data/"
    shell:
        """
        datasets download virus genome taxon sars-cov-2 --filename {output} --annotated --complete-only --refseq
        """

# snakemake --snakefile Snakefile.py --cores 1 unzip_covid_dataset_annotated_complete_refseq
rule unzip_covid_dataset_annotated_complete_refseq:
    input:
        "data/sars-cov-2.annotated.complete.refseq.zip"
    output:
        "data/sars-cov-2.annotated.complete.refseq"
    resources:
        tmpdir="data/"
    shell:
        """
        unzip -d {output} {input}
        """
        
################################################################################################
# BY TAXON ID GENBANK
# Collecting 92 genome accessions [================================================] 100% 92/92
# Found 429 files for rehydration
################################################################################################

# snakemake --snakefile Snakefile.py --cores 1 download_taxon_2697049_genbank_dehydrated
rule download_taxon_2697049_genbank_dehydrated:
    output:
        "data/taxon_2697049.genbank.zip"
    resources:
        tmpdir="data/"
    shell:
        """
        datasets download genome taxon 2697049 --dehydrated --assembly-source genbank --filename {output}
        """

# snakemake --snakefile Snakefile.py --cores 1 unzip_taxon_2697049_genbank
rule unzip_taxon_2697049_genbank:
    input:
        "data/taxon_2697049.genbank.zip"
    resources:
        tmpdir="data/"
    output:
        directory("data/taxon_2697049.genbank")
    shell:
        """
        unzip -d {output} {input}
        """
        
# snakemake --snakefile Snakefile.py --cores 1 rehydrate_taxon_2697049_genbank
rule rehydrate_taxon_2697049_genbank:
    input:
        directory("data/taxon_2697049.genbank")
    resources:
        tmpdir="data/"
    shell:
        """
        datasets rehydrate --directory data/taxon_2697049.genbank
        """

# ################################################################################################
# # BY TAXON ID
# # this was effectively the same as genbank, let's just use that
# Collecting 93 genome accessions [================================================] 100% 93/93
# Found 434 files for rehydration
# ################################################################################################

# # snakemake --snakefile Snakefile.py --cores 1 download_taxon_2697049_dehydrated
# rule download_taxon_2697049_dehydrated:
#     output:
#         "data/taxon_2697049.zip"
#     resources:
#         tmpdir="data/"
#     shell:
#         """
#         datasets download genome taxon 2697049 --dehydrated --filename {output}
#         """

# # snakemake --snakefile Snakefile.py --cores 1 unzip_taxon_2697049
# rule unzip_taxon_2697049:
#     input:
#         "data/taxon_2697049.zip"
#     resources:
#         tmpdir="data/"
#     output:
#         directory("data/taxon_2697049")
#     shell:
#         """
#         unzip -d {output} {input}
#         """
        
# # snakemake --snakefile Snakefile.py --cores 1 rehydrate_taxon_2697049
# rule rehydrate_taxon_2697049:
#     input:
#         directory("data/taxon_2697049")
#     resources:
#         tmpdir="data/"
#     shell:
#         """
#         datasets rehydrate --directory {input}
#         """

# snakemake --snakefile Snakefile.py --cores 1 unzip_taxon_2697049 && snakemake --snakefile Snakefile.py --cores 1 rehydrate_taxon_2697049
        
        
# rule download_ncbi_tax_dump:
#     output:
#     shell:


# rule build_genome_from_dataset:
#     input:
#         "dataset_directory:"
#     shell:
#         """
#         # julia --project="./Project.toml" -e 'import Pkg; Pkg.instantiate()'
#         # julia --project="./Project.toml" -e 'import Mycelia; println(pathof(Mycelia))'
#         # papermill path_of_mycelia/build_pangenome.ipynb --data_directory={input}
#         """

# dataformat tsv virus-genome --inputfile data/sars-cov-2.annotated.complete.refseq/ncbi_dataset/data/data_report.jsonl

# this doesn't work for biosample
# dataformat tsv virus-genome --inputfile data/sars-cov-2.annotated.complete.refseq/ncbi_dataset/data/biosample.jsonl

# dataformat tsv genome --inputfile human/ncbi_dataset/data/assembly_data_report.jsonl

# download SRA runs
# https://www.ncbi.nlm.nih.gov/sra/?term=txid2697049%5BOrganism:noexp%5D%20NOT%200[Mbases

# snakemake --snakefile Snakefile.py --cores 1 rehydrate_taxon_2697049_genbank
rule download_metadata_from_google_drive:
    output:
        "metadata/sequences.csv.gz"
    resources:
        tmpdir="data/"
    shell:
        """
        rclone copy google_drive:Projects/sars-cov2-pangenome-analysis/metadata ./metadata
        """
        
rule unzip_sequences_table:
    shell:
        """
        gzip -d ./metadata/sequences.csv.gz
        """
        
# TODO:
# downloading over API using the accessions in the above sequences table is working
# join the sequences table with the kmers
# genome metadata = node
# 