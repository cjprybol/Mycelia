import os

################################################################################################
# DOWNLOAD KRAKEN DATABASES
################################################################################################

# archaea, bacteria, viral, plasmid, human1, UniVec_Core
# Standard plus protozoa & fungi
# Standard plus protozoa, fungi & plant
# capped at 8
# snakemake --snakefile Snakefile.py --cores 1 download_kraken_pluspfp8_db

out_directory = config["out_directory"]
kraken_pluspfp8 = "k2_pluspfp_08gb_20220607"

# snakemake --use-conda  --snakefile snakemake/download-databases.snakemake.py --cores 1 download_kraken_pluspfp8_db --config out_directory=data/kraken-databases
rule download_kraken_pluspfp8_db:
    conda:
        "../environment.yml"
    output:
        directory(out_directory + "/" + kraken_pluspfp8)
    shell:
        """
        wget --no-clobber --directory-prefix {out_directory} https://genome-idx.s3.amazonaws.com/kraken/{kraken_pluspfp8}.tar.gz
        mkdir -p {output}
        tar -xvzf {out_directory}/{kraken_pluspfp8}.tar.gz --directory {output}
        """

# done!
# capped at 16
kraken_pluspfp16 = "k2_pluspfp_16gb_20220607"

# snakemake --use-conda  --snakefile snakemake/download-databases.snakemake.py --cores 1 download_kraken_pluspfp16_db --config out_directory=data/kraken-databases
rule download_kraken_pluspfp16_db:
    conda:
        "../environment.yml"
    output:
        directory(out_directory + "/" + kraken_pluspfp16)
    shell:
        """
        wget --no-clobber --directory-prefix {out_directory} https://genome-idx.s3.amazonaws.com/kraken/{kraken_pluspfp16}.tar.gz
        mkdir -p {output}
        tar -xvzf {out_directory}/{kraken_pluspfp16}.tar.gz --directory {output}
        """
        
# # TODO
# # snakemake --use-conda  --snakefile snakemake/classify-reads.snakemake.py --cores 1 download_kraken_pluspfp_db --config out_directory=data/kraken-databases
# rule download_kraken_pluspfp_db:
#     output:
#         directory("data/kraken-databases/k2_pluspfp_20220607")
#     shell:
#         """
#         wget --no-clobber --directory-prefix data/kraken-databases/ https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20220607.tar.gz
#         mkdir -p {output}
#         tar -xvzf data/kraken-databases/k2_pluspfp_20220607.tar.gz --directory {output}
#         """

########################################################################################
# DOWNLOAD CENTRIFUGE DATABASES
# couldn't get to run
########################################################################################

# # https://benlangmead.github.io/aws-indexes/centrifuge
# # NCBI: nucleotide non-redundant sequences 	March, 2018
# #
# https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz

# couldn't get this to run
# https://benlangmead.github.io/aws-indexes/centrifuge
# NCBI: nucleotide non-redundant sequences 	March, 2018
#
# https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
# snakemake --snakefile Snakefile.py --cores 1 download_centrifuge_nt
# done!
# rule download_centrifuge_nt:
#     output:
#         directory("data/centrifuge-databases/nt_2018_3_3")
#     shell:
#         """
#         mkdir -p {output}
#         wget --no-clobber --directory-prefix data/centrifuge-databases/ https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz
#         tar -xvzf data/centrifuge-databases/nt_2018_3_3.tar.gz --directory {output}
#         """

########################################################################################
# DOWNLOAD KAIJU DATABASES
# https://github.com/bioinformatics-centre/kaiju
# dropped because of non-inclusive database
########################################################################################

# mkdir data/kaijudb
# cd data/kaijudb
# kaiju-makedb -s refseq