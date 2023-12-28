#!/bin/bash
#SBATCH --mail-user=cameron.prybol@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# LBL
#SBATCH --job-name=prodigal
#SBATCH --partition=lr3
#SBATCH --qos=lr_normal

# FASTA=$1
# # conda create -n prodigal -c bioconda prodigal
# PRODIGAL_BASE="${FASTA}_prodigal/$(basename ${FASTA})"

# conda run --no-capture-output -n prodigal \
#     prodigal \
#         -m \
#         -p meta \
#         -i ${FASTA} \
#         -f gff \
#         -o ${FASTA}.prodigal.gff \
#         -a ${FASTA}.prodigal.faa \
#         -d ${FASTA}.prodigal.fna \
#         -s ${FASTA}.prodigal.all_potential_gene_scores.txt
