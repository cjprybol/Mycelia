#!/bin/bash
# SCG3 setup
#SBATCH --mail-user=cameron.prybol@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --partition=batch
#SBATCH --account=mpsnyder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=bwa-index-human-genome
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
CPU=1

# sbatch 0.1.bwa-index.sh
# conda create -n bwa-mem2 -c bioconda bwa-mem2

conda run --no-capture-output -n bwa-mem2 bwa-mem2 index $HOME/workspace/Mycelia/projects/ME-CFS/reference_data/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna
