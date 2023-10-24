#!/bin/bash
source sbatch-config.sh
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name=qc-short-reads
CPU=1
OUTPUT_DIRECTORY=$1
FORWARD=$2
REVERSE=$3

# sbatch qc-short-reads.sh $OUTPUT_DIRECTORY $FORWARD $REVERSE
# check progress
# squeue -u $USER 

if [ $(conda env list | grep "trim_galore" | wc -l) -eq 0 ]; then
    conda create -n trim_galore -c bioconda trim-galore
fi
conda run -n trim_galore trim_galore --suppress_warn --cores $CPU --output_dir $OUTPUT_DIRECTORY --paired $FORWARD $REVERSE