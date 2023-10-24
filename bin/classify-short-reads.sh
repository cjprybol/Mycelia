#!/bin/bash
source sbatch-config.sh
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --job-name=classify-short-reads
CPU=2
OUTPUT_DIRECTORY=$1
FORWARD=$2
REVERSE=$3

# sbatch qc-short-reads.sh $OUTPUT_DIRECTORY $FORWARD $REVERSE
# check progress
# squeue -u $USER 