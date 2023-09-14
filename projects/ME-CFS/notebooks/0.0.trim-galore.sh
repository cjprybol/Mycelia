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
#SBATCH --job-name=trim-galore
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
CPU=1

# conda create -n trim_galore -c bioconda trim-galore
OUTPUT_DIRECTORY=$1
FORWARD=$2
REVERSE=$3

# sbatch trim-galore.sh $OUTPUT_DIRECTORY $FORWARD $REVERSE
# squeue -u $USER

conda run -n trim_galore trim_galore --suppress_warn --cores $CPU --output_dir $OUTPUT_DIRECTORY --paired $FORWARD $REVERSE

