#!/bin/bash
#SBATCH --mail-user=cameron.prybol@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --constraint=cpu
#SBATCH --qos=shared
# we had a timeout with 12, so doubling to 24 (max)
#SBATCH --time=24:00:00
##SBATCH --time-min=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# we had a timeout with 8, so doubling to 16
#SBATCH --cpus-per-task=16
# just shy of 2Gb of RAM per thread
# 256 threads/logical CPUs, 128 physical cores per node
# 512 Gb of RAM per node
ASSEMBLED_FASTA=$1
OUTDIR=$2
GENOMAD_DB_DIR=$3
# this should match cpus-per-task above
N_THREADS=16
conda run --no-capture-output -n genomad genomad end-to-end --threads $N_THREADS --splits $N_THREADS $ASSEMBLED_FASTA $OUTDIR $GENOMAD_DB_DIR