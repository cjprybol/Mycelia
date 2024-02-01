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
#SBATCH --job-name=megahit
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
CPU=16
# ~ 8Gb RAM / CPU

# PARTITION    NODEL TIMELIMIT  AVAIL      NODES(A/I/O/T)      
# batch*       cfx12 14-00:00:0 up         77/0/4/81           
# interactive  smsh1 1-00:00:00 up         4/0/0/4 

OUTPUT_DIRECTORY=$1
FORWARD=$2
REVERSE=$3

# sbatch 1.0.megahit.sh $OUTPUT_DIRECTORY $FORWARD $REVERSE
# squeue -u $USER

# module load megahit
# megahit --num-cpu-threads $CPU -o $OUTPUT_DIRECTORY -1 $FORWARD -2 $REVERSE 
conda run --live-stream --no-capture-output -n megahit megahit --num-cpu-threads $CPU -o $OUTPUT_DIRECTORY -1 $FORWARD -2 $REVERSE 