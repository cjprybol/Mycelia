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
#SBATCH --job-name=map-reads-to-human
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
CPU=8

FASTA=$1
FORWARD=$2
REVERSE=$3
OUTFILE=$4

# sbatch trim-galore.sh $OUTPUT_DIRECTORY $FORWARD $REVERSE
# squeue -u $USER
# conda create -n bwa-mem2 -c bioconda bwa-mem2
# conda create -n samtools -c bioconda samtools

conda run --live-stream --no-capture-output -n bwa-mem2 bwa-mem2 \
    mem \
    -t $CPU \
    ${FASTA} \
    ${FORWARD} \
    ${REVERSE} \
    | conda run --live-stream --no-capture-output -n samtools samtools view -bh \
    | conda run --live-stream --no-capture-output -n samtools samtools sort \
    > $OUTFILE
    
    
    

