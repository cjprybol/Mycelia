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
#SBATCH --job-name=mmseqs_easy_search
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=16G
CPU=16
# ~ 8Gb RAM / CPU

## # LBL
## SBATCH --mail-user=cameron.prybol@gmail.com
## SBATCH --mail-type=ALL
## SBATCH --error=%x-%j.err
## SBATCH --output=%x-%j.out
## SBATCH --ntasks=1
## SBATCH --cpus-per-task=16
## SBATCH --job-name=mmseqs_easy_search
## SBATCH --partition=lr3
## SBATCH --qos=lr_normal

# NERSC
# #SBATCH --mail-user=cameron.prybol@gmail.com
# #SBATCH --mail-type=ALL
# #SBATCH --error=%x-%j.err
# #SBATCH --output=%x-%j.out
# #SBATCH --constraint=cpu
# #SBATCH --qos=shared
# #SBATCH --time=24:00:00
# ##SBATCH --time-min=01:00:00
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=128
# N_THREADS=128
# # just shy of 2Gb of RAM per thread
# # 256 threads/logical CPUs, 128 physical cores per node
# # 512 Gb of RAM per node

FASTA=$1
DATABASE=$2
OUTPATH=$3

# conda create -n mmseqs2 -c bioconda mmseqs2
# mkdir -p $HOME/workspace/mmseqs
# conda run --no-capture-output -n mmseqs2 mmseqs databases UniRef50 $HOME/workspace/mmseqs/UniRef50 $HOME/workspace/mmseqs/tmp

# kept getting OOO errors, so just running defaults
# --start-sens 1 \
# -s 7 \
# --sens-steps 3 \
# --threads $CPU \

# sbatch mmseqs-easy-search.sh $FASTA $DATABASE $OUTPATH
conda run --no-capture-output -n mmseqs2 \
mmseqs easy-search \
$FASTA \
$DATABASE \
${OUTPATH} \
${OUTPATH}_TMP \
--format-mode 4 \
--format-output query,qheader,target,theader,pident,fident,nident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,taxid

rm -rf ${OUTPATH}_TMP
            
           
            
            