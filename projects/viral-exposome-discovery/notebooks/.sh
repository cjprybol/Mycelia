#!/bin/bash
#SBATCH --mail-user=cameron.prybol@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=test
#SBATCH --partition=lr3
#SBATCH --qos=lr_normal

# I can only use lr3?
##SBATCH --partition=lr
##SBATCH --qos=lr_normal

# didn't work
##SBATCH --partition=cm1
##SBATCH --qos=cm1_normal

# didn't work
##SBATCH --partition=cf1
##SBATCH --qos=cf_lowprio+

# didn't work
##SBATCH --partition=es1
##SBATCH --qos=es_lowprio+

# partitions memory cores 16-64
# lr3 64
# lr4 64
# lr5 64/128
# lr6 96/128
# lr7 256
# lr_bigmem 1584
# cm1 256, 48
# cf1 192, 64
# es1 96/512, 8/64

# >= 4Gb/core

# sinfo – view the current status of the queues, e.g., sinfo
# sinfo -o "%12P %5N %10l %10a %20F"
# Allocated
# Idle
# Other
# Total
# PARTITION    NODEL TIMELIMIT  AVAIL      NODES(A/I/O/T)      
# cm1          n0000 infinite   up         2/11/1/14
# es1          n0000 infinite   up         10/30/7/47
# cf1          n0000 infinite   up         2/62/8/72
# lr3          n0000 infinite   up         156/1/175/332       
# lr4          n0000 infinite   up         23/0/118/141        
# lr5          n0000 infinite   up         145/0/47/192        
# lr6          n0000 infinite   up         326/18/66/410       
# lr7          n0000 infinite   up         60/32/8/100         
# lr_bigmem    n0272 infinite   up         0/2/0/2          


# $ sacctmgr show associations user=$USER | awk '{print $4}' | tail -n+3

# sbatch – submit a job to the batch queue system, e.g., sbatch myjob.sh
# sbatch lbl-slurm-test.sh "this is a test string"

# squeue – check the current jobs in the batch queue system, e.g., squeue
# squeue -u $USER

# sacct -j –format=JobID,JobName,MaxRSS,Elapsed

# srun – to run interactive jobs, e.g., srun –pty bash
# scancel – cancel a job, e.g., scancel 123

echo "$@"