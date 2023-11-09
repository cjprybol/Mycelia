#!/bin/bash
#SBATCH --mail-user=cameron.prybol@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --time=12:00:00
##SBATCH --time-min=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# In Slurm a "cpu" corresponds to a hyperthread. There are 2 cpus per physical core.
# just shy of 2Gb of RAM per thread
# 256 threads, 128 physical cores per node
# 512 Gb of RAM per node
podman-hpc run --rm --volume /global:/global quay.io/biocontainers/virsorter:2.2.4--pyhdfd78af_0 virsorter run --db-dir $1 -w $2 -i $3 --min-length 1 -j 4 all