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