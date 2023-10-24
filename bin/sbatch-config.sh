#!/bin/bash

current_host=$(hostname)
echo "running on ${current_host}"

if [[ $current_host == *"stanford.edu" ]]; then
# if echo "$current_host" | grep -q "stanford\.edu$"; then
    echo "using slurm parameters for Stanford"
    source sbatch-scg3.sh
elif [[ "${string}" == *"login" ]]; then
# elif echo "$current_host" | grep -q "^login"; then
    echo "using slurm parameters for NERSC"
    source sbatch-nersc.sh
else
    echo "are you on Lawrencium?"
    source sbatch-lawrencium.sh
fi