#!/bin/bash
# Step 5 in the pipeline
#
#SBATCH --job-name=cdhit
#SBATCH --output=cdhit.out
#
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1
#
module load cd-hit
#
cd /nfs/scistore08/kondrgrp/sovchinn/PEGP/projects/$1/
pwd

srun cd-hit -i filtered_hits.faa -o clustered90.faa -c 0.9 -n 5 -M 16000
