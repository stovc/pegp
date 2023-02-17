#!/bin/bash
#
#SBATCH --job-name=blastp
#SBATCH --output=blastp.out
#
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1
#

srun bash scripts/01_blastp.sh $1 $2
echo 1 $? >> projects/$1/exit_log.txt
