#!/bin/bash
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

srun bash scripts/05_cluster.sh $1
echo 5 $? >> projects/$1/exit_log.txt
