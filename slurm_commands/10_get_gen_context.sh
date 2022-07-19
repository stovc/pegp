#!/bin/bash
#
#SBATCH --job-name=genome_context
#SBATCH --output=genome_context_out.txt
#
#SBATCH --time=60:00:00
#
#SBATCH --mem=128G
#64 was not sufficient
#
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
export OMP_NUM_THREADS=1
#
#run the respective binary through SLURM's srun

ml load python/3.9.5
pip install pandas

cd ../scripts
srun --cpu_bind=verbose python3 10_get_gen_context.py $1 $2
