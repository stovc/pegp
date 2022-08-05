#!/bin/bash
#
#SBATCH --job-name=genome_context
#SBATCH --output=mk_blastp_df_out.txt
#
#SBATCH --time=2:00:00   # HOW LONG DO I NEED?
#
#SBATCH --mem=128G         # 64 was not sufficient
#
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
export OMP_NUM_THREADS=1
#
#run the respective binary through SLURM's srun

ml load python/3.9.5
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate /nfs/scistore08/kondrgrp/sovchinn/anaconda3/envs/pg

srun --cpu_bind=verbose python3 scripts/02_mk_blastp_df.py $1 $2
echo 2 $? >> projects/$1/exit_log.txt
