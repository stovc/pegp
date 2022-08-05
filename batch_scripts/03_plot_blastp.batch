#!/bin/bash
#
#SBATCH --job-name=genome_context
#SBATCH --output=plot_blastp_out.txt
#
#SBATCH --time=2:00:00   # HOW LONG DO I NEED?
#
#SBATCH --mem=32G
#
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
export OMP_NUM_THREADS=1
#
#run the respective binary through SLURM's srun

source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate /nfs/scistore08/kondrgrp/sovchinn/anaconda3/envs/pg

srun --cpu_bind=verbose python3 scripts/03_plot_blastp.py $1 $3 $4 $5
echo 3 $? >> projects/$1/exit_log.txt
