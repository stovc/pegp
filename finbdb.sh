#!/bin/bash
#
#SBATCH --job-name=finbdb
#SBATCH --output=finbdb_out.txt
#
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#
#SBATCH --no-requeue
#
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#for single-CPU jobs make sure that they use a single thread
export OMP_NUM_THREADS=1

ml load python/3.9.5
pip install pandas
pip install biopython

cd /nfs/scistore08/kondrgrp/sovchinn/PEGP
srun --cpu_bind=verbose python3 FinishBlastDatabaseInput.py base1