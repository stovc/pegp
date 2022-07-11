#!/bin/bash
#
#-------------------------------------------------------------
#example script for running a single-CPU serial job via SLURM
#-------------------------------------------------------------
#
#SBATCH --job-name=mkbdb
#SBATCH --output=mkbdb_out.txt
#
#Define the number of hours the job should run. 
#Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=120:00:00
#
#Define the amount of RAM used by your job in GigaBytes
#SBATCH --mem=16G
#
#SBATCH --no-requeue
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#for single-CPU jobs make sure that they use a single thread
export OMP_NUM_THREADS=1
#
#run the respective binary through SLURM's srun

ml load python/3.9.5
pip install pandas
pip install biopython

cd /nfs/scistore08/kondrgrp/sovchinn/PEGP
srun --cpu_bind=verbose python3 MakeBlastDatabaseInput.py base1 archaea bacteria viruses