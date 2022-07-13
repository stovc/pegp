#!/bin/bash
#
#SBATCH --job-name=blastp
#SBATCH --output=blastp.out
#
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1
#
module load ncbi-blast/2.11.0+
#
cd /nfs/scistore08/kondrgrp/sovchinn/PEGP/
echo "databases/$2/protein.faa"
echo "projects/$1/input.faa"
echo "projects/$1/blastp.xml"

srun blastp -query projects/$1/input.faa -db databases/$2/protein.faa -out projects/$1/blastp.xml -outfmt 5 -max_target_seqs 100000
