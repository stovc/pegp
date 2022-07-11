#!/bin/bash
#
#SBATCH --job-name=mk_bdb
#SBATCH --output=mk_bdb.txt
#
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1
#
#load the respective software module you intend to use
module load ncbi-blast/2.11.0+
#
cd /nfs/scistore08/kondrgrp/sovchinn/PEGP/databases/base1
srun makeblastdb -in protein.faa -blastdb_version 5 -taxid_map taxid_map.txt -dbtype prot -parse_seqids -title "base1"