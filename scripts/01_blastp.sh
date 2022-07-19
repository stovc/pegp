#!/bin/bash
module load ncbi-blast/2.11.0+
#
echo "databases/$2/protein.faa"
echo "projects/$1/input.faa"
echo "projects/$1/blastp.xml"

blastp -query projects/$1/input.faa -db databases/$2/protein.faa -out projects/$1/blastp.xml -outfmt 5 -max_target_seqs 100000
echo 1 $? >> exit_log.txt
