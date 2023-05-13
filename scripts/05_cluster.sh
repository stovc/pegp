#!/bin/bash
# cluster sequences into clusters of 90% and more similarity and pick representatives using cd-hit
# input - filtered_hits.faa
# output:
# - clustered90.faa - fasta file with representative sequences
# - clustered90.faa - description of all computed clusters

# arg
PROJECT=$1

# log step number and started status
echo 5 started >> projects/"$PROJECT"/exit_log.txt

# cluster sequences
cd-hit -i projects/"$PROJECT"/filtered_hits.faa -o projects/"$PROJECT"/clustered90.faa -c 0.9 -n 5 -M 16000

# log clustering program version
cd-hit | cat > projects/"$PROJECT"/clustering_log.txt

# write step number and exit status to the exit log
echo 5 $? >> projects/$1/exit_log.txt
