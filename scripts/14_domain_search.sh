#!/bin/bash
# Search protein domains in proteins from "clustered_full.faa" from project folder specified in $1 against the Pfam-A.hmm database
# Input:
#    - "clustered_full.faa" - fasta with full aa sequences of the hits
#    - "databases/Pfam/Pfam-A.hmm" - Pfam, target database for domain search
# Output:
#    - "hmmscan.tbl" - domain search report

# args
PROJECT=$1
DOMAIN_DB="databases/Pfam/Pfam-A.hmm"

echo 14 started >> projects/"$PROJECT"/exit_log.txt
hmmscan --domtblout projects/"$PROJECT"/hmmscan.tbl --cut_ga "$DOMAIN_DB" projects/"$PROJECT"/clustered_full.faa

# write step number and exit status to the exit log
echo 14 $? >> projects/$1/exit_log.txt
