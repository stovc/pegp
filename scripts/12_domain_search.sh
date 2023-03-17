#!/bin/bash
# Search protein domains in proteins from "clustered_full.faa" from project folder specified in $1 against against the Pfam-A.hmm database
# output - domain search report "hmmscan.tbl"

# args
PROJECT=$1
DOMAIN_DB="databases/Pfam/Pfam-A.hmm"

echo 12 started >> projects/"$PROJECT"/exit_log.txt
hmmscan --domtblout projects/"$PROJECT"/hmmscan.tbl --cut_ga "$DOMAIN_DB" projects/"$PROJECT"/clustered_full.faa

# write step number and exit status to the exit log
echo 12 $? >> projects/$1/exit_log.txt
