#!/bin/bash
# blast "input.faa" from project folder specified in $1 against a database specified in $2
# input - "input.faa" - search query fasta
# output - "blastp.xml" - BLASTP search report
# e-value cutoff - 0.05; out format - xml; maximum reported hits - as much as possible

# args
PROJECT=$1
DATABASE=$2

# log step number and started status
echo 1 started >> projects/"$PROJECT"/exit_log.txt

#blast
hmmsearch -o projects/"$PROJECT"/hits.txt projects/"$PROJECT"/input.faa databases/"$DATABASE"/protein.faa

# write step number and exit status to the exit log
echo 1 $? >> projects/$PROJECT/exit_log.txt
