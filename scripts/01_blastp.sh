#!/bin/bash
# blast input.faa from project folder specified in $1 against a database specified in $2

# log step number and started status
echo 1 started >> projects/$1/exit_log.txt

#blast
blastp -query projects/$1/input.faa -db databases/$2/protein.faa -out projects/$1/blastp.xml -evalue 0.05 -outfmt 5 -max_target_seqs 100000

# write step number and exit status to the exit log
echo 1 $? >> projects/$1/exit_log.txt
