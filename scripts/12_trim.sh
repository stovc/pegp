#!/bin/bash
# trim the alignment by removing columns that contain fraction of gaps lower than $THRESHOLD
# Input:
#    - "aligned.fa";
#    - [trimming threshold]
# Output:
#    - "trimed.fa" - trimmed alignment
#    - trim.html - triming report

# arg
PROJECT=$1
THRESHOLD=$3

# log step number and started status
echo 12 started >> projects/"$PROJECT"/exit_log.txt

# trim
trimal -in projects/"$PROJECT"/aligned.fa -out projects/"$PROJECT"/trimed.fa -htmlout projects/"$PROJECT"/trim.html -gt "$THRESHOLD"

# write step number and exit status to the exit log
echo 12 $? >> projects/"$1"/exit_log.txt
