#!/bin/bash
# trim the alignment
# remove columns that contain >= 50% of gaps

# arg
PROJECT=$1

# log step number and started status
echo 11 started >> projects/$PROJECT/exit_log.txt

# trim
trimal -in projects/"$PROJECT"/aligned.fa -out projects/"$PROJECT"/trimed50.fa -htmlout projects/"$PROJECT"/trim50.html -gt 0.5

# write step number and exit status to the exit log
echo 11 $? >> projects/$1/exit_log.txt
