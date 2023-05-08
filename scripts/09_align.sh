#!/bin/bash
# align clustered and genome-consistent hits by MAFFT in the L-ins-i mode
# input - clustered90-gc.faa
# output - aligned.fa

# args
PROJECT=$1

# log step number and started status
echo 9 started >> projects/"$PROJECT"/exit_log.txt

# align
linsi --thread -1 projects/"$PROJECT"/clustered90-gc.faa > projects/"$PROJECT"/aligned.fa

# write step number and exit status to the exit log
echo 9 $? >> projects/"$PROJECT"/exit_log.txt
