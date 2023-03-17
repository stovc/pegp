#!/bin/bash
# Construct maximum likelihood tree with 1000 ultrafast bootstrap replicates
# Input - trimed.fa
# Output - tree/tree.treefile - ML tree
#          tree/tree.treefile - consensus tree
#          tree/tree.iqtree - iqtree output

# args
PROJECT=$1

# log step number and started status
echo 15 started >> projects/"$PROJECT"/exit_log.txt

# make directory if does not exist
mkdir -p projects/"$PROJECT"/tree

# align
iqtree -s projects/"$PROJECT"/trimed.fa -B 1000 -T AUTO --prefix projects/"$PROJECT"/tree/tree

# write step number and exit status to the exit log
echo 15 $? >> projects/"$PROJECT"/exit_log.txt
