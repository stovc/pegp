#!/bin/bash
# this step is reserved for testing which substitution model suits
# the alignment better for maximum likelihood inference
# currently does nothing

# log step number and started status
echo 13 started >> projects/$1/exit_log.txt

echo this step currently does nothing

# write step number and exit status to the exit log
echo 13 $? >> projects/$1/exit_log.txt
