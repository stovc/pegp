#!/bin/bash
# This step is reserved for testing which substitution model suits
# the alignment better for maximum likelihood inference.
# It currently does nothing.

# log step number and started status
echo 15 started >> projects/$1/exit_log.txt

echo this step currently does nothing

# write step number and exit status to the exit log
echo 15 $? >> projects/$1/exit_log.txt
