#!/bin/bash

echo 5 started >> projects/$1/exit_log.txt
cd-hit -i projects/$1/filtered_hits.faa -o projects/$1/clustered90.faa -c 0.9 -n 5 -M 16000
