#!/bin/bash

# Run ncbi-blast+ makeblastdb program with proteins from databases/$1/annotation.csv to make a blast database

INPUT=databases/$1/protein.faa
TAXID_MAP=databases/$1/taxid_map.txt

makeblastdb -in "$INPUT" -blastdb_version 5 -taxid_map "$TAXID_MAP" -dbtype prot -parse_seqids -title "$1"
