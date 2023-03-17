#!/bin/bash
# download and unzip Pfam and prepare it for hmmsearch

# download HMM profiles of the Pfam database
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -O databases/Pfam/Pfam-A.hmm.gz

# unzip
gunzip databases/Pfam/Pfam-A.hmm.gz

# make auxiliary binaries for hmmsearch
hmmpress databases/Pfam/Pfam-A.hmm
