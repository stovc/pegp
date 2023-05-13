"""Filter HMMER hits by parameters: query coverage, e-value.

- input:
    - "hits.txt" - hmmer report. needed to extract aligned parts of the hits
    - "hits_df.csv" - csv with annotations to subset
    - parameters with filtering thresholds: [query coverage], [e-value]

- output:
    - 'filtered_hits.faa' - fasta file with filtered hits containing aligned parts of the proteins

- add 'filtered' column to hits_df.csv. Value is True, if the hit successfully passed filtering
"""

import sys
from pathlib import Path
from Bio import Seq, SeqIO
from Bio import SeqFeature
from Bio import SearchIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gc
import subprocess
import traceback
import numpy as np
import pandas as pd

if __name__ == '__main__':
    try:
        project = sys.argv[1]
        coverage_threshold = float(sys.argv[3])
        evalue_threshold = float(sys.argv[4])

        # construct input and output paths
        in_path = Path('projects') / project / 'hits.txt'
        df_path = Path('projects') / project / 'hits_df.csv'
        out_faa_path = Path('projects') / project / 'filtered_hits.faa'

        # logging to exit log
        exitlog_path = Path('projects') / project / 'exit_log.txt'

        # logging start to exit log
        with open(exitlog_path, 'a') as outfile:
            outfile.write('4 started\n')
        search_result = SearchIO.read(in_path, "hmmer3-text")  # this is currently xml output of blast. to be custom later

        df = pd.read_csv(df_path, index_col=0)

        out = open(out_faa_path, 'w')

        #  filter df by thresholds
        df['filtered'] = [df['query_coverage'][i] > coverage_threshold and 
                           df['evalue'][i] < evalue_threshold 
                           for i in range(len(df['query_coverage']))
                         ]

        # retrieve aligned sequences of filtered hits
        ids_to_retrieve = list(df[df.filtered == True].index.values)
        for hit in search_result:
            hsp_no = 0
            for hsp in hit:
                hsp_no += 1

                ID = hit.id + str(hsp_no)

                if ID in ids_to_retrieve:
                    sequence = str(hsp.hit.seq)
                    sequence = sequence.replace('-', '')  # remove gaps
                    out.write('>' + ID + '\n' + sequence + '\n')

        df.to_csv(df_path)

        out.close()
        
        filtering_log_path = Path('projects') / project / 'filtering_log.txt'
        with open(filtering_log_path, 'w') as filtering_log:
            filtering_log.write("Coverage threshold: " + str(coverage_threshold) + '\n')
            filtering_log.write("e-value threshold: " + str(evalue_threshold))

    except Exception as e:
        ecx_type = str(type(e))
        print(f"RAISE {ecx_type} EXCEPTION")

        with open(exitlog_path, 'a') as outfile:
            outfile.write('4 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 to the exit log
        with open(exitlog_path, 'a') as outfile:
            outfile.write('4 0\n')
