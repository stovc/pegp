"""Filter BLAST hits by parameters: identity, overlap, e-value.

- Step 4 in the pipeline
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
        print('args')
        for i in sys.argv:
            print(i)
        project = sys.argv[1]
        ident_threshold = float(sys.argv[3])
        coverage_threshold = float(sys.argv[4])
        evalue_threshold = float(sys.argv[5])

        # construct input and output paths
        in_path = Path('projects') / project / 'blastp.xml'
        df_path = Path('projects') / project / 'blastp_df.csv'
        out_df_path = Path('projects') / project / 'filtered_hits.csv'
        out_faa_path = Path('projects') / project / 'filtered_hits.faa'

        # logging to exit log
        exitlog_path = Path('projects') / project / 'exit_log.txt'

        # logging start to exit log
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '4 started'], stdout=outfile)
        search_result = SearchIO.read(in_path, "blast-xml")  # this is currently xml output of blast. to be custom later

        df = pd.read_csv(df_path, index_col=0)

        out = open(out_faa_path, 'w')

        #  filter df by thresholds
        df = df[df.identity > ident_threshold]
        df = df[df.query_coverage > coverage_threshold]
        df = df[df.evalue < evalue_threshold]

        # retrieve aligned sequences of filtered hits
        ids_to_retrieve = list(df.index.values)
        for hit in search_result:
            hsp_no = 0
            for hsp in hit:
                hsp_no += 1

                ID = hit.id + str(hsp_no)

                if ID in ids_to_retrieve:
                    sequence = str(hsp.hit.seq)
                    sequence = sequence.replace('-', '')  # remove gaps
                    out.write('>' + ID + '\n' + sequence + '\n')

        df.to_csv(out_df_path)

        out.close()
    except Exception as e:
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '4 ' + ecx_type], stdout=outfile)

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 to the exit log
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '4 0'], stdout=outfile)
