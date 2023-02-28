"""Filter BLAST hits by parameters: identity, overlap, e-value.

- Step 4 in the pipeline
"""

import sys, getopt
from pathlib import Path
from Bio import Seq, SeqIO
from Bio import SeqFeature
from Bio.Blast import NCBIXML
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
        ident_threshold = float(sys.argv[3])  # 0.15
        coverage_threshold = float(sys.argv[4])  # 115
        evalue_threshold = float(sys.argv[5])  # 0.05

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

    proteins_to_filter = []

    n = 0
    for blast_record in blast_records:
        query_length = blast_record.query_length
        for alignment in blast_record.alignments:
            hsp_no = 0
            for hsp in alignment.hsps:
                hsp_no += 1
                n += 1
                print(n)

                ID = alignment.accession + str(hsp_no)
                identity = hsp.identities / hsp.align_length
                overlap = hsp.align_length
                length = alignment.length
                evalue = hsp.expect
                sequence = hsp.sbjct.replace('-', '')

                if identity > IDENT_THRESHOLD and overlap > OVERLAP_THRESHOLD and evalue < EVALUE_THRESHOLD:
                    out.write('>' + ID + '\n' + sequence + '\n')

    df = df[df.identity > IDENT_THRESHOLD]
    df = df[df.overlap > OVERLAP_THRESHOLD]
    df = df[df.evalue < EVALUE_THRESHOLD]
    df = df.drop(columns=['query', 'evalue', 'overlap', 'identity'])
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
