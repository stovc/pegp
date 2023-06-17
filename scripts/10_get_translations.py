"""Retrieves full translations of `clustered_gc.faa` sequences from `hits_df.csv`, which passed 
filtration, clusterization and genome concistency operation.

- Full translations are needed for the domain search,
while the truncated ones are needed for the alignment and tree construction

Input:
    - "hits_df.csv"
    - "translation.csv"
Output:
    - "clustered_full.faa"
"""

import sys
from pathlib import Path
import pandas as pd
import subprocess
import traceback

if __name__ == '__main__':
    try:
        # arguments - project name and database
        project = sys.argv[1]
        database = sys.argv[2]

        # logging to exit log
        exitlog_path = Path('projects') / project / 'exit_log.txt'

        with open(exitlog_path, 'a') as outfile:
            outfile.write('10 started\n')

        # construct input and output paths
        in_path = Path('projects') / project / 'hits_df.csv'
        data_path = Path('databases') / database / 'translation.csv'

        out_path = Path('projects') / project / 'clustered_full.faa'

        translation_db = pd.read_csv(data_path, index_col=0, names=['translation'])
        clustered_db = pd.read_csv(in_path, usecols=['ID', 'filtered_clustered'], index_col=0)
        clustered_db = clustered_db[clustered_db.filtered_clustered == True]

        ids = list(clustered_db.index)

        with open(out_path, 'w') as f:
            for id in ids:
                lcs = id[:-1]
                f.write('>' + id + '\n')  # get a separate record for each hsp to be compatible with tree-annotation
                f.write(translation_db.loc[lcs, 'translation'] + '\n')
    except Exception as e:
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            outfile.write('10 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 to the exit log
        with open(exitlog_path, 'a') as outfile:
            outfile.write('10 0\n')
