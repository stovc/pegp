"""Retrieves full translations of `clustered_gc.faa` sequences from `filtered_clustered.csv`.

- `clustered-gc.faa` contains truncated translations as they result from a homology search

- Full translations are needed for the domain search,
while the truncated ones are needed for the alignment and tree construction

- Step 7a in the pipeline
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
            subprocess.run(["echo", '9 started'], stdout=outfile)

        # construct input and output paths
        in_path = Path('projects') / project / 'filtered_clustered.csv'
        data_path = Path('databases') / database / 'translation.csv'

        out_path = Path('projects') / project / 'clustered_full.faa'

        translation_db = pd.read_csv(data_path, index_col=0, names=['translation'])
        clustered_db = pd.read_csv(in_path, index_col=0)  # TODO: to get index, there's no need to read all columns

        ids = list(clustered_db.index)

        with open(out_path, 'w') as f:
            for id in ids:
                protID = id[:-1]
                f.write('>' + protID + '\n')
                f.write(translation_db.loc[protID, 'translation'] + '\n')
    except Exception as e:
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '9 ' + ecx_type], stdout=outfile)

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 to the exit log
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '9 0'], stdout=outfile)
