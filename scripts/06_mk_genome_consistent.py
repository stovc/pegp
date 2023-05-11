"""Ensure that all filtered hits that are present with clustered hits in the same genome are kept for the analysis.

- `clustered.faa` file is updated with sequences from `filtered_hits.faa`
if they are present in the same genome according to `filtered_hits.csv`

Input:
    - clustered.faa
    - filtered_hits.faa
    - hits_df.csv

Output:
    - clustered-gc.faa file

Add 'filtered_clustered' column into hits_df.csv. Value is True, if hit passed filtration 
and genome consistency operation
"""

import sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import subprocess
import traceback

if __name__ == '__main__':
    try:
        # script argument - project
        project = sys.argv[1]

        # construct input and output paths
        in_faa_path = Path('projects') / project / 'clustered90.faa'
        in_unclustered_faa_path = Path('projects') / project / 'filtered_hits.faa'
        df_path = Path('projects') / project / 'hits_df.csv'

        out_faa_path = Path('projects') / project / 'clustered90-gc.faa'

        # logging to exit log
        exitlog_path = Path('projects') / project / 'exit_log.txt'

        with open(exitlog_path, 'a') as outfile:
            outfile.write('6 started\n')

        df = pd.read_csv(df_path, index_col=0)

        clustered_ids = []
        for seq_record in SeqIO.parse(in_faa_path, "fasta"):
            clustered_ids.append(seq_record.id)

        # subset by list of IDs -> subset the column of genomes -> drop duplicates -> to list
        filtered_df = df[df.filtered == True]
        genomes = filtered_df[filtered_df.index.isin(clustered_ids)]['assembly'].drop_duplicates().to_list()
        retained_ids = list(df[df['assembly'].isin(genomes)].index)
        df['filtered_clustered'] = [df['filtered'][i] and df.index[i] in retained_ids
                                    for i in range(len(df['filtered']))
                                    ]

        with open(out_faa_path, 'w') as f:
            for seq_record in SeqIO.parse(in_unclustered_faa_path, "fasta"):
                if seq_record.id in retained_ids:
                    f.write('>' + seq_record.id + '\n')
                    f.write(str(seq_record.seq) + '\n')

        df.to_csv(df_path)
    except Exception as e:
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            outfile.write('6 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 to the exit log
        with open(exitlog_path, 'a') as outfile:
            outfile.write('6 0\n')
