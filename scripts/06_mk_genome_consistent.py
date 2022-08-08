"""Ensure that all filtered hits that are present with clustered hits in the same genome are kept for the analysis.

- `clustered.faa` file is updated with sequences from `filtered_fits.faa`
if they are present in the same genome according to `filtered_hits_data.csv`

- Produces `clustered-gc.faa` file

- Step 6 in the pipeline
"""

import sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import subprocess

if __name__ == '__main__':
    # script argument - project
    project = sys.argv[1]

    # construct input and output paths
    in_faa_path = Path('projects') / project / 'clustered90.faa'
    in_unclustered_faa_path = Path('projects') / project / 'filtered_hits.faa'
    in_df_path = Path('projects') / project / 'filtered_hits_data.csv'

    out_faa_path = Path('projects') / project / 'clustered90-gc.faa'
    out_df_path = Path('projects') / project / 'filtered_clustered.csv'

    # logging to exit log
    exitlog_path = Path('projects') / project / 'exit_log.txt'

    with open(exitlog_path, 'a') as outfile:
        subprocess.run(["echo", '6 started'], stdout=outfile)

    df = pd.read_csv(in_df_path, index_col=0)

    clustered_ids = []
    for seq_record in SeqIO.parse(in_faa_path, "fasta"):
        clustered_ids.append(seq_record.id)

    # subset by list of IDs -> subset the column of genomes -> drop duplicates -> to list
    genomes = df[df.index.isin(clustered_ids)]['assembly'].drop_duplicates().to_list()
    retained_ids = list(df[df['assembly'].isin(genomes)].index)
    df = df[df.index.isin(retained_ids)]

    with open(out_faa_path, 'w') as f:
        for seq_record in SeqIO.parse(in_unclustered_faa_path, "fasta"):
            if seq_record.id in retained_ids:
                f.write('>' + seq_record.id + '\n')
                f.write(str(seq_record.seq) + '\n')

    df.to_csv(out_df_path)
