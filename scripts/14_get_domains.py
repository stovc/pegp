"""Make a .csv file with domain data from hmmscan tabular output.

- Step 9a in the pipeline
"""

import sys
from pathlib import Path
import pandas as pd

# argument - project name
project = sys.argv[1]
#database = sys.argv[2]

# construct input and output paths
in_hmmscan_path = Path('projects') / project / 'hmmscan.tbl'
in_data_path = Path('projects') / project / 'filtered_clustered.csv'

out_path = Path('projects') / project / 'domains.csv'

#names = ['target', 'accession', 'tlen', 'query', 'qaccession', 'qlen', 'evalue', 'escore', 'ebias', '#',
#         'of', 'c-Evalue', 'i-Evalue', 'iscore', 'ibias', 'hmm from', 'hmm to', 'start', 'stop',
#         'env from', 'env to', 'acc', 'description']

names = ['domain', 'id', 'evalue', 'start', 'stop']

hmmscan_df = pd.read_csv(in_hmmscan_path,
                         delim_whitespace=True,
                         skiprows=3, skipfooter=10,
                         usecols=[0, 3, 6, 17, 18],
                         index_col=1,
                         names=names)
prot_series = pd.read_csv(in_data_path, usecols=[0,2], index_col=0)
prot_series = prot_series['targ_dom_pos']
print(prot_series)
print(type(prot_series))
grouped = hmmscan_df.groupby('id')

out_df = pd.DataFrame()

for name, group in grouped:
    for row in group.itertuples():
        print(type(row))
        print(row)
        mask = (row.evalue < group.evalue) & (row.start < group.stop) & (row.stop > group.start)
        print(mask)
        overlaps = bool(sum(mask))
        print(overlaps)
        if not overlaps:
            to_append = pd.DataFrame(row).T
            print(to_append)
            out_df = pd.concat([out_df, to_append], ignore_index=True)
            print(out_df)

ids = prot_series.index
for i in ids:
    position = prot_series[i]
    to_append = {0: i, 1: '.', 2: 0, 3: position, 4: position + 1}
    print(to_append)
    out_df = out_df.append(to_append, ignore_index=True)

out_df = out_df.drop(columns=[2])
out_df = out_df.rename(columns={0: 'id', 1: 'domain', 3: 'start', 4: 'end'})
print(out_df)
out_df = out_df.set_index('id')
out_df.to_csv(out_path)
