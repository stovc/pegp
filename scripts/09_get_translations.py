"""Retrieves full translations of `clustered_gc.faa` sequences from `filtered_clustered.csv`.

- `clustered-gc.faa` contains truncated translations as they result from a homology search

- Full translations are needed for the domain search,
while the truncated ones are needed for the alignment and tree construction

- Step 7a in the pipeline
"""

import sys
from pathlib import Path
import pandas as pd

# arguments - project name and database
project = sys.argv[1]
database = sys.argv[2]

# construct input and output paths
in_path = Path('projects') / project / 'filtered_clustered.csv'
data_path = Path('databases') / database / 'translation.csv'

out_path = Path('projects') / project / 'clustered_full.faa'

translation_db = pd.read_csv(data_path, index_col=0, names=['translation'])
clustered_db = pd.read_csv(in_path, index_col=0)

ids = list(clustered_db.index)

with open(out_path, 'w') as f:
    for i in ids:
        f.write('>' + i + '\n')
        f.write(translation_db.loc[i, 'translation'] + '\n')
