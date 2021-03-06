"""Construct a .csv containing cluster representative information for non-clustered proteins.

- input: cdhit 'clustered90.faa.clstr' output file
- output: cluster_dict.csv

```csv
id,representative
[unclustered protein id],[cluster representative protein id]
AB00047H1T,AB00047H1T
```

- Step 6a in the pipeline
"""

import sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO

# script argument - project
project = sys.argv[1]

# construct input and output paths
in_clstr_path = Path('projects') / project / 'clustered90.faa.clstr'

out_df_path = Path('projects') / project / 'cluster_dict.csv'

# open data
clstr = open(in_clstr_path, 'r')
out = open(out_df_path, 'w')

representative = ''
species_set = ''
species = []

d = {}

clstr.readline()

out.write('id,representative\n')

for line in clstr:
    print(line)
    if '*' in line:
        representative = line.split(' ')[-2][1:-3]       # e.g., '5\t191aa, >AB000MSVVJ... at 99.48%' -> 'AB000MSVVJ'
    if '%' in line:
        spec = line.split(' ')[-3][1:-3]               # e.g., '3\t191aa, >AB000MT09Q... at 99.48%' -> AB000MT09Q
        species.append(spec)
    if 'Cluster' in line:
        out.write(f'{representative},{representative}\n')
        for sp in species:
            out.write(f'{sp},{representative}\n')
        species = []

out.close()
clstr.close()
