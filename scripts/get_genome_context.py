"""Get a dataframe describing the genome context of the hits.

- input: `filtered_clustered.csv`, `annotation.csv`
- output:  `genome_context.csv`

- Step 7b in the pipeline
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# arguments - project name and database
project = sys.argv[1]
database = sys.argv[2]

# construct input and output paths
in_path = Path('projects') / project / 'filtered_clustered.csv'
data_path = Path('databases') / database / 'annotation.csv'
out_path = Path('projects') / project / 'genome_context.csv'

annotation_db = pd.read_csv(data_path, index_col=0)
clustered_db = pd.read_csv(in_path)

ids = list(clustered_db.index)

with open(out_path, 'w') as f:
    f.write('molecule;gene;target;start;end;strand\n')
    n = 0
    for row in clustered_db.itertuples():
        n += 1
        print(n)

        feature_id = row.ID

        gene = str(row.gene)
        if gene == 'nan':
            gene = 'p: ' + str(row.product)

        target = 'T'
        start = row.start
        end = row.end
        sm = row.strand  # strand modifier

        if sm == 1:
            to_write = f'{feature_id};{gene};{target};{sm*start};{sm*end};1\n'
        elif sm == -1:
            to_write = f'{feature_id};{gene};{target};{sm * end};{sm * start};1\n'

        f.write(to_write)
        print(to_write)

        # writing dummy for alignment
        to_write = f'{feature_id};.;F;{sm*(start+(end-start)/2)};{sm*(start+(end-start)/2)};1\n'
        f.write(to_write)

        context = row.context
        context = context.split(';')

        for feature_id in context:
            rec = annotation_db.loc[feature_id, ]

            # assign gene
            gene = str(rec['gene'])
            if gene == 'nan':
                gene = str(rec['product'])
            if gene == 'nan':
                gene = str(rec['feature'])
            if gene == 'nan':
                gene = 'FTR'

            target = 'F'
            start = int(rec['start']) * sm
            end = int(rec['end']) * sm
            strand = int(rec['strand']) * sm


            if sm == 1:
                to_write = f'{row.ID};{gene};{target};{start};{end};{strand}\n'
            elif sm == -1:
                to_write = f'{row.ID};{gene};{target};{end};{start};{strand}\n'

            f.write(to_write)
            print(to_write)
