"""Get a list of all taxids in the database from the `annotation.csv` file

- out `taxids.txt`

looks like
```
398511
511145
456320
306537
...
```

TODO: I should implement `get_all_taxids.py` in the `mk_blast_db.py` script

"""

import pandas as pd

df = pd.read_csv('annotation.csv')

taxids = df['taxid'].unique()

out = open('../taxids.txt', 'w')
for i in taxids:
    out.write(f'{i}\n')

out.close()
