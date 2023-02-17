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

TODO: The job of `get_taxids.py` should be in the `mk_db.py` script

"""

import sys
from pathlib import Path
import pandas as pd

database = sys.argv[1]

database_path = Path('../databases/') / database
data_path = database_path / 'annotation.csv'
output_path = database_path / 'taxids.txt'

df = pd.read_csv(data_path)
print(df)

taxids = df['taxid'].unique()
print(type(taxids))
print(taxids)

out = open(output_path, 'w')
for i in taxids:
    out.write(f'{i}\n')

out.close()
