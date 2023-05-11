"""Make a .csv file with domain data from hmmscan tabular output.

Input:
    - "hmmscan.tbl" - hmmscan output
    - "hits_df.csv" - annotations
Output:
    - "domains.csv" - predicted domains with annotations
"""

import sys
import subprocess
import traceback
from pathlib import Path
import pandas as pd

if __name__ == '__main__':
    try:
        # argument - project name
        project = sys.argv[1]

        # logging to exit log
        exitlog_path = Path('projects') / project / 'exit_log.txt'

        with open(exitlog_path, 'a') as outfile:
            outfile.write('15 started\n')

        # construct input and output paths


        out_path = Path('projects') / project / 'domains.csv'

        names = ['domain', 'id', 'evalue', 'start', 'stop']

        # hmmscan data
        in_hmmscan_path = Path('projects') / project / 'hmmscan.tbl'
        hmmscan_df = pd.read_csv(in_hmmscan_path,
                                 delim_whitespace=True,
                                 skiprows=3, skipfooter=10,
                                 usecols=[0, 3, 6, 17, 18],
                                 index_col=1,
                                 names=names,
                                 engine='python')

        grouped = hmmscan_df.groupby('id')  # each group consists of all domains found in particular protein

        # protein data
        in_data_path = Path('projects') / project / 'hits_df.csv'
        prot_data = pd.read_csv(in_data_path, usecols=['ID', 'length', 'targ_dom_pos', 'filtered_clustered'], index_col=0)
        prot_data = prot_data[prot_data.filtered_clustered == True]
        ids = prot_data.index
        target_domain_positions = prot_data['targ_dom_pos']
        protein_lengths = prot_data['length']

        out_df = pd.DataFrame()

        for name, group in grouped:
            for row in group.itertuples():
                mask = (row.evalue < group.evalue) & (row.start < group.stop) & (row.stop > group.start)
                overlaps = bool(sum(mask))
                if not overlaps:
                    to_append = pd.DataFrame(row).T
                    out_df = pd.concat([out_df, to_append], ignore_index=True)

        # add dummy annotations to align domain plot
        for i in ids:
            # append target domain positions
            position = target_domain_positions[i]
            to_append = {0: i, 1: '.', 2: 0, 3: position, 4: position + 1}   # apparently I can remove this column 2 --- test!
            out_df = out_df.append(to_append, ignore_index=True)  # TODO: change append to concat
            # append whole proteins
            protein_length = protein_lengths[i]
            to_append = {0: i, 1: '_', 2: 0, 3: 0, 4: protein_length}
            out_df = out_df.append(to_append, ignore_index=True)  # TODO: change append to concat

        out_df = out_df.drop(columns=[2])
        out_df = out_df.rename(columns={0: 'id', 1: 'domain', 3: 'start', 4: 'end'})
        out_df = out_df.set_index('id')
        out_df.to_csv(out_path)
    except Exception as e:
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            outfile.write('15 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 to the exit log
        with open(exitlog_path, 'a') as outfile:
            outfile.write('15 0\n')
