"""Make the dataframe of BLAST hits with annotations.

- input `blastp_hits.xml`
- output `blastp_hits_data.csv`

- Step 2 in the pipeline
"""

import sys, getopt
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from Bio.Blast import NCBIXML
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from ete3 import NCBITaxa


def get_lineage(taxid):
    lineage = ncbi.get_lineage(int(taxid))  # get the lineage as a list of taxids
    names = ncbi.get_taxid_translator(lineage)  # make a dictionary {taxid: taxon name}
    ranks = ncbi.get_rank(lineage)  # make a dictionary {taxid: taxon rank}

    lineage_dict = {ranks[i]: names[i] for i in lineage}  # dict {rank name: taxon name}

    out = []
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for i in ranks:
        taxon = lineage_dict.get(i)
        out.append(taxon)  # write organism
    return out


def get_replicon(replicon):
    if replicon == 'main':
        replicon_type = 'main'
        replicon_name = 'main'
    else:
        print(replicon)
        replicon_type = replicon.split('_')[0]
        replicon_name = ''.join(replicon.split('_')[1:])
    return replicon_type, replicon_name


project = sys.argv[1]
database = sys.argv[2]
data_path = Path('databases') / database / 'annotation.csv'  # will be stored in analysis configs

# NCBI taxonomy
ncbi = NCBITaxa()

# construct input and output paths
in_path = Path('projects') / project / 'blastp.xml'
out_path = Path('projects') / project / 'blastp_df.csv'

# open input and output
result_handle = open(in_path, 'r')
prot_df = pd.read_csv(data_path, index_col=0)

out_df = pd.DataFrame(columns=['ID', 'query', 'evalue', 'overlap', 'identity', 'length', 'targ_dom_pos'])

genomes = {}
found_hits = []

blast_records = NCBIXML.parse(result_handle)

n = 0
for blast_record in blast_records:
    query_length = blast_record.query_length
    query_id = blast_record.query_id.split('_')[1]  # Query_1 -> 1
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:

            n += 1
            print(n)

            ID = alignment.accession

            new_row = {'ID': ID, 'query': query_id,
                       'evalue': hsp.expect, 'overlap': hsp.align_length,
                       'identity': hsp.identities / hsp.align_length, 'length': alignment.length,
                       'targ_dom_pos': (hsp.sbjct_start + (len(hsp.sbjct) // 2))}
            out_df = out_df.append(new_row, ignore_index=True)

out_df = pd.merge(out_df, prot_df, how='left', on='ID')

out_df['superkingdom'], out_df['phylum'], out_df['class'], out_df['order'], \
out_df['family'], out_df['genus'], out_df['species'] = zip(*map(get_lineage, out_df['taxid']))

out_df['replicon_type'], out_df['replicon_name'] = zip(*map(get_replicon, out_df['replicon']))

out_df = out_df.drop('replicon', axis = 1)

out_df.to_csv(out_path, index=False)
