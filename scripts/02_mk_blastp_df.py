"""Make the dataframe of BLAST hits with annotations.

- input `blastp_hits.xml`
- output `blastp_hits_data.csv`

- Step 2 in the pipeline
"""

import subprocess
import sys
from pathlib import Path
import traceback
from Bio import SearchIO
import pandas as pd
from ete3 import NCBITaxa


def get_lineage(taxid):
    """Get NCBI taxid of a species
    Return the lineage of the species as a list of taxon names
    The lineage includes taxonomic ranks according to the `ranks` variable
    """

    lineage_taxid = ncbi.get_lineage(int(taxid))  # get the lineage as a list of taxids
    names = ncbi.get_taxid_translator(lineage_taxid)  # make a dictionary {taxid: taxon name}
    ranks = ncbi.get_rank(lineage_taxid)  # make a dictionary {taxid: taxon rank}

    lineage_dict = {ranks[i]: names[i] for i in lineage}  # dict {rank name: taxon name}

    lineage_names = []
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for i in ranks:
        taxon = lineage_dict.get(i)
        lineage_names.append(taxon)  # write organism
    return lineage_names


def get_replicon(replicon):
    """e.g.
    in - 'chromosome_II'
    out - ('chromosome', 'II')

    but 'main' --> ('main', 'main')

    TODO: this should be done in database_building/mk_db.py
    """

    if replicon == 'main':
        replicon_type = 'main'
        replicon_name = 'main'
    else:
        print(replicon)
        replicon_type = replicon.split('_')[0]
        replicon_name = ''.join(replicon.split('_')[1:])
    return replicon_type, replicon_name


if __name__ == '__main__':
    try:
        # args
        project = sys.argv[1]
        database = sys.argv[2]

        # NCBI taxonomy
        ncbi = NCBITaxa()

        # construct paths
        in_path = Path('projects') / project / 'blastp.xml'
        data_path = Path('databases') / database / 'annotation.csv'  # database will be stored in analysis configs
        out_path = Path('projects') / project / 'blastp_df.csv'
        exitlog_path = Path('projects') / project / 'exit_log.txt'

        # logging start to exit log
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '2 started'], stdout=outfile)  # TODO: make it through write

    # open input and output
    result_handle = open(in_path, 'r')
    prot_df = pd.read_csv(data_path, index_col=0)
    # print(prot_df)
    # print(prot_df.index)
    # print(prot_df.index.name)
    # print(prot_df.index.values)
    prot_df['protID'] = prot_df.index

        out_df = pd.DataFrame(
            columns=['ID', 'protID', 'evalue', 'query_coverage', 'identity', 'length', 'targ_dom_pos'])

        # parse homology search output and make dataframe of hits
        search_result = SearchIO.read(in_path, "blast-xml")  # this is currently xml output of blast. to be custom later
        n = 0
        for hit in search_result:
            hsp_no = 0
            for hsp in hit:
                hsp_no += 1
                n += 1
                print(n)

                # row to be added to the dataframe
                print('query span', hsp.query_span)
                print('full query', search_result.seq_len)

                new_row = pd.DataFrame(
                    {
                        'ID': hit.id + str(hsp_no),  # hit id = protein_id + hsp number
                        'protID': hit.id,  # id of the protein the hsp belongs to
                        'evalue': hsp.evalue,
                        'query_coverage': 100 * hsp.query_span / search_result.seq_len,
                        'identity': hsp.ident_num / hsp.aln_span,
                        'length': hit.seq_len,
                        'targ_dom_pos': (hsp.hit_start + hsp.query_span // 2)  # coordinate of the center of the query
                    },
                    index=[0]
                )
                out_df = pd.concat([out_df, new_row])

        print(out_df)

        # add annotations from `annotation.csv` to the hit dataframe
        out_df = pd.merge(out_df, prot_df, how='left', on='protID')

        print(out_df)

        # add columns with belonging to taxa on different taxonomic levels based on taxid
        out_df['superkingdom'], out_df['phylum'], out_df['class'], out_df['order'], \
            out_df['family'], out_df['genus'], out_df['species'] = zip(*map(get_lineage, out_df['taxid']))

        # add columns with replicon type and name based on the value in the dataframe
        out_df['replicon_type'], out_df['replicon_name'] = zip(*map(get_replicon, out_df['replicon']))
        out_df = out_df.drop('replicon', axis=1)

        out_df.to_csv(out_path, index=False)

    except Exception as e:
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '2 ' + ecx_type], stdout=outfile)

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 nto the exit log
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '2 0'], stdout=outfile)
