"""Make the dataframe of BLAST hits with annotations.

- input:
    - `blastp_hits.xml` - BLASTP search output (from `project` specified in arguments)
    - `annotation.csv` - .csv with annotations for each potential hit in the `database` specified in arguments

- output `blastp_hits_data.csv` - dataframe containing hit ids, BLAST output properties (identity, query coverage,
hit length), and additional annotations from the database

TODO: hit length can be precomputed and be a part of the database
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
    The lineage includes taxonomic ranks according to the `filter_ranks` variable

    Example:
        input: 511145
        output: ['Bacteria', 'Pseudomonadota','Gammaproteobacteria','Enterobacterales','Enterobacteriaceae',
        'Escherichia','Escherichia coli']
    """
    # taxonomic ranks to keep in the lineage
    filter_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    # get the lineage as a list of taxids
    lineage_taxid = ncbi.get_lineage(int(taxid))
    # make the lineage as dictionary {taxid: taxon name}
    lineage_names = ncbi.get_taxid_translator(lineage_taxid)
    # make a dictionary {taxid: taxon rank}
    lineage_ranks = ncbi.get_rank(lineage_taxid)
    # dict {rank name: taxon name}
    lineage_dict = {lineage_ranks[i]: lineage_names[i] for i in lineage_names}

    # filter lineage_dict by filter_ranks
    lineage_names = []
    for i in filter_ranks:
        taxon = lineage_dict.get(i)
        lineage_names.append(taxon)  # write organism
    return lineage_names


def get_replicon(replicon):
    """Get replicon type and name by replicon annotation of the source feature in a Genbank file.
    e.g.
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

        # NCBI taxonomy database
        ncbi = NCBITaxa()

        # log start to exit log
        exitlog_path = Path('projects') / project / 'exit_log.txt'
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '2 started'], stdout=outfile)  # TODO: make it through write

        # read input dataframe
        data_path = Path('databases') / database / 'annotation.csv'  # database will be stored in analysis configs
        prot_df = pd.read_csv(data_path, index_col=0)
        prot_df['protID'] = prot_df.index

        # create output dataframe
        out_df = pd.DataFrame(
            columns=['ID', 'protID', 'evalue', 'query_coverage', 'identity', 'length', 'targ_dom_pos'])

        # open homology search output to be parsed
        search_result_path = Path('projects') / project / 'blastp.xml'
        search_result = SearchIO.read(search_result_path, "blast-xml")  # it is xml output of blast. to be custom later

        # parse homology search output (iterate hits and their HSPs) and make dataframe of hits
        for hit in search_result:
            hsp_no = 0
            for hsp in hit:
                hsp_no += 1

                # row to be added to the dataframe
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

        # add annotations from `annotation.csv` to the hit dataframe
        out_df = pd.merge(out_df, prot_df, how='left', on='protID')

        # add columns with taxa of different taxonomic levels based on taxid
        out_df['superkingdom'], out_df['phylum'], out_df['class'], out_df['order'], \
            out_df['family'], out_df['genus'], out_df['species'] = zip(*map(get_lineage, out_df['taxid']))

        # add columns with replicon type and name based on the value in the dataframe
        out_df['replicon_type'], out_df['replicon_name'] = zip(*map(get_replicon, out_df['replicon']))
        out_df = out_df.drop('replicon', axis=1)

        # export output dataframe to csv
        out_df_path = Path('projects') / project / 'blastp_df.csv'
        out_df.to_csv(out_df_path, index=False)

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
