"""Make the dataframe of BLAST hits with annotations.

- input:
    - `hits.txt` - hmm search output (from `project` specified in arguments)
    - `annotation.csv` - .csv with annotations for each potential hit in the `database` specified in arguments

- output `blastp_hits_data.csv` - dataframe containing hit ids, hmmer output properties (identity, query coverage,
hit length), and additional annotations from the database

TODO: hit length can be precomputed and be a part of the database

FIXME: not blast anymore
"""

import subprocess
import sys
import math
from pathlib import Path
import traceback
from Bio import SearchIO
import pandas as pd


def get_lineage(gtdb_taxonomy):
    """Parse GTDB taxonomy string and return lineage names.

    Example:
        input:
            d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhizobiales;
            f__Xanthobacteraceae;g__Azorhizobium;s__Azorhizobium caulinodans

        output:
            [
                "Bacteria",
                "Pseudomonadota",
                "Alphaproteobacteria",
                "Rhizobiales",
                "Xanthobacteraceae",
                "Azorhizobium",
                "Azorhizobium caulinodans",
            ]
    """

    rank_prefixes = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }

    output_ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    lineage = {rank: None for rank in output_ranks}

    if pd.isna(gtdb_taxonomy):
        return [lineage[rank] for rank in output_ranks]

    for taxon in str(gtdb_taxonomy).split(";"):
        if "__" not in taxon:
            continue

        prefix, name = taxon.split("__", 1)
        rank = rank_prefixes.get(prefix)

        if rank is not None:
            lineage[rank] = name if name else None

    return [lineage[rank] for rank in output_ranks]


if __name__ == '__main__':
    try:
        # args
        project = sys.argv[1]
        database = sys.argv[2]

        # log start to exit log
        exitlog_path = Path('projects') / project / 'exit_log.txt'
        with open(exitlog_path, 'a') as outfile:
            outfile.write('2 started\n')

        # read input dataframe
        data_path = Path('databases') / database / 'annotation.csv'  # database will be stored in analysis configs
        prot_df = pd.read_csv(data_path, index_col=0)
        # prot_df['lcs'] = prot_df.index

        # create output dataframe
        out_df = pd.DataFrame(
            columns=['hsp', 'lcs', 'evalue', 'lg_evalue', 'query_coverage', 'targ_dom_pos'])

        # open homology search output to be parsed
        search_result_path = Path('projects') / project / 'hits.txt'
        search_result = SearchIO.read(search_result_path, "hmmer3-text")

        # parse homology search output (iterate hits and their HSPs) and make dataframe of hits
        for hit in search_result:
            hsp_no = 0
            for hsp in hit:
                hsp_no += 1

                # row to be added to the dataframe
                new_row = pd.DataFrame(
                    {
                        'hsp': hit.id + str(hsp_no),  # hit id = protein_id + hsp number
                        'lcs': hit.id,  # id of the protein the hsp belongs to
                        'evalue': hsp.evalue,
                        'lg_evalue': math.log10(hsp.evalue),
                        'query_coverage': 100 * hsp.query_span / search_result.seq_len,
                        'targ_dom_pos': (hsp.hit_start + hsp.query_span // 2)  # coordinate of the center of the query
                    },
                    index=[0]
                )
                out_df = pd.concat([out_df, new_row])

        # add annotations from `annotation.csv` to the hit dataframe
        out_df = pd.merge(out_df, prot_df, how='left', on='lcs')

        # add columns with taxa of different taxonomic levels based on GTDB taxonomy
        out_df['domain'], out_df['phylum'], out_df['class'], out_df['order'], \
            out_df['family'], out_df['genus'], out_df['species'] = zip(
            *map(get_lineage, out_df['gtdb_taxonomy'])
        )

        # export output dataframe to csv
        out_df_path = Path('projects') / project / 'hits_df.csv'
        out_df.to_csv(out_df_path, index=False)

    except Exception as e:
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            outfile.write('2 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 nto the exit log
        with open(exitlog_path, 'a') as outfile:
            outfile.write('2 0\n')
