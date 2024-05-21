"""Make a homology search database with metadata.

Generate species trees as .nwk and corresponding annotation data as .csv from `taxid.txt` related to the database.
The trees differ by the taxonomic rank to which they are collapsed

Input
-----
1) Path to a folder with genomes
2) Path to  metadata file
3) Path to the output database file

Output
------
1) Homology search database with metadata:
- Protein fasta file with all CDSs extracted from the genomes
- CSV file with all metadata for the CDSs and other features
- CSV file with upstream nucleotide sequences, nucleotide sequences, downstream nucleotide sequences
- CSV file with translations

2) Phylogenetic trees of the genomes collapsed to various taxonomic ranks.
[full, genus, family, order, class, phylum]."""

from ete3 import NCBITaxa, PhyloNode

from argparse import ArgumentParser, Namespace
import logging
import os
import pandas as pd
from pathlib import Path
import sys

import constants as c
from genome_parsing import GenomeParserConfig, parse_genome
from phylogeny import export_tree, export_tree_annotation, prune_tree
from utils import concat_files_from_folder, delete_folder, make_id_iterator


def create_logger(log_path) -> logging.Logger:
    """Create a logger."""
    logger = logging.getLogger(__name__)

    # Set the logging level of the logger to the lowest level needed
    logger.setLevel(logging.DEBUG)

    # Handler for stdout (logs INFO and above)
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)

    # Handler for file (logs DEBUG and above)
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)

    # Create a formatter and set it for both handlers
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    stdout_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # Add the handlers to the logger
    logger.addHandler(stdout_handler)
    logger.addHandler(file_handler)

    return logger


def make_all_paths_list(genomes_path):
    file_names = os.listdir(genomes_path)
    genome_file_names = [i for i in file_names if c.GENOME_EXTENSION in i]
    all_genome_paths = [genomes_path / i for i in genome_file_names]

    return all_genome_paths


def make_completed_paths_list(database_path):
    """Make a list paths of genomes completed in previous runs."""
    with open(database_path / 'completed_paths.txt', 'r') as f:
        completed_paths = f.readlines()
    completed_paths = [i.strip() for i in completed_paths]
    completed_paths = [Path(i) for i in completed_paths]

    return completed_paths


def make_uncompleted_genome_path_list(all_genome_paths, completed_genome_paths):
    genome_paths = [i for i in all_genome_paths if i not in completed_genome_paths]
    genome_paths.sort()

    return genome_paths


def parse_arguments():
    # expected arguments: 1) path to genomes; 2) path to metadata; 3) path to the database

    parser = ArgumentParser(prog='PegpDatabaseBuilding',
                            description='Make a homology search database with metadata',
                            epilog='LET\'S SEE HOW IT WORKS'
                            )

    # parser.add_argument('-h', '--help', help='Show this help message', action='help')
    parser.add_argument('genomes',
                        help='Path to the directory containing genomes')

    parser.add_argument('metadata',
                        help='Path to the file containing metadata')

    parser.add_argument('database',
                        help='Path to the directory with the output database. If it does not exist, it will be created')

    parser.add_argument('--no-context',
                        action='store_true', default=False,
                        help='Do not compute genome context of the features',)

    parser.add_argument('--no-ncbi-update',
                        action='store_true', default=False,
                        help='Do not update NCBI Taxonomy database')

    args: Namespace = parser.parse_args()

    return args


if __name__ == '__main__':
    # parse arguments
    args = parse_arguments()

    genomes_path = Path(args.genomes)
    database_path = Path(args.database)
    enable_genome_context = not args.no_context

    # make temporary folders and an assembly progress file
    if not os.path.exists(database_path):
        os.makedirs(database_path)
        for i in ['protein', 'sequence', 'translation', 'annotation']:
            os.makedirs(database_path / i)
        open(database_path / 'completed_paths.txt', 'a').close()

    # create logger
    log_path = database_path / 'database_assembly.log'
    logger = create_logger(log_path)
    logger.info(f'Database assembly. database: {args.database}, genomes: {args.genomes}')

    # generate a list of all genome paths to be processed
    all_genome_paths = make_all_paths_list(genomes_path=genomes_path)
    completed_genome_paths = make_completed_paths_list(database_path=database_path)
    genome_paths = make_uncompleted_genome_path_list(all_genome_paths, completed_genome_paths)

    completed_paths = open(database_path / 'completed_paths.txt', 'a')

    logger.info(f'total genomes {len(all_genome_paths)}, completed genomes {len(completed_genome_paths)}\n')

    # configure genome parser
    id_iterator = make_id_iterator(prefix=c.ID_PREFIX, symbols=c.ID_SUFFIX_SYMBOLS, suffix_length=c.ID_SUFFIX_LENGTH)
    genome_metadata_df = pd.read_csv(args.metadata, index_col=0, sep='\t', usecols=c.GENOME_METADATA_COLUMNS_OF_INTEREST)
    genome_parser_config = GenomeParserConfig(database_path, genome_metadata_df, id_iterator, enable_genome_context)

    # iterate genomes
    iteration = 1
    for genome_path in genome_paths:
        logger.info(f'{iteration} {genome_path}')
        parse_genome(genome_path, config=genome_parser_config)
        completed_paths.write(f'{genome_path}\n')
        iteration += 1

    concat_files_from_folder(folder=(database_path / 'protein'), extension='faa')

    FOLDERS_TO_CONCATENATE_CSV = ['annotation', 'sequence', 'translation']
    for folder in FOLDERS_TO_CONCATENATE_CSV:
        concat_files_from_folder(folder=(database_path / folder), extension='csv')

    # REMOVE FOLDERS
    # for folder in FOLDERS_TO_CONCATENATE_CSV + ['protein']:
    #     delete_folder(database_path / folder)

    # MAKE ORG TREE
    # GET TAXIDS
    data_path = database_path / 'annotation.csv'
    df = pd.read_csv(data_path)

    # df = df[df.gtdb_taxonomy.notnull()]

    taxids = df['taxid'].unique()

    # NCBI taxonomy database object
    ncbi = NCBITaxa()

    if not args.no_ncbi_update:
        ncbi.update_taxonomy_database()

    # get tree topology as PhyloNode object
    tree = ncbi.get_topology(taxids, intermediate_nodes=True)

    tree_full = prune_tree(tree, ['leaf', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'], database=ncbi)
    tree_genus = prune_tree(tree, ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'], database=ncbi)
    tree_family = prune_tree(tree, ['family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'], database=ncbi)
    tree_order = prune_tree(tree, ['order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'], database=ncbi)
    tree_class = prune_tree(tree, ['class', 'phylum', 'superkingdom', 'kingdom', 'root'], database=ncbi)
    tree_phylum = prune_tree(tree, ['phylum', 'superkingdom', 'kingdom', 'root'], database=ncbi)

    # EXPORT TREES AND ANNOTATIONS
    # construct path for the output folder and make it
    out_folder = database_path / 'org_trees'
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # export trees pruned to different taxonomic levels
    export_tree(tree_full, out_folder / 'org_tree_full.nwk')
    export_tree(tree_genus, out_folder / 'org_tree_genus.nwk')
    export_tree(tree_family, out_folder / 'org_tree_family.nwk')
    export_tree(tree_order, out_folder / 'org_tree_order.nwk')
    export_tree(tree_class, out_folder / 'org_tree_class.nwk')
    export_tree(tree_phylum, out_folder / 'org_tree_phylum.nwk')

    # export .csv annotations for the trees pruned to different taxonomic levels
    export_tree_annotation(tree_full, out_folder / 'org_tree_full_data.csv', database=ncbi)
