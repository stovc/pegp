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
- CSV file with upstream nucleotide sequences
- CSV file with nucleotide sequences
- CSV file with downstream nucleotide sequences
- CSV file with translations

2) Phylogenetic trees of the genomes collapsed to various taxonomic ranks.
[full, genus, family, order, class, phylum]."""

from Bio import SeqIO
from ete3 import NCBITaxa, PhyloNode

from argparse import ArgumentParser, Namespace
import itertools
import os
import pandas as pd
from pathlib import Path

from utils import circular_slice, get_first, concat_files_from_folder, delete_folder, distance_between_sequences
from phylogeny import export_tree, export_tree_annotation, prune_tree
import config as c


def add_genome_context(seq_record, seq_record_data):
    if seq_record.annotations['topology'] == 'circular':
        circular_length = len(seq_record.seq)
    else:
        circular_length = None

    context = ''

    short_data = seq_record_data[['lcs', 'start', 'end']]
    short_data = short_data.assign(context='')

    length = len(short_data.index)

    # i is the number of the iterated element the context is being inferred for
    for i in range(length):
        # iterate j elements to the right from the i-th element
        short = False
        within_window = True
        j = i + 1
        while within_window:
            # jump to start if reached the end
            if j >= length:
                if circular_length is not None:
                    j = 0
                else:
                    break

            # starts and ends of the i-th and j-th features
            start_i = short_data.iat[i, 1]
            end_i = short_data.iat[i, 2]
            start_j = short_data.iat[j, 1]
            end_j = short_data.iat[j, 2]

            if j == i:
                within_window = False
                short = True
            else:
                if distance_between_sequences(start_i, end_i, start_j, end_j, circular_length) <= c.CONTEXT_WINDOW_SIZE:
                    context += short_data.iat[j, 0] + ';'
                else:
                    within_window = False
            j += 1

        # iterate j elements to the left from the i-th element IF IT THE REPLICON ISN'T TOO SHORT
        if not short:
            outside_window = False
            j = i - 1
            while not outside_window:

                start1 = short_data.iat[i, 1]
                end1 = short_data.iat[i, 2]
                start2 = short_data.iat[j, 1]
                end2 = short_data.iat[j, 2]

                if distance_between_sequences(start1, end1, start2, end2, circular_length) <= c.CONTEXT_WINDOW_SIZE:
                    context += seq_record_data.iat[j, 0] + ';'
                else:
                    outside_window = True
                j -= 1

                if (j < 0) and (circular_length is None):
                    outside_window = True

        context = context[:-1]  # remove last comma
        short_data.iat[i, 3] = context
        context = ''

    seq_record_data_with_context = pd.merge(seq_record_data, short_data[['lcs', 'context']], on='lcs')

    return seq_record_data_with_context


def get_feature_type(feature):
    if feature.type == 'repeat_region':
        feature_type = get_first(feature.qualifiers, 'rpt_family')
    elif feature.type == 'regulatory':
        feature_type = get_first(feature.qualifiers, 'regulatory_class')
    elif feature.type == 'ncRNA':
        feature_type = get_first(feature.qualifiers, 'ncRNA_class')
    else:
        feature_type = feature.type

    return feature_type


def get_replicon_data_from_source(feature):
    plasmid = feature.qualifiers.get('plasmid')
    log.write(f'Plasmid: {plasmid}\n')
    segment = feature.qualifiers.get('segment')
    log.write(f'Segment: {segment}\n')
    chromosome = feature.qualifiers.get('chromosome')
    log.write(f'Chromosome: {chromosome}\n')

    if plasmid is not None:
        replicon_type = 'plasmid'
        replicon = plasmid[0]
    elif segment is not None:
        replicon_type = 'segment'
        replicon = segment[0]
    elif chromosome is not None:
        replicon_type = 'chromosome'
        replicon = chromosome[0]
    else:
        replicon_type = 'chromosome'
        replicon = 'main'

    return [replicon_type, replicon]


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

    parser = ArgumentParser()

    # parser.add_argument('-h', '--help', help='Show this help message', action='help')
    parser.add_argument('genomes', help='Path to the directory containing genomes')
    parser.add_argument('metadata', help='Path to the file containing metadata')

    parser.add_argument('database',
                        help='Path to the directory with the output database. If it does not exist, it will be created')

    parser.add_argument('--nocontext', help='Do not compute genome context of the features',
                        action='store_true', default=False)

    args: Namespace = parser.parse_args()

    return args


def parse_genome(genome_path):
    for seq_record in SeqIO.parse(genome_path, 'genbank'):

        seq_record_data = []
        accession = genome_path.name
        accession = accession.split('_')[0] + accession.split('_')[1]  # genome name, e.g. GCF_000701165.1
        log.write('Assembly: ' + accession + ' ')

        accession = accession[:3] + '_' + accession[3:]

        # open output files
        # fasta file with all protein sequences
        protein_output = open(database_path / 'protein' / f'{accession}', 'w')
        # other annotations discribed
        upstream_out = open(database_path / 'upstream' / f'{accession}', 'w')
        sequence_out = open(database_path / 'sequence' / f'{accession}', 'w')
        downstream_out = open(database_path / 'downstream' / f'{accession}', 'w')
        translation_out = open(database_path / 'translation' / f'{accession}', 'w')
        annotation_out = open(database_path / 'annotation' / f'{accession}', 'w')

        definition = seq_record.description  # Bifidobacterium breve strain NRBB18 chromosome, complete genome
        log.write('Definition: ' + definition + '\n')

        taxid = 'NONE'
        replicon = 'NONE'
        replicon_metadata = ['NONE', 'NONE']

        genome_metadata = [accession] + genome_metadata_df.loc[[accession]].values.flatten().tolist()

        for feature in seq_record.features:  # extract taxid and replicon from the source feature
            # "source" feature is single per seq_record/replicon and contains replicon related information
            # this feature is not written to the database; its properties are assigned to every database entry

            # EXTRACT REPLICON METADATA

            if feature.type == 'source':
                replicon_metadata = get_replicon_data_from_source(feature)

            # EXTRACT FEATURE METADATA

            elif feature.type not in ['source', 'gene']:
                # generate lcs id
                lcs = next(id_iterator)
                lcs = ''.join(lcs)
                lcs = id_prefix + lcs

                feature_type = get_feature_type(feature)

                # coordinates
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                coordinates = [start, end, strand]

                # sequence
                seq_record_seq = str(seq_record.seq)
                upstream = circular_slice(seq_record_seq, start - c.UTR_WINDOW_SIZE, start)
                sequence = seq_record_seq[start:end]
                downstream = circular_slice(seq_record_seq, end, end + c.UTR_WINDOW_SIZE)
                translation = get_first(feature.qualifiers, 'translation')

                # annotations
                gene = get_first(feature.qualifiers, 'gene')
                product = get_first(feature.qualifiers, 'product')

                # annotations i dont use
                # protein_id = get_first(feature.qualifiers, 'protein_id')

                # WRITE METADATA
                if translation is not None:
                    protein_output.write(f'>{lcs}\n{translation}\n')  # write a fasta record with translation
                    translation_out.write(f'{lcs},{translation}\n')
                    protein_length = len(translation)
                else:
                    protein_length = None

                upstream_out.write(f'{lcs},{upstream}\n')
                sequence_out.write(f'{lcs},{sequence}\n')
                downstream_out.write(f'{lcs},{downstream}\n')

                feature_metadata = [feature_type, gene, product, start, end, strand, protein_length]
                annotation = [lcs] + genome_metadata + replicon_metadata + feature_metadata
                seq_record_data.append(annotation)

        # make a dataframe with metadata
        columns = ['lcs', 'assembly'] + genome_metadata_df.columns.to_list() + ['replicon_type', 'replicon'] + \
                  ['feature_type', 'gene', 'product', 'start', 'end', 'strand', 'protein_length']
        seq_record_data = pd.DataFrame(seq_record_data, columns=columns)

        # INFER GENOME CONTEXT
        if not args.nocontext:
            seq_record_data = add_genome_context(seq_record, seq_record_data)

        annotation_out.write(seq_record_data.to_csv(index=False))

        protein_output.close()
        annotation_out.close()
        upstream_out.close()
        sequence_out.close()
        downstream_out.close()
        translation_out.close()

    completed_paths.write(str(genome_path) + '\n')


if __name__ == '__main__':
    args = parse_arguments()

    genomes_path = Path(args.genomes)
    database_path = Path(args.database)

    # make temporary folders and an assembly progress file
    if not os.path.exists(database_path):
        os.makedirs(database_path)

        for i in ['protein', 'upstream', 'sequence', 'downstream', 'translation', 'annotation']:
            os.makedirs(database_path / i)

        open(database_path / 'completed_paths.txt', 'a').close()

    log = open(database_path / 'database_assembly_log.txt', 'a')  # log
    print(f'Database assembly. database: {args.database}, genomes: {args.genomes}')
    log.write(f'Database assembly. db: {args.database}, folders: {args.genomes}\n')

    # generate a list of all genome paths to be processed
    all_genome_paths = make_all_paths_list(genomes_path=genomes_path)
    completed_genome_paths = make_completed_paths_list(database_path=database_path)

    # log msg
    print(f'total genomes {len(all_genome_paths)}, completed genomes {len(completed_genome_paths)}')
    log.write(f'total genomes {len(all_genome_paths)}, completed genomes {len(completed_genome_paths)}\n')

    genome_paths = make_uncompleted_genome_path_list(all_genome_paths, completed_genome_paths)
    completed_paths = open(database_path / 'completed_paths.txt', 'a')

    # generated IDs will consist of a prefix and an iterated part
    id_prefix = 'AB'
    id_iterator = itertools.product(c.HEADER_SYMBOLS, repeat=8)

    genome_metadata_df = pd.read_csv(args.metadata, index_col=0, sep='\t', usecols=c.GENOME_METADATA_COLUMNS_OF_INTEREST)

    iteration = 1
    # iterate genomes
    for genome_path in genome_paths:
        log.write(f'{iteration} {genome_path}\n')
        print(f'{iteration} {genome_path}')
        iteration += 1
        # iterate seq. records (replicons) in a genbank file

        parse_genome(genome_path)

    log.close()

    concat_files_from_folder(folder=(database_path / 'protein'), extension='faa')

    FOLDERS_TO_CONCATENATE_CSV = ['annotation', 'upstream', 'sequence', 'downstream', 'translation']
    for folder in FOLDERS_TO_CONCATENATE_CSV:
        concat_files_from_folder(folder=(database_path / folder), extension='csv')

    # REMOVE FOLDERS

    for folder in FOLDERS_TO_CONCATENATE_CSV + ['protein']:
        delete_folder(database_path / folder)

    # MAKE ORG TREE
    # GET TAXIDS
    data_path = database_path / 'annotation.csv'
    df = pd.read_csv(data_path)

    # df = df[df.gtdb_taxonomy.notnull()]

    taxids = df['taxid'].unique()

    # NCBI taxonomy database object
    ncbi = NCBITaxa()
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
