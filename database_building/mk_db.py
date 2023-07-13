"""Make homology search database input and annotation files.

Iterate all Genebank genomes (.gbk) in the input directory in order to
prepare
- an _input file for construction of the homology search database_
- _associated annotation files_

Input - [database_name] [path to the folder with genomes]

For each genome, iterate all *sequence features* and make the following files in `databases/[database_name]`:
- 'protein.faa' with translations of protein coding genes (*this file will be used to generate BLAST database in a further step*)

and other files containing interesting information about our proteins:
- 'annotation.csv' with the most useful _annotations_, described below

other less frequently used and more heavy-weighted annotations are written into separate .csv files
- 'upstream.csv' with 200 bp upstream nucleotide sequences
- 'downstream.csv' with 200 bp downstream nucleotide sequences
- 'sequence.csv' with the nucleotide sequence of the genes
- 'translation.csv' with translations protein coding genes

In all generated files, each record is preceded by an ID that is consistent across the files

_Annotations_ are in the annotation.csv file and include:
- Identity
    - ID (generated for the created database)
    - locus tag
    - assembly
    - accession
    - protein_id
- Position in the genome
    - replicon (identity of the chromosome or the plasmid)
    - start
    - stop
    - strand
- Taxonomy information
    - taxid (NCBI taxonomy ID)
- Gene annotation
    - feature type
    - gene
    - product

Generate species trees as .nwk and corresponding annotation data as .csv from `taxid.txt` related to the database.
The trees differ by the taxonomic rank to which they are collapsed

input - taxid.txt
output - org_tree_[TAXONOMIC RANK].nwk + 'org_tree_[TAXONOMIC RANK]_data.csv' where
the taxonomic ranks are [full, genus, family, order, class, phylum]

"""

import os, sys
from pathlib import Path
from argparse import ArgumentParser, Namespace
import string
import itertools
from Bio import Seq, SeqIO
from Bio import SeqFeature
import pandas as pd
from ete3 import NCBITaxa, PhyloNode


def distance(start1, end1, start2, end2, circular_length=None):
    """Return distance between two SeqFeatures. Supports circular nucleic acid
    Nucleic acid is supposed to be linear unless circular_length argument is passed.

    Arguments:
    start1, end1 -- start and end nucleotides of the 1st SeqFeature
    start2, end2 -- start and end nucleotides of the 2nd SeqFeature
    circular_length -- length of a circular nucleic acid; nucleic acid supposed linear if None
    """
    if circular_length is None:
        x = min(abs(start1 - end2),
                abs(end1 - start2))
    else:
        x = min(abs(start1 - end2),
                abs(end1 - start2),
                abs(circular_length - end1 + start2),
                abs(circular_length - end2 + start1))
    return x


def circular_slice(s, left, right):
    """Slice a string as circular."""
    if left < 0:
        return s[left:] + s[:right]
    elif right > len(s):
        return s[left:] + s[:right - len(s)]
    else:
        return s[left:right]


def get_first(dict_arg, key):
    """Return first element of a dict item if possible. If not subscriptable, return the item."""
    get = dict_arg.get(key)
    try:
        return get[0]
    except:
        return get


def concatenate(folder, extension):
    """Concatenates all files of 'folder' into a single file of 'extension' extension."""
    print('folder', folder)
    record_list = os.listdir(database_path / folder)
    record_list.sort()

    with open(database_path / f'{folder}.{extension}', 'w') as outfile:
        for f in record_list:
            with open(database_path / folder / f, 'r') as infile:
                outfile.write(infile.read())


def concatenate_csv(folder):
    """Concatenates all files of 'folder' into a single csv file."""
    print('folder', folder)
    record_list = os.listdir(database_path / folder)
    record_list.sort()

    df_concat = pd.concat([pd.read_csv(database_path / folder / f) for f in record_list], ignore_index=True)
    df_concat.to_csv(database_path / f'{folder}.csv', index=False)


# ORG TREE
def prune_tree(tree: PhyloNode, keep: list) -> PhyloNode:
    """Remove nodes not listed in `keep`
    If `keep` contains 'leaf', tips of the tree are not removed.
    Return prunned tree."""
    tree2 = tree.copy()

    # list of nodes to be discarded
    to_prune = []

    # iterate nodes and add their names into to_prune list if their rank is not in to prune
    for node in tree2.traverse():
        rank = get_rank(node.name)
        if 'leaf' in keep and node.is_leaf():  # node is removed if it is leaf and `keep` doe not contain `leaf`
            to_prune.append(node.name)
        if rank in keep:
            to_prune.append(node.name)

    tree2.prune(to_prune)

    return tree2


def export_tree(tree: PhyloNode, path) -> None:
    """Export tree as .nwk to specified path."""
    nwk_string = tree.write(format=1)
    with open(path, 'w') as out_file:
        out_file.write(nwk_string)
    return None


def export_annotation(tree: PhyloNode, path) -> None:
    """Export a csv annotation for a tree to specified path.
    The csv annotation contains [taxid,name,rank] for each node of the tree"""

    with open(path, 'w') as out_file:
        out_file.write('taxid;name;rank\n')
        for node in tree.traverse():
            taxid = node.name
            name = get_taxid_name(taxid)
            rank = get_rank(taxid)

            if node.is_leaf():    # all leaves get the "species" rank. it is done for simplicity
                rank = 'species'  # it is now compatible with the  R scripts. TODO: ranking at the strain level

            if name == 'root':
                rank = 'root'
            out_file.write(f'{taxid};{name};{rank}\n')
    return None


def get_taxid_name(taxid: int) -> str:
    """Return name of taxid.
    Return 'missing' if taxid is missing in the database"""
    name = ncbi.get_taxid_translator([taxid])
    name = list(name.values())
    if len(name) == 1:
        name = name[0]
    else:
        name = 'missing'
    return name


def get_rank(taxid: int) -> str:
    """Return rank of taxid."""
    rank = ncbi.get_rank([taxid])
    rank = list(rank.values())
    if len(rank) == 1:
        rank = rank[0]
    else:
        rank = 'missing'
    return rank


# constants
GENOMES_LOCATION = Path("genomes/")      # folder containing genome collections to construct a database from
DATABASES_LOCATION = Path("../databases/")  # folder containing databases
GENOME_EXTENSION = '.gbff'               # used to filter out genome files
SYMBOLS = string.digits + string.ascii_uppercase  # symbols used for generating headers
UTR_WINDOW = 200                         # window for recording 3' and 5' UTRs
CONTEXT_WINDOW = 10000                   # window for recording genomic context

# parse arguments
# expected arguments: 1) path to genomes; 2) path to metadata; 3) path to the database

parser = ArgumentParser()

parser.add_argument('-h', '--help', help='Show this help message', action='help')
parser.add_argument('genomes', help='Path to the directory containing genomes')
parser.add_argument('metadata', help='Path to the file containing metadata')

parser.add_argument('database',
                    help='Path to the directory with the output database. If doesn\'t exist, will be created')

parser.add_argument('--nocontext', help='Do not compute genome context of the features',
                    action='store_true', default=False)

args: Namespace = parser.parse_args()

database_path = Path(args.database)


# make temporary folders and an assembly progress file
if not os.path.exists(database_path):
    os.makedirs(database_path)
    for i in ['protein', 'upstream', 'sequence', 'downstream', 'translation', 'annotation']:
        os.makedirs(database_path / i)
    open(database_path / 'completed_paths.txt', 'a').close()

log = open(database_path / 'database_assembly_log.txt', 'a')  # log
print(f'Database assembly. db: {database_name}, folders: {genome_folders}')
log.write(f'Database assembly. db: {database_name}, folders: {genome_folders}\n')

# generate a list of all paths of genomes to be processed
genome_paths = []
for genome_folder in genome_folders:  # iterate genome folders
    genome_paths_in_folder = os.listdir(GENOMES_LOCATION / genome_folder)
    genome_paths_in_folder = [i for i in genome_paths_in_folder if GENOME_SIGNATURE in i]
    genome_paths_in_folder = [GENOMES_LOCATION / genome_folder / i for i in genome_paths_in_folder]
    genome_paths += genome_paths_in_folder

# generate a list paths of genomes completed in previous runs
with open(database_path / 'completed_paths.txt', 'r') as f:
    completed_paths = f.readlines()
completed_paths = [i.strip() for i in completed_paths]
completed_paths = [Path(i) for i in completed_paths]

# infer genomes to be completed
print(f'total genomes {len(genome_paths)}, completed genomes {len(completed_paths)}')
log.write(f'total genomes {len(genome_paths)}, completed genomes {len(completed_paths)}\n')
genome_paths = [i for i in genome_paths if i not in completed_paths]
genome_paths.sort()
completed_paths = open(database_path / 'completed_paths.txt', 'a')

# generated IDs will consist of a prefix and an iterated part
id_prefix = 'AB'
id_iterator = itertools.product(SYMBOLS, repeat=8)

# open metadata
metadata = pd.read_csv(Path(metadata_path) / 'metadata.csv', index_col='accession')

# make the list of columns
columns = ['lcs', 'assembly'] + metadata.columns.to_list() + ['replicon_type', 'replicon'] + \
          ['feature_type', 'gene', 'product', 'start', 'end', 'strand', 'protein_length']

iteration = 1
# iterate genomes
for genome_path in genome_paths:
    log.write(f'{iteration} {genome_path}\n')
    print(f'{iteration} {genome_path}')
    iteration += 1
    # iterate seq. records (replicons) in a genbank file
    for seq_record in SeqIO.parse(genome_path, 'genbank'):

        seq_record_data = []
        accession = genome_path.name
        accession = accession.split('_')[0] + accession.split('_')[1]  # genome name, e.g. GCF_000701165.1
        log.write('Assembly: ' + accession + ' ')

        accession = accession[:3] + '_' + accession[3:]

        # open output files
        ## fasta file with all protein sequences
        protein_output = open(database_path / 'protein' / f'{accession}', 'w')
        ## other annotations discribed
        upstream_out = open(database_path / 'upstream' / f'{accession}', 'w')
        sequence_out = open(database_path / 'sequence' / f'{accession}', 'w')
        downstream_out = open(database_path / 'downstream' / f'{accession}', 'w')
        translation_out = open(database_path / 'translation' / f'{accession}', 'w')
        annotation_out = open(database_path / 'annotation' / f'{accession}', 'w')

        definition = seq_record.description  # Bifidobacterium breve strain NRBB18 chromosome, complete genome
        log.write('Definition: ' + definition + '\n')

        if seq_record.annotations['topology'] == 'circular':
            circular_length = len(seq_record.seq)
        else:
            circular_length = None

        taxid = 'NONE'
        replicon = 'NONE'
        replicon_metadata = ['NONE', 'NONE']
        for feature in seq_record.features:  # extract taxid and replicon from the source feature
            # "source" feature is single per seq_record/replicon and contains replicon related information
            # this feature is not written to the database; its properties are assigned to every database entry

            # EXTRACT GENOME METADATA
            genome_metadata = [accession] + metadata.loc[[accession]].values.flatten().tolist()

            # EXTRACT REPLICON METADATA

            if feature.type == 'source':
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

                replicon_metadata = [replicon_type, replicon]

            # EXTRACT FEATURE METADATA

            elif feature.type not in ['source', 'gene']:
                # generate lcs id
                lcs = next(id_iterator)
                lcs = ''.join(lcs)
                lcs = id_prefix + lcs

                # feature type
                if feature.type == 'repeat_region':
                    feature_type = get_first(feature.qualifiers, 'rpt_family')
                elif feature.type == 'regulatory':
                    feature_type = get_first(feature.qualifiers, 'regulatory_class')
                elif feature.type == 'ncRNA':
                    feature_type = get_first(feature.qualifiers, 'ncRNA_class')
                else:
                    feature_type = feature.type

                # coordinates
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                coordinates = [start, end, strand]

                # sequence
                seq_record_seq = str(seq_record.seq)
                upstream = circular_slice(seq_record_seq, start - UTR_WINDOW, start)
                sequence = seq_record_seq[start:end]
                downstream = circular_slice(seq_record_seq, end, end + UTR_WINDOW)
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

        seq_record_data = pd.DataFrame(seq_record_data, columns=columns)

        """
        # INFER GENOME CONTEXT
        context = ''

        short_data = seq_record_data[['lcs', 'start', 'end']]
        short_data = short_data.assign(context='')

        length = len(short_data.index)
        for i in range(length):
            short = False
            done = False
            j = i
            while not done:
                if j + 1 < length:
                    j += 1
                else:
                    j = 0
                
                start1 = short_data.iat[i, 1]
                end1 = short_data.iat[i, 2]
                start2 = short_data.iat[j, 1]
                end2 = short_data.iat[j, 2]

                if distance(start1, end1, start2, end2, circular_length) <= CONTEXT_WINDOW:
                    context += short_data.iat[j, 0] + ';'
                else:
                    done = True

                if j == i:
                    done = True
                    short = True

            if not short:
                done = False
                j = i - 1
                while not done:
                    start1 = short_data.iat[i, 1]
                    end1 = short_data.iat[i, 2]
                    start2 = short_data.iat[j, 1]
                    end2 = short_data.iat[j, 2]

                    if distance(start1, end1, start2, end2, circular_length) <= CONTEXT_WINDOW:
                        context += seq_record_data.iat[j, 0] + ';'
                    else:
                        done = True
                    j -= 1

            context = context[:-1]  # remove last comma
            short_data.iat[i, 3] = context
            context = ''
        
        seq_record_data = pd.merge(seq_record_data, short_data[['ID', 'context']], on='lcs')
        """
        annotation_out.write(seq_record_data.to_csv(index=False))
        annotation_out.close()
        upstream_out.close()
        sequence_out.close()
        downstream_out.close()
        translation_out.close()
        protein_output.close()
    completed_paths.write(str(genome_path) + '\n')

log.close()

concatenate('protein', 'faa')

concatenate('upstream', 'csv')
concatenate('sequence', 'csv')
concatenate('downstream', 'csv')
concatenate('translation', 'csv')

FOLDERS_TO_CONCATENATE_CSV = ['annotation']
for folder in FOLDERS_TO_CONCATENATE_CSV:
    concatenate_csv(folder)

# GETTING TAXIDS
data_path = database_path / 'annotation.csv'
df = pd.read_csv(data_path)

# df = df[df.gtdb_taxonomy.notnull()]

taxids = df['taxid'].unique()

# MAKE ORG TREE

# NCBI taxonomy database object
ncbi = NCBITaxa()

# get tree topology as PhyloNode object
tree = ncbi.get_topology(taxids, intermediate_nodes=True)

tree_full = prune_tree(tree, ['leaf', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_genus = prune_tree(tree, ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_family = prune_tree(tree, ['family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_order = prune_tree(tree, ['order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_class = prune_tree(tree, ['class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_phylum = prune_tree(tree, ['phylum', 'superkingdom', 'kingdom', 'root'])

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
export_annotation(tree_full, out_folder / 'org_tree_full_data.csv')
