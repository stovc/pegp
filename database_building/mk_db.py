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
"""

import os, sys
from pathlib import Path
import string
import itertools
from Bio import Seq, SeqIO
from Bio import SeqFeature
import pandas as pd


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


# constants
GENOMES_LOCATION = Path("genomes/")      # folder containing genome collections to construct a database from
DATABASES_LOCATION = Path("../databases/")  # folder containing databases
GENOME_SIGNATURE = '.gbff'               # used to filter out genome files
SYMBOLS = string.digits + string.ascii_uppercase  # symbols used for generating headers
UTR_WINDOW = 200                         # window for recording 3' and 5' UTRs
CONTEXT_WINDOW = 10000                   # window for recording genomic context
COLUMNS = ['ID', 'locus_tag',
           'assembly', 'accession',
           'start', 'end', 'strand',
           'taxid', 'replicon',
           'feature', 'gene', 'product', "length", 'protein_id']  # columns in the annotation dataframe

# parse arguments
# expected arguments: 1) database name; 2-N) genome folders to extract the genomes from
arguments = sys.argv
database_name = arguments[1]
genome_folders = arguments[2:]
database_path = Path(DATABASES_LOCATION) / database_name

# make temporary folders and an assembly progress file
if not os.path.exists(database_path):
    os.makedirs(database_path)
    for i in ['protein', 'taxid_map', 'upstream', 'sequence', 'downstream', 'translation', 'annotation']:
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
f = open(database_path / 'completed_paths.txt', 'r')
completed_paths = f.readlines()
f.close()
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

iteration = 1
# iterate genomes
for genome_path in genome_paths:
    log.write(f'{iteration} {genome_path}\n')
    print(f'{iteration} {genome_path}')
    iteration += 1
    # iterate seq. records (replicons) in a genbank file
    for seq_record in SeqIO.parse(genome_path, 'genbank'):

        seq_record_data = []
        assembly = genome_path.name
        assembly = assembly.split('_')[0] + assembly.split('_')[1]  # genome name, e.g. GCF_000701165.1
        log.write('Assembly: ' + assembly + ' ')

        accession = seq_record.id  # genome id in NCBI, e.g. NZ_CP023193.1
        log.write('Accession: ' + accession + ' ')

        # open output files
        ## fasta file with all protein sequences
        protein_output = open(database_path / 'protein' / f'{assembly}-{accession}', 'w')
        ## taxid map file with taxids for BLAST database
        map_output = open(database_path / 'taxid_map' / f'{assembly}-{accession}', 'w')
        ## other annotations discribed
        upstream_out = open(database_path / 'upstream' / f'{assembly}-{accession}', 'w')
        sequence_out = open(database_path / 'sequence' / f'{assembly}-{accession}', 'w')
        downstream_out = open(database_path / 'downstream' / f'{assembly}-{accession}', 'w')
        translation_out = open(database_path / 'translation' / f'{assembly}-{accession}', 'w')
        annotation_out = open(database_path / 'annotation' / f'{assembly}-{accession}', 'w')

        definition = seq_record.description  # Bifidobacterium breve strain NRBB18 chromosome, complete genome
        log.write('Definition: ' + definition + '\n')

        if seq_record.annotations['topology'] == 'circular':
            circular_length = len(seq_record.seq)
        else:
            circular_length = None

        taxid = 'NONE'
        replicon = 'NONE'
        for feature in seq_record.features:  # extract taxid and replicon from the source feature
            if feature.type == 'source':

                db_xref_list = feature.qualifiers.get('db_xref')  # extract taxid from the source feature
                for db_xref in db_xref_list:
                    if db_xref.split(':')[0] == 'taxon':  # 'taxon:1685'
                        taxid = db_xref.split(':')[1]     # 1685
                log.write('Taxid: ' + taxid + '\n')

                plasmid = feature.qualifiers.get('plasmid')
                log.write(f'Plasmid: {plasmid}\n')
                segment = feature.qualifiers.get('segment')
                log.write(f'Segment: {segment}\n')
                chromosome = feature.qualifiers.get('chromosome')
                log.write(f'Chromosome: {chromosome}\n')

                if plasmid is not None:
                    replicon = 'plasmid_' + plasmid[0]
                elif segment is not None:
                    replicon = 'segment_' + segment[0]
                elif chromosome is not None:
                    replicon = 'chromosome_' + chromosome[0]
                else:
                    replicon = 'main'

            elif feature.type not in ['source', 'gene']:
                ID = next(id_iterator)
                ID = ''.join(ID)
                ID = id_prefix + ID

                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand

                seq_record_seq = str(seq_record.seq)
                upstream = circular_slice(seq_record_seq, start - UTR_WINDOW, start)
                sequence = seq_record_seq[start:end]
                downstream = circular_slice(seq_record_seq, end, end + UTR_WINDOW)

                # define feature type
                if feature.type == 'repeat_region':
                    feature_type = get_first(feature.qualifiers, 'rpt_family')
                elif feature.type == 'regulatory':
                    feature_type = get_first(feature.qualifiers, 'regulatory_class')
                elif feature.type == 'ncRNA':
                    feature_type = get_first(feature.qualifiers, 'ncRNA_class')
                else:
                    feature_type = feature.type

                locus_tag = get_first(feature.qualifiers, 'locus_tag')
                gene = get_first(feature.qualifiers, 'gene')
                product = get_first(feature.qualifiers, 'product')
                translation = get_first(feature.qualifiers, 'translation')
                protein_id = get_first(feature.qualifiers, 'protein_id')

                if translation is not None:
                    protein_output.write(f'>{ID}\n{translation}\n')  # write a fasta record with translation
                    map_output.write(ID + ' ' + taxid + '\n')
                    translation_out.write(f'{ID},{translation}\n')
                    protein_length = len(translation)
                else:
                    protein_length = None

                upstream_out.write(f'{ID},{upstream}\n')
                sequence_out.write(f'{ID},{sequence}\n')
                downstream_out.write(f'{ID},{downstream}\n')

                annotation = [ID, locus_tag,
                              assembly, accession, start, end, strand,
                              taxid, replicon,
                              feature_type, gene,
                              product, protein_length, protein_id]
                seq_record_data.append(annotation)

        seq_record_data = pd.DataFrame(seq_record_data, columns=COLUMNS)
        context = ''

        short_data = seq_record_data[['ID', 'start', 'end']]
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
        seq_record_data = pd.merge(seq_record_data, short_data[['ID', 'context']], on='ID')
        annotation_out.write(seq_record_data.to_csv(index=False))
        annotation_out.close()
        upstream_out.close()
        sequence_out.close()
        downstream_out.close()
        translation_out.close()
        protein_output.close()
        map_output.close()
    completed_paths.write(str(genome_path) + '\n')

log.close()

concatenate('protein', 'faa')
concatenate('taxid_map', 'txt')

concatenate('upstream', 'csv')
concatenate('sequence', 'csv')
concatenate('downstream', 'csv')
concatenate('translation', 'csv')

FOLDERS_TO_CONCATENATE_CSV = ['annotation']
for folder in FOLDERS_TO_CONCATENATE_CSV:
    concatenate_csv(folder)
