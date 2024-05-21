from Bio import SeqIO
import numpy as np
import pandas as pd

import csv
import logging
from pathlib import Path

from utils import circular_slice, get_first, distance_between_sequences
import constants as c

logger = logging.getLogger(__name__)

class GenomeParserConfig:
    def __init__(self, database_path, genome_metadata_df, id_iterator, enable_genome_context):
        self.enable_genome_context = enable_genome_context
        self.database_path = database_path
        self.genome_metadata_df = genome_metadata_df
        self.id_iterator = id_iterator
        self.enable_genome_context = enable_genome_context


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
    logger.info(f'Plasmid: {plasmid}\n')
    segment = feature.qualifiers.get('segment')
    logger.info(f'Segment: {segment}\n')
    chromosome = feature.qualifiers.get('chromosome')
    logger.info(f'Chromosome: {chromosome}\n')

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


def parse_genome(genome_path: Path, config: GenomeParserConfig) -> None:
    """Parse features in a genbank genome file. Save protein sequences and metadata to files."""
    logger.info(f'Parsing genome: {genome_path}\n')

    filename = genome_path.name                    # e.g., GCF_000012525.1_ASM1252v1_genomic.gbff
    accession = '_'.join(filename.split('_')[:2])  # e.g., GCF_000701165.1

    logger.info(f'Assembly: {accession}\n')

    annotation_path = Path(config.database_path / 'annotation' / f'{accession}')

    columns_genome = config.genome_metadata_df.columns.to_list()
    annotation_columns = c.COLUMNS_ID + columns_genome + c.COLUMNS_REPLICON + c.COLUMNS_FEATURE

    for seq_record in SeqIO.parse(genome_path, 'genbank'):
        metadata_main = []
        metadata_sequence = [['lcs', 'upstream', 'sequence', 'downstream']]
        metadata_translation = [['lcs', 'translation']]

        # open output files
        protein_fasta_file = open(config.database_path / 'protein' / f'{accession}', 'w')

        replicon_metadata = ['NONE', 'NONE']

        genome_metadata = [accession] + config.genome_metadata_df.loc[[accession]].values.flatten().tolist()

        for feature in seq_record.features:  # extract taxid and replicon from the source feature
            # "source" feature is single per seq_record/replicon and contains replicon related information
            # this feature is not written to the database; its properties are assigned to every database entry

            # EXTRACT REPLICON METADATA

            if feature.type == 'source':
                replicon_metadata = get_replicon_data_from_source(feature)

            # EXTRACT FEATURE METADATA

            elif feature.type not in ['source', 'gene']:
                # generate lcs id
                lcs = next(config.id_iterator)

                feature_type = get_feature_type(feature)

                # annotations
                gene = get_first(feature.qualifiers, 'gene')
                product = get_first(feature.qualifiers, 'product')

                # coordinates
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand

                # sequence
                seq_record_seq = str(seq_record.seq)
                upstream = circular_slice(seq_record_seq, start - c.UTR_WINDOW_SIZE, start)
                sequence = seq_record_seq[start:end]
                downstream = circular_slice(seq_record_seq, end, end + c.UTR_WINDOW_SIZE)
                translation = get_first(feature.qualifiers, 'translation')

                if len(sequence) > 100000:  # WORKAROUND
                    sequence = sequence[:100000]

                # write protein sequence data
                if translation is not None:
                    protein_fasta_file.write(f'>{lcs}\n{translation}\n')  # write a fasta record with translation
                    metadata_translation.append([lcs, translation])

                    protein_length = len(translation)
                else:
                    protein_length = np.nan

                # WRITE METADATA
                # Write main metadata
                feature_metadata = [feature_type, gene, product, start, end, strand, protein_length]
                metadata_main_row = [lcs] + genome_metadata + replicon_metadata + feature_metadata
                metadata_main.append(metadata_main_row)

                # write sequence metadata
                metadata_sequence.append([lcs, upstream, sequence, downstream])

        # make a dataframe with metadata
        metadata_main_df = pd.DataFrame(metadata_main, columns=annotation_columns)

        # INFER GENOME CONTEXT
        if config.enable_genome_context:
            metadata_main_df = add_genome_context(seq_record, metadata_main_df)

        # SAVE METADATA TO FILES
        metadata_main_df.to_csv(annotation_path, index=False)

        def export_list_to_csv(list_to_export, path):
            with open(path, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(list_to_export)

        sequence_csv_path = config.database_path / 'sequence' / f'{accession}'
        export_list_to_csv(metadata_sequence, sequence_csv_path)

        translation_csv_path = config.database_path / 'translation' / f'{accession}'
        export_list_to_csv(metadata_translation, translation_csv_path)

        protein_fasta_file.close()

    return None
