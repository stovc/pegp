import string

# constants
GENOME_EXTENSION = '.gbff'               # used to filter out genome files
HEADER_SYMBOLS = string.digits + string.ascii_uppercase  # symbols used for generating headers
UTR_WINDOW_SIZE = 200                         # window for recording 3' and 5' UTRs
CONTEXT_WINDOW_SIZE = 10000                   # window for recording genomic context
GENOME_METADATA_COLUMNS_OF_INTEREST = ['refseq_category', 'taxid', 'organism_name',
                                       'genome_size', 'gc_percent', 'replicon_count', 'total_gene_count']

