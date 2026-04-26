import string

# constants
GENOME_EXTENSION = '.gbff'               # used to filter out genome files
ID_PREFIX = 'lc_'
ID_SUFFIX_SYMBOLS = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'  # symbols used for generating headers
ID_SUFFIX_LENGTH = 10

UTR_WINDOW_SIZE = 200                         # window for recording 3' and 5' UTRs
CONTEXT_WINDOW_SIZE = 10000                   # window for recording genomic context
GENOME_METADATA_COLUMNS_OF_INTEREST = ['accession', 'gtdb_taxonomy']
COLUMNS_ID = ['lcs', 'assembly']
COLUMNS_REPLICON = ['replicon_type', 'replicon']
COLUMNS_FEATURE = ['feature_type', 'gene', 'product', 'start', 'end', 'strand', 'protein_length']

