WORK IN PROGRESS

# pegp

PEGP (protein evolutionary genomics pipeline). A software tool for
- protein homology search
- phylogenetic tree construction with annotations (taxonomy, domain architecture, genome context)
- reconstruction of phyletic pattern with identification of paralogs

## Dependencies

### Software

- CD-HIT
  - https://github.com/weizhongli/cdhit/wiki/2.-Installation - installation
- MAFFT
  - https://mafft.cbrc.jp/alignment/software/linux.html - installation (recommend using .deb for installation)
- trimAl v1.2 `mamba install -c bioconda trimal`
- hmmer `mamba install -c bioconda hmmer`
- iqtree v2.x `mamba install -c bioconda iqtree`

### Python packages

- python 3.10
- biopython `mamba install biopython`
- ete3 `mamba install -c conda-forge ete3`
- matplotlib `mamba install matplotlib`
- numpy, pandas, seaborn `mamba install seaborn`
- prettytable `mamba install prettytable`
- colorama `mamba install colorama`
- reportlab `mamba install conda-forge::reportlab`

## Installation
- `git clone https://github.com/stovc/pegp`
- Run `bash utility_scripts/mk_pfam.sh` from the root of the repository

## Running

#### Download genomes

1. `bash prepare_genomes/download_genome_metadata.sh`
2. `python prepare_genomes/filter_gtdb_metadata.py --taxonomic-rank genus` (filter outputs need to be concatenated manually at the moment to meta.tsv at the moment)
3. `bash prepare_genomes/download_genomes.sh --out-dir data/genomes/r232_rs_complete_max-gen1_min-phy3`

#### Build a database
`python database_building/make_db GENOMES_FOLDER_PATH METADATA_PATH OUTPUT_DATABASE_PATH`
sample:
`python database_building/make_db.py test/genomes test/genomes/metadata.tsv databases/test`

#### Run an analysis

run `python3 pegp.py`

###### Create a project
For a test project:

- `n project test/hmms/clpP_TIGR00493.1.HMM`
- dir `reports` needs to be created manually at the moment
###### Running the analysis steps

1. `s project 1`
2. `s project 2`
3. `s project 3 `
4. `s project 4 70 0.00001` - query coverage percent, e-value cutoff
5. `s project 5 0.9` - cd-hit "-c", clustering threshold
6. `s project 6`
7. `s project 7`
8. `s project 8`
9. `s project 9`
10. `s project 10`
11. `s project 11`
12. `s project 12 0.5` - trimal "-gt", max fraction of gaps allowed
13. `s project 13`
14. `s project 14`
15. `s project 15`
16. `s project 16`
17. `s project 17`
