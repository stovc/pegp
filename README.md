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
  - https://mafft.cbrc.jp/alignment/software/linux.html - installation (recommend usinf .deb for installation)
- trimAl v1.2 `conda install -c bioconda trimal`
- hmmer `conda install -c bioconda hmmer`
- iqtree v2.x `conda install -c bioconda iqtree`

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
- Run `bash utility_scripts/mk_pfam.sh` from the root of the repo

## Running

#### Download genomes

`bash prepare_genomes/download_genome_metadata.sh`
`python3 prepare_genomes/filter_gtdb_metadata.py`
`bash prepare_genomes/download_genomes.sh`

#### Build a database
`python3 database_building/make_db GENOMES_FOLDER_PATH METADATA_PATH OUTPUT_DATABASE_PATH`
sample:
`python3 database_building/make_db.py test/genomes test/genomes/metadata.tsv databases/test`

#### Run an analysis

run `python3 pegp.py`

###### Create a project
For a test project:

n project test/hmms/clpP_TIGR00493.1.HMM

###### Running the analysis steps

1. s project 1
2. s project 2
3. s project 3 
4. s project 4 70 0.05
5. s project 5 0.9
6. s project 6
7. s project 7
8. s project 8
9. s project 9
10. s project 10
