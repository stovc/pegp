WORK IN PROGRESS

# pegp

PEGP (protein evolutionary genomics pipeline). A software tool for
- protein homology search
- phylogenetic tree construction with annotations (taxonomy, domain architecture, genome context)
- reconstruction of phyletic pattern with identification of paralogs

## Dependencies
### Software
- ncbi-blast+
  - https://www.ncbi.nlm.nih.gov/books/NBK52640/ - installation and setup
- CD-HIT
  - https://github.com/weizhongli/cdhit/wiki/2.-Installation - installation
- MAFFT
  - https://mafft.cbrc.jp/alignment/software/linux.html - installation (recommend usinf .deb for installation)
- trimAl v1.2
  - http://trimal.cgenomics.org/downloads - installation
- hmmer
  - `conda install -c bioconda hmmer`
- iqtree v2.x
  - `conda install -c bioconda iqtree`

### Python packages
- Python 3.9
- BioPython
- ete3=3.1.2
- numpy
- pandas
- matplotlib
- seaborn
- prettytable
- colorama

## Installation
- `git clone https://github.com/stovc/pegp`
- Run `bash utility_scripts/mk_pfam.sh` from the root of the repo
