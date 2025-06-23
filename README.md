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
