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
- biopython `conda install biopython`
- ete3 `conda install -c etetoolkit ete3`
- matplotlib `conda install matplotlib`
- numpy, pandas, seaborn `conda install seaborn`
- prettytable `conda install prettytable`
- colorama `conda install colorama`

## Installation
- `git clone https://github.com/stovc/pegp`
- Run `bash utility_scripts/mk_pfam.sh` from the root of the repo
