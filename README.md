WORK IN PROGRESS

# PEGP

PEGP (protein evolutionary genomics pipeline). A software tool for analysis of evolutionary history of protein families on large evolutionary scales.

- protein homology search
- phylogenetic tree construction with annotations (taxonomy, domain architecture, genome context)
- reconstruction of phyletic pattern with identification of paralogs

## System requirements

Linux and macOS are supported.

For pipeline:

 - Python 3 (3.9 or higher)

For tree visualisation:
 - R (4.2.2)

## Dependencies

### Software

- CD-HIT v4.8.1
`conda install -c bioconda cd-hit`
- MAFFT v7.310
 `conda install -c bioconda mafft`
- trimAl v1.2 
    `conda install -c bioconda trimal`
- hmmer v3.3.2
    `conda install -c bioconda hmmer`
- iqtree v2.x 
    `conda install -c bioconda iqtree`

### Python packages

- biopython v1.79
`conda install biopython`
- ete3 v3.1.3
`conda install -c etetoolkit ete3`
- matplotlib v3.7.0
`conda install matplotlib`
- numpy v1.23.5, pandas v1.5.3, seaborn v0.12.2 
`conda install numpy pandas seaborn
- prettytable v3.6.0
`conda install prettytable`
- colorama v0.4.6
`conda install colorama`

### R pakages
 - tidytree
 - dplyr
 - ggplot2

Visualization of phylogenetic trees

- ggtree
- Gggenes

To help install packages and additional scripts of R, please, go to the link: https://github.com/stovc/tree-annotation. There is a link to the repository dedicated to the visualization of trees in R with a description of the dependencies installation.

## Installation
`git clone https://github.com/stovc/pegp`

## Database assembly
- Download GenBank (.gbff) genomes of interest into a directory 
- Build a database from genomes by running 
`cd database_building/`
`python3 mk_db.py [database_name] genomes_family/`

- Extract all taxids (taxonomic IDs of species our genomes belong to in NCBI Taxonomy database) from `annotation.csv` by running
`python3 get_taxids.py [database_name]`

- Build organism trees collapsed to different taxonomic ranks and corresponding dataframes by running
`python3 mk_org_tree.py [database_name]`
- Build blast database by running from the root of the repository:
`bash database_building/makeblastdb.sh [database_name]`

- If you want to run a test analysis, build HMM-profile by running `bash utility_scripts/mk_pfam.sh` from the root of the repo 

## Preparation of the query
Prepare your query protein domain.
- Download (or copy to clipboard) an amino acid sequence of the protein of interest from https://www.ncbi.nlm.nih.gov/gene/ 
- Predict its domain architecture at https://www.ebi.ac.uk/interpro/

- Subset the amino acid sequence corresponding to the domain of interest according to Pfam 

- Create a fasta amino acid file (.faa) with your domain
- Create a query HMM-profile from fasta file


## Run tool
Run the pipeline manager
`python3 pegp.py`


### Create a project
Create a new project by typing 
`n [project_name] [query_path]`

### Run the analysis
To run a step, type
`s [project_name] [step_number]`

Steps 4 (hits filtering) and 11 (triming) requere additional arguments:
`s [project_name] 4 [coverage_threshpld] [e-value_threshold]`
`s [project_name] 11 [threshold]`

After every step, a table with steps status will be printed:
![](https://hackmd.io/_uploads/B1oW6AOHh.png)

Only steps with status "Ready" should be run. 

All results of the scripts, including reports, are stored in `projects/[project_name]` folder.

Example of graph from the report (hits 3d plot of identity, query_coverage, and length colored by taxon):
![](https://hackmd.io/_uploads/S1xD00uBn.jpg)

Finally, in result tool running user get:
1. phylogenetic tree in NEWICK format
2. annotation.csv (contain information about genes of different species)
3. upstream an downstream sequences of genes 
4. taxon data for tree building (csv format)

Example of phylogenetic tree reconstruction (you can learn more by going to the R repository: https://github.com/stovc/tree-annotation):
![](https://hackmd.io/_uploads/SkgRb1tS2.png)

For assessment:

The scripts (written together with the supervisor) and example of data received by Evgeny Kovalenko student are located in a separated directory `Evgeny_work`
