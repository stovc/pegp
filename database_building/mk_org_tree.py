"""Generate species trees as .nwk and corresponding annotation data as .csv from `taxid.txt` related to the database.
The trees differ by the taxonomic rank to which they are collapsed

input - taxid.txt
output - org_tree_[TAXONOMIC RANK].nwk + 'org_tree_[TAXONOMIC RANK]_data.csv' where
the taxonomic ranks are [full, genus, family, order, class, phylum]

TODO: The job of `mk_org_tree.py` should be done in `mk_db.py`
"""

import os
import sys
from pathlib import Path
from ete3 import NCBITaxa, PhyloNode


def prune_tree(tree: PhyloNode, keep: list) -> PhyloNode:
    """Remove nodes not listed in `keep`
    If `keep` contains 'leaf', tips of the tree are not removed.
    Return prunned tree."""
    tree2 = tree.copy()

    # list of nodes to be discarded
    to_prune = []

    # iterate nodes and add their names into to_prune list if their rank is not in to prune
    for node in tree2.traverse():
        rank = get_rank(node.name)
        if 'leaf' in keep and node.is_leaf():  # node is removed if it is leaf and `keep` doe not contain `leaf`
            to_prune.append(node.name)
        if rank in keep:
            to_prune.append(node.name)

    tree2.prune(to_prune)

    return tree2


def export_tree(tree: PhyloNode, path) -> None:
    """Export tree as .nwk to specified path."""
    nwk_string = tree.write(format=1)
    with open(path, 'w') as out_file:
        out_file.write(nwk_string)
    return None


def export_annotation(tree: PhyloNode, path) -> None:
    """Export annotation .csv annotation for tree to specified path.
    The annotation .csv contains [taxid,name,rank] for each node of the tree"""

    with open(path, 'w') as out_file:
        out_file.write('taxid;name;rank\n')
        for node in tree.traverse():
            taxid = node.name
            name = get_taxid_name(taxid)
            rank = get_rank(taxid)

            if node.is_leaf():    # all leaves get the "species" rank. it is done for simplicity
                rank = 'species'  # it is now compatible with the  R scripts. TODO: ranking at the strain level

            if name == 'root':
                rank = 'root'
            out_file.write(f'{taxid};{name};{rank}\n')
    return None


def get_taxid_name(taxid: int) -> str:
    """Return name of taxid.
    Return 'missing' if taxid is missing in the database"""
    name = ncbi.get_taxid_translator([taxid])
    name = list(name.values())
    if len(name) == 1:
        name = name[0]
    else:
        name = 'missing'
    return name


def get_rank(taxid: int) -> str:
    """Return rank of taxid."""
    rank = ncbi.get_rank([taxid])
    rank = list(rank.values())
    if len(rank) == 1:
        rank = rank[0]
    else:
        rank = 'missing'
    return rank


if __name__ == '__main__':

    database = sys.argv[1]
    database_path = Path('../databases/') / database

    # NCBI taxonomy database object
    ncbi = NCBITaxa()

    # read taxids from txt to a list
    input_path = database_path / 'taxids.txt'
    with open(input_path, 'r') as input_file:
        taxids = input_file.readlines()

    # get tree topology as PhyloNode object
    tree = ncbi.get_topology(taxids, intermediate_nodes=True)

    tree_full = prune_tree(tree, ['leaf', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
    tree_genus = prune_tree(tree, ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
    tree_family = prune_tree(tree, ['family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
    tree_order = prune_tree(tree, ['order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
    tree_class = prune_tree(tree, ['class', 'phylum', 'superkingdom', 'kingdom', 'root'])
    tree_phylum = prune_tree(tree, ['phylum', 'superkingdom', 'kingdom', 'root'])

    # EXPORT TREES AND ANNOTATIONS
    # construct path for the output folder and make it
    out_folder = database_path / 'org_trees'
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # export trees pruned to different taxonomic levels
    export_tree(tree_full, out_folder / 'org_tree_full.nwk')
    export_tree(tree_genus, out_folder / 'org_tree_genus.nwk')
    export_tree(tree_family, out_folder / 'org_tree_family.nwk')
    export_tree(tree_order, out_folder / 'org_tree_order.nwk')
    export_tree(tree_class, out_folder / 'org_tree_class.nwk')
    export_tree(tree_phylum, out_folder / 'org_tree_phylum.nwk')

    # export .csv annotations for the trees pruned to different taxonomic levels
    export_annotation(tree_full, out_folder / 'org_tree_full_data.csv')
