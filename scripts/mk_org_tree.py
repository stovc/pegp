"""Generates species trees and corresponding annotation data from `taxid.txt` related to the database.
The trees are collapsed to different taxonomic levels
"""

from ete3 import NCBITaxa
ncbi = NCBITaxa()
from pathlib import Path


FOLDER = Path('../org_tree')


def get_rank(taxid):
    """Return rank of taxid"""
    rank = ncbi.get_rank([taxid])
    rank = list(rank.values())
    if len(rank) == 1:
        rank = rank[0]
    else:
        rank = 'missing'
    return rank


def get_name(taxid):
    """Return name of taxid"""
    name = ncbi.get_taxid_translator([taxid])
    name = list(name.values())
    if len(name) == 1:
        name = name[0]
    else:
        name = 'missing'
    return name


def collapse_leaf(node):
    keep = ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root']
    rank = get_rank(node.name)

    if rank in keep:
       return True
    else:
       return False


def cleanup_tree(tree, keep):
    tree2 = tree.copy()

    to_delete = []
    for n in tree2.traverse():
        rank = get_rank(n.name)
        if 'leaf' not in keep and n.is_leaf():
            to_delete.append(n)
        if rank not in keep and not n.is_leaf():
            to_delete.append(n)

    for node in to_delete:
        node.delete()

    return tree2


def prune_tree(tree, keep):
    tree2 = tree.copy()

    to_prune = []
    for n in tree2.traverse():
        rank = get_rank(n.name)
        if 'leaf' in keep and n.is_leaf():
            to_prune.append(n.name)
        if rank in keep:
            to_prune.append(n.name)

    tree2.prune(to_prune)

    return tree2

def export_tree(tree, name):
    nwk_string = tree.write(format=1)
    out = open(FOLDER / name, 'w')
    out.write(nwk_string)
    out.close()

    return None


def export_annotation(tree, name):
    out = open(FOLDER / name, 'w')
    out.write('taxid;name;rank\n')
    for n in tree.traverse():
        rank = get_rank(n.name)
        taxid = n.name
        name = get_name(taxid)
        rank = get_rank(taxid)
        if rank not in ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom']:
            rank = 'species'
        if name == 'root':
            rank = 'root'
        out.write(f'{taxid};{name};{rank}\n')
    out.close()

    return None


inp = open('../taxids.txt')
taxids = inp.readlines()
inp.close()

tree = ncbi.get_topology(taxids)
print('constructed tree from taxids')

tree_full = prune_tree(tree, ['leaf', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_genus = prune_tree(tree, ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_family = prune_tree(tree, ['family', 'order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_order = prune_tree(tree, ['order', 'class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_class = prune_tree(tree, ['class', 'phylum', 'superkingdom', 'kingdom', 'root'])
tree_phylum = prune_tree(tree, ['phylum', 'superkingdom', 'kingdom', 'root'])
print('prunned trees')

export_tree(tree_full, 'org_tree_full.nwk')
export_tree(tree_genus, 'org_tree_genus.nwk')
export_tree(tree_family, 'org_tree_family.nwk')
export_tree(tree_order, 'org_tree_order.nwk')
export_tree(tree_class, 'org_tree_class.nwk')
export_tree(tree_phylum, 'org_tree_phylum.nwk')
print('exported trees')

export_annotation(tree_full, 'org_tree_full_data.csv')
export_annotation(tree_genus, 'org_tree_genus_data.csv')
export_annotation(tree_family, 'org_tree_family_data.csv')
export_annotation(tree_order, 'org_tree_order_data.csv')
export_annotation(tree_class, 'org_tree_class_data.csv')
export_annotation(tree_phylum, 'org_tree_phylum_data.csv')
print('exported annotation data')


