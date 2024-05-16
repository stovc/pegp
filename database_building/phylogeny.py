from ete3 import NCBITaxa, PhyloNode


def export_tree(tree: PhyloNode, path) -> None:
    """Export tree as .nwk to specified path."""
    nwk_string = tree.write(format=1)
    with open(path, 'w') as out_file:
        out_file.write(nwk_string)
    return None


def export_tree_annotation(tree: PhyloNode, path, database) -> None:
    """Export a csv annotation for a tree to specified path.
    The csv annotation contains [taxid,name,rank] for each node of the tree"""

    with open(path, 'w') as out_file:
        out_file.write('taxid;name;rank\n')
        for node in tree.traverse():
            taxid = node.name
            name = get_name_of_taxid(taxid, database)
            rank = get_rank_of_taxid(taxid, database)

            if node.is_leaf():    # all leaves get the "species" rank. it is done for simplicity
                rank = 'species'  # it is now compatible with the  R scripts. TODO: ranking at the strain level

            if name == 'root':
                rank = 'root'
            out_file.write(f'{taxid};{name};{rank}\n')
    return None


def get_name_of_taxid(taxid: int, database) -> str:
    """Return name of taxid.
    Return 'missing' if taxid is missing in the database"""
    name = database.get_taxid_translator([taxid])
    name = list(name.values())
    if len(name) == 1:
        name = name[0]
    else:
        name = 'missing'
    return name


def get_rank_of_taxid(taxid: int, database) -> str:
    """Return rank of taxid."""
    rank = database.get_rank([taxid])
    rank = list(rank.values())
    if len(rank) == 1:
        rank = rank[0]
    else:
        rank = 'missing'
    return rank


def prune_tree(tree: PhyloNode, keep: list, database) -> PhyloNode:
    """Remove nodes not listed in `keep`
    If `keep` contains 'leaf', tips of the tree are not removed.
    Return prunned tree."""
    tree2 = tree.copy()

    # list of nodes to be discarded
    to_prune = []

    # iterate nodes and add their names into to_prune list if their rank is not in to prune
    for node in tree2.traverse():
        rank = get_rank_of_taxid(node.name, database)
        if 'leaf' in keep and node.is_leaf():  # node is removed if it is leaf and `keep` doe not contain `leaf`
            to_prune.append(node.name)
        if rank in keep:
            to_prune.append(node.name)

    tree2.prune(to_prune)

    return tree2
