"""
Extra functionality for msprime which we need here. 
Most if not all of this file can be removed once 
https://github.com/jeromekelleher/msprime/issues/354
is addressed
"""
import numpy as np

def treestring(name, tree):
    return "TREE " + name + " = [&R] " + tree.newick(precision=14)[:-1] + ":0;\n"

def header(n_tips, node_labels):
    #convert back to 0-based tip labels, or use string node labels (escaped with single quotes)
    tip_map = [
        str(i + 1) + " '{}'".format(node_labels[i].replace("'","''") if i in node_labels else i) 
        for i in range(n_tips)]
    return "#NEXUS\nBEGIN TREES;\nTRANSLATE\n{};\n".format(",\n".join(tip_map))
        
def footer():
    return "END;\n"

def write_nexus_tree(tree, treefile, node_labels={}):
    """
    Writes out a single SparseTree from a tree sequence to a nexus tree file
    which allows tips to be relabelled so that they either correspond to the
    tskit numbers (0 based) or to a user specified set of names
    """
    treefile.write(header(tree.num_samples(), node_labels))
    treefile.write(treestring(str(tree.interval[1]), tree))
    treefile.write(footer())
    

def write_nexus_trees(
        ts, treefile, node_labels={}, tree_labels_between_variants=False):
    """
    Writes out all the trees in this tree sequence to a single nexus file
    (see https://doi.org/10.2307/2413497 for file format information)

    The node labels, if given, are expected to be a dict that maps integers
    from 0 (the tips) onto strings. As per the Nexus specification, these strings
    will be quoted in the nexus file, and single quotes escaped by doubling, so
    there is no need to do this yourself. Any node without a label will be allocated
    the tip number in the tree sequence.

    The names of each tree in the treefile are meaningful. They give the
    positions along the sequence where we switch to a new tree. These
    positions are half-open intervals, that is a tree named '3' will *not*
    include position 3, but will include position 2.999999. By default this means 
    if a tree sequence of (say) length 5 has been inferred from variants, say at
    position 3 and 4, the first tree will be named '4' (i.e. covering all positions
    up to but not including 4) and the second tree will be named '5' (the max
    sequence length). Thus variants at position 4 will fall into the tree labelled 5
    not the tree labelled 4 - this seems slightly odd.

    If tree_labels_between_variants is True, we change labelling so that the
    positions are halfway between the two nearest variants
    i.e. if there are variants with different trees at position 3 & 4, then
    the first tree has the name '3.5' and the second '5'. Thus we are assured
    that variant 3 lies in the tree labelled 3.5 (covering 0-3.5), and the
    variant 4 lies in the tree labelled 5.

    Note that inferred trees that have been simplified to include only a subset
    of tips may well have breakpoints (switches between trees) that do not 
    occur at the position of a variant on the subsampled trees.
    """
    treefile.write(header(ts.num_samples, node_labels))

    if tree_labels_between_variants:
        variant_pos = np.array([m.position for m in ts.mutations()])
        pos_between_vars = np.concatenate([[0],np.diff(variant_pos)/2+variant_pos[:-1],
                                          [ts.get_sequence_length()]])
        for t in ts.trees():
            # index by the average position between the two nearest variants
            av_upper_pos = pos_between_vars[np.searchsorted(variant_pos,t.get_interval()[1])]
            assert av_upper_pos <= t.get_interval()[1]
            assert t.num_roots == 1, \
                "Couldn't write Nexus to {} as Newick at interval {} has more than one root".format(
                treefile.name, t.get_interval())
            treefile.write(treestring(str(av_upper_pos), t))
    else:
        for t in ts.trees():
            # index by rightmost genome position
            assert t.num_roots == 1, \
                "Couldn't write Nexus to {} as Newick at interval {} has more than one root".format(
                treefile.name, t.get_interval())
            treefile.write(treestring(str(t.get_interval()[1]), t))
    treefile.write(footer())

def save_nexus_trees(ts, fn, **kwargs):
    """
    Same as write_nexus_trees() only use a file name not a file handle
    """
    with open(fn, "w+") as out:
        write_nexus_trees(ts, out, **kwargs)

def save_nexus_tree(tree, fn, **kwargs):
    """
    Same as write_nexus_tree() only use a file name not a file handle
    """
    with open(fn, "w+") as out:
        write_nexus_tree(tree, out, **kwargs)
