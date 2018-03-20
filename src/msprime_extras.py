"""
Extra functionality for msprime which we need here. 
Most if not all of this file can be removed once 
https://github.com/jeromekelleher/msprime/issues/354
is addressed
"""
import logging
import numpy as np

def write_nexus_trees(
        ts, treefile, tree_labels_between_variants=False):
    """
    Writes out all the trees in this tree sequence to a single nexus file.
    The names of each tree in the treefile are meaningful. They give the
    positions along the sequence where we switch to a new tree.

    The positions are half-open intervals, that is a tree named '3' will *not*
    include position 3, but will include position 2.999999. That means if a
    treesequence of (say) length 5 has been inferred from variants, say at
    position 3 and 4, the first tree will be named '4' (i.e. covering all positions
    up to but not including 4) and the second tree will be named with the sequence
    length '5'. Thus variants at position 4 will fall into the tree labelled 5
    not the tree labelled 4 - this seems slightly odd.

    If tree_labels_between_variants is True, we change this so that instead the
    positions are halfway between the two nearest variants
    i.e. if there are variants with different trees at position 3 & 4, then
    the first tree has the name '3.5' and the second '5'. Thus we are assured
    that variant 3 lies in the tree labelled 3.5 (covering 0-3.5), and the
    variant 4 lies in the tree labelled 5.

    Note that inferred trees that have been simplified to include only a subset
    of tips may well have breakpoints (switches between trees) that do not 
    occur at the position of a variant on the subsampled trees
    """
    print("#NEXUS\nBEGIN TREES;", file=treefile)
    tip_map = [
        "{} {}".format(i + 1, i) #convert back to 0-based tip labels
        for i in range(ts.get_sample_size())]
    print("TRANSLATE\n{};".format(",\n".join(tip_map)), file=treefile)
    variant_index = 0

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
            print("TREE " + str(av_upper_pos) + " = [&R] " + t.newick(precision=14,time_scale=(1/4))[:-1] + ":0;", file=treefile)
    else:
        for t in ts.trees():
            # index by rightmost genome position
            assert t.num_roots == 1, \
                "Couldn't write Nexus to {} as Newick at interval {} has more than one root".format(
                treefile.name, t.get_interval())
            print("TREE " + str(t.get_interval()[1]) + " = [&R] " + t.newick(precision=14,time_scale=(1/4))[:-1] + ":0;", file=treefile)
    print("END;", file=treefile)
