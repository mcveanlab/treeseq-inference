"""
Extra functionality for msprime which we need here.
"""
import logging

import numpy as np

import msprime
# TMP - remove once proper creation API is in place
import _msprime

def sparse_tree_to_newick(st, precision, Ne):
    """
    Converts the specified sparse tree to an ms-compatible Newick tree.
    """
    branch_lengths = {}
    root = st.get_root()
    stack = [root]
    while len(stack) > 0:
        node = stack.pop()
        if st.is_internal(node):
            for child in st.get_children(node):
                stack.append(child)
                length = (st.get_time(node) - st.get_time(child)) / (4 * Ne)
                s = "{0:.{1}f}".format(length, precision)
                branch_lengths[child] = s
    return _build_newick(root, root, st, branch_lengths)


def _build_newick(node, root, tree, branch_lengths):
    if tree.is_leaf(node):
        s = "{0}:{1}".format(node + 1, branch_lengths[node])
    else:
        subtrees = []
        for c in tree.get_children(node):
            s = _build_newick(c, root, tree, branch_lengths)
            subtrees.append(s)
        if node == root:
            # The root node is treated differently
            s = "({});".format(",".join(subtrees))
        else:
            s = "({}):{}".format(",".join(subtrees), branch_lengths[node])
    return s


def newick_trees(ts, precision=3, Ne=1):
    """
    Returns an iterator over the trees in the specified tree sequence.

    This is the same as the equivalent function in the msprime library,
    except it handles non binary trees and is much slower.
    """
    for t in ts.trees():
        l, r = t.interval
        yield r - l, sparse_tree_to_newick(t, precision, Ne)


def discretise_mutations(ts):
    """
    Takes the specified tree sequence returns a new tree sequence which discretises
    the mutations such that all positions are integers.

    Raises a ValueError if this is not possible without violating the underlying
    topology.
    """
    logging.warning("The current implementation of the discretise_mutations() function" +
        " creates sequences which do not fit the standard evolutionary model well." +
        " That causes inference problems for probabalistic inference programs" +
        " such as ARGweaver.\n" + "It is not recommended to test ARGweaver inference" +
        " on this dataset, especially for high mutation rates")
    breakpoints = sorted(set([r.left for r in ts.records()] + [ts.sequence_length]))
    old_mutations = list(ts.mutations())
    positions = []
    new_mutations = []
    k = 0
    for j in range(len(breakpoints) - 1):
        left, right = breakpoints[j], breakpoints[j + 1]
        interval_mutations = []
        while k < ts.num_mutations and old_mutations[k].position < right:
            interval_mutations.append(old_mutations[k])
            k += 1
        if len(interval_mutations) > right - left:
            raise ValueError(
                "Cannot discretise {} mutations in interval ({}, {})".format(
                    len(interval_mutations), left, right))
        # Start the new mutations at the beginning of the interval and place them
        # one by one at unit intervals. This will lead to very strange looking
        # spatial patterns, but will do until we implement a proper finite sites
        # mutation model in msprime.
        x = left
        for mut in interval_mutations:
            new_mutations.append(
                msprime.Mutation(position=x, node=mut.node, index=mut.index))
            x += 1
    records = list(ts.records())
    # Ugly!!! Fix this once the new creation API is in place.
    samples = [msprime.Sample(time=0, population=0) for _ in range(ts.sample_size)]
    ll_ts = _msprime.TreeSequence()
    ll_ts.load_records(samples=samples, records=records)
    ll_ts.set_mutations(new_mutations)
    return msprime.TreeSequence(ll_ts)


def write_nexus_trees(
        ts, treefile, tree_labels_between_variants=False, zero_based_tip_numbers=True):
    """
    Writes out all the trees in this tree sequence to a single nexus file.
    The names of each tree in the treefile are meaningful. They give the 
    positions along the sequence where we switch to a new tree.
    
    The positions are half-open intervals, that is a tree named '3' will *not* 
    include position 3, but will include position 2.999999. That means if a
    treesequence of (say) legnth 5 has been inferred from variants, say at 
    position 3 and 4, the first tree will be named '4' (i.e. covering all positions 
    up to but not including 4) and the second tree will be named with the sequence 
    length '5'. Thus variants at position 4 will fall into the tree labelled 5 
    not the tree labelled 4 - this seems slightly odd.
    
    If tree_labels_between_variants is True, we change this so that the 
    positions are instead, halfway between the two nearest variants
    i.e. if there are variants with different trees at position 3 & 4, then
    the first tree has the name '3.5' and the second '5'. Thus we are assured
    that variant 3 lies in the tree labelled 3.5 (covering 0-3.5), and the 
    variant 4 lies in the tree labelled 5.
    
    """
    print("#NEXUS\nBEGIN TREES;", file=treefile)
    increment = 0 if zero_based_tip_numbers else 1
    tip_map = [
        "{} {}".format(i + increment, i + increment)
        for i in range(ts.get_sample_size())]
    print("TRANSLATE\n{};".format(",\n".join(tip_map)), file=treefile)
    variant_index = 0

    assert zero_based_tip_numbers==False, 'At the moment, we can only output 1-based tip labels'
    if tree_labels_between_variants:
        variant_pos = np.array([m.position for m in ts.mutations()])
        pos_between_vars = np.concatenate([[0],np.diff(variant_pos)/2+variant_pos[:-1],
                                          [ts.get_sequence_length()]])
        for t, (_, newick) in zip(ts.trees(), newick_trees(ts)):
            # index by the average position between the two nearest variants
            av_upper_pos = pos_between_vars[np.searchsorted(variant_pos,t.get_interval()[1])]
            assert av_upper_pos <= t.get_interval()[1]
            print("TREE " + str(av_upper_pos) + " = " + newick, file=treefile)
    else:
        for t, (_, newick) in zip(ts.trees(), newick_trees(ts)):
            # index by rightmost genome position
            print("TREE " + str(t.get_interval()[1]) + " = " + newick, file=treefile)
    print("END;", file=treefile)
