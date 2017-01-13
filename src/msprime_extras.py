"""
Extra functionality for msprime which we need here.
"""

import numpy as np


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


# JK - I've been staring at this for 10 minutes, an I can't figure out what the
# index_trees_by_variants argument does here. It seems
# to be two mutually exclusive arguments, which give different behaviours? We have to
# simplify this. Maybe this should be two different functions? What do we use this
# behaviour for??

def write_nexus_trees(
        ts, treefile, index_trees_by_variants=False, zero_based_tip_numbers=True):
    """
    trees in the file correspond to the upper breakpoint position along the
    genome. If index_trees_by_variants == True then each tree name is instead
    the number of variants observed along the sequence so far (i.e. the
    uppermost variant index + 1). If index_tree_by_variants is a list, it is
    taken to be a list of variant positions to use as positions at which trees
    are output (with adjacent identical trees output as a single tree and
    indexed by the uppermost variant position, in the same manner as for
    index_trees_by_variants = True).

    Using the variant position to index trees allows trees to be compared for
    batches of variants rather than via chromosome position, allowing fair
    comparison between inferred trees.
    """
    print("#NEXUS\nBEGIN TREES;", file=treefile)
    increment = 0 if zero_based_tip_numbers else 1
    tip_map = [
        "{} {}".format(i + increment, i + increment)
        for i in range(ts.get_sample_size())]
    print("TRANSLATE\n{};".format(",\n".join(tip_map)), file=treefile)
    variant_index = 0
    if hasattr(index_trees_by_variants, "__len__"):
        # index_trees_by_variants contains an array of variant positions, which can
        # be used to index the trees
        assert min(index_trees_by_variants) >= 0, \
            "The variant positions passed in must all be greater than or equal to 0"
        assert max(index_trees_by_variants) < ts.get_sequence_length(), \
            "The variant positions passed in must all be less than {}".format(
                        ts.get_sequence_length())
        positions = np.sort(np.array(index_trees_by_variants))
    else:
        positions = None

    # TODO: should do ts.newick_trees(zero_based_tip_numbers)
    for t, (_, newick) in zip(ts.trees(), newick_trees(ts)):
        if not index_trees_by_variants:
            # index by rightmost genome position
            print("TREE " + str(t.get_interval()[1]) + " = " + newick, file=treefile)
        else:
            # index by rightmost variant number
            if positions is None:  # use the built-in variant positions
                n = t.get_num_mutations()
            else:
                l, r = t.get_interval()
                n = np.count_nonzero((positions >= l) & (positions < r))
            if n:
                # only print a tree if is is covered by at least 1 variant
                variant_index += n
                print("TREE " + str(variant_index) + " = " + newick, file=treefile)
    print("END;", file=treefile)
