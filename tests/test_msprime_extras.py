"""
Tests for the msprime extras module.
"""
import unittest
import tempfile

import dendropy

import msprime
import msprime_extras


class TestNewickTrees(unittest.TestCase):
    """
    Tests involving newick tree generation.
    """
    def get_binary_tree_sequence(self):
        return msprime.simulate(sample_size=10, recombination_rate=10, random_seed=1)

    def get_nonbinary_tree_sequence(self):
        events = [msprime.SimpleBottleneck(time=1, proportion=0.75)]
        ts = msprime.simulate(
            sample_size=10, recombination_rate=10, random_seed=1,
            demographic_events=events)
        non_binary = False
        for record in ts.records():
            if len(record.children) > 2:
                non_binary = True
        assert non_binary
        return ts

    def test_write_nexus_trees(self):
        ts = self.get_binary_tree_sequence()
        with tempfile.TemporaryFile("w+") as f:
            msprime_extras.write_nexus_trees(ts, f)
            f.seek(0)
            tree_list = dendropy.TreeList.get(file=f, schema="nexus")
            self.assertEqual(ts.num_trees, len(tree_list))
            # Hard to do more than this, as there is no sense of equality
            # in dendropy.

    def test_binary_newick(self):
        ts = self.get_binary_tree_sequence()
        l1 = list(ts.newick_trees())
        l2 = list(msprime_extras.newick_trees(ts))
        self.assertEqual(l1, l2)

    def test_nonbinary_newick(self):
        ts = self.get_nonbinary_tree_sequence()
        newick_trees = [newick for _, newick in msprime_extras.newick_trees(
            ts, precision=14, Ne=0.25)]
        self.assertEqual(len(newick_trees), ts.num_trees)
        for t, newick in zip(ts.trees(), newick_trees):
            tree = dendropy.Tree.get(data=newick, schema="newick")
            num_children_1 = [len(x.child_nodes()) for x in tree.preorder_node_iter()]
            num_children_2 = [len(t.children(u)) for u in t.nodes()]
            self.assertEqual(num_children_1, num_children_2)
            for msp_node, dendropy_node in zip(t.nodes(), tree.preorder_node_iter()):
                if t.parent(msp_node) != msprime.NULL_NODE:
                    branch_length_1 = dendropy_node.edge.length
                    branch_length_2 = t.branch_length(msp_node)
                    self.assertAlmostEqual(branch_length_1, branch_length_2)
