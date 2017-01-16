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

class TestDiscretisePositions(unittest.TestCase):
    """
    Tests for discretising the positions of mutations so that we can
    convert easily to other formats.
    """
    def verify_discretise(self, ts):
        ts_new = msprime_extras.discretise_mutations(ts)
        self.assertEqual(list(ts.records()), list(ts_new.records()))
        breakpoints = sorted(set([r.left for r in ts.records()] + [ts.sequence_length]))
        for j in range(len(breakpoints) - 1):
            left, right = breakpoints[j], breakpoints[j + 1]
            old_mutations = [
                mut for mut in ts.mutations() if left <= mut.position < right]
            new_mutations = [
                mut for mut in ts_new.mutations() if left <= mut.position < right]
            self.assertEqual(len(old_mutations), len(new_mutations))
            self.assertEqual(
                [mut.node for mut in old_mutations],
                [mut.node for mut in new_mutations])
            self.assertEqual(
                [mut.index for mut in old_mutations],
                [mut.index for mut in new_mutations])
            self.assertEqual(
                len(new_mutations),
                len(set(mut.position for mut in new_mutations)))
            for mut in new_mutations:
                self.assertTrue(int(mut.position) == mut.position)

    def test_no_recombination(self):
        ts = msprime.simulate(10, length=100, mutation_rate=0.05, random_seed=1)
        self.assertGreater(ts.num_mutations, 5)
        self.verify_discretise(ts)

    def test_nominal_case(self):
        recomb_map = msprime.RecombinationMap.uniform_map(
            length=3000, rate=0.005, num_loci=3000)
        ts = msprime.simulate(
            10, recombination_map=recomb_map, mutation_rate=0.005, random_seed=1)
        self.assertGreater(ts.num_records, 20)
        self.assertGreater(ts.num_mutations, 10)
        self.assertLess(ts.num_mutations, 500)
        self.verify_discretise(ts)

    def test_unit_interval(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        self.assertGreater(ts.num_mutations, 1)
        self.assertRaises(ValueError, msprime_extras.discretise_mutations, ts)

    def test_dense_mutations(self):
        ts = msprime.simulate(10, length=100, mutation_rate=0.5, random_seed=1)
        self.assertGreater(ts.num_mutations, 10)
        self.assertRaises(ValueError, msprime_extras.discretise_mutations, ts)
