"""
Tests for the argmetrics code. Compare the results obtained using the
module with those from ETE and dendropy.
"""

import unittest
import tempfile

import msprime
import dendropy

import msprime_extras
import ARG_metrics


def get_single_binary_example(random_seed=1, length=1):
    return msprime.simulate(10, length=length, random_seed=random_seed)


def get_binary_example(random_seed=1, sample_size=10, length=10):
    ts = msprime.simulate(
        sample_size, recombination_rate=0.5, length=length, random_seed=random_seed)
    assert 5 < ts.num_trees < 100
    return ts


def get_single_nonbinary_example(random_seed=1):
    demographic_events = [
        msprime.SimpleBottleneck(time=0.5, proportion=0.9)]
    ts = msprime.simulate(
        20, random_seed=random_seed, demographic_events=demographic_events)
    nonbinary = False
    for e in ts.edgesets():
        if len(e.children) > 2:
            nonbinary = True
    assert nonbinary
    return ts


def get_nonbinary_example(random_seed=1, sample_size=20, length=10):
    demographic_events = [
        msprime.SimpleBottleneck(time=0.5, proportion=0.9)]
    ts = msprime.simulate(
        sample_size, random_seed=random_seed, demographic_events=demographic_events,
        recombination_rate=2)
    nonbinary = False
    for e in ts.edgesets():
        if len(e.children) > 2:
            nonbinary = True
    assert nonbinary
    assert 5 < ts.num_trees < 100
    return ts


class TestSingleTree(unittest.TestCase):
    """
    Tests on the simplest case of a single tree.
    """
    def verify_rf_metrics(self, ts1, ts2):
        """
        Verifies that RF metrics computed using dendropy and ARGmetrics for the
        specified tree sequence are equal, and returns the value.
        """
        with tempfile.NamedTemporaryFile("w+") as nexus_file_1, \
                tempfile.NamedTemporaryFile("w+") as nexus_file_2:
            msprime_extras.write_nexus_trees(
                ts1, nexus_file_1, zero_based_tip_numbers=False)
            msprime_extras.write_nexus_trees(
                ts2, nexus_file_2, zero_based_tip_numbers=False)
            nexus_file_1.flush()
            nexus_file_2.flush()
            nexus_file_1.seek(0)
            nexus_file_2.seek(0)
            tns = dendropy.TaxonNamespace()
            t1 = dendropy.Tree.get(
                file=nexus_file_1, schema="nexus", taxon_namespace=tns)
            t2 = dendropy.Tree.get(
                file=nexus_file_2, schema="nexus", taxon_namespace=tns)
            rf_distance_1 = dendropy.calculate.treecompare.symmetric_difference(t1, t2)
            metrics = ARG_metrics.get_metrics(nexus_file_1.name, nexus_file_2.name)
            rf_distance_2 = metrics["RFunrooted"]
        self.assertEqual(rf_distance_1, rf_distance_2)
        return rf_distance_1

    def test_same_binary_tree(self):
        ts1 = get_single_binary_example(1)
        ts2 = get_single_binary_example(1)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertEqual(d, 0)

    def test_binary_tree(self):
        ts1 = get_single_binary_example(1)
        ts2 = get_single_binary_example(2)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertGreater(d, 0)

    def test_same_nonbinary_tree(self):
        ts1 = get_single_nonbinary_example(1)
        ts2 = get_single_nonbinary_example(1)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertEqual(d, 0)

    def test_nonbinary_tree(self):
        ts1 = get_single_nonbinary_example(1)
        ts2 = get_single_nonbinary_example(2)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertGreater(d, 0)


class TestMultipleTrees(unittest.TestCase):
    """
    Tests on the situation where we have many trees.
    """

    def get_mean_rf_distance(self, ts1, ts2):
        """
        Returns the mean distance between the trees in the specified tree sequences.
        """
        assert ts1.sample_size == ts2.sample_size
        assert ts1.sequence_length == ts2.sequence_length
        trees1 = []
        intervals1 = []
        trees2 = []
        intervals2 = []
        tns = dendropy.TaxonNamespace()
        for t in ts1.trees():
            newick = msprime_extras.sparse_tree_to_newick(t, 1, 6)
            dt = dendropy.Tree.get(data=newick, schema="newick", taxon_namespace=tns)
            trees1.append(dt)
            intervals1.append(t.interval)
        assert len(trees1) == ts1.num_trees
        for t in ts2.trees():
            newick = msprime_extras.sparse_tree_to_newick(t, 1, 6)
            dt = dendropy.Tree.get(data=newick, schema="newick", taxon_namespace=tns)
            trees2.append(dt)
            intervals2.append(t.interval)
        assert len(trees2) == ts2.num_trees
        j1 = 0
        j2 = 0
        total_distance = 0
        total_metric = 0
        # I haven't tested this algorithm thoroughly, so there might be corner cases
        # not handled correctly. However, the total_distance assert below should
        # catch the problem if it occurs.
        while j1 < len(trees1) and j2 < len(trees2):
            # Each iteration of this loop considers one overlapping interval and
            # increments the counters.
            l1, r1 = intervals1[j1]
            l2, r2 = intervals2[j2]
            l = max(l1, l2)
            r = min(r1, r2)
            rf_distance = dendropy.calculate.treecompare.symmetric_difference(
                    trees1[j1], trees2[j2])
            total_metric += rf_distance * (r - l)
            total_distance += r - l
            if r1 <= r2:
                j1 += 1
            if r1 >= r2:
                j2 += 1
        self.assertAlmostEqual(total_distance, ts1.sequence_length)
        return total_metric / total_distance

    def verify_rf_metrics(self, ts1, ts2):
        """
        Verifies that RF metrics computed using dendropy and ARGmetrics for the
        specified tree sequence are equal, and returns the value.
        """
        with tempfile.NamedTemporaryFile("w+") as nexus_file_1, \
                tempfile.NamedTemporaryFile("w+") as nexus_file_2:
            msprime_extras.write_nexus_trees(
                ts1, nexus_file_1, zero_based_tip_numbers=False)
            msprime_extras.write_nexus_trees(
                ts2, nexus_file_2, zero_based_tip_numbers=False)
            nexus_file_1.flush()
            nexus_file_2.flush()
            metrics = ARG_metrics.get_metrics(nexus_file_1.name, nexus_file_2.name)
            rf_distance_2 = metrics["RFunrooted"]
        rf_distance_1 = self.get_mean_rf_distance(ts1, ts2)
        self.assertAlmostEqual(rf_distance_1, rf_distance_2)
        return rf_distance_1

    def test_binary_trees(self):
        ts1 = get_binary_example(1)
        ts2 = get_binary_example(3)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertGreater(d, 0)

    def test_same_tree_sequence(self):
        ts1 = get_binary_example(1)
        ts2 = get_binary_example(1)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertEqual(d, 0)

    def test_single_tree_with_sequence(self):
        ts1 = get_binary_example(1, length=10)
        ts2 = get_single_binary_example(3, length=10)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertGreater(d, 0)

    def test_two_single_trees(self):
        ts1 = get_single_binary_example(5, length=10)
        ts2 = get_single_binary_example(3, length=10)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertGreater(d, 0)

    def test_nonbinary_trees(self):
        ts1 = get_nonbinary_example(1)
        ts2 = get_nonbinary_example(3)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertGreater(d, 0)

    def test_binary_nonbinary_trees(self):
        ts1 = get_binary_example(1, sample_size=20)
        ts2 = get_nonbinary_example(3, sample_size=20)
        d = self.verify_rf_metrics(ts1, ts2)
        self.assertGreater(d, 0)
