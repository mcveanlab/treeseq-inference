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


def get_single_binary_example(random_seed=1):
    return msprime.simulate(10, random_seed=random_seed)


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


class TestSingleTree(unittest.TestCase):
    """
    Tests on the simplest case of a single tree.
    """
    def setUp(self):
        self.nexus_file_1 = tempfile.NamedTemporaryFile("w+")
        self.nexus_file_2 = tempfile.NamedTemporaryFile("w+")

    def tearDown(self):
        self.nexus_file_1.close()
        self.nexus_file_2.close()

    def verify_rf_metrics(self, ts1, ts2):
        """
        Verifies that RF metrics computed using dendropy and ARGmetrics for the
        specified tree sequence are equal, and returns the value.
        """
        msprime_extras.write_nexus_trees(
            ts1, self.nexus_file_1, zero_based_tip_numbers=False)
        msprime_extras.write_nexus_trees(
            ts2, self.nexus_file_2, zero_based_tip_numbers=False)
        self.nexus_file_1.flush()
        self.nexus_file_2.flush()
        self.nexus_file_1.seek(0)
        self.nexus_file_2.seek(0)
        tns = dendropy.TaxonNamespace()
        t1 = dendropy.Tree.get(
            file=self.nexus_file_1, schema="nexus", taxon_namespace=tns)
        t2 = dendropy.Tree.get(
            file=self.nexus_file_2, schema="nexus", taxon_namespace=tns)
        rf_distance_1 = dendropy.calculate.treecompare.symmetric_difference(t1, t2)
        metrics = ARG_metrics.get_metrics(self.nexus_file_1.name, self.nexus_file_2.name)
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
