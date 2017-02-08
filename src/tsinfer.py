"""
Implementation of the Li and Stephens algorithm for inferring a tree
sequence.

Python 3 only.
"""

import sys
import os
import time
import collections
import threading
import tempfile
import resource
import subprocess
import multiprocessing

import numpy as np

import msprime
import _msprime
import _tsinfer

if sys.version_info[0] < 3:
    raise Exception("Python 3 you idiot!")


class ReferencePanel(object):
    """
    Class representing the reference panel for inferring a tree sequence
    from observed data.
    """
    def __init__(self, samples, sites, sequence_length):
        self._ll_reference_panel = _tsinfer.ReferencePanel(samples)
        self.sites = sites
        self.sequence_length = sequence_length
        assert len(self.sites) == self._ll_reference_panel.num_sites

    def infer_paths(self, rho, num_workers=None):
        N = self._ll_reference_panel.num_haplotypes
        n = self._ll_reference_panel.num_samples
        m = self._ll_reference_panel.num_sites
        P = np.zeros((N, m), dtype=np.uint32)
        P[-1,:] = -1

        work = []
        for j in range(N - 2, n - 1, -1):
            work.append((j, N - j - 1))
        for j in range(n):
            work.append((j, N - n))
        work_index = 0
        lock = threading.Lock()

        def worker():
            nonlocal work_index
            threader = _tsinfer.Threader(self._ll_reference_panel)
            while True:
                with lock:
                    if work_index >= len(work):
                        break
                    haplotype_index, panel_size = work[work_index]
                    work_index += 1
                threader.run(haplotype_index, panel_size, rho, P[haplotype_index])

        threads = []
        if num_workers is None:
            num_workers = multiprocessing.cpu_count()
        for _ in range(num_workers):
            t = threading.Thread(target=worker)
            t.start()
            threads.append(t)
        for t in threads:
            t.join()
        return P

    def convert_records(self, P):
        N = self._ll_reference_panel.num_haplotypes
        n = self._ll_reference_panel.num_samples
        m = self._ll_reference_panel.num_sites
        H = self._ll_reference_panel.get_haplotypes()
        sites = self.sites

        assert N, m == P.shape
        assert H.shape == P.shape
        C = [[[] for _ in range(m)] for _ in range(N)]

        for j in range(N - 1):
            for l in range(m):
                C[P[j][l]][l].append(j)

        # print("children = ")
        # for j in range(N):
        #     print(j, "\t", C[j])

        mutations = []
        records = []
        # Create the mutations by finding the oldest 1 in each locus.
        # print("mutations = ")
        for l in range(m):
            u = np.where(H[:,l] == 1)[0][0]
            # u is a sample with this mutations. Follow its path upwards until
            # we find the oldest node with the mutation.
            while H[u, l] == 1:
                v = u
                u = P[u][l]
            mutations.append((sites[l], (v,)))
        assert len(mutations) == m
        for u in range(n, N):
            row = C[u]
            last_c = row[0]
            left = 0
            for l in range(1, m):
                if row[l] != last_c:
                    if len(last_c) > 0:
                        records.append(msprime.CoalescenceRecord(
                            left=sites[left], right=sites[l], node=u, children=tuple(last_c),
                            time=u, population=0))
                    left = l
                    last_c = row[l]
            if len(last_c) > 0:
                records.append(msprime.CoalescenceRecord(
                    left=sites[left], right=self.sequence_length, node=u,
                    children=tuple(last_c), time=u, population=0))
        records.sort(key=lambda r: r.time)
        # print("records = ")
        # for r in records:
        #     print(r)
        # for m in mutations:
        #     print(m)
        samples = [(0, 0) for _ in range(n)]
        ll_ts = _msprime.TreeSequence()
        ll_ts.load_records(
            samples=samples, coalescence_records=records, mutations=mutations)
        assert ll_ts.get_sample_size() == n
        assert ll_ts.get_num_mutations() == m
        ts = msprime.TreeSequence(ll_ts)
        assert ts.num_mutations == m
        assert len(list(ts.mutations())) == m
        return ts


def check_paths(H, P, n):
    """
    Checks that the specified paths through the haplotypes make sense.
    """
    N, m = P.shape
    assert H.shape == P.shape
    i = np.arange(m)
    for j in range(N - 1):
        p = P[j]
        h = H[j]
        # hp = np.array([H[p[l], l] for l in range(n)])
        hp = H[p, i]
        # print("path:")
        # print(p)
        # print(h)
        # print(hp)
        for l in np.where(hp != h)[0]:
            assert hp[l] == 0

