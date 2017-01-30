"""
One-off test to see if subsetting works.
"""
import sys
import os
import tempfile
import subprocess

import numpy as np

curr_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(1,os.path.join(curr_dir,'..','msprime'))
import msprime
import msprime_extras
import tsinfer

tsinfer_executable = os.path.join(curr_dir,'run_tsinfer.py')

msprime.TreeSequence.write_nexus_trees = msprime_extras.write_nexus_trees
#R tree metrics assume tips are numbered from 1 not 0
tree_tip_labels_start_at_0 = False

def make_errors(v, p):
    """
    Flip each bit with probability p.
    """
    m = v.shape[0]
    mask = np.random.random(m) < p
    return np.logical_xor(v.astype(bool), mask).astype(int)

def generate_samples(ts, error_p):
    """
    Returns samples with a bits flipped with a specified probability.

    Rejects any variants that result in a fixed column.
    """
    S = np.zeros((ts.sample_size, ts.num_mutations), dtype="u1")
    for variant in ts.variants():
        done = False
        # Reject any columns that have no 1s or no zeros
        while not done:
            S[:,variant.index] = make_errors(variant.genotypes, error_p)
            s = np.sum(S[:, variant.index])
            done = 0 < s < ts.sample_size
    return S


def main():
    length=5000
    Ne= 5000
    rho = 1.5e-8
    mu = 2.5e-8
    ts = msprime.simulate(sample_size= 10000, Ne=Ne, length=length, recombination_rate=rho, mutation_rate=mu, random_seed= 224899943)
    small_ts = ts.simplify(list(range(10)))
    print("simulated with {} mutations".format(ts.get_num_mutations()))
    with open("tmp/sim.nex", "w+") as out:
        #tree metrics assume tips are numbered from 1 not 0
        small_ts.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
    
    for subsample in [10,20,50,100,200,500,1000, 5000, 10000]:
        ts_sub = ts.simplify(list(range(subsample)))
        S = generate_samples(ts_sub, 0)
        sample_fn="tmp/sub{}".format(subsample) + ".npy"
        positions_fn="tmp/sub{}".format(subsample) + ".pos.npy"
        nex_fn="tmp/sub{}".format(subsample) + ".nex"
        tmp_fn="tmp/tmp.hdf5"
        np.save(sample_fn, S)
        np.save(positions_fn, np.array([v.position for v in ts_sub.variants()]))
        cmd = ["python3", tsinfer_executable, sample_fn, positions_fn, tmp_fn,
            "--length", str(int(length)), "--recombination-rate", str(rho)]
        print("running {}".format(" ".join(cmd)))
        subprocess.call(cmd)
        ts_simplified = msprime.load(tmp_fn)
        ts_inferred = ts_simplified.simplify(list(range(10)))
        with open(nex_fn, "w+") as out:
            #tree metrics assume tips are numbered from 1 not 0
            ts_inferred.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
        
if __name__ == "__main__":
    main()
