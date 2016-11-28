#!/usr/bin/env python3

"""
Generate data to run fastARG on.
"""
import os

import numpy as np
import sys
sys.path.insert(0, '../msprime/') # look at the local version of ms
import msprime

def main():
    sample_sizes = list(map(int, np.linspace(500, 5000, 10)))
    tmpdir = "./"
    for n in sample_sizes:
        ts = msprime.simulate(
            sample_size=n, Ne=1e4, mutation_rate=2e-8, recombination_rate=2e-8,
            length=1e6)
        assert ts.get_sample_size() ==  n
        print("Finished simulating for {}; mutations = {}; trees = {}".format(
            n, ts.get_num_mutations(), ts.get_num_trees()))
        outfile = os.path.join(tmpdir, "{}.hdf5".format(n))
        ts.dump(outfile)

        fastarg_file = os.path.join(tmpdir, "{}.haps".format(n))
        with open(fastarg_file, "w") as f:
            for j, v in enumerate(ts.variants(as_bytes=True)):
                print(j, v.genotypes.decode(), sep="\t", file=f)

if __name__ == "__main__":
    main()
