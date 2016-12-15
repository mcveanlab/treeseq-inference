#!/usr/bin/env python3

"""
Generate data to run fastARG on.
"""
import os

import numpy as np
import sys
sys.path.insert(0, '../msprime/') # look at the local version of ms
import msprime

def main(treefile_format = "nex"):
    sample_sizes = list(map(int, np.linspace(500, 5000, 10)))
    sample_sizes = [500] #just to test with a small sample size
    tmpdir = "../test_files"
    for n in sample_sizes:
        ts = msprime.simulate(
            sample_size=n, Ne=1e4, mutation_rate=2e-8, recombination_rate=2e-8,
            length=1e6)
        assert ts.get_sample_size() ==  n
        print("Finished simulating for {}; mutations = {}; trees = {}".format(
            n, ts.get_num_mutations(), ts.get_num_trees()))
        outfile = os.path.join(tmpdir, "{}.hdf5".format(n))
        ts.dump(outfile)
        if treefile_format.lower().startswith("nex"):
            with open(os.path.join(tmpdir, "{}.nex".format(n)), "w+") as treefile:
                print("#NEXUS\nBEGIN TREES;", file = treefile)
                print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i,i) for i in range(ts.get_sample_size())])), file = treefile)
                endpoint = 0
                for l,t in ts.newick_trees():
                    endpoint = endpoint+l
                    print("TREE " + str(endpoint) + " = " + t, file=treefile) 
                print("END;", file = treefile)
        else:
            with open(os.path.join(tmpdir, "{}.nwk".format(n)), "w+") as treefile:
                endpoint = 0
                for l,t in ts.newick_trees():
                    endpoint = endpoint+l
                    print(str(endpoint) + t, file=treefile) #put the treename first, to read into a multiphy object
        fastarg_file = os.path.join(tmpdir, "{}.haps".format(n))
        with open(fastarg_file, "w") as f:
            for j, v in enumerate(ts.variants(as_bytes=True)):
                print(j, v.genotypes.decode(), sep="\t", file=f)

if __name__ == "__main__":
    main()
