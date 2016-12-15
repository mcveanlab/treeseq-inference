#!/usr/bin/env python3

"""
Generate treefiles to test treecmp
"""
import os
import sys
import numpy as np
import os.path
sys.path.insert(1,os.path.join(sys.path[0],'..','fastARG')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(sys.path[0],'..','msprime')) # use the local copy of msprime in preference to the global one
import msprime
from warnings import warn
from tempfile import NamedTemporaryFile
from msprime_fastARG import *

def print_nexus(ts, treefile):
    import sys
    
    print("#NEXUS\nBEGIN TREES;", file = treefile)
    print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i,i) for i in range(ts.get_sample_size())])), file = treefile)
    breakpoint = 0
    trees = 0
    epsilon = 1e-8
    for l,t in ts.newick_trees():
        trees += 1
        breakpoint = breakpoint +l
        print("TREE " + str(breakpoint) + " = " + t, file=treefile)         
    print("END;", file = treefile)
    if treefile.name != "<stdout>":
        print("output {} trees".format(trees))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generate treefiles to test treecmp.r')
    args = parser.parse_args()
    import csv
    import time
    import random
    import filecmp
    seed_generator = random.Random()
    seed_generator.seed(123)
    n_simulation_replicates = 1
    n_fastARG_replicates = 1
    sample_size_sims = [{'sample_size':ss,'length':int(1e6)} for ss in map(int, np.linspace(500, 5000, 10))]
    sample_size_sims = [{'sample_size':500,'length':int(1e6)}]
    length_sims = [{'sample_size':1000,'length':l} for l in map(int, np.linspace(1e6, 1e8, 10))]
    length_sims = []
    for params in sample_size_sims + length_sims:
        for sim_replicates in range(n_simulation_replicates):
            sim_seed = seed_generator.randint(1,2**32-1)
            n = int(params['sample_size'])
            l = int(params['length'])
            print("Simulating n_haps={} seq_len={} (seed is {})".format(n, l, sim_seed))
            
            ts = msprime.simulate(
                sample_size=n, Ne=1e4, mutation_rate=2e-8, recombination_rate=2e-8,
                length=l, random_seed=sim_seed)
            assert ts.get_sample_size() ==  n
            with open("{}-{}.nex".format(n,l), "w+") as nex:
                print_nexus(ts, nex)
            print("Finished simulating for {}; mutations = {}; trees = {}".format(n, ts.get_num_mutations(), ts.get_num_trees()))
            for inference_replicates in range(n_fastARG_replicates):
                with NamedTemporaryFile("w+") as fa_in, \
                     NamedTemporaryFile("w+") as fa_out, \
                     NamedTemporaryFile("w+") as tree, \
                     NamedTemporaryFile("w+") as fa_revised, \
                     NamedTemporaryFile("w+") as muts:
                    msprime_to_fastARG_in(ts, fa_in)
                    inf_seed = seed_generator.randint(1,2**32-1)
                    start_time = time.time()
                    run_fastARG("../fastARG/fastARG", fa_in, fa_out, seed=inf_seed)
                    end_time = time.time()
                    root_seq = fastARG_root_seq(fa_out)
                    fastARG_out_to_msprime_txts(fa_out, variant_positions_from_fastARGin(fa_in), tree, muts, l)
                    ts_new = msprime_txts_to_fastARG_in_revised(tree, muts, root_seq, fa_revised)
                    #quick check
                    if filecmp.cmp(fa_in.name, fa_revised.name, shallow=False) == False:
                        warn("Initial fastARG input file differs from processed fastARG file")
                    with open("{}-{}.fa.nex".format(n,l), "w+") as nex:
                        print_nexus(ts_new, nex)
