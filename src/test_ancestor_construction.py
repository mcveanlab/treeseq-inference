#!/usr/bin/env python3.5
description = """
Run some simulations with different methods of ancestor reconstruction and output how 
well we do in temrs of number of trees, edges, and overall tree metrics.
"""
import sys
import os
import re
import tempfile
import subprocess
import random

import numpy as np
import pandas as pd

sys.path.insert(1,os.path.join(sys.path[0],'..','msprime'))
sys.path.insert(1,os.path.join(sys.path[0],'..','tsinfer'))
import msprime
import msprime_extras
import tsinfer
import ARG_metrics
from dev import single_real_ancestor_injection


def print_treeseq(ts, positions):
    for t in ts.trees():
        interval = t.get_interval()
        tree_variants = np.logical_and(interval[0] <= positions,  positions < interval[1])
        span = "".join(["-" if pos else " " for pos in tree_variants])
        print(re.sub(r"-(?= |$)", "<", re.sub(r"((?<= )|^)-", '>', span)))
        tree = t.draw(format="unicode")
        pos = np.where(tree_variants)[0][0]
        print("\n".join([" "*(pos) + line for line in tree.splitlines()]))


def main(verbosity, only_one_method):
    #monkey patch a nexus-writing routine in here, so we can test tree metrics
    msprime.TreeSequence.write_nexus_trees = msprime_extras.write_nexus_trees
    method_range = slice(0,4) if only_one_method==None else slice(only_one_method, only_one_method+1)

    method = "C"
    print_trees=True if verbosity else False
    path_compression = True
    msprime_args = dict(sample_size=4, recombination_rate=0.35, model="smc_prime")
    rng1 = random.Random(123)
    rng2 = random.Random(12)
    if not print_trees:
        print("seed", "trees","sites","edges:", sep="\t", end="\t")
        print("\t".join(["tsinfer_norm", "known_anc_orig", "known_anc_jerome", "known_anc_yan"][method_range]))

    #msprime_args.update(dict(random_seed = 58547, sample_size = 6))
    #msprime_args.update(dict(random_seed = 23900, sample_size = 5))
    msprime_args.update(dict(random_seed = 58763, sample_size = 4, model="smc_prime"))

    reps = 1 if 'random_seed' in msprime_args else 10000
    for i in range(reps):
        if 'random_seed' not in msprime_args or i>0:
            msprime_args['random_seed'] = rng1.randint(1, 100000)
        ts, full_inferred_ts, orig_anc_ts, jk_anc_ts, hyw_anc_ts = \
            single_real_ancestor_injection(method, path_compression, simplify=True, **msprime_args)
        #for e in ts.edges():
        #    print(e)
        #for m in ts.mutations():
        #    print(m)

        #print(ts.num_trees, ts.num_sites, ts.num_edges, sep="\t", end="\t")
        positions = [v.position for v in ts.variants()]

        if print_trees:
            print(list(enumerate(positions)))
            for h in ts.haplotypes():
                print(h)
            print("Real tree sequence")
            print_treeseq(ts, np.array(positions))
        with tempfile.NamedTemporaryFile("w+") as original_nexus:
            ts.write_nexus_trees(original_nexus)
            original_nexus.flush()
            for i, i_ts in list(enumerate((full_inferred_ts, orig_anc_ts, jk_anc_ts, hyw_anc_ts)))[method_range]:
                with tempfile.NamedTemporaryFile("w+") as inferred_nexus:
                    i_ts.write_nexus_trees(inferred_nexus)
                    inferred_nexus.flush()
                    metrics = []
                    for polytomy_reps in range(1):
                        np.set_printoptions(precision=12)
                        #print(ARG_metrics.get_full_metrics(original_nexus.name, inferred_nexus.name, variant_positions=positions))
                        metrics.append(ARG_metrics.get_metrics(
                            original_nexus.name, inferred_nexus.name,
                            #randomly_resolve_inferred = rng2.randint(1, 2**31),
                            variant_positions=positions,
                            ))
                    metrics = pd.DataFrame(metrics).mean().to_dict()
                    
                    setattr(i_ts,"treestat", metrics['KCrooted'])
                    subprocess.call(["cp", original_nexus.name, "./orig.nex"])
                    subprocess.call(["cp", inferred_nexus.name, "./inferred{}.nex".format(i)])
                    if print_trees:
                        print("Ancestor construction method " + str(i), metrics)
                        print_treeseq(i_ts, np.array(positions))
                        print("_"*80)
        print(msprime_args['random_seed'], ts.num_trees,ts.num_sites,ts.num_edges, sep="\t", end="\t")
        inferred = np.array([[x.num_edges, x.num_trees, x.treestat] \
            for x in (full_inferred_ts, orig_anc_ts, jk_anc_ts, hyw_anc_ts)[method_range]])
        print("\t".join(["{:>2.0f} {:>2.0f}/{:>2.0f} {:.4f}".format(x[1],x[0], ts.num_edges, x[2]) + \
            ("*" if x[2]==min(inferred[:,2]) and np.sum(np.isclose(inferred[:,2], min(inferred[:,2]), atol=1e-6))==1 else " ")\
            for x in inferred]))
        if hyw_anc_ts.treestat > 0.0005 and hyw_anc_ts.num_trees < 5:
            sys.exit()


        
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=description, add_help=False)
    parser.add_argument("-v","--verbosity", action="count",
            help="Verbosity level", default=0)    
    parser.add_argument("-a","--ancestor_method", default=None, type=int, choices=range(4),
            help="Only try one method of ancestor reconstruction")
    args = parser.parse_args()

    main(args.verbosity, args.ancestor_method)
