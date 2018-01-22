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


def main(verbosity):
    #monkey patch a nexus-writing routine in here, so we can test tree metrics
    msprime.TreeSequence.write_nexus_trees = msprime_extras.write_nexus_trees

    method = "C"
    print_trees=True if verbosity else False
    path_compression = True
    rng1 = random.Random(1234)
    rng2 = random.Random(12)
    if not print_trees:
        print("seed", "trees","sites","edges:", "tsinfer_norm", "known_anc_orig", 
            "known_anc_jerome", "known_anc_yan", sep="\t")
    for i in range(100):
        r_int = rng1.randint(1, 100000)
        #r_int = 5190
        ts, full_inferred_ts, orig_anc_ts, jk_anc_ts, hyw_anc_ts = \
            single_real_ancestor_injection(method, path_compression, r_int, simplify=True)

        #print(ts.num_trees, ts.num_sites, ts.num_edges, sep="\t", end="\t")
        positions = [v.position for v in ts.variants()]

        if print_trees:
            for h in ts.haplotypes():
                print(h)
            print("Real tree sequence")
            print_treeseq(ts, np.array(positions))
        with tempfile.NamedTemporaryFile("w+") as original_nexus:
            ts.write_nexus_trees(original_nexus, tree_labels_between_variants=True)
            original_nexus.flush()
            for i, i_ts in enumerate((full_inferred_ts, orig_anc_ts, jk_anc_ts, hyw_anc_ts)):
                with tempfile.NamedTemporaryFile("w+") as inferred_nexus:
                    i_ts.write_nexus_trees(inferred_nexus, tree_labels_between_variants=True)
                    inferred_nexus.flush()
                    metrics = []
                    for polytomy_reps in range(1):
                        metrics.append(ARG_metrics.get_metrics(
                            original_nexus.name, inferred_nexus.name, variant_positions = positions,
                            #randomly_resolve_inferred = rng2.randint(1, 2**31)
                            ))
                    metrics = pd.DataFrame(metrics).mean().to_dict()
                    
                    setattr(i_ts,"treestat", metrics['KCrooted'])
                    subprocess.call(["cp", original_nexus.name, "./orig.nex"])
                    subprocess.call(["cp", inferred_nexus.name, "./inferred{}.nex".format(i)])
                    if print_trees:
                        print("Ancestor construction method " + str(i), metrics)
                        print_treeseq(i_ts, np.array(positions))
                        print("_"*80)
        print(r_int, ts.num_trees,ts.num_sites,ts.num_edges, sep="\t", end="\t")
        inferred = np.array([[x.num_edges, x.num_trees, x.treestat] \
            for x in (full_inferred_ts, orig_anc_ts, jk_anc_ts, hyw_anc_ts)])
        print("\t".join(["{:>2.0f} {:>2.0f}/{:>2.0f} {:.2f}".format(x[1],x[0], ts.num_edges, x[2]) + \
            ("*" if x[2]==min(inferred[:,2]) and np.sum(inferred[:,2]==min(inferred[:,2]))==1 else " ")\
            for x in inferred]))


        
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=description, add_help=False)
    parser.add_argument("-v","--verbosity", action="count",
            help="Verbosity level", default=0)    
    args = parser.parse_args()

    main(args.verbosity)
