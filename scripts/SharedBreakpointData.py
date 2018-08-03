#!/usr/bin/env python3
"""
Simply output the number of shared breakpoints seen in a tree sequence
"""
import os
import argparse
import collections

import numpy as np
import msprime

def main():
    parser = argparse.ArgumentParser(
        description=__doc__)
    parser.add_argument(
        "tree_sequence_file",
        help="The input tree sequence file.")
    args = parser.parse_args()


    ts = msprime.load(args.tree_sequence_file)
    srb = collections.defaultdict(int)
    for (left, right), edges_out, edges_in in ts.edge_diffs():
        if len(edges_out) and len(edges_in):
            child_data = collections.defaultdict(list)
            for edge in edges_out:
                child_data[edge.child].append(edge.parent)
            for edge in edges_in:
                child_data[edge.child].append(edge.parent)
            for c, parents in child_data.items():
                assert 0 < len(parents) < 3
                if len(parents) == 2:
                    srb[(left, parents[0],parents[1])] += 1
            #for k,v in srb.items():
            #    if v>1:
            #        for c, ps in child_data.items():
            #            #print(tuple(ps), k[1:])
            #            #print(tuple(ps), k[1:])
            #            if tuple(ps) == k[1:]:
            #                print(k[0], c, ps)

    count = np.bincount(np.array(list(srb.values()), dtype=np.int))
    
    print("Sharing counts for recombination breakpoints"
    print(np.array([count[1:], np.arange(1,len(count))]).transpose())

if __name__ == "__main__":
    main()
