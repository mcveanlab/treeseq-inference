"""
Extra functions for tsinfer, in particular, for resolving polytomies

We might want to merge this at some point into the existing tsinfer file

Python 3 only.
"""

import sys
import os
import time
import collections
import threading
import resource
import logging
import itertools

import numpy as np

import msprime
import tsinfer


def get_nodes_and_mutations(ts):
    """
    Extract the nodes & edgesets from a treesequence.
    Will need recoding when JK changes the calling convention
    """
    nodes = msprime.NodeTable()
    mutations = msprime.MutationsTable()
    ts.dump_tables(nodes=nodes, mutations=mutations) 
    return nodes, mutations

def resolve_polytomies(ts, polytomy_func):
    """
    polytomy_func should take a set of edge records, and an edgesets and a nodes object
    to be added to.
    """
    new_edgesets = msprime.EdgesetTable()
    nodes, mutations = get_nodes_and_mutations(ts)
    edge_records = [[]] #store the edge records per parent, split into contiguous blocks
    for e in ts.edgesets():    #assume records are in order
        if len(edge_records[0]==0) or e.parent == records[0][0].parent:
            if e.right==edge_records[-1][-1].left:
                #contiguous with the last record
                edge_records[-1].append(e)
            else:
                #this is the same parent, but not contiguous
                edge_records.append([e])
        else:
            #submit records for polytomy resolution - may require new nodes to be created
            polytomy_func(edge_records, new_edgesets, nodes)
            edge_records = [[e]]
    if edge_records:
        #last loop
        polytomy_func(edge_records, nodes, new_edgeset)

    return msprime.load_tables(nodes=nodes, edgesets=new_edgesets, mutations=mutations)


def runs_of_ones(bits):
  for bit, group in itertools.groupby(bits):
    if bit: yield sum(group)

def resolve_polytomy_jk(old_edgesets_with_shared_parent, new_edgesets, nodes):
    """
    JK version - break into contigous lengths. Then for each contigous length, if any in the 
    old_edgesets contain >2 children, find a subset of children that is common to the largest span
    and create new node BELOW the parent N, which groups these together.
    
    Since we are doing this over all edgesets that share the same parent, we need not worry 
    that the new node is created below N: it will be shared by all the relevant sections of genome.
    
    The advantage to creating a node below N is that elsewhere on the tree, the links to N will remain
    valid (e.g. placement of mutations, N as a child of another node, etc)
    
    Question: should 'largest span' be measured in genome units or in number of C records shared. If we 
    measure in C records, they might (approximately) correspond to the number of recombination events
    that have occurred elsewhere but have *not* separated these two pairs. However, because of the 
    simplify routines
    suggests the latter, since this corresponds to the number of recombination events. But we can try both
    
    """
    for contig in old_edgesets_with_shared_parent:
        pairlengths = {} 
        #for each possible pair, make an array of 0/1s such that array[i] indicates
        # whether this pair is present in edgeset i.        
        for i, edgeset in enumerate(contig):
            for pairs in itertools.combinations(edgeset.children)
                if pairs not in pairlengths:
                    pairlengths[pairs] = numpy.zeroes(len(contig), dtype=bool)
                pairlengths[pairs][i] = True
        for pairs, bitarray in pairlengths.items()
            max(runs_of_ones(bitarray))
    
        #make a new node joining the pair with the longests shared span
        curr_node = old_edgesets_with_shared_parent[0][0].parent
        curr_time = nodes[curr_node].time
        nodes.add_row(time=curr_time-epsilon, population=2)
        new_edgesets.add_row(left=0, right=1, parent=4, children=(0, 1, 3))
        longest_pair[0]
    
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Test polytomy resolving')
    parser.add_argument('--sample_size', '-n', type=int, default=5, help='the sample size if an .trees file is not given')
    parser.add_argument('--effective_population_size', '-Ne', type=float, default=5000, help='the effective population size')
    parser.add_argument('--sequence_length', '-l', type=float, default=55000, help='the sequence length')
    parser.add_argument('--recombination_rate', '-rho', type=float, default=2.5e-8, help='the recombination rate')
    parser.add_argument('--mutation_rate', '-mu', type=float, default=5e-8, help='the mutation rate')
    parser.add_argument('--random_seed', '-seed', type=int, default=1234, help='a random seed for msprime & AW simulation')
    parser.add_argument('--verbosity', '-v', action='count', default=0)
    args = parser.parse_args()
    log_level = logging.WARNING
    if args.verbosity == 1:
        log_level = logging.INFO
    if args.verbosity >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(
        format='%(asctime)s %(message)s', level=log_level, stream=sys.stdout)
    ts = msprime.simulate(
            sample_size = args.sample_size, 
            Ne=args.effective_population_size, 
            length=args.sequence_length,
            recombination_rate=args.recombination_rate, 
            mutation_rate=args.mutation_rate, 
            random_seed=args.random_seed)
    resolve_jk(ts, args.random_seed, args.mutation_rate*2)