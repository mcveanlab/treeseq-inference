#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert fastARG format to msprime hdf5.

The node numbers output by fastARG seem to be strictly increasing, so can simply be used as node ages within msprime

"""    

import os
import sys
import re
import numpy as np
import csv
from warnings import warn

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert fastARG format to msprime hdf5')
    parser.add_argument('infile', type=argparse.FileType('r', encoding='UTF-8'), help='an output file from fastARG (https://github.com/lh3/fastARG )')
    parser.add_argument('outfile', help='The name (without suffix) of an output file, or "-" for stdout. If saving to a file, the suffix ".msprime" will be added to create the tree output file (format explained at http://msprime.readthedocs.io/en/stable/api.html#msprime.load_txt), and the suffix ".msmut" will be added to create the mutations file')
    args = parser.parse_args()


    coalescence_points={}
    #intermediate storage format - allows fastARG records to be split up
    #e.g. cr[X] = child1:[left,right], child2:[left,right],...

    mutations={}
    mutation_nodes=set() #this is just to check that there aren't duplicate nodes with the same mutations - a fastARG restriction    
    
    for line_num, fields in enumerate(csv.reader(args.infile, delimiter='\t')):
        if   fields[0]=='E':
            srand48seed = int(fields[1])
            
        elif fields[0]=='N':
            haplotypes, loci = int(fields[1]), int(fields[2])
            
        elif fields[0]=='C' or fields[0]=='R':
            #coalescence or recombination events - 
            #coalescence events should come in pairs
            #format is curr_node, child_node, seq_start_inclusive, seq_end_noninclusive, num_mutations, mut_loc1, mut_loc2, ....
            curr_node, child_node, left, right, n_mutations = [int(i) for i in fields[1:6]]
            if curr_node<child_node:
                warn("Line {} has a ancestor node_id less than the id of its children, so node ID cannot be used as a proxy for age".format(line_num))
                sys.exit()
            #one record for each node
            if curr_node not in coalescence_points:
                coalescence_points[curr_node] = {}
            if child_node in coalescence_points[curr_node]:
                warn("Child node {} already exists for node {} (line {})".format(child_node, curr_node, line_num))
            coalescence_points[curr_node][child_node]=[left, right]
            #mutations in msprime are placed as ancestral to a target node, rather than descending from a child node (as in fastARG)
            #we should check that (a) mutation at the same locus does not occur twice 
            # (b) the same child node is not used more than once (NB this should be a fastARG restriction)
            if n_mutations:
                if child_node in mutation_nodes:
                    warn("Node {} already has some mutations: some more are being added from line {}.".format(child_node, line_num))
                mutation_nodes.add(child_node)
                for m in [float(pos) for pos in fields[6:(6+n_mutations)]]:
                    if m in mutations:
                        warn("Duplicate mutations at a single locus: {}. One is from node {}.".format(pos, child_node))
                    else:
                        mutations[m]=child_node
                

        elif fields[0]=='S':
            #sequence at root
            root_node = int(fields[1])
            base_sequence = fields[2]
            
        else:
            warn("Bad line format - the fastARG file has a line that does not begin with E, N, C, R, or S:\n{}".format("\t".join(fields)))
    
    #Done reading in
    
    if min(coalescence_points.keys()) != haplotypes:
            warn("The smallest internal node id is expected to be the same as the number of haplotypes, but it is not ({} vs {})".format(min(coalescence_points.keys()), haplotypes))
    
    tree_filehandle = sys.stdout if args.outfile=="-" else open(args.outfile + ".msprime", 'w')
    mutations_filehandle = sys.stdout if args.outfile=="-" else open(args.outfile + ".msmut", 'w')

    for node,children in sorted(coalescence_points.items()): #sort by node number == time
        #look at the break points for all the child sequences, and break up into that number of records
        breaks = set()
        for leftright in children.values():
            breaks.update(leftright)
        breaks = sorted(breaks)
        for i in range(1,len(breaks)):
            leftbreak = breaks[i-1]    
            rightbreak = breaks[i]
            #NB - here we could try to input the values straight into an msprime python structure,
            #but until this feature is implemented in msprime, we simply output to a correctly formatted text file
            tree_filehandle.write("{}\t{}\t{}\t".format(leftbreak, rightbreak,node))
            tree_filehandle.write(",".join([str(cnode) for cnode, cspan in sorted(children.items()) if cspan[0]<rightbreak and cspan[1]>leftbreak]))
            tree_filehandle.write("\t{}\t{}\n".format(node, 0))
    
    for pos, node in sorted(mutations.items()):
        mutations_filehandle.write("{}\t{}\n".format(pos, node))

    if tree_filehandle is not sys.stdout:
        tree_filehandle.close()
    if mutations_filehandle is not sys.stdout:
        mutations_filehandle.close()                