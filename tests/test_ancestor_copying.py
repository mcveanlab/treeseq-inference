#!/usr/bin/env python3
"""
Use msprime to create a simulation, and reconstruct the ancestral haplotype matrix, 
H from the simulation
Also keep track of which parts are from which ancestors (i.e. the true copying matrix P),
so that we can compare how we do in *infering* a copying matrix from the ancestral 
haplotype matrix
"""
import sys
import os
import operator

import numpy as np
script_path = __file__ if "__file__" in locals() else "./dummy.py"
# use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime'))
# use the local copy of tsinfer in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','tsinfer')) 

import msprime
import tsinfer
import msprime_to_inference_matrices

def print_matrix(M, break_at):
    for row in range(M.shape[0]):
        if row == break_at:
            print()
        print("{:>3}  ".format(row) + \
            " ".join(["{:>2}".format(x if x!=-1 else '*') for x in M[row,:]]))


rho = 4
ts = msprime.simulate(
    6, mutation_rate=2, recombination_rate=rho, random_seed=6)


#make haplotype (h) & sequence (p) matrices
h, p = msprime_to_inference_matrices.make_ancestral_matrices(ts)

#simplify down
is_singleton = msprime_to_inference_matrices.singleton_sites(ts)
is_unused_row = msprime_to_inference_matrices.unused_ancestors(h[:,~is_singleton])
H = h[:,~is_singleton][~is_unused_row,:]
P, row_map = msprime_to_inference_matrices.relabel_copy_matrix(p[:,~is_singleton], ~is_unused_row)
    
#check if any mutations are no longer covered
mutated_nodes = [m.node for m in ts.mutations()]
assert not np.any(row_map[mutated_nodes] < 0)

print("True haplotype matrix (singletons & blank rows removed)")
print_matrix(H, ts.sample_size)
print("True copying matrix "
    "(no singletons, blank rows removed and nodes renumbered accordingly)")
print_matrix(P, ts.sample_size)


#an alternative simplification
inf_order = [(m.node, np.sum(v.genotypes), v.site.index) for v in ts.variants() for m in v.site.mutations]
inf_order.sort(key=operator.itemgetter(1), reverse=True)
#exclude ancestors of singletons
freq_ordered_mutation_nodes = np.array([i[0] for i in inf_order if i[1] > 1], dtype=np.int)
#add the samples
keep_nodes = np.append(freq_ordered_mutation_nodes, range(ts.sample_size))
H = h[keep_nodes,:]
P, row_map = msprime_to_inference_matrices.relabel_copy_matrix(p, keep_nodes)

print("True haplotype matrix (nodes without mutations removed)")
print_matrix(H, H.shape[0]-ts.sample_size)
print("True copying matrix "
    "(nodes without mutations removed and nodes renumbered accordingly)")
print_matrix(P, P.shape[0]-ts.sample_size)

#We want to construct a new ancestor for each non-sample node
# afterwards we could delete ancestors to create polytomies
#  (if there are no mutations on a branch)
# or create new ancestors 
#  (if there are multiple mutations or separated lengths of genome)
