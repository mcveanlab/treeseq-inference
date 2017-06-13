#!/usr/bin/env python3
"""
Use msprime to create a simulation, and reconstruct the ancestral haplotype matrix, H from the simulation
Also keep track of which parts are from which ancestors (i.e. the true copying matrix P), so that we can compare
how we do in *infering* a copying matrix from the ancestral haplotype matrix
"""
import sys
import os
import string

import numpy as np
script_path = __file__ if "__file__" in locals() else "./dummy.py"
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','tsinfer')) # use the local copy of tsinfer in preference to the global one
    
import msprime
import tsinfer

rho = 2
ts = msprime.simulate(
    5, mutation_rate=5, recombination_rate=rho, random_seed=6)

#We want to construct a new ancestor for each non-sample node
# afterwards we could delete ancestors to create polytomies (if there are no mutations on a branch)
# or create new ancestors (if there are multiple mutations or separated lengths of genome)


#make the array - should be a row for each node
column_mapping = np.zeros(ts.num_sites)
H = -np.ones((ts.num_nodes, ts.num_sites), dtype=np.int)
root = np.zeros(ts.num_sites, dtype=np.int)
mutations = {k:[] for k in range(ts.num_nodes)}
for v in ts.variants(as_bytes=False):
    column_mapping[v.index] = v.position
    H[0:ts.sample_size,v.site.index] = v.genotypes
    root[v.site.index] = int(v.site.ancestral_state)
    mutation_state = v.site.ancestral_state
    for m in v.site.mutations:
        mutations[m.node].append((m.site, mutation_state,m.derived_state)) # (site, prev state, new state)
        mutation_state = m.derived_state


for es in ts.edgesets():
    #these are in youngest -> oldest order, so by iterating in the order given,
    # we should always be able to fill out parents from their children
    mask = (es.left <= column_mapping) & (column_mapping < es.right)
    previous_child = -1
    for child in es.children:
        genotype = -np.ones(ts.num_sites, dtype=np.int)
        genotype[mask] = H[child,mask]
        #add mutations. These are in oldest -> youngest order, so as we are working from child->parent, we need to reverse
        for (site, prev_state, new_state) in reversed(mutations[child]):
            if mask[site]: #only set mutations on sites contained in this edgeset
                assert genotype[site] == int(new_state), "state of child {} for pos {} is not {}".format(child, site, new_state)
                genotype[site] = int(prev_state)
        
        if previous_child != -1:
            assert np.array_equal(genotype[mask],H[es.parent,mask]), "children {} and {} disagree".format(previous_child, child)
        
        H[es.parent,mask] = genotype[mask]
        previous_child = child

with open("ancestral_haplotypes.txt", "w") as text_file:
    for row in range(H.shape[0]):
        print("{:>3}  ".format(row)+ " ".join([(str(x) if x!=-1 else '*') for x in H[row,:]]), file=text_file)

trees = ts.newick_trees()
next(trees)

'''

S = np.zeros((ts.sample_size, ts.num_sites), dtype="u1")
for variant in ts.variants():
    S[:, variant.index] = variant.genotypes

sites = [mut.position for mut in ts.mutations()]
panel = tsinfer.ReferencePanel(S, sites, ts.sequence_length, rho=rho)



P, mutations = panel.infer_paths(num_workers=1)
'''