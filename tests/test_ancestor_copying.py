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
import string

import numpy as np
script_path = __file__ if "__file__" in locals() else "./dummy.py"
# use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime'))
# use the local copy of tsinfer in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','tsinfer')) 

import msprime
import tsinfer

rho = 2
ts = msprime.simulate(
    4, mutation_rate=1, recombination_rate=rho, random_seed=6)

#We want to construct a new ancestor for each non-sample node
# afterwards we could delete ancestors to create polytomies
#  (if there are no mutations on a branch)
# or create new ancestors 
#  (if there are multiple mutations or separated lengths of genome)

#make the array - should be a row for each node
column_mapping = np.zeros(ts.num_sites)
haplotype_matrix = -np.ones((ts.num_nodes, ts.num_sites), dtype=np.int)
copying_matrix   = -np.ones((ts.num_nodes, ts.num_sites), dtype=np.int)
root = np.zeros(ts.num_sites, dtype=np.int)
mutations = {k:[] for k in range(ts.num_nodes)}
for v in ts.variants(as_bytes=False):
    column_mapping[v.index] = v.position
    haplotype_matrix[0:ts.sample_size,v.site.index] = v.genotypes
    root[v.site.index] = int(v.site.ancestral_state)
    mutation_state = v.site.ancestral_state
    for m in v.site.mutations:
        # (site, prev state, new state)
        mutations[m.node].append((m.site, mutation_state,m.derived_state))
        mutation_state = m.derived_state


for es in ts.edgesets():
    # these are in youngest -> oldest order, so by iterating in the order given,
    # we should always be able to fill out parents from their children
    mask = (es.left <= column_mapping) & (column_mapping < es.right)
    previous_child = -1
    for child in es.children:
        #fill in the copying matrix
        copying_matrix[child,mask] = es.parent
        #create row for the haplotype matrix
        genotype = -np.ones(ts.num_sites, dtype=np.int)
        genotype[mask] = haplotype_matrix[child,mask]
        # add mutations. These are in oldest -> youngest order, so as we are working from
        # child->parent, we need to reverse
        for site, prev_state, new_state in reversed(mutations[child]):
            if mask[site]: #only set mutations on sites contained in this edgeset
                assert genotype[site] == int(new_state)
                genotype[site] = int(prev_state)
        
        if previous_child != -1:
            assert np.array_equal(genotype[mask], haplotype_matrix[es.parent,mask]), \
                "children {} and {} disagree".format(previous_child, child)
        
        haplotype_matrix[es.parent,mask] = genotype[mask]
        previous_child = child


#Check that the ancestral state matches the oldest variant for each column
for c in range(ts.num_sites):
    col = haplotype_matrix[:,c]
    filled = col[col != -1]
    assert filled[-1] == root[c], \
        "the ancestral state does not match the oldest variant for column {}".format(c)

#identify singletons, so we can exclude singleton sites (columns)
is_singleton = np.array([
    np.count_nonzero(v.genotypes != int(v.site.ancestral_state)) <= 1 \
        for v in ts.variants(as_bytes=False) \
    ], dtype=np.bool)

#identify nodes that are not used (or only used for singletons)
#these are cases where the node is used in the Treeseq, but since the
#genome span does not have any variable sites, it is not filled in the H matrix
is_unused_node = np.all(haplotype_matrix[:,~is_singleton] == -1, 1)

H = haplotype_matrix[:,~is_singleton][~is_unused_node,:]

#copying matrix is more tricky, as we need to renumber contents if we remove rows
bad_parent_id = -2 #used to relabel refs to rows we have removed (not -1 which is N/A)
#make the index one longer than the rows to allow mapping -1 -> -1
node_relabelling = np.full(copying_matrix.shape[0]+1, bad_parent_id, dtype=np.int)
node_relabelling[-1] = -1
node_relabelling[~is_unused_node]=np.arange(np.count_nonzero(~is_unused_node))
P = node_relabelling[copying_matrix[:,~is_singleton][~is_unused_node,:]]
assert not np.any(P == bad_parent_id)
assert not np.any(P >= P.shape[0])

print("True haplotype matrix (singletons & blank rows removed)")
for row in range(H.shape[0]):
    if row == ts.sample_size:
        print()
    print("{:>3}  ".format(row) + \
        " ".join(["{:>2}".format(x if x!=-1 else '*') for x in H[row,:]]))

print("True copying matrix "
    "(no singletons, blank rows removed and nodes renumbered accordingly)")
for row in range(P.shape[0]):
    if row == ts.sample_size:
        print()
    print("{:>3}  ".format(row) + \
        " ".join(["{:>2}".format(x if x!=-1 else '*') for x in P[row,:]]))


'''

S = np.zeros((ts.sample_size, ts.num_sites), dtype="u1")
for variant in ts.variants():
    S[:, variant.index] = variant.genotypes

sites = [mut.position for mut in ts.mutations()]
panel = tsinfer.ReferencePanel(S, sites, ts.sequence_length, rho=rho)



P, mutations = panel.infer_paths(num_workers=1)
'''
