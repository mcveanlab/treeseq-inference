#!/usr/bin/env python3
"""
Take an msprime simulation and output the number of nodes with no mutations above them
"""
import sys
import os
import argparse
from collections import defaultdict

import numpy as np

script_path = __file__ if "__file__" in locals() else "./dummy.py"
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','tsinfer')) # use the local copy of msprime in preference to the global one
import msprime
import tsinfer

def has_mutation(mutations, node, start_end):
    for mutation_location in mutations.get(node, []):
        if mutation_location >= start_end[0] and mutation_location < start_end[1]:
            return True
    return False
    
def update_edges(edges, out_edges, in_edges):
    #update edges
    for node, edge in out_edges.items():
        del edges[node]
    for node, edge in in_edges.items():
        edges[node] = edge

def generate_samples(ts):
    S = np.zeros((ts.sample_size, ts.num_mutations), dtype="i1")
    for variant in ts.variants():
        S[:,variant.index] = variant.genotypes
    pos = np.array([v.position for v in ts.variants()])
    return S, pos

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Output information on ARG nodes.", 
    )
    parser.add_argument('-v', '--verbosity', action="count")
    args = parser.parse_args()


    for l in range(50000, 1500000, 50000):
        recombination_rate=1e-8
        mutation_rate=1e-8
        ts = msprime.simulate(20, 
            Ne=50000, length=l, recombination_rate=recombination_rate, mutation_rate=mutation_rate, 
            random_seed=6)
        nodes_with_mutations = {m.node:True for m in ts.mutations()}
        nodes_with_recombinations = {}
        nodes_with_informative_recombinations = {}
        #save mutations - note that this doesn't account for back-mutations
        mutations = defaultdict(list)
        for m in ts.mutations():
            mutations[m.node].append(m.position)
        diff_iterator = ts.edge_diffs()
        tree_iterator = ts.trees()
        start_end, o_rec, i_rec = next(diff_iterator) #first diff 
        #keep our own record of the live edgesets, which will allow us to reconstruct the tree
        edges = {i.child:i for i in i_rec} #save edges indexed by node number (should be unique)
        for ((location, end_location), o_rec, i_rec), tree in zip(diff_iterator, tree_iterator):
            updated = False
            out_edges = {o.child:o for o in o_rec}
            in_edges = {i.child:i for i in i_rec}
            assert len(o_rec) == len(i_rec), (len(o_rec), len(i_rec))
            #should either have 2, 3, or 4 changed edges
            assert len(o_rec) <= 4 and len(o_rec) > 1, len(o_rec)
            if len(o_rec) >= 3:
                #find the deleted parent node
                out_parent = set([r.parent for r in o_rec]) - set([r.parent for r in i_rec])
                #find the created parent node
                in_parent = set([r.parent for r in i_rec]) - set([r.parent for r in o_rec])
                #symmetric difference - these should be the deleted and newly introduced parents
                changed_parents = out_parent | in_parent
                #the moved node (above which the recombination has occurred) is the one that shares these parents
                moved_node = set([r.child for r in o_rec if r.parent in changed_parents]) & set([r.child for r in i_rec if r.parent in changed_parents])
                #
                assert len(moved_node) <= 2
                if len(moved_node) == 1:
                    try:
                        if args.verbosity:
                            print("moved {}, out_parent {}, in parent {}".format(moved_node, out_parent, in_parent))
                            print(edges.keys())
                            print(out_edges)
                            print(in_edges)
                        moved_node = moved_node.pop()
                        assert(len(out_parent)==1)
                        assert(len(in_parent)==1)
                            
                        nodes_with_recombinations[moved_node] = True
                        #is this node informative (requires the two branches on each side to acquire different mutation patterns)
                        #could we possibly check 2 different ways, i.e. by looking at the halotypes themselves (finding unshared (unnested) mutations)
                        #or by going up the ancestry, gradually reducing in size until we hit the same ancestor. The ancestry method is coded up here
                        #but it is an underestimate, as we lose some information in cases where msprime throws away parts of the ancestry (where the
                        #trees change rooting).
                        out_edge = out_edges[moved_node]
                        in_edge = in_edges[moved_node]
                        
                        #find MRCA
                        ancestors=set()
                        focal_node = out_edge.parent
                        while focal_node is not None:
                            ancestors.add(focal_node)
                            focal_node = edges[focal_node].parent if focal_node in edges else None
                        if args.verbosity:                        
                            print("ancestors: {}".format(ancestors))
                        focal_node = in_edge.parent
                        while focal_node not in ancestors:
                            if focal_node in in_edges:
                                #check within new edges first, in case this is a replacement
                                focal_node = in_edges[focal_node].parent
                            else:
                                if focal_node in edges:
                                    #NB - this could include edges which will not be present after node mutation
                                    focal_node = edges[focal_node].parent
                                else:
                                    focal_node = None
                                    break;
                        MRCA = focal_node
                        if args.verbosity:
                            print("MRCA = {}".format(MRCA))
                        #(if MRCA is none, we just go as far up the tree as possible)
                        
                        #now look for mutations on LHS
                        unrecombined_LH_boundary = out_edge.left
                        LHS_mutation = has_mutation(mutations, out_edge.child, (unrecombined_LH_boundary, location))
                        if not LHS_mutation:
                            #look in parents of this node, all the way to the MRCA, but only if unrecombined
                            focal_node = out_edge.parent
                            while focal_node != MRCA:
                                if focal_node not in edges:
                                    break;
                                    #we dan't know anything more about this route, so have to bail: caution, this makes the informative stat an underestimate
                                unrecombined_LH_boundary = max(unrecombined_LH_boundary, edges[focal_node].left)
                                #loop over mutations for this node
                                if has_mutation(mutations, focal_node, (unrecombined_LH_boundary, location)):
                                    LHS_mutation = True
                                    break;
                                focal_node = edges[focal_node].parent
                                
                        update_edges(edges, out_edges, in_edges)
                        updated = True
                        
                        unrecombined_RH_boundary = in_edge.right
                        RHS_mutation = has_mutation(mutations, in_edge.child, (location, unrecombined_RH_boundary))
                        if not RHS_mutation:
                            #look in parents of this node, all the way to the MRCA, but only if unrecombined
                            focal_node = in_edge.parent
                            while focal_node != MRCA:
                                if focal_node not in edges:
                                    break;
                                    #we dan't know anything more about this route, so have to bail: caution, this makes the informative stat an underestimate
                                unrecombined_RH_boundary = min(unrecombined_RH_boundary, edges[focal_node].right)
                                #loop over mutations for this node
                                if has_mutation(mutations, focal_node, (location, unrecombined_RH_boundary)):
                                    RHS_mutation = True
                                    break;
                                focal_node = edges[focal_node].parent
        
                        
                        if LHS_mutation and RHS_mutation:
                            nodes_with_informative_recombinations[moved_node] = True
                            if args.verbosity:
                                print("informative SRB on node {}".format(moved_node))
                        else:
                            if args.verbosity:
                                print("uninformative SRB on node {}".format(moved_node))
                    except:
                        print(tree.draw(format="unicode"))
                        print(next(tree_iterator).draw(format="unicode"))
                        raise
    
                elif len(moved_node) == 2:
                    #two children have simultaneously swapped from one specific parent to another specific parent
                    #we can't tell which node the recombination occurred above, but it may not matter, as the
                    #information on the split is the same
                    assert set([r.child for r in o_rec if r.parent in changed_parents]) == set([r.child for r in i_rec if r.parent in changed_parents])
            elif len(o_rec) == 2:
                #two children have simultaneously swapped from one specific parent to another specific parent
                #we can't tell which node the recombination occurred above, but it may not matter, as the
                #information on the split is the same            
                assert set([r.child for r in o_rec]) == set([r.child for r in i_rec])
                assert o_rec[0].parent == o_rec[0].parent
                assert i_rec[0].parent == i_rec[0].parent
    
            if not updated:
                update_edges(edges, out_edges, in_edges)
                updated = True

        #now compare to inference
        S, pos = generate_samples(ts)
        inferred_ts1 = tsinfer.infer(S, pos, l, recombination_rate, resolve_shared_recombinations=False).simplify()
        inferred_ts2 = tsinfer.infer(S, pos, l, recombination_rate, resolve_shared_recombinations=True).simplify()
        
        resolvable_nodes = len({**nodes_with_mutations, **nodes_with_informative_recombinations})
                
        print("For {:.3f}Mb {} / {} resolvable nodes, {} ({:.2f}%) via mutations (always resolvable in msprime), {} ({:.2f}%) via recombinations (of {} or {:.2f}% nodes with breakpoints above). There are {} ({:.2f}%) unresolvable nodes. Recombinations have added {} ({}%)".format(
            l/1e6, resolvable_nodes, ts.get_num_nodes(), 
            len(nodes_with_mutations), len(nodes_with_mutations)/ts.get_num_nodes()*100,
            len(nodes_with_informative_recombinations), len(nodes_with_informative_recombinations)/ts.get_num_nodes()*100,
            len(nodes_with_recombinations), len(nodes_with_recombinations)/ts.get_num_nodes()*100,
            ts.get_num_nodes() - resolvable_nodes, 100-(resolvable_nodes/ts.get_num_nodes()*100),
            resolvable_nodes - len(nodes_with_mutations), (resolvable_nodes - len(nodes_with_mutations))/ts.get_num_nodes()*100))
        print("From tsinfer: {} nodes inferred without using shared breakpoints, {} nodes using SRBs (an increase of {} nodes)".format(
            inferred_ts1.get_num_nodes(), inferred_ts2.get_num_nodes(), inferred_ts2.get_num_nodes()-inferred_ts1.get_num_nodes()))
        
        
