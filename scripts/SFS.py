import msprime
import collections
import math

recombination_rate=1e-8
mutation_rate=1e-8
ts = msprime.load("../ftprime/sweepfile0.8.hdf5")

check_method = True

trees = ts.trees(leaf_counts=True)
diffs = ts.edge_diffs()

SFS = collections.Counter()
for i, ((start, end), out_edges, in_edges) in enumerate(diffs):
    try:
        if i:
            check_nodes = set([edge.child for edge in out_edges] + [edge.parent for edge in in_edges])
            removed_nodes = set()
            for node in check_nodes:
                #should be no out_edges at the start
                #need to go up from this node and modify
                while node not in removed_nodes and node != msprime.NULL_NODE:
                    removed_nodes.add(node)
                    node = tree.get_parent(node)
            for node in removed_nodes:
                try:
                    SFS[tree.get_num_leaves(node)] -= tree.get_branch_length(node)
                except ValueError:
                    if tree.get_parent(node) == msprime.NULL_NODE:
                        #the root node will not have a branch length, and can be ignored anyway
                        pass
                    else:
                        raise
        tree = next(trees)
        check_nodes = set([edge.child for edge in in_edges] + [edge.parent for edge in out_edges])
        added_nodes = set()
        for node in check_nodes:
            #need to go up from this node and modify
            while node not in added_nodes and node != msprime.NULL_NODE:
                added_nodes.add(node)
                node = tree.get_parent(node)
        for n in added_nodes:
            try:
                SFS[tree.get_num_leaves(n)] += tree.get_branch_length(n)
            except ValueError:
                if tree.get_parent(n) == msprime.NULL_NODE:
                    #the root node will not have a branch length, and can be ignored anyway
                    pass
                else:
                    raise
        if check_method:
            SFS2 = collections.Counter()
            for node in tree.nodes():
                if tree.get_parent(node) != msprime.NULL_NODE:
                    SFS2[tree.get_num_leaves(node)] += tree.get_branch_length(node)
            for k in SFS.keys():
                assert math.isclose(SFS[k],SFS2[k], abs_tol=1e-9), (i, k, SFS[k],SFS2[k])
        print("Done tree {} from {} to {} ({} root{})".format(i, start, end, tree.num_roots, 's' if tree.num_roots > 1 else ''))
        
    except:
        print(tree.num_roots)
        print([n for n in tree.nodes() if tree.get_parent(n) == 12557])
        print(out_edges)
        print(in_edges)
        #print(tree.draw(format="unicode"))
        raise