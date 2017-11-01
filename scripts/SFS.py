import msprime
import collections
import math

check_method = True

recombination_rate=1e-8
mutation_rate=1e-8
ts = msprime.simulate(10000, 
    Ne=10000, length=1000000, recombination_rate=recombination_rate, mutation_rate=mutation_rate, 
    random_seed=6)
trees = ts.trees(leaf_counts=True)
diffs = ts.edge_diffs()
SFS = collections.Counter()

#first iteration introduces all the edges
diff = next(diffs)
tree = next(trees)
for node in tree.nodes():
    try:
        SFS[tree.get_num_leaves(node)] += tree.get_branch_length(node)
    except ValueError:
        #can ignore root node (no branch length, mutations above it are fixed anyway) but not other nodes
        if tree.get_parent(node) != msprime.NULL_NODE:
            raise
            
for i, ((start, end), out_edges, in_edges) in enumerate(diffs):
    try:
        #REMOVE OUT EDGES: must check out nodes + ancestors *and* ancestors of in_nodes
        #as an in node parent with >2 children will not appear in out_nodes 
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
                #can ignore root node (no branch length, mutations above it are fixed anyway) but not other nodes
                if tree.get_parent(node) != msprime.NULL_NODE:
                    raise
        
        tree = next(trees)
        #ADD IN EDGES: must check out nodes + ancestors *and* ancestors of out_nodes
        #as an out node parent with >2 children will not appear in in_nodes 
        check_nodes = set([edge.child for edge in in_edges] + [edge.parent for edge in out_edges])
        added_nodes = set()
        for node in check_nodes:
            #need to go up from this node and modify
            while node not in added_nodes and node != msprime.NULL_NODE:
                added_nodes.add(node)
                node = tree.get_parent(node)
        for node in added_nodes:
            try:
                SFS[tree.get_num_leaves(node)] += tree.get_branch_length(node)
            except ValueError:
                #can ignore root node (no branch length, mutations above it are fixed anyway) but not other nodes
                if tree.get_parent(node) != msprime.NULL_NODE:
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