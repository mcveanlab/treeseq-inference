def write_nexus_trees(ts, treefile, index_trees_by_variants=False, zero_based_tip_numbers=True):
    """
    if index_trees_by_variants == False (the default), then the names of the trees in the file correspond to the
    upper breakpoint position along the genome. If index_trees_by_variants == True then each tree name is instead
    the number of variants observed along the sequence so far (i.e. the uppermost variant index + 1). 
    If index_tree_by_variants is a list, it is taken to be a list of variant positions to use as positions at which
    trees are output (with adjacent identical trees output as a single tree and indexed by the uppermost variant position,
    in the same manner as for index_trees_by_variants = True). 
    
    Using the variant position to index trees allows trees to be compared for batches of variants rather than 
    via chromosome position, allowing fair comparison between inferred trees.
    """
    import sys
    import numpy as np
    
    print("#NEXUS\nBEGIN TREES;", file = treefile)
    increment = 0 if zero_based_tip_numbers else 1
    print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i+increment,i+increment) for i in range(ts.get_sample_size())])), file = treefile)
    variant_index = 0
    epsilon = 1e-8
    if hasattr(index_trees_by_variants, "__len__"):
        #index_trees_by_variants contains an array of variant positions, which can be used to index the trees
        assert min(index_trees_by_variants) >= 0, "The variant positions passed in must all be greater than or equal to 0"
        assert max(index_trees_by_variants) < ts.get_sequence_length(), "The variant positions passed in must all be less than {}".format(ts.get_sequence_length())
        positions = np.sort(np.array(index_trees_by_variants))
    else:
        positions = None
        
    for t, (_, newick) in zip(ts.trees(), ts.newick_trees()): #TO DO: should do ts.newick_trees(zero_based_tip_numbers)
        if index_trees_by_variants == False:
            #index by rightmost genome position
            print("TREE " + str(t.get_interval()[1]) + " = " + newick, file=treefile)
        else:
            #index by rightmost variant number
            if positions is None: #use the built-in variant positions
                n = t.get_num_mutations()
            else:
                l, r = t.get_interval()
                n = np.count_nonzero((positions >= l) & (positions < r))
            if n:
                #only print a tree if is is covered by at least 1 variant
                variant_index += n
                print("TREE " + str(variant_index) + " = " + newick, file=treefile)
    print("END;", file = treefile)
