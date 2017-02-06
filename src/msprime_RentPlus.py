#!/usr/bin/env python3
"""Various functions to convert msprime simulation files to RentPlus input format, and from RentPlus output to msprime."""
import sys
from math import ceil
import numpy as np
import os.path
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','msprime')) # use the local copy of msprime in preference to the global one
import msprime
import logging
    
def variant_matrix_to_RentPlus_in(var_matrix, var_positions, seq_length, RentPlus_filehandle, infinite_sites=True):
    """
    Takes an variant matrix, and outputs a file in .dat format, suitable for input 
    into RentPlus (see https://github.com/SajadMirzaei/RentPlus)
    
    We could either use infinite sites, with a floating point number between 0 and 1, or
    integer sites. If integer sites, note that RentPlus does not claim to deal well with
    recurrent mutations.
    
    if infinite_sites==False, then use a basic discretising function, which simply rounds 
    upwards to the nearest int, ANDing the results if 2 or more variants end up at the same 
    integer position.
    
    if infinite_sites==True, then we must convert all the positions to lie between 0 and 1
    
    Note that with integer positions rentplus uses position coordinates (1,N) - i.e. [0,N). 
    That compares to msprime which uses (0..N-1) - i.e. (0,N].
    
    Compared e.g. with fastARG, the matrix is transposed, so that rows are positions,
    this is a pain when trying to merge adjacent genotypes
    """
    n_variants, n_samples = var_matrix.shape
    assert len(var_matrix)==n_variants
    if infinite_sites:
        #normalize to between 0 and 1
        print(" ".join([float(p)/seq_length for p in var_positions]), file=RentPlus_filehandle)
        for v in var_matrix:
            print("".join(v), file=RentPlus_filehandle)
    else:
        #compress runlengths of integers - we need to do this on the first line and also 
        prev_position = 0
        delim=""
        for pos in var_positions:
            if int(ceil(pos)) != prev_position:
                RentPlus_filehandle.write(delim + str(int(ceil(pos))))
                delim=" "
            prev_position = int(ceil(pos))

        prev_position = 0
        ANDed_genotype = None
        for pos, genotype in zip(var_positions, var_matrix):
            if int(ceil(pos)) != prev_position:
                #this is a new position. Print the genotype at the old position, and then reset everything
                if prev_position:
                    print("".join(ANDed_genotype), file=RentPlus_filehandle)
                ANDed_genotype = genotype
            else:
                ANDed_genotype =  np.logical_and(ANDed_genotype, genotype)
            prev_position = int(ceil(pos))
        if ANDed_genotype is not None: # print out the last site
            print("".join(ANDed_genotype), file=RentPlus_filehandle)

    RentPlus_filehandle.flush()
    RentPlus_filehandle.seek(0)
        
    
def RentPlus_trees_to_nexus(trees_filename, outfilehandle, seq_length, num_tips, zero_based_tip_numbers=True):
    """
    RentPlus outputs a tree for every position, so there is a lot of duplication. We can merge duplicate lines
    in this file as long as we keep track of which variants correspond to which tree, and take the average 
    of the variant positions in the file
    
    We generally treat the tree numbers as 
    setting zero_based_tip_numbers=True is not supported
    """
    #by default, RentPlus produces 1-based tip numbers - if we want 0-base, we need to convert them
    def tip_subtract_1(tree):
        import re
        return re.sub(r'(?<=[(,)])(\d+)',lambda m: str(int(m.group(1))-1), tree)
    increment = 0 if zero_based_tip_numbers else 1
    with open(trees_filename, 'rt+') as RentPlusTrees:
        print("#NEXUS\nBEGIN TREES;", file = outfilehandle)
        print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i+increment,i+increment) for i in range(num_tips)])), 
            file = outfilehandle)
        oldline = 0,''
        for line in RentPlusTrees:
            pos, tree = line.rstrip().split(None, 1)
            if tree:
                #RentPlus has many repeated tree lines. We only need to print out 1
                if tree != oldline[1]:
                    if oldline[1]:
                        print("TREE " + str((float(oldline[0])+float(pos))/2) + " = " +
                            tip_subtract_1(tree) if zero_based_tip_numbers else tree, 
                            end = "\n" if tree.endswith(';') else ";\n")
                    oldline = pos, tree
        if oldline[1]:
            print("TREE " + str((float(oldline[0])+float(pos))/2) + " = " +
                tip_subtract_1(tree) if zero_based_tip_numbers else tree, 
                end = "\n" if tree.endswith(';') else ";\n")
        print("END;", file = outfilehandle)
    outfilehandle.flush()