#!/usr/bin/env python3
"""Various functions to convert tree sequences to RentPlus input format, and from RentPlus output to tree seqs."""
import sys
from math import ceil
import os.path
import logging

import numpy as np
    
def samples_to_RentPlus_in(sample_data, RentPlus_filehandle, infinite_sites=True):
    """
    Takes a SampleData object, and outputs a file in .dat format, suitable for input 
    into RentPlus (see https://github.com/SajadMirzaei/RentPlus)
    
    We could either use infinite sites, with a floating point number between 0 and 1, or
    integer sites. If integer sites, note that RentPlus does not claim to deal well with
    recurrent mutations.
    
    if infinite_sites==False, then use a basic discretising function, which simply rounds 
    upwards to the nearest int, ANDing the results if 2 or more variants end up at the same 
    integer position.
    
    if infinite_sites==True, then we must convert all the positions to lie between 0 and 1
    
    Note that with integer positions rentplus uses position coordinates (1,N) - i.e. [0,N). 
    That compares to tree sequences which use (0..N-1) - i.e. (0,N].
    
    Compared e.g. with fastARG, the matrix is transposed, so that rows are positions,
    this is a pain when trying to merge adjacent genotypes
    """
    position = sample_data.sites_position[:] #decompress all in one go to avoid sequential unpacking
    if infinite_sites:
        unique_positions = np.unique(position)
        assert unique_positions.shape[0] == sample_data.num_sites, \
            'If infinite sites is assumed, positions need to be unique, but there are dups'
        # normalize to between 0 and 1, as required by RentPlus for non-integer positions
        print(" ".join([str(float(p)/sample_data.sequence_length) for p in unique_positions]),
            file=RentPlus_filehandle)
        for id, haplotype in sample_data.haplotypes():
            print((haplotype + ord('0')).tobytes().decode(),file=RentPlus_filehandle)
    else:
        #compress runlengths of integers
        unique_positions = np.unique(np.ceil(position).astype(np.uint64))
        output = np.zeros((sample_data.num_samples, unique_positions.shape[0]), dtype=np.uint8)
        prev_position = 0
        index=0
        ANDed_genotype = None
        for id, genotype in sample_data.genotypes():
            if int(ceil(position[id])) != prev_position:
                #this is a new position. Save the genotype at the old position, and then reset everything
                if ANDed_genotype is not None:
                    output[:,index]=ANDed_genotype
                    index += 1
                ANDed_genotype = genotype
            else:
                ANDed_genotype =  np.logical_and(ANDed_genotype, genotype)
            prev_position = int(ceil(position[id]))
        if ANDed_genotype is not None: # print out the last site
            output[:,index]=ANDed_genotype
        
        assert index+1 == unique_positions.shape[0], \
            "Saved {} columns, should be {}".format(index+1, unique_positions.shape[0])
        np.savetxt(RentPlus_filehandle, output, "%u", delimiter="", comments='',
            header=" ".join([str(p) for p in unique_positions]))

    RentPlus_filehandle.flush()
    RentPlus_filehandle.seek(0)
        
    
def RentPlus_trees_to_nexus(trees_filename, outfilehandle, seq_length, num_tips):
    """
    RentPlus outputs a tree for every position, so there is a lot of duplication. 
    We can merge duplicate lines in this file as long as we keep track of which 
    variants correspond to which tree, and take the highest position in the file
    """
    with open(trees_filename, 'rt+') as RentPlusTrees:
        print("#NEXUS\nBEGIN TREES;", file = outfilehandle)
        # RentPlus creates 1-based tip numbers from a set of sequences, so we convert
        # back to 0-based by using the Nexus TRANSLATE functionality
        print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i+1,i) for i in range(num_tips)])), 
            file = outfilehandle)
        buffered_tree = tree = ''
        for line in RentPlusTrees:
            pos, tree = line.rstrip().split(None, 1)
            if tree:
                # RentPlus has many repeated tree lines. We only need to print out 1
                if tree != buffered_tree:
                    if buffered_tree:
                        # Print the previous (buffered) tree with the new position 
                        # marking where we switch *off* this tree into the next
                        print("TREE", pos, "=", buffered_tree, 
                            sep=" ",
                            end = "\n" if buffered_tree.endswith(';') else ";\n", 
                            file = outfilehandle)
                    buffered_tree = tree
        #print out the last tree
        if tree:
            print("TREE", str(seq_length), "=", tree, 
                sep=" ",
                end = "\n" if tree.endswith(';') else ";\n", 
                file = outfilehandle)
        print("END;", file = outfilehandle)
    outfilehandle.flush()