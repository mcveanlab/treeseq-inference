#!/usr/bin/env python3
"""Various functions to convert tree sequences to RentPlus input format, and from RentPlus output to tree seqs."""
import sys
from math import ceil
import numpy as np
import os.path
import logging
    
def samples_to_RentPlus_in(sample_data, RentPlus_filehandle, continuous_positions=True):
    """
    Takes a SampleData object, and outputs a file in .dat format, suitable for input 
    into RentPlus (see https://github.com/SajadMirzaei/RentPlus)
    
    We could either use infinite sites, with a floating point number between 0 and 1, or
    integer sites. If integer sites, note that RentPlus does not claim to deal well with
    recurrent mutations.
        
    Note that with integer positions rentplus uses position coordinates (1,N) - i.e. [0,N). 
    That compares to tree sequences which use (0..N-1) - i.e. (0,N].
    
    Compared e.g. with fastARG, the matrix is transposed, so that rows are positions,
    this is a pain when trying to merge adjacent genotypes
    """
    position = sample_data.sites_position[:] #decompress all in one go to avoid sequential unpacking
    assert sample_data.num_sites == np.unique(position).shape[0], 'Site positions need to be unique, but there are dups'
    if continuous_positions:
        #normalize to between 0 and 1
        print(" ".join([str(float(p)/sample_data.sequence_length) for p in position]), 
            file=RentPlus_filehandle)
        for id, genotypes in sample_data.genotypes():
            print((genotypes.astype("i1") + ord('0')).tobytes().decode(),file=RentPlus_filehandle)
    else:
        output = np.zeros((sample_data.num_samples, position.shape[0]), dtype=np.uint8)
        for i, genotypes in sample_data.genotypes():
            output[:,i]=genotypes
        np.savetxt(RentPlus_filehandle, output, "%u", delimiter="", comments='',
            header=" ".join([str(int(p+1)) for p in position]))

    RentPlus_filehandle.flush()
    RentPlus_filehandle.seek(0)
        
    
def RentPlus_trees_to_nexus(trees_filename, outfilehandle, seq_length, num_tips):
    """
    RentPlus outputs a tree for every position, so there is a lot of duplication. We can merge duplicate lines
    in this file as long as we keep track of which variants correspond to which tree, and take the average 
    of the variant positions in the file
    """
    with open(trees_filename, 'rt+') as RentPlusTrees:
        print("#NEXUS\nBEGIN TREES;", file = outfilehandle)
        #RentPlus creates 1-based tip numbers from a set of sequences, so we convert back to 0-based
        #by using the Nexus TRANSLATE functionality
        print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i+1,i) for i in range(num_tips)])), 
            file = outfilehandle)
        oldline = 0,''
        for line in RentPlusTrees:
            pos, tree = line.rstrip().split(None, 1)
            if tree:
                #RentPlus has many repeated tree lines. We only need to print out 1
                if tree != oldline[1]:
                    if oldline[1]:
                        print("TREE", str((float(oldline[0])+float(pos))/2), "=", tree,
                            sep=" ",
                            end = "\n" if tree.endswith(';') else ";\n", 
                            file = outfilehandle)
                    oldline = pos, tree
        if oldline[1]:
            print("TREE", str(seq_length), "=", tree, end = "\n" if tree.endswith(';') else ";\n", 
                file = outfilehandle)
        print("END;", file = outfilehandle)
    outfilehandle.flush()