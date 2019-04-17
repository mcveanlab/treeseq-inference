#!/usr/bin/env python3
"""
Various functions to convert sample data to argentum (tcPBWT) input format, and from
argentum output to newick."""
import sys
from math import ceil
import os.path
import logging

import numpy as np
    
def samples_to_argentum_in(sample_data, argentum_in_filehandle):
    """
    Takes a SampleData object, and outputs a file in .binary format, suitable for input 
    into argentum (see https://github.com/vlshchur/stable)
    
    This simply has one line per site with 0 (ancestral) and 1 (derived) states for each
    sample concatenated on a line. The first line is ignored, so we use it to store site
    positions, 
    """
    np.savetxt(
        argentum_in_filehandle, sample_data.sites_position[:], 
        newline=" ", header="argentum_input=sites:", delimiter="", comments="#")
    argentum_in_filehandle.write("\n")
    for id, genotype in sample_data.genotypes(): 
        np.savetxt(argentum_in_filehandle, genotype, fmt="%i", delimiter="", newline="")
        argentum_in_filehandle.write("\n")
    argentum_in_filehandle.flush()
    argentum_in_filehandle.seek(0)
        
def variant_positions_from_argentum_in(argentum_in_filename):
    with open(argentum_in_filename, "rt") as argentum_in_filehandle:
        # omit the first 
        assert argentum_in_filehandle.read(1) == "#"
        while argentum_in_filehandle.read(1) != ' ':
            pass #skip to first char after a space
        try:
            line = argentum_in_filehandle.readline().rstrip()
            positions = np.fromstring(line, sep=" ")
        except:
            logging.warning("Could not convert the first line to a set of floating point values:\n {}".format(line))
        return(positions)
 
def argentum_out_to_nexus(argentum_out, variant_positions, seq_length, outfilehandle):
    """
    argentum outputs a tree for every position, so there is a lot of duplication. 
    We can merge duplicate lines in this file as long as we keep track of which variants
    correspond to which tree.
    """
    with open(argentum_out, "rt") as argentum_out_fh:
        print("#NEXUS\nBEGIN TREES;", file = outfilehandle)
        buffered_tree = tree = ''
        site = 0
        for line in argentum_out_fh:
            # any line in an argentum output beginning with a brace is a tree
            if line.startswith("("):
                tree = line.rstrip()
                if tree != buffered_tree:
                    # argentum has many repeated tree lines. We only need to print out 1
                    if buffered_tree != '':
                        # Print the previous (buffered) tree with the new position 
                        # marking where we switch *off* this tree into the next
                        print("TREE", variant_positions[site], "=", buffered_tree, 
                            sep=" ",
                            end = "\n" if buffered_tree.endswith(';') else ";\n", 
                            file = outfilehandle)
                    buffered_tree = tree
                site += 1
        #print out the last tree
        if tree != '':
            print("TREE", str(seq_length), "=", tree, 
                sep=" ",
                end = "\n" if tree.endswith(';') else ";\n", 
                file = outfilehandle)
        print("END;", file = outfilehandle)
        outfilehandle.flush()