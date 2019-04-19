#!/usr/bin/env python3
"""
Various functions to convert sample data to argentum (tcPBWT) input format, and from
argentum output to newick."""
import sys
from math import ceil
import os.path
import logging

import numpy as np
    
def samples_to_argentum_in(sample_data, argentum_in_filehandle, positions_filehandle):
    """
    Takes a SampleData object, and outputs a file in .binary format, suitable for input 
    into argentum (see https://github.com/vlshchur/argentum)
    
    This simply has one line per site with 0 (ancestral) and 1 (derived) states for each
    sample concatenated on a line. 
    """
    for id, genotype in sample_data.genotypes(): 
        np.savetxt(argentum_in_filehandle, genotype, fmt="%i", delimiter="", newline="")
        argentum_in_filehandle.write("\n")
    argentum_in_filehandle.flush()
    argentum_in_filehandle.seek(0)
    np.savetxt(
        positions_filehandle, sample_data.sites_position[:], newline=" ", delimiter="")
    positions_filehandle.flush()
    positions_filehandle.seek(0)
        
def variant_positions_from_fn(positions_filename):
    with open(positions_filename, "rt") as positions_fh:
        for line in positions_fh:
            try:
                positions = np.fromstring(line.rstrip(), sep=" ")
            except:
                logging.warning("Could not convert the first line to a set of floating point values:\n {}".format(line))
            return(positions)


def planar_order_to_newick(order_string, height_string, branch_lengths=True):
    orders = order_string.split(",")
    heights = np.fromstring(height_string, sep=",")
    while len(heights):
        # cluster tips together, starting at the smallest height
        target_height = np.min(heights)
        runs_of_min_height = (heights == target_height)
        # make sure all runs of ones are well-bounded
        bounded = np.hstack(([0], runs_of_min_height, [0]))
        # get 1 at run starts and -1 at run ends
        diffs = np.diff(bounded)
        removed = 0
        nwk_height = ":{}".format(target_height) if branch_lengths else ""
        for start, end in zip(np.where(diffs > 0)[0], np.where(diffs < 0)[0]):
            # merge adjacent tips of the same min_height together,
            # comma-separated, within a brace (i.e. newick format)
            start -= removed
            end -= removed
            orders[end] = "({})".format(
                               ",".join([o+nwk_height for o in orders[start:(end+1)]]))
            del orders[start:end]
            removed += end-start
        heights = heights[~runs_of_min_height]
    return(orders[0] + ";")

def argentum_out_to_nexus(argentum_out, variant_positions, seq_length, outfilehandle):
    """
    argentum outputs a pair or "planar order" lines for every position, so there is a lot
    of duplication. 
    We can merge duplicate lines in this file as long as we keep track of which variants
    correspond to which tree.
    """
    with open(argentum_out, "rt") as argentum_out_fh:
        print("#NEXUS\nBEGIN TREES;", file = outfilehandle)
        buffered_planar_order = ("", "")
        height = order = ''
        site = 0
        for order in argentum_out_fh:
            if not order[0].isdigit():
                continue # This is not a planar order line
            order = order.rstrip()
            height = next(argentum_out_fh).rstrip()
            if buffered_planar_order != (order, height):
                # argentum has many repeated tree lines. We only need to print out 1
                if buffered_planar_order[0] != '':
                    # Print the previous (buffered) tree with the new position 
                    # marking where we switch *off* this tree into the next
                    newick_tree = planar_order_to_newick(*buffered_planar_order)
                    print("TREE", variant_positions[site], "=", newick_tree, 
                        sep=" ",
                        end = "\n" if newick_tree.endswith(';') else ";\n", 
                        file = outfilehandle)
                buffered_planar_order = (order, height)
            site += 1
        #print out the last tree
        if site != len(variant_positions):
            raise ValueError("argentum bug hit: {} trees but {} sites"
                .format(site, len(variant_positions)))
        if height != '':
            newick_tree = planar_order_to_newick(*buffered_planar_order)
            print("TREE", str(seq_length), "=", newick_tree, 
                sep=" ",
                end = "\n" if newick_tree.endswith(';') else ";\n", 
                file = outfilehandle)
        print("END;", file = outfilehandle)
        outfilehandle.flush()