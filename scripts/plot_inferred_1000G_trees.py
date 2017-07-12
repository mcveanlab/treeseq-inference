#!/usr/bin/env python3
"""
Create an svg file for each tree in a 1000G inferred tree sequence. Convert to a movie using something like

> convert 'tree*.svg' tree.mpg

or use ffmpeg if `convert` is too slow/memory intensive

> mogrify -format png 'tree*.svg'
> ffmpeg -r 25 -i tree%05d.png -pix_fmt yuv420p tree.mp4


"""
import sys
import os

import csv
import h5py
script_path = __file__ if "__file__" in locals() else "./dummy.py"
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime')) # use the local copy of msprime in preference to the global one
import msprime


def percolate_unambiguous_colours(tree, colours, node):
    """
    Percolate colours up the tree where all children have the same colour    
    """
    if tree.is_internal(node):
        children = tree.get_children(node)
        for c in children:
            percolate_unambiguous_colours(tree, colours, c)
        if node not in colours:
            unique_colours = set([colours.get(c) for c in children]) #could be None
            if len(unique_colours) == 1:
                single_colour = unique_colours.pop()
                if single_colour is not None:
                    colours[node] = single_colour

def colour_children(tree, colours):
    ret_colours = {}
    for node in tree.nodes():
        if node in colours:
            for c in tree.get_children(node):
                ret_colours[c]=colours[node]
    return ret_colours

#from http://www.internationalgenome.org
pop_cols = {"African":"yellow", "American": "red", "East Asian":"green","European":"blue", "South Asian":"purple"}
sample_cols = {}
with open("/Volumes/SDdisk/1000G/igsr_samples.tsv") as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter='\t')
    for row in reader:
        if 'phase 1' in row['Data collections']:
            sample_cols[row['Sample name']]= pop_cols[row['Superpopulation name']]

ts = msprime.load("../tmp/1000G22_infer2000.hdf5").simplify()
#sample names are not stored in msprime files, so we need to get them from the base file
#sample names in the base file have an extra character ('a' or 'b') appended
with h5py.File("/Volumes/SDdisk/1000G/1000G22.hdf5", "r") as f:
    data = f['data']
    samples = data['samples']
    sample_colours = {i:sample_cols[str(d, 'utf8')[:-1]] for i,d in enumerate(data['samples']) if str(d, 'utf8')[-1] in 'ab'}

for i, t in enumerate(ts.trees()):
    if i==0:
        start = t.get_interval()[1]
    branch_colours = {},{}
    node_colours=sample_colours.copy()
    percolate_unambiguous_colours(t, node_colours, t.get_root(), )
    branch_colours = colour_children(t, node_colours)
    file = "/Volumes/SDdisk/1000G/tree{:05d}.svg".format(i)
    svg=t.draw(file,
        1400,850,
        branch_colours=branch_colours, node_colours=node_colours)
    print("written {} spanning genome interval {} ({} base pairs). Total={:.3f}Mb".format(
        file, t.get_interval(), t.get_length(),(t.get_interval()[1]-start)/1e6))    
