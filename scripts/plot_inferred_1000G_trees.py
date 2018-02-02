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
import argparse

import csv
import h5py
script_path = __file__ if "__file__" in locals() else "./dummy.py"
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime')) # use the local copy of msprime in preference to the global one
import msprime

parser = argparse.ArgumentParser(description='Plot a large set of 1000 genomes trees from a treeseq.')
parser.add_argument('msprime_infile', default="1000G/1000G22_infer2000.hdf5",
                    help='an msprime hdf5 file inferred from 1000 genomes data (phase 1)')
parser.add_argument('base_file', default="1000G/1000G22.hdf5",
                    help='the equivalent hdf5 file from which the msprime file was inferred,' +
                        " only needed because we can't currently store sample names in the msprime file")
parser.add_argument('1000genomes_tsv', default="1000G/igsr_samples.tsv",
                    help='a path to the igsr_samples.tsv file')
parser.add_argument('--outdir', '-d', default="1000G",
                    help='the directory in which to save svg files')
args = parser.parse_args()

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

ts = msprime.load(args.msprime_infile).simplify()
#sample names are not stored in msprime files, so we need to get them from the base file
#sample names in the base file have an extra character ('a' or 'b') appended
with h5py.File(args.base_file, "r") as f:
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
    file = os.path.join(args.outdir,"tree{:05d}.svg".format(i))
    svg=t.draw(file,
        1400,850,
        show_mutation_labels=True, show_internal_node_labels=False,
        show_leaf_node_labels=False,
        branch_colours=branch_colours, node_colours=node_colours)
    print("written {} spanning genome interval {} ({} base pairs). Total={:.3f}Mb".format(
        file, t.get_interval(), t.get_length(),(t.get_interval()[1]-start)/1e6))    
