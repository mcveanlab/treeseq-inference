#!/usr/bin/env python3
"""
Take an .hdf5 file from vcf2tsinfer.py and run tsinfer on it.

"""
import sys
import os
import argparse

import h5py
import numpy as np

script_path = __file__ if "__file__" in locals() else "./dummy.py"
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','tsinfer')) # use the local copy of tsinfer in preference to the global one

import tsinfer

parser = argparse.ArgumentParser(description='Take an .hdf5 file from vcf2tsinfer.py and run tsinfer on it.')
parser.add_argument('infile', 
                    help='an hdf5 file output by vcf2tsinfer.py')
parser.add_argument('outfile',
                    help='the output file, in hdf5 treesequence format (see msprime)')
parser.add_argument('--restrict_sites', '-s', type=int, default=1000,
                    help='only use s sites from the total (set to 0 to use all)')
parser.add_argument('--start_at_site', '-start', type=int, default=0,
                    help='omit this many sites at the start')

args = parser.parse_args()

select = slice(args.start_at_site, args.restrict_sites or None) #if 0, select all
with h5py.File(args.infile, "r") as f:
    data = f['data']
    inferred_ts = tsinfer.infer(
        samples   = data['variants'][()][:,select], 
        positions = data['position'][select], 
        length    = max(data['position'][select])+1,
        recombination_rate = 1e-8, 
        error_rate = 1e-50,
        num_threads=10,
        show_progress=True)
    inferred_ts.dump(args.outfile)