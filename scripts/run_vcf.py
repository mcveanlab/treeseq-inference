#!/usr/bin/env python3
"""
Take an .hdf5 file from vcf2tsinfer.py and run tsinfer on it.

"""
import sys
import os

import h5py
import numpy as np

script_path = __file__ if "__file__" in locals() else "./dummy.py"
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','tsinfer')) # use the local copy of tsinfer in preference to the global one

import tsinfer


#To do - use argparse module 

n_variants = 1000 #take the first 1000 variants

with h5py.File(sys.argv[1], "r") as f:
    data = f['data']
    inferred_ts = tsinfer.infer(
        samples   = data['variants'][()][:,:n_variants], 
        positions = data['position'][:n_variants], 
        length    = max(data['position'][:n_variants])+1,
        recombination_rate = 1e-8, 
        error_rate = 1e-50,
        num_threads=10,
        show_progress=True)
    inferred_ts.dump(sys.argv[2])