#!/usr/bin/env python3
"""
Tests for the illustration code.
"""
import sys
import os
import string

import numpy as np
script_path = __file__ if "__file__" in locals() else "./dummy.py"
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','msprime')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(script_path)),'..','tsinfer')) # use the local copy of tsinfer in preference to the global one
    

import msprime
import tsinfer

rho = 2
ts = msprime.simulate(
    5, mutation_rate=5, recombination_rate=rho, random_seed=6)

S = np.zeros((ts.sample_size, ts.num_sites), dtype="u1")
for variant in ts.variants():
    S[:, variant.index] = variant.genotypes

sites = [mut.position for mut in ts.mutations()]
panel = tsinfer.ReferencePanel(S, sites, ts.sequence_length, rho=rho)
P, mutations = panel.infer_paths(num_workers=1)



labels = np.arange(P.shape[0])[~exclude]
corresponding_letters = [(v,string.ascii_lowercase[i]) for i,v in enumerate(labels)]
corresponding_letters.append((-1,"*"))

plot = tsinfer.Illustrator(panel, new_P)
plot.run(0, 'plot{:02d}'.format(0), corresponding_letters)
#for s in range(ts.num_sites):
#    plot.run(s, 'plot{:02d}'.format(s))
        