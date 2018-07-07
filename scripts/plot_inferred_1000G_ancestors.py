#!/usr/bin/env python3
"""
Plot ancestor characteristics for 1000G
"""
import sys
import os
import argparse

import numpy as np
import pandas as pd
import matplotlib as mp
# Force matplotlib to not use any Xwindows backend.
mp.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(1,os.path.join(sys.path[0],'..','tsinfer'))
import tsinfer

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

parser = argparse.ArgumentParser(description='Plot ancestors generated from the 1000G data.')
parser.add_argument('samples_file',
                    help='a path to the 1000G.samples file')
args = parser.parse_args()

sample_data = tsinfer.load(args.samples_file)
anc = tsinfer.generate_ancestors(
    sample_data, 
    progress_monitor=tsinfer.cli.ProgressMonitor(generate_ancestors=True))
frequency = anc.ancestors_time[:]
lengths = anc.ancestors_end[:]-anc.ancestors_start[:]

df = pd.DataFrame(
    {'l':lengths, 'f':frequency, 'nsites': [len(x) for x in anc.ancestors_focal_sites[:]]}
    )

lengths_by_anc_time = df.iloc[df['nsites'].nonzero()].groupby('f', sort=True).mean()


plt.semilogy(lengths_by_anc_time.index, lengths_by_anc_time['l'], label="Estimated ancestors")
plt.ylim(0, np.max(lengths_by_anc_time['l'][0:lengths_by_anc_time.shape[0]//2])) #set to max of younger ancestors (we have a few very long old ancestors
plt.xlabel("Ancestors_time (=freq, youngest to oldest)")
plt.ylabel("Length")
plt.legend()
plt.savefig("{}_length-time.png".format(samples_file.replace(".samples","")))
plt.clf()
