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
parser.add_argument('infile',
                    help='a path to the 1000G.samples file or 1000G.ancestors file. If a sample file, the ancestors file is created and saved')
parser.add_argument("--length-scale", "-X", choices=['linear', 'log'], default="linear",
    help='Should we transform lengths when plotting')
args = parser.parse_args()


data = tsinfer.load(args.infile)
fname = os.path.splitext(args.infile)[0]
if data.__class__ == tsinfer.formats.SampleData:
    #must create ancestors
    anc = tsinfer.generate_ancestors(
        data, 
        path=fname + ".ancestors",
        progress_monitor=tsinfer.cli.ProgressMonitor(generate_ancestors=True))
else:
    anc=data

frequency = anc.ancestors_time[:] + 1
positions = np.append(anc.sites_position[:], anc.sequence_length)
lengths_by_pos = (positions[anc.ancestors_end[:]]-positions[anc.ancestors_start[:]])/1000


df = pd.DataFrame({
    'l':lengths_by_pos,
    'f':frequency, 
    'nsites': [len(x) for x in anc.ancestors_focal_sites[:]]})

df_all = pd.DataFrame(
    {'lengths_per_site': np.repeat(df.l.values, df.nsites.values),
    'f': np.repeat(df.f.values, df.nsites.values),
    't':1
    }).sort_values(by=['f'])
df_all['pos'] = range(df_all.shape[0])
df_all['mean_pos'] = np.repeat(df_all.groupby('f').mean().pos.values, df_all.groupby('f').sum().t.values)
df_all['width'] = np.repeat(df_all.groupby('f').sum().t.values, df_all.groupby('f').sum().t.values)

mean_by_anc_time = df.iloc[df['nsites'].nonzero()].groupby('f', sort=True).mean()
median_by_anc_time = df.iloc[df['nsites'].nonzero()].groupby('f', sort=True).median()
sum_by_anc_time = df.iloc[df['nsites'].nonzero()].groupby('f', sort=True).sum()

fpos = np.insert(np.cumsum(sum_by_anc_time['nsites']).values, 0, 0)

fig = plt.figure(figsize=(10,10), dpi=100)
x_jittered = df_all.mean_pos.values+np.random.uniform(-df_all.width.values*9/20, df_all.width.values*9/20, len(df_all.mean_pos.values))
#plot with jitter
plt.scatter(x_jittered, df_all.lengths_per_site.values, marker='.', s=72./fig.dpi, alpha=0.05, color="black")
max_mean_y = np.max(mean_by_anc_time['l'][0:mean_by_anc_time.shape[0]//2]) #only use younger ancestors (we have a few very long old ancestors)
if args.length_scale == "log":
    plt.ylim(0.1, max_mean_y**2) 
    plt.yscale("log")
else:
    plt.ylim(1, max_mean_y*2) 

ax = plt.gca()
ax.step(fpos[:-1], mean_by_anc_time.l, label="Mean", where='post', color="orange")
ax.step(fpos[:-1], median_by_anc_time.l, label="Median", where='post', color="darkorange", linestyle=":")
ax.set_xlim(xmin=0)
ax.tick_params(axis='x', which="major", length=0)
ax.set_xticklabels('', minor=True)
ax.set_xticks(fpos[:-1], minor=True)
ax.set_xticks(fpos[:-1]+np.diff(fpos)/2)
ax.set_xticklabels(np.where(
    np.isin(mean_by_anc_time.index, np.array([1,2,3,4,5,6,10,50,1000, 5000])), 
    mean_by_anc_time.index,
    ""))
plt.xlabel("Prevalance (frequency), as ancestors_time+1 (~youngest to oldest)")
plt.ylabel("Length" + ("(kb)" if args.physical_length else "(# sites)"))
plt.legend(loc='upper center')
plt.savefig("{}_{}_length-time.png".format(
    os.path.basename(fname),"kb" if args.physical_length else "sites"), dpi=100)
plt.clf()