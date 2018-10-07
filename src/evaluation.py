#!/usr/bin/env python3
"""
Code to run simulations, inference methods and generate all plots
in the paper.

Run as e.g.

./evaluation.py setup metrics_by_mutation_rate -P
./evaluation.py infer metrics_by_mutation_rate -P -p 30 -t 8 #this may take a long time
./evaluation.py figure kc_rooted_by_mutation_rate

"""

import argparse
import collections
import filecmp
import glob
import itertools
import json
import logging
import multiprocessing
import os
import sys
import random
import shutil
import signal
import statistics
import subprocess
import tempfile
import time
import math
import gzip

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.backends.backend_pdf
import pandas as pd
import tqdm

import msprime
import _msprime
import tsinfer

# The following are all contained in files local to this directory
import ts_extras
import ts_fastARG
import ts_ARGweaver
import ts_RentPlus
import ARG_metrics

fastARG_executable = os.path.join('tools','fastARG','fastARG')
ARGweaver_executable = os.path.join('tools','argweaver','bin','arg-sample')
smc2arg_executable = os.path.join('tools','argweaver','bin','smc2arg')
RentPlus_executable = os.path.join('tools','RentPlus','RentPlus.jar')
tsinfer_executable = os.path.join('src','run_tsinfer.py')

#monkey-patch nexus saving/writing into msprime/tskit
msprime.TreeSequence.write_nexus_trees = ts_extras.write_nexus_trees
msprime.TreeSequence.save_nexus_trees = ts_extras.save_nexus_trees

FASTARG = "fastARG"
ARGWEAVER = "ARGweaver"
RENTPLUS = "RentPlus"
TSINFER = "tsinfer"

#names for optional columns in the dataset
ERROR_COLNAME = 'error_param'
SUBSAMPLE_COLNAME = 'subsample_size'
SIMTOOL_COLNAME = 'sim_tool' #if column does not exist, default simulation tool = 'msprime'
MUTATION_SEED_COLNAME = 'mut_seed'
SELECTION_COEFF_COLNAME = 'selection_coefficient'
DOMINANCE_COEFF_COLNAME = 'dominance_coefficient'
SELECTED_FREQ_COLNAME = 'output_frequency'
SELECTED_POSTGEN_COLNAME = 'output_after_generations'

#bit flags for metrics, first bit indicates which locations are used, 
# second bit indicates whether to break polytomies randomly
METRICS_LOCATION_ALL      = 0 #metrics are averaged over all positions in the genome
METRICS_LOCATION_VARIANTS = 2**0 #metrics only averaged over variant sites
METRICS_POLYTOMIES_LEAVE  = 0
METRICS_POLYTOMIES_BREAK  = 2**1

if sys.version_info[0] < 3:
    raise Exception("Python 3 only")

save_stats = dict(
    cpu = "cputime",
    mem =  "memory",
    n_edges = "edges",
    ts_filesize = "ts_filesize"
)

def latex_float(f):
    """
    Return an exponential number in nice LaTeX form. 
    In titles etc for plots this needs to be encased in $ $ signs, and r'' strings used
    """
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

def nanblank(val):
    """hack around a horrible pandas syntax, which puts nan instead of blank strings"""
    return "" if pd.isnull(val) else val

def ts_has_non_singleton_variants(ts):
    """
    Check if there are any doubletons or higher in a treeseq
    which have not been fixed
    (required for meaningful inference)
    """
    for v in ts.variants():
        if 1 < np.sum(v.genotypes) < ts.num_samples:
            return True
    return False


def make_errors_genotype_model(g, error_probs):
    """
    Given an empirically estimated error probability matrix, resample for a particular
    variant. Determine variant frequency and true genotype (g0, g1, or g2),
    then return observed genotype based on row in error_probs with nearest
    frequency. Treat each pair of alleles as a diploid individual. 
    """
    w = np.copy(g)
    
    #Make diploid (iterate each pair of alleles)
    genos=[(w[i],w[i+1]) for i in range(0,w.shape[0],2)]
    
    #Record the true genotypes
    g0 = [i for i, x in enumerate(genos) if x == (0,0)]
    g1a = [i for i, x in enumerate(genos) if x == (1,0)]
    g1b = [i for i, x in enumerate(genos) if x == (0,1)]
    g2 = [i for i, x in enumerate(genos) if x == (1,1)]
    

    for idx in g0:
        result=[(0,0),(1,0),(1,1)][np.random.choice(3,
            p=error_probs[['p00','p01','p02']].values[0])]
        if result == (1,0):
            genos[idx]=[(0,1),(1,0)][np.random.choice(2)]
        else:
            genos[idx] = result
    for idx in g1a:
        genos[idx]=[(0,0),(1,0),(1,1)][np.random.choice(3,
            p=error_probs[['p10','p11','p12']].values[0])]
    for idx in g1b:
        genos[idx]=[(0,0),(0,1),(1,1)][np.random.choice(3,
            p=error_probs[['p10','p11','p12']].values[0])]
    for idx in g2:
        result=[(0,0),(1,0),(1,1)][np.random.choice(3,
            p=error_probs[['p20','p21','p22']].values[0])]
        if result == (1,0):
            genos[idx]=[(0,1),(1,0)][np.random.choice(2)]
        else:
            genos[idx] = result
            
    return(np.array(sum(genos, ())))


def mk_sim_name(sample_size, Ne, length, recombination_rate, mutation_rate, seed, 
    mut_seed=None, directory=None, tool="msprime",
    s=None, h=None, freq=None, post_gens=None):
    """
    Create a filename for saving an msprime simulation (without extension)
    Other functions add error rates and/or a subsample sizes to the name
    """
    #format mut rate & recomb rate to print non-exponential notation without
    # trailing zeroes. 12 dp should be ample for these rates
    rho = "{:.12f}".format(float(recombination_rate)).rstrip('0')
    mu = "{:.12f}".format(float(mutation_rate)).rstrip('0')
    file = "{}-n{}_Ne{}_l{}_rho{}_mu{}-gs{}".format(tool, int(sample_size), \
        Ne if isinstance(Ne, str) else float(Ne), int(length), rho, mu, int(seed))
    if mut_seed is not None:
        file += "_ms{}".format(int(mut_seed))
    if s is not None:
        file += "_s{}".format(s)
    if h is not None:
        file += "_h{}".format(h)
    if freq is not None:
        file += "_f{}".format(freq)
        if post_gens is not None:
            file += "+{}".format(post_gens)
    if directory is None:
        return file
    else:
        return os.path.join(directory,file)

def mk_sim_name_from_row(row, directory=None, error_col=ERROR_COLNAME, subsample_col=SUBSAMPLE_COLNAME):
    """
    If error_col and subsample_col are None, this is the same as mk_sim_name()
    but filled out using data from a row. If error_col is a string which exists as
    a column name in the row then add_error_param_to_name() is also called, using
    the error rate specified in that column.
    """
    #fill out optional colnames (these might not exist)
    tool = row.get(SIMTOOL_COLNAME, 'msprime') #default to 'msprime' if tool not specified in row
    mut_seed = row.get(MUTATION_SEED_COLNAME)
    s = row.get(SELECTION_COEFF_COLNAME)
    h = row.get(DOMINANCE_COEFF_COLNAME)
    freq = row.get(SELECTED_FREQ_COLNAME)
    post_gens = row.get(SELECTED_POSTGEN_COLNAME)

    name = mk_sim_name(row.sample_size, row.Ne, row.length, row.recombination_rate,
        row.mutation_rate, row.seed, mut_seed, directory,
        tool=tool, s=s, h=h, freq=freq, post_gens = post_gens)
    #in some cases the original simulation (or results from the simulation) have been
    # modified and we need a name that reflects that modification
    if subsample_col and subsample_col in row:
        name = add_subsample_param_to_name(name, row[subsample_col])
    if error_col and error_col in row:
        name = add_error_param_to_name(name, row[error_col])
    return(name)

def add_subsample_param_to_name(sim_name, subsample_size=None):
    """
    Mark a filename as containing only a subset of the samples of the full sim
    """
    if subsample_size is not None and not pd.isnull(subsample_size):
        if sim_name.endswith("+") or sim_name.endswith("-"):
            #this is the first param
            return sim_name + "max{}".format(int(subsample_size))
        else:
            return sim_name + "_max{}".format(int(subsample_size))
    else:
        return sim_name

def add_error_param_to_name(sim_name, error_param=None):
    """
    Append the error param to the simulated tree sequence filename.
    Only relevant for files downstream of the step where sequence error is added
    """
    if error_param is not None and not pd.isnull(error_param):
        if sim_name.endswith("+") or sim_name.endswith("-"):
            #this is the first param
            return sim_name + "_err{}".format(error_param)
        else:
            #this is the first param
            return sim_name + "err{}".format(error_param)
    else:
        return sim_name

def construct_fastarg_name(sim_name, seed, directory=None):
    """
    Returns a fastARG inference filename (without file extension),
    based on a simulation name
    """
    d,f = os.path.split(sim_name)
    return os.path.join(d,'+'.join(['fastarg', f, "fs"+str(int(seed))]))

def fastarg_name_from_row(row, sim_dir):
    """
    return the fa name based on an sim specified by row
    """
    return construct_fastarg_name(mk_sim_name_from_row(row, sim_dir),
                                  seed=row.seed)


def construct_argweaver_name(sim_name, burnin, ntimes, seed, iteration_number=None):
    """
    Returns an ARGweaver inference filename (without file extension),
    based on a simulation name. The iteration number (used in .smc and .nex output)
    is usually added by the ARGweaver `arg-sample` program,
    in the format .10, .20, etc. (we append an 'i' too, giving
    'i.10', 'i.100', etc
    """
    d,f = os.path.split(sim_name)
    suffix = "burn"+str(int(burnin)) + "nt"+str(int(ntimes)) + "ws"+str(int(seed))
    if iteration_number is not None:
        suffix += "_i."+ str(int(iteration_number))
    return os.path.join(d,'+'.join(['aweaver', f, suffix]))

def argweaver_names_from_row(row, sim_dir):
    """
    return multiple argweaver names based on an sim specified by row
    there is one name per argweaver iteration listed in row.ARGweaver_iterations
    """
    return [construct_argweaver_name(mk_sim_name_from_row(row, sim_dir),
                                     burnin=row.ARGweaver_burnin, ntimes=ARGweaver_ntimes,
                                     seed=row.seed, iteration_number=it)
                for it in nanblank(row.ARGweaver_iterations).split(",") if it]

def construct_rentplus_name(sim_name):
    """
    Returns an RentPlus inference filename (without file extension),
    based on a simulation name.
    """
    d,f = os.path.split(sim_name)
    return os.path.join(d,'+'.join(['rentpls', f, ""]))

def rentplus_name_from_ts_row(row, sim_dir):
    """
    return the rentplus name based on a simulated tree sequence specified by row
    """
    return construct_rentplus_name(mk_sim_name_from_row(row, sim_dir))

def construct_tsinfer_name(sim_name, subsample_size, error_param=None):
    """
    Returns a TSinfer filename.
    If the file is a subset of the original, this can be added to the
    basename in this function, or later using the
    add_subsample_param_to_name() routine.
    """
    d,f = os.path.split(sim_name)
    suffix = "" if error_param is None else "err{}".format(error_param)
    name = os.path.join(d,'+'.join(['tsinfer', f, suffix]))
    if subsample_size is not None and not pd.isnull(subsample_size):
        name = add_subsample_param_to_name(name, subsample_size)
    return name

def time_cmd(cmd, stdout=sys.stdout):
    """
    Runs the specified command line (a list suitable for subprocess.call)
    and writes the stdout to the specified file object.
    """
    if sys.platform == 'darwin':
        #on OS X, install gtime using `brew install gnu-time`
        time_cmd = "/usr/local/bin/gtime"
    else:
        time_cmd = "/usr/bin/time"
    full_cmd = [time_cmd, "-f%M %S %U"] + cmd
    with tempfile.TemporaryFile() as stderr:
        exit_status = subprocess.call(full_cmd, stderr=stderr, stdout=stdout)
        stderr.seek(0)
        if exit_status != 0:
            raise ValueError(
                "Error running '{}': status={}:stderr{}".format(
                    " ".join(cmd), exit_status, stderr.read()))

        split = stderr.readlines()[-1].split()
        # From the time man page:
        # M: Maximum resident set size of the process during its lifetime,
        #    in Kilobytes.
        # S: Total number of CPU-seconds used by the system on behalf of
        #    the process (in kernel mode), in seconds.
        # U: Total number of CPU-seconds that the process used directly
        #    (in user mode), in seconds.
        max_memory = int(split[0]) * 1024
        system_time = float(split[1])
        user_time = float(split[2])
    return user_time + system_time, max_memory

def ARGmetric_params_from_row(row):
    """
    Create some ARGmetric params (see ARGmetrics.py) from a row
    We hack the same polytomy seed as inference seed

    This stupidly has to be defined at the top level, since it is
    part of the function passed in to a multiprocessing Pool, and
    hence needs to be 'pickle'able :(
    """
    return {'make_bin_seed':row.seed, 'reps':row.tsinfer_biforce_reps}


class InferenceRunner(object):
    """
    Class responsible for running a single inference tool and returning results for
    the dataframe. Should create results files that are bespoke for each tool, and
    also for each tool, convert these into nexus files that can be used for comparing
    metrics.
    metric_params is a list of possible parameter combinations for which to run metrics.
    If it is an empty array, it indicates metric calcs should be turned off
        """
    def __init__(
            self, tool, row, simulations_dir, num_threads,
            metric_params, polytomy_reps):
        self.tool = tool
        self.row = row
        self.num_threads = num_threads
        self.metric_params = metric_params
        self.polytomy_reps = polytomy_reps
        #the original simulation file does not include error or subsampling
        self.orig_sim_fn = mk_sim_name_from_row(row, simulations_dir, subsample_col=None, error_col=None)
        #the known tree file to compare inference against: is subset but no err
        self.cmp_fn = mk_sim_name_from_row(
            row, simulations_dir, subsample_col='restrict_sample_size_comparison', error_col=None)
        #sample_fn is used for sample data from the simulation, so must include error,
        #and (if used) the subsample size
        self.sample_fn = mk_sim_name_from_row(row, simulations_dir)
        # This should be set by the run_inference methods. Append ".nex" to
        # get the nexus files, or ".hdf5" to get the TS files
        self.inferred_filenames = None

    def run(self, metrics_only):
        logging.debug("parameters = {}".format(self.row.to_dict()))
        if self.tool == TSINFER:
            ret = self.__run_tsinfer(skip_infer = metrics_only)
        elif self.tool == FASTARG:
            ret = self.__run_fastARG(skip_infer = metrics_only)
        elif self.tool == ARGWEAVER:
            ret = self.__run_ARGweaver(skip_infer = metrics_only)
        elif self.tool == RENTPLUS:
            ret = self.__run_RentPlus(skip_infer = metrics_only)
        else:
            raise ValueError("unknown tool {}".format(self.tool))
        if self.inferred_filenames is None:
            logging.info("No inferred tree files so metrics skipped for {} row {} = {}".format(
                self.tool, int(self.row[0]), ret))
        else:
            for metric in self.metric_params:
            #NB Jerome thinks it may be clearer to have get_metrics() return a single set of metrics
            #rather than an average over multiple inferred nexus files, and do the averaging in python
                if metric & METRICS_LOCATION_VARIANTS:
                    #get positions from the samples store, for use in metric calcs
                    try:
                        positions = tsinfer.SampleData.load(path=self.cmp_fn + ".samples").sites_position[:].tolist()
                    except FileNotFoundError:
                        #no such file exists, could be a case with no sites
                        continue
                else:
                    positions = None
                source_nexus_file = self.cmp_fn + ".nex"
                inferred_nexus_files = [fn + ".nex" for fn in self.inferred_filenames]
                #here we should create a separate set of metrics for tsinfer with and without polytomy breaking
                #we should check if it is TSINFER, and then prepend '' for the default metric and ''
                if metric & METRICS_POLYTOMIES_BREAK:
                    metrics = []
                    for i in range(self.polytomy_reps):
                        metrics.append(ARG_metrics.get_metrics(
                            source_nexus_file, inferred_nexus_files, variant_positions = positions,
                            randomly_resolve_inferred = int(self.row.seed)+i*11))
                    metrics = pd.DataFrame(metrics).mean().to_dict()
                else:
                    metrics = ARG_metrics.get_metrics(
                        source_nexus_file, inferred_nexus_files, variant_positions = positions)
                ret.update({(str(metric) + "_" + k):v for k,v in metrics.items()})
        logging.debug("returning infer results for {} row {} = {}".format(
            self.tool, int(self.row[0]), ret))
        #flatten the return value and prepend the tool name
        row = {(self.tool + "_" + k):v for k,v in ret.items()}
        return row

    def __run_tsinfer(self, skip_infer=False):
        #default to no subsampling
        subsample_size = getattr(self.row,'subsample_size', None)
        restrict_sample_size_comparison = getattr(self.row,'restrict_sample_size_comparison', None)
        #construct filenames - these can be used even if inference does not occur
        samples_fn = self.sample_fn + ".samples"
        out_fn = construct_tsinfer_name(self.sample_fn,
            restrict_sample_size_comparison)
        self.inferred_filenames = [out_fn]
        if skip_infer:
            return {}
        #Now perform the inference
        time = memory = fs = counts = inferred_ts = None
        logging.debug("loading samples for ts inference from {}".format(
            samples_fn))
        try:
            inferred_ts, time, memory = self.run_tsinfer(
                samples_fn, self.row.length, self.num_threads,
                #uncomment below to inject real ancestors - will need adjusting for subsampling
                #inject_real_ancestors_from_ts_fn = self.orig_sim_fn + ".hdf5",
                )
            if restrict_sample_size_comparison is not None:
                if len(self.metric_params):
                    with open(construct_tsinfer_name(self.sample_fn, None) + ".nex", "w+") as out:
                        tree_labels_between_variants=(True if subsample_size is None else False)
                        inferred_ts.write_nexus_trees(
                            out, tree_labels_between_variants=tree_labels_between_variants)
                inferred_ts = inferred_ts.simplify(list(range(restrict_sample_size_comparison)))
                
            inferred_ts.dump(out_fn + ".inferred.hdf5")
            fs = os.path.getsize(out_fn + ".inferred.hdf5")
            if len(self.metric_params):
                with open(self.inferred_filenames[0] + ".nex", "w+") as out:
                    #For subsampled trees, the locations of trees along the
                    #genome (the breakpoints) may not correspond to variant positions
                    #on the subsampled trees: indeed, there may be multiple trees
                    #between two adjacent variants, so we cannot reposition breakpoints
                    #between the nearest left and right variant positions
                    tree_labels_between_variants=(True if subsample_size is None else False)
                    inferred_ts.write_nexus_trees(
                        out, tree_labels_between_variants=tree_labels_between_variants)
            unique, counts = np.unique(np.array([e.parent for e in inferred_ts.edges()], dtype="u8"), return_counts=True)
        except ValueError as e:
            if "No inference sites" in str(e):
                logging.warning("No inference sites in {}. Skipping".format(samples_fn))
                self.inferred_filenames = None
            elif "No samples file" in  str(e):
                logging.warning("No samples file at {}. Skipping".format(samples_fn))
                self.inferred_filenames = None
            else:
                raise            
        return  {
            save_stats['cpu']: time,
            save_stats['mem']: memory,
            save_stats['n_edges']: None if counts is None else np.sum(counts),
            save_stats['ts_filesize']: fs,
            'mean_polytomy': None if counts is None else np.mean(counts),
            'var_polytomy': None if counts is None else np.var(counts),
            'max_polytomy': None if counts is None else np.max(counts)
        }

    def __run_fastARG(self, skip_infer=False):
        inference_seed = self.row.seed  # TODO do we need to specify this separately?
        if skip_infer:
            return {}
        infile = self.sample_fn + ".hap"
        time = memory = edges = fs = None
        logging.debug("reading: {} for fastARG inference".format(infile))
        try:
            inferred_ts, time, memory = self.run_fastarg(infile, self.row.length, inference_seed)
            edges = inferred_ts.num_edges
            self.inferred_filenames = [construct_fastarg_name(self.sample_fn, inference_seed)]
            inferred_ts.dump(self.inferred_filenames[0] + ".hdf5")
            fs = os.path.getsize(self.inferred_filenames[0] + ".hdf5")
            if len(self.metric_params):
                for fn in self.inferred_filenames:
                    with open(fn + ".nex", "w+") as out:
                        inferred_ts.write_nexus_trees(out)
        except ValueError as e:
            if "0 samples;" in str(e):
                logging.warning("No samples in {}. Skipping".format(infile))
                self.inferred_filenames = None
            else:
                raise
        return {
            save_stats['cpu']: time,
            save_stats['mem']: memory,
            save_stats['n_edges']: edges,
            save_stats['ts_filesize']: fs,
        }

    def __run_RentPlus(self, skip_infer=False):
        self.inferred_filenames = [construct_rentplus_name(self.sample_fn)]
        if skip_infer:
            return {}
        infile = self.sample_fn + ".dat"
        time = memory = None
        logging.debug("reading: {} for RentPlus inference".format(infile))
        try:
            treefile, num_tips, time, memory = self.run_rentplus(infile, self.row.length)
            if len(self.metric_params):
                for fn in self.inferred_filenames:
                    with open(fn + ".nex", "w+") as out:
                        ts_RentPlus.RentPlus_trees_to_nexus(treefile, out, self.row.length, num_tips)
        except ValueError as e:
            self.inferred_filenames = None
            #allow the RentPlus error described at https://github.com/mcveanlab/treeseq-inference/issues/45
            if "main.Main.findCompatibleRegion(Main.java:660)" in str(e):
                logging.warning("RENTplus bug hit"\
                    " (https://github.com/SajadMirzaei/RentPlus/issues/5),"\
                    " aborting this replicate.")
            elif "main.Main.makeDistanceMatrices(Main.java:228)":
                logging.warning("RENTplus bug hit"\
                    " (https://github.com/SajadMirzaei/RentPlus/issues/6),"\
                    " aborting this replicate.")
            else:
                raise
        return {
            save_stats['cpu']: time,
            save_stats['mem']: memory,
        }


    def __run_ARGweaver(self, skip_infer=False):
        """
        The default values are used if there are no columns which specify the
        number of burn in iterations or number of discretised timesteps
        """
        #set some AW params
        n_out_samples=10
        sample_step=20
        #these ones can be overridden
        default_ARGweaver_burnin=1000
        default_ARGweaver_ntimes=20

        #this is a hack when we specify a complex model which does not have an explicit Ne
        #just set the Ne to something likely (HACK!!)
        ARGweaver_Ne = 10000 if isinstance(self.row.Ne, str) else str(self.row.Ne)

        inference_seed = self.row.seed  # TODO do we need to specify this separately?
        burnin = getattr(self.row,'ARGweaver_burnin', default_ARGweaver_burnin)
        n_timesteps = getattr(self.row,'ARGweaver_ntimes', default_ARGweaver_ntimes)
        self.inferred_filenames = []
        iteration_ids = []
        out_fn = construct_argweaver_name(self.sample_fn, burnin, n_timesteps, inference_seed)
        if skip_infer:
            #must get the iteration IDs off the row
            if self.row.ARGweaver_iterations:
                for it in self.row.ARGweaver_iterations.split(","):
                    self.inferred_filenames.append(
                        construct_argweaver_name(self.sample_fn, burnin, n_timesteps, inference_seed, it))
            return {}
        infile = self.sample_fn + ".sites"
        time = memory = None
        filesizes = []
        edges = []
        iteration_ids = []
        logging.debug("reading: {} for ARGweaver inference".format(infile))
        try:
            iteration_ids, stats_file, time, memory = self.run_argweaver(
                infile, ARGweaver_Ne, self.row.recombination_rate, self.row.mutation_rate,
                out_fn, inference_seed, n_out_samples, sample_step, burnin, n_timesteps,
                verbose = logging.getLogger().isEnabledFor(logging.DEBUG))
            for it in iteration_ids:
                base = construct_argweaver_name(self.sample_fn, burnin, n_timesteps, inference_seed, it)
                self.inferred_filenames.append(base)
                if skip_infer==False:
                    if len(self.metric_params):
                        with open(base + ".nex", "w+") as out:
                            ts_ARGweaver.ARGweaver_smc_to_nexus(base+".smc.gz", out)
                    try:
                        with open(base+".TSnodes", "w+") as ts_nodes, \
                                open(base+".TSedges", "w+") as ts_edges:
                            ts_ARGweaver.ARGweaver_smc_to_ts_txts(
                                smc2arg_executable, base, ts_nodes, ts_edges)
                            inferred_ts = msprime.load_text(
                                nodes=ts_nodes, edges=ts_edges).simplify()
                            inferred_ts.dump(base + ".hdf5")
                            filesizes.append(os.path.getsize(base + ".hdf5"))
                            edges.append(inferred_ts.num_edges)
                    except ts_ARGweaver.CyclicalARGError as e:
                        logging.warning("Cyclical ARG Exception when converting {}: {}".format(
                            base + ".msp", e))
        except ValueError as e:
            if 'src/argweaver/sample_thread.cpp:517:' in str(e):
                logging.warning("Hit argweaver bug " \
                "https://github.com/mcveanlab/treeseq-inference/issues/25" \
                " for {}. Skipping".format(out_fn))
            elif "Assertion `trans[path[i]] != 0.0' failed" in str(e):
                logging.warning("Hit argweaver bug " \
                "https://github.com/mcveanlab/treeseq-inference/issues/42" \
                " for {}. Skipping".format(out_fn))
            else:
                raise

        return {
            save_stats['cpu']: time,
            save_stats['mem']: memory,
            save_stats['n_edges']: statistics.mean(edges) if len(edges) else None,
            save_stats['ts_filesize']: statistics.mean(filesizes) if len(filesizes) else None,
            "iterations": ",".join(iteration_ids),
        }

    @staticmethod
    def run_tsinfer(sample_fn, length,
        num_threads=1, inject_real_ancestors_from_ts_fn=None, rho=None, error_probability=None):
            with tempfile.NamedTemporaryFile("w+") as ts_out:
                cmd = [sys.executable, tsinfer_executable, sample_fn, "--length", str(int(length))]
                cmd += ["--threads", str(num_threads), ts_out.name]
                if inject_real_ancestors_from_ts_fn:
                    logging.debug("Injecting real ancestors constructed from {}".format(
                        inject_real_ancestors_from_ts_fn))
                    cmd.extend(["--inject-real-ancestors-from-ts", inject_real_ancestors_from_ts_fn])
                cpu_time, memory_use = time_cmd(cmd)
                ts_simplified = msprime.load(ts_out.name)
            return ts_simplified, cpu_time, memory_use

    @staticmethod
    def run_fastarg(file_name, seq_length, seed):
        with tempfile.NamedTemporaryFile("w+") as fa_out, \
                tempfile.NamedTemporaryFile("w+") as tree, \
                tempfile.NamedTemporaryFile("w+") as muts:
            cmd = ts_fastARG.get_cmd(fastARG_executable, file_name, seed)
            cpu_time, memory_use = time_cmd(cmd, fa_out)
            logging.debug("ran fastarg for seq length {} [{} s]: '{}'".format(seq_length, cpu_time, cmd))
            var_pos = ts_fastARG.variant_positions_from_fastARGin_name(file_name)
            inferred_ts = ts_fastARG.fastARG_out_to_ts(fa_out, var_pos, seq_len=seq_length)
            return inferred_ts, cpu_time, memory_use

    @staticmethod
    def run_rentplus(file_name, seq_length):
        """
        runs RentPlus, returning the output filename, the total CPU, & max mem
        must check here if we are using 0..1 positions (infinite sites) or integers
        """
        haplotype_lines = 0
        integer_positions = True
        with open(file_name, "r+") as infile:
            for pos in next(infile).split():
                try:
                    dummy = int(pos)
                except ValueError:
                    integer_positions = False
            for line in infile:
                if line.rstrip():
                    haplotype_lines += 1
        cmd = ["java", "-jar", RentPlus_executable]
        cmd += [file_name] if integer_positions else ['-l', seq_length, file_name]
        with tempfile.NamedTemporaryFile("w+") as script_output:
            cpu_time, memory_use = time_cmd(cmd, script_output)
        logging.debug("ran RentPlus for {} haplotypes with seq length {} [{} s]: '{}'".format(
            haplotype_lines, seq_length, cpu_time, cmd))
        #we cannot back-convert RentPlus output to treeseq form - just return the trees file
        assert os.path.isfile(file_name + ".trees"), 'No trees file created when running Rent+'
        #we might also want to look at TMRCAs (file_name + '.Tmrcas')
        return file_name + ".trees", haplotype_lines, cpu_time, memory_use

    @staticmethod
    def run_argweaver(
            sites_file, Ne, recombination_rate, mutation_rate, path_prefix, seed,
            MSMC_samples, sample_step, burnin_iterations, ntimes, verbose=False):
        """
        this produces a whole load of .smc files labelled <path_prefix>i.0.smc,
        <path_prefix>i.10.smc, etc, which we then convert into sts format

        Returns the iteration numbers ('0', '10', '20' etc), the name of the
        statistics file, the total CPU time, and the max memory usage.
        """
        cpu_time = []
        memory_use = []
        burn_prefix = None
        exe = [ARGweaver_executable, '--sites', sites_file.name if hasattr(sites_file, "name") else sites_file,
               '--popsize', str(Ne),
               '--recombrate', str(recombination_rate),
               '--mutrate', str(mutation_rate),
               '--ntimes', str(ntimes),
               '--overwrite']
        if not verbose:
            exe += ['--quiet']
        if seed is not None:
            exe += ['--randseed', str(int(seed))]
        if burnin_iterations > 0:
            burn_in = str(int(burnin_iterations))
            burn_prefix = path_prefix+"_burn"
            logging.info("== Burning in ARGweaver MCMC using {} steps ==".format(burn_in))
            logging.debug("== ARGweaver burnin command is {} ==".format(" ".join(exe)))
            c, m = time_cmd(exe + ['--iters', burn_in,
                                   '--sample-step', burn_in,
                                   '--output', burn_prefix])
            cpu_time.append(c)
            memory_use.append(m)
            #if burn_in, read from the burn in arg file, rather than the original .sites
            exe += ['--arg', burn_prefix+"."+ burn_in +".smc.gz"]
        else:
            exe += ['--sites', sites_file]

        new_prefix = path_prefix + "_i" #we append a '_i' to mark iteration number
        iterations = int(sample_step * (MSMC_samples-1))
        exe += ['--output', new_prefix]
        exe += ['--iters', str(iterations)]
        exe += ['--sample-step', str(int(sample_step))]
        logging.info("== Running ARGweaver for {} steps to collect {} samples ==".format( \
            int(iterations), MSMC_samples))
        logging.debug("== ARGweaver command is {} ==".format(" ".join(exe)))
        c, m = time_cmd(exe)
        cpu_time.append(c)
        memory_use.append(m)

        smc_prefix = new_prefix + "." #the arg-sample program adds .iteration_num
        saved_iterations = [f[len(smc_prefix):-7] for f in glob.glob(smc_prefix + "*" + ".smc.gz")]
        new_stats_file_name = path_prefix+".stats"

        #concatenate all the stats together
        with open(new_stats_file_name, "w+") as stats:
            if burn_prefix:
                shutil.copyfileobj(open(burn_prefix + ".stats"), stats)
                print("\n", file=stats)
            shutil.copyfileobj(open(new_prefix + ".stats"), stats)
        #To Do: translate these to treesequence objects, so that e.g. edges can be calculated
        #as https://github.com/mdrasmus/argweaver/issues/20 is now closed
        return saved_iterations, new_stats_file_name, sum(cpu_time), max(memory_use)


def infer_worker(work):
    """
    Entry point for running a single inference task in a worker process.
    """
    tool, row, sims_dir, n_threads, metric_params, metrics_only, polytomy_reps = work
    runner = InferenceRunner(tool, row, sims_dir, n_threads, metric_params, polytomy_reps)
    result = runner.run(metrics_only)
    result[tool + '_completed'] = True
    return int(row[0]), tool, result


class Dataset(object):
    """
    A dataset is some collection of simulations and associated data.
    Results for datasets are stored in a general data_dir
    """
    data_dir = os.path.join(sys.path[0],"..","data")
    """
    Each dataset has a unique name. This is used as the prefix for the data
    file and raw_data_dir directory. Within this, replicate instances of datasets
    (each with a different RNG seed) are saved under the seed number
    """
    name = None

    """
    Two defaults that can be overridden on the command line
    """
    default_replicates = 10
    default_seed = 123

    """
    Which tools to run for this dataset, and which tree metric parameters to run for each tool
    (tree metrics calculate distances from the source tree seq to the inferred tree seq(s)).
    One tool can have multiple tree metric parameters, each of which gets calculated.
    tools_and_metrics can be overridden for particular datasets, and if we do not want to run 
    any metrics we can set the value for a tool to [].
    """
    tools_and_metrics = {
        ARGWEAVER : [METRICS_LOCATION_ALL | METRICS_POLYTOMIES_LEAVE],
        FASTARG:    [METRICS_LOCATION_ALL | METRICS_POLYTOMIES_LEAVE],
        RENTPLUS:   [METRICS_LOCATION_ALL | METRICS_POLYTOMIES_LEAVE],
        #tsinfer regularly produces polytomies, so need to check both methods here
        TSINFER:    [METRICS_LOCATION_ALL | METRICS_POLYTOMIES_LEAVE, METRICS_LOCATION_ALL | METRICS_POLYTOMIES_BREAK],
    }
    
    """
    If any metrics are calculated by randomly resolving polytomies (by setting
    one of the tools_and_metrics values to METRICS_POLYTOMIES_BREAK, this parameter gives
    the number of times the same tree will replicate the polytomy breaking process
    (the final metric will be a simple mean of these replicates)
    """
    random_resolve_polytomy_replicates = 10

    """
    These are the basic columns that record the simulation used (can be overridden)
    """
    sim_cols = [
        "replicate", "sample_size", "Ne", "length", "recombination_rate", "mutation_rate",
        ERROR_COLNAME, "edges", "n_trees", "n_sites", "seed"]

    extra_sim_cols = []

    """
    For a tidier csv file, we can exclude any of the save_stats values or ARGmetrics 
    columns using a regexp, e.g. 'rooted$'
    """
    exclude_colnames_matching = []

    """
    We store empirically determined error stats (e.g. from 1000 genomes platinum analysis
    in a separate csv file in the data_dir. We can use this filename as the label in the
    ERROR_COLNAME column in the results file, so we know which error matrix we have used
    """
    error_filename = "EmpiricalErrorPlatinum1000G.csv"

    def __init__(self):
        self.data_file = os.path.abspath(os.path.join(self.data_dir, self.name)) + ".csv"
        self.raw_data_dir = os.path.join(self.data_dir, "raw__NOBACKUP__", self.name)
        self.simulations_dir = os.path.join(self.raw_data_dir, "simulations")
        self.last_data_write_time = time.time()
        self.metric_prefixes = [t + "_" + str(bits) \
            for t, metrics in self.tools_and_metrics.items() \
                for bits in metrics]

    def load_data(self):
        self.data = pd.read_csv(self.data_file)

    def dump_data(self, write_index=False, force_flush=True):
        """
        Dumps data if it hasn't been written in the last 30 seconds. If force is true,
        write it out anyway.
        """
        now = time.time()
        if force_flush or (now - self.last_data_write_time) > 30:
            logging.info("Flushing data file")
            self.data.to_csv(self.data_file, index=write_index)
            self.last_data_write_time = time.time()

    #
    # Main entry points.
    #
    def setup(self, args):
        """
        Creates the dataframe and storage directories and then runs the initial
        simulations.
        """
        if os.path.exists(self.simulations_dir):
            shutil.rmtree(self.simulations_dir)
            logging.info("Deleting dir {}".format(self.simulations_dir))
        os.makedirs(self.simulations_dir)
        try:
            self.error_matrix = pd.read_csv(
                os.path.join(self.data_dir, self.error_filename))
        except FileNotFoundError:
            pass # Allow no error file in case we are not using error: catch issues later
        self.verbosity = args.verbosity
        logging.info("Creating dir {}".format(self.simulations_dir))
        self.data = self.run_simulations(
            args.replicates, args.seed, args.progress, args.processes)
        for t in self.tools_and_metrics.keys():
            col = t + "_completed"
            if col not in self.data:
                self.data[col] = False
        # Other result columns are added later during the infer step.
        self.dump_data(write_index=True)

    def run_simulations(self, replicates=None, seed=None, show_progress=False, num_processes=1):
        """
        Called from setup(): requires us to define the dictionaries
        self.between_sim_params and self.within_sim_params.
        if self.filter_between_sim_params is given then it is a function passed to
        filter to only leave in certain param combinations
        """
        rng = random.Random(seed or self.default_seed)
        n_reps = replicates or self.default_replicates
        if hasattr(self,'filter_between_sim_params'):
            logging.info("Filtering out some parameter values")
            param_iter, count_param_iter = itertools.tee(
                filter(lambda params: self.filter_between_sim_params(params),
                    itertools.product(
                        #always call with replicate number first, then all other params
                        range(n_reps), *self.between_sim_params.values())))
        else:
            param_iter, count_param_iter = itertools.tee(
                itertools.product(
                    #always call with replicate number first, then all other params
                    range(n_reps), *self.between_sim_params.values()))

        num_sims = sum(1 for _ in count_param_iter)
        num_rows = num_sims * np.prod([len(a) for a in self.within_sim_params.values()])

        self.progress = tqdm.tqdm(total=num_rows) if show_progress else None
        #Predefine a set of seeds so that we get a different seed for each
        #thread (NB if we filter out sims, this may be unnecessarily long)
        seeds = [rng.randint(1, 2**31) for i in range(num_sims)]
        combined_iterator = enumerate(zip(seeds, param_iter))

        def save_result(data, values_by_row):
            for i,d in values_by_row.items():
                for colname, val in d.items():
                    data.iloc[i][colname]=val
                if self.progress is not None:
                    self.progress.update()

        data = pd.DataFrame(index=np.arange(0, num_rows), columns=self.sim_cols + self.extra_sim_cols)

        if num_processes > 1:
            logging.info("Setting up using multiprocessing ({} processes)".format(num_processes))
            with multiprocessing.Pool(processes=num_processes, maxtasksperchild=2) as pool:
                for result in pool.imap_unordered(self.run_single_sim, combined_iterator):
                    save_result(data, result)
        else:
            # When we have only one process it's easier to keep everything in the same
            # process for debugging.
            logging.info("Setting up using a single process")
            for result in map(self.run_single_sim, combined_iterator):
                save_result(data, result)
        return data

    def run_single_sim(self, runtime_information):
        """
        Run the routine self.single_sim() - could be as a separate process
        """
        i, (seed_plus_params) = runtime_information
        rng_seed, params = seed_plus_params
        sim_params = dict(zip(['replicate'] + list(self.between_sim_params.keys()), params))
        #sims may have multiple rows, e.g. one for each error rate
        #so the row numbers will start at multiples of this
        initial_row_id = i * np.prod([len(x) for x in self.within_sim_params.values()])
        return self.single_sim(initial_row_id, sim_params, random.Random(rng_seed))

    def save_within_sim_data(self, base_row_id, ts, base_fn, base_params, iterate_over_dict=None):
        """
        save trees, and variant matrices with error, and return the values to be save in the csv file
        iterate_over_dict can be used to pass a subset of values to iterate over, rather than the default
        of self.within_sim_params
        """

        if iterate_over_dict is None:
            iterate_over_dict = self.within_sim_params
        if any(self.tools_and_metrics.values()):
            if base_params.get('restrict_sample_size_comparison'):
                #we are not comparing the trees with all samples, but just the first n
                n = base_params['restrict_sample_size_comparison']
                cmp_fn = add_subsample_param_to_name(base_fn, n)
                small_ts = ts.simplify(list(range(n)))
                small_ts.save_nexus_trees(cmp_fn +".nex")
                #We must also get the locations of variants out of this file, so we can compare
                #metrics fairly (i.e. higher resolution inference with more samples doesn't
                #have an advantage in getting more information to locate breakpoints
                #Note that we might accidentally create a TS with no valid sites here
                if small_ts.num_sites:
                    self.generate_samples(small_ts, cmp_fn)
            else:
                ts.save_nexus_trees(base_fn +".nex")
        return_value = {}
        for params in itertools.product(*iterate_over_dict.values()):
            #may be iterating here over error rates
            keyed_params = dict(zip(iterate_over_dict.keys(), params))
            keyed_params.update(base_params)
            row = return_value[base_row_id] = {}
            base_row_id += 1
            for k,v in keyed_params.items():
                row[k] = v
            #add some stats from the ts
            row['edges'] = ts.num_edges
            row['n_trees'] = ts.num_trees
            row['n_sites'] = ts.num_sites
            subsample = keyed_params.get(SUBSAMPLE_COLNAME)
            if subsample:
                self.save_variant_matrices(
                    ts.simplify(list(range(subsample))),
                    add_subsample_param_to_name(base_fn, subsample),
                    keyed_params.get(ERROR_COLNAME) or 0,
                    infinite_sites=False)
            else:
                self.save_variant_matrices(
                    ts, base_fn, 
                    keyed_params.get(ERROR_COLNAME) or 0,
                    infinite_sites=False)
        return return_value

    def single_sim(self, row_id, sim_params, rng):
        """
        Return a dictionary whose keys are the row numbers which
        will be inserted into the csv file, and whose values are
        dictionaries keyed by column name
        """
        raise NotImplementedError(
            "A dataset should implement a single_sim method that runs" \
            " a simulation, e.g. by calling self.single_neutral_simulation")

    def generate_samples(self, ts, filename, error_param=0):
        """
        Generate a samples file from a simulated ts based on the empirically estimated 
        error matrix saved in self.error_matrix.
        Reject any variants that result in a fixed column. 
        """
        record_rate = logging.getLogger().isEnabledFor(logging.INFO)
        n_variants = bits_flipped = 0
        assert ts.num_sites != 0
        sample_data = tsinfer.SampleData(path=filename+".samples", sequence_length=ts.sequence_length)
        for v in ts.variants():
            n_variants += 1
    
            try:
                error_param = float(error_param)
                if error_param == 0:
                    genotypes=v.genotypes
                else:
                    raise NotImplementedError
            except ValueError:
                # Error_param is not a number => is a error file
                # First record the allele frequency
                m = v.genotypes.shape[0]
                frequency = np.sum(v.genotypes) / m
                try:
                    #Find closest row in error matrix file
                    closest_row = (self.error_matrix['freq']-frequency).abs().argsort()[:1]
                    closest_freq = self.error_matrix.iloc[closest_row]
                except AttributeError:
                    logging.warning("Empirical error matrix not available at: {}".format(
                        self.error_filename))
                    raise
                genotypes = make_errors_genotype_model(v.genotypes,closest_freq)
                if record_rate:
                    bits_flipped += np.sum(np.logical_xor(genotypes, v.genotypes))

            sample_data.add_site(
                position=v.site.position, alleles=v.alleles,
                genotypes=genotypes)

        if error_param != 0:
            logging.info("Error = {} used with actual error rate = {}".format(
                error_param,
                bits_flipped/(n_variants*ts.sample_size)) if record_rate else "")
      
        sample_data.finalise()
        return sample_data
    
    def infer(
            self, num_processes, num_threads, force=False, metrics_only=False,
            specific_tool=None, specific_row=None, flush_all=False, show_progress=False):
        """
        Runs the main inference processes and stores results in the dataframe.
        can 'force' all rows to be (re)run, or specify a specific row to run.
        """
        self.load_data()
        tools = self.tools_and_metrics.keys()
        if specific_tool is not None:
            if specific_tool not in self.tools_and_metrics:
                raise ValueError("Tool '{}' not recognised: options = {}".format(
                    specific_tool, list(self.tools_and_metrics.keys())))
            tools = [specific_tool]
        row_ids = self.data.index
        if specific_row is not None:
            if specific_row < 0 or specific_row > len(self.data.index):
                raise ValueError("Row {} out of bounds".format(specific_row))
            row_ids = [specific_row]
        work = []
        tool_work_total = {tool: 0 for tool in tools}
        for row_id in row_ids:
            row = self.data.iloc[row_id]
            for tool in tools:
                #special case for ARGweaver, to allow multiple runs for the same sim
                if getattr(row,"only_AW",False) and tool != ARGWEAVER:
                    logging.info("Data row {} is set to skip {} inference".format(
                        row_id,tool) + " (probably to avoid duplicate effort)")
                    continue
                # Only run for a particular tool if it has not already completed
                # of if --force is specified.
                if metrics_only and not row[tool + "_completed"]:
                    logging.info(
                        "Data row {} has not completed {} inference: skipping metric calculations".format(
                            row_id, tool))
                elif not metrics_only and not force and row[tool + "_completed"]:
                    logging.info(
                        "Data row {} is filled out for {} inference: skipping".format(
                            row_id, tool))
                else:
                    work.append((
                        tool, row, self.simulations_dir, num_threads,
                        self.tools_and_metrics[tool], metrics_only,
                        self.random_resolve_polytomy_replicates))
                    tool_work_total[tool] += 1
        logging.info(
            "running {} {} (max {} tools over {} of {} rows) with {} "
            "processes and {} threads".format(
                len(work), "metric calculations" if metrics_only else "inference trials",
                len(tools), int(np.ceil(len(work)/len(self.tools_and_metrics))),
                len(self.data.index), num_processes, num_threads))

        # Randomise the order that work is done in so that we get results for all parts
        # of the plots through rather than get complete results for the initial data
        # points.
        random.shuffle(work)
        tool_work_completed = {tool: 0 for tool in tools}
        if show_progress:
            width = max(len(tool) for tool in tools)
            progress = {
                tool:tqdm.tqdm(
                    desc="{:>{}}".format(tool, width),
                    total=tool_work_total[tool]) for tool in tools}

        def store_result(row_id, tool, results):
            tool_work_completed[tool] += 1
            logging.info("{} {}/{} completed for {}".format(
                "Metric calculation" if metrics_only else "Inference",
                tool_work_completed[tool], tool_work_total[tool], tool))
            for k, v in results.items():
                if any([re.search(regexp, k) for regexp in self.exclude_colnames_matching]):
                    pass
                else:
                    self.data.ix[row_id, k] = v
            self.dump_data(force_flush=flush_all)
            # Update the progress meters
            if show_progress:
                progress[tool].update()

        if num_processes > 1:
            with multiprocessing.Pool(processes=num_processes) as pool:
                for result in pool.imap_unordered(infer_worker, work):
                    store_result(*result)
        else:
            # When we have only one process it's easier to keep everything in the same
            # process for debugging.
            for result in map(infer_worker, work):
                store_result(*result)
        self.dump_data(force_flush=True)

        if show_progress:
            for p in progress.values():
                p.close()
    #
    # Utilities for running simulations and saving files.
    #
    def single_neutral_simulation(self, sample_size, Ne, length, recombination_rate, mutation_rate, seed,
        mut_seed=None, replicate=None, **kwargs):
        """
        The standard way to run one msprime simulation for a set of parameter
        values. Saves the output to an .hdf5 file, and also saves variant files
        for use in fastARG (a .hap file, in the format specified by
        https://github.com/lh3/fastARG#input-format) ARGweaver (a .sites file:
        http://mdrasmus.github.io/argweaver/doc/#sec-file-sites) tsinfer
        (currently a numpy array containing the variant matrix)

        mutation_seed allows the same
        ancestry to be simulated (if the same genealogy_seed is given) but have
        different mutations thrown onto the trees (even with different
        mutation_rates). If None, then different vaues of mu will create different
        simulations.

        If Ne is a string, it is an identified used for the population model
        which will be saved into the filename & results file, while the simulation
        will be called with Ne=1

        replicate is useful to pass in for error messages etc.

        Returns a tuple of treesequence, filename (without file type extension)
        """
        logging.info(
            "Running neutral simulation for "
            "n = {}, l = {:.2g}, Ne = {}, rho = {}, mu = {}".format(
                sample_size, length, Ne, recombination_rate, mutation_rate))
        # Since we want to have a finite site model, we force the recombination map
        # to have exactly l loci with a recombination rate of rho between them.
        recombination_map = msprime.RecombinationMap.uniform_map(length, recombination_rate, length)
        rng1 = random.Random(seed)
        sim_seed = rng1.randint(1, 2**31)
        if "population_configurations" in kwargs:
            Ne_used = 1
            n_used = None
        else:
            Ne_used = Ne
            n_used = sample_size
        if mut_seed is None:
            ts = msprime.simulate(
                n_used, Ne_used, recombination_map=recombination_map, mutation_rate=mutation_rate,
                random_seed=sim_seed, **kwargs)
        else:
            #run with no mutations (should give same result regardless of mutation rate)
            ts = msprime.simulate(
                n_used, Ne_used, recombination_map=recombination_map, mutation_rate=0,
                random_seed=sim_seed, **kwargs)
            #now add the mutations
            rng2 = msprime.RandomGenerator(mut_seed)
            muts = msprime.MutationTable()
            tables = ts.dump_tables()
            mutgen = msprime.MutationGenerator(rng2, mutation_rate)
            mutgen.generate(tables.nodes, tables.edges, tables.sites, muts)
            msprime.sort_tables(
                nodes=tables.nodes, edges=tables.edges, sites=tables.sites, mutations=muts)
            ts = msprime.load_tables(nodes=nodes, edges=edges, sites=sites, mutations=muts)

        logging.info(
            "Neutral simulation done; {} sites, {} trees".format(ts.num_sites, ts.num_trees))
        sim_fn = mk_sim_name(sample_size, Ne, length, recombination_rate, mutation_rate, seed, mut_seed, self.simulations_dir)
        logging.debug("writing {}.hdf5".format(sim_fn))
        ts.dump(sim_fn+".hdf5")

        # Make sure that there is *some* information in this simulation that can be used
        # to infer a ts, otherwise it's pointless
        if ts.get_num_mutations() == 0:
            raise ValueError("No mutations present")
        if ts_has_non_singleton_variants(ts) == False:
            raise ValueError("No non-singleton variants present ({} singletons)".format(
                sum([np.sum(v.genotypes)==1 for v in ts.variants()])))

        return ts, sim_fn

    def single_simulation_with_human_demography(self, sample_size, sim_name, length,
        recombination_rate, mutation_rate, seed, mut_seed=None, **kwargs):
        """
        Run an msprime simulation with a rough approximation of recent human demography
        using the  Gutenkunst et al., (2009) model used by a number of other simulators
        (e.g. msms), which is an elaboration of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1310645/

        (see http://msprime.readthedocs.io/en/stable/tutorial.html#demography)

        the sample_size, n, is divided into 3 roughly equal-sized groups from africa, europe, & china
        (YRI, CEU, and CHB)
        """
        def out_of_africa(nYRI, nCEU, nCHB):
            """
            Straight from the msprime docs, but return the setup params for passing to the
            simulate() function
            """
            # First we set out the maximum likelihood values of the various parameters
            # given in Table 1.
            N_A = 7300
            N_B = 2100
            N_AF = 12300
            N_EU0 = 1000
            N_AS0 = 510
            # Times are provided in years, so we convert into generations.
            generation_time = 25
            T_AF = 220e3 / generation_time
            T_B = 140e3 / generation_time
            T_EU_AS = 21.2e3 / generation_time
            # We need to work out the starting (diploid) population sizes based on
            # the growth rates provided for these two populations
            r_EU = 0.004
            r_AS = 0.0055
            N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
            N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
            # Migration rates during the various epochs.
            m_AF_B = 25e-5
            m_AF_EU = 3e-5
            m_AF_AS = 1.9e-5
            m_EU_AS = 9.6e-5
            # Population IDs correspond to their indexes in the population
            # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
            # initially.
            population_configurations = [
                msprime.PopulationConfiguration(
                    sample_size= nYRI, initial_size=N_AF),
                msprime.PopulationConfiguration(
                    sample_size= nCEU, initial_size=N_EU, growth_rate=r_EU),
                msprime.PopulationConfiguration(
                    sample_size= nCHB, initial_size=N_AS, growth_rate=r_AS)
            ]
            migration_matrix = [
                [      0, m_AF_EU, m_AF_AS],
                [m_AF_EU,       0, m_EU_AS],
                [m_AF_AS, m_EU_AS,       0],
            ]
            demographic_events = [
                # CEU and CHB merge into B with rate changes at T_EU_AS
                msprime.MassMigration(
                    time=T_EU_AS, source=2, destination=1, proportion=1.0),
                msprime.MigrationRateChange(time=T_EU_AS, rate=0),
                msprime.MigrationRateChange(
                    time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
                msprime.MigrationRateChange(
                    time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
                msprime.PopulationParametersChange(
                    time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
                # Population B merges into YRI at T_B
                msprime.MassMigration(
                    time=T_B, source=1, destination=0, proportion=1.0),
                # Size changes to N_A at T_AF
                msprime.PopulationParametersChange(
                    time=T_AF, initial_size=N_A, population_id=0)
            ]

            return dict(
                population_configurations=population_configurations,
                migration_matrix=migration_matrix,
                demographic_events=demographic_events
                )

        asian_samples = sample_size//3
        afro_european_samples = sample_size - asian_samples
        european_samples = afro_european_samples//2
        african_samples = afro_european_samples - european_samples
        assert african_samples + european_samples + asian_samples == sample_size

        kwargs.update(out_of_africa(african_samples, european_samples, asian_samples))
        return self.single_neutral_simulation(sample_size, sim_name, length,
            recombination_rate, mutation_rate, seed, mut_seed, **kwargs)

    def single_simulation_with_selective_sweep(self, sample_size, Ne, length, recombination_rate, mutation_rate,
        selection_coefficient, dominance_coefficient, stop_freqs, seed, mut_seed=None, replicate=None):
        """
        Run a forward simulation with a selective sweep for a set of parameter values.
        using simuPOP.

        Forward simulations are slow, especially for large chromosomes & population sizes.
        So to save time, we can get multiple results from the same simulation, by saving
        populations at different times


        if stop is a single number it is taken as a stop frequency, and
        the simulation is stopped when the selected variant reaches that frequency. If it is a
        tuple, it is taken as (freq, generations), and the simulation is stopped at the
        specified number of generations after that frequency has been reached (e.g. if
        stop = (1.0, 200), the simulation is stopped 200 generations after fixation.

        Convert the output to multiple .hdf5 files using ftprime. Other
        details as for single_neutral_simulation()
        Returns an iterator over (treesequence, filename, output_freq) tuples
        (without file type extension)

        """
        #Don't `import` simulate_sweep at top of file as it fires up a simupop instance
        # which is not always needed
        from selective_sweep import simulate_sweep
        logging.info(
            "Running forward simulation with selection for " + \
            "n = {}, Ne = {}, l = {:.2g}, rho = {}, mu = {}, s = {}".format(
                sample_size, Ne, length, recombination_rate, mutation_rate, selection_coefficient))
        sim_fn = mk_sim_name(sample_size, Ne, length, recombination_rate, mutation_rate, seed, mut_seed, self.simulations_dir,
            tool="ftprime", s=selection_coefficient, h=dominance_coefficient) + "_f" #freq + post_gens added by simulate_sweep()


        saved_files = simulate_sweep(
            popsize = Ne,
            chrom_length = length,
            recomb_rate = recombination_rate,
            mut_rate = mutation_rate,
            selection_coef = selection_coefficient,
            dominance_coef = dominance_coefficient,
            output_at_freqs = [(s,0) if isinstance(s, (str, float, int)) else (s[0],s[1]) for s in stop_freqs],
            nsamples = int(sample_size/2), #sample n/2 diploid individuals from the population
            max_generations=int(1e8), #bail after ridiculous numbers of generations
            mutations_after_simulation=True,
            treefile_prefix=sim_fn,
            seed=seed)

        expected_suffix = ".hdf5"
        for outfreq, fn in saved_files.items():
            assert fn.endswith(expected_suffix)
            ts = msprime.load(fn)
            logging.info(
                "Selective simulation from '{}'; {} sites ({} mutations), {} trees".format(
                    fn, ts.num_sites,ts.get_num_mutations(), ts.num_trees))

            # Simulations with a selective sweep might well have reduced diversity
            # so allow this, but log a warning (not sure if this will break e.g. ARGweaver)
            if ts.get_num_mutations() == 0:
                logging.warning("No mutations present")
            if ts_has_non_singleton_variants(ts) == False:
                logging.warning("No non-singleton variants present ({} singletons) for output at {}".format(
                    sum([np.sum(v.genotypes)==1 for v in ts.variants()]), outfreq))
            yield ts, fn[:-len(expected_suffix)], outfreq



    def save_variant_matrices(self, ts, filename, error_param=0, infinite_sites=True):
        """
        Make sample data from a tree sequence
        """
        if infinite_sites:
            #for infinite sites, assume we have discretised mutations to ints
            if not all(p.is_integer() for p in pos):
                raise ValueError("Variant positions are not all integers")
        filename = add_error_param_to_name(filename, error_param)
        if ts.num_sites == 0:
            logging.warning("No sites to save for {}".format(filename))
        else:
            logging.debug("Saving samples to {}".format(filename))
            s = self.generate_samples(ts, filename, error_param)
            if FASTARG in self.tools_and_metrics:
                logging.debug("writing samples to {}.hap for fastARG".format(filename))
                with open(filename+".hap", "w+") as file_in:
                    ts_fastARG.samples_to_fastARG_in(s, file_in)
            if ARGWEAVER in self.tools_and_metrics:
                logging.debug("writing samples to {}.sites for ARGweaver".format(filename))
                with open(filename+".sites", "w+") as file_in:
                    ts_ARGweaver.samples_to_ARGweaver_in(s, file_in, infinite_sites=infinite_sites)
            if RENTPLUS in self.tools_and_metrics:
                logging.debug("writing samples to {}.dat for RentPlus".format(filename))
                with open(filename+".dat", "wb+") as file_in:
                    ts_RentPlus.samples_to_RentPlus_in(s, file_in, infinite_sites=infinite_sites)



class AllToolsDataset(Dataset):
    """
    ARG inference using all tools (requires small sample
    sizes to allow ARGweaver to run). The x-axis is centred around regions 
    of biological interest (the mutation rate is symmetrical in log space
    either side of mu/rho = 1). 
    
    For symmetry around e.g. mu = 1e-8, we can pick a lower bound x and set
    the upper bound as (1e-8 ** 2)/x
    """
    name = "all_tools"
    default_replicates = 100
    default_seed = 123

    def geomspace_around(middle, bound, **kwargs):
        """
        return a list of numbers evenly spaced in log space around a central value
        """
        if bound < middle:
            return np.geomspace(bound, middle**2/bound, **kwargs)
        else:
            return np.geomspace(middle**2/bound, bound, **kwargs)
    
    #params that change BETWEEN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    between_sim_params = {
        'Ne': [5000],
        'mutation_rate': geomspace_around(1e-8, 5e-7, num=7), #gives upper bound of 5e-7 and lower of 2e-10
        'sample_size':   [16],
        'length':        [1000000], #1Mb ensures that low mutation rates ~2e-10 still have some variants
        'recombination_rate': [1e-8],
    }
    #params that change WITHIN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    within_sim_params = {
        ERROR_COLNAME : [0, Dataset.error_filename.replace(".csv","")],
    }
    
    def single_sim(self, row_id, sim_params, rng):
        while True:
            sim_params['seed'] = rng.randint(1, 2**31)
            try:
                # Run the simulation until we get an acceptable one
                ts, fn = self.single_neutral_simulation(**sim_params)
                break
            except ValueError as e: #No non-singleton variants
                logging.warning(e)
        row_data = self.save_within_sim_data(row_id, ts, fn, sim_params)
        return row_data


class AllToolsAccuracyDataset(AllToolsDataset):
    """
    ARG inference using all tools, with mutation rate increasing 
    to very high levels (relative to recombination). We expect accuracy
    (measured by tree metrics) to increase (metrics to decrease to 0)
    as mutation rate increases (and more data exists for inference).
    """
    name = "all_tools_accuracy"


    #params that change BETWEEN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    between_sim_params = {
        'Ne': [5000],
        'mutation_rate': np.geomspace(0.5e-8, 3.5e-6, num=7),
        'sample_size':   [16],
        'length':        [100000],  #should be enough for ~ 50 trees
        'recombination_rate': [1e-8],
    }

    #params that change WITHIN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    within_sim_params = {
        ERROR_COLNAME : [0, AllToolsDataset.error_filename.replace(".csv","")],
    }

class AllToolsPerformanceDataset(AllToolsDataset):
    """
    ARG inference using all tools, with various sample sizes so that we 
    can inspect the scaling properties of the algorithms. 
    """
    name = "all_tools_performance"
    default_replicates = 10

    fixed_sample_size = 10
    fixed_length = 2 * 10**5
    num_points = 4
    #Ensure sample sizes are even so we can output diploid VCF.
    sample_sizes = [10, 20, 40, 65, 100]
    lengths = np.linspace(fixed_length / 5, 4 * fixed_length, num_points).astype(int)

    #params that change BETWEEN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    between_sim_params = {
        'Ne':            [5000],
        'mutation_rate': [1e-8],
        #Ensure sample sizes are even so we can output diploid VCF.
        'sample_size':   np.unique(np.append(sample_sizes, fixed_sample_size)),
        'length':        np.unique(np.append(lengths, fixed_length)),
        'recombination_rate': [1e-8],
    }

    #params that change WITHIN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    within_sim_params = {
        ERROR_COLNAME : [0],
    }

    def filter_between_sim_params(self, params):
        """
        We want to only use cases where the sample_size OR the length are at the fixed values
        """
        keyed_params = dict(zip(['replicates']+list(self.between_sim_params.keys()), params))
        return (keyed_params['sample_size']==self.fixed_sample_size) or \
            (keyed_params['length']==self.fixed_length)



class TsinferPerformanceDataset(AllToolsPerformanceDataset):
    """
    The performance of tsinfer in terms of CPU time, memory usage and
    compression rates for large scale datasets. Contains data for
    the two main dimensions of sample size and sequence length.
    """
    name = "tsinfer_performance"
    tools_and_metrics = {TSINFER:[]} #run tsinfer, but switch metrics off (too slow)

    default_replicates = 10
    fixed_sample_size = 50000
    fixed_length = 50 * 10**6
    num_points = 20
    #Ensure sample sizes are even so we can output diploid VCF.
    sample_sizes = np.linspace(0, fixed_sample_size, num_points + 1)[1:].astype(int) * 2
    lengths = np.linspace(fixed_length / 10, 2 * fixed_length, num_points).astype(int)

    #params that change BETWEEN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    between_sim_params = {
        'Ne':            [5000],
        'mutation_rate': [1e-8],
        #Ensure sample sizes are even so we can output diploid VCF.
        'sample_size':   np.unique(np.append(sample_sizes, fixed_sample_size)),
        'length':        np.unique(np.append(lengths, fixed_length)),
        'recombination_rate': [1e-8/2, 1e-8, 1e-8*2],
    }

    #params that change WITHIN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    within_sim_params = {
        ERROR_COLNAME : [0],
    }

    extra_sim_cols = ["ts_filesize", "vcf_filesize", "vcfgz_filesize"]


    def single_sim(self, row_id, sim_params, rng):
        """
        Modify to include file sizes, which is useful for large dataset comparison
        Note that the vcf files may take up *substantial* amounts of space
        """
        while True:
            sim_params['seed'] = rng.randint(1, 2**31)
            try:
                # Run the simulation until we get an acceptable one
                ts, fn = self.single_neutral_simulation(**sim_params)
                break
            except ValueError as e: #No non-singleton variants
                logging.warning(e)
        vcf_filename = fn + ".vcf"
        with open(vcf_filename, "w") as vcf_file:
            ts.write_vcf(vcf_file, 2)
        vcfgz_filename = vcf_filename + ".gz"
        vcf_filesize = os.path.getsize(vcf_filename)
        # gzip removes the original by default, so might as well save some space.
        subprocess.check_call(["gzip", vcf_filename])
        vcfgz_filesize = os.path.getsize(vcfgz_filename)
        # only need the file size, so can delete the .gz file too
        os.remove(vcfgz_filename)
        row_data = self.save_within_sim_data(row_id, ts, fn, dict(sim_params,
            ts_filesize = os.path.getsize(fn + ".hdf5"),
            vcf_filesize = vcf_filesize,
            vcfgz_filesize = vcfgz_filesize
            ))
        return row_data


class FastargTsinferComparisonDataset(AllToolsPerformanceDataset):
    """
    The performance of the various programs in terms of running time and memory usage
    """
    name = "fastarg_tsinfer_comparison"
    default_replicates = 20
    default_seed = 1000
    tools_and_metrics = {FASTARG:[], TSINFER:[]} # Everything else is too slow
    fixed_sample_size = 10000
    fixed_length = 5 * 10**6
    num_points = 20
    sample_sizes = np.linspace(10, 2*fixed_sample_size, num_points).astype(int)
    lengths = np.linspace(fixed_length / 10, 2 * fixed_length, num_points).astype(int)
    #params that change BETWEEN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    between_sim_params = {
        'Ne':            [5000],
        'mutation_rate': [1e-8],
        'recombination_rate': [1e-8],
        'sample_size':   np.unique(np.append(sample_sizes, fixed_sample_size)),
        'length':        np.unique(np.append(lengths, fixed_length)),
    }
    #params that change WITHIN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    within_sim_params = {
        ERROR_COLNAME : [0],
    }

### SUPPLEMENTARY MATERIAL
### The following classes are used for figures in the supplementary material
###

class SubsamplingDataset(Dataset):
    """
    Accuracy of ARG inference of a fixed subset of samples (measured by various statistics)
    as the population sample size increases.
    """
    name = "subsampling"
    #to make this a fair comparison, we need to calculate only at the specific variant sites
    #because different sample sizes will be based on different original variant positions
    #which are then cut down to the subsampled ones.
    tools_and_metrics = {
        TSINFER: [METRICS_LOCATION_VARIANTS| METRICS_POLYTOMIES_LEAVE, METRICS_LOCATION_VARIANTS | METRICS_POLYTOMIES_BREAK]
    }
    default_replicates = 100
    default_seed = 123

    #params that change BETWEEN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    between_sim_params = {
        'sample_size':        [1000], #the maximum sample size - we will trim this down
        'length':             [10000, 100000, 1000000],
        'Ne':                 [5000],
        'mutation_rate':      [1e-8],
        'recombination_rate': [1e-8],
         #we always simplify down to the first n samples for comparison purposes
        'restrict_sample_size_comparison': [6],
    }

    within_sim_params = {
        'tsinfer_srb' : [True], #, False], #should we use shared recombinations ("path compression")
        SUBSAMPLE_COLNAME:  [12, 50, 100, 500, 1000], #we infer based on this many samples
        ERROR_COLNAME: [0, Dataset.error_filename.replace(".csv","")]
    }


    extra_sim_cols = [SUBSAMPLE_COLNAME, 'restrict_sample_size_comparison']

    def single_sim(self, row_id, sim_params, rng):
        assert max(self.within_sim_params[SUBSAMPLE_COLNAME]) <= sim_params['sample_size']
        base_row_id = row_id
        while True:
            sim_params['seed'] = rng.randint(1, 2**31)
            try:
                # Run the simulation until we get an acceptable one
                base_ts, fn = self.single_neutral_simulation(
                    **{k:v for k,v in sim_params.items() if k != 'restrict_sample_size_comparison'})
                break
            except ValueError as e: #No non-singleton variants
                logging.warning(e)
        return_data = {}
        #it's worth saving the full tree file for other analysis
        base_ts.save_nexus_trees(fn + ".nex")
        #this will loop over the subsamples and errors, saving samples files for inference
        return_data = self.save_within_sim_data(base_row_id, base_ts, fn, sim_params)
        return return_data

class AllToolsAccuracyWithDemographyDataset(Dataset):
    """
    Accuracy of ARG inference (measured by various statistics)
    tending to fully accurate as mutation rate increases, but with a
    demographic (out of africa) rather than a neutral model.
    
    This is faster than the non-demography simulation, as the 
    population expansions etc mean that the Ne is slightly lower
    """

    name = "all_tools_accuracy_with_demography"

    default_replicates = 100
    default_seed = 123

    #params that change BETWEEN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    between_sim_params = {
        'mutation_rate': np.geomspace(0.5e-8, 3.5e-6, num=7),
        'sample_size':   [15], #will be split across the 3 human sub pops
        'length':        [100000],
        'recombination_rate': [1e-8],
        'sim_name': ["Gutenkunst.out.of.africa"]
    }

    #params that change WITHIN simulations. Keys should correspond
    # to column names in the csv file. Values should all be arrays.
    within_sim_params = {
        ERROR_COLNAME : [0, Dataset.error_filename.replace(".csv","")],
    }

    def single_sim(self, row_id, sim_params, rng):

        while True:
            sim_params['seed'] = rng.randint(1, 2**31)
            try:
                # Run the simulation until we get an acceptable one
                ts, fn = self.single_simulation_with_human_demography(**sim_params)
                break
            except ValueError as e: #No non-singleton variants
                logging.warning(e)
        row_data = self.save_within_sim_data(row_id, ts, fn, sim_params)
        #set Ne to store the name 'Gutenkunst.out.of.africa' rather than a meaningless Ne param
        #this will be used in the filenames etc.
        for row in row_data.values():
            row['Ne'] = row['sim_name']
        return row_data

class AllToolsAccuracyWithSelectiveSweepDataset(Dataset):
    """
    Accuracy of ARG inference (measured by various statistics)
    tending to fully accurate as mutation rate increases, but with a
    selective sweep superimposed.

    The time since / post sweep is likely to be a major factor in
    misfitting the model. We judge this by stopping the sweep at a
    specific frequency, or at a given number of generations post fixation.
    """

    name = "all_tools_accuracy_with_selective_sweep"

    default_replicates = 50
    default_seed = 123

    between_sim_params = {
        'mutation_rate': np.geomspace(0.5e-8, 3.5e-6, num=5),
        'sample_size': [15],
        'Ne':     [5000],
        'length': [100000],
        'recombination_rate': [1e-8],
        SELECTION_COEFF_COLNAME: [0.1],
        DOMINANCE_COEFF_COLNAME: [0.5],
    }

    within_sim_params = {
        #NB - these are strings because they are output as part of the filename
        'stop_freqs': ['0.2', '0.5', '0.8', '1.0', ('1.0', 200), ('1.0', 1000)], #frequencies when file is saved.
        ERROR_COLNAME : [0, Dataset.error_filename.replace(".csv","")],
    }

    extra_sim_cols = [SIMTOOL_COLNAME, SELECTION_COEFF_COLNAME,
            DOMINANCE_COEFF_COLNAME, SELECTED_FREQ_COLNAME, SELECTED_POSTGEN_COLNAME]

    def single_sim(self, row_id, sim_params, rng):
        while True:
            base_row_id = row_id
            return_data = {}
            sim_params['seed'] = rng.randint(1, 2**31)
            try:
                #we have multiple rows per simulation for results after different generations post-fixation
                #these are returned in an iterator by the single_simulation_with_selective_sweep() method
                #so we pop this set of permutations off the auto-iterated list
                within_sim_params = self.within_sim_params.copy()
                stop_freqs = within_sim_params.pop('stop_freqs')
                # Pass stopping freqs into the sim and loop through the returned freqs
                # Reject this entire set if we got no mutations at any point, or all mutations are fixed
                for ts, fn, output_info in self.single_simulation_with_selective_sweep(stop_freqs = stop_freqs, **sim_params):
                    #create new parameters to save to the csv file
                    params = {SIMTOOL_COLNAME:    'ftprime',
                        SELECTED_FREQ_COLNAME:    output_info[0],
                        SELECTED_POSTGEN_COLNAME: output_info[1] if len(output_info)>1 else 0}
                    params.update(sim_params)
                    row_data = self.save_within_sim_data(base_row_id, ts, fn, params, iterate_over_dict = within_sim_params)
                    base_row_id += len(row_data)
                    assert all(k not in return_data for k in row_data)
                    return_data.update(row_data)
                break #not rejected, can finish
            except ValueError as e: #No non-singleton variants
                logging.warning(e)
            except (_msprime.LibraryError, MemoryError) as e:
                #we sometimes run out of memory here
                logging.warning("Error when running `single_sim()`, probably in single_simulation_with_selective_sweep({})"\
                    .format(",".join(["{}={}".format(k,v) for k,v in sim_params.items()])))
                raise
        return return_data

### ADDITIONAL DEBUG
### The following classes are used for drilling down into particular features
### of the data, e.g. why we do better than AW for high mutation rates,
### some debugging routines, etc.
###


class ARGweaverParamChangesDataset(Dataset):
    """
    Accuracy of ARGweaver inference (measured by various tree statistics)
    as we adjust burn-in time and number of time slices.
    This is an attempt to explain why ARGweaver can do badly e.g. for high mutation rates

    You can check that the timeslices really *are* having an effect by looking at the unique
    times (field 3) in the .arg files within raw__NOBACKUP__/argweaver_param_changes e.g.
    cut -f 3 data/raw__NOBACKUP__/argweaver_param_changes/simulations/<filename>.arg | sort | uniq
    """
    name = "argweaver_param_changes"
    tools = {
        ARGWEAVER:[METRICS_LOCATION_ALL | METRICS_POLYTOMIES_LEAVE], 
        TSINFER:  [METRICS_LOCATION_ALL | METRICS_POLYTOMIES_LEAVE, METRICS_LOCATION_ALL | METRICS_POLYTOMIES_BREAK]}
    default_replicates = 40
    default_seed = 123

    between_sim_params = {
        'mutation_rate': np.geomspace(1e-6, 1e-3, num=6),
        'sample_size': [6],
        'Ne':     [5000],
        'length': [10**4],
        'recombination_rate': [1e-8],
    }

    within_sim_params = {
        'AW_burnin_steps': [1000,2000,5000], #by default, bin/arg-sim has no burnin time
        'AW_num_timepoints': [20,60,200], #by default, bin/arg-sim has n=20
        ERROR_COLNAME : [0],
    }

    extra_sim_cols = ["only_AW", "ARGweaver_burnin", "ARGweaver_ntimes"]

    def single_sim(self, row_id, sim_params, rng):
        return_data = {}
        while True:
            sim_params['seed'] = rng.randint(1, 2**31)
            try:
                # Run the simulation
                ts, fn = self.single_neutral_simulation(**sim_params)
                break;
            except ValueError as e: #No non-singleton variants
                logging.warning(e)
        within_sim_params = self.within_sim_params.copy()
        AW_burnin_steps = within_sim_params.pop('AW_burnin_steps')
        AW_num_timepoints = within_sim_params.pop('AW_num_timepoints')
        for i, burnin in enumerate(AW_burnin_steps):
            for j, n_timesteps in enumerate(AW_num_timepoints):
                #only run other infers for the first row in this set of AW run parameters
                sim_params['only_AW'] = (0 if i==0 and j==0 else 1)
                sim_params['ARGweaver_burnin'] = burnin
                sim_params['ARGweaver_ntimes'] = n_timesteps
                row_data = self.save_within_sim_data(row_id, ts, fn, sim_params, iterate_over_dict = within_sim_params)
                row_id += len(row_data)
                assert all(k not in return_data for k in row_data)
                return_data.update(row_data)
        return return_data


######################################
#
# Data summaries
#  Here we take the CSV data files 
#  produced by setup() and infer() and
#  condense them into small CSV files
#  for plotting
#
######################################


######################################
#
# Figures
#
######################################

        

class Figure(object):
    """
    Superclass of all figures. Each figure depends on a dataset.
    """
    datasetClass = None
    name = None
    figures_dir = os.path.join(sys.path[0],"..","figures")

    tools_format = collections.OrderedDict([
        (ARGWEAVER, {"mark":"o", "col":"green"}),
        (RENTPLUS,  {"mark":"^", "col":"magenta"}),
        (FASTARG,   {"mark":"s", "col":"red"}),
        (TSINFER,   {"mark":"*", "col":"blue"}),
    ])
    metrics_params_format = collections.OrderedDict([
        (METRICS_LOCATION_ALL      | METRICS_POLYTOMIES_LEAVE, {"linestyle":"-"}),
        (METRICS_LOCATION_ALL      | METRICS_POLYTOMIES_BREAK, {"linestyle":"--"}),
        (METRICS_LOCATION_VARIANTS | METRICS_POLYTOMIES_LEAVE, {"linestyle":"-."}),
        (METRICS_LOCATION_VARIANTS | METRICS_POLYTOMIES_BREAK, {"linestyle":":"}),
    ])
    
    metric_titles = {
        "RFrooted": "Robinson-Foulds metric",
        "RFunrooted": "Robinson-Foulds metric (unrooted)",
        "SPRunrooted": "estimated SPR difference (unrooted)",
        "pathunrooted": "Path difference (unrooted)",
        "KCrooted": "Kendall-Colijn metric",
    }
    """
    Each figure has a unique name. This is used as the identifier and the
    file name for the output plots.
    """

    def __init__(self):
        self.dataset = self.datasetClass()
        self.dataset.load_data()
        self.data_file = os.path.abspath(os.path.join(self.dataset.data_dir, self.name)) + ".csv"
        
        self.tools = collections.OrderedDict(
            [(tool,a) for tool,a in self.tools_format.items() if tool in self.dataset.tools_and_metrics]
        )
        
        self.tools_and_metrics_params = collections.OrderedDict(
            [(tool + "_" + str(metric_param), dict(a, **b)) 
                for tool,a in self.tools_format.items()
                    for metric_param,b in self.metrics_params_format.items()
                        if tool in self.dataset.tools_and_metrics and metric_param in self.dataset.tools_and_metrics[tool]
            ]
        )

    def error_label(self, error, label_for_no_error = "No error"):
        """
        Make a nice label explaining this error parameter
        """
        try:
            error = float(error)
            return "Error rate = {}".format(error) if error else label_for_no_error
        except (ValueError, TypeError):
            return "{} error".format(error) if error else label_for_no_error

    def df_wide_to_long(self, input_df, prefixes_to_split, split_names=['tool']):
        """
        Convert an input dataframe to long format, on the basis of columns beginning
        with a set of prefixes (e.g. tool names)
        
        The column names can be of the form
        
            "tsinfer_per site_breaking polytomies_rooted_RF_treedist"
        
        and values separated by underscores become separate stacking variables, here

        "tsinfer", "per site", "breaking polytomies", "rooted", "RF", "treedist"

        the names to use for each variable is given by the split_names array
        """
        # set the non-tool-specific columns to indexes
        df = input_df.set_index([c for c in input_df.columns
            if not c.split("_",1)[0] in prefixes_to_split])
        row_index_depth = len(df.index.names)

        # stack the remaining columns by tool name (and possibly extra params, separated
        # by underscores, after the tool name, e.g. tsinfer_3_RF_rooted
        df.columns = pd.MultiIndex.from_tuples(
            [col.split('_', len(split_names)) for col in df.columns])
        #first part of name ("ARGweaver_" etc) is the one to stack
        df = df.stack(level=list(range(len(split_names))))

        # nix any row index levels with the same name as a stacked column (e.g. 'completed')
        for colname in df.columns:
            if colname:
                try:
                    df.index = df.index.droplevel([colname])
                    row_index_depth -= 1
                except KeyError:
                    pass
        
        # set the name of the new tool column 
        row_index_names = list(df.index.names)
        row_index_names[row_index_depth:(row_index_depth+len(split_names))] = split_names
        df.index.names = row_index_names

        # reinstate the unchanged columns (duplicates the params etc)
        return df.reset_index()
    
    def convert_treemetric_colname(self, colname, allowed_prefixes):
        """
        for a column name
        """
        colsplit = colname.split("_")
        if colsplit[0] in allowed_prefixes \
            and len(colsplit) > 1 \
            and colsplit[-1].endswith("rooted"):
            # This is a column giving an average tree metric
            if colsplit[1].isdigit():
                # This is a bit flag for metrics params:
                bitflag = int(colsplit[1])
                # (first bit is location, second is polytomies)
                bitnames = [
                    "per variant" if (METRICS_LOCATION_VARIANTS & bitflag) else "per site",
                    "breaking polytomies" if (METRICS_POLYTOMIES_BREAK & bitflag) else ""
                ]
                colsplit = colsplit[0:1] + bitnames + colsplit[2:]
            if colsplit[-1].endswith("unrooted"):
                # remove the "unrooted" suffix, and stick it as an earlier param
                colsplit.append(colsplit[-1][:-len("unrooted")])
                colsplit[-2] = "unrooted"
            elif colname.endswith("rooted"):
                colsplit.append(colsplit[-1][:-len("rooted")])
                colsplit[-2] = "rooted"
            colsplit.append("treedist")
        return "_".join(colsplit)
            
    def mean_se(self, df, groupby_columns):
        g = df.groupby(groupby_columns)
        df =pd.concat([g.mean().add_suffix("_mean"), g.sem().add_suffix("_se")], axis=1)
        return df.reset_index()
        
    def savefig(self, *figures):
        filename = os.path.join(self.figures_dir, "{}.pdf".format(self.name))
        with matplotlib.backends.backend_pdf.PdfPages(filename) as pdf:
            for fig in figures:
                pdf.savefig(fig)

    def plot(self):
        raise NotImplementedError()


    def summarize(self):
        """
        Save a summarized csv file of the results, appropriate for a particular plot.
        The default is simply to save a long (not wide) version of all the columns.
        """
        toolnames = self.dataset.tools_and_metrics.keys()
        summary_dataset = self.df_wide_to_long(self.dataset.data, toolnames)
        summary_dataset.to_csv(self.data_file)

class MetricsAllToolsFigure(Figure):
    """
    Simple figure that shows all the metrics at the same time.
    Assumes at most 2 sample sizes
    """
    datasetClass = AllToolsDataset
    name = "metrics_all_tools"
    
    
    fillstyles = ['full', 'none']

    def summarize(self):
        param_cols = ['Ne', 'length', 'sample_size',
            'recombination_rate', 'mutation_rate', ERROR_COLNAME]
        toolnames = self.dataset.tools_and_metrics.keys()

        summary_df = self.dataset.data
        
        # Turn the column names containing tree metric distances to names containing
        # rational descriptions of the parameters associated with each metric
        summary_df.columns = [
            self.convert_treemetric_colname(cn, toolnames) for cn in summary_df.columns]
        metric_param_names = ["tool",'averaging', 'polytomy_treatment', 'rooting', 'metric']

        # Remove unused columns (any not in param_cols, or a tree distance measure)
        response_cols = [cn for cn in summary_df.columns if cn.endswith("treedist")]
        summary_df = summary_df.loc[:,param_cols+response_cols]
        # Convert metric params to 
        summary_df = self.df_wide_to_long(summary_df, toolnames, metric_param_names)
        summary_df = self.mean_se(summary_df, param_cols+metric_param_names)
        summary_df.to_csv(self.data_file)

    def plot(self):
        from matplotlib.legend_handler import HandlerTuple
        df = self.dataset.data
        error_params = df[ERROR_COLNAME].unique()
        sample_sizes = df.sample_size.unique()
        metrics = ARG_metrics.get_metric_names()
        topology_only_metrics = [m for m in metrics if not m.startswith('w')]
        fig, axes = pyplot.subplots(len(topology_only_metrics),
            len(error_params), figsize=(4*len(error_params), 15), sharey='row')
        for j, metric in enumerate(topology_only_metrics):
            for k, error in enumerate(error_params):
                ax = axes[j][k] if len(error_params)>1 else axes[j]
                ax.set_xscale('log')
                if j == 0:
                    ax.set_title(self.error_label(error))
                if k == 0:
                    ax.set_ylabel(getattr(self, 'y_axis_label', self.metric_titles[metric]))
                if j == len(topology_only_metrics) - 1:
                    ax.set_xlabel("Mutation rate")
                for n, fillstyle in zip(sample_sizes, self.fillstyles):
                    df_s = df[np.logical_and(df.sample_size == n, df[ERROR_COLNAME] == error)]
                    group = df_s.groupby(["mutation_rate"])
                    mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
                    for tool_and_metrics_param, setting in self.tools_and_metrics_params.items():
                        if metric.startswith("RF") \
                            and tool_and_metrics_param.startswith(TSINFER) \
                            and (int(tool_and_metrics_param.rsplit("_",1)[1]) & METRICS_POLYTOMIES_BREAK == 0):
                                # RF metrics are not well behaved for trees with polytomies (i.e. tsinfer trees)
                                # so we should omit tsinfer cases where METRICS_POLYTOMIES_BREAK == 0
                            continue
                        colname = tool_and_metrics_param + "_" + metric
                        if colname in df.columns:
                            ax.errorbar(
                                [m['mu'] for m in mean_sem],
                                [m['mean'][tool_and_metrics_param + "_" + metric] for m in mean_sem],
                                yerr=[m['sem'][tool_and_metrics_param + "_" + metric] for m in mean_sem] \
                                     if getattr(self, 'error_bars', False) else None,
                                linestyle=setting["linestyle"],
                                fillstyle=fillstyle,
                                color=setting["col"],
                                marker=None if len(sample_sizes)==1 else setting['mark'],
                                elinewidth=1)
                ax.set_ylim(ymin=0)
                ax.axvline(x=df.recombination_rate.unique()[0], 
                    color = 'gray', zorder=-1, linestyle=":", linewidth=1)
                ax.text(df.recombination_rate.unique()[0], ax.get_ylim()[1]/40, r'$\mu=\rho$',
                    va = "bottom",  ha="right", color = 'gray', rotation=90)
                ax.get_yaxis().set_label_coords(-0.15,0.5)

        artists = [
            pyplot.Line2D((0,1),(0,0), color= setting["col"],
                linewidth=2, linestyle=setting["linestyle"],
                marker = None if len(sample_sizes)==1 else setting['mark'])
            for tool,setting in self.tools_and_metrics_params.items()]
        tool_labels = [(l.replace("_0", "").replace("_2"," breaking polytomies") if l.startswith(TSINFER) else l.replace("_0", "")) 
            for l in self.tools_and_metrics_params.keys()]
        axes[0][0].legend(
            artists, tool_labels, numpoints=1, labelspacing=0.1)
            
            #numpoints=3, loc="upper right")
            # bbox_to_anchor=(0.0, 0.1))
        # ax = pyplot.gca().add_artist(first_legend)
        maintitle = "Average tree distance for neutral simulations"
        if len(sample_sizes)>1:
            artists = [
                tuple([
                    pyplot.Line2D(
                    (0,0),(0,0), color=setting['col'], fillstyle=fillstyle, marker=setting['mark'], linestyle='None')
                    for tool, setting in self.tools_and_metrics_params.items()])
                for fillstyle in self.fillstyles]
            axes[0][-1].legend(
                artists, ["Sample size = {}".format(n) for n in [15,20]],
                loc="upper right", numpoints=1, handler_map={tuple: HandlerTuple(ndivide=None, pad=2)})
        else:
            maintitle += " of {} samples".format(sample_sizes[0])
        maintitle += " spanning "  + ",".join(["{:g}kb\n".format(x/1e3) for x in df.length.unique()])
        if any([isinstance(Ne, str) for Ne in df.Ne.unique()]):
            maintitle += "(Model{}=".format('s' if len(df.Ne.unique()) > 1 else "") \
                + ",".join(["“{}”".format(x.replace("."," ")) for x in df.Ne.unique()])
        else:
            maintitle += "($N_e$=" + ",".join(["{}".format(x) for x in df.Ne.unique()])
        maintitle += r"; $\rho="+ ",".join(["{}".format(latex_float(x)) for x in df.recombination_rate.unique()])
        maintitle += "$)"
        fig.suptitle(maintitle, fontsize=16)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        self.savefig(fig)

class MetricsAllToolsAccuracyFigure(MetricsAllToolsFigure):
    """
    Show the metrics tending to 0 as mutation rate increases
    """
    datasetClass = AllToolsAccuracyDataset
    name = "metrics_all_tools_accuracy"
    error_bars = True
    
class MetricsAllToolsAccuracyDemographyFigure(MetricsAllToolsFigure):
    """
    Simple figure that shows all the metrics at the same time for
    a genome under a more complex demographic model (the Gutenkunst 
    Out Of Africa model), as mutation rate increases to high values
    """
    datasetClass = AllToolsAccuracyWithDemographyDataset
    name = "metrics_all_tools_accuracy_demography"


class MetricAllToolsFigure(Figure):
    """
    Superclass of the metric all tools figure. Each subclass should be a
    single figure for a particular metric.
    """
    datasetClass = AllToolsDataset
    fillstyles = ['full', 'none']

    def plot(self):
        df = self.dataset.data
        metric = self.metric
        error_params = df[ERROR_COLNAME].unique()
        sample_sizes = df.sample_size.unique()

        fig, axes = pyplot.subplots(1, len(error_params),
            figsize=getattr(self,'figsize',(12, 6)), sharey=True)
        lines = []
        for k, error in enumerate(error_params):
            ax = axes[k] if len(error_params)>1 else axes
            ax.set_title(self.error_label(error))
            ax.set_xlabel("Mutation rate")
            ax.set_xscale('log')
            if k == 0:
                ax.set_ylim(self.ylim)
                ax.set_ylabel(getattr(self, 'y_axis_label', self.metric_titles[metric]))
            for n, fillstyle in zip(sample_sizes, self.fillstyles):
                df_s = df[np.logical_and(df.sample_size == n, df[ERROR_COLNAME] == error)]
                group = df_s.groupby(["mutation_rate"])
                mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
                for tool_and_metrics_param,setting in self.tools_and_metrics_params.items():
                    if metric.startswith("RF") \
                        and tool_and_metrics_param.startswith(TSINFER) \
                        and (int(tool_and_metrics_param.rsplit("_",1)[1]) & METRICS_POLYTOMIES_BREAK == 0):
                            # RF metrics are not well behaved for trees with polytomies (i.e. tsinfer trees)
                            # so we should omit tsinfer cases where METRICS_POLYTOMIES_BREAK == 0
                        continue
                    ax.errorbar(
                        [m['mu'] for m in mean_sem],
                        [m['mean'][tool_and_metrics_param + "_" + metric] for m in mean_sem],
                        yerr=[m['sem'][tool_and_metrics_param + "_" + metric] for m in mean_sem] \
                             if getattr(self, 'error_bars', False) else None,
                        linestyle=setting["linestyle"],
                        fillstyle=fillstyle,
                        color=setting["col"],
                        marker=setting["mark"],
                        elinewidth=1)
            ax.set_ylim(ymin=0)
            ax.axvline(x=df.recombination_rate.unique()[0], 
                color = 'gray', zorder=-1, linestyle=":", linewidth=1)
            ax.text(df.recombination_rate.unique()[0], ax.get_ylim()[1]/40, r'$\mu=\rho$',
                va = "bottom",  ha="right", color = 'gray', rotation=90)

        # Create legends from custom artists
        artists = [
            pyplot.Line2D((0,1),(0,0), color= setting["col"], fillstyle=self.fillstyles[0],
                marker= setting["mark"], linestyle=setting["linestyle"])
            for tool,setting in self.tools_and_metrics_params.items()]
        tool_labels = [(l.replace("_0", "").replace("_2"," breaking polytomies") if l.startswith(TSINFER) else l.replace("_0", "")) 
            for l in self.tools_and_metrics_params.keys()]
        first_legend = axes[0].legend(
            artists, tool_labels, numpoints=1, labelspacing=0.1, loc="upper right")
            # bbox_to_anchor=(0.0, 0.1))
        # ax = pyplot.gca().add_artist(first_legend)
        if len(sample_sizes)>1:
            artists = [
                pyplot.Line2D(
                    (0,0),(0,0), color="black", fillstyle=fillstyle, linewidth=2)
                for n, fillstyle in zip(sample_sizes, self.fillstyles)]
            axes[-1].legend(
                artists, ["Sample size = {}".format(n) for n, fillstyle in zip(sample_sizes, self.fillstyles)],
                loc="upper right")
        self.savefig(fig)

class MetricAllToolsAccuracyFigure(MetricAllToolsFigure):
    """
    Superclass of the metric all tools figure. Each subclass should be a
    single figure for a particular metric.
    """
    datasetClass = AllToolsAccuracyDataset

class RFRootedAllToolsFigure(MetricAllToolsFigure):
    name = "rf_rooted_all_tools"
    metric = "RFrooted"
    ylim = None


class KCRootedAllToolsFigure(MetricAllToolsFigure):
    name = "kc_rooted_all_tools"
    metric = "KCrooted"
    ylim = (0, 24)
    error_bars = True

class KCRootedAllToolsBasicFigure(MetricAllToolsFigure):
    """For the publication: remove the polytomy-breaking line, and make symbols small and solid"""
    name = "kc_rooted_all_tools_basic"
    metric = "KCrooted"
    ylim = (0, 24)
    fillstyles = ['full']
    figsize = (9,4.5)
    error_bars = True
    y_axis_label="Average distance from true trees (mean $\pm$ s.e.)"
    def __init__(self):
        super().__init__()
        self.tools_and_metrics_params = {k:v for k,v in self.tools_and_metrics_params.items() \
            if k.endswith("_" + str(METRICS_LOCATION_ALL | METRICS_POLYTOMIES_LEAVE))}        

class RFRootedAllToolsAccuracyFigure(MetricAllToolsAccuracyFigure):
    name = "rf_rooted_all_tools_accuracy"
    metric = "RFrooted"
    ylim = None


class KCRootedAllToolsAccuracyFigure(MetricAllToolsAccuracyFigure):
    name = "kc_rooted_all_tools_accuracy"
    metric = "KCrooted"
    ylim = (0, 40)
    error_bars = True

class MetricAllToolsAccuracyDemographyFigure(MetricAllToolsFigure):
    """
    Superclass of the metric all tools figure for demographic simulations. 
    Each subclass should be a single figure for a particular metric.
    """
    datasetClass = AllToolsAccuracyWithDemographyDataset

class KCRootedAllToolsAccuracyDemographyFigure(MetricAllToolsAccuracyDemographyFigure):
    name = "kc_rooted_all_tools_accuracy_demography"
    metric = "KCrooted"
    ylim = (0, 20)
    error_bars = True

class CputimeAllToolsFigure(Figure):
    """
    This figure is useful because we can only really get the CPU times
    for all four methods in the same scale for these tiny examples.
    We can show that ARGWeaver and RentPlus are much slower than tsinfer
    and FastARG here and compare tsinfer and FastARG more thoroughly
    in a dedicated figure.
    """
    name = "cputime_all_tools"
    datasetClass = AllToolsPerformanceDataset
    error_bars = True

    def plot(self):
        df = self.dataset.data
        sample_sizes = df.sample_size.unique()

        fig, (ax_hi, ax_lo) = pyplot.subplots(2, 1, sharex=True)
        lines = []
        noerror = [e for e in df[ERROR_COLNAME].unique() if self.error_label(e, None) is None]
        assert len(noerror) == 1
        noerror = noerror[0]
        ax_lo.set_xlabel("Sample Size")
        ax_hi.set_ylabel("CPU time (hours)")
        #ax_lo.set_xlim(sample_sizes.min(), sample_sizes.max())

        # zoom-in / limit the view to different portions of the data
        #ax_hi.set_ylim(0.6, 2.7)  # outliers only
        ax_lo.set_ylim(0, 0.005)  # most of the data
        #ax_hi.set_ylim(0.01, 3)  # outliers only
        #ax_lo.set_ylim(0, 0.002)  # most of the data
        
        # hide the spines between ax and ax2
        ax_hi.spines['bottom'].set_visible(False)
        ax_lo.spines['top'].set_visible(False)
        ax_hi.xaxis.tick_top()
        ax_hi.tick_params(labeltop='off')  # don't put tick labels at the top
        ax_lo.xaxis.tick_bottom()

        df_s = df[np.logical_and(
            df[ERROR_COLNAME] == noerror, 
            df['length'] == self.dataset.fixed_length)]
        group = df_s.groupby(["sample_size"])
        mean_sem = [{'sz':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
        for tool, setting in self.tools.items():
            for ax in (ax_lo, ax_hi):
                ax.errorbar(
                    [m['sz'] for m in mean_sem],
                    [m['mean'][tool + "_" + "cputime"]/60/60 for m in mean_sem],
                    yerr=[m['sem'][tool + "_" + "cputime"]/60/60 for m in mean_sem] \
                         if getattr(self, 'error_bars', False) else None,
                    color=setting["col"],
                    marker=setting["mark"],
                    elinewidth=1,
                    label=tool)
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax_hi.transAxes, color='k', clip_on=False)
        ax_hi.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax_hi.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
        
        kwargs.update(transform=ax_lo.transAxes)  # switch to the bottom axes
        ax_lo.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax_lo.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        ax_hi.legend(loc="upper left")
        self.savefig(fig)

class MetricsAllToolsAccuracySweepFigure(Figure):
    """
    Simple figure that shows all the metrics at the same time for
    a genome under a selective sweep
    """
    datasetClass = AllToolsAccuracyWithSelectiveSweepDataset
    name = "metrics_all_tools_accuracy_sweep"
    error_bars = True
    fillstyles = ['full', 'none']

    def plot(self):
        df = self.dataset.data
        error_params = df[ERROR_COLNAME].unique()
        output_freqs = df[[SELECTED_FREQ_COLNAME, SELECTED_POSTGEN_COLNAME]].drop_duplicates()
        sample_sizes = df.sample_size.unique()

        metrics = ARG_metrics.get_metric_names()
        topology_only_metrics = [m for m in metrics if not m.startswith('w')] #don't use the weighted metrics
        figures = [] #save figures for putting onto a multi-page pdf
        for error in error_params:
            fig, axes = pyplot.subplots(len(topology_only_metrics),
                len(output_freqs), figsize=(20, 20))
            fig.suptitle(self.error_label(error))
            for j, metric in enumerate(topology_only_metrics):
                for k, output_data in enumerate(output_freqs.itertuples()):
                    freq = output_data.output_frequency
                    gens = output_data.output_after_generations
                    ax = axes[j][k]
                    ax.set_xscale('log')
                    if j == 0:
                        ax.set_title("Swept variant @ freq {}{}".format(
                            freq, "+{} gens".format(int(gens)) if gens else ""))
                    if k == 0:
                        ax.set_ylabel(getattr(self, 'y_axis_label', self.metric_titles[metric]))
                    if j == len(topology_only_metrics) - 1:
                        ax.set_xlabel("Mutation rate")
                    if np.isclose(freq, 1.0) and not gens:
                        # This is *at* fixation - set the plot background colour
                        ax.set_facecolor('0.9')
                    for n, fillstyle in zip(sample_sizes, self.fillstyles):
                        df_s = df[np.logical_and.reduce((
                            df.sample_size == n,
                            df[ERROR_COLNAME] == error,
                            df[SELECTED_FREQ_COLNAME] == freq,
                            df[SELECTED_POSTGEN_COLNAME] == gens))]
                        group = df_s.groupby(["mutation_rate"])
                        #NB pandas.DataFrame.mean and pandas.DataFrame.sem have skipna=True by default
                        mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
                        for tool_and_metrics_param, setting in self.tools_and_metrics_params.items():
                            if metric.startswith("RF") \
                                and tool_and_metrics_param.startswith(TSINFER) \
                                and (int(tool_and_metrics_param.rsplit("_",1)[1]) & METRICS_POLYTOMIES_BREAK == 0):
                                    # RF metrics are not well behaved for trees with polytomies (i.e. tsinfer trees)
                                    # so we should omit tsinfer cases where METRICS_POLYTOMIES_BREAK == 0
                                continue
                            if all(np.isnan(m['mean'][tool_and_metrics_param + "_" + metric]) for m in mean_sem):
                                #don't plot if all NAs
                                continue
                            ax.errorbar(
                                [m['mu'] for m in mean_sem],
                                [m['mean'][tool_and_metrics_param + "_" + metric] for m in mean_sem],
                                yerr=[m['sem'][tool_and_metrics_param + "_" + metric] for m in mean_sem] \
                                    if getattr(self, 'error_bars', False) else None,
                                color=setting["col"],
                                linestyle=setting["linestyle"],
                                marker=setting["mark"],
                                fillstyle=fillstyle,
                                elinewidth=1)
                    ax.set_ylim(ymin=0)
                    ax.get_yaxis().set_label_coords(-0.15,0.5)
                    ax.axvline(x=df.recombination_rate.unique()[0], 
                        color = 'gray', zorder=-1, linestyle=":", linewidth=1)
                    ax.text(df.recombination_rate.unique()[0], ax.get_ylim()[1]/40, r'$\mu=\rho$',
                        va = "bottom",  ha="right", color = 'gray', rotation=90)

            # Create legends from custom artists
            artists = [
                pyplot.Line2D((0,1),(0,0), color= setting["col"], fillstyle=self.fillstyles[0],
                    marker= setting["mark"], linestyle=setting["linestyle"])
                for tool,setting in self.tools_and_metrics_params.items()]
            tool_labels = [(l.replace("_0", "").replace("_2"," breaking polytomies") if l.startswith(TSINFER) else l.replace("_0", "")) 
                for l in self.tools_and_metrics_params.keys()]
            first_legend = axes[0][0].legend(
                artists, tool_labels, numpoints=1, labelspacing=0.1, loc="upper right")
            if len(sample_sizes)>1:
                artists = [
                    pyplot.Line2D(
                        (0,0),(0,0), color="black", fillstyle=fillstyle, linewidth=2)
                    for n, fillstyle in zip(sample_sizes, self.fillstyles)]
                axes[0][-1].legend(
                    artists, ["Sample size = {}".format(n) for n, fillstyle in zip(sample_sizes, self.fillstyles)],
                    loc="upper right")
            figures.append(fig)
        self.savefig(*figures)



class MetricAllToolsAccuracySweepFigure(Figure):
    """
    Superclass of the metric all tools figure for simulations with selection. 
    Each subclass should be a single figure for a particular metric.
    """
    datasetClass = AllToolsAccuracyWithSelectiveSweepDataset
    datasetClass = AllToolsAccuracyWithSelectiveSweepDataset
    name = "metric_all_tools_accuracy_sweep"
    error_bars = True
    fillstyles = ['full', 'none']

    def plot(self):
        df = self.dataset.data
        metric = self.metric
        error_params = df[ERROR_COLNAME].unique()
        output_freqs = df[[SELECTED_FREQ_COLNAME, SELECTED_POSTGEN_COLNAME]].drop_duplicates()
        sample_sizes = df.sample_size.unique()

        fig, axes = pyplot.subplots(len(output_freqs), len(error_params),
            figsize=(7*len(error_params), 3*len(output_freqs)))
        for j, output_data in enumerate(output_freqs.itertuples()):
            for k, error in enumerate(error_params):
                freq = output_data.output_frequency
                gens = output_data.output_after_generations
                ax = axes[j][k]
                ax.set_xscale('log')
                if j == 0:
                    ax.set_title(self.error_label(error))
                if j == len(output_freqs) - 1:
                    ax.set_xlabel("Neutral mutation rate")
                if k == 0:
                    ax.set_ylabel(getattr(self, 'y_axis_label', self.metric_titles[metric].replace("Kendall-Colijn", "KC")) + 
                        " @ {}{}".format(
                            "fixation " if np.isclose(freq, 1.0) else "freq {}".format(freq),
                            "+{} gens".format(int(gens)) if gens else ""))
                if np.isclose(freq, 1.0) and not gens:
                    # This is *at* fixation - set the plot background colour
                    ax.set_facecolor('0.9')
                for n, fillstyle in zip(sample_sizes, self.fillstyles):
                    df_s = df[np.logical_and.reduce((
                        df.sample_size == n,
                        df[ERROR_COLNAME] == error,
                        df[SELECTED_FREQ_COLNAME] == freq,
                        df[SELECTED_POSTGEN_COLNAME] == gens))]
                    group = df_s.groupby(["mutation_rate"])
                    #NB pandas.DataFrame.mean and pandas.DataFrame.sem have skipna=True by default
                    mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
                    for tool_and_metrics_param, setting in self.tools_and_metrics_params.items():
                        if metric.startswith("RF") \
                            and tool_and_metrics_param.startswith(TSINFER) \
                            and (int(tool_and_metrics_param.rsplit("_",1)[1]) & METRICS_POLYTOMIES_BREAK == 0):
                                # RF metrics are not well behaved for trees with polytomies (i.e. tsinfer trees)
                                # so we should omit tsinfer cases where METRICS_POLYTOMIES_BREAK == 0
                            continue
                        if all(np.isnan(m['mean'][tool_and_metrics_param + "_" + metric]) for m in mean_sem):
                            #don't plot if all NAs
                            continue
                        ax.errorbar(
                            [m['mu'] for m in mean_sem],
                            [m['mean'][tool_and_metrics_param + "_" + metric] for m in mean_sem],
                            yerr=[m['sem'][tool_and_metrics_param + "_" + metric] for m in mean_sem] \
                                if getattr(self, 'error_bars', False) else None,
                            color=setting["col"],
                            linestyle=setting["linestyle"],
                            marker=setting["mark"],
                            fillstyle=fillstyle,
                            elinewidth=1)
                ax.set_ylim(getattr(self,'ylim',0))
                ax.axvline(x=df.recombination_rate.unique()[0], 
                    color = 'gray', zorder=-1, linestyle=":", linewidth=1)
                ax.text(df.recombination_rate.unique()[0], ax.get_ylim()[1]/40, r'$\mu=\rho$',
                    va = "bottom",  ha="right", color = 'gray', rotation=90)

        # Create legends from custom artists
        artists = [
            pyplot.Line2D((0,1),(0,0), color= setting["col"], fillstyle=self.fillstyles[0],
                marker= setting["mark"], linestyle=setting["linestyle"])
            for tool,setting in self.tools_and_metrics_params.items()]
        tool_labels = [(l.replace("_0", "").replace("_2"," breaking polytomies") if l.startswith(TSINFER) else l.replace("_0", "")) 
            for l in self.tools_and_metrics_params.keys()]
        first_legend = axes[0][0].legend(
            artists, tool_labels, numpoints=1, labelspacing=0.1, loc="upper right")
        if len(sample_sizes)>1:
            artists = [
                pyplot.Line2D(
                    (0,0),(0,0), color="black", fillstyle=fillstyle, linewidth=2)
                for n, fillstyle in zip(sample_sizes, self.fillstyles)]
            axes[0][-1].legend(
                artists, ["Sample size = {}".format(n) for n, fillstyle in zip(sample_sizes, self.fillstyles)],
                loc="upper right")
        self.savefig(fig)

class KCRootedAllToolsAccuracySweepFigure(MetricAllToolsAccuracySweepFigure):
    name = "kc_rooted_all_tools_accuracy_sweep"
    metric = "KCrooted"
    ylim = (0, 30)
    error_bars = True


class MetricsSubsamplingFigure(Figure):
    """
    Figure that shows whether increasing sample size helps with the accuracy of
    reconstructing the ARG for a fixed subsample. We only use tsinfer for this.
    """
    datasetClass = SubsamplingDataset
    name = "metrics_subsampling"
    error_bars=True

    def plot(self):
        figures = [] #save figures for putting onto a multi-page pdf
        df = self.dataset.data
        lengths = df.length.unique()
        restrict_sample_size_comparison = df.restrict_sample_size_comparison.unique()
        assert len(restrict_sample_size_comparison) == 1
        col="black"
        metrics = ARG_metrics.get_metric_names()
        topology_only_metrics = [m for m in metrics if not m.startswith('w')]
        for error in df[ERROR_COLNAME].unique():
            fig, axes = pyplot.subplots(len(topology_only_metrics),
                len(lengths), figsize=(12, 16))
            for j, metric in enumerate(topology_only_metrics):
                for k, length in enumerate(lengths):
                    ax = axes[j][k]
                    if j == 0:
                        ax.set_title("Sequence length = {} Kb".format(int(length/1000)))
                    if k == 0:
                        ax.set_ylabel(getattr(self, 'y_axis_label', metric + " metric"))
                    if j == len(topology_only_metrics) - 1:
                        ax.set_xlabel("Original sample size")
                    for tool_and_metrics_param, setting in self.tools_and_metrics_params.items():
                        colname = tool_and_metrics_param + "_" + metric
                        if colname in df.columns \
                            and metric.startswith("RF") \
                            and tool_and_metrics_param.startswith(TSINFER) \
                            and (int(tool_and_metrics_param.rsplit("_",1)[1]) & METRICS_POLYTOMIES_BREAK == 0):
                                # RF metrics are not well behaved for trees with polytomies (i.e. tsinfer trees)
                                # so we should omit tsinfer cases where METRICS_POLYTOMIES_BREAK == 0
                                continue
                        dfp = df[np.logical_and(df.length == length, df[ERROR_COLNAME] == error)]
                        group = dfp.groupby(["subsample_size"])
                        mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
                        #NB pandas.DataFrame.mean and pandas.DataFrame.sem have skipna=True by default
                        #print(mean_sem)
                        if getattr(self, 'error_bars', None):
                            yerr=[m['sem'][colname] for m in mean_sem]
                        else:
                            yerr = None
                        ax.errorbar(
                            [m['mu'] for m in mean_sem],
                            [m['mean'][colname] for m in mean_sem],
                            yerr=yerr,
                            color= setting['col'],
                            linestyle=setting['linestyle'],
                            marker=setting['mark'],
                            #elinewidth=1
                            )
            artists = [
                pyplot.Line2D((0,1),(0,0), color= setting["col"],
                    marker= setting["mark"], linestyle=setting["linestyle"])
                for tool,setting in self.tools_and_metrics_params.items()]
            tool_labels = [(l.replace("_1", "").replace("_3"," breaking polytomies") if l.startswith(TSINFER) else l.replace("_0", "")) 
                for l in self.tools_and_metrics_params.keys()]
            first_legend = axes[0][-1].legend(
                artists, tool_labels, numpoints=1, labelspacing=0.1, loc="lower right")
            fig.suptitle('ARG metric for trees subsampled down to {} tips'.format(
                restrict_sample_size_comparison[0], self.error_label(error)))
            figures.append(fig)
        self.savefig(*figures)


class MetricARGweaverParametersFigure(Figure):
    """
    See the effect of burnin time and number of timeslices on the accuracy of ARGweaver
    as compared to TSinfer, when looking at tree metrics.
    """
    datasetClass = ARGweaverParamChangesDataset

    def plot(self):
        df = self.dataset.data
        metric = self.metric
        tools = self.dataset.tools
        AW_burnin = df.ARGweaver_burnin.unique()
        AW_discrete_timeslices = df.ARGweaver_ntimes.unique()

        # TODO move this into the superclass so that we have consistent styling.
        linestyles = [":","-.", "-"]
        fig, axes = pyplot.subplots(1, 3,
            figsize=getattr(self,'figsize',(12, 6)), sharey=True)
        lines = []
        for k, ntimes in enumerate(AW_discrete_timeslices):
            ax = axes[k]
            ax.set_title("ARGweaver timeslices = {}".format(ntimes))
            ax.set_xlabel("Mutation rate")
            ax.set_xscale('log')
            if k == 0:
                ax.set_ylabel(getattr(self, 'y_axis_label', metric + " metric"))

            #get the only tsinfer data for this selection (regardless of ntimes)
            tool = TSINFER
            df_s = df[np.logical_not(df[tool + "_" + metric].isnull())]
            group = df_s.groupby(["mutation_rate"])
            mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
            if getattr(self, 'error_bars', None):
                yerr=[m['sem'][tool + "_" + metric] for m in mean_sem]
            else:
                yerr = None
            ax.errorbar(
                [m['mu'] for m in mean_sem],
                [m['mean'][tool + "_" + metric] for m in mean_sem],
                yerr=yerr,
                linestyle=linestyles[0],
                color=self.tools[tool]["col"],
                marker=self.tools[tool]["mark"],
                elinewidth=1)


            for n, linestyle in zip(AW_burnin, linestyles):
                tool = ARGWEAVER
                df_s = df[np.logical_and(df.ARGweaver_burnin == n, df.ARGweaver_ntimes == ntimes)]
                group = df_s.groupby(["mutation_rate"])
                mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
                if getattr(self, 'error_bars', None):
                    yerr=[m['sem'][tool + "_" + metric] for m in mean_sem]
                else:
                    yerr = None
                ax.errorbar(
                    [m['mu'] for m in mean_sem],
                    [m['mean'][tool + "_" + metric] for m in mean_sem],
                    yerr=yerr,
                    linestyle=linestyle,
                    color=self.tools[tool]["col"],
                    marker=self.tools[tool]["mark"],
                    elinewidth=1)

        axes[0].set_ylim(self.ylim)

        # Create legends from custom artists
        artists = [
            pyplot.Line2D((0,1),(0,0), color=self.tools[tool]["col"],
                marker=self.tools[tool]["mark"], linestyle='')
            for tool in tools]
        first_legend = axes[0].legend(
            artists, tools, numpoints=3, loc="upper center")
            # bbox_to_anchor=(0.0, 0.1))
        # ax = pyplot.gca().add_artist(first_legend)
        artists = [
            pyplot.Line2D(
                (0,0),(0,0), color="black", linestyle=linestyle, linewidth=2)
            for linestyle in linestyles]
        axes[-1].legend(
            artists, ["Burn in = {} gens".format(n) for n in AW_burnin],
            loc="upper center")
        self.savefig(fig)

class KCRootedARGweaverParametersFigure(MetricARGweaverParametersFigure):
    name = "kc_rooted_argweaver_params"
    metric = "KCrooted"
    ylim = (0, 4)
    error_bars = True

class RFRootedARGweaverParametersFigure(MetricARGweaverParametersFigure):
    name = "rf_rooted_argweaver_params"
    metric = "RFrooted"
    ylim = None
    #ylim = (0, 4)
    error_bars = True


class TsinferPerformanceLengthSamplesFigure(Figure):
    """
    Superclass for the performance metrics figures. Each of these figures
    has two panels; one for scaling by sequence length and the other
    for scaling by sample size. Different lines are given for each
    of the different combinations of tsinfer parameters
    """
    datasetClass = TsinferPerformanceDataset
    plotted_column = None
    y_axis_label = None

    def plot(self):
        df = self.dataset.data
        # Rescale the length to MB
        df.length /= 10**6
        # Set statistics to the ratio of observed over expected
        source_colour = "red"
        inferred_colour = self.tools[TSINFER]["col"]
        recombination_linestyles = [':', '-', '--']
        recombination_rates = df.recombination_rate.unique()
        mutation_rates = self.datasetClass.between_sim_params['mutation_rate']
        assert len(recombination_linestyles) >= len(recombination_rates)
        assert len(mutation_rates) ==1
        mu = mutation_rates[0]
        #inferred_markers =    {False:'s',True:'b'}
        fig, (ax1, ax2) = pyplot.subplots(1, 2, figsize=(12, 6), sharey=True)
        ax1.set_title("Fixed number of chromosomes ({})".format(self.datasetClass.fixed_sample_size))
        ax1.set_xlabel("Sequence length (MB)")
        ax1.set_ylabel(self.y_axis_label)
        for rho_index, rho in enumerate(recombination_rates):
            dfp = df[np.logical_and.reduce((
                df.sample_size == self.datasetClass.fixed_sample_size,
                df.recombination_rate == rho))]
            group = dfp.groupby(["length"])
            #NB pandas.DataFrame.mean and pandas.DataFrame.sem have skipna=True by default
            mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
            ax1.errorbar(
                [m['mu'] for m in mean_sem],
                [m['mean'][self.plotted_column] for m in mean_sem],
                yerr=[m['sem'][self.plotted_column] for m in mean_sem] if getattr(self, 'error_bars', False) else None,
                linestyle= recombination_linestyles[rho_index],
                color=inferred_colour,
                #elinewidth=1
                )


        ax2.set_title("Fixed sequence length ({:g} Mb)".format(self.datasetClass.fixed_length / 10**6))
        ax2.set_xlabel("Sample size")
        ax2.set_ylabel(self.y_axis_label)
        for rho_index, rho in enumerate(recombination_rates):
            dfp = df[np.logical_and.reduce((
                df.length == self.datasetClass.fixed_length / 10**6,
                df.recombination_rate == rho))]
            group = dfp.groupby(["sample_size"])
            #NB pandas.DataFrame.mean and pandas.DataFrame.sem have skipna=True by default
            mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
            ax2.errorbar(
                [m['mu'] for m in mean_sem],
                [m['mean'][self.plotted_column] for m in mean_sem],
                yerr=[m['sem'][self.plotted_column] for m in mean_sem] if getattr(self, 'error_bars', False) else None,
                linestyle= recombination_linestyles[rho_index],
                color=inferred_colour,
                #elinewidth=1
                )
        params = [
            pyplot.Line2D(
                (0,0),(0,0), color= inferred_colour,
                linestyle=recombination_linestyles[rho_index], linewidth=2)
            for rho_index, rho in enumerate(recombination_rates)]
        ax1.legend(
            params, [r"$\rho$ = {}".format("$\mu$" if rho==mu else r"{:g}$\mu$".format(rho/mu) if rho>mu else r"$\mu$/{:g}".format(mu/rho))
                for rho_index, rho in enumerate(recombination_rates)],
            loc="lower right", fontsize=10, title="Relative rate of\nrecombination")
    
        pyplot.suptitle(r'Tsinfer large dataset performance for $\mu$=${}$'.format(latex_float(mu)))
        self.savefig(fig)


class TsinferEdgesPerformanceFigure(TsinferPerformanceLengthSamplesFigure):
    name = "tsinfer_edges_ln"
    plotted_column = "metric"
    y_axis_label = "inferred_edges / real_edges"
    def plot(self):
        self.dataset.data[self.plotted_column] = (
            self.dataset.data["tsinfer_edges"] / self.dataset.data["edges"])
        TsinferPerformanceLengthSamplesFigure.plot(self)


class TsinferFileSizePerformanceFigure(TsinferPerformanceLengthSamplesFigure):
    name = "tsinfer_filesize_ln"
    plotted_column = "metric"
    y_axis_label = "inferred_filesize / real_filesize"
    y_axis_label = "File size relative to original"
    def plot(self):
        self.dataset.data[self.plotted_column] = (
            self.dataset.data["tsinfer_ts_filesize"] / self.dataset.data["ts_filesize"])
        TsinferPerformanceLengthSamplesFigure.plot(self)


class CompressionPerformanceFigure(TsinferPerformanceLengthSamplesFigure):
    name = "tsinfer_compression_ln"
    plotted_column = "metric"
    y_axis_label = "inferred_filesize / real_filesize"
    y_axis_label = "Compression factor relative to vcf.gz"
    def plot(self):
        self.dataset.data[self.plotted_column] = (
             self.dataset.data["vcfgz_filesize"] / self.dataset.data["tsinfer_ts_filesize"])
        TsinferPerformanceLengthSamplesFigure.plot(self)


class TsinferPerformanceSizesSamplesFigure(Figure):
    """
    Class for the performance metrics against sites as well as length
    (the first is the same as the LHS TsinferEdgesPerformanceFigure,
    but the second is the same using sites instead)
    """
    name="tsinfer_edges_sn"
    datasetClass = TsinferPerformanceDataset
    plotted_column = None
    y_axis_label = None
    error_bars = True

    def plot(self):
        df = self.dataset.data
        self.dataset.data[self.plotted_column] = self.dataset.data["tsinfer_edges"] / self.dataset.data["edges"]
        # Rescale the length to MB
        df.length /= 10**6
        # Set statistics to the ratio of observed over expected
        source_colour = "red"
        inferred_colour = self.tools[TSINFER]["col"]
        recombination_linestyles = [':', '-', '--']
        recombination_rates = df.recombination_rate.unique()
        mutation_rates = self.datasetClass.between_sim_params['mutation_rate']
        assert len(recombination_linestyles) >= len(recombination_rates)
        assert len(mutation_rates) ==1
        mu = mutation_rates[0]
        inferred_linestyles = {False:':',True:'-'}
        fig, (ax1, ax2) = pyplot.subplots(1, 2, figsize=(12, 6), sharey=True)
        ax1.set_title("Fixed number of chromosomes ({})".format(self.datasetClass.fixed_sample_size))
        ax1.set_xlabel("Sequence length (MB)")
        ax1.set_ylabel(self.y_axis_label)
        for shared_breakpoint in df.tsinfer_srb.unique():
            dfp = df[np.logical_and.reduce((
                df.sample_size == self.datasetClass.fixed_sample_size,
                df.tsinfer_srb == shared_breakpoint))]
            group = dfp.groupby(["length"])
            #NB pandas.DataFrame.mean and pandas.DataFrame.sem have skipna=True by default
            mean_sem = [{'mu':g, 'mean':data.mean(), 'sem':data.sem()} for g, data in group]
            ax1.errorbar(
                [m['mu'] for m in mean_sem],
                [m['mean'][self.plotted_column] for m in mean_sem],
                yerr=[m['sem'][self.plotted_column] for m in mean_sem] if getattr(self, 'error_bars') else None,
                linestyle=inferred_linestyles[shared_breakpoint],
                color=inferred_colour,
                )

        ax2.set_title("Fixed number of chromosomes ({})".format(self.datasetClass.fixed_sample_size))
        ax2.set_xlabel("Number of variable sites")
        ax2.set_ylabel(self.y_axis_label)
        for shared_breakpoint in df.tsinfer_srb.unique():
            dfp = df[np.logical_and.reduce((
                df.length == self.datasetClass.fixed_length / 10**6,
                df.tsinfer_srb == shared_breakpoint))]
            ax2.plot(
                dfp['n_sites'],
                dfp[self.plotted_column],
                color=inferred_colour,
                linestyle="",marker="o")

        pyplot.suptitle('Tsinfer large dataset performance for mu={}'.format(mu))
        self.savefig(fig)


class FastargTsinferComparisonFigure(Figure):
    """
    Superclass for the program comparison figures (comparing tsinfer with fastarg)
    Each figure has two panels; one for scaling by sequence length and the other
    for scaling by sample size.
    """
    datasetClass = FastargTsinferComparisonDataset

    def plot(self):
        df = self.dataset.data
        tools = self.dataset.tools_and_metrics.keys()

        # Rescale the length to Mb
        length_scale = 10**6
        df.length /= length_scale
        # Scale time to hours
        time_scale = 3600
        for tool in tools:
            df[tool + "_cputime"] /= time_scale
        # Scale memory to GiB
        for tool in tools:
            df[tool + "_memory"] /= 1024 * 1024 * 1024

        fig, (ax1, ax2) = pyplot.subplots(1, 2, sharey=True, figsize=(8, 5.5))

        dfp = df[df.sample_size == self.datasetClass.fixed_sample_size]
        group = dfp.groupby(["length"])
        group_mean = group.mean()

        for tool in tools:
            col = tool + "_" + self.plotted_column
            ax1.plot(group_mean[col], label=tool, color=self.tools[tool]["col"])
        ax1.legend(
            loc="upper left", numpoints=1, fontsize="small")

        ax1.set_xlabel("Length (Mb) for fixed sample of {}".format(self.datasetClass.fixed_sample_size))
        ax1.set_ylabel(self.y_label)

        dfp = df[df.length == self.datasetClass.fixed_length / length_scale]
        group = dfp.groupby(["sample_size"])
        group_mean = group.mean()

        for tool in tools:
            col = tool + "_" + self.plotted_column
            ax2.plot(group_mean[col], label=tool, color=self.tools[tool]["col"])

        ax2.set_xlabel("Sample size for fixed length of {} Mb".format(
            self.datasetClass.fixed_length / length_scale))

        # ax1.set_xlim(-5, 105)
        # ax1.set_ylim(-5, 250)
        # ax2.set_xlim(-5, 105)
        self.tweak(ax1, ax2)

        fig.tight_layout()

        # fig.text(0.19, 0.97, "Sample size = 1000")
        # fig.text(0.60, 0.97, "Sequence length = 50Mb")
        # pyplot.savefig("plots/simulators.pdf")
        self.savefig(fig)


class FastargTsinferComparisonTimeFigure(FastargTsinferComparisonFigure):
    name = "fastarg_tsinfer_comparison_time"
    plotted_column = "cputime"
    y_label = "CPU time (hours)"

    def tweak(self, ax1, ax2):
        ax1.set_ylim(0, 2.5)


class FastargTsinferComparisonMemoryFigure(FastargTsinferComparisonFigure):
    name = "fastarg_tsinfer_comparison_memory"
    plotted_column = "memory"
    y_label = "Memory (GiB)"

    def tweak(self, ax1, ax2):
        pass


def run_setup(cls, args):
    if args.processes > 1 and args.progress == True:
        raise ValueError("Due to limitations of the Python multiprocessing module" \
            " you can't show a tqdm progress-bar when setting up simulations in parallel." \
            " Try again without the -P (--progress) flag.")
    f = cls()
    f.setup(args)

def run_infer(cls, args):
    logging.info("Inferring {}".format(cls.name))
    f = cls()
    f.infer(
        args.processes, args.threads, force=args.force, metrics_only=args.metrics_only,
        specific_tool=args.tool, specific_row=args.row,
        flush_all=args.flush_all, show_progress=args.progress)

def run_plot(cls, args):
    f = cls()
    f.plot()

def run_summarize(cls, args):
    f = cls()
    f.summarize()


def main():
    def get_subclasses(cls):
        for subclass in cls.__subclasses__():
            yield from get_subclasses(subclass)
            yield subclass
    datasets = list(get_subclasses(Dataset))
    figures = list(get_subclasses(Figure))
    name_map = dict([(d.name, d) for d in datasets + figures])
    parser = argparse.ArgumentParser(
        description="Set up base data, generate inferred datasets, process datasets and plot figures.")
    parser.add_argument('--verbosity', '-v', action='count', default=0)
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    subparser = subparsers.add_parser('setup',
        help="Run simulations, outputting true histories & genome sequences for analysis" +
            "(created files will overwrite previous runs of the same name)")
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier', choices=sorted([d.name for d in datasets if d.name]))
    subparser.add_argument(
        "--processes", '-p', type=int, default=1,
        help="number of worker processes, e.g. 40")
    subparser.add_argument(
         '--replicates', '-r', type=int, help="number of replicates")
    subparser.add_argument(
         '--seed', '-s', type=int, help="use a non-default RNG seed")
    subparser.add_argument(
         '--hack_finite_sites', action='store_true',
         help="Mutations at the same (integer) location are superimposed, not shifted along")
    subparser.add_argument(
         '--progress',  "-P", action='store_true',
         help="Show a progress bar.", )
    subparser.set_defaults(func=run_setup)

    subparser = subparsers.add_parser('infer',
        help="Infer genealogical histories from previously generated sequence data," +
            " and generate statistics for comparison")
    subparser.add_argument(
        "--processes", '-p', type=int, default=1,
        help="number of worker processes, e.g. 40")
    subparser.add_argument(
        "--threads", '-t', type=int, default=1,
        help="number of threads per worker process (for supporting tools), e.g. 8")
    subparser.add_argument(
        "--tool", '-T', default=None,
        help="Only run this specific tool")
    subparser.add_argument(
        "--row", '-r', type=int, default=None,
        help="Only run for a specific row")
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier', choices=sorted([d.name for d in datasets if d.name]))
    subparser.add_argument(
         '--force',  "-f", action='store_true',
         help="redo all the inferences, even if we have already filled out some", )
    subparser.add_argument(
         '--metrics-only',  "-m", action='store_true',
         help="skip the inference step & just re-calculate metrics for already-run sims", )
    subparser.add_argument(
         '--progress',  "-P", action='store_true',
         help="Show a progress bar.", )
    subparser.add_argument(
         '--flush-all',  "-F", action='store_true',
         help="flush the result file after every result.", )
    subparser.set_defaults(func=run_infer)

    subparser = subparsers.add_parser('summarize',
        help="Create data files for later plotting")
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the figure identifier', choices=sorted([f.name for f in figures if f.name]))
    subparser.set_defaults(func=run_summarize)

    subparser = subparsers.add_parser('figure')
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the figure identifier', choices=sorted([f.name for f in figures if f.name]))
    subparser.set_defaults(func=run_plot)

    args = parser.parse_args()
    log_level = logging.WARNING
    if args.verbosity == 1:
        log_level = logging.INFO
    if args.verbosity >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(
        format='%(asctime)s %(message)s', level=log_level, stream=sys.stdout)

    # Create a new process group and become the leader.
    os.setpgrp()
    k = args.name[0]
    try:
        if k == "all":
            classes = datasets
            if args.func == run_plot:
                classes = figures
            for name, cls in name_map.items():
                if cls in classes:
                    args.func(cls, args)
        else:
            try:
                cls = name_map[k]
            except KeyError as e:
                e.args = (e.args[0] + ". Select from datasets={} or figures={}".format(
                        [d.name for d in datasets], [f.name for f in figures]),)
                raise
            args.func(cls, args)
    except KeyboardInterrupt:
        print("Interrupted! Trying to kill subprocesses")
        os.killpg(0, signal.SIGINT)

if __name__ == "__main__":
    main()
