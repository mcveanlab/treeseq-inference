"""
Code to run simulations, inference methods and generate all plots
in the paper.
"""

import argparse
import logging
import os.path
import shutil
import sys
import time
import glob
import collections
import itertools
import tempfile
import subprocess
import filecmp
import random
import multiprocessing
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import pandas as pd

# import the local copy of msprime in preference to the global one
curr_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(1,os.path.join(curr_dir,'..','msprime'))
import msprime
import msprime_extras
import msprime_fastARG
import msprime_ARGweaver
import tsinfer
import ARG_metrics

fastARG_executable = os.path.join(curr_dir,'..','fastARG','fastARG')
ARGweaver_executable = os.path.join(curr_dir,'..','argweaver','bin','arg-sample')
tsinfer_executable = os.path.join(curr_dir,'run_tsinfer.py')
#monkey-patch write.nexus into msprime
msprime.TreeSequence.write_nexus_trees = msprime_extras.write_nexus_trees

#R tree metrics assume tips are numbered from 1 not 0
tree_tip_labels_start_at_0 = False

if sys.version_info[0] < 3:
    raise Exception("Python 3 only")


def cpu_time_colname(tool):
    return tool + "_cputime"

def memory_colname(tool):
    return tool + "_memory"

def metric_colnames(metrics_for):
    metric_names = list(ARG_metrics.get_ARG_metrics())
    return ["{}_{}".format(f, metric) for f in metrics_for for metric in metric_names]

def nanblank(val):
    """hack around a horrible pandas syntax, which puts nan instead of blank strings"""
    return "" if pd.isnull(val) else val

def always_true(*pargs):
    """
    A func that returns True for any input value
    """
    return True
    
def make_errors(v, p):
    """
    Flip each bit with probability p.
    """
    m = v.shape[0]
    mask = np.random.random(m) < p
    return np.logical_xor(v.astype(bool), mask).astype(int)

def generate_samples(ts, error_p):
    """
    Returns samples with a bits flipped with a specified probability.

    Rejects any variants that result in a fixed column.
    """
    S = np.zeros((ts.sample_size, ts.num_mutations), dtype="u1")
    for variant in ts.variants():
        done = False
        # Reject any columns that have no 1s or no zeros
        while not done:
            S[:,variant.index] = make_errors(variant.genotypes, error_p)
            s = np.sum(S[:, variant.index])
            done = 0 < s < ts.sample_size
    return S

def msprime_name(n, Ne, l, rho, mu, genealogy_seed, mut_seed, directory=None):
    """
    Create a filename for saving an msprime simulation (without extension)
    Other functions add error rates and/or a subsample sizes to the name
    """
    #format mut rate & recomb rate to print non-exponential notation without
    # trailing zeroes 12 dp should be ample for these rates
    rho = "{:.12f}".format(float(rho)).rstrip('0')
    mu = "{:.12f}".format(float(mu)).rstrip('0')
    file = "msprime-n{}_Ne{}_l{}_rho{}_mu{}-gs{}_ms{}".format(int(n), float(Ne), int(l), rho, \
        mu, int(genealogy_seed), int(mut_seed))
    if directory is None:
        return file
    else:
        return os.path.join(directory,file)

def msprime_name_from_row(row, directory=None, error_col=None, subsample_col=None):
    """
    If error_col and subsample_col are None, this is the same as msprime_name() 
    but filled out using data from a row. If error_col is a string which exists as
    a column name in the row then add_error_param_to_name() is also called, using 
    the error rate specified in that column. Alternatively (e.g. if error_col is a 
    number) then it is treated as the error rate to add via add_error_param_to_name(). 
    The same goes for subsample_col.
    """
    name = msprime_name(row.sample_size, row.Ne, row.length, row.recombination_rate,
        row.mutation_rate, row.seed, row.seed, directory)
    if subsample_col is not None:
        if isinstance(subsample_col, str):
            if subsample_col in row:
                name = add_subsample_param_to_name(name, row[subsample_col])
        else:
            name = add_subsample_param_to_name(name, subsample_col)
    if error_col is not None:
        if isinstance(error_col, str):
            if error_col in row:
                name = add_error_param_to_name(name, row[error_col])
        else:
            name = add_error_param_to_name(name, error_col)
    return(name)

def add_subsample_param_to_name(sim_name, subsample_size):
    """
    Mark a filename as containing only a subset of the samples of the full sim
    Can be used on msprime output files but also e.g. tsinfer output files
    """
    if sim_name.endswith("+") or sim_name.endswith("-"):
        #this is the first param
        return sim_name + "max{}".format(int(subsample_size))
    else:
        return sim_name + "_max{}".format(int(subsample_size))

def add_error_param_to_name(sim_name, error_rate):
    """
    Append the error param to the msprime simulation filename.
    Only relevant for files downstream of the step where sequence error is added
    """
    if sim_name.endswith("+") or sim_name.endswith("-"):
        #this is the first param
        return sim_name + "_err{}".format(float(error_rate))
    else:
        #this is the first param
        return sim_name + "err{}".format(float(error_rate))


def construct_fastarg_name(sim_name, seed, directory=None):
    """
    Returns a fastARG inference filename (without file extension),
    based on a simulation name
    """
    d,f = os.path.split(sim_name)
    return os.path.join(d,'+'.join(['fastarg', f, "fs"+str(int(seed))]))

def construct_argweaver_name(sim_name, seed, iteration_number=None):
    """
    Returns an ARGweaver inference filename (without file extension),
    based on a simulation name. The iteration number (used in .smc and .nex output)
    is usually added by the ARGweaver `arg-sample` program,
    in the format .10, .20, etc. (we append an 'i' too, giving
    'i.10', 'i.100', etc
    """
    d,f = os.path.split(sim_name)
    suffix = "ws"+str(int(seed))
    if iteration_number is not None:
        suffix += "_i."+ str(int(iteration_number))
    return os.path.join(d,'+'.join(['aweaver', f, suffix]))

def construct_tsinfer_name(sim_name, subsample_size=None):
    """
    Returns an MSprime Li & Stevens inference filename.
    If the file is a subset of the original, this can be added to the
    basename in this function, or later using the 
    add_subsample_param_to_name() routine.
    """
    d,f = os.path.split(sim_name)
    name = os.path.join(d,'+'.join(['tsinfer', f, ""]))
    if subsample_size is not None:
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


class InferenceRunner(object):
    """
    Class responsible for running a single inference tool and returning results for
    the dataframe. Should create results files that are bespoke for each tool, and
    also for each tool, convert these into nexus files that can be used for comparing
    metrics.
    """
    def __init__(self, tool, row, simulations_dir, num_threads):
        self.tool = tool
        self.row = row
        self.num_threads = num_threads
        self.base_fn = msprime_name_from_row(row, simulations_dir, 
            'error_rate', 'subsample')

    def run(self):
        logging.info("Row {}: running {} inference".format(int(self.row[0]),self.tool))
        logging.debug("parameters = {}".format(self.row.to_dict()))
        if self.tool == "tsinfer":
            ret = self.__run_tsinfer()
        elif self.tool == "fastARG":
            ret = self.__run_fastARG()
        elif self.tool == "ARGweaver":
            ret = self.__run_ARGweaver()
        else:
            raise KeyError("unknown tool {}".format(self.tool))
        logging.debug("returning infer results for {} row {} = {}".format(
            self.tool, int(self.row[0]), ret))
        return ret


    def __run_tsinfer(self):

        samples_fn = self.base_fn + ".npy"
        logging.debug("reading: variant matrix {} for msprime inference".format(samples_fn))
        positions_fn = self.base_fn + ".pos.npy"
        logging.debug("reading: positions {} for msprime inference".format(positions_fn))
        scaled_recombination_rate = 4 * self.row.recombination_rate * self.row.Ne
        inferred_ts, time, memory = self.run_tsinfer(
            samples_fn, positions_fn, self.row.length, scaled_recombination_rate,
            num_threads=self.num_threads)
        if 'tsinfer_subset' in self.row:
            inferred_ts = inferred_ts.simplify(list(range(self.row.tsinfer_subset)))
            out_fn = construct_tsinfer_name(self.base_fn, self.row.tsinfer_subset)
        else:
            out_fn = construct_tsinfer_name(self.base_fn)
            
        with open(out_fn +".nex", "w+") as out:
            #tree metrics assume tips are numbered from 1 not 0
            inferred_ts.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
        return  {
            cpu_time_colname(self.tool): time,
            memory_colname(self.tool): memory
        }

    def __run_fastARG(self):
        inference_seed = self.row.seed  # TODO do we need to specify this separately?
        infile = self.base_fn + ".hap"
        out_fn = construct_fastarg_name(self.base_fn, inference_seed)
        logging.debug("reading: {} for fastARG inference".format(infile))
        inferred_ts, time, memory = self.run_fastarg(infile, self.row.length, inference_seed)
        with open(out_fn +".nex", "w+") as out:
            inferred_ts.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
        return {
            cpu_time_colname(self.tool): time,
            memory_colname(self.tool): memory
        }

    def __run_ARGweaver(self):
        inference_seed = self.row.seed  # TODO do we need to specify this separately?
        infile = self.base_fn + ".sites"
        out_fn = construct_argweaver_name(self.base_fn, inference_seed)
        logging.debug("reading: {} for ARGweaver inference".format(infile))
        iteration_ids, stats_file, time, memory = self.run_argweaver(
            infile, self.row.Ne, self.row.recombination_rate, self.row.mutation_rate,
            out_fn, inference_seed, self.row.aw_n_out_samples,
            self.row.aw_iter_out_freq, self.row.aw_burnin_iters)
        #now must convert all of the .smc files to .nex format
        for it in iteration_ids:
            base = construct_argweaver_name(self.base_fn, inference_seed, it)
            with open(base + ".nex", "w+") as out:
                msprime_ARGweaver.ARGweaver_smc_to_nexus(
                    base+".smc.gz", out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
        results = {
            cpu_time_colname(self.tool): time,
            memory_colname(self.tool): memory,
            "ARGweaver_iterations": ",".join(iteration_ids)
        }
        return results

    @staticmethod
    def run_tsinfer(sample_fn, positions_fn, length, rho, num_threads=1):
        with tempfile.NamedTemporaryFile("w+") as ts_out:
            cmd = [
                sys.executable, tsinfer_executable, sample_fn, positions_fn,
                "--length", str(int(length)), "--recombination-rate", str(rho),
                "--threads", str(num_threads), ts_out.name]
            cpu_time, memory_use = time_cmd(cmd)
            ts_simplified = msprime.load(ts_out.name)
        return ts_simplified, cpu_time, memory_use

    @staticmethod
    def run_fastarg(file_name, seq_length, seed):
        with tempfile.NamedTemporaryFile("w+") as fa_out, \
                tempfile.NamedTemporaryFile("w+") as tree, \
                tempfile.NamedTemporaryFile("w+") as muts, \
                tempfile.NamedTemporaryFile("w+") as fa_revised:
            cmd = msprime_fastARG.get_cmd(fastARG_executable, file_name, seed)
            cpu_time, memory_use = time_cmd(cmd, fa_out)
            logging.debug("ran fastarg [{} s]: '{}'".format(cpu_time, cmd))
            var_pos = msprime_fastARG.variant_positions_from_fastARGin_name(file_name)
            root_seq = msprime_fastARG.fastARG_root_seq(fa_out)
            msprime_fastARG.fastARG_out_to_msprime_txts(
                    fa_out, var_pos, tree, muts, seq_len=seq_length)
            inferred_ts = msprime_fastARG.msprime_txts_to_fastARG_in_revised(
                    tree, muts, root_seq, fa_revised, simplify=True)
            try:
                assert filecmp.cmp(file_name, fa_revised.name, shallow=False), \
                    "{} and {} differ".format(file_name, fa_revised.name)
            except AssertionError:
                warn("File '{}' copied to 'bad.hap' for debugging".format(
                    fa_revised.name), file=sys.stderr)
                shutil.copyfile(fa_revised.name, os.path.join('bad.hap'))
                raise
            return inferred_ts, cpu_time, memory_use

    @staticmethod
    def run_argweaver(
            sites_file, Ne, recombination_rate, mutation_rate, path_prefix, seed,
            n_samples, sampling_freq, burnin_iterations=0, quiet=True):
        """
        this produces a whole load of .smc files labelled <path_prefix>i.0.smc,
        <path_prefix>i.10.smc, etc. the iteration numbers ('i.0', 'i.10', etc)
        are returned by this function

        TO DO: if verbosity < 0 (logging level == warning) we should set quiet = TRUE

        """
        new_prefix = path_prefix + "_i" #we append a '_i' to mark iteration number
        before = time.clock()
        # TODO change this to get a command line for argweaver which we then run
        # using time_cmd.
        msprime_ARGweaver.run_ARGweaver(Ne=Ne,
            mut_rate     = mutation_rate,
            recomb_rate  = recombination_rate,
            executable   = ARGweaver_executable,
            ARGweaver_in = sites_file,
            sample_step  = sampling_freq,
            iterations   = sampling_freq * (n_samples-1),
            #e.g. for 3 samples at sampling_freq 10, run for 10*(3-1) iterations
            #which gives 3 samples: t=0, t=10, & t=20
            rand_seed    = int(seed),
            burn_in      = int(burnin_iterations),
            quiet        = quiet,
            out_prefix   = new_prefix)
        cpu_time = time.clock() - before
        memory_use = 0 #To DO
        #here we need to look for the .smc files output by argweaver
        smc_prefix = new_prefix + "." #the arg-sample program adds .iteration_num
        iterations = [f[len(smc_prefix):-7] for f in glob.glob(smc_prefix + "*" + ".smc.gz")]
        new_stats_file_name = path_prefix+".stats"
        os.rename(new_prefix + ".stats", new_stats_file_name)
        #cannot translate these to msprime ts objects, as smc2arg does not work
        #see https://github.com/mdrasmus/argweaver/issues/20
        return iterations, new_stats_file_name, cpu_time, memory_use



def infer_worker(work):
    """
    Entry point for running a single inference task in a worker process.
    """
    tool, row, simulations_dir, num_threads = work
    runner = InferenceRunner(tool, row, simulations_dir, num_threads)
    return int(row[0]), runner.run()

class MetricsRunner(object):
    """
    Class responsible for firing off the rpy2 code to calculate metrics that
    compare an original simulation file against the result from a set of tools
    (one class instantiated per row).
    Results are returned that can be incorporated into the main dataframe.
    """
    def __init__(self, row, nexus_dir, num_threads):
        self.row = row
        self.nexus_dir = nexus_dir
        self.num_threads = num_threads
        self.tools=collections.OrderedDict()
        # the original file against which others should be compared is the one 
        # without errors injected, but could potentially be a subsetted one
        self.simulation_comparison_fn=msprime_name_from_row(row, self.nexus_dir, 
            error_col=None, subsample_col='base_subsample_size')

    def add_tool(self, toolname, filenames, reps=None, make_bin_seed=None):
        """
        Adds a tool for this simulation. Some simulations wil create
        multiple output files, results from which can be averaged.
        Other simulations may have a single file, but the metric calculation 
        may need to be run multiple times (e.g. if polytomies/multifurcations
        are resolved at random on each successive calculation of a metric.
        If so, the 'reps' parameter gives a number of replicates, and the 
        result returned is an average over different reps.
        
        Columns named `toolname`-`metric` should exist in the row
        """
        self.tools[toolname] = {'nexus':filenames, 'reps':reps, 'make_bin_seed': make_bin_seed}
    

    def run(self):
        if self.tools:
            return ARG_metrics.get_ARG_metrics(self.simulation_comparison_fn + ".nex", 
                threads=self.num_threads, **self.tools)
        else:
            #called with nothing to compare!
            return None

def metric_worker(work):
    """
    Entry point for running a single set of metric calculations in a worker process.
    This is called multiple times for each worker process.
    
    Each row of the dataset compares a single original file with the result from 
    multiple tools. To save having to read the same original file multiple times
    we do the calculation for all tools withing a single instance of metric_worker.
    """
    metrics, row, simulations_dir, num_threads = work
    runner = MetricsRunner(row, simulations_dir, num_threads)
    #logging.debug("Running metrics on row")
    for tool, nex_params in metrics.items():
        try:
            nexus, params = nex_params(row, simulations_dir)
        except ValueError:
            nexus, = nex_params(row, simulations_dir)
            params={}
        #now add the .nex extension
        if isinstance(nexus, str):
            nexus += ".nex"
            if not os.path.isfile(nexus):
                nexus = None
        else:
            nexus = [fn + ".nex" for fn in nexus if os.path.isfile(fn + ".nex")]
        if nexus:
            runner.add_tool(tool, nexus, **params)
        else:
            logging.debug("Skipping metrics for " + tool +
                " (row " + str(row[0]) + ") " +
                "because of missing nexus file(s)")
    return int(row[0]), runner.run()

class Dataset(object):
    """
    A dataset is some collection of simulations and associated data.
    """
    name = None
    """
    Each dataset has a unique name. This is used as the prefix for the data
    file and raw_data_dir directory. Within this, replicate instances of datasets
    (each with a different RNG seed) are saved under the seed number
    """

    data_dir = "data"

    #the tools dict contains functions that are called on each row
    #that define when to use each tool. This allows us to run only
    #some tools for some of the simulations. The defaults here can 
    #be overridden in each dataset class, to run only a subset
    tools = {
        "ARGweaver": always_true,
        "fastARG"  : always_true,
        "tsinfer"  : always_true}
    #the metrics_for dict defines what metrics we calculate
    #and how to calculate them, given a row in a df
    #it is assumed that each will be a function that
    # is called with a row of data and a simulation dir 
    #The function should output a tuple whose first value 
    #is the interence file name, and whose (optional) second 
    #value give further parameters to the metrics function
    metrics_for = {
        "fastARG": lambda row, sim_dir: (
            construct_fastarg_name(msprime_name_from_row(row, sim_dir, 'error_rate', 'subsample'),
                                   seed=row.seed),
            ),
        "Aweaver": lambda row, sim_dir: ([
            construct_argweaver_name(msprime_name_from_row(row, sim_dir, 'error_rate', 'subsample'), 
                                     seed=row.seed, iteration_number=it)
                for it in nanblank(row.ARGweaver_iterations).split(",") if it],
            ),
        "tsipoly": lambda row, sim_dir: (
            construct_tsinfer_name(msprime_name_from_row(row, sim_dir, 'error_rate', 'subsample')),
            ),
        "tsibiny": lambda row, sim_dir: (
            construct_tsinfer_name(msprime_name_from_row(row, sim_dir, 'error_rate', 'subsample')),
            #hack same polytomy seed as inference seed
            {'make_bin_seed':row.seed, 'reps':row.tsinfer_biforce_reps} 
            )
    }


    def __init__(self):
        self.data_path = os.path.abspath(
            os.path.join(self.data_dir, "{}".format(self.name)))
        self.data_file = self.data_path + "_data.csv"
        self.param_file = self.data_path + "_setup.json"
        self.raw_data_dir = os.path.join(self.data_dir, "raw__NOBACKUP__", self.name)
        self.simulations_dir = os.path.join(self.raw_data_dir, "simulations")

    def load_data(self):
        self.data = pd.read_csv(self.data_file)

    def dump_data(self, write_index=False):
        self.data.to_csv(self.data_file, index=write_index)

    def dump_setup(self, arg_dict):
        with open(self.param_file , "w+") as setup:
            json.dump(arg_dict, setup, sort_keys=True, indent=2)

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
        logging.info("Creating dir {}".format(self.simulations_dir))
        self.data = self.run_simulations(args.replicates, args.seed)
        self.verbosity = args.verbosity
        # Add the result columns
        extra_cols = collections.OrderedDict()
        for tool in sorted(self.tools.keys()):
            extra_cols[cpu_time_colname(tool)]=np.NaN
            extra_cols[memory_colname(tool)]=np.NaN
        # We need to store more information in the case of ARGweaver, since
        # each ARGweaver run produces a whole set of iterations. We join these
        # together in a single column
        if "ARGweaver" in self.tools:
            extra_cols["ARGweaver_iterations"]=""
        
        #add the columns for the ARG metrics
        extra_cols.update([(k,np.NaN) for k in metric_colnames(self.metrics_for.keys())])
        for col, default in extra_cols.items():
            self.data[col] = default
        self.dump_data(write_index=True)
        self.dump_setup({k:v for k,v in vars(args).items() if k != "func"})

    def infer(self, num_processes, num_threads, force=False):
        """
        Runs the main inference processes and stores results in the dataframe.
        """
        self.load_data()
        work = []
        for i in self.data.index:
            row = self.data.iloc[i]
            for tool in [tool for tool,func in self.tools.items() if func(self, row)]:
                # All values that are unset should be NaN, so we only run those that
                # haven't been filled in already. This allows us to stop and start the
                # infer processes without having to start from scratch each time.
                if force or pd.isnull(row[cpu_time_colname(tool)]):
                    work.append((tool, row, self.simulations_dir, num_threads))
                else:
                    logging.info("Data row {} is filled out for {} inference: skipping".format(i, tool))
        logging.info("running {} inference trials (max {} tools over {} of {} rows) with {} processes and {} threads".format(
            len(work), len(self.tools), int(np.ceil(len(work)/len(self.tools))), 
            len(self.data.index), num_processes, num_threads))
        if num_processes > 1:
            with multiprocessing.Pool(processes=num_processes) as pool:
                for row_id, updated in pool.imap_unordered(infer_worker, work):
                    for k, v in updated.items():
                        self.data.ix[row_id, k] = v
                    self.dump_data()
        else:
            # When we have only one process it's easier to keep everything in the same
            # process for debugging.
            for row_id, updated in map(infer_worker, work):
                for k, v in updated.items():
                    self.data.ix[row_id, k] = v
                self.dump_data()


    #
    # Utilities for running simulations and saving files.
    #
    def single_simulation(self, n, Ne, l, rho, mu, seed, mut_seed=None, 
        discretise_mutations=True):
        """
        The standard way to run one msprime simulation for a set of parameter
        values. Saves the output to an .hdf5 file, and also saves variant files
        for use in fastARG (a .hap file, in the format specified by
        https://github.com/lh3/fastARG#input-format) ARGweaver (a .sites file:
        http://mdrasmus.github.io/argweaver/doc/#sec-file-sites) tsinfer
        (currently a numpy array containing the variant matrix)

        mutation_seed is not yet implemented, but should allow the same
        ancestry to be simulated (if the same genealogy_seed is given) but have
        different mutations thrown onto the trees (even with different
        mutation_rates)

        Returns a tuple of treesequence, filename (without file type extension)
        """
        logging.info(
            "Running simulation for n = {}, l = {}, rho={}, mu = {} seed={}".format(
                n, l, rho, mu, seed))
        # Since we want to have a finite site model, we force the recombination map
        # to have exactly l loci with a recombination rate of rho between them.
        recombination_map = msprime.RecombinationMap.uniform_map(l, rho, l)
        # We need to rejection sample any instances that we can't discretise under
        # the current model. The simplest way to do this is to have a local RNG
        # which we seed with the specified seed.
        rng = random.Random(seed)
        # TODO replace this with a proper finite site mutation model in msprime.
        done = False
        while not done:
            sim_seed = rng.randint(1, 2**31)
            ts = msprime.simulate(
                n, Ne, recombination_map=recombination_map, mutation_rate=mu,
                random_seed=sim_seed)
            if discretise_mutations:
                try:
                    ts = msprime_extras.discretise_mutations(ts)
                    done = True
                except ValueError as ve:
                    logging.info("Rejecting simulation: seed={}: {}".format(sim_seed, ve))
            else:
                done=True
        # Here we might want to iterate over mutation rates for the same
        # genealogy, setting a different mut_seed so that we can see for
        # ourselves the effect of mutation rate variation on a single topology
        # but for the moment, we don't bother, and simply write
        # mut_seed==genealogy_seed
        sim_fn = msprime_name(n, Ne, l, rho, mu, seed, seed, self.simulations_dir)
        logging.debug("writing {}.hdf5".format(sim_fn))
        ts.dump(sim_fn+".hdf5", zlib_compression=True)
        return ts, sim_fn

    @staticmethod
    def save_variant_matrices(ts, fname, error_rate=None, integer_positions=True):
        if error_rate is None:
            S = generate_samples(ts, 0)
            filename = fname
        else:
            S = generate_samples(ts, error_rate)
            filename = add_error_param_to_name(fname, error_rate)
        outfile = filename + ".npy"
        logging.debug("writing variant matrix to {} for tsinfer".format(outfile))
        np.save(outfile, S)
        outfile = filename + ".pos.npy"
        pos = np.array([v.position for v in ts.variants()])
        logging.debug("writing variant positions to {} for tsinfer".format(outfile))
        np.save(outfile, pos)
        if integer_positions:
            assert all(p.is_integer() for p in pos), \
                "Variant positions are not all integers in {}".format()
        logging.debug("writing variant matrix to {}.hap for fastARG".format(filename))
        with open(filename+".hap", "w+") as fastarg_in:
            msprime_fastARG.variant_matrix_to_fastARG_in(np.transpose(S), pos, fastarg_in)
        logging.debug("writing variant matrix to {}.sites for ARGweaver".format(filename))
        with open(filename+".sites", "w+") as argweaver_in:
            if error_rate is None and integer_positions == False:
                msprime_ARGweaver.msprime_to_ARGweaver_in(ts, argweaver_in)
            else:
                msprime_ARGweaver.variant_matrix_to_ARGweaver_in(np.transpose(S), pos,
                    ts.get_sequence_length(), argweaver_in)

    def process(self, num_processes, num_threads, force=False):
        """
        Runs the main metric calculating processes and stores results in the dataframe.
        Should be able to cope with missing nexus files, e.g. if inference only run
        for a subset of tools
        """
        self.load_data()
        work = []
        metric_cols = metric_colnames(self.metrics_for.keys())
        for i in self.data.index:
            # Any row without the metrics columns unset should be NaN, so we only run those that
            # haven't been filled in already. This allows us to stop and start the
            # infer processes without having to start from scratch each time.
            if force or np.all(pd.isnull(self.data.ix[i, metric_cols])):
                work.append((self.metrics_for, self.data.iloc[i], self.simulations_dir, num_threads))
            else:
                logging.info("Data row {} has metrics (partially) filled: skipping".format(i))
        logging.info("running metrics on {} rows with {} processes and {} threads".format(
            len(work), num_processes, num_threads))
        if num_processes > 1:
            with multiprocessing.Pool(processes=num_processes) as pool:
                for row_id, dataframe in pool.imap_unordered(metric_worker, work):
                    try:
                        for rowname, row in dataframe.iterrows():
                            colnames = ["{}_{}".format(rowname,col) for col in row.index]
                            self.data.ix[row_id, colnames] = tuple(row)
                        self.dump_data()
                    except AttributeError:
                        logging.debug("No dataframe returned from metric calculation")
        else:
            # When we have only one process it's easier to keep everything in the same
            # process for debugging.
            for row_id, dataframe in map(metric_worker, work):
                try:
                    for rowname, row in dataframe.iterrows():
                        logging.debug("got {} for {}".format(tuple(row), rowname))
                        colnames = ["{}_{}".format(rowname,col) for col in row.index]
                        self.data.ix[row_id, colnames] = tuple(row)
                    self.dump_data()
                except AttributeError:
                    logging.debug("No dataframe returned from metric calculation")


class BasicTestDataset(Dataset):
    """
    This attempts to replicate the simulations and inferences previously carried out
    in test_treecmp.py - plotting fastarg against argweaver inference
    """
    name = "basic_test"
    tools = {
        "fastARG":always_true,
        "ARGweaver":always_true}
    default_replicates = 20
    default_seed = 123

    def __init__(self):
        super().__init__()
        #remove ts inferences
        for k in list(self.metrics_for.keys()):
            if k.startswith('ts'):
                del self.metrics_for[k]
        
    def setup(self, args):
        self.JeromesDiscretise = not args.hack_finite_sites
        super().setup(args)
    
    def run_simulations(self, replicates, seed):
        # TODO there is a fair amount of shared code here between this and the
        # MetricsByMutationRateDataset. Factor these out once a few more datasets
        # have been added and the common patterns are clear.
        if replicates is None:
            replicates = self.default_replicates
        if seed is None:
            seed = self.default_seed
        rng = random.Random(seed)
        seeds = set()
        cols = [
            "sample_size", "Ne", "length", "recombination_rate", "mutation_rate",
            "replicate", "seed", "aw_burnin_iters", "aw_n_out_samples", "aw_iter_out_freq"]
        sample_size= 8
        Ne = 1e4
        length=int(1e4)
        recombination_rate = 2e-8
        mutation_rates = [2e-8, 5e-8, 1e-7, 5e-7, 1e-6, 3e-6, 5e-6]
        aw_burnin_iters = 100
        aw_n_out_samples = 20
        aw_iter_out_freq = 10
        num_rows = replicates * len(mutation_rates)
        data = pd.DataFrame(index=np.arange(0, num_rows), columns=cols)
        row_id = 0
        for mutation_rate in mutation_rates:
            for replicate in range(replicates):
                done = False
                while not done:
                    replicate_seed = rng.randint(1, 2**31)
                    if replicate_seed not in seeds:
                        seeds.add(replicate_seed)
                        done = True
                # Run the simulation
                ts, fn = self.single_simulation(
                    sample_size, Ne, length, recombination_rate, mutation_rate,
                    replicate_seed, replicate_seed, discretise_mutations=self.JeromesDiscretise)
                with open(fn +".nex", "w+") as out:
                    ts.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
                # Add the rows for each of the error rates in this replicate
                row = data.iloc[row_id]
                row_id += 1
                row.sample_size = sample_size
                row.recombination_rate = recombination_rate
                row.mutation_rate = mutation_rate
                row.length = length
                row.Ne = Ne
                row.seed = replicate_seed
                row.replicate = replicate
                row.aw_n_out_samples = aw_n_out_samples
                row.aw_burnin_iters = aw_burnin_iters
                row.aw_iter_out_freq = aw_iter_out_freq
                self.save_variant_matrices(ts, fn, error_rate=None, 
                    integer_positions=self.JeromesDiscretise)
        return data

class NumRecordsBySampleSizeDataset(Dataset):
    """
    Information on the number of coalescence records inferred by tsinfer
    and FastARG for various sample sizes, under 3 different error rates
    """
    name = "num_records_by_sample_size"
    #the tools dict contains functions that are called on each row
    #that define when to use each tool
    tools = {
        "fastARG":always_true,
        "tsinfer":always_true}
    
    default_replicates = 10
    default_seed = 123
        
    def run_simulations(self, replicates, seed):
        # TODO there is a fair amount of shared code here between this and the
        # MetricsByMutationRateDataset. Factor these out once a few more datasets
        # have been added and the common patterns are clear.
        if replicates is None:
            replicates = self.default_replicates
        if seed is None:
            seed = self.default_seed
        rng = random.Random(seed)
        seeds = set()
        cols = [
            "sample_size", "Ne", "length", "recombination_rate", "mutation_rate",
            "error_rate", "replicate", "seed"]
        sample_sizes = np.linspace(10, 500, num=10).astype(int)
        error_rates = [0, 0.1, 0.01]
        recombination_rate = 2.5e-8
        mutation_rate = 1.5e-8
        length = 50000
        Ne = 5000
        num_rows = replicates * len(sample_sizes) * len(error_rates)
        data = pd.DataFrame(index=np.arange(0, num_rows), columns=cols)
        row_id = 0
        for sample_size in sample_sizes:
            for replicate in range(replicates):
                done = False
                while not done:
                    replicate_seed = rng.randint(1, 2**31)
                    if replicate_seed not in seeds:
                        seeds.add(replicate_seed)
                        done = True
                # Run the simulation
                ts, fn = self.single_simulation(
                    sample_size, Ne, length, recombination_rate, mutation_rate,
                    replicate_seed, replicate_seed)
                with open(fn +".nex", "w+") as out:
                    ts.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
                # Add the rows for each of the error rates in this replicate
                for error_rate in error_rates:
                    row = data.iloc[row_id]
                    row_id += 1
                    row.sample_size = sample_size
                    row.recombination_rate = recombination_rate
                    row.mutation_rate = mutation_rate
                    row.length = length
                    row.Ne = Ne
                    row.seed = replicate_seed
                    row.error_rate = error_rate
                    row.replicate = replicate
                    self.save_variant_matrices(ts, fn, error_rate)
        return data



class MetricsByMutationRateDataset(Dataset):
    """
    Dataset for Figure 1
    Accuracy of ARG inference (measured by various statistics)
    tending to fully accurate as mutation rate increases
    """
    name = "metrics_by_mutation_rate"
    
    default_replicates = 10
    default_seed = 123

    
    def run_simulations(self, replicates, seed):
        if replicates is None:
            replicates = self.default_replicates
        if seed is None:
            seed = self.default_seed
        rng = random.Random(seed)
        seeds = set()
        cols = [
            "sample_size", "Ne", "length", "recombination_rate", "mutation_rate",
            "error_rate", "replicate", "seed", "aw_burnin_iters",
            "aw_n_out_samples", "aw_iter_out_freq", "tsinfer_biforce_reps"]
        # Variable parameters
        mutation_rates = np.logspace(-8, -5, num=6)[:-2] * 1.5
        error_rates = [0, 0.01, 0.1]

        # Fixed parameters
        sample_size = 12
        Ne = 5000
        length = 50000
        recombination_rate = 2.5e-8
        ## argweaver params: aw_n_out_samples will be produced, every argweaver_iter_out_freq
        aw_burnin_iters = 5000
        aw_n_out_samples = 100
        aw_iter_out_freq = 10
        # TMP for development
        ## tsinfer params: number of times to randomly resolve into bifurcating (binary) trees
        tsinfer_biforce_reps = 20
        num_rows = replicates * len(mutation_rates) * len(error_rates)
        data = pd.DataFrame(index=np.arange(0, num_rows), columns=cols)
        row_id = 0
        for mutation_rate in mutation_rates:
            for replicate in range(replicates):
                done = False
                while not done:
                    replicate_seed = rng.randint(1, 2**31)
                    if replicate_seed not in seeds:
                        seeds.add(replicate_seed)
                        done = True
                # Run the simulation
                ts, fn = self.single_simulation(
                    sample_size, Ne, length, recombination_rate, mutation_rate,
                    replicate_seed, replicate_seed)
                with open(fn +".nex", "w+") as out:
                    ts.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
                # Add the rows for each of the error rates in this replicate
                for error_rate in error_rates:
                    row = data.iloc[row_id]
                    row_id += 1
                    row.sample_size = sample_size
                    row.recombination_rate = recombination_rate
                    row.mutation_rate = mutation_rate
                    row.length = length
                    row.Ne = Ne
                    row.seed = replicate_seed
                    row.error_rate = error_rate
                    row.replicate = replicate
                    row.aw_n_out_samples = aw_n_out_samples
                    row.aw_burnin_iters = aw_burnin_iters
                    row.aw_iter_out_freq = aw_iter_out_freq
                    row.tsinfer_biforce_reps = tsinfer_biforce_reps
                    self.save_variant_matrices(ts, fn, error_rate)
        return data

class SampleSizeEffectOnSubsetDataset(Dataset):
    """
    Dataset for Figure 3
    Dataset testing the effect on a subset of samples when extra genomes are
    added to the dataset used for inference. The hope is that the extra samples
    allowed by our tsinference method will more than compensate for biases in the
    inference method over more statistically rigorous methods (e.g. ARGweaver).

    We should expect this to be the case, as the benefit of adding more
    intermediate branches to trees (e.g. to halt long-branch attraction
    artifacts) is well documented in the phylogenetic literature, under
    the general heading of 'taxon sampling'.

    Method: take a large simulated dataset (thousands?), with fixed, realistic
    mutation and recombination rates, over a large span of genome. Use the msprime
    simplify() function to cut it down to a subset of n genomes (n ~ 10 or 20).
    Call this subset S_n_base. Run all the inference methods on this base subset
    to get inference outputs (eventually saved as .nex files). Then take larger
    and larger subsets of the original simulation - S_N for N=20, 40, 100, etc. -
    and run tsinfer (but not the other methods) on these larger subsets, ensuring
    that after each ts inference attempt, the resulting TreeSequence is subset
    (simplify()ed) back down to cover only the first n samples. Use ARGmetrics to
    compare these secondarily-simplified ARGs with the true (original, simulated)
    S_n subset. We would hope that the ARG accuracy metrics tend to 0 as the tsinfer
    subset size N goes up.
    """
    name = "sample_size_effect_on_subset"
    base_subsample_size = 10
    tools = {
        "tsinfer":   always_true,
        #only use fastARG & ARGweaver tools on small subsamples
        "fastARG":   lambda self, row: row.subsample==self.base_subsample_size, 
        "ARGweaver": lambda self, row: row.subsample==self.base_subsample_size
        }
    default_replicates = 10
    default_seed = 123

    def __init__(self):
        super().__init__()
        #override the built-in tsinference metrics, so that we use the subsetted files
        self.metrics_for["tsipoly"] = lambda row, sim_dir: (
            construct_tsinfer_name(msprime_name_from_row(row, sim_dir, 'error_rate', 'subsample'),
                                   self.base_subsample_size),
            )
        self.metrics_for["tsibiny"] = lambda row, sim_dir: (
            construct_tsinfer_name(msprime_name_from_row(row, sim_dir, 'error_rate', 'subsample'),
                                   self.base_subsample_size),
            #hack same polytomy seed as inference seed
            {'make_bin_seed':row.seed, 'reps':row.tsinfer_biforce_reps} 
            )
        
    def run_simulations(self, replicates, seed):
        if replicates is None:
            replicates = self.default_replicates
        if seed is None:
            seed = self.default_seed
        rng = random.Random(seed)
        seeds = set()
        cols = [
            "sample_size", "subsample", "Ne", "length", "recombination_rate", "mutation_rate",
            "error_rate", "replicate", "seed", "base_subsample_size", "aw_burnin_iters",
            "aw_n_out_samples", "aw_iter_out_freq", "tsinfer_biforce_reps"]
        # Variable parameters
        error_rates = [0, 0.01, 0.1]
        subsamples  = [self.base_subsample_size, 20, 50, 100, 200, 500, 1000]
        # Fixed parameters
        mutation_rate = 1.5e-8
        sample_size = 1000
        Ne = 5000
        length = 50000
        recombination_rate = 2.5e-8
        ## argweaver params: aw_n_out_samples will be produced, every argweaver_iter_out_freq
        aw_burnin_iters = 5000
        aw_n_out_samples = 100
        aw_iter_out_freq = 10
        # TMP for development
        ## tsinfer params: number of times to randomly resolve into bifurcating (binary) trees
        tsinfer_biforce_reps = 20
        num_rows = replicates * len(error_rates) * len(subsamples)
        data = pd.DataFrame(index=np.arange(0, num_rows), columns=cols)
        row_id = 0
        for replicate in range(replicates):
            done = False
            while not done:
                replicate_seed = rng.randint(1, 2**31)
                if replicate_seed not in seeds:
                    seeds.add(replicate_seed)
                    done = True
            # Run the simulation
            ts, fn = self.single_simulation(
                sample_size, Ne, length, recombination_rate, mutation_rate,
                replicate_seed, replicate_seed)
            with open(fn +".nex", "w+") as out:
                ts.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
            # Add the rows for each of the error rates in this replicate
            for subsample in subsamples:
                subfn = add_subsample_param_to_name(fn, subsample)
                ts_sub = ts.simplify(list(range(subsample)))
                with open(subfn +".nex", "w+") as out:
                    ts_sub.write_nexus_trees(out, zero_based_tip_numbers=tree_tip_labels_start_at_0)
                for error_rate in error_rates:
                    row = data.iloc[row_id]
                    row_id += 1
                    row.sample_size = sample_size
                    row.subsample = subsample
                    row.recombination_rate = recombination_rate
                    row.mutation_rate = mutation_rate
                    row.length = length
                    row.Ne = Ne
                    row.seed = replicate_seed
                    row.error_rate = error_rate
                    row.replicate = replicate
                    row.tsinfer_subset = self.base_subsample_size
                    row.aw_n_out_samples = aw_n_out_samples
                    row.aw_burnin_iters = aw_burnin_iters
                    row.aw_iter_out_freq = aw_iter_out_freq
                    row.tsinfer_biforce_reps = tsinfer_biforce_reps
                    self.save_variant_matrices(ts_sub, subfn, error_rate)
        return data

class Figure(object):
    """
    Superclass of all figures. Each figure depends on a dataset.
    """
    dataset = None
    name = None
    """
    Each figure has a unique name. This is used as the identifier and the
    file name for the output plots.
    """
    pdf_width_inches = 10
    pdf_height_inches = 7

    def __init__(self):
        self.filepath = self.dataset.data_path + "+" + self.name

    def R_plot(self, Rcmds):
        """
        Take the R commands used to generate a plot (as a single string
        or an array of lines) and evaluate them in R, saving the result 
        to a pdf file.
        """
        if isinstance(Rcmds, str):
            Rcmds = [Rcmds]
        Rcmds.insert(0, "if (!interactive()) pdf('{}',width={}, height={})".format(
            self.filepath + ".pdf", self.pdf_width_inches, self.pdf_height_inches))
        Rcmds.append("if (!interactive()) dev.off()")
        script = self.filepath + ".R"
        with open(script, "w+") as source:
            for line in Rcmds:
                print(line, file=source)
        subprocess.call(['R', 'CMD', 'BATCH', '--no-save', '--no-restore', script])
        logging.info("Plot file saved to {}. Code for generating plots interactively at {}".format(
            self.filepath + ".pdf", script))

    def R_plot_data(self, Rcmds, Rdata_cmd = None):
        """
        Rdata_cmd is an R command to create a 'data' object in the R
        session - if None, it defaults to reading the dataset csv file
        """
        if isinstance(Rcmds, str):
            Rcmds = [Rcmds]
        if Rdata_cmd is None:
            Rdata_cmd = "data <- read.csv('%s')".format(self.dataset.data_file)
        self.R_plot([Rdata_cmd] + Rcmds)
        
    def plot(self):
        raise NotImplementedError()

class BasicARGweaverVSfastARG(Figure):
    dataset = BasicTestDataset()
    name = "aw_vs_fa"
    
    def plot(self):
        self.R_plot_data(\
"""
datamean <- aggregate(subset(data, select=-ARGweaver_iterations), list(data$mutation_rate), mean)
tools <- c('fastARG', 'Aweaver')
metrics <- c('RFrooted', 'RFunrooted', 'wRFrooted', 'wRFunrooted', 'SPRunrooted', 'pathunrooted')
cols <- c(fastARG='red', Aweaver='green')
layout(matrix(1:6,2,3))
sapply(metrics, function(m) {
    colnames = paste(tools, m, sep='_')
    matplot(data$mutation_rate, data[, colnames], type='p', pch=c(1,2), col=cols[tools], main=m, log='x')
    matlines(datamean$mutation_rate, datamean[, colnames], type='l', lty=1, col=cols[tools])
    mtext(names(cols), line=c(-1.2,-2), adj=1, cex=0.7, col=cols)
})
"""
            )
                
def run_setup(cls, args):
    f = cls()
    f.setup(args)

def run_infer(cls, args):
    logging.info("Inferring {}".format(cls.name))
    f = cls()
    f.infer(args.processes, args.threads, args.force)


def run_process(cls, args):
    logging.info("Processing {}".format(cls.name))
    f = cls()
    f.process(args.processes, args.threads, args.force)


def run_plot(cls, args):
    f = cls()
    f.plot()


def main():
    datasets = [
        BasicTestDataset,
        NumRecordsBySampleSizeDataset,
        MetricsByMutationRateDataset,
        SampleSizeEffectOnSubsetDataset,
    ]
    figures = [
        BasicARGweaverVSfastARG
    ]
    name_map = dict([(d.name, d) for d in datasets + figures])
    parser = argparse.ArgumentParser(
        description="Set up base data, generate inferred datasets, process datasets and plot figures.")
    parser.add_argument('--verbosity', '-v', action='count', default=0)
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    subparser = subparsers.add_parser('setup')
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
    subparser.add_argument(
         '--replicates', '-r', type=int, help="number of replicates")
    subparser.add_argument(
         '--seed', '-s', type=int, help="use a non-default RNG seed")
    subparser.add_argument(
         '--hack_finite_sites', action='store_true',
         help="Mutations at the same (integer) location are superimposed, not shifted along")
    subparser.set_defaults(func=run_setup)

    subparser = subparsers.add_parser('infer')
    subparser.add_argument(
        "--processes", '-p', type=int, default=1,
        help="number of worker processes")
    subparser.add_argument(
        "--threads", '-t', type=int, default=1,
        help="number of threads per worker process (for supporting tools)")
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
    subparser.add_argument(
         '--data-file', '-f', type=str,
         help="which CSV file to use for existing data, if not the default", )
    subparser.add_argument(
         '--force',  action='store_true', 
         help="redo all the inferences, even if we have already filled out some", )
    subparser.set_defaults(func=run_infer)

    subparser = subparsers.add_parser('process')
    subparser.add_argument(
        "--processes", '-p', type=int, default=1,
        help="number of worker processes")
    subparser.add_argument(
        "--threads", '-t', type=int, default=1,
        help="number of threads per worker process (for supporting tools)")
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
    subparser.add_argument(
         '--data-file', '-f', type=str,
         help="which CSV file to use for existing data, if not the default")
    subparser.add_argument(
         '--force',  action='store_true', 
         help="redo all the metrics, even if we have already filled out some", )
    subparser.set_defaults(func=run_process)

    subparser = subparsers.add_parser('figure')
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the figure identifier')
    subparser.set_defaults(func=run_plot)

    args = parser.parse_args()
    log_level = logging.WARNING
    if args.verbosity == 1:
        log_level = logging.INFO
    if args.verbosity >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(
        format='%(asctime)s %(message)s', level=log_level, stream=sys.stdout)

    k = args.name[0]
    if k == "all":
        classes = datasets
        if args.func == run_plot:
            classes = figures
        for name, cls in name_map.items():
            if cls in classes:
                args.func(cls, args)
    else:
        cls = name_map[k]
        args.func(cls, args)

if __name__ == "__main__":
    main()
