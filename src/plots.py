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
#monkey-patch write.nexus into msprime
msprime.TreeSequence.write_nexus_trees = msprime_extras.write_nexus_trees

if sys.version_info[0] < 3:
    raise Exception("Python 3 only")


def cpu_time_colname(tool):
    return tool + "_cpu_time"

def memory_colname(tool):
    return tool + "_memory"

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

def add_subsample_param_to_name(sim_name, subsample_size):
    """
    Mark a filename as containing only a subset of the samples of the full sim
    Can be used on msprime output files but also e.g. tsinfer output files
    """
    if sim_name.endswith("+") or sim_name.endswith("-"):
        #this is the first param
        return sim_basename + "max{}".format(int(subsample_size))
    else:
        return sim_basename + "_max{}".format(int(subsample_size))

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

def construct_tsinfer_name(sim_name):
    """
    Returns an MSprime Li & Stevens inference filename.
    If the file is a subset of the original, this can be added to the
    basename later using the add_subsample_param_to_name() routine.
    """
    d,f = os.path.split(sim_name)
    return os.path.join(d,'+'.join(['tsinfer', f, ""]))

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


class Figure(object):
    """
    Superclass of all figures. Each figure depends on a dataset.
    """
    name = None
    """
    Each figure has a unique name. This is used as the identifier and the
    file name for the output plots.
    """

    def plot(self):
        raise NotImplementedError()


class InferenceRunner(object):
    """
    Class responsible for running a single inference tool and returning results for
    the dataframe.
    """
    def __init__(self, tool, row, simulations_dir, num_threads):
        self.tool = tool
        self.row = row
        self.num_threads = num_threads
        self.sim_fn = msprime_name(
            row.sample_size, row.Ne, row.length, row.recombination_rate,
            row.mutation_rate, row.seed, row.seed, simulations_dir)
        self.err_fn = add_error_param_to_name(self.sim_fn, row.error_rate)

    def run(self):
        logging.info("running {} inference on row {}".format(self.tool, int(self.row[0])))
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

        samples_fn = self.err_fn + ".npy"
        logging.debug("reading: variant matrix {} for msprime inference".format(samples_fn))
        positions_fn = self.err_fn + ".pos.npy"
        logging.debug("reading: positions {} for msprime inference".format(positions_fn))
        out_fn = construct_tsinfer_name(self.err_fn)
        scaled_recombination_rate = 4 * self.row.recombination_rate * self.row.Ne
        inferred_ts, time, memory = self.run_tsinfer(
            samples_fn, positions_fn, self.row.length, scaled_recombination_rate,
            num_threads=self.num_threads)
        with open(out_fn +".nex", "w+") as out:
            inferred_ts.write_nexus_trees(out)
        return  {
            cpu_time_colname(self.tool): time,
            memory_colname(self.tool): memory
        }

    def __run_fastARG(self):
        inference_seed = self.row.seed  # TODO do we need to specify this separately?
        infile = self.err_fn + ".hap"
        out_fn = construct_fastarg_name(self.err_fn, inference_seed)
        logging.debug("reading: {} for fastARG inference".format(infile))
        inferred_ts, time, memory = self.run_fastarg(infile, self.row.length, inference_seed)
        with open(out_fn +".nex", "w+") as out:
            inferred_ts.write_nexus_trees(out)
        return {
            cpu_time_colname(self.tool): time,
            memory_colname(self.tool): memory
        }

    def __run_ARGweaver(self):
        inference_seed = self.row.seed  # TODO do we need to specify this separately?
        infile = self.err_fn + ".sites"
        out_fn = construct_argweaver_name(self.err_fn, inference_seed)
        logging.debug("reading: {} for ARGweaver inference".format(infile))
        iteration_ids, stats_file, time, memory = self.run_argweaver(
            infile, self.row.Ne, self.row.recombination_rate, self.row.mutation_rate,
            out_fn, inference_seed, self.row.aw_n_out_samples,
            self.row.aw_iter_out_freq, self.row.aw_burnin_iters)
        #now must convert all of the .smc files to .nex format
        for it in iteration_ids:
            base = construct_argweaver_name(self.err_fn, inference_seed, it)
            with open(base + ".nex", "w+") as out:
                msprime_ARGweaver.ARGweaver_smc_to_nexus(
                    base+".smc.gz", out, zero_based_tip_numbers=False)
        results = {
            cpu_time_colname(self.tool): time,
            memory_colname(self.tool): memory,
            "ARGWeaver_iterations": ",".join(iteration_ids)
        }
        return results

    @staticmethod
    def run_tsinfer(sample_fn, positions_fn, length, rho, num_threads=1):
        with tempfile.NamedTemporaryFile("w+") as ts_out:
            cmd = [
                sys.executable, "src/run_tsinfer.py", sample_fn, positions_fn,
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
            n_samples, sampling_freq, burnin_iterations=0):
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
            quiet        = True,
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

    def __init__(self):
        self.data_file = os.path.abspath(
            os.path.join(self.data_dir, "{}.csv".format(self.name)))
        self.raw_data_dir = os.path.join(self.data_dir, "raw__NOBACKUP__", self.name)
        self.simulations_dir = os.path.join(self.raw_data_dir, "simulations")

    def load_data(self):
        self.data = pd.read_csv(self.data_file)

    def dump_data(self, write_index=False):
        self.data.to_csv(self.data_file, index=write_index)

    #
    # Main entry points.
    #
    def setup(self, replicates=None, seed=None):
        """
        Creates the dataframe and storage directories and then runs the initial
        simulations.
        """
        if os.path.exists(self.simulations_dir):
            shutil.rmtree(self.simulations_dir)
            logging.info("Deleting dir {}".format(self.simulations_dir))
        os.makedirs(self.simulations_dir)
        logging.info("Creating dir {}".format(self.simulations_dir))
        self.data = self.run_simulations(replicates, seed)
        # Add the result columns
        tool_cols = []
        for tool in self.tools:
            tool_cols.append(cpu_time_colname(tool))
            tool_cols.append(memory_colname(tool))
        # We need to store more information in the case of ARGweaver, since
        # each ARGweaver run produces a whole set of iterations. We join these
        # together in a single column
        if "ARGweaver" in self.tools:
            tool_cols.append("ARGweaver_iterations")
        #add the columns for the ARG metrics
        metric_colnames = list(ARG_metrics.get_ARG_metrics())
        for col in tool_cols:
            for metric_name in metric_colnames:
                self.data[col] = np.NaN
        self.dump_data(write_index=True)

    def infer(self, num_processes, num_threads):
        """
        Runs the main inference processes and stores results in the dataframe.
        """
        logging.info("running infer with {} processes and {} threads".format(
            num_processes, num_threads))
        self.load_data()
        work = []
        for i in self.data.index:
            for tool in self.tools:
                # All values that are unset should be NaN, so we only run those that
                # haven't been filled in already. This allows us to stop and start the
                # infer processes without having to start from scratch each time.
                if np.isnan(self.data.ix[i, cpu_time_colname(tool)]):
                    work.append((tool, self.data.iloc[i], self.simulations_dir, num_threads))
                else:
                    logging.info("Skipping row {} for {}".format(i, tool))
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
    def single_simulation(self, n, Ne, l, rho, mu, seed, mut_seed=None, subsample=None):
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
            try:
                ts = msprime_extras.discretise_mutations(ts)
                done = True
            except ValueError as ve:
                logging.info("Rejecting simulation: seed={}: {}".format(sim_seed, ve))
        # Here we might want to iterate over mutation rates for the same
        # genealogy, setting a different mut_seed so that we can see for
        # ourselves the effect of mutation rate variation on a single topology
        # but for the moment, we don't bother, and simply write
        # mut_seed==genealogy_seed
        sim_fn = msprime_name(n, Ne, l, rho, mu, seed, seed, self.simulations_dir)
        logging.debug("writing {}.hdf5".format(sim_fn))
        ts.dump(sim_fn+".hdf5", zlib_compression=True)
        if subsample is not None:
            ts = ts.simplify(list(range(subsample)))
            sim_fn = add_subsample_param_to_name(sim_fn, subsample)
            ts.dump(sim_fn +".hdf5", zlib_compression=True)
        return ts, sim_fn

    @staticmethod
    def save_variant_matrices(ts, fname, error_rate):
        S = generate_samples(ts, error_rate)
        err_filename = add_error_param_to_name(fname, error_rate)
        outfile = err_filename + ".npy"
        logging.debug("writing variant matrix to {} for tsinfer".format(outfile))
        np.save(outfile, S)
        outfile = err_filename + ".pos.npy"
        pos = np.array([v.position for v in ts.variants()])
        logging.debug("writing variant positions to {} for tsinfer".format(outfile))
        np.save(outfile, pos)
        assert all(p.is_integer() for p in pos), \
            "Variant positions are not all integers in {}".format()
        logging.debug("writing variant matrix to {}.hap for fastARG".format(err_filename))
        with open(err_filename+".hap", "w+") as fastarg_in:
            msprime_fastARG.variant_matrix_to_fastARG_in(np.transpose(S), pos, fastarg_in)
        logging.debug("writing variant matrix to {}.sites for ARGweaver".format(err_filename))
        with open(err_filename+".sites", "w+") as argweaver_in:
            msprime_ARGweaver.variant_matrix_to_ARGweaver_in(np.transpose(S), pos,
                ts.get_sequence_length(), argweaver_in)

    def process(self):
        """
        processes the inferred ARGs
        """
        raise NotImplementedError()

class NumRecordsBySampleSizeDataset(Dataset):
    """
    Information on the number of coalescence records inferred by tsinfer
    and FastARG for various sample sizes, under 3 different error rates
    """
    name = "num_records_by_sample_size"
    tools = ["fastARG", "tsinfer"]
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
            "error_rate", "replicate", "seed", "error_rate"]
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
                    ts.write_nexus_trees(out)
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
    tools = ["fastARG", "tsinfer", "ARGweaver"]
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
            "error_rate", "replicate", "seed", "error_rate", "aw_burnin_iters",
            "aw_n_out_samples", "aw_iter_out_freq", "tsinfer_biforce_reps"]
        # Variable parameters
        mutation_rates = np.logspace(-8, -5, num=10)[:-1] * 1.5
        error_rates = [0, 0.01, 0.1]

        # Fixed parameters
        sample_size = 12
        Ne = 5000
        length = 50000
        recombination_rate = 2.5e-8
        ## argweaver params: aw_n_out_samples will be produced, every argweaver_iter_out_freq
        # aw_burnin_iters = 5000
        # aw_n_out_samples = 100
        # TMP for development
        aw_burnin_iters = 5
        aw_n_out_samples = 10
        aw_iter_out_freq = 10
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
                    ts.write_nexus_trees(out)
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

    def process(self):
        """
        Extracts metrics from the .nex files using R, and saves the results
        in the csv file under fastARG_RFunrooted, etc etc.
        The tsinfer routine has two possible sets of stats: for nonbinary trees
        and for randomly resolved strictly bifurcating trees (with a random seed given)
        These are stored in 'msipoly_RFunrooted' and
        'msibifu_RFunrooted', and the random seed used (as passed to R via
        genome_trees_dist_forcebin_b) is stored in 'msibifu_seed'
        """
        # TODO update this to use a similar pattern the to the 'infer' method.
        # All columns should be added to the file at setup() time, and we can
        # run the argmetrics code in a subprocess.

        for i in self.data.index:
            d = self.data.iloc[i]
            sim_fn=msprime_name(d.sample_size, d.Ne, d.length, d.recombination_rate,
                                d.mutation_rate, d.seed, d.seed, self.simulations_dir)
            base_fn=add_error_param_to_name(sim_fn, d.error_rate)
            polytomy_resolution_seed = d.inference_seed #hack: use same seed as in inference step
            toolfiles = {
                "fastARG":construct_fastarg_name(base_fn, d.inference_seed)+".nex",
                "ARGweaver":{'nexus':construct_argweaver_name(base_fn, d.inference_seed, it)+".nex" \
                    for it in d['ARGweaver_iterations'].split(",")},
                "MSIpoly":construct_tsinfer_name(base_fn)+".nex",
                "MSIbifu":{'nexus':construct_tsinfer_name(base_fn)+".nex",
                           'reps':d.tsinfer_biforce_reps,
                           'seed':polytomy_resolution_seed}
            }
            logging.info("Processing ARG to extract metrics: mu = {}, err = {}.".format(
                d.mutation_rate, d.error_rate))
            m = ARG_metrics.get_ARG_metrics(sim_fn + ".nex", **toolfiles)
            self.data.loc[i,"MSIbifu_seed"] = polytomy_resolution_seed
            #Now add all the metrics to the data file
            for prefix, metrics in m.items():
                colnames = ["_".join([prefix,k]) for k in sorted(metrics.keys())]
                add_columns(self.data, colnames)
                self.data.loc[i, colnames] = metrics.values()

             # Save each row so we can use the information while it's being built
            self.data.to_csv(self.data_file)

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
    that after each inference attempt, the resulting TreeSequence is subset
    (simplify()ed) back down to cover only the first n samples. Use ARGmetrics to
    compare these secondarily-simplified ARGs with the true (original, simulated)
    S_n subset. We would hope that the ARG accuracy metrics tend to 0 as the tsinfer
    subset size N goes up.
    """
    name = "sample_size_effect_on_subset"


def run_setup(cls, args):
    f = cls()
    f.setup(args.replicates, args.seed)

def run_infer(cls, args):
    logging.info("Inferring {}".format(cls.name))
    f = cls()
    f.infer(args.processes, args.threads)


def run_process(cls, args):
    logging.info("Processing {}".format(cls.name))
    f = cls()
    f.process()


def run_plot(cls, args):
    f = cls()
    f.plot()


def main():
    datasets = [
        NumRecordsBySampleSizeDataset,
        MetricsByMutationRateDataset,
    ]
    figures = [
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
    subparser.set_defaults(func=run_infer)

    subparser = subparsers.add_parser('process')
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
    subparser.add_argument(
         '--data-file', '-f', type=str,
         help="which CSV file to use for existing data, if not the default")
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
