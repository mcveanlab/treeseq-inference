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
import collections
import itertools
import tempfile
import subprocess
import filecmp

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
import tsinf

fastARG_executable = os.path.join(curr_dir,'..','fastARG','fastARG')
ARGweaver_executable = os.path.join(curr_dir,'..','argweaver','bin','arg-sample')
#monkey-patch write.nexus into msprime
msprime.TreeSequence.write_nexus_trees = msprime_extras.write_nexus_trees

if sys.version_info[0] < 3:
    raise Exception("Python 3 only")

def expand_grid(data_dict):
    """from last example in http://pandas.pydata.org/pandas-docs/stable/cookbook.html"""
    rows = itertools.product(*data_dict.values())
    return pd.DataFrame.from_records(rows, columns=data_dict.keys())

def add_columns(dataframe, colnames):
    """
    add the named columns to a data frame unless they already exist
    """
    for cn in colnames:
        if cn not in dataframe.columns:
            dataframe[cn]=np.NaN

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
    Can be used on msprime output files but also e.g. msinfer output files
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
    return os.path.join(d,'+'.join(['fastarg', f, "fs"+str(seed)]))

def construct_argweaver_name(sim_name, seed):
    """
    Returns an ARGweaver inference filename (without file extension),
    based on a simulation name. The iteration number (used in .smc and .nex output)
    is not added here, but is added instead by the ARGweaver `arg-sample` program,
    in the format .10, .20, etc.
    """
    d,f = os.path.split(sim_name)
    return os.path.join(d,'+'.join(['aweaver', f, "ws"+str(seed)]))

def construct_msinfer_name(sim_name):
    """
    Returns an MSprime Li & Stevens inference filename.
    If the file is a subset of the original, this can be added to the
    basename later using the add_subsample_param_to_name() routine.
    """
    d,f = os.path.split(sim_name)
    return os.path.join(d,'+'.join(['msinfer', f, ""]))


def time_cmd(cmd, stdout):
    full_cmd = ["/usr/bin/time", "-f%M %S %U"] + cmd
    with tempfile.TemporaryFile() as stderr:
        exit_status = subprocess.call(full_cmd, stderr=stderr, stdout=stdout)
        stderr.seek(0)
        if exit_status != 0:
            raise ValueError(
                "Error running '{}': status={}:stderr{}".format(
                    cmd, exit_status, stderr.read()))

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

class Dataset(object):
    """
    A dataset is some collection of simulations which are run using
    the generate() method and stored in the raw_data_path directory.
    The process() method then processes the raw data and outputs
    the results into the data file.
    """
    name = None
    """
    Each dataset has a unique name. This is used as the prefix for the data
    file and raw_data_dir directory. Within this, replicate instances of datasets
    (each with a different RNG seed) are saved under the seed number
    """

    data_dir = "data"

    def __init__(self, seed):
        """
        If we initialize a dataset with the same seed, it should exactly
        duplicate the previous set of simulations (overwriting the old files)

        If we initialize a dataset with a different seed, it should create a new
        set of files which re-run the simulations w/ different throws of the dice.
        """
        self.instance = seed
        self.data_file = os.path.abspath(os.path.join(self.data_dir,
            "{}_{}.csv".format(self.name, self.instance)))
        self.raw_data_dir = os.path.join(self.data_dir, "raw__NOBACKUP__", self.name, str(self.instance))
        if not os.path.exists(self.raw_data_dir):
            logging.info("Making raw data dir {}".format(self.raw_data_dir))
        os.makedirs(self.raw_data_dir, exist_ok=True)

    def setup(self):
        """
        creates the simulations to use as a base
        """
        raise NotImplementedError()

    def generate(self):
        """
        generates the inferences
        """
        raise NotImplementedError()

    def process(self):
        """
        processes the inferred ARGs
        """
        raise NotImplementedError()

    def get_seeds(self, n, upper=1e7):
        """
        returns a set of n unique seeds suitable for seeding RNGs
        """
        np.random.seed(self.instance)
        seeds = np.random.randint(upper, size=n)
        seeds = np.unique(seeds)
        while len(seeds) != n:
            seeds = np.append(np.random.randint(upper,size = n-len(seeds)))
            seeds = np.unique(seeds)
        return seeds

    def single_simulation(self, n, Ne, l, rho, mu, seed, mut_seed=None, subsample=None):
        """
        The standard way to run one msprime simulation for a set of parameter values.
        Saves the output to an .hdf5 file, and also saves variant files for use in
        fastARG (a .hap file, in the format specified by https://github.com/lh3/fastARG#input-format)
        ARGweaver (a .sites file: http://mdrasmus.github.io/argweaver/doc/#sec-file-sites)
        msinfer (currently a numpy array containing the variant matrix)

        mutation_seed is not yet implemented, but should allow the same ancestry to be simulated (if the same
        genealogy_seed is given) but have different mutations thrown onto the trees (even with different mutation_rates)

        Returns a tuple of treesequence, filename (without file type extension)
        """
        ts = msprime.simulate(n, Ne, l, recombination_rate=rho, mutation_rate=mu, random_seed=int(seed))
        #here we might want to iterate over mutation rates for the same genealogy, setting a different mut_seed
        #but only so that we can see for ourselves the effect of mutation rate variation on a single topology
        #for the moment, we don't bother, and have mut_seed==genealogy_seed
        sim_fn = msprime_name(n, Ne, l, rho, mu, seed, seed, self.simulations_dir)
        logging.debug("writing {}.hdf5".format(sim_fn))
        ts.dump(sim_fn+".hdf5", zlib_compression=True)
        if subsample is not None:
            ts = ts.simplify(list(range(subsample)))
            sim_fn = add_subsample_param_to_name(sim_fn, subsample)
            ts.dump(sim_fn +".hdf5", zlib_compression=True)
        return ts, sim_fn

    @staticmethod
    def save_variant_matrices(ts, fname, error_rates):
        for error_rate in error_rates:
            S = generate_samples(ts, error_rate)
            err_filename = add_error_param_to_name(fname, error_rate)
            logging.debug("writing variant matrix to {}.npy for msinfer".format(err_filename))
            np.save(err_filename+".npy", S)

            pos = (v.position for v in ts.variants())
            logging.debug("writing variant matrix to {}.hap for fastARG".format(err_filename))
            with open(err_filename+".hap", "w+") as fastarg_in:
                msprime_fastARG.variant_matrix_to_fastARG_in(S, pos, fastarg_in)
            logging.debug("writing variant matrix to {}.sites for ARGweaver".format(err_filename))
            with open(err_filename+".sites", "w+") as argweaver_in:
                msprime_ARGweaver.variant_matrix_to_ARGweaver_in(S, pos, argweaver_in)

    @staticmethod
    def run_tsinf(S, rho):
        before = time.clock()
        panel = tsinf.ReferencePanel(S)
        P = panel.infer_paths(rho, num_workers=4)
        ts_new = panel.convert_records(P)
        ts_simplified = ts_new.simplify()
        cpu_time = time.clock() - before
        memory_use = 0 #To DO
        return ts_simplified, cpu_time, memory_use

    @staticmethod
    def run_fastarg(file_name, seed):
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
                    fa_out, var_pos, tree, muts, status_to=None)
            inferred_ts = msprime_fastARG.msprime_txts_to_fastARG_in_revised(
                    tree, muts, root_seq, fa_revised, simplify=True, status_to=None)
            assert filecmp.cmp(file_name, fa_revised.name, shallow=False)
            return inferred_ts, cpu_time, memory_use


class NumRecordsBySampleSizeDataset(Dataset):
    """
    Information on the number of coalescence records inferred by tsinf
    and FastARG for various sample sizes, under 3 different error rates
    """
    name = "num_records_by_sample_size"

    def __init__(self, seed=123, replicates=10):
        """
        Everything is done via a Data Frame which contains the initial params and is then used to
        store the output values.
        """
        super(NumRecordsBySampleSizeDataset, self).__init__(seed)
        self.simulations_dir = os.path.join(self.raw_data_dir, "simulations")

        #set up a pd Data Frame containing all the params
        params = collections.OrderedDict()
        params['sample_size']       = np.linspace(10, 500, num=10).astype(int)
        params['Ne']                 = 5000,
        params['length']             = 50000,
        params['recombination_rate'] = 2.5e-8,
        params['mutation_rate']      = 1.5e-8,
        params['replicate']         = np.arange(replicates)
        sim_params = list(params.keys())
        #now add some other params, where multiple vals exists for one simulation
        params['error_rate']         = [0, 0.01, 0.1]
        params['tool']               = ["fastARG", "msinfer"]

        #all combinations of params in a DataFrame
        self.data = expand_grid(params)
        #assign a unique index for each simulation
        self.data['sim'] = pd.Categorical(self.data[sim_params].astype(str).apply("".join, 1)).codes
        #set unique seeds for each sim
        self.data['seed'] = self.get_seeds(max(self.data.sim)+1)[self.data.sim]

    def setup(self):
        """
        Generates the initial simulations from which we can infer data.
        Should be quite fast, but will generate lots of files, and probably
        take up a fair bit of disk space.
        """
        # clear out any old simulations to avoid confusion.
        if os.path.exists(self.simulations_dir):
            shutil.rmtree(self.simulations_dir)
        os.makedirs(self.simulations_dir)
        for s in range(max(self.data.sim)+1):
            sim = self.data[self.data.sim == s]
            logging.info("Running simulation for n = {}".format(pd.unique(sim.sample_size)))
            #note to Jerome - I don't use the num_replicates option in simulate(), although it is presumably more
            # efficient, as then I can't replicate each msprime simulation independently, given a RNG seed.
            ts, fn = self.single_simulation(pd.unique(sim.sample_size),
                                            pd.unique(sim.Ne),
                                            pd.unique(sim.length),
                                            pd.unique(sim.recombination_rate),
                                            pd.unique(sim.mutation_rate),
                                            pd.unique(sim.seed),
                                            pd.unique(sim.seed)) #same mutation_seed as genealogy_seed
            self.save_variant_matrices(ts, fn, pd.unique(sim.error_rate))

    def generate(self, save_trees=False):
        """
        This runs the inference methods. It will be the most time consuming bit.
        The final files saved will be nexus files containing multiple trees.
        In the case of fastARG and msinfer methods, we could also save hdf5 files for reference

        Note that we should be able to kill a python instance after 'generate()' has run,
        and fire up another instance in which we run 'process()'.

        """
        add_columns(self.data, ['source_records','inferred_records','cpu_time','memory'])
        for i in self.data.index:
            d = self.data.iloc[i]
            sim_fn = msprime_name(d.sample_size, d.Ne, d.length, d.recombination_rate,
                                  d.mutation_rate, d.seed, d.seed, self.simulations_dir)
            err_fn = add_error_param_to_name(sim_fn, d.error_rate)
            ts = msprime.load(sim_fn+".hdf5")
            inferred_ts = None
            assert ts.sample_size == d.sample_size
            if d.tool == 'msinfer':
                infile = err_fn + ".npy"
                out_fn = construct_msinfer_name(err_fn)
                logging.info("processing msinfer inference for n = {}".format(d.sample_size))
                logging.debug("reading: {} for msprime inference".format(infile))
                S = np.load(infile)
                assert S.shape == (ts.sample_size, ts.num_mutations)
                inferred_ts, time, memory = self.run_tsinf(S, 4*d.recombination_rate*d.Ne)
                inferred_ts.dump(out_fn +".hdf5", zlib_compression=True)
            elif d.tool == 'fastARG':
                infile = err_fn + ".hap"
                out_fn = construct_fastarg_name(err_fn, d.seed)
                logging.info("processing fastARG inference for n = {}".format(d.sample_size))
                logging.debug("reading: {} for msprime inference".format(infile))
                inferred_ts, time, memory = self.run_fastarg(infile, d.seed)
            else:
                raise ValueError("Unknown tool: {}".format(tool))
            self.data.loc[i,['source_records','inferred_records','cpu_time','memory']] = \
                (ts.get_num_records(), inferred_ts.get_num_records(),time, memory)

            # Save each row so we can use the information while it's being built
            self.data.to_csv(self.data_file)

class MetricsByMutationRate(Dataset):
    """
    Dataset for Figure 1
    Accuracy of ARG inference (measured by various statistics)
    tending to fully accurate as mutation rate increases
    """
    name = "num_records_by_sample_size"

    def __init__(self, seed=111):
        """
        Everything is done via a Data Frame which contains the initial params and is then used to
        store the output values.
        """
        super(NumRecordsBySampleSizeDataset, self).__init__(seed)
        self.simulations_dir = os.path.join(self.raw_data_dir, "simulations")

        #set up a pd Data Frame containing all the params
        params = collections.OrderedDict()
        params['sample_size']       = 12
        params['Ne']                 = 5000,
        params['length']             = 50000,
        params['recombination_rate'] = 2.5e-8,
        params['mutation_rate']      = np.linspace(1.5e-8, 500, num=10)
        params['replicate']         = np.arange(100)
        sim_params = list(params.keys())
        #now add some other params, where multiple vals exists for one simulation
        params['error_rate']         = [0, 0.01, 0.1]
        params['tool']               = ["fastARG","ARGweaver", "msinfer"]

        #all combinations of params in a DataFrame
        self.data = expand_grid(params)
        #assign a unique index for each simulation
        self.data['sim'] = pd.Categorical(self.data[sim_params].astype(str).apply("".join, 1)).codes
        #set unique seeds for each sim
        self.data['seed'] = self.get_seeds(max(self.data.sim)+1)[self.data.sim]





def run_setup(cls, extra_args):
    logging.info("Setting up base data for {}{}".format(cls.name,
        " with extra arguments {}".format(extra_args) if extra_args else ""))
    f = cls(**extra_args)
    f.setup()

def run_generate(cls, extra_args):
    logging.info("Generating {}".format(cls.name))
    f = cls()
    f.generate()


def run_process(cls, extra_args):
    logging.info("Processing {}".format(cls.name))
    f = cls()
    f.process()


def run_plot(cls, extra_args):
    f = cls()
    f.plot()


def main():
    datasets = [
        NumRecordsBySampleSizeDataset,
    ]
    figures = [
    ]
    name_map = dict([(d.name, d) for d in datasets + figures])
    parser = argparse.ArgumentParser(
        description= "Set up base data, generate inferred datasets, process datasets and plot figures.")
    parser.add_argument('--verbosity', '-v', action='count', default=0)
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    subparser = subparsers.add_parser('setup')
    # TODO: something like this will be useful to set up smaller runs for
    # testing purposes and to control the number of processes used.
    # subparser.add_argument(
    #     "--processes", '-p', help="number of processes",
    #     type=int, default=1)
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
    subparser.add_argument(
         '--replicates', '-r', help="number of replicates", type=int)
    subparser.add_argument(
         '--seed', '-s', help="use a non-default RNG seed", type=int)
    subparser.set_defaults(func=run_setup)

    subparser = subparsers.add_parser('generate')
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
    subparser.set_defaults(func=run_generate)

    subparser = subparsers.add_parser('process')
    subparser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
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

    base_args = ['name','verbosity','command','func']
    extra_args = {k:v for k,v in vars(args).items() if k not in base_args and v is not None}
    k = args.name[0]
    if k == "all":
        classes = datasets
        if args.func == run_plot:
            classes = figures
        for name, cls in name_map.items():
            if cls in classes:
                args.func(cls, extra_args)
    else:
        cls = name_map[k]
        args.func(cls, extra_args)

if __name__ == "__main__":
    main()
