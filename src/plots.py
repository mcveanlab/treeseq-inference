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

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import pandas as pd
import seaborn as sns

# TODO force import of local version
import msprime
import tsinf

if sys.version_info[0] < 3:
    raise Exception("Python 3 only")


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
    Each dataset a unique name. This is used as the prefix for the data
    file and raw_data_dir directory.
    """

    data_dir = "data"

    def __init__(self):
        self.data_file = os.path.join(self.data_dir, "{}.csv".format(self.name))
        self.raw_data_dir = os.path.join(self.data_dir, "raw__NOBACKUP__", self.name)
        if not os.path.exists(self.raw_data_dir):
            logging.info("Making raw data dir {}".format(self.raw_data_dir))
        os.makedirs(self.raw_data_dir, exist_ok=True)

    def generate(self):
        raise NotImplementedError()

    def process(self):
        raise NotImplementedError()

class NumRecordsBySampleSizeDataset(Dataset):
    """
    Information on the number of coalescence records inferred by tsinf
    and FastARG for various sample sizes.
    """
    name = "num_records_by_sample_size"

    def __init__(self):
        super(NumRecordsBySampleSizeDataset, self).__init__()
        self.simulations_dir = os.path.join(self.raw_data_dir, "simulations")
        self.tree_sequence_filename_pattern = os.path.join(
            self.simulations_dir, "n={}_rep={}.hdf5")
        self.samples_filename_pattern = os.path.join(
            self.simulations_dir, "n={}_rep={}_error={}.npy")
        self.sample_sizes = np.linspace(10, 500, num=10).astype(int)
        self.error_rates = [0, 0.01, 0.1]
        self.mutation_rate = 1.5
        self.recombination_rate = 2.5
        self.length = 50
        self.num_replicates = 10

    def run_tsinf(self, S):
        before = time.clock()
        panel = tsinf.ReferencePanel(S)
        P = panel.infer_paths(self.recombination_rate, num_workers=4)
        ts_new = panel.convert_records(P)
        ts_simplified = ts_new.simplify()
        cpu_time = time.clock() - before
        return ts_simplified.get_num_records(), cpu_time, 0

    def generate(self):
        # clear out any old simulations to avoid confusion.
        if os.path.exists(self.simulations_dir):
            shutil.rmtree(self.simulations_dir)
        os.makedirs(self.simulations_dir)
        for n in self.sample_sizes:
            logging.info("Running simulation for n = {}".format(n))
            replicates = msprime.simulate(
                n, Ne=0.5, mutation_rate=self.mutation_rate, length=self.length,
                recombination_rate=self.recombination_rate,
                num_replicates=self.num_replicates)
            for j, ts in enumerate(replicates):
                filename = self.tree_sequence_filename_pattern.format(n, j)
                logging.debug("writing {}".format(filename))
                ts.dump(filename, zlib_compression=True)
                for error_rate in self.error_rates:
                    S = generate_samples(ts, error_rate)
                    filename = self.samples_filename_pattern.format(n, j, error_rate)
                    logging.debug("writing {}".format(filename))
                    np.save(filename, S)

    def process(self):
        df = pd.DataFrame(columns=(
            "tool", "num_samples", "error_rate", "replicate", "source_records",
            "inferred_records", "cpu_time", "memory"))
        row_index = 0
        for n in self.sample_sizes:
            logging.info("processing for n = {}".format(n))
            for j in range(self.num_replicates):
                filename = self.tree_sequence_filename_pattern.format(n, j)
                ts = msprime.load(filename)
                assert ts.sample_size == n
                for error_rate in self.error_rates:
                    filename = self.samples_filename_pattern.format(n, j, error_rate)
                    logging.debug("reading: {}".format(filename))
                    S = np.load(filename)
                    assert S.shape == (ts.sample_size, ts.num_mutations)
                    inferred_records, time, memory = self.run_tsinf(S)
                    df.loc[row_index] = (
                        "tsinf", n, error_rate, j, ts.get_num_records(),
                        inferred_records, time, memory)
                    row_index += 1
                    # Save each row so we can use the information while it's being built
                    df.to_csv(self.data_file)


def run_generate(cls, args):
    logging.info("Generating {}".format(cls.name))
    f = cls()
    f.generate()


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
    ]
    figures = [
    ]
    name_map = dict([(d.name, d) for d in datasets + figures])
    parser = argparse.ArgumentParser(
        description= "Generate datasets, process raw data and generate figures.")
    parser.add_argument('--verbose', '-v', action='count', default=0)
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    generate_parser = subparsers.add_parser('generate')
    # TODO: something like this will be useful to run smaller runs for
    # testing purposes and to control the number of processes used.
    # generate_parser.add_argument(
    #     '-n', help="number of replicates", type=int, default=-1)
    # generate_parser.add_argument(
    #     "--processes", '-p', help="number of processes",
    #     type=int, default=1)
    generate_parser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the dataset identifier')
    generate_parser.set_defaults(func=run_generate)

    process_parser = subparsers.add_parser('process')
    process_parser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the simulation identifier')
    process_parser.set_defaults(func=run_process)

    figure_parser = subparsers.add_parser('figure')
    figure_parser.add_argument(
        'name', metavar='NAME', type=str, nargs=1,
        help='the figure dentifier')
    figure_parser.set_defaults(func=run_plot)

    args = parser.parse_args()
    log_level = logging.WARNING
    if args.verbose == 1:
        log_level = logging.INFO
    if args.verbose >= 2:
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
