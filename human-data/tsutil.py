"""
Various utilities for manipulating tree sequences and running tsinfer.
"""
import argparse

import msprime
import tsinfer
import daiquiri
import numpy as np

def run_simplify(args):
    print("simplify", args.input, args.output)
    ts = msprime.load(args.input)
    ts = ts.simplify()
    ts.dump(args.output)


def run_augment(sample_data, ancestors_ts, subset, num_threads):
    progress_monitor = tsinfer.cli.ProgressMonitor(enabled=True, augment_ancestors=True)
    return tsinfer.augment_ancestors(
        sample_data, ancestors_ts, subset, num_threads=num_threads,
        progress_monitor=progress_monitor)

def run_match_samples(sample_data, ancestors_ts, num_threads):
    progress_monitor = tsinfer.cli.ProgressMonitor(enabled=True, match_samples=True)
    return tsinfer.match_samples(
        sample_data, ancestors_ts, num_threads=num_threads,
        simplify=False, progress_monitor=progress_monitor)


def run_sequential_augment(args):

    base = ".".join(args.input.split(".")[:-1])

    sample_data = tsinfer.load(args.input)
    num_samples = sample_data.num_samples
    ancestors_ts = msprime.load(base + ".ancestors.trees")

    samples = np.arange(num_samples, dtype=int)
    n = 2
    offset = -1
    while n < num_samples // 4:
        augmented_file = base + ".augmented_{}.ancestors.trees".format(n)
        final_file = base + ".augmented_{}.unsimplified.trees".format(n)
        print("RUNNING", augmented_file)
        subset = np.linspace(offset, num_samples - offset, n + 2, dtype=int)[1: -1]
        ancestors_ts = run_augment(sample_data, ancestors_ts, subset, args.num_threads)
        ancestors_ts.dump(augmented_file)
        n *= 2
        offset += 1

    final_ts = run_match_samples(sample_data, ancestors_ts, args.num_threads)
    final_ts.dump(final_file)


def main():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = "command"

    subparser = subparsers.add_parser("simplify")
    subparser.add_argument(
        "input", type=str, help="Input tree sequence")
    subparser.add_argument(
        "output", type=str, help="Input tree sequence")
    subparser.set_defaults(func=run_simplify)

    subparser = subparsers.add_parser("sequential-augment")
    subparser.add_argument(
        "input", type=str, help="Input tree sequence")
    subparser.add_argument("--num-threads", type=int, default=0)
    subparser.set_defaults(func=run_sequential_augment)


    daiquiri.setup(level="INFO")

    args = parser.parse_args()
    args.func(args)
     
main()
