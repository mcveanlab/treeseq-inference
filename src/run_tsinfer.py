"""
Simple CLI to run tsinf on the command line.
"""
import sys
import os
import argparse

import numpy as np

# use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(sys.path[0],'..','msprime'))
sys.path.insert(1,os.path.join(sys.path[0],'..','tsinfer'))
import msprime
import tsinfer

def main():

    description = """Simple CLI wrapper for tsinfer
        msprime version: {}
        tsinfer version: {}""".format(msprime.__version__, tsinfer.__version__)
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbosity', '-v', action='count', default=0)
    parser.add_argument(
        "samples",
        help="The observed haplotypes as a numpy array file")
    parser.add_argument(
        "positions",
        help="The positions of sites as a numpy array file")
    parser.add_argument(
        "output",
        help="The path to write the output file to")
    parser.add_argument(
        "-srb", "--shared-recombinations", action='store_true',
        help="Use shared recombinations to break polytomies")
    parser.add_argument(
        "-sl", "--shared-lengths", action='store_true',
        help="Use shorter shared lengths of sequence to break polytomies")
    parser.add_argument(
        "-l", "--length", default=None, type=int,
        help="The total sequence length")
    parser.add_argument(
        "-r", "--recombination-rate", default=1, type=float,
        help="The scaled recombination rate.")
    parser.add_argument(
        "-e", "--error-probability", default=0, type=float,
        help="The probability of observing an error")
    parser.add_argument(
        "-t", "--threads", default=1, type=int,
        help="The number of worker threads to use")

    args = parser.parse_args()
    S = np.load(args.samples)
    pos = np.load(args.positions)
    # We need to transpose this now as
    genotypes = S.astype(np.uint8).T
    ts = tsinfer.infer(
        genotypes, pos, args.length,
        args.recombination_rate,
        args.error_probability,
        num_threads=args.threads)
    ts.dump(args.output)

    # # TODO add command line arg here for when we're comparing run time performance.
    # # Quickly verify that we get the sample output.
    # Sp = np.zeros(S.shape)
    # for j, h in enumerate(ts.haplotypes()):
    #     Sp[j] = np.fromstring(h, np.uint8) - ord('0')
    # assert np.all(Sp == S)

if __name__ == "__main__":
    main()
