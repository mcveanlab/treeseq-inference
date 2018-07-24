"""
Simple CLI to run tsinf on the command line.
"""
import sys
import os
import argparse
import logging

import numpy as np

# use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(sys.path[0],'..','msprime'))
sys.path.insert(1,os.path.join(sys.path[0],'..','tsinfer'))
import msprime
import tsinfer
import tsinfer.eval_util as eval_util
import tsinfer.formats as formats

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
        help="The samples file name, as saved by tsinfer.SampleData.initialise()")
    parser.add_argument(
        "output",
        help="The path to write the output file to")
    parser.add_argument(
        "-srb", "--shared-recombinations", action='store_true',
        help="Use shared recombinations (path compression) to break polytomies")
    parser.add_argument(
        "-l", "--length", default=None, type=int,
        help="The total sequence length")
    parser.add_argument(
        "-r", "--recombination-rate", default=None, type=float,
        help="The scaled recombination rate (now ignored)")
    parser.add_argument(
        "-e", "--error-probability", default=None, type=float,
        help="The probability of observing an error (now ignored)")
    parser.add_argument(
        "-t", "--threads", default=1, type=int,
        help="The number of worker threads to use")
    parser.add_argument(
        "-m", "--method", default="C", choices=['C','P'],
        help="Which implementation to use, [C] (faster) or [P]ython (more debuggable)")
    parser.add_argument(
        "--inject-real-ancestors-from-ts", default=None,
        help="Instead of inferring ancestors, construct known ones from this tree sequence file path")
    parser.add_argument(
        "-V", "--version", action='version', version=description)
    

    args = parser.parse_args()
        
    if args.recombination_rate is not None:
        logging.warning("TSinfer now simply ignores recombination rate. You can omit this parameter")
    if args.error_probability is not None:
        logging.warning("TSinfer now simply ignores error probabilities. You can omit this parameter")
    engine = tsinfer.PY_ENGINE if args.method == "P" else tsinfer.C_ENGINE

    if not os.path.isfile(args.samples):
        raise ValueError("No samples file")
    sample_data = tsinfer.load(args.samples)
    if all(False for _ in sample_data.genotypes(inference_sites=True)):
        raise ValueError("No inference sites")
    if args.inject_real_ancestors_from_ts is not None:
        ancestor_data = tsinfer.AncestorData.initialise(sample_data, compressor=None)
        orig_ts = msprime.load(args.inject_real_ancestors_from_ts)
        eval_util.build_simulated_ancestors(sample_data, ancestor_data, orig_ts)
        ancestor_data.finalise()
        ancestors_ts = tsinfer.match_ancestors(
            sample_data, ancestor_data, path_compression=args.shared_recombinations, engine=engine)
        ts = tsinfer.match_samples(
            sample_data, ancestors_ts, path_compression=args.shared_recombinations, engine=engine, simplify=True)
    else:
        ts = tsinfer.infer(
            sample_data, num_threads=args.threads, path_compression=args.shared_recombinations, engine=engine)
    ts.dump(args.output)

    # # TODO add command line arg here for when we're comparing run time performance.
    # # Quickly verify that we get the sample output.
    # Sp = np.zeros(S.shape)
    # for j, h in enumerate(ts.haplotypes()):
    #     Sp[j] = np.fromstring(h, np.uint8) - ord('0')
    # assert np.all(Sp == S)

if __name__ == "__main__":
    main()
