"""
Various utilities for manipulating tree sequences and running tsinfer.
"""
import argparse

import msprime
import tsinfer


def run_simplify(args):
    print("simplify", args.input, args.output)
    ts = msprime.load(args.input)
    ts = ts.simplify()
    ts.dump(args.output)


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

    args = parser.parse_args()
    args.func(args)

main()
