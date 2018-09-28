"""
Quick script to benchmark bgen reading. Uses the local simplebgen
module which is a thin wrapper around the bgen library.
"""
import sys
import os
import argparse
import time

import simplebgen

parser = argparse.ArgumentParser()

parser.add_argument("input", type=str, help="Input bgen file")
parser.add_argument(
    "--num-variants", type=int, default=None,
    help="Number of variants to benchmark genotypes decoding performance on")

args = parser.parse_args()

bg = simplebgen.BgenReader(args.input)
print("PID = ", os.getpid())
print("num_samples  = ", bg.num_samples)
print("num_variants = ", bg.num_variants)
num_variants = 0
before = time.perf_counter()
for j in range(bg.num_variants):
    if num_variants == args.num_variants:
        break
    P = bg.get_probabilities(j)
    num_variants += 1

duration = time.perf_counter() - before
total_genotypes = (2 * bg.num_samples * num_variants) / 10**6
print("Iterated over {} variants in {:.2f}s @ {:.2f} M genotypes/s".format(
    num_variants, duration, total_genotypes / duration))
