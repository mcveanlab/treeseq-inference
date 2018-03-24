"""
Script to produce the figure showing how much space we need
in principle to store the ancestry of 10 billion humans.
"""
import time
import os
import os.path
import argparse

import numpy as np
import msprime
import scipy.optimize as optimize


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot

data_prefix = "tmp__NOBACKUP__"
mutation_rate = 1e-8
recombination_rate = 1e-8
length = 100 * 10**6
Ne = 10**4

def run_simulation(sample_size):
    before = time.perf_counter()
    ts = msprime.simulate(
        sample_size=sample_size, Ne=Ne, mutation_rate=mutation_rate,
        length=length, recombination_rate=recombination_rate)
    duration = time.perf_counter() - before
    print("Simulated {} in {} hours".format(sample_size, duration / 3600))
    ts.dump(os.path.join(data_prefix, "{}.uncompressed.hdf5".format(sample_size)))
    ts.dump(
        os.path.join(data_prefix, "{}.compressed.hdf5".format(sample_size)),
        zlib_compression=True)

def run_simulate():
    for k in range(1, 8):
        run_simulation(10**k)

def run_process():
    sample_size = []
    uncompressed = []
    compressed = []
    rho = 4 * Ne * recombination_rate * length
    mu = 4 * Ne * mutation_rate * length
    for k in range(1, 8):
        n = 10**k
        filename = os.path.join(data_prefix, "{}.uncompressed.hdf5".format(n))
        if not os.path.exists(filename):
            break
        ts = msprime.load(filename)
        sample_size.append(n)
        uncompressed.append(os.path.getsize(filename))
        filename = os.path.join(data_prefix, "{}.compressed.hdf5".format(n))
        compressed.append(os.path.getsize(filename))
    GB = 1024**3
    sample_size = np.array(sample_size)
    compressed = np.array(compressed) / GB
    uncompressed = np.array(uncompressed) / GB

    df = pd.DataFrame({
        "sample_size": sample_size,
        "compressed": compressed,
        "uncompressed": uncompressed})



def run_plot():

    def func(n, a, b, c):
        z = rho * np.log(n)
        return a * (n + z) + b * z + c

    popt, pcov = optimize.curve_fit(func, sample_size[:-1], uncompressed[:-1])
    print("a = ", popt[0], "b = ", popt[1], "c = ", popt[2])

    projected_n = 10**np.arange(1, 11)
    fitted = func(projected_n, *popt)
    pyplot.loglog(projected_n, fitted, "o", ls="--", label="Model")
    print(
        "predicted file size for 10 billion samples is ",
        func(10**10, *popt), "gigabytes")

    pyplot.loglog(sample_size, uncompressed, label="Uncompressed file size")
    # pyplot.loglog(sample_size, compressed)

    pyplot.xlabel("Number of chromosomes")
    pyplot.ylabel("File size (GiB)")
    pyplot.legend()
    pyplot.tight_fit()
    pyplot.savefig("figures/ancestry_of_everyone.pdf")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the plot showing the files size storing everyone")
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    subparser = subparsers.add_parser("simulate")
    subparser.set_defaults(func=run_simulate)

    subparser = subparsers.add_parser("process")
    subparser.set_defaults(func=run_process)

    args = parser.parse_args()
    args.func()
