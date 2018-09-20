"""
Script to produce the figure showing how much space we need
in principle to store the ancestry of 10 billion humans.
"""
import time
import os
import os.path
import argparse
import subprocess

import numpy as np
import msprime
import scipy.optimize as optimize
import scipy.stats as stats
import pandas as pd
import humanize

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot

datafile = "data/storing_everyone.csv"
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
    trees_file = os.path.join(data_prefix, "{}.trees".format(sample_size))
    ts.dump(trees_file)
    subprocess.check_call(["gzip", "-k", trees_file])


def run_simulate():
    for k in range(1, 8):
        run_simulation(10**k)


def run_benchmark():
    print("msprime version:", msprime.__version__)

    before = time.perf_counter()
    filename = os.path.join(data_prefix, "{}.trees".format(10**7))
    ts = msprime.load(filename)
    duration = time.perf_counter() - before
    print("loaded {} tree sequense in {:.2f}s".format(
        humanize.naturalsize(os.path.getsize(filename)), duration))

    size = ts.num_samples * ts.num_sites
    print("Total size of genotype matrix = ", humanize.naturalsize(size))

    before = time.perf_counter()
    j = 0
    for tree in ts.trees(sample_counts=False, sample_lists=False):
        j += 1
    assert j == ts.num_trees
    duration = time.perf_counter() - before
    print("Iterated over {} trees in {:.2f}s".format(ts.num_trees, duration))

    before = time.perf_counter()
    samples = np.arange(10**6)
    freq = np.zeros(ts.num_sites)
    for tree in ts.trees(tracked_samples=samples):
        for site in tree.sites():
            node = site.mutations[0].node
            freq[site.id] = tree.num_tracked_samples(node)
    duration = time.perf_counter() - before
    print("Computed {} allele frequencies in {:.2f}s".format(ts.num_sites, duration))


def run_process():
    sample_size = []
    uncompressed = []
    compressed = []
    vcf = []
    bcf = []
    for k in range(1, 8):
        n = 10**k
        filename = os.path.join(data_prefix, "{}.trees".format(n))
        if not os.path.exists(filename):
            break
        sample_size.append(n)
        uncompressed.append(os.path.getsize(filename))
        ts = msprime.load(filename)
        filename += ".gz"
        compressed.append(os.path.getsize(filename))
        if k < 7:
            filename = os.path.join(data_prefix, "{}.vcf".format(n))
            with open(filename, "w") as vcf_file:
                ts.write_vcf(vcf_file, 2)
            print("Wrote ", filename)
            vcf.append(os.path.getsize(filename))

            bcf_filename = os.path.join(data_prefix, "{}.bcf".format(n))
            subprocess.check_call(["bcftools view -O b {} > {}".format(
                filename, bcf_filename)], shell=True)
            bcf.append(os.path.getsize(bcf_filename))
        else:
            vcf.append(0)
            bcf.append(0)


    GB = 1024**3
    sample_size = np.array(sample_size)
    compressed = np.array(compressed) / GB
    uncompressed = np.array(uncompressed) / GB
    vcf = np.array(vcf) / GB
    bcf = np.array(bcf) / GB

    df = pd.DataFrame({
        "sample_size": sample_size,
        "compressed": compressed,
        "uncompressed": uncompressed,
        "vcf": vcf, "bcf": bcf})
    print(df)
    df.to_csv(datafile)


def run_plot():
    GB = 1024**3
    # Fit the model for the observed data. Wek
    rho = 4 * Ne * recombination_rate * length

    def additive_model(n, a, b, c):
        # The model is based on the proof that we need O(n + rho log n) space
        # to store edge. We then need a further O(theta log n) space to
        # store the mutations. Because rho == mu, we don't need a separate term
        # for rho, and represent them both by z. Thus, a here is the per-edge
        # cost, b is the per site cost, and c is a constant for the overheads.
        z = rho * np.log(n)
        return a * (n + z) + b * z + c

    def mulplicative_model(n, a, b, c):
        # This is the number of sites.
        z = rho * np.log(n)
        return a * n * z + b * z + c

    df = pd.read_csv(datafile)
    df = df[df.sample_size > 100]
    projected_n = 10**np.arange(3, 11)

    msp_fit_params, _ = optimize.curve_fit(
        additive_model, df.sample_size[:-1], df.uncompressed[:-1])
    msp_fit = additive_model(projected_n, *msp_fit_params)
    mspz_fit_params, _ = optimize.curve_fit(
        additive_model, df.sample_size[:-1], df.compressed[:-1])
    mspz_fit = additive_model(projected_n, *mspz_fit_params)

    index = df.sample_size < 10**6
    vcf_fit_params, _ = optimize.curve_fit(
        mulplicative_model, df.sample_size[index], df.vcf[index])
    vcf_fit = mulplicative_model(projected_n, *vcf_fit_params)

    bcf_fit_params, _ = optimize.curve_fit(
        mulplicative_model, df.sample_size[index], df.bcf[index])
    bcf_fit = mulplicative_model(projected_n, *bcf_fit_params)
    bcf_fit_params, _ = optimize.curve_fit(
        mulplicative_model, df.sample_size[index], df.bcf[index])

    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)
    xytext = (18, 0)

    index = df.vcf > 0
    line, = ax1.loglog(df.sample_size[index], df.vcf[index], "^", label="VCF")
    ax1.loglog(projected_n, vcf_fit, "--", color=line.get_color())
    ax1.annotate(
        humanize.naturalsize(vcf_fit[-1] * GB, binary=True, format="%d"),
        textcoords="offset points", xytext=xytext,
        xy=(projected_n[-1], vcf_fit[-1]), xycoords="data")

    line, = ax1.loglog(df.sample_size[index], df.bcf[index], "s", label="BCF")
    ax1.loglog(projected_n, bcf_fit, "--", color=line.get_color())
    ax1.annotate(
        humanize.naturalsize(bcf_fit[-1] * GB, binary=True, format="%d"),
        textcoords="offset points", xytext=xytext,
        xy=(projected_n[-1], bcf_fit[-1]), xycoords="data")

    line, = ax1.loglog(
        df.sample_size, df.uncompressed, "o", label="msprime uncompressed")
    ax1.loglog(projected_n, msp_fit, "--", color=line.get_color())
    ax1.annotate(
        humanize.naturalsize(msp_fit[-1] * GB, binary=True, format="%d"),
        textcoords="offset points", xytext=xytext,
        xy=(projected_n[-1], msp_fit[-1]), xycoords="data")

    line, = ax1.loglog(
        df.sample_size, df.compressed, "*", label="msprime compressed")
    ax1.loglog(projected_n, mspz_fit, "--", color=line.get_color())
    ax1.annotate(
        humanize.naturalsize(mspz_fit[-1] * GB, binary=True, format="%d"),
        textcoords="offset points", xytext=xytext,
        xy=(projected_n[-1], mspz_fit[-1]), xycoords="data")

    ax1.set_xlabel("Number of chromosomes")
    ax1.set_ylabel("File size (GiB)")
    pyplot.legend()
    # pyplot.tight_layout()
    pyplot.savefig("figures/storing_everyone.pdf")


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

    subparser = subparsers.add_parser("benchmark")
    subparser.set_defaults(func=run_benchmark)

    subparser = subparsers.add_parser("plot")
    subparser.set_defaults(func=run_plot)

    args = parser.parse_args()
    args.func()
