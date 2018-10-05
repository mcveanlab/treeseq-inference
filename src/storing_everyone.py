"""
Script to produce the figure showing how much space we need
in principle to store the ancestry of 10 billion humans.
"""
import time
import os
import os.path
import argparse
import subprocess
import io

import numpy as np
import msprime
import scipy.optimize as optimize
import pandas as pd
import humanize
import cyvcf2
# Used for newick
from Bio import Phylo

datafile = "data/storing_everyone.csv"
data_prefix = "data/raw__NOBACKUP__/storing_everyone"
mutation_rate = 1e-8
recombination_rate = 1e-8
length = 100 * 10**6
Ne = 10**4

if not os.path.exists(data_prefix):
    os.mkdir(data_prefix)


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


def benchmark_bcf(ts):
    total_sites = ts.num_sites
    num_sites = 10**4
    vcf_filename = os.path.join(data_prefix, "large-subset.vcf")
    bcf_filename = os.path.join(data_prefix, "large-subset.bcf")
    if not os.path.exists(vcf_filename):

        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        for site in ts.sites():
            if site.id == num_sites:
                break
            site_id = tables.sites.add_row(
                site.position, ancestral_state=site.ancestral_state,
                metadata=site.metadata)
            for mutation in site.mutations:
                tables.mutations.add_row(
                    site_id, node=mutation.node, parent=mutation.parent,
                    derived_state=mutation.derived_state,
                    metadata=mutation.metadata)
        ts = tables.tree_sequence()
        print("Subsetted tree sequence")

        with open(vcf_filename, "w") as vcf_file:
            ts.write_vcf(vcf_file, 2)
        print("Wrote ", vcf_filename)
        subprocess.check_call(["bcftools view -O b {} > {}".format(
            vcf_filename, bcf_filename)], shell=True)
        print("Wrote ", bcf_filename)

    before = time.perf_counter()
    records = cyvcf2.VCF(bcf_filename)
    duration = time.perf_counter() - before
    print("Read BCF header in {:.2f} seconds".format(duration))
    before = time.perf_counter()
    count = 0
    for record in records:
        count += 1
    assert count == num_sites
    duration = time.perf_counter() - before
    print("Read {} BCF records in {:.2f} seconds".format(count, duration))

    estimated_time = total_sites * (duration / num_sites)
    print("Esimated time to read {} BCF records = {:.2f} hours".format(
        total_sites, estimated_time / 3600))


def run_benchmark_newick(ts, num_trees):

    total_length = 0
    total_duration = 0
    for tree in ts.trees():
        if num_trees == tree.index:
            break
        ns = tree.newick()
        handle = io.StringIO(ns)
        before = time.perf_counter()
        Phylo.read(handle, "newick")
        total_duration += time.perf_counter() - before
        total_length += len(ns)
    else:
        raise ValueError("not enough trees!")

    expected_size = (total_length / num_trees) * ts.num_trees
    mean_time = total_duration / num_trees
    print("Expected size of newick file: {:.2f}TiB".format(expected_size / 1024**6))
    print("Expected time to read {} hours @ {:.2f} secs each".format(
        (mean_time * ts.num_trees) / 3600, mean_time))


def run_benchmark():
    print("msprime version:", msprime.__version__)

    before = time.perf_counter()
    filename = os.path.join(data_prefix, "{}.trees".format(10**7))
    ts = msprime.load(filename)
    duration = time.perf_counter() - before
    print("loaded {} tree sequence in {:.2f}s".format(
        humanize.naturalsize(os.path.getsize(filename), binary=True), duration))

    run_benchmark_newick(ts, 2)

    size = ts.num_samples * ts.num_sites
    print("Total size of genotype matrix = ", humanize.naturalsize(size, binary=True))

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

    benchmark_bcf(ts)


def run_convert_files():
    for k in range(1, 8):
        n = 10**k
        filename = os.path.join(data_prefix, "{}.trees".format(n))
        if not os.path.exists(filename):
            break
        ts = msprime.load(filename)
        filename += ".gz"
        if k < 7:
            filename = os.path.join(data_prefix, "{}.vcf".format(n))
            with open(filename, "w") as vcf_file:
                ts.write_vcf(vcf_file, 2)
            print("Wrote ", filename)
            gz_filename = filename + ".gz"
            subprocess.check_call("gzip -c {} > {}".format(filename, gz_filename), shell=True)
            print("Wrote ", gz_filename)


def run_make_data():
    sample_size = 10**np.arange(1, 11)
    uncompressed = np.zeros(sample_size.shape)
    compressed = np.zeros(sample_size.shape)
    vcf = np.zeros(sample_size.shape)
    vcfz = np.zeros(sample_size.shape)

    GB = 1024**3
    for j, n in enumerate(sample_size):
        files = [
            (uncompressed, os.path.join(data_prefix, "{}.trees".format(n))),
            (compressed, os.path.join(data_prefix, "{}.trees.gz".format(n))),
            (vcf, os.path.join(data_prefix, "{}.vcf".format(n))),
            (vcfz, os.path.join(data_prefix, "{}.vcf.gz".format(n)))]
        for array, filename in files:
            if os.path.exists(filename):
                array[j] = os.path.getsize(filename) / GB

    # Fit the model for the observed data.
    rho = 4 * Ne * recombination_rate * length

    def additive_model(n, a, b, c):
        # The model is based on the proof that we need O(n + rho log n) space
        # to store edge. We then need a further O(theta log n) space to
        # store the mutations. Because rho == mu, we don't need a separate term
        # for rho, and represent them both by z. Thus, a here is the per-edge
        # cost, b is the per site cost, and c is a constant for the overheads.
        z = rho * np.log(n)
        return a * (n + z) + b * z + c

    def mulplicative_model(n, a, b):
        # Fit a simple exponential function.
        return a * np.power(n, b)

    # Keep on data point back from all these fits so that we can see how the
    # fit is doing.
    index = uncompressed > 0
    tsk_fit_params, _ = optimize.curve_fit(
        additive_model, sample_size[index][:-1], uncompressed[index][:-1])
    tsk_fit = additive_model(sample_size, *tsk_fit_params)
    tskz_fit_params, _ = optimize.curve_fit(
        additive_model, sample_size[index][:-1], compressed[index][:-1])
    tskz_fit = additive_model(sample_size, *tskz_fit_params)

    index = (vcf > 0) & (sample_size > 10)
    vcf_fit_params, _ = optimize.curve_fit(
        mulplicative_model, sample_size[index][:-1], vcf[index][:-1])
    vcf_fit = mulplicative_model(sample_size, *vcf_fit_params)
    vcfz_fit_params, _ = optimize.curve_fit(
        mulplicative_model, sample_size[index][:-1], vcfz[index][:-1])
    vcfz_fit = mulplicative_model(sample_size, *vcfz_fit_params)

    df = pd.DataFrame({
        "sample_size": sample_size,
        "compressed": compressed,
        "uncompressed": uncompressed,
        "vcf": vcf,
        "vcfz": vcfz,
        "vcf_fit": vcf_fit,
        "vcfz_fit": vcfz_fit,
        "tsk_fit": tsk_fit,
        "tskz_fit": tskz_fit})
    df.to_csv(datafile)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run the plot showing the files size storing everyone")
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    subparser = subparsers.add_parser("simulate")
    subparser.set_defaults(func=run_simulate)

    subparser = subparsers.add_parser("convert-files")
    subparser.set_defaults(func=run_convert_files)

    subparser = subparsers.add_parser("make-data")
    subparser.set_defaults(func=run_make_data)

    subparser = subparsers.add_parser("benchmark")
    subparser.set_defaults(func=run_benchmark)

    args = parser.parse_args()
    args.func()
