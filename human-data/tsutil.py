"""
Various utilities for manipulating tree sequences and running tsinfer.
"""
import argparse
import subprocess
import time
import collections
import json
import sys

import msprime
import tsinfer
import daiquiri
import numpy as np
import tqdm


def run_simplify(args):
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
        final_file = base + ".augmented_{}.nosimplify.trees".format(n)
        print("RUNNING", augmented_file)
        subset = np.linspace(offset, num_samples - offset, n + 2, dtype=int)[1: -1]
        ancestors_ts = run_augment(sample_data, ancestors_ts, subset, args.num_threads)
        ancestors_ts.dump(augmented_file)
        n *= 2
        offset += 1

    final_ts = run_match_samples(sample_data, ancestors_ts, args.num_threads)
    final_ts.dump(final_file)

def run_benchmark(args):

    before = time.perf_counter()
    ts = msprime.load(args.input)
    duration = time.perf_counter() - before
    print("Loaded in {:.2f}s".format(duration))

    before = time.perf_counter()
    j = 0
    for tree in ts.trees(sample_counts=False):
        j += 1
    assert j == ts.num_trees
    duration = time.perf_counter() - before
    print("Iterated over trees in {:.2f}s".format(duration))

    before = time.perf_counter()
    num_variants = 0
    for var in ts.variants():
        if num_variants == args.num_variants:
            break
        num_variants += 1
    duration = time.perf_counter() - before
    total_genotypes = (ts.num_samples * num_variants) / 10**6
    print("Iterated over {} variants in {:.2f}s @ {:.2f} M genotypes/s".format(
        num_variants, duration, total_genotypes / duration))


def run_combine_ukbb_1kg(args):
    ukbb_samples_file = "ukbb_{}.samples".format(args.chromosome)
    tg_ts_file = "1kg_{}.nosimplify.trees".format(args.chromosome)
    ancestors_ts_file = "1kg_ukbb_{}.ancestors.trees".format(args.chromosome)
    samples_file = "1kg_ukbb_{}.samples".format(args.chromosome)

    ukbb_samples = tsinfer.load(ukbb_samples_file)
    tg_ts = msprime.load(tg_ts_file)
    print("Loaded ts:", tg_ts.num_nodes, tg_ts.num_edges)
     
    # Subset the sites down to the UKBB sites.
    tables = tg_ts.dump_tables()
    ukbb_sites = set(ukbb_samples.sites_position[:])
    ancestors_sites = set(tables.sites.position[:])
    intersecting_sites = ancestors_sites & ukbb_sites
    # Create a map from position to alleles so we can make sure they are 
    # compatible between the two datasets.
    tg_alleles = {}
    for site in tg_ts.sites():
        if site.position in intersecting_sites:
            assert len(site.mutations) == 1
            tg_alleles[site.position] = [
                site.ancestral_state, site.mutations[0].derived_state]
    
    for pos, alleles in zip(ukbb_samples.sites_position[:], ukbb_samples.sites_alleles[:]):
        if pos in tg_alleles and alleles != tg_alleles[pos]:
            print(
                "Removing site with incompatible alleles:", pos, alleles, 
                tg_alleles[pos])
            intersecting_sites.remove(pos)

    print("Intersecting sites = ", len(intersecting_sites))
    tables.sites.clear()
    tables.mutations.clear()
    for site in tg_ts.sites():
        if site.position in intersecting_sites:
            # The ancestors tree sequence requires that sites be 0/1
            site_id = tables.sites.add_row(position=site.position, ancestral_state="0")
            assert len(site.mutations) == 1
            mutation = site.mutations[0]
            tables.mutations.add_row(site=site_id, node=mutation.node, derived_state="1")

    # Update the nodes to mark everything as a sample. Don't bother trying to keep 
    # the metadata or flags, as we can recover these from the original.
    tables.nodes.set_columns(
        flags=np.ones_like(tables.nodes.flags),
        time=tables.nodes.time + 1)
    # Reduce this to the site topology now to make things as quick as possible.
    tables.simplify(reduce_to_site_topology=True, filter_sites=False)
    tg_ancestors_ts = tables.tree_sequence()
    print("Reduced to :", tg_ancestors_ts.num_nodes, tg_ancestors_ts.num_edges)
    tg_ancestors_ts.dump(ancestors_ts_file)

    # Now create a new samples file to get rid of the missing sites.
    git_hash = subprocess.check_output(["git", "rev-parse", "HEAD"])
    git_provenance = {
        "repo": "git@github.com:mcveanlab/treeseq-inference.git",
        "hash": git_hash.decode().strip(),
        "dir": "human-data",
        "notes:": (
            "Use the Makefile to download and process the upstream data files")}
    
    with tsinfer.SampleData(
            path=samples_file, num_flush_threads=4,
            sequence_length=ukbb_samples.sequence_length) as samples:

        for _ in tqdm.tqdm(range(ukbb_samples.num_individuals)):
            samples.add_individual(ploidy=2)
        for variant in tqdm.tqdm(ukbb_samples.variants(), total=ukbb_samples.num_sites):
            if variant.site.position in intersecting_sites:
                samples.add_site(
                    position=variant.site.position, alleles=variant.alleles,
                    genotypes=variant.genotypes, metadata=variant.site.metadata)

        for timestamp, record in ukbb_samples.provenances():
            samples.add_provenance(timestamp, record)
        samples.record_provenance(
            command=sys.argv[0], args=sys.argv[1:], git=git_provenance)
        
    print(samples)

def run_compute_1kg_ancestry(args):
    ts = msprime.load(args.input)
    # This is specialised for the 1000 genomes metadata format. We could compute the 
    # ancestry for all populations, but it's a good bit slower.
    superpops = collections.defaultdict(list)
    for population in ts.populations():
        md = json.loads(population.metadata.decode())
        superpops[md["super_population"]].extend(ts.samples(population.id))
    A = tsinfer.mean_sample_ancestry(ts, list(superpops.values()), show_progress=True)
    np.save(args.output, A)


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

    subparser = subparsers.add_parser("combine-ukbb-1kg")
    subparser.add_argument("chromosome", type=str, help="chromosome stem")
    subparser.set_defaults(func=run_combine_ukbb_1kg)

    subparser = subparsers.add_parser("benchmark")
    subparser.add_argument(
        "input", type=str, help="Input tree sequence")
    subparser.add_argument(
        "--num-variants", type=int, default=None,
        help="Number of variants to benchmark genotypes decoding performance on")
    subparser.set_defaults(func=run_benchmark)

    subparser = subparsers.add_parser("compute-1kg-ancestry")
    subparser.add_argument(
        "input", type=str, help="Input tree sequence")
    subparser.add_argument(
        "output", type=str, help="Filename to write numpy array to.")
    subparser.set_defaults(func=run_compute_1kg_ancestry)

    daiquiri.setup(level="INFO")

    args = parser.parse_args()
    args.func(args)

main()
