"""
Scripts used to analyse data in the human_data directory and produce the
data files used for plotting.
"""
import argparse
import os.path
import json

import numpy as np
import pandas as pd
import msprime

data_prefix = "human-data"


def get_sgdp_sample_edges():
    filename = os.path.join(data_prefix, "sgdp_chr20.nosimplify.trees")

    ts = msprime.load(filename)
    population_name = []
    population_region = []
    for pop in ts.populations():
        md = json.loads(pop.metadata.decode())
        population_name.append(md["name"])
        population_region.append(md["region"])

    tables = ts.tables
    child_counts = np.bincount(tables.edges.child)
    datasets = []
    samples = []
    strands = []
    populations = []
    regions = []
    sample_edges = []
    for ind in ts.individuals():
        md = json.loads(ind.metadata.decode())
        for j, node_id in enumerate(ind.nodes):
            node = ts.node(node_id)
            samples.append(md["sgdp_id"])
            strands.append(j)
            populations.append(population_name[node.population])
            regions.append(population_region[node.population])
            sample_edges.append(child_counts[node_id])
            datasets.append("sgdp")

    df = pd.DataFrame({
        "dataset": datasets,
        "sample": samples,
        "strand": strands,
        "population": populations,
        "region": regions,
        "sample_edges": sample_edges})
    return df


def get_1kg_sample_edges():
    filename = os.path.join(data_prefix, "1kg_chr20.nosimplify.trees")

    ts = msprime.load(filename)
    population_name = []
    population_region = []
    for pop in ts.populations():
        md = json.loads(pop.metadata.decode())
        population_name.append(md["name"])
        population_region.append(md["super_population"])

    tables = ts.tables
    child_counts = np.bincount(tables.edges.child)
    datasets = []
    samples = []
    strands = []
    populations = []
    regions = []
    sample_edges = []
    for ind in ts.individuals():
        md = json.loads(ind.metadata.decode())
        for j, node_id in enumerate(ind.nodes):
            node = ts.node(node_id)
            samples.append(md["individual_id"])
            strands.append(j)
            populations.append(population_name[node.population])
            regions.append(population_region[node.population])
            sample_edges.append(child_counts[node_id])
            datasets.append("1kg")

    df = pd.DataFrame({
        "dataset": datasets,
        "sample": samples,
        "strand": strands,
        "population": populations,
        "region": regions,
        "sample_edges": sample_edges})
    return df


def process_sample_edges():
    """
    Processes data from the SGDP and 1KG data files to produce a data frame
    containing the number of sample edges for every sample.
    """
    df_sgdp = get_sgdp_sample_edges()
    df_1kg = get_1kg_sample_edges()
    df_all = pd.concat([df_sgdp, df_1kg], ignore_index=True)
    datafile = "data/sample_edges.csv"
    df_all.to_csv(datafile)


def main():
    name_map = {
        "sample_edges": process_sample_edges
    }

    parser = argparse.ArgumentParser(
        description="Process the human data and make data files for plotting.")
    parser.add_argument(
        "name", type=str, help="figure name", choices=list(name_map.keys()))

    args = parser.parse_args()
    name_map[args.name]()


if __name__ == "__main__":
    main()
