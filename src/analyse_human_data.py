"""
Scripts used to analyse data in the human_data directory and produce the
data files used for plotting.
"""
import argparse
import os.path
import json
import collections

import numpy as np
import pandas as pd

import sys
# FIXME!!! mean_descendants not yet merged to msprime.
sys.path.insert(0, "/gpfs1/well/mcvean/ukbb12788/jk/msprime")

import msprime
import tqdm

data_prefix = "human-data"

def print_sample_edge_stats(ts):
    """ 
    Print out some basic stats about the sample edges in the specified tree 
    sequence.
    """
    tables = ts.tables
    child_counts = np.bincount(tables.edges.child)
    samples = ts.samples()
    all_counts = child_counts[samples]
    print("mean sample count    = ", np.mean(all_counts))

    index = tables.nodes.flags[tables.edges.child] == msprime.NODE_IS_SAMPLE
    length = tables.edges.right[index]- tables.edges.left[index]
    print("mean length          = ", np.mean(length))
    print("median length        = ", np.median(length))

    n50 = np.zeros(ts.num_samples)
    edges = tables.edges
    child = edges.child
    left = edges.left
    right = edges.right
    for j, sample in enumerate(samples):
        index = child == sample
        length = right[index]- left[index]
        length = np.sort(length)[::-1]
        cumulative = np.cumsum(length)
        # N50 is the first length such that the cumulative sum up to that point
        # is >= L / 2.
        n50[j] = length[cumulative >= ts.sequence_length / 2][0]
    print("Average N50          = ", np.mean(n50))


def get_sgdp_sample_edges():
    filename = os.path.join(data_prefix, "sgdp_chr20.nosimplify.trees")

    ts = msprime.load(filename)
    print("SGDP")
    print_sample_edge_stats(ts)
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
    print("TGP")
    print_sample_edge_stats(ts)
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


def process_hg01933_parent_ancestry():
    filename = os.path.join(data_prefix, "1kg_chr20.snipped.trees")
    ts = msprime.load(filename)
    tables = ts.tables
    region_sample_set_map = collections.defaultdict(list)
    for population in ts.populations():
        md = json.loads(population.metadata.decode())
        region = md["super_population"]
        region_sample_set_map[region].extend(list(ts.samples(population=population.id)))
    regions = list(region_sample_set_map.keys())
    region_sample_sets = [region_sample_set_map[k] for k in regions]

    D = ts.mean_descendants(region_sample_sets)   

    def parent_gnn(sample_id):
        index = tables.edges.child == sample_id
        left = tables.edges.left[index]
        right = tables.edges.right[index]
        parent = tables.edges.parent[index]
        k = parent.shape[0]
        length = (right - left).reshape((k, 1))
        D_parent = D[parent]
        total = np.sum(D_parent, axis=1).reshape(k, 1)
        P_gnn = D_parent / total
        df = pd.DataFrame({region: P_gnn[:,j] for j, region in enumerate(regions)})
        df["left"] = left
        df["right"] = right
        df = df.sort_values("left")
        return df.set_index("left")

    for ind in ts.individuals():
        md = json.loads(ind.metadata.decode())
        if md["individual_id"] == "HG01933":
            for j, node in enumerate(ind.nodes):
                df = parent_gnn(node)
                df.to_csv("data/HG01933_parent_ancestry_{}.csv".format(j))


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


def process_sample_edge_outliers():
    """
    Runs the analysis for finding the sample edge outliers.
    """
    filename = os.path.join(data_prefix, "1kg_chr20.nosimplify.trees")
    ts = msprime.load(filename)

    # construct the dictionary mapping individual names to their metadata
    tables = ts.tables
    individual_name_map = {}
    for individual in ts.individuals():    
        metadata = json.loads(individual.metadata.decode())
        name = metadata["individual_id"]
        individual_name_map[name] = individual

    #construct a dictionary linking individual's names to their number of breakpoints within 100bp of each other
    close_breakpoints = dict()
    child = tables.edges.child
    left = tables.edges.left

    for key, individual in tqdm.tqdm(individual_name_map.items()):
        index_0 = child == individual.nodes[0]
        left_0 = left[index_0]
        index_1 = child == individual.nodes[1]
        left_1 = left[index_1]
        close_100 = 0
        for breakpoint in left_0:
            close_100 += len(left_1[(left_1 >= breakpoint-100) & (left_1 <= breakpoint+100)])
        close_breakpoints[key] = close_100

    print("Average = ", np.mean(list(close_breakpoints.values())))
    for ind in ["NA20289", "HG02789"]:
        print(ind, ":", close_breakpoints[ind])



def process_1kg_ukbb_gnn():
    source_file = os.path.join(data_prefix, "1kg_ukbb_chr20.snipped.trees.gnn.csv")
    df = pd.read_csv(source_file)

    # Use TGP populations here to make sure we don't leak any metadata.
    tgp_populations = [
        'CHB', 'JPT', 'CHS', 'CDX', 'KHV',
        'CEU', 'TSI', 'FIN', 'GBR', 'IBS',
        'YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB',
        'MXL', 'PUR', 'CLM', 'PEL',
        'GIH', 'PJL', 'BEB', 'STU', 'ITU']
    # Overall GNN by 1KG population
    dfg = df.groupby(df.ethnicity).mean()
    dfg = dfg[tgp_populations]
    datafile = "data/1kg_ukbb_ethnicity.csv"
    dfg.to_csv(datafile)

    # Subset down to the british ethnicity.
    df = df[df.ethnicity == "British"]
    print("British subset = ", len(df))
    dfg = df.groupby(df.centre).mean()
    dfg = dfg[tgp_populations]
    datafile = "data/1kg_ukbb_british_centre.csv"
    dfg.to_csv(datafile)

def process_ukbb_ukbb_gnn():
    source_file = os.path.join(data_prefix, "ukbb_chr20.augmented_131072.snipped.trees.gnn.csv")
    df = pd.read_csv(source_file)
    
    # Subset down to the british ethnicity.
    df = df[df.ethnicity == "British"]
    print("British subset = ", len(df))
    dfg = df.groupby(df.centre).mean()
    dfg = dfg[df.centre.unique()]
    datafile = "data/ukbb_ukbb_british_centre.csv"
    dfg.to_csv(datafile)

def main():
    name_map = {
        "sample_edges": process_sample_edges,
        "sample_edge_outliers": process_sample_edge_outliers,
        "1kg_ukbb_gnn": process_1kg_ukbb_gnn,
        "ukbb_ukbb_gnn": process_ukbb_ukbb_gnn,
        "hg01933_parent_ancestry": process_hg01933_parent_ancestry,
    }

    parser = argparse.ArgumentParser(
        description="Process the human data and make data files for plotting.")
    parser.add_argument(
        "name", type=str, help="figure name", choices=list(name_map.keys()))

    args = parser.parse_args()
    name_map[args.name]()


if __name__ == "__main__":
    main()
