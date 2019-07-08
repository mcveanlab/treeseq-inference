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

import tskit
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

    ts = tskit.load(filename)
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

    ts = tskit.load(filename)
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


def process_hg01933_local_gnn():
    filename = os.path.join(data_prefix, "1kg_chr20.snipped.trees")
    ts = tskit.load(filename)
    region_sample_set_map = collections.defaultdict(list)
    for population in ts.populations():
        md = json.loads(population.metadata.decode())
        region = md["super_population"]
        region_sample_set_map[region].extend(list(ts.samples(
            population=population.id)))
    regions = list(region_sample_set_map.keys())
    region_sample_sets = [region_sample_set_map[k] for k in regions]

    def local_gnn(ts, focal, reference_sets):
        reference_set_map = np.zeros(ts.num_nodes, dtype=int) - 1
        for k, reference_set in enumerate(reference_sets):
            for u in reference_set:
                if reference_set_map[u] != -1:
                    raise ValueError("Duplicate value in reference sets")
                reference_set_map[u] = k

        K = len(reference_sets)
        A = np.zeros((len(focal), ts.num_trees, K))
        lefts = np.zeros(ts.num_trees, dtype=float)
        rights = np.zeros(ts.num_trees, dtype=float)
        parent = np.zeros(ts.num_nodes, dtype=int) - 1
        sample_count = np.zeros((ts.num_nodes, K), dtype=int)

        # Set the intitial conditions.
        for j in range(K):
            sample_count[reference_sets[j], j] = 1

        for t, ((left, right),edges_out, edges_in) in enumerate(ts.edge_diffs()):
            for edge in edges_out:
                parent[edge.child] = -1
                v = edge.parent
                while v != -1:
                    sample_count[v] -= sample_count[edge.child]
                    v = parent[v]
            for edge in edges_in:
                parent[edge.child] = edge.parent
                v = edge.parent
                while v != -1:
                    sample_count[v] += sample_count[edge.child]
                    v = parent[v]

            # Process this tree.
            for j, u in enumerate(focal):
                focal_reference_set = reference_set_map[u]
                p = parent[u]
                lefts[t] = left
                rights[t] = right
                while p != tskit.NULL:
                    total = np.sum(sample_count[p])
                    if total > 1:
                        break
                    p = parent[p]
                if p != tskit.NULL:
                    scale = 1 / (total - int(focal_reference_set != -1))
                    for k, reference_set in enumerate(reference_sets):
                        n = sample_count[p, k] - int(focal_reference_set == k)
                        A[j, t, k] = n * scale
        return (A, lefts, rights)

    for ind in ts.individuals():
        md = json.loads(ind.metadata.decode())
        if md["individual_id"] == "HG01933":
            for j, node in enumerate(ind.nodes):
                A, left, right = local_gnn(ts, [node], region_sample_sets)
                df = pd.DataFrame(data=A[0], columns=regions)
                df["left"] = left
                df["right"] = right
                # Remove rows with no difference in GNN to next row
                keep_rows = ~(df.iloc[:, 0:5].diff(axis=0) == 0).all(axis=1)
                df = df[keep_rows]
                df.to_csv("data/HG01933_local_gnn_{}.csv".format(j))


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
    ts = tskit.load(filename)

    # construct the dictionary mapping individual names to their metadata
    tables = ts.tables
    individual_name_map = {}
    for individual in ts.individuals():    
        metadata = json.loads(individual.metadata.decode())
        name = metadata["individual_id"]
        individual_name_map[name] = individual

    # construct a dictionary linking individual's names to their number of 
    # breakpoints within 100bp of each other
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
            close_100 += len(left_1[(left_1 >= breakpoint - 100) & (left_1 <= breakpoint + 100)])
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

def process_ukbb_1kg_duplicates():
    source_file = os.path.join(data_prefix, "1kg_ukbb_chr20.nosimplify.trees")
    ts = tskit.load(source_file)
    print("loaded")
    tables = ts.tables
    child_counts = np.bincount(tables.edges.child)
    samples = ts.samples()
    all_counts = child_counts[ts.samples()]

    # First find all samples with < 50 edges.
    candidates = samples[np.where(all_counts < 50)]

    # Now pull out those that are consecutive (i.e. from the same individual)
    c1 = candidates[np.where(candidates[:-1] == (candidates[1:] - 1))]

    print("Found", c1.shape[0], "matches")
     
    edge_child = tables.edges.child
    edge_parent = tables.edges.parent
    edge_left = tables.edges.left
    edge_right = tables.edges.right

    for c in c1:
        node = ts.node(c)
        ind = ts.individual(node.individual)
        print("Individual", ind.id, ind.nodes)
        for u in [c, c + 1]:
            print("\tChild", u)
            index = edge_child == u
            p = edge_parent[index]
            length = edge_right[index] - edge_left[index]
            # Compute the total length covered by all the parents
            counter = collections.Counter()
            for parent, L in zip(p, length):
                counter[parent] += L / ts.sequence_length

            for parent, L in counter.most_common(2):
                node = ts.node(parent)
                print("\t\tparent=", parent, "ind=", node.individual, "len=", L)
                if L > 0.9:
                    break

def main():
    name_map = {
        "sample_edges": process_sample_edges,
        "sample_edge_outliers": process_sample_edge_outliers,
        "1kg_ukbb_gnn": process_1kg_ukbb_gnn,
        "ukbb_ukbb_gnn": process_ukbb_ukbb_gnn,
        "hg01933_local_gnn": process_hg01933_local_gnn,
        "ukbb_1kg_duplicates": process_ukbb_1kg_duplicates,
    }

    parser = argparse.ArgumentParser(
        description="Process the human data and make data files for plotting.")
    parser.add_argument(
        "name", type=str, help="figure name", choices=list(name_map.keys()))

    args = parser.parse_args()
    name_map[args.name]()


if __name__ == "__main__":
    main()
