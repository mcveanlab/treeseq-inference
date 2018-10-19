"""
Script to produce the data for plot showing accuracy of frequency as a proxy
for allele age (supplementary figure s1).
"""
import random
import itertools
import logging
import argparse
import multiprocessing

import pandas as pd
import numpy as np
import tqdm

import msprime
import tsinfer


error_matrix=pd.read_csv("data/EmpiricalErrorPlatinum1000G.csv")

def make_errors_genotype_model(g, error_probs):
    """
    Given an empirically estimated error probability matrix, resample for a particular
    variant. Determine variant frequency and true genotype (g0, g1, or g2),
    then return observed genotype based on row in error_probs with nearest
    frequency. Treat each pair of alleles as a diploid individual.
    """
    w = np.copy(g)

    # Make diploid (iterate each pair of alleles)
    genos = [(w[i], w[i+1]) for i in range(0, w.shape[0], 2)]

    # Record the true genotypes
    g0 = [i for i, x in enumerate(genos) if x == (0, 0)]
    g1a = [i for i, x in enumerate(genos) if x == (1, 0)]
    g1b = [i for i, x in enumerate(genos) if x == (0, 1)]
    g2 = [i for i, x in enumerate(genos) if x == (1, 1)]

    for idx in g0:
        result = [(0, 0), (1, 0), (1, 1)][
            np.random.choice(3, p=error_probs[['p00', 'p01', 'p02']].values[0])]
        if result == (1, 0):
            genos[idx] = [(0, 1), (1, 0)][np.random.choice(2)]
        else:
            genos[idx] = result
    for idx in g1a:
        genos[idx] = [(0, 0), (1, 0), (1, 1)][
            np.random.choice(3, p=error_probs[['p10', 'p11', 'p12']].values[0])]
    for idx in g1b:
        genos[idx] = [(0, 0), (0, 1), (1, 1)][
            np.random.choice(3, p=error_probs[['p10', 'p11', 'p12']].values[0])]
    for idx in g2:
        result = [(0, 0), (1, 0), (1, 1)][
            np.random.choice(3, p=error_probs[['p20', 'p21', 'p22']].values[0])]
        if result == (1, 0):
            genos[idx] = [(0, 1), (1, 0)][np.random.choice(2)]
        else:
            genos[idx] = result

    return np.array(sum(genos, ()))    
    

def generate_samples_empirical(ts):

    """
    Generate a samples file from a simulated ts based on an empirically estimated error matrix
    Reject any variants that result in a fixed column. 
    """
    assert ts.num_sites != 0
    sample_data = tsinfer.SampleData(sequence_length=ts.sequence_length)

    for v in ts.variants():
        #Record the allele frequency
        m = v.genotypes.shape[0]
        frequency = np.sum(v.genotypes) / m

        #Find closest row in error matrix file
        closest_freq = error_matrix.iloc[(error_matrix['freq']-frequency).abs().argsort()[:1]]

        #make new genotypes with error
        # Reject any columns that have no 1s or no zeros.
        # Unless the original also has them, as occasionally we have
        # some sims (e.g. under selection) where a variant is fixed
        genotypes = make_errors_genotype_model(v.genotypes,closest_freq)
        
        sample_data.add_site(
            position=v.site.position, alleles=v.alleles,
            genotypes=genotypes)

    sample_data.finalise()
    return sample_data

def generate_samples(ts):
    """
    Generate a samples file from a simulated ts
    Samples may have bits flipped with a specified probability.
    (reject any variants that result in a fixed column)
    """

    assert ts.num_sites != 0
    
    sample_data = tsinfer.SampleData(sequence_length=ts.sequence_length)
    for v in ts.variants():
        sample_data.add_site(
            position=v.site.position, alleles=v.alleles,
            genotypes=v.genotypes)

    sample_data.finalise()
    return sample_data


def run_comparisons(params):
    sample_size, Ne, length, rec_rate, mut_rate, seed = params
    n_bins = 20
    bin_size = length//n_bins
    assert bin_size == length/n_bins
    total = np.zeros(n_bins, dtype=np.int)
    agree = np.zeros(n_bins, dtype=np.int)
    total_error = np.zeros(n_bins, dtype=np.int)
    agree_error = np.zeros(n_bins, dtype=np.int)

    simulation = msprime.simulate(
        sample_size=sample_size, Ne=Ne, length=length, recombination_rate=rec_rate, 
        mutation_rate=mut_rate, random_seed=seed)
    sample_data = generate_samples(simulation)
    sample_data_error = generate_samples_empirical(simulation)
    num_mutations = len(sample_data.sites_genotypes[:])
    variant_ages = [(mutation.id,simulation.node(mutation.node).time) for mutation in list(simulation.mutations())]
    variant_frequencies = [len(geno[geno==1])/len(geno) for geno in sample_data.sites_genotypes[:]]
    variant_frequencies_error = [len(geno[geno==1])/len(geno) for geno in sample_data_error.sites_genotypes[:]]

    mutation_positions=[mutation.position for mutation in list(simulation.mutations())]

    for combo in list(itertools.combinations(range(0,num_mutations),2)):
        age1 = variant_ages[combo[0]]
        age2 = variant_ages[combo[1]]

        freq1 = variant_frequencies[combo[0]]
        freq2 = variant_frequencies[combo[1]]

        freq1_error = variant_frequencies_error[combo[0]]
        freq2_error = variant_frequencies_error[combo[1]]
        
        distance = int(abs(mutation_positions[combo[0]] - mutation_positions[combo[1]]))
        
        if (age1[1] != age2[1]):
            total[distance//bin_size] += 1

            if (freq1 != freq2):
                if ((age1[1] < age2[1]) == (freq1 < freq2)):
                    agree[distance//bin_size] += 1

            else:
                if (random.choice([True, False])):
                    agree[distance//bin_size] += 1

            if (freq1_error != freq2_error):
                if ((age1[1] < age2[1]) == (freq1_error < freq2_error)):
                    agree_error[distance//bin_size] += 1

            else:
                if (random.choice([True, False])):
                    agree_error[distance//bin_size] += 1

    return total, agree, agree_error



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
         '--replicates', '-r', type=int, default=4000, help="number of replicates")
    parser.add_argument(
         '--seed', '-s', type=int, default=123, help="use a non-default RNG seed")
    parser.add_argument(
        "--processes", '-p', type=int, default=1,
        help="number of worker processes, e.g. 40")
    parser.add_argument(
         '--progress',  "-P", action='store_true',
         help="Show a progress bar.", )

    args = parser.parse_args()
    
    np.random.seed(args.seed)
    random_seeds = np.random.randint(1 , 2**32, size = args.replicates)
    run_params = zip(
        itertools.repeat(50),     # sample_size
        itertools.repeat(5000),   # Ne 
        itertools.repeat(100000), # length
        itertools.repeat(1e-8), # recomb_rate
        itertools.repeat(1e-8), # mut_rate
        random_seeds # seed
        )
    
    df = pd.DataFrame.from_dict({"Agree":np.zeros(20), "Total":np.zeros(20), "ErrorAgree":np.zeros(20)})
    
    if args.processes > 1:
        logging.info("Setting up using multiprocessing ({} processes)".format(args.processes))
        with multiprocessing.Pool(processes=args.processes, maxtasksperchild=2) as pool:
            for total, agree, agree_error in tqdm.tqdm(
                pool.imap_unordered(run_comparisons, run_params), total=args.replicates, disable=not args.progress):
                df.Total += total
                df.Agree += agree
                df.ErrorAgree += agree_error
                
    else:
        # When we have only one process it's easier to keep everything in the same
        # process for debugging.
        logging.info("Setting up using a single process")
        for total, agree, agree_error in tqdm.tqdm(
            map(run_comparisons, run_params), total=args.replicates, disable=not args.progress):
            df.Total += total
            df.Agree += agree
            df.ErrorAgree += agree_error

    df.to_csv("data/frequency_distance_accuracy_singletons.csv")



