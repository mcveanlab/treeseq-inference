"""
Script to produce the data for plot showing accuracy of frequency as a proxy
for allele age (supplementary figure s1).
"""

import msprime

import pandas as pd
import random

random_seed = 13

results_total = {}
results_agree = {}
results_total_error = {}
results_agree_error = {}

for i in range(0,100):
    results_total[i] = 0
    results_agree[i] = 0
    results_total_error[i] = 0
    results_agree_error[i] = 0

def make_errors_genotype_model(g, error_matrix):
    """
    Given an empirically estimated error matrix, resample for a particular
    variant. Given a true genotype of g0, g1, or g2, return observed genotype
    depending on the error_matrix. For example, given a variant with genotype
    g0, return g0 with probability 99.942%, g1 with probability 0.041%, and
    g2 with probability 0.017%. Treat each pair of alleles as a diploid 
    individual. 
    """

    w = np.copy(g)
    #Make diploid (iterate each pair of alleles)
    genos=[(w[i],w[i+1]) for i in range(0,w.shape[0],2)]
    #Record the true genotypes
    g0 = [i for i, x in enumerate(genos) if x == (0,0)]
    g1a = [i for i, x in enumerate(genos) if x == (1,0)]
    g1b = [i for i, x in enumerate(genos) if x == (0,1)]
    g2 = [i for i, x in enumerate(genos) if x == (1,1)]

    for idx in g0:
        result=random.choices(
            [(0,0),(1,0),(1,1)], weights=error_matrix[['p00','p01','p02']].values[0])
        if result == (1,0):
            genos[idx]=random.choices([(0,1),(1,0)])
        else:
            genos[idx] = result

    for idx in g1a:
        genos[idx]=random.choices(
            [(0,0),(1,0),(1,1)], weights=error_matrix[['p10','p11','p12']].values[0])

    for idx in g1b:
        genos[idx]=random.choices(
            [(0,0),(0,1),(1,1)], weights=error_matrix[['p10','p11','p12']].values[0])

    for idx in g2:
        result=random.choices(
            [(0,0),(1,0),(1,1)], weights=error_matrix[['p20','p21','p22']].values[0])

        if result == (1,0):
            genos[idx]=random.choices([(0,1),(1,0)])

        else:
            genos[idx] = result
        
    return(np.array(sum([chromo for tup in genos for chromo in tup], ())))
    
    

def generate_samples_empirical(ts, filename, error_matrix):

    """
    Generate a samples file from a simulated ts based on an empirically estimated error matrix
    Reject any variants that result in a fixed column. 
    """

    error_matrix=pd.read_csv(error_matrix)
    record_rate = logging.getLogger().isEnabledFor(logging.INFO)
    n_variants = bits_flipped = 0
    assert ts.num_sites != 0
    sample_data = tsinfer.SampleData(
        path="../data/"+filename + "/"+filename+".samples",
        sequence_length=ts.sequence_length)
    for v in ts.variants():
        n_variants += 1
        #Record the allele frequency
        m = v.genotypes.shape[0]
        frequency = np.sum(v.genotypes) / m

        #Find closest row in error matrix file
        closest_freq = error_matrix.iloc[(error_matrix['freq']-frequency).abs().argsort()[:1]]

        #make new genotypes with error
        # Reject any columns that have no 1s or no zeros.
        # Unless the original also has them, as occasionally we have
        # some sims (e.g. under selection) where a variant is fixed
        while True:
            genotypes = make_errors_genotype_model(v.genotypes,closest_freq)
            s = np.sum(genotypes)
            if 0 < s < ts.sample_size:
                break
            if s == np.sum(v.genotypes):
                break
        if record_rate:
            bits_flipped += np.sum(np.logical_xor(genotypes, v.genotypes))
        sample_data.add_site(
            position=v.site.position, alleles=v.alleles,
            genotypes=genotypes)
    logging.info(": actual error rate = {}".format(
            bits_flipped/(n_variants*ts.sample_size)) if record_rate else "")
    sample_data.finalise()
    return sample_data

def generate_samples(ts, filename, real_error_rate=0,no_error=False):
    """
    Generate a samples file from a simulated ts
    Samples may have bits flipped with a specified probability.
    (reject any variants that result in a fixed column)
    """

    assert ts.num_sites != 0
    
    sample_data = tsinfer.SampleData(path="../data/"+filename + "/"+filename+".samples", sequence_length=ts.sequence_length)
    for v in ts.variants():
        sample_data.add_site(
            position=v.site.position, alleles=v.alleles,
            genotypes=v.genotypes)

    sample_data.finalise()
    return sample_data


def run_comparisons(sample_size,Ne,length,rec_rate,mut_rate,seed):
    simulation = msprime.simulate(sample_size=sample_size, Ne=Ne, length=length, recombination_rate=rec_rate, mutation_rate=mut_rate, random_seed=seed)
    sample_data = generate_samples(simulation,"samples_file",0,True)
    sample_data_error = generate_samples_empirical(simulation,"samples_file","/home/wilderwohns/tsinfer_paper/data/result.profile.1000g.platinum.empirical.csv")
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
        
        distance = abs(mutation_positions[combo[0]] - mutation_positions[combo[1]])
        
        if (age1[1] != age2[1]):
            if (freq1 != 1/sample_size) and (freq2 != 1/sample_size) and (freq1 != (1-sample_size)/sample_size) and (freq2 != (1-sample_size)/sample_size):
                results_total[distance//1000] += 1
                #results.loc[distance//1000,"Total"] += 1
                if freq1 < 0.02:
                    print(freq1)
                if (freq1 != freq2):
                    if ((age1[1] < age2[1]) == (freq1 < freq2)):
                        results_agree[distance//1000] += 1
                        #results.loc[distance//1000,"Correct"] += 1

                else:
                    if (random.choice([True, False])):
                        results_agree[distance//1000] += 1
                        #results.loc[distance//1000,"Correct"] += 1
                         
            if (freq1_error != 1/sample_size) and (freq2_error != 1/sample_size) and (freq1_error != (1-sample_size)/sample_size) and (freq2_error != (1-sample_size)/sample_size):  
                results_total_error[distance//1000] += 1
                #results.loc[distance//1000,"Total"] += 1
                if (freq1_error != freq2_error):
                    if ((age1[1] < age2[1]) == (freq1_error < freq2_error)):
                        results_agree_error[distance//1000] += 1
                        #results.loc[distance//1000,"Correct"] += 1

                else:
                    if (random.choice([True, False])):
                        results_agree_error[distance//1000] += 1
                        #results.loc[distance//1000,"Correct"] += 1







if __name__ == "__main__":
	for i in tqdm(range(0,1000)): 
	    cur_seed = random.seed(a=random_seed)
	    run_comparisons(50,10000,100000,2e-8,2e-8,cur_seed)
	   


