"""
Simulations for the Out-of-Africa simulations of the sample edges 
statistic.
"""
import math

import msprime
import numpy as np
import pandas as pd


def out_of_africa():
    """
    Exactly the same as https://msprime.readthedocs.io/en/stable/tutorial.html 
    except with different sample sizes
    """
    # Based on 1000G sample sizes.
    AF_sample_size = 1322
    AS_sample_size = 1986
    EU_sample_size = 1006

    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
 
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    return dict(
        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=AF_sample_size, initial_size=N_AF),
            msprime.PopulationConfiguration(
                sample_size=EU_sample_size, initial_size=N_EU, growth_rate=r_EU),
            msprime.PopulationConfiguration(
                sample_size=AS_sample_size, initial_size=N_AS, growth_rate=r_AS)
        ],
        migration_matrix = [
            [      0, m_AF_EU, m_AF_AS],
            [m_AF_EU,       0, m_EU_AS],
            [m_AF_AS, m_EU_AS,       0],
        ],
        demographic_events = [
            # CEU and CHB merge into B with rate changes at T_EU_AS
            msprime.MassMigration(
                time=T_EU_AS, source=2, destination=1, proportion=1.0),
            msprime.MigrationRateChange(time=T_EU_AS, rate=0),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
            # Population B merges into YRI at T_B
            msprime.MassMigration(
                time=T_B, source=1, destination=0, proportion=1.0),
            # Size changes to N_A at T_AF
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]
    )

def main():
    # simulate the length of chromosome 20 using out of africa model
    length = 64444167
    ts = msprime.simulate(
        **out_of_africa(), length=length, recombination_rate=2e-8, mutation_rate=2e-8,
        random_seed=12345)

    population_name = ["African", "European", "Asian"]
    tables = ts.tables
    child_counts = np.bincount(tables.edges.child) 

    sample_edges = []
    population = []
    for u in ts.samples():
        population.append(population_name[ts.node(u).population])
        sample_edges.append(child_counts[u])

    df = pd.DataFrame({"population": population, "sample_edges": sample_edges})
    df.to_csv("data/ooa_sample_edges_sim.csv")


if __name__ == "__main__":
    main()
