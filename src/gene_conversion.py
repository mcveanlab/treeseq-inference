import subprocess
import random
import string
import logging
import argparse

import numpy as np

import msprime, pyslim

eidos_cmd = string.Template("""

initialize() {
    $set_random_seed_cmd
	initializeTreeSeq(checkCoalescence=T);
	initializeMutationRate(0);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, $length-1);
	initializeRecombinationRate($recombination_rate);
    initializeGeneConversion($gene_conversion_proportion, $gene_conversion_mean_length);
}

1 late() {
	sim.addSubpop("p0", $popsize);
}

1: late() {
if (sim.treeSeqCoalesced()) {
	catn("Simulation finished in " + sim.generation + " generations");
	sim.treeSeqOutput("$file_prefix" + ".trees");
	sim.simulationFinished();}
}

$max_generations {
	catn("NO COALESCENCE BY GENERATION $max_generations");
}
""")

def mutate_simplify(ts, mu, samples, seed):
    subsamples = ts.samples()[samples]
    return msprime.mutate(ts, mu, random_seed=seed, keep=True).simplify(subsamples)


def WF_GC_simulation(
    popsize, chrom_length, recomb_rate, mut_rate, gc_rate, gc_mean_len, nsamples,
    file_prefix = "gc_sim", max_generations=1e9, seed=None, slimname="slim"):
    """
    Carry out a simulation of a wright-fisher population with gene conversion
    
    nsamples is number of haploid (genome) samples
    """
    full_pop_file_prefix = file_prefix + "_all"
    eidos_seed_cmd = "" if seed is None else "setSeed({});".format(int(seed))
    cmd = eidos_cmd.substitute(
        set_random_seed_cmd = eidos_seed_cmd,
        file_prefix = full_pop_file_prefix,
        popsize = popsize,
        length = chrom_length,
        recombination_rate = recomb_rate + gc_rate,
        gene_conversion_proportion = gc_rate / (recomb_rate + gc_rate),
        gene_conversion_mean_length = gc_mean_len,
        max_generations = int(max_generations),
    )

    logging.info("Simulating")
    process = subprocess.Popen(slimname,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
        universal_newlines=True)
    process.stdin.write(cmd)
    process.stdin.close()
    suppress_header_line = 3
    for line in iter(process.stdout.readline, ''):
        print(line)
        if suppress_header_line:
            if suppress_header_line != 3:
                suppress_header_line -= 1 #decrement
            if line.startswith("// Starting run at generation"):
                suppress_header_line = 2 #this is the penultimate line of the header
        else:
            logging.info(line.rstrip())
    ts = pyslim.load(full_pop_file_prefix + ".trees")
    samp = np.random.choice(ts.num_samples, nsamples, replace=False)
    ts = mutate_simplify(ts, mut_rate, samp, seed)
    logging.info("Finished mutating")
    return ts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run a simulation with gene conversion.")
    parser.add_argument(
        "output",
        help="The path to write the output file to")
    parser.add_argument(
        "--popsize", "-N", type=int, default=5000, 
        help="Ne")
    parser.add_argument(
        "--seed", "-s", type=int, default=None, 
        help="run a single simulation with this seed")
    parser.add_argument(
        "--slim_binary", "-S", type=str, default="slim", 
        help="Path to the slim binary")
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    ts = WF_GC_simulation(
        popsize=args.popsize, chrom_length=100000, recomb_rate=1e-8, mut_rate=1e-8, gc_rate=1e-7,
        gc_mean_len = 500, nsamples=16, file_prefix = args.output, seed = args.seed, slimname=args.slim_binary)
    ts.dump(args.output + ".trees")