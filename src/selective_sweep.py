#!/usr/bin/env python3.5
description = '''
Use ftprime to simulate a chromosome with a neutral mutations and a single advantageous variant swept to a 
user-specified frequency (or set of frequencies) using simuPOP. Sampled tree sequences will be
written at each user-specified point, and the variant reintroduced on a single individual if lost
from the population. The simulation ends when the condition to produce the final output file is met.
'''

import gzip
import sys
import math
import time
import random
import logging
from collections import OrderedDict

# some simupop options involving mutation type
import simuOpt
if logging.getLogger().isEnabledFor(logging.INFO):
    simuOpt.setOptions(alleleType='mutant')
else:
    simuOpt.setOptions(alleleType='mutant', quiet=True)

import simuPOP as sim
from ftprime import RecombCollector
import msprime

REPORTING_STEP = 50

def fileopt(fname,opts):
    '''Return the file referred to by fname, open with options opts;
    if fname is "-" return stdin/stdout; if fname ends with .gz run it through gzip.
    '''
    if fname == "-":
        if opts == "r":
            fobj = sys.stdin
        elif opts == "w":
            fobj = sys.stdout
        else:
            print("Something not right here.")
    elif fname[len(fname)-3:len(fname)]==".gz":
        fobj = gzip.open(fname,opts)
    else:
        fobj = open(fname,opts)
    return fobj

# Check on frequency of selected allele (run every generation)
def check_freqs_reached(pop, selected_locus, output_frequency_reached, countdown):
    """
    output_freqs is a mutable dictionary of (frequency,generations):g 
    giving the number of generations after a given frequency at which
    to output files. If g is None, we have already met this criterion
    """
    curr_freq = pop.dvars().alleleFreq[selected_locus][1]
    for output in [o for o,reached in output_frequency_reached.items() if not reached]:
        if curr_freq >= float(output[0]):
            countdown[output] = output[1] if len(output)>1 else 0 #if no "post_gen", output 0 gens after freq
            output_frequency_reached[output] = True
    return True

# output if necessary    
def output_to_msprime(pop, recomb_collector, treefile_prefix, nsamples, mut_rate, 
    output_freqs, countdown, mutations_after_simulation):
    """
    Check if any files need outputting this generation, if so, output and add names to output_freqs
    """
    for output in list(countdown.keys()):
        if countdown[output] == 0:
            assert output_freqs[output]==True
            #save msprime file as a sample from pop
            diploid_samples = random.sample(pop.indInfo("ind_id"), nsamples)
            ts = recomb_collector.tree_sequence(diploid_samples)
            
            logging.info("Loaded into tree sequence!")
            logging.info("----------")
            
            fn = treefile_prefix + output[0] + \
                ("+{}".format(output[1]) if len(output)>1 else "") + ".hdf5"
            if mutations_after_simulation:
                ts = msprime.mutate(
                    ts, mut_rate, random_seed=random.randint(1, 2**32 - 1), keep=True)
                ts.dump(fn)
            else:
                #to do
                assert False, \
                    "Injecting simuPOP mutations into the treesequence is not yet supported in ftprime"
                ts.dump(fn)
            
            logging.info("Written out samples to {} ({} variable sites)".format(fn, ts.get_num_mutations()))
            logging.info("----------")

            
            output_freqs[output] = fn #hacky: save so that we can return the file names
            del countdown[output]
        else:
            #count down until output
            countdown[output] -= 1
    return any(reached == False for reached in output_freqs.values()) or len(countdown) > 0

class FixedFitness:
    def __init__(self, s, h):
        # mean is alpha/beta
        self.s = s
        self.h = h
    def __call__(self, loc, alleles):
        # needn't return fitness for alleles=(0,0) as simupop knows that's 1
        if 0 in alleles:
            return 1. + self.h*self.s
        else:
            return 1. + self.s

def add_mutant_to_RC_ARG(simuPOP_population, ind, pos, chrom, val, recomb_collector):
    """
    Take an index into the specified population, get the unique 
    individual ID corresponding to that index, convert it to
    the ID used in the ARG stored in a recomb_collector instance,
    and mutate it to the derived value
    """
    ind_id = simuPOP_population.individual(ind).info('ind_id')
    haploid_id = recomb_collector.i2c(ind_id, chrom)
    logging.info("Gen: {:4d}. Mutation (re)introduced at locus {} in simupop individual {:g} chrom {} = ftprime id {}".format(
        simuPOP_population.dvars().gen, pos, ind_id, chrom, haploid_id))
    recomb_collector.args.add_mutation(position=pos, node=haploid_id, derived_state=b'1', ancestral_state=b'0')
    return True

def simulate_sweep(popsize, chrom_length, recomb_rate, mut_rate, selection_coef, dominance_coef, 
    output_at_freqs, nsamples, simplify_interval=500, generations_before_sweep=0, max_generations=-1, 
    mutations_after_simulation=True, treefile_prefix="sweepfile", seed=None, verbosity = 0):
    """
    Carry out a simulation of a selective sweep, and save msprime-format files at frequencies
    specified by output_at_freqs, which is a list of (frequency, post_generation) tuples
    Note that some of these files may have fixed variants (i.e. mutations above the root node)
    """
    
    if seed is not None:
        sim.setRNG(seed=seed)
        random.seed(seed)

    # locations of the loci along the chromosome?
    # hard code defaults for simupop:
    # >The default positions are 1, 2, 3, 4, ... on each
    # >chromosome.
    locus_position = list(range(0, chrom_length))
    
    # which loci are under selection?
    selected_locus = math.ceil(chrom_length / 2)
    neutral_loci = list(set(range(1,chrom_length)) - set([selected_locus]))
        
    output_countdown = {} #will contain the countdown generations before output
    output_frequency_reached = OrderedDict() # (freq,gens) tuples showing if these
    for params in output_at_freqs:
        output_frequency_reached[(params[0],int(params[1])) if len(params)>1 else (params[0],)]=False


    pop = sim.Population(
            size=popsize,
            loci=[chrom_length],
            lociPos=locus_position,
            infoFields=['ind_id','fitness'])
    
    
    # set up recomb collector
    # NB: we have to simulate an initial tree sequence, but we can put neutral mutations on
    # *after* we have run the simulation, since regardless of selective forces,
    # the probability of neutral mutations is simply proportional to branch length 
    
    id_tagger = sim.IdTagger()
    id_tagger.apply(pop)
    first_gen = pop.indInfo("ind_id")
    length = max(locus_position)
    # Since we want to have a finite site model, we force the recombination map
    # to have exactly `length` loci with a fixed recombination rate between them.
    rcmb_map = msprime.RecombinationMap.uniform_map(length, recomb_rate, length)
    if mutations_after_simulation:
        init_ts = msprime.simulate(2*len(first_gen), Ne=popsize, recombination_map=rcmb_map)
    else:
        init_ts = msprime.simulate(2*len(first_gen), Ne=popsize, recombination_map=rcmb_map, mutation_rate=mut_rate)
    
    haploid_labels = [(k,p) for k in first_gen
                            for p in (0,1)]
    node_ids = {x:j for x, j in zip(haploid_labels, init_ts.samples())}
    
    rc = RecombCollector(ts=init_ts, node_ids=node_ids,
                         locus_position=locus_position)
    
    if mutations_after_simulation:
        mutate = []
        # initially, population is monogenic
        init_geno=[sim.InitGenotype(freq=1.0)]
    else:
        #create a mutator to inject mutations in forward time
        mutate = [sim.SNPMutator(u=neut_mut_rate,v=0,loci=neutral_loci)]
        assert False, "Need to implement conversion of msprime haplotypes -> simuPOP initial haplotypes"

    if verbosity:
        logger = logging.getLogger('')
        report = [sim.PyEval(r"'Gen: %4d - report. Focal allele at freq %f' % (gen, alleleFreq[{}][1])".format(
            selected_locus), step=REPORTING_STEP, output=logger.info)]
    else:
        report = []
        
    pop.evolve(
        initOps=[
            sim.InitSex(),
        ]+init_geno,
        preOps= mutate + [
            sim.PyOperator(lambda pop: rc.increment_time() or True),
            sim.PyMlSelector(FixedFitness(selection_coef, dominance_coef),
                loci=selected_locus),
        ],
        matingScheme=sim.RandomMating(
            ops=[
                id_tagger,
                sim.Recombinator(rates=recomb_rate, output=rc.collect_recombs,
                                 infoFields="ind_id"),
            ] ),
        postOps= [
            sim.PyOperator(lambda pop: rc.simplify(pop.indInfo("ind_id")) or True,
                           step=simplify_interval),
            sim.Stat(alleleFreq=selected_locus),
            ##reintroduce if lost (http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec6.html#manually-introduced-mutations-pointmutator)
            sim.IfElse('gen > {} and alleleNum[{}][1] == 0'.format(generations_before_sweep, selected_locus), ifOps=[
                sim.PointMutator(inds=0, loci=selected_locus, allele=1),
                #add this mutation to the ARG being tracked in rc
                sim.PyOperator(lambda pop: add_mutant_to_RC_ARG(
                    pop, ind=0, pos=selected_locus, chrom=0, val=1, recomb_collector=rc))
            ]),
            ## check frequencies, and output when necessary
            sim.PyOperator(lambda pop: check_freqs_reached(
                pop, selected_locus, output_frequency_reached, output_countdown)),
            sim.PyOperator(lambda pop: output_to_msprime(
                pop, rc, treefile_prefix, nsamples, mut_rate, output_frequency_reached, output_countdown, 
                mutations_after_simulation)),
        ] + report,
        gen = max_generations
    )
    
    logging.info("Done simulating!")
    logging.info("----------")
    
    del pop
    del rc
        
    return output_frequency_reached #these should now contain file names

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description=description, add_help=False)
    parser.add_argument('--help', action='help', help='show this help message and exit')
    parser.add_argument("-N","--popsize", dest="popsize", type=int,
            help="Size of the population", default=5000)
    parser.add_argument("-r","--recomb_rate", dest="recomb_rate", type=float,
            help="Recombination rate", default=1e-7)
    parser.add_argument("-L","--length", dest="chrom_length", type=int,
            help="Number of bp in the chromosome", default=1000)
    parser.add_argument("-U","--neut_mut_rate", dest="neut_mut_rate", type=float,
            help="Neutral mutation rate", default=1e-7)
    parser.add_argument("-s","--selection_coefficient", dest="selection_coefficient", type=float,
            help="Selective advantage, s, of the homozygote (homozygote fitness = 1+s)", default=0.1)
    parser.add_argument("-h","--dominance_coefficient", dest="dominance_coefficient", type=float, 
            help="Dominance coefficient, h (heterozygote fitness = 1 + h*s)", default=0.5)
    parser.add_argument("-k","--nsamples", dest="nsamples", type=int,
            help="Number of *diploid* samples, total", default=100)
    parser.add_argument("--treefile_prefix","-t", type=str, dest="treefile_prefix",
            help="Prefix used when saving treefiles (will have freq+gens .hdf5 appended)", default="sweepfile")
    parser.add_argument("-B","--generations_before_sweep", dest="generations_before_sweep", type=int, default=0,
            help="Start introducing the selective variant this many generations into the simulation")
    parser.add_argument("-M","--max_generations", dest="max_generations", type=int, default=int(1e6),
            help="Abort the simulation after this many generations if the highest target frequency has not been reached")
    parser.add_argument("-G", "--gc", dest="simplify_interval", type=int,
            help="Interval between simplify steps.", default=500)
    parser.add_argument("-g","--logfile", dest="logfile", type=str,
            help="Name of log file (or '-' for stdout)", default="-")
    parser.add_argument("-v","--verbosity", action="count",
            help="Verbosity level", default=0)
    parser.add_argument("--track_mutations", dest="track_mutations", action="store_true",
            help=("If multiple treesequence files are produced, keep track of mutations forwards in time "
            "rather than adding mutations independently to each sample. This is slower, but produces " 
            "files which can be treated as successive samples of the same overall population"))
    parser.add_argument("--seed", "-d", dest="seed", type=int, help="random seed (if None, use simuPOP default)", default=None)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-of","--output_frequency", dest="output_frequency", 
            action='append', nargs="+", type=str, required=True, metavar=('FREQ', 'GENS'),
            help=("Output an msprime tree sequence file once the frequency of the selected variant has reached this value."
            "A second number may also be given, delaying output for that number of generations after the specified frequency "
            "is reached. Multiple uses of this parameter are allowed, with files output at each specified step. "
            "e.g. `-of 0.5 -of 0.8 -of 1.0 200` results in a file saved when the target variant attains frequency=0.5, "
            "again at freq 0.8, and another 200 generations after fixation (freq = 1.0)."
            "Once all specified files have been output, the simulation will stop."))        
    
    args = parser.parse_args()

    # do we need mutations added in forwards time (slower) or can they be added in retrospect
    mutations_after_simulation = args.track_mutations == False or len(args.output_frequency) == 1

    log_level = logging.WARNING
    if args.verbosity == 1:
        log_level = logging.INFO
    if args.verbosity >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(
        format='%(asctime)s %(message)s', level=log_level, stream=fileopt(args.logfile, "w"))
    
    logging.info("Options:")
    logging.info(str(args))
    logging.info("----------")
        
    saved_files = simulate_sweep(args.popsize, args.chrom_length, args.recomb_rate, args.neut_mut_rate, 
        args.selection_coefficient, args.dominance_coefficient, args.output_frequency, args.nsamples, 
        args.simplify_interval, generations_before_sweep=args.generations_before_sweep,
        max_generations=args.max_generations,  mutations_after_simulation=mutations_after_simulation,
        treefile_prefix=args.treefile_prefix, seed=args.seed, verbosity = args.verbosity)

    logging.info("All done, msprime files for specified output timepoints saved in: {}".format(saved_files))
    if args.verbosity > 1:
        #check this has all worked
        for freq,fn in saved_files.items():
            logging.info("Checking variants in file {}".format(fn))
            ts = msprime.load(fn)
            for v in ts.variants():
                logging.info(" variant pos: {} present in {}/{} samples ({:.2f}%)".format(
                    v.position, sum(v.genotypes), len(v.genotypes), sum(v.genotypes)/len(v.genotypes)*100))


if __name__ == '__main__':
    main()