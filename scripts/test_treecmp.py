#!/usr/bin/env python3

"""
Generate treefiles to test treecmp.r routine in R
"""
import os
import sys
import numpy as np
import os.path
sys.path.insert(1,os.path.join(sys.path[0],'..','fastARG')) # use the local copy of fastARG
sys.path.insert(1,os.path.join(sys.path[0],'..','msprime')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(sys.path[0],'..','scripts'))
import msprime
from warnings import warn
from tempfile import NamedTemporaryFile, TemporaryDirectory
from msprime_fastARG import *
from msprime_ARGweaver import *


def write_nexus_trees(ts, treefile, index_trees_by_variants=False, zero_based_tip_numbers=True):
    """
    if index_trees_by_variants == False (the default), then the names of the trees in the file correspond to the
    upper breakpoint position along the genome. If index_trees_by_variants == True then each tree name is instead
    the number of variants observed along the sequence so far (i.e. the uppermost variant index + 1). 
    If index_tree_by_variants is a list, it is taken to be a list of variant positions to use as positions at which
    trees are output (with adjacent identical trees output as a single tree and indexed by the uppermost variant position,
    in the same manner as for index_trees_by_variants = True). 
    
    Using the variant position to index trees allows trees to be compared for batches of variants rather than 
    via chromosome position, allowing fair comparison between inferred trees.
    """
    import sys
    import numpy as np
    
    print("#NEXUS\nBEGIN TREES;", file = treefile)
    increment = 0 if zero_based_tip_numbers else 1
    print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i+increment,i+increment) for i in range(ts.get_sample_size())])), file = treefile)
    variant_index = 0
    epsilon = 1e-8
    if hasattr(index_trees_by_variants, "__len__"):
        #index_trees_by_variants contains an array of variant positions, which can be used to index the trees
        assert min(index_trees_by_variants) >= 0, "The variant positions passed in must all be greater than or equal to 0"
        assert max(index_trees_by_variants) < ts.get_sequence_length(), "The variant positions passed in must all be less than {}".format(ts.get_sequence_length())
        positions = np.sort(np.array(index_trees_by_variants))
    else:
        positions = None
        
    for t, (_, newick) in zip(ts.trees(), ts.newick_trees()): #TO DO: should do ts.newick_trees(zero_based_tip_numbers)
        if index_trees_by_variants == False:
            #index by rightmost genome position
            print("TREE " + str(t.get_interval()[1]) + " = " + newick, file=treefile)
        else:
            #index by rightmost variant number
            if positions is None: #use the built-in variant positions
                n = t.get_num_mutations()
            else:
                l, r = t.get_interval()
                n = np.count_nonzero((positions >= l) & (positions < r))
            if n:
                #only print a tree if is is covered by at least 1 variant
                variant_index += n
                print("TREE " + str(variant_index) + " = " + newick, file=treefile)
    print("END;", file = treefile)

def write_variant_positions(ts, pos_file):
    for v in ts.variants():
        print(v.position, file=pos_file)

def test_fixed(nexus_dir, coalescence_records, n_mutations, n_sim_replicates=1):
    import random
    import filecmp
    import os
    from subprocess import call
    from itertools import cycle
    from warnings import warn
    orig_filenames = {"gpos":"original_trees_by_genome_pos.nex", "mpos":"original_trees_%d_muts_by_mut_idx.nex"}
    fa_filenames = {"gpos":"fastarg_trees_%d_muts_by_genome_pos_rep%d.nex", "mpos":"fastarg_trees_%d_muts_by_mut_idx_rep%d.nex"}
    orig_filenames = {k:os.path.join(nexus_dir,v) for k,v in orig_filenames.items()}
    fa_filenames = {k:os.path.join(nexus_dir,v) for k,v in fa_filenames.items()}
    random.seed(1234)
    n_fastARG_replicates = 1 #should be the same each time
    mut_params = []
    rep_params = []
    with NamedTemporaryFile("w+") as tree_in:
        tree_in.write(coalescence_records)
        tree_in.flush()
        ts = msprime.load_txt(tree_in.name)
        with open(orig_filenames["gpos"], "w+") as nex:
            write_nexus_trees(ts, nex, index_trees_by_variants=False, zero_based_tip_numbers=False) #ape in R requires 1-based numbers in .nex files
        for n_muts in n_mutations:
            for sim_replicate in range(n_sim_replicates):
                muts = [None] * n_muts
                print("Placing %d mutations on the example tree" % n_muts, file=sys.stderr)
                #mutations are placed above internal nodes (but not the root).
                #so we can allocate mutations simply to child IDs in each coalescent record
                
                #for testing purposes, we allocate mutations to sequential coalescent records
                #so that we cover all the possible mutation places
                record_iterator = cycle(ts.records())
                while n_muts != 0:
                    rec = next(record_iterator)
                    for c in rec.children:
                        n_muts -= 1
                        muts[n_muts] = (random.uniform(rec.left, rec.right), c)
                        if n_muts==0:
                            break
                    
                ts.set_mutations(muts)
                with open(orig_filenames["mpos"]% len(muts), "w+") as nex:
                    write_nexus_trees(ts, nex, index_trees_by_variants=True, zero_based_tip_numbers=False)
                for inference_replicates in range(n_fastARG_replicates):
                        with NamedTemporaryFile("w+") as fa_in, \
                             NamedTemporaryFile("w+") as fa_out, \
                             NamedTemporaryFile("w+") as tree, \
                             NamedTemporaryFile("w+") as fa_revised, \
                             NamedTemporaryFile("w+") as mutation_file:
                            #create the fastarg input file
                            msprime_to_fastARG_in(ts, fa_in)
                            #run fastarg
                            run_fastARG("../fastARG/fastARG", fa_in, fa_out, seed=1234*(inference_replicates+2), status_to=None)
                            #turn the fastARG output into msprime input format
                            fastARG_out_to_msprime_txts(fa_out, variant_positions_from_fastARGin(fa_in), tree, mutation_file, ts.get_sequence_length(), status_to=None)
                            #read in the msprime input
                            ts_new = msprime_txts_to_fastARG_in_revised(tree, mutation_file, fastARG_root_seq(fa_out), fa_revised, status_to=None)
                            #quick check that the generated haplotypes are the same
                            if filecmp.cmp(fa_in.name, fa_revised.name, shallow=False) == False:
                                warn("Initial fastARG input file differs from processed fastARG file")
                            else:
                                #output them for inspection
                                call(["cp", fa_in.name, "%d_muts.haps" % len(muts)])
                            with open(fa_filenames["gpos"] % (len(muts), sim_replicate), "w+") as nex:
                                write_nexus_trees(ts_new, nex, index_trees_by_variants=False, zero_based_tip_numbers=False)
                            with open(fa_filenames["mpos"] % (len(muts), sim_replicate), "w+") as nex:
                                write_nexus_trees(ts_new, nex, index_trees_by_variants=True, zero_based_tip_numbers=False)
                            mut_params.append(len(muts))
                            rep_params.append(sim_replicate)
    print("output R commands to stdout:\n\n", file=sys.stderr)
    r_cmds = []
    r_cmds.append("orig_g <- read.nexus('%s', force.multi=TRUE)" % os.path.abspath(orig_filenames["gpos"]))
    r_cmds.append("muts <- c({})".format(",".join([str(m) for m in mut_params])))
    r_cmds.append("reps <- c({})".format(",".join([str(r) for r in rep_params])))
    r_cmds.append("names(muts) <- sprintf('%dmuts%d', muts, reps)")
    r_cmds.append("gpos_metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(orig_g, read.nexus(sprintf('{full_filename}',m,r), force.multi=TRUE)), muts, reps)))".format(full_filename=os.path.abspath(fa_filenames["gpos"])))
    r_cmds.append("mpos_metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(read.nexus(sprintf('{full_orig_filename}',m), force.multi=TRUE), read.nexus(sprintf('{full_new_filename}',m,r), force.multi=TRUE)), muts, reps)))".format(full_orig_filename=os.path.abspath(orig_filenames["mpos"]), full_new_filename=os.path.abspath(fa_filenames["mpos"])))
    r_cmds.append("gpos_rooted_metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(orig_g, read.nexus(sprintf('{full_filename}',m,r), force.multi=TRUE), rooted=TRUE), muts, reps)))".format(full_filename=os.path.abspath(fa_filenames["gpos"])))
    r_cmds.append("mpos_rooted_metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(read.nexus(sprintf('{full_orig_filename}',m), force.multi=TRUE), read.nexus(sprintf('{full_new_filename}',m,r), force.multi=TRUE), rooted=TRUE), muts, reps)))".format(full_orig_filename=os.path.abspath(orig_filenames["mpos"]), full_new_filename=os.path.abspath(fa_filenames["mpos"])))
    print("\n".join(r_cmds))

def msprime_filename(n, Ne, l, rc, mu, tree_seed, mut_seed, subset=None):
    if subset is None:
        prefix = "msprime"
    else:
        prefix="ms{0:02d}".format(subset)
    #some jiggery pokery needed to print non-exponential notation without trailing zeroes
    return(prefix + "-n{}_Ne{}_l{}_rc{}_mu{}-{}_{}".format(n, Ne, l, "{:.12f}".format(rc).rstrip('0'), "{:.12f}".format(mu).rstrip('0'), tree_seed, mut_seed))
    
def construct_fastarg_basename(seed, sim_basename, base_dir=None):
    """returns a fastARG filename but without file extension"""
    if base_dir is not None:
        import os
        return(os.path.join(base_dir, '+'.join(['fastarg', sim_basename, seed])))
    else:
        return('+'.join(['fastarg', sim_basename, str(seed)]))
        
def construct_argweaver_basename(seed, sim_basename, base_dir=None):
    """returns an ARGweaver filename but without .iter number or file extension"""
    if base_dir is not None:
        import os
        return(os.path.join(base_dir, '+'.join(['aweaver', sim_basename, seed])))
    else:
        return('+'.join(['aweaver', sim_basename, str(seed)]))

def construct_ls_basename(sim_basename, subset=None, base_dir=None):
    """returns an ARGweaver filename but without .iter number or file extension"""
    if subset is None:
        prefix = "ls{0:02d}".format(subset)
    else:
        prefix = "lsinfer"
            
    if base_dir is not None:
        import os
        return(os.path.join(base_dir, '+'.join([prefix, sim_basename, ""])))
    else:
        return('+'.join([prefix, sim_basename, ""]))
 
def seed_set(n, rng):
    """
    returns a set of n seeds suitable for seeding RNGs
    """
    seeds=set()
    while len(seeds) < n:
        seeds.add(rng.randint(1,999999))
    return(seeds)
    
def test_sim(nexus_dir, mut_rates=[2e-8], sample_size=8, length=1e4, Ne=1e4, recombination_rate=2e-8, rand_seed=1234, n_fastARG_replicates = 1):
    import random
    import filecmp
    import os
    from subprocess import call
    from itertools import cycle
    from warnings import warn
    import shutil
    import re
    import csv
    n_fastARG_replicates = 1 #should be the same each time
    fa_mut_params = []
    fa_rep_params = []
    aw_mut_params = []
    aw_rep_params = []
    aw_samp_params = []
    aw_likelihood = []
    random.seed(rand_seed)
    mutseeds = seed_set(len(mut_rates), random)
    ts = msprime.simulate(sample_size=sample_size, Ne=Ne, mutation_rate=None, recombination_rate=recombination_rate, length=length, random_seed=rand_seed)
    for mu, mut_seed in zip(mut_rates, mutseeds):
        print("Throwing mutations onto the ARG at rate {} ".format(mu), file=sys.stderr)
        rng = msprime.RandomGenerator(mut_seed)
        ts.generate_mutations(mu, rng)
        simname = msprime_filename(sample_size, Ne, length, recombination_rate, mu, rand_seed, mut_seed)
        with open(os.path.join(nexus_dir,simname+".nex"), "w+") as nex, \
             open(os.path.join(nexus_dir, simname + ".vpos"), "w+") as v_pos:
            write_nexus_trees(ts, nex, index_trees_by_variants=False, zero_based_tip_numbers=False)
            write_variant_positions(ts, v_pos)
        #create seeds for inference params
        for inference_seed in seed_set(n_fastARG_replicates, random):
            with TemporaryDirectory() as tmp:
                fastarg_base = construct_fastarg_basename(inference_seed, simname)
                with open(os.path.join(tmp, fastarg_base + ".fava"), "w+") as fa_in, \
                     open(os.path.join(tmp, fastarg_base + ".farg"), "w+") as fa_out, \
                     open(os.path.join(tmp, fastarg_base + ".mscr"), "w+") as tree, \
                     open(os.path.join(tmp, fastarg_base + ".new.fava"), "w+") as fa_revised, \
                     open(os.path.join(tmp, fastarg_base + ".msmu"), "w+")  as mutation_file, \
                     open(os.path.join(nexus_dir, fastarg_base + ".nex"), "w+") as nex_g:
                    #create the fastarg input file
                    msprime_to_fastARG_in(ts, fa_in)
                    #run fastarg
                    run_fastARG("../fastARG/fastARG", fa_in, fa_out, seed=inference_seed, status_to=None)
                    #turn the fastARG output into msprime input format
                    fastARG_out_to_msprime_txts(fa_out, variant_positions_from_fastARGin(fa_in), tree, mutation_file, seq_len=ts.get_sequence_length(), status_to=None)
                    #read in the msprime input
                    ts_fa = msprime_txts_to_fastARG_in_revised(tree, mutation_file, fastARG_root_seq(fa_out), fa_revised, status_to=None)
                    #quick check that the generated haplotypes are the same
                    if filecmp.cmp(fa_in.name, fa_revised.name, shallow=False) == False:
                        warn("Initial fastARG input file differs from processed fastARG file")
                    write_nexus_trees(ts_fa, nex_g, index_trees_by_variants=False, zero_based_tip_numbers=False)
                argweaver_base = construct_argweaver_basename(inference_seed, simname)
                with open(os.path.join(tmp, argweaver_base+".sites"), "w+") as aw_in:
                    msprime_to_ARGweaver_in(ts, aw_in)
                    aw_prefix=os.path.join(tmp, argweaver_base)
                    run_ARGweaver(Ne=Ne, mut_rate=mu, recomb_rate=recombination_rate, executable="../argweaver/bin/arg-sample", \
                                  rand_seed=inference_seed, quiet=True, out_prefix=aw_prefix, iterations=100, sample_step=10, \
                                  ARGweaver_in_filehandle=aw_in, ARGweaver_out_dir=tmp)
                    #convert all the generated smc files to nexus
                    for smc_file in os.listdir(tmp):
                        if smc_file.endswith(".smc.gz"):
                            prefix = smc_file.replace(".smc.gz", "")
                            with open(os.path.join(nexus_dir,prefix+".nex"), "w+") as nex_g:
                                ARGweaver_smc_to_nexus(os.path.join(tmp,smc_file), nex_g, zero_based_tip_numbers=False)
                    #copy the stats file so we can get at the likelihood etc.
                    shutil.copyfile(os.path.join(tmp, aw_prefix + ".stats"), os.path.join(nexus_dir, argweaver_base+".stats"))
                    #lik = {}
                    #with open(os.path.join(tmp, aw_prefix + ".stats"), "r+") as stats:
                    #    for line in csv.DictReader(stats, delimiter='\t'):
                    #        if line['stage'] == 'resample':
                    #            lik[line['iter']] = float(line['joint']) #TO DO - check with Gil if we should be using 'joint' or 'likelihood' or what
                    #convert to .nex files
                        
def tmp():
    print("output R commands to stdout:\n\n", file=sys.stderr)
    tree_stat = "RF" #which tree statistic to measure
    r_cmds = []
    r_cmds.append(r"fa.muts <- c({})".format(",".join([str(m) for m in fa_mut_params])))
    r_cmds.append(r"fa.reps <- c({})".format(",".join([str(r) for r in fa_rep_params])))
    r_cmds.append(r"names(fa.muts) <- sprintf('%dmuts%d', fa.muts, fa.reps)")
    r_cmds.append(r"aw.muts <- c({})".format(",".join([str(m) for m in aw_mut_params])))
    r_cmds.append(r"aw.reps <- c({})".format(",".join([str(r) for r in aw_rep_params])))
    r_cmds.append(r"aw.samp <- c({})".format(",".join(["'{}'".format(r) for r in aw_samp_params])))
    r_cmds.append(r"aw.loglik <- c({})".format(",".join(["{}".format(r) for r in aw_likelihood])))
    r_cmds.append(r"names(aw.muts) <- sprintf('%dmuts%d%s', aw.muts, aw.reps, aw.samp)")
    r_cmds.append(r"re.match <- '(\\d+)muts(\\d+)(.*)'")
    r_cmds.append(r"layout(matrix(1:8,4,2, byrow=TRUE))")    
    r_cmds.append(r"for (rooted in c(FALSE, TRUE)) {")
    r_cmds.append(r" for (filenames in list(genome.pos=c('%s','%s', '%s', 'using absolute genome positions'), mut.pos=c('%s','%s', '%s', 'using variant positions'))) {" % ( \
        os.path.abspath(orig_filenames["gpos"]), \
        os.path.abspath(fa_filenames["gpos"]), \
        os.path.abspath(aw_filenames["gpos"]), \
        os.path.abspath(orig_filenames["mpos"]), \
        os.path.abspath(fa_filenames["mpos"]), \
        os.path.abspath(aw_filenames["mpos"])))
    r_cmds.append(r"  fa.metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(read.nexus(sprintf(filenames[1],m,r), force.multi=TRUE), read.nexus(sprintf(filenames[2],m,r), force.multi=TRUE), rooted=rooted), fa.muts, fa.reps)))")
    r_cmds.append(r"  aw.metrics <- as.data.frame(t(mapply(function(m, r, s) genome.trees.dist(read.nexus(sprintf(filenames[1],m,r), force.multi=TRUE), read.nexus(sprintf(filenames[3],m,r,s), force.multi=TRUE), rooted=rooted), aw.muts, aw.reps, aw.samp)))")
    r_cmds.append(r"  mutation.rate <- unstack(fa.metrics, 1/as.numeric(sub(re.match,'\\1',rownames(fa.metrics)))~as.numeric(sub(re.match,'\\2',rownames(fa.metrics))))")
    r_cmds.append(r"  %s.metric <- unstack(fa.metrics, %s~as.numeric(sub(re.match,'\\2',rownames(fa.metrics))))" % (tree_stat, tree_stat))
    r_cmds.append(r"  #ARGweaver results are more tricky as there are multiple results per simulation")
    r_cmds.append(r"  sims <- sub(re.match,'\\1muts\\2',rownames(aw.metrics))")
    r_cmds.append(r"  aw.metrics.summary <- do.call(rbind,by(aw.metrics, factor(sims, levels=unique(sims)), function(x) data.frame(%s.mean=mean(x$%s), %s.sd=sd(x$%s))))" % (tree_stat, tree_stat, tree_stat, tree_stat))
    r_cmds.append(r"  aw.mutation.rate <- unstack(aw.metrics.summary, 1/as.numeric(sub(re.match,'\\1',rownames(aw.metrics.summary)))~as.numeric(sub(re.match,'\\2',rownames(aw.metrics.summary))))")
    r_cmds.append(r"  aw.%s.metric <- unstack(aw.metrics.summary, %s.mean~as.numeric(sub(re.match,'\\2',rownames(aw.metrics.summary))))" % (tree_stat, tree_stat))
    r_cmds.append(r"  aw.%s.sd <- unstack(aw.metrics.summary, %s.sd~as.numeric(sub(re.match,'\\2',rownames(aw.metrics.summary))))" % (tree_stat, tree_stat))

    r_cmds.append(r"  matplot(mutation.rate, %s.metric, type='l', lty=1, log='x', main=ifelse(rooted,'Rooted metric', 'Unrooted metric'))" % tree_stat)
    r_cmds.append(r"  abline(h=0,lty=3)")
    r_cmds.append(r"  mtext(filenames[4], line=0.5, cex=0.9)")
    r_cmds.append(r"  matplot(aw.mutation.rate, aw.%s.metric, type='l', lty=1, log='x', main=ifelse(rooted,'Rooted metric', 'Unrooted metric'))" % tree_stat)
    r_cmds.append(r"  #arrows(aw.mutation.rate, aw.%s.metric, type='l', lty=1, log='x', main=ifelse(rooted,'Rooted metric', 'Unrooted metric'))" % tree_stat)
    r_cmds.append(r"  abline(h=0,lty=3)")
    r_cmds.append(r"  mtext(filenames[4], line=0.5, cex=0.9)")
    r_cmds.append(r"}}")
    
    print("\n".join(r_cmds))

 
if __name__ == "__main__":
    import argparse
    import random
    parser = argparse.ArgumentParser(description='Generate treefiles to test treecmp.r')
    parser.add_argument('--sim', "-s", type=int, default=None, help='run this many replicate simulations, rather than using a fixed coalescence record structure')
    parser.add_argument('--save_nexus_dir', "-d", default='nexus_dir', help='where to save all the nexus files')
    args = parser.parse_args()
    if args.sim:
        mut_rates = [2e-8, 5e-8, 1e-7, 5e-7, 1e-6, 3e-6, 5e-6]
        seeds = seed_set(args.sim, random)
        for reps, seed in enumerate(seeds):
            print("Simulating an ARG for replicate %d" % reps, file=sys.stderr)
            test_sim(args.save_nexus_dir, mut_rates, rand_seed=seed)

    else:
        #test using the example in Fig 4 of the original 2015 msprime paper, which
        test_fixed(args.save_nexus_dir, """\
        2       10      4       2,3         0.071    0
        0       2       5       1,3         0.090    0
        2       10      5       1,4         0.090    0
        0       7       6       0,5         0.170    0
        7       10      7       0,5         0.202    0
        0       2       8       2,6         0.253    0""", n_mutations=range(3,200,2), n_sim_replicates=1)