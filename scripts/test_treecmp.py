#!/usr/bin/env python3

"""
Generate treefiles to test treecmp.r routine in R
"""
import os
import sys
import numpy as np
import os.path
sys.path.insert(1,os.path.join(sys.path[0],'..','fastARG')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(sys.path[0],'..','msprime')) # use the local copy of msprime in preference to the global one
sys.path.insert(1,os.path.join(sys.path[0],'..','scripts')) # use the local copy of msprime in preference to the global one
import msprime
from warnings import warn
from tempfile import NamedTemporaryFile
from msprime_fastARG import *
#from msprime_ARGweaver import *


def write_nexus_trees(ts, treefile, index_trees_by_variant_number=False):
    """
    if index_trees_by_variant_number == False (the default), then the names of the trees in the file correspond to the
    upper breakpoint position along the genome. If index_trees_by_variant_number == True then each tree name is instead
    the number of variants observed along the sequence so far (i.e. the uppermost variant index + 1). Using the variant
    number allows trees to be compared for batches of variants rather than chromosome position, allowing fair comparison 
    between inferred trees.
    """
    import sys
    
    print("#NEXUS\nBEGIN TREES;", file = treefile)
    print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(i,i) for i in range(ts.get_sample_size())])), file = treefile)
    trees = 0
    variant_index = 0
    epsilon = 1e-8
    for t, (_, newick) in zip(ts.trees(), ts.newick_trees()):
        trees += 1
        if index_trees_by_variant_number:
            #index by rightmost variant number
            n = t.get_num_mutations()
            if n:
                variant_index += n
                print("TREE " + str(variant_index) + " = " + newick, file=treefile)
        else:
            #index by rightmost genome position
            print("TREE " + str(t.get_interval()[1]) + " = " + newick, file=treefile)
    print("END;", file = treefile)

def test_fixed(coalescence_records, n_mutations, n_sim_replicates=1):
    import random
    import filecmp
    import os
    from subprocess import call
    from itertools import cycle
    from warnings import warn
    nexus_dir = "nex_files"
    orig_filenames = {"gpos":"original_trees_by_genome_pos.nex", "mpos":"original_trees_%d_muts_by_mut_idx.nex"}
    new_filenames = {"gpos":"fastarg_trees_%d_muts_by_genome_pos_rep%d.nex", "mpos":"fastarg_trees_%d_muts_by_mut_idx_rep%d.nex"}
    orig_filenames = {k:os.path.join(nexus_dir,v) for k,v in orig_filenames.items()}
    new_filenames = {k:os.path.join(nexus_dir,v) for k,v in new_filenames.items()}
    random.seed(1234)
    n_fastARG_replicates = 1 #should be the same each time
    mut_params = []
    rep_params = []
    with NamedTemporaryFile("w+") as tree_in:
        tree_in.write(coalescence_records)
        tree_in.flush()
        ts = msprime.load_txt(tree_in.name)
        with open(orig_filenames["gpos"], "w+") as nex:
            write_nexus_trees(ts, nex, index_trees_by_variant_number=False)
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
                    write_nexus_trees(ts, nex, index_trees_by_variant_number=True)
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
                            with open(new_filenames["gpos"] % (len(muts), sim_replicate), "w+") as nex:
                                write_nexus_trees(ts_new, nex, index_trees_by_variant_number=False)
                            with open(new_filenames["mpos"] % (len(muts), sim_replicate), "w+") as nex:
                                write_nexus_trees(ts_new, nex, index_trees_by_variant_number=True)
                            mut_params.append(len(muts))
                            rep_params.append(sim_replicate)
    print("output R commands to stdout:\n\n", file=sys.stderr)
    r_cmds = []
    r_cmds.append("orig_g <- read.nexus('%s', force.multi=TRUE)" % os.path.abspath(orig_filenames["gpos"]))
    r_cmds.append("muts <- c({})".format(",".join([str(m) for m in mut_params])))
    r_cmds.append("reps <- c({})".format(",".join([str(r) for r in rep_params])))
    r_cmds.append("names(muts) <- sprintf('%dmuts%d', muts, reps)")
    r_cmds.append("gpos_metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(orig_g, read.nexus(sprintf('{full_filename}',m,r), force.multi=TRUE)), muts, reps)))".format(full_filename=os.path.abspath(new_filenames["gpos"])))
    r_cmds.append("mpos_metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(read.nexus(sprintf('{full_orig_filename}',m), force.multi=TRUE), read.nexus(sprintf('{full_new_filename}',m,r), force.multi=TRUE)), muts, reps)))".format(full_orig_filename=os.path.abspath(orig_filenames["mpos"]), full_new_filename=os.path.abspath(new_filenames["mpos"])))
    r_cmds.append("gpos_rooted_metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(orig_g, read.nexus(sprintf('{full_filename}',m,r), force.multi=TRUE), rooted=TRUE), muts, reps)))".format(full_filename=os.path.abspath(new_filenames["gpos"])))
    r_cmds.append("mpos_rooted_metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(read.nexus(sprintf('{full_orig_filename}',m), force.multi=TRUE), read.nexus(sprintf('{full_new_filename}',m,r), force.multi=TRUE), rooted=TRUE), muts, reps)))".format(full_orig_filename=os.path.abspath(orig_filenames["mpos"]), full_new_filename=os.path.abspath(new_filenames["mpos"])))
    print("\n".join(r_cmds))

def test_sim(recipr_mut_rates=[50000000], n_sim_replicates=1, sample_size=8, length=1e4, Ne=1e4, recombination_rate=2e-8):
    import random
    import filecmp
    import os
    from subprocess import call
    from itertools import cycle
    from warnings import warn
    nexus_dir = "nex_files"
    orig_filenames = {"gpos":"original_trees_%d_muts_by_genome_pos_rep%d.nex", "mpos":"original_trees_%d_muts_by_mut_idx_rep%d.nex"}
    new_filenames = {"gpos":"fastarg_trees_%d_muts_by_genome_pos_rep%d.nex", "mpos":"fastarg_trees_%d_muts_by_mut_idx_rep%d.nex"}
    orig_filenames = {k:os.path.join(nexus_dir,v) for k,v in orig_filenames.items()}
    new_filenames = {k:os.path.join(nexus_dir,v) for k,v in new_filenames.items()}
    rng = msprime.RandomGenerator(1234)
    n_fastARG_replicates = 1 #should be the same each time
    mut_params = []
    rep_params = []
    for sim_replicate in range(n_sim_replicates):
        print("Simulating an ARG for replicate %d" % sim_replicate, file=sys.stderr)
        ts = msprime.simulate(sample_size=sample_size, Ne=Ne, mutation_rate=None, recombination_rate=recombination_rate, length=length)
        for r_mut in recipr_mut_rates:
            print("Throwing mutations onto the ARG at rate {} ".format(1.0/r_mut), file=sys.stderr)
            ts.generate_mutations(1.0/r_mut, rng)
            with open(orig_filenames["gpos"] % (r_mut, sim_replicate), "w+") as nex:
                write_nexus_trees(ts, nex, index_trees_by_variant_number=False)
            with open(orig_filenames["mpos"] % (r_mut, sim_replicate), "w+") as nex:
                write_nexus_trees(ts, nex, index_trees_by_variant_number=True)
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
                        fastARG_out_to_msprime_txts(fa_out, variant_positions_from_fastARGin(fa_in), tree, mutation_file, seq_len=ts.get_sequence_length(), status_to=None)
                        #read in the msprime input
                        ts_new = msprime_txts_to_fastARG_in_revised(tree, mutation_file, fastARG_root_seq(fa_out), fa_revised, status_to=None)
                        #quick check that the generated haplotypes are the same
                        if filecmp.cmp(fa_in.name, fa_revised.name, shallow=False) == False:
                            warn("Initial fastARG input file differs from processed fastARG file")
                        with open(new_filenames["gpos"] % (r_mut, sim_replicate), "w+") as nex:
                            write_nexus_trees(ts_new, nex, index_trees_by_variant_number=False)
                        with open(new_filenames["mpos"] % (r_mut, sim_replicate), "w+") as nex:
                            write_nexus_trees(ts_new, nex, index_trees_by_variant_number=True)
                        mut_params.append(r_mut)
                        rep_params.append(sim_replicate)
    print("output R commands to stdout:\n\n", file=sys.stderr)
    r_cmds = []
    r_cmds.append("muts <- c({})".format(",".join([str(m) for m in mut_params])))
    r_cmds.append("reps <- c({})".format(",".join([str(r) for r in rep_params])))
    r_cmds.append("names(muts) <- sprintf('%dmuts%d', muts, reps)")
    r_cmds.append("layout(matrix(1:4,2,2))")    
    r_cmds.append("for (rooted in c(FALSE, TRUE)) {")
    r_cmds.append(" for (filenames in list(genome.pos=c('%s','%s', 'using absolute genome positions'), mut.pos=c('%s','%s', 'using variant positions'))) {" % ( \
        os.path.abspath(orig_filenames["gpos"]), \
        os.path.abspath(new_filenames["gpos"]), \
        os.path.abspath(orig_filenames["mpos"]), \
        os.path.abspath(new_filenames["mpos"])))
    r_cmds.append("  metrics <- as.data.frame(t(mapply(function(m, r) genome.trees.dist(read.nexus(sprintf(filenames[1],m,r), force.multi=TRUE), read.nexus(sprintf(filenames[2],m,r), force.multi=TRUE), rooted=rooted), muts, reps)))")
    r_cmds.append(r"  mut.unstack <- unstack(metrics, 1/as.numeric(sub('muts\\d+','',rownames(metrics)))~as.numeric(sub('\\d+muts','',rownames(metrics))))")
    r_cmds.append(r"  met.unstack <- unstack(metrics, RF~as.numeric(sub('\\d+muts','',rownames(metrics))))")
    r_cmds.append("  matplot(mut.unstack, met.unstack, type='l', lty=1, log='x', main=ifelse(rooted,'Rooted metric', 'Unrooted metric'))")
    r_cmds.append("  mtext(filenames[3], line=0.5, cex=0.9)")
    r_cmds.append("}}")
    
    print("\n".join(r_cmds))

 
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generate treefiles to test treecmp.r')
    parser.add_argument('--sim', "-s", action='store_true', help='run a set of simulations, rather than use a fixed coalescence record structure')
    args = parser.parse_args()
    if args.sim:
        mut_rates = [1e-8, 2e-8, 5e-8, 1e-7, 5e-7, 1e-6, 3e-6, 5e-6]
        test_sim(recipr_mut_rates=[int(1.0/x) for x in mut_rates], n_sim_replicates=10)
    else:
        #test using the example in Fig 4 of the original 2015 msprime paper, which
        test_fixed("""\
        2       10      4       2,3         0.071    0
        0       2       5       1,3         0.090    0
        2       10      5       1,4         0.090    0
        0       7       6       0,5         0.170    0
        7       10      7       0,5         0.202    0
        0       2       8       2,6         0.253    0""", n_mutations=range(3,200,2), n_sim_replicates=1)