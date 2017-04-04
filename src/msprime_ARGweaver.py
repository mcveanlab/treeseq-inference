#!/usr/bin/env python3
"""Various functions to convert msprime simulation file to ARGweaver input format, and from .arg files to msprime input.

When run as a script, takes an msprime simulation in hdf5 format, saves to ARGweaver input format (haplotype sequences), 
runs ARGweaver inference on it to make .smc files, converts the .smc ARGweaver output files to msprime input files, 
reads these into msprime, and checks that the trees from msprime are the same/. Dues to a bug in smc2arg 
(https://github.com/mdrasmus/argweaver/issues/20) these are *NOT* the same, but this tests tree balance statistics 
(which are tip-label agnostic) and (if sample <= 5) looks at all possible tip permutations to see if the trees are
essentially the same but with the labels lost.

E.g. to look ove many small trees (seq length 1 million bases), try

python3 ./src/msprime_ARGweaver.py tmp/AWtest -v -l 1000000 


"""
import sys
from math import ceil
import numpy as np
import os.path
sys.path.insert(1,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','msprime')) # use the local copy of msprime in preference to the global one
import msprime
import logging

def msprime_hdf5_to_ARGweaver_in(msprime_hdf5, ARGweaver_filehandle):
    """
    take an hdf5 file, and convert it into an input file suitable for ARGweaver
    Returns the simulation parameters (Ne, mu, r) used to create the hdf5 file
    """
    logging.info("== Saving to ARGweaver input format ==")
    try:
        ts = msprime.load(msprime_hdf5.name) #msprime_hdf5 is a fh
    except AttributeError:
        ts = msprime.load(msprime_hdf5)
    msprime_to_ARGweaver_in(ts, ARGweaver_filehandle)
    #here we should extract the /provenance information from the hdf5 file and return {'Ne':XXX, 'mutation_rate':XXX, 'recombination_rate':XXX}
    #but this information is currently not encoded in the hdf5 file (listed as TODO)

    return {'Ne':None, 'mutation_rate':None, 'recombination_rate':None}
    
def msprime_to_ARGweaver_in(ts, ARGweaver_filehandle):
    """
    Takes an msprime TreeSequence, and outputs a file in .sites format, suitable for input 
    into ARGweaver (see http://mdrasmus.github.io/argweaver/doc/#sec-file-sites)
    The documentation (http://mdrasmus.github.io/argweaver/doc/#sec-prog-arg-sample) 
    states that the only mutation model is Jukes-Cantor (i.e. equal mutation between 
    all bases). Assuming adjacent sites are treated independently, we convert variant 
    format (0,1) to sequence format (A, T, G, C) by simply converting 0->A and 1->T
    
    Msprime simulations assume infinite sites by allowing mutations to occur at 
    floating-point positions along a sequence. ARGweaver has discrete sites instead. 
    This routine implements a basic discretising function, which simply rounds upwards
    to the nearest int, ANDing the results if 2 or more variants end up at the same 
    integer position.

    Note that ARGweaver uses position coordinates (1,N) - i.e. [0,N). 
    That compares to msprime which uses (0..N-1) - i.e. (0,N].
    """
    simple_ts = ts.simplify()
    print("\t".join(["NAMES"]+[str(x) for x in range(simple_ts.get_sample_size())]), file=ARGweaver_filehandle)
    print("\t".join(["REGION", "chr", "1", str(int(simple_ts.get_sequence_length()))]), file=ARGweaver_filehandle)
    genotypes = None
    position = 0
    for v in simple_ts.variants():
        if int(ceil(v.position)) != position:
            #this is a new position. Print the genotype at the old position, and then reset everything
            if position:
                print(position, "".join(np.where(genotypes==0,"A","T")), sep="\t", file=ARGweaver_filehandle)
            genotypes = v.genotypes
            position = int(ceil(v.position))
        else:
            genotypes =  np.logical_and(genotypes, v.genotypes)
    if position:
        print(position, "".join(np.where(genotypes==0,"A","T")), sep="\t", file=ARGweaver_filehandle)

    ARGweaver_filehandle.flush()
    ARGweaver_filehandle.seek(0)

def variant_matrix_to_ARGweaver_in(var_matrix, var_positions, seq_length, ARGweaver_filehandle, infinite_sites=True):
    """
    Takes an variant matrix, and outputs a file in .sites format, suitable for input 
    into ARGweaver (see http://mdrasmus.github.io/argweaver/doc/#sec-file-sites)
    The documentation (http://mdrasmus.github.io/argweaver/doc/#sec-prog-arg-sample) 
    states that the only mutation model is Jukes-Cantor (i.e. equal mutation between 
    all bases). Assuming adjacent sites are treated independently, we convert variant 
    format (0,1) to sequence format (A, T, G, C) by simply converting 0->A and 1->T
    
    if infinite_sites==False, then use a basic discretising function, which simply rounds 
    upwards to the nearest int, ANDing the results if 2 or more variants end up at the same 
    integer position.
    
    Note that ARGweaver uses position coordinates (1,N) - i.e. [0,N). 
    That compares to msprime which uses (0..N-1) - i.e. (0,N].
    """
    n_variants, n_samples = var_matrix.shape
    assert len(var_matrix)==n_variants
    print("\t".join(["NAMES"]+[str(x) for x in range(n_samples)]), file=ARGweaver_filehandle)
    print("\t".join(["REGION", "chr", "1", str(seq_length)]), file=ARGweaver_filehandle)
    genotypes = None
    if infinite_sites:
        for pos, v in zip(var_positions, var_matrix):
            print(pos+1, "".join(np.where(v==0,"A","T")), sep="\t", file=ARGweaver_filehandle)
    else:
        prev_position = 0
        ANDed_genotype = None
        for pos, genotype in zip(var_positions, var_matrix):
            if int(ceil(pos)) != prev_position:
                #this is a new position. Print the genotype at the old position, and then reset everything
                if prev_position:
                    print(prev_position, "".join(np.where(ANDed_genotype==0,"A","T")), sep="\t", 
                        file=ARGweaver_filehandle)
                ANDed_genotype = genotype
            else:
                ANDed_genotype =  np.logical_and(ANDed_genotype, genotype)
            prev_position = int(ceil(pos))
        if ANDed_genotype is not None: # print out the last site
            print(prev_position, "".join(np.where(ANDed_genotype ==0,"A","T")), sep="\t", file=ARGweaver_filehandle)

    ARGweaver_filehandle.flush()
    ARGweaver_filehandle.seek(0)
        
    
def ARGweaver_smc_to_msprime_txts(smc2bin_executable, prefix, nodes_fh, edgesets_fh, override_assertions=False):
    """
    convert the ARGweaver smc representation to coalescence records format
    """
    assert override_assertions, "smc2arg is currently broken and should not be used." + \
        "See https://github.com/mdrasmus/argweaver/issues/20"
    from subprocess import call
    logging.debug("== Converting the ARGweaver smc output file '{}' to .arg format using '{}' ==".format(\
        prefix + ".smc.gz", smc2bin_executable))
    call([smc2bin_executable, prefix + ".smc.gz", prefix + ".arg"])
    with open(prefix + ".arg", "r+") as arg_fh:
        return(ARGweaver_arg_to_msprime_txts(arg_fh, nodes_fh, edgesets_fh))

def ARGweaver_arg_to_msprime_txts(ARGweaver_arg_filehandle, nodes_fh, edgesets_fh):
    """
    convert the ARGweaver arg representation to coalescence records format
    
    We need to split ARGweaver records that extend over the whole genome into sections that cover just that 
    coalescence point.
    
    returns the mapping of ARGweaver node names to msprime node names
    """
    import re
    import csv
    logging.debug("== Converting .arg output to msprime ==")
    ARG_nodes={} #cr[X] = child1:[left,right], child2:[left,right],... : serves as intermediate ARG storage 
    ARG_node_times={} #node_name => time
    node_names={} #map of ARGweaver names -> numbers
    tips = set()
    root_node = None

    #first row gives start and end
    ARGweaver_arg_filehandle.seek(0)
    firstline = next(ARGweaver_arg_filehandle)
    m = re.match(r'^start=(\d+)\s+end=(\d+)\s*$', firstline)
    if m:
        start=float(m.group(1))
        end=float(m.group(2))
    else:
        raise ValueError("Could not find start and end positions in .arg file")

    try:
        for line_num, fields in enumerate(csv.DictReader(ARGweaver_arg_filehandle, delimiter='\t')):
            assert (fields['name'] not in ARG_node_times), "duplicate node names identified: line {}".format(line_num)
            #HACK: make sure that parent nodes are strictly older than children. 
            #This assumes that parents always have a higher node number
            ARG_node_times[fields['name']] = float(fields['age'])
            #we save info about nodes when looking at their children, so we should save info into parent nodes
            if fields['parents'] == '':
                assert(root_node == None)
                root_node = fields['name']
                #don't need to record anything here, as we will grab details of the root when looking at children
            else:
                if fields['event']=='recomb':
                    #each recombination event has multiple parents
                    
                    for second_parent, parent in enumerate(fields['parents'].split(",")):
                        if parent not in ARG_nodes:
                            ARG_nodes[parent]={}
                        ARG_nodes[parent][fields['name']]=[(float(fields['pos']) if second_parent else start), (end if second_parent else float(fields['pos']))]
                else:
                    #these should all have one parent
                    if fields['parents'] not in ARG_nodes:
                        ARG_nodes[fields['parents']]={}
                    ARG_nodes[fields['parents']][fields['name']]=[start,end]
                    
                    if fields['event']=='gene':
                        node_names[fields['name']] = len(node_names)
                        tips.add(fields['name'])
        #now relabel the internal nodes
        for key in ARG_nodes:
            node_names[key]=len(node_names)
        
        
        #recursive hack to make times strictly decreasing, using depth-first topological sorting algorithm
        def set_child_times(node_name, node_order, temporary_marks=set()):
            if node_name in ARG_nodes:
                assert node_name not in temporary_marks, "ARG has a cycle in it, around node {}. This should not be possible. Aborting this conversion!".format(node_name)
                if node_name not in node_order:
                    temporary_marks.add(node_name)
                    for child_name in ARG_nodes[node_name]:
                        set_child_times(child_name, node_order, temporary_marks)
                    node_order.append(node_name)
                    temporary_marks.remove(node_name)
                    
        node_order = [] #contains the internal nodes, such that parent is always after child
        set_child_times(root_node, node_order)
    
        max_epsilon = len(node_order)
        for epsilon, nm in enumerate(node_order):
            ARG_node_times[nm] += 0.001 * (epsilon+1) / max_epsilon
        
        cols = {
            'nodes':     ["id","is_sample","time"],
            'edgesets':  ["left","right","parent","children"]
        }
        lines = {k:"\t".join(["{%s}" % s for s in v]) for k,v in cols.items()}
        for k,v in cols.items():
            #print the header lines
            print("\t".join(v), file=locals()[k+'_fh'])


        for node_name in sorted(node_names, key=node_names.get): #sort by id
            print(lines['nodes'].format(id=node_names[node_name], 
                is_sample=int(node_name in tips), 
                time=ARG_node_times[node_name]),
                file=nodes_fh)


        for node_name in sorted(ARG_node_times, key=ARG_node_times.get): #sort by time
            #look at the break points for all the child sequences, and break up into that number of records
            try:
                children = ARG_nodes[node_name]
                assert all([ARG_node_times[child]<ARG_node_times[node_name] for child in children]), "ARGweaver node {} has an adjusted time of {} but its children ({}) are not strictly younger (times @ {}).".format(node_name, str(ARG_node_times[node_name]), ", ".join(children), ", ".join([str(ARG_node_times[c]) for c in children]))
                breaks = set()
                for leftright in children.values():
                    breaks.update(leftright)
                breaks = sorted(breaks)
                for i in range(1,len(breaks)):
                    leftbreak = breaks[i-1]    
                    rightbreak = breaks[i]
                    #NB - here we could try to input the values straight into an msprime python structure,
                    #but until this feature is implemented in msprime, we simply output to a correctly formatted text file
                    print(lines['edgesets'].format(left=leftbreak, right=rightbreak, 
                        parent=node_names[node_name],
                        children=",".join([str(x) for x in sorted([node_names[cnode] for cnode, cspan in children.items() if cspan[0]<rightbreak and cspan[1]>leftbreak])])),
                        file=edgesets_fh)
            except KeyError:
                #these should all be the tips
                assert node_name in tips, "The node {} is not a parent of any other node, but is not a tip either".format(node_name)
    except AssertionError as e:
        import sys
        import shutil
        import os
        working_dir = os.path.dirname(os.path.realpath(ARGweaver_arg_filehandle.name))
        shutil.make_archive('bad', 'zip', working_dir)
        raise type(e)(str(e) + "\n" +
                  "A copy of the directory '{}' has been compressed to 'bad.zip' for debugging purposes\n".format(ARGweaver_arg_filehandle.name) +
                  "You should inspect the file '{}' within that archive".format(os.path.basename(ARGweaver_arg_filehandle.name))).with_traceback(sys.exc_info()[2])
    
    for k in cols.keys():
        locals()[k+'_fh'].flush()
        locals()[k+'_fh'].seek(0)
    return(node_names)
    
def ARGweaver_smc_to_nexus(smc_filename, outfilehandle, zero_based_tip_numbers=True):
    """
    setting zero_based_tip_numbers=False changes the tip labelling to be more like the standard NEXUS format, with tips labelled from 1..N not 0..(N-1)
    """
    import gzip
    import re
    increment = 0 if zero_based_tip_numbers else 1
    with (gzip.open(smc_filename, 'rt+') if smc_filename.endswith(".gz") else open(smc_filename, 'rt+')) as smc:
        print("#NEXUS\nBEGIN TREES;", file = outfilehandle)
        for line in smc:
            if line.startswith("NAMES"):
                print("TRANSLATE\n{};".format(",\n".join(["{} {}".format(int(i)+increment,int(i)+increment) for i in line.split() if i != "NAMES"])), file = outfilehandle)
            elif line.startswith("TREE"):
                if not zero_based_tip_numbers:
                    #hack the labels in the tree, by looking for digits precceded by '(' ')' or ','
                    line = re.sub(r'(?<=[(,)])(\d+)',lambda m: str(int(m.group(1))+1), line)
                #rightmost sequence position (X) is correct ( < X )
                print(re.sub(r'^TREE\s+\d+\s+(\d+)\s+',lambda m: "TREE "+ m.group(1) + " = ", line.rstrip()), file = outfilehandle)
        print("END;", file = outfilehandle)


def main(args):
    import os
    import itertools
    import subprocess
    from dendropy import TreeList, calculate
    import msprime_extras
    def msprime_txts_to_hdf5(msprime_nodes, msprime_edgesets, hdf5_outname=None):
        import shutil
        import msprime
        logging.info("== Converting new msprime ARG as hdf5 ===")
        try:
            ts = msprime.load_text(nodes=msprime_nodes, edgesets=msprime_edgesets).simplify()
        except:
            logging.warning("Can't load the texts file properly. Saved copied to 'bad.nodes' & 'bad.edgesets' for inspection")
            shutil.copyfile(msprime_nodes.name, "bad.nodes")
            shutil.copyfile(msprime_edgesets.name, "bad.edgesets")
            raise
        logging.info("== loaded {}, {}===".format(msprime_nodes.name, msprime_edgesets.name))
        try:
            simple_ts = ts.simplify()
        except:
            ts.dump("bad.hdf5")
            logging.warning("Can't simplify. HDF5 file dumped to 'bad.hdf5'")
            raise
        if hdf5_outname:
            simple_ts.dump(hdf5_outname)
        return(simple_ts)

    msprime.TreeSequence.write_nexus_trees = msprime_extras.write_nexus_trees
    iterations = 20
    full_prefix = os.path.join(args.outputdir, os.path.splitext(os.path.basename(args.hdf5file))[0])
    with open(full_prefix+".sites", "w+") as aw_in:
        msprime_hdf5_to_ARGweaver_in(args.hdf5file, aw_in)
        cmd = [os.path.join(args.ARGweaver_executable_dir, args.ARGweaver_sample_executable), 
            '--sites', aw_in.name,
            '--popsize', str(args.effective_population_size),
            '--recombrate', str(args.recombination_rate),
            '--mutrate', str(args.mutation_rate),
            '--overwrite',
            '--randseed', str(int(args.random_seed)),
            '--iters', str(iterations),
            '--sample-step', str(iterations),
            '--output', full_prefix]
        assert os.stat(aw_in.name).st_size > 0,  "Initial .sites file is empty"
        logging.debug("running '{}'".format(" ".join(cmd)))
        subprocess.call(cmd)
        smc = full_prefix + "." + str(iterations) + ".smc.gz"
        assert os.path.isfile(smc),  "No output file names {}".format(smc)
        smc_nex = smc.replace(".smc.gz", ".nex")
        with open(smc_nex, "w+") as smc_nex_out:
            ARGweaver_smc_to_nexus(smc, smc_nex_out, zero_based_tip_numbers=False)
        arg_nex=smc.replace(".smc.gz", ".msp_nex")
        with open(smc.replace(".smc.gz", ".TSnodes"), "w+") as nodes, \
            open(smc.replace(".smc.gz", ".TSedgesets"), "w+") as edgesets, \
            open(arg_nex, "w+") as msp_nex:
            ARGweaver_smc_to_msprime_txts(
                os.path.join(args.ARGweaver_executable_dir, args.ARGweaver_smc2arg_executable), 
                smc.replace(".smc.gz", ""),
                nodes, edgesets,
                override_assertions=True)
                
            ts = msprime_txts_to_hdf5(nodes, edgesets)
            ts.write_nexus_trees(msp_nex, zero_based_tip_numbers=False)
        smc_trees = TreeList.get(path=smc_nex, schema="nexus")
        arg_trees = TreeList.get(path=arg_nex, schema='nexus') 
        #zero_based_tip_numbers assumed False)
        #Check the smc trees against the msprime-imported equivalents
        #NB, the ARGweaver output does not specify where mutations occur on the ARG, so we cannot
        #reconstruct the sequences implied by this ARG for testing purposes, and thus cannot compare
        #the original sequences with the reconstructed ones
            
        assert len(smc_trees)==len(arg_trees), "number of trees in original and msprime-processed files is not the same"
        assert [int(float(t.label)) for t in smc_trees] == [int(float(t.label)) for t in arg_trees], "names are different"
        tot=0
        for (smc_tree, arg_tree) in zip(smc_trees, arg_trees):
            tot+=abs(calculate.treemeasure.colless_tree_imbalance(smc_tree, None)-
                calculate.treemeasure.colless_tree_imbalance(arg_tree, None))
        print("sum of tree balance differences is {} (should be 0) over {} trees of {} tips".format(
            tot, len(smc_trees), len(smc_trees.taxon_namespace)))
        
        if ts.get_sample_size() <= 5:
            stats={}
            print("Testing all permutations of tips")
            for i, perm in enumerate(itertools.permutations(range(1, ts.get_sample_size()+1))):
                arg_trees = TreeList.get(path=arg_nex, schema='nexus')
                for taxon in arg_trees.taxon_namespace:
                    taxon.label = perm[int(taxon.label)-1]
                test_trees = arg_trees.migrate_taxon_namespace(smc_trees.taxon_namespace)
                tot=0
                for (smc_tree, arg_tree) in zip(smc_trees, arg_trees):
                    tot+=calculate.treecompare.symmetric_difference(smc_tree, arg_tree)
                stats[i] = tot, perm
            
            for i in sorted(stats, key=stats.get, reverse=True):
                print("Permutation {}, sum stat = {} over {} trees".format(
                    (i, stats[i][1]) if stats[i][0]==0 else i, stats[i][0], len(smc_trees)))

if __name__ == "__main__":
    import argparse
    import filecmp
    import os
    parser = argparse.ArgumentParser(description='Check ARGweaver imports by running an msprime simulation to create an ARGweaver import file, inferring some args from it in smc format, converting the .smc format to .arg format, reading the .arg into msprime, and comparing the nexus output trees with the trees in the .smc file. This testing process requires the dendropy library')
    parser.add_argument('--hdf5file', type=argparse.FileType('r', encoding='UTF-8'), default=None, 
        help='an msprime hdf5 file. If none, simulate one with defaults')
    parser.add_argument('--ARGweaver_executable_dir', '-d', 
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','argweaver/bin/'), 
        help='the path to the directory containing the ARGweaver executables')
    parser.add_argument('--ARGweaver_sample_executable', '-x', default="arg-sample", help='the name of the ARGweaver executable')
    parser.add_argument('--ARGweaver_smc2arg_executable', '-s', default="smc2arg", help='the name of the ARGweaver executable')
    parser.add_argument('--sample_size', '-n', type=int, default=5, help='the sample size if an hdf5 file is not given')
    parser.add_argument('--effective_population_size', '-Ne', type=float, default=5000, help='the effective population size if an hdf5 file is not given')
    parser.add_argument('--sequence_length', '-l', type=float, default=55000, help='the sequence length if an hdf5 file is not given')
    parser.add_argument('--recombination_rate', '-rho', type=float, default=2.5e-8, help='the recombination rate if an hdf5 file is not given')
    parser.add_argument('--mutation_rate', '-mu', type=float, default=5e-8, help='the mutation rate if an hdf5 file is not given')
    parser.add_argument('--random_seed', '-seed', type=int, default=1234, help='a random seed for msprime & AW simulation')
    parser.add_argument('outputdir', nargs="?", default=None, help='the directory in which to store the intermediate files. If None, files are saved under temporary names')
    parser.add_argument('--verbosity', '-v', action='count', default=0)
    args = parser.parse_args()
    log_level = logging.WARNING
    if args.verbosity == 1:
        log_level = logging.INFO
    if args.verbosity >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(
        format='%(asctime)s %(message)s', level=log_level, stream=sys.stdout)

    if args.hdf5file is None:
        logging.info("Running a new simulation with n {}, Ne {}, l {}, rho {}, mu {}".format(
        args.sample_size, args.effective_population_size, args.sequence_length,
        args.recombination_rate, args.mutation_rate))
        ts = msprime.simulate(
            sample_size = args.sample_size, 
            Ne=args.effective_population_size, 
            length=args.sequence_length,
            recombination_rate=args.recombination_rate, 
            mutation_rate=args.mutation_rate, 
            random_seed=args.random_seed)
    else:
        logging.warning("Loading a user-specified simulation file: WARNING, argweaver may end up being run with different parameters from the simulation")
    if args.outputdir == None:
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as aw_out_dir:
            logging.info("Saving everything to temporary files (temporarily stored in {})".format(aw_out_dir))
            args.outputdir = aw_out_dir
            if args.hdf5file is None:
                args.hdf5file = os.path.join(aw_out_dir, "sim.hdf5")
                ts.dump(args.hdf5file, zlib_compression=True)
            main(args)
    else:
        if not os.path.isdir(args.outputdir):
            logging.info("Output dir {} does not exist: creating it".format(args.outputdir))
            os.mkdir(args.outputdir)
        if len(os.listdir(args.outputdir)) > 0:
            logging.info("Output dir {} already contains files: deleting them".format(args.outputdir))
            import shutil
            shutil.rmtree(args.outputdir) 
            os.mkdir(args.outputdir)
        if args.hdf5file is None:
            args.hdf5file = os.path.join(args.outputdir, "sim.hdf5")
            ts.dump(args.hdf5file, zlib_compression=True)
        else:
            args.hdf5file = args.hdf5file.name
        main(args)
