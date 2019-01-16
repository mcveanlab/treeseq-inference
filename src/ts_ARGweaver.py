#!/usr/bin/env python3
"""
Various functions to convert a ts file to ARGweaver input format,
and from .arg files to tree seq input.

When run as a script, takes an msprime simulation in .trees format, saves to
ARGweaver input format (haplotype sequences), runs ARGweaver inference on it to
make .smc files, converts the .smc ARGweaver output files to ts input
files, reads these in, and checks that the msprime ts is the
same. We also test tree balance statistics (which are tip-label agnostic) and (if
sample <= 5) looks at all possible tip permutations to see if the trees are
essentially the same but with the labels lost.

E.g. to look over many small trees (seq length 1 million bases), try

python3 ./src/ts_ARGweaver.py tmp/AWtest -v -l 1000000

"""
import sys
import subprocess
import logging
import math
import re
import gzip
import csv
import os.path

import numpy as np

import msprime

class CyclicalARGError(Exception):
    """
    Exception raised when ARG Weaver generates a cyclical ARG. This is a bug in
    ARGWeaver, so there's nothing we can do about it other than catch the
    error and abort the conversion.

    See https://github.com/mdrasmus/argweaver/issues/19
    """

def tsfile_to_ARGweaver_in(trees, ARGweaver_filehandle):
    """
    take a .trees file, and convert it into an input file suitable for ARGweaver
    Returns the simulation parameters (Ne, mu, r) used to create the .trees file
    """
    logging.info("== Saving to ARGweaver input format ==")
    try:
        ts = msprime.load(trees.name) #trees is a fh
    except AttributeError:
        ts = msprime.load(trees)
    ts_to_ARGweaver_in(ts, ARGweaver_filehandle)
    #here we should extract the /provenance information from the .trees file and return 
    # {'Ne':XXX, 'mutation_rate':XXX, 'recombination_rate':XXX}
    #but this information is currently not encoded in the .trees file (listed as TODO)

    return {'Ne':None, 'mutation_rate':None, 'recombination_rate':None}

def ts_to_ARGweaver_in(ts, ARGweaver_filehandle):
    """
    Takes a TreeSequence, and outputs a file in .sites format, suitable for input
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
    That compares to tree sequences which use (0..N-1) - i.e. (0,N].
    """
    simple_ts = ts.simplify()
    print("\t".join(["NAMES"]+[str(x) for x in range(simple_ts.get_sample_size())]), file=ARGweaver_filehandle)
    print("\t".join(["REGION", "chr", "1", str(int(simple_ts.get_sequence_length()))]), file=ARGweaver_filehandle)
    genotypes = None
    position = 0
    for v in simple_ts.variants():
        if int(math.ceil(v.position)) != position:
            #this is a new position. Print the genotype at the old position, and then reset everything
            if position:
                print(position, "".join(np.where(genotypes==0,"A","T")), sep="\t", file=ARGweaver_filehandle)
            genotypes = v.genotypes
            position = int(math.ceil(v.position))
        else:
            genotypes =  np.logical_and(genotypes, v.genotypes)
    if position:
        print(position, "".join(np.where(genotypes==0,"A","T")), sep="\t", file=ARGweaver_filehandle)

    ARGweaver_filehandle.flush()
    ARGweaver_filehandle.seek(0)

def samples_to_ARGweaver_in(sample_data, ARGweaver_filehandle, infinite_sites=True):
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

    Note that ARGweaver uses position coordinates (1,N) - i.e. (0,N].
    That compares to tree sequences which use (0..N-1) - i.e. [0,N).
    """
    print("\t".join(["NAMES"]+[str(x) for x in range(sample_data.num_samples)]), file=ARGweaver_filehandle)
    print("\t".join(["REGION", "chr", "1", str(sample_data.sequence_length)]), file=ARGweaver_filehandle)
    
    position = sample_data.sites_position[:] #decompress all in one go to avoid sequential unpacking
    if infinite_sites:
        for id, genotype in sample_data.genotypes():
            print(position[id]+1, "".join(np.where(genotype==0,"A","T")),
                sep="\t", file=ARGweaver_filehandle)
    else:
        prev_position = 0
        ANDed_genotype = None
        for id, genotype in sample_data.genotypes():
            if int(math.ceil(position[id])) != prev_position:
                #this is a new position. Print the genotype at the old position, and then reset everything
                if prev_position:
                    print(prev_position, "".join(np.where(ANDed_genotype==0,"A","T")), sep="\t",
                        file=ARGweaver_filehandle)
                ANDed_genotype = genotype
            else:
                ANDed_genotype =  np.logical_and(ANDed_genotype, genotype)
            prev_position = int(math.ceil(position[id]))
        if ANDed_genotype is not None: # print out the last site
            print(prev_position, "".join(np.where(ANDed_genotype ==0,"A","T")), sep="\t", file=ARGweaver_filehandle)

    ARGweaver_filehandle.flush()
    ARGweaver_filehandle.seek(0)


def ARGweaver_smc_to_ts_txts(smc2bin_executable, prefix, nodes_fh, edges_fh):
    """
    convert the ARGweaver smc representation to tree sequence text format
    """
    logging.debug(
        "== Converting the ARGweaver smc output file '{}' to .arg format using '{}' ==".format(
            prefix + ".smc.gz", smc2bin_executable))
    subprocess.call([smc2bin_executable, prefix + ".smc.gz", prefix + ".arg"])
    with open(prefix + ".arg", "r+") as arg_fh:
        return ARGweaver_arg_to_ts_txts(arg_fh, nodes_fh, edges_fh)

def ARGweaver_arg_to_ts_txts(ARGweaver_arg_filehandle, nodes_fh, edges_fh):
    """
    convert the ARGweaver arg representation to tree sequence tables

    We need to split ARGweaver records that extend over the whole genome into sections
    that cover just that coalescence point.

    returns the mapping of ARGweaver node names to TS node names
    """
    logging.debug("== Converting .arg output to tree seq ==")
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

    for line_num, fields in enumerate(csv.DictReader(ARGweaver_arg_filehandle, delimiter='\t')):
        assert (fields['name'] not in ARG_node_times), \
                "duplicate node names identified: line {}".format(line_num)
        #HACK: make sure that parent nodes are strictly older than children.
        #This assumes that parents always have a higher node number
        ARG_node_times[fields['name']] = float(fields['age'])
        #we save info about nodes when looking at their children, so we
        # should save info into parent nodes
        if fields['parents'] == '':
            assert(root_node == None)
            root_node = fields['name']
            #don't need to record anything here, as we will grab details of the
            # root when looking at children
        else:
            if fields['event']=='recomb':
                #each recombination event has multiple parents
                for second_parent, parent in enumerate(fields['parents'].split(",")):
                    if parent not in ARG_nodes:
                        ARG_nodes[parent]={}
                    ARG_nodes[parent][fields['name']]=[
                        (float(fields['pos']) if second_parent else start),
                        (end if second_parent else float(fields['pos']))]
            else:
                #these should all have one parent
                if fields['parents'] not in ARG_nodes:
                    ARG_nodes[fields['parents']]={}
                ARG_nodes[fields['parents']][fields['name']]=[start,end]

                if fields['event']=='gene':
                    #we should trust the labels from 
                    node_names[fields['name']] = int(fields['name'])
                    tips.add(fields['name'])
    #now relabel the internal nodes
    for key in ARG_nodes:
        node_names[key]=len(node_names)

    #recursive hack to make times strictly decreasing, using depth-first topological
    # sorting algorithm
    def set_child_times(node_name, node_order, temporary_marks=set()):
        if node_name in ARG_nodes:
            if node_name in temporary_marks:
                raise CyclicalARGError(
                    "ARG has a cycle in it, around node {}. This should not be possible."
                    "Aborting this conversion!".format(node_name))
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

    print("id\tis_sample\ttime", file=nodes_fh)
    for node_name in sorted(node_names, key=node_names.get): #sort by id
        print("{id}\t{is_sample}\t{time}".format(
            id=node_names[node_name],
            is_sample=int(node_name in tips),
            time=ARG_node_times[node_name]),
            file=nodes_fh)

    print("left\tright\tparent\tchild", file=edges_fh)
    for node_name in sorted(ARG_node_times, key=ARG_node_times.get): #sort by time
        # look at the break points for all the child sequences, and break up
        # into that number of records
        try:
            children = ARG_nodes[node_name]
            assert all([ARG_node_times[child] < ARG_node_times[node_name] for child in children])
            breaks = set()
            for leftright in children.values():
                breaks.update(leftright)
            breaks = sorted(breaks)
            for i in range(1,len(breaks)):
                leftbreak = breaks[i-1]
                rightbreak = breaks[i]
                #The read_text function allows `child` to be a comma-separated list of children
                children_str = ",".join(map(str, sorted([
                    node_names[cnode] for cnode, cspan in children.items()
                        if cspan[0]<rightbreak and cspan[1]>leftbreak])))
                print("{left}\t{right}\t{parent}\t{children}".format(
                    left=leftbreak, right=rightbreak, parent=node_names[node_name],
                    children=children_str), file=edges_fh)
        except KeyError:
            #these should all be the tips
            assert node_name in tips, (
                "The node {} is not a parent of any other node, but is not a tip "
                "either".format(node_name))
    nodes_fh.flush()
    nodes_fh.seek(0)
    edges_fh.flush()
    edges_fh.seek(0)
    return node_names

def ARGweaver_smc_to_nexus(smc_filename, outfilehandle):
    """
    ARGweaver always exports smc trees with tips labelled from 0..N-1. 
    Whereas Nexus format expects 1..N, so we must always relabel 
    them. The true labels should be on the NAMES line
    """
    with (gzip.open(smc_filename, 'rt+') if smc_filename.endswith(".gz") else open(smc_filename, 'rt+')) as smc:
        print("#NEXUS\nBEGIN TREES;", file = outfilehandle)
        for line in smc:
            if line.startswith("NAMES"):
                names = [n for n in line.split() if n != "NAMES"]
                print("TRANSLATE\n{};".format(",\n".join([
                    "{} {}".format(int(i+1),n)
                    for i,n in enumerate(names)])), file=outfilehandle)
            elif line.startswith("TREE"):
                #hack the labels in the tree to increment each by 1, by looking for digits preceeded by '(' ')' or ','
                line = re.sub(r'(?<=[(,)])(\d+)',lambda m: str(int(m.group(1))+1), line)
                #rightmost sequence position (X) is correct ( < X )
                print(
                    re.sub(
                        r'^TREE\s+\d+\s+(\d+)\s+',
                        lambda m: "TREE "+ m.group(1) + " = ",
                        line.rstrip()),
                    file = outfilehandle)
        print("END;", file = outfilehandle)


def main(args):
    import os
    import itertools
    import subprocess
    from dendropy import TreeList
    from dendropy.calculate import treecompare
    import ts_extras
    def ts_txts_to_trees(ts_nodes, ts_edges, trees_outname=None):
        import shutil
        import msprime
        logging.info("== Converting new ts ARG to .trees ===")
        try:
            ts = msprime.load_text(nodes=ts_nodes, edges=ts_edges)
        except:
            logging.warning("Can't load the texts file properly. Saved copied to 'bad.nodes' & 'bad.edges' for inspection")
            shutil.copyfile(ts_nodes.name, "bad.nodes")
            shutil.copyfile(ts_edges.name, "bad.edges")
            raise
        logging.info("== loaded {}, {}===".format(ts_nodes.name, ts_edges.name))
        try:
            simple_ts = ts.simplify()
        except:
            ts.dump("bad.trees")
            logging.warning("Can't simplify. .trees file dumped to 'bad.trees'")
            raise
        if trees_outname:
            simple_ts.dump(trees_outname)
        return(simple_ts)

    msprime.TreeSequence.write_nexus_trees = ts_extras.write_nexus_trees
    iterations = 20
    full_prefix = os.path.join(args.outputdir, os.path.splitext(os.path.basename(args.trees_file))[0])
    with open(full_prefix+".sites", "w+") as aw_in:
        tsfile_to_ARGweaver_in(args.trees_file, aw_in)
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
        #now check that the smc file produced can be converted to nodes
        smc = full_prefix + "." + str(iterations) + ".smc.gz"
        assert os.path.isfile(smc),  "No output file names {}".format(smc)
        smc_nex = smc.replace(".smc.gz", ".nex")
        with open(smc_nex, "w+") as smc_nex_out:
            ARGweaver_smc_to_nexus(smc, smc_nex_out)
        arg_nex=smc.replace(".smc.gz", ".ts_nex")
        with open(smc.replace(".smc.gz", ".TSnodes"), "w+") as nodes, \
            open(smc.replace(".smc.gz", ".TSedges"), "w+") as edges, \
            open(arg_nex, "w+") as ts_nex:
            ARGweaver_smc_to_ts_txts(
                os.path.join(args.ARGweaver_executable_dir, args.ARGweaver_smc2arg_executable),
                smc.replace(".smc.gz", ""),
                nodes, edges)

            ts = ts_txts_to_trees(nodes, edges)
            ts.write_nexus_trees(ts_nex)
        
        smc_trees = TreeList.get(path=smc_nex, schema="nexus")
        arg_trees = TreeList.get(path=arg_nex, schema="nexus", 
            taxon_namespace=smc_trees[0].taxon_namespace)
        #zero_based_tip_numbers assumed False)
        #Check the smc trees against the ts-imported equivalents
        #NB, the ARGweaver output does not specify where mutations occur on the ARG, so we cannot
        #reconstruct the sequences implied by this ARG for testing purposes, and thus cannot compare
        #the original sequences with the reconstructed ones

        assert len(smc_trees)==len(arg_trees)
        assert [int(float(t.label)) for t in smc_trees] == [int(float(t.label)) for t in arg_trees]
        for i, (smc_tree, arg_tree) in enumerate(zip(smc_trees, arg_trees)):
            if treecompare.symmetric_difference(smc_tree, arg_tree) == 0:
                print("✓ Tree " + str(i+1) + 
                    " in AW SMC file is identical to that produced by SMC->ARG->STS")
            else:
                raise Exception("Tree {} differs\n".format(i+1) + \
                    smc_tree.label + " (smc) = " + smc_tree.as_string(schema="newick", 
                        suppress_edge_lengths=True, 
                        suppress_internal_node_labels = True, 
                        suppress_rooting = True) + \
                    arg_tree.label + " (arg) = " + arg_tree.as_string(schema="newick", 
                        suppress_edge_lengths=True, 
                        suppress_internal_node_labels = True,
                        suppress_rooting = True))
                
            
if __name__ == "__main__":
    import argparse
    import filecmp
    import os
    parser = argparse.ArgumentParser(description='Check ARGweaver imports by running an msprime simulation to create an ARGweaver import file, inferring some args from it in smc format, converting the .smc format to .arg format, reading the .arg into msprime, and comparing the nexus output trees with the trees in the .smc file. This testing process requires the dendropy library')
    parser.add_argument('--trees_file', type=argparse.FileType('r', encoding='UTF-8'), default=None,
        help='an msprime .trees file. If none, simulate one with defaults')
    parser.add_argument('--ARGweaver_executable_dir', '-d',
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','argweaver/bin/'),
        help='the path to the directory containing the ARGweaver executables')
    parser.add_argument('--ARGweaver_sample_executable', '-x', default="arg-sample", help='the name of the ARGweaver executable')
    parser.add_argument('--ARGweaver_smc2arg_executable', '-s', default="smc2arg", help='the name of the ARGweaver executable')
    parser.add_argument('--sample_size', '-n', type=int, default=5, help='the sample size if a .trees file is not given')
    parser.add_argument('--effective_population_size', '-Ne', type=float, default=5000, help='the effective population size if a .trees file is not given')
    parser.add_argument('--sequence_length', '-l', type=float, default=55000, help='the sequence length if a .trees file is not given')
    parser.add_argument('--recombination_rate', '-rho', type=float, default=2.5e-8, help='the recombination rate if a. trees file is not given')
    parser.add_argument('--mutation_rate', '-mu', type=float, default=5e-8, help='the mutation rate if a .trees file is not given')
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

    if args.trees_file is None:
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
            if args.trees_file is None:
                args.trees_file = os.path.join(aw_out_dir, "sim.trees")
                ts.dump(args.trees_file, zlib_compression=True)
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
        if args.trees_file is None:
            args.trees_file = os.path.join(args.outputdir, "sim.trees")
            ts.dump(args.trees_file, zlib_compression=True)
        else:
            args.trees_file = args.trees_file.name
        main(args)
