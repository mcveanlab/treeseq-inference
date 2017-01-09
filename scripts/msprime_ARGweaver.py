#!/usr/bin/env python3
"""Various functions to convert msprime simulation file to ARGweaver input format, and from .arg files to msprime input.

When run as a script, takes an msprime simulation in hdf5 format, saves to ARGweaver input format (haplotype sequences), runs ARGweaver inference on it, converts the ARGweaver output file (arg and mutations) to msprime input files, reads these into msprime outputs the haplotype sequences again, and checks that the haplotype sequences are the same.

E.g. for 8.hdf files produced by generate_data.py

./msprime_ARGweaver.py ../test_files/8.hdf5 -d ../argweaver/bin/

"""
import sys
import os.path
sys.path.insert(1,os.path.join(sys.path[0],'..','msprime')) # use the local copy of msprime in preference to the global one
import msprime
from warnings import warn

def msprime_hdf5_to_ARGweaver_in(msprime_hdf5, ARGweaver_filehandle, status_to=sys.stdout):
    """
    take an hdf5 file, and convert it into an input file suitable for ARGweaver
    Returns the simulation parameters (Ne, mu, r) used to create the hdf5 file
    """
    if status_to:
        print("== Saving to ARGweaver input format ==", file=status_to)
    ts = msprime.load(msprime_hdf5.name)
    msprime_to_ARGweaver_in(ts, ARGweaver_filehandle)
    #here we should extract the /provenance information from the hdf5 file and return {'Ne':XXX, 'mutation_rate':XXX, 'recombination_rate':XXX}
    #but this information is currently not encoded in the hdf5 file (listed as TODO)
    #so here we just hack round it for the time being (use the params in 
    return {'Ne':1e4, 'mutation_rate':2e-8, 'recombination_rate':2e-8}
    
def msprime_to_ARGweaver_in(ts, ARGweaver_filehandle):
    """
    Takes an msprime TreeSequence object, and outputs a file in .sites format, suitable for input into ARGweaver
    (see http://mdrasmus.github.io/argweaver/doc/#sec-file-sites)
    Note that the documentation (http://mdrasmus.github.io/argweaver/doc/#sec-prog-arg-sample) states that the only mutation
    model is Jukes-Cantor (i.e. equal mutation between all bases). Assuming adjacent sites are treated independently, we
    convert variant format (0,1) to sequence format (A, T, G, C) by simply converting 0->A and 1->T
    
    Also note that the msprime simulations assume infinite sites by allowing mutations to occur at floating-point
    positions long a sequence. These need to be converted to integer positions along the genome for use in ARGweaver
    
    """
    from math import ceil
    import numpy as np
    simple_ts = ts.simplify()
    print("\t".join(["NAMES"]+[str(x) for x in range(simple_ts.get_sample_size())]), file=ARGweaver_filehandle)
    print("\t".join(["REGION", "chr", "1", str(int(simple_ts.get_sequence_length()))]), file=ARGweaver_filehandle)
    genotypes = None
    pos = 0
    for j, v in enumerate(simple_ts.variants()):
        if int(ceil(v.position)) != pos:
            #this is a new position. Print the genotype at the old position, and then reset everything
            if pos:
                print(pos, "".join(np.where(genotypes==0,"A","T")), sep="\t", file=ARGweaver_filehandle)
            genotypes = v.genotypes
            pos = int(ceil(v.position))
        else:
            genotypes =  np.logical_and(genotypes, v.genotypes)
    if pos:
        print(pos, "".join(np.where(genotypes==0,"A","T")), sep="\t", file=ARGweaver_filehandle)

    ARGweaver_filehandle.flush()
    ARGweaver_filehandle.seek(0)

def run_ARGweaver(Ne, mut_rate, recomb_rate, executable, ARGweaver_in_filehandle, ARGweaver_out_dir, out_prefix="aw", seed=None, iterations=None, sample_step=None, status_to=sys.stdout, quiet=False, rand_seed=None):
    import os
    from subprocess import call
    ARGweaver_out_dir = os.path.join(ARGweaver_out_dir,out_prefix)
    exe = [str(executable), '--output', ARGweaver_out_dir, '--popsize', str(Ne), '--mutrate', str(mut_rate), '--recombrate', str(recomb_rate), '--sites', ARGweaver_in_filehandle.name, '--overwrite']
    if quiet:
        exe.extend(['--quiet'])
    if rand_seed is not None:
        exe.extend(['--randseed', str(rand_seed)])
    if iterations is not None:
        exe.extend(['--iters', str(iterations)])
    if sample_step is not None:
        exe.extend(['--sample-step', str(sample_step)])
    if status_to:
        print("== Running ARGweaver as `{}` ==".format(" ".join(exe)), file=status_to)
    call(exe)
    
    
def ARGweaver_smc_to_msprime_txts(smc2bin_executable, prefix, tree_filehandle, status_to=sys.stdout):
    """
    convert the ARGweaver smc representation to coalescence records format
    """
    assert False, "smc2arg is currently broken and should not be used. See https://github.com/mdrasmus/argweaver/issues/20"
    from subprocess import call
    if status_to:
        print("== Converting the ARGweaver smc output file '{}' to .arg format using '{}' ==".format(prefix + ".smc.gz", smc2bin_executable), file=status_to)
    call([smc2bin_executable, prefix + ".smc.gz", prefix + ".arg"])
    with open(prefix + ".arg", "r+") as arg_fh:
        return(ARGweaver_arg_to_msprime_txts(arg_fh, tree_filehandle, status_to))

def ARGweaver_arg_to_msprime_txts(ARGweaver_arg_filehandle, tree_filehandle, status_to=sys.stdout):
    """
    convert the ARGweaver arg representation to coalescence records format
    
    We need to split ARGweaver records that extend over the whole genome into sections that cover just that 
    coalescence point.
    
    returns the mapping of ARGweaver node names to msprime node names
    """
    import re
    import csv
    if status_to:
        print("== Converting .arg output to msprime ==", file=status_to)
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
                    tree_filehandle.write("{}\t{}\t{}\t".format(leftbreak, rightbreak, node_names[node_name]))
                    tree_filehandle.write(",".join([str(x) for x in sorted([node_names[cnode] for cnode, cspan in children.items() if cspan[0]<rightbreak and cspan[1]>leftbreak])]))
                    tree_filehandle.write("\t{}\t{}\n".format(ARG_node_times[node_name], 0))
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
    
    tree_filehandle.flush()
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

def msprime_txts_to_hdf5(tree_filehandle, hdf5_outname=None, status_to=sys.stdout):
    from warnings import warn
    import shutil
    import msprime
    if status_to:
        print("== Converting new msprime ARG as hdf5 ===", file = status_to)
    try:
        ts = msprime.load_txt(tree_filehandle.name)
    except:
        warn("Can't load the txt file properly. Saved a copy to 'bad.msprime' for inspection")
        shutil.copyfile(tree_filehandle.name, "bad.msprime")
        raise
    if status_to:
        print("== loaded {} ===".format(tree_filehandle.name), file=status_to)
    try:
        simple_ts = ts.simplify()
    except:
        ts.dump("bad.hdf5")
        warn("Can't simplify. HDF5 file dumped to 'bad.hdf5'")
        raise
    if hdf5_outname:
        simple_ts.dump(hdf5_outname)
    return(simple_ts)


def main(directory, hdf5_filehandle):
    import os
    from dendropy import TreeList, calculate
    full_prefix = os.path.join(directory, os.path.splitext(os.path.basename(hdf5_filehandle.name))[0])
    with open(full_prefix+".sites", "w+") as aw_in, open(full_prefix+".msp_recs", "w+") as ms_rec, open(full_prefix+".msp_muts", "w+") as ms_mut:
        params = msprime_hdf5_to_ARGweaver_in(hdf5_filehandle, aw_in)
        ts = msprime.load(hdf5_filehandle.name)
        ts.write_records(ms_rec)
        ts.write_mutations(ms_mut)
        run_ARGweaver(Ne=params['Ne'], 
                      mut_rate=params['mutation_rate'],
                      recomb_rate=params['recombination_rate'],
                      executable= os.path.join(args.ARGweaver_executable_dir, args.ARGweaver_sample_executable),
                      ARGweaver_in_filehandle=aw_in,
                      iterations=20,
                      rand_seed=1234,
                      ARGweaver_out_dir=directory)
        if os.stat(aw_in.name).st_size == 0:
            warn("Initial .sites file is empty")            
        for smc_file in os.listdir(directory):
            if smc_file.endswith(".smc.gz"):
                smc = os.path.join(directory,smc_file)
                with open(smc.replace(".smc.gz", ".msp_recs"), "w+") as tree, \
                     open(smc.replace(".smc.gz", ".msp_muts"), "w+") as muts, \
                     open(smc.replace(".smc.gz", ".msp_nex"), "w+") as msp_nex, \
                     open(smc.replace(".smc.gz", ".nex"), "w+") as smc_nex:
                    ARGweaver_smc_to_msprime_txts(os.path.join(args.ARGweaver_executable_dir, args.ARGweaver_smc2arg_executable), smc.replace(".smc.gz", ""), tree)
                    ts = msprime_txts_to_hdf5(tree, full_prefix + ".hdf5")
                    ARGweaver_smc_to_nexus(smc, smc_nex, zero_based_tip_numbers=False)
                    print("#NEXUS\nBEGIN TREES;\nTRANSLATE\n{};".format(",\n".join(["{} {}".format(i+1,i+1) for i in range(ts.get_sample_size())])), file=msp_nex)
                    variant_index = 0
                    for (l, nwk) in ts.newick_trees():
                        variant_index += l
                        print("TREE "+str(variant_index)+" = "+nwk, file=msp_nex)
                    print("END;",file=msp_nex)
                    smc_nex.flush()
                    smc_nex.seek(0)
                    msp_nex.flush()
                    msp_nex.seek(0)
                    msp_trees = TreeList.get(file=msp_nex, schema='nexus') #zero_based_tip_numbers assumed False)
                    smc_trees = TreeList.get(file=smc_nex, schema="nexus", taxon_namespace=msp_trees.taxon_namespace)
                    #Check the smc trees against the msprime-imported equivalents
                    #NB, the ARGweaver output does not specify where mutations occur on the ARG, so we cannot
                    #reconstruct the sequences implied by this ARG for testing purposes, and thus cannot compare
                    #the original sequences with the reconstructed ones
                    assert len(smc_trees)==len(msp_trees), "number of trees in original and msprime-processed files is not the same"
                    print("The following RF and wRF statistics should be zero or near zero")
                    for (smc_tree, msp_tree) in zip(smc_trees, msp_trees):
                        print("Tree up to smc position {}, ms prime position {}: RF={} wRF={}".format(smc_tree._label, msp_tree._label, calculate.treecompare.symmetric_difference(smc_tree, msp_tree), calculate.treecompare.weighted_robinson_foulds_distance(smc_tree, msp_tree)))

if __name__ == "__main__":
    import argparse
    import filecmp
    import os
    from warnings import warn
    parser = argparse.ArgumentParser(description='Check ARGweaver imports by running an msprime simulation to create an ARGweaver import file, inferring some args from it in smc format, converting the .smc format to .arg format, reading the .arg into msprime, and comparing the nexus output trees with the trees in the .smc file. This testing process requires the dendropy library')
    parser.add_argument('hdf5file', type=argparse.FileType('r', encoding='UTF-8'), help='an msprime hdf5 file')
    parser.add_argument('--ARGweaver_executable_dir', '-d', default="../argweaver/bin/", help='the path to the directory containing the ARGweaver executables')
    parser.add_argument('--ARGweaver_sample_executable', '-x', default="arg-sample", help='the name of the ARGweaver executable')
    parser.add_argument('--ARGweaver_smc2arg_executable', '-s', default="smc2arg", help='the name of the ARGweaver executable')
    parser.add_argument('outputdir', nargs="?", default=None, help='the directory in which to store the intermediate files. If None, files are saved under temporary names')
    args = parser.parse_args()
    
    if args.outputdir == None:
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as aw_out_dir:
            print("Saving everything to temporary files (temporarily stored in {})".format(aw_out_dir))
            main(aw_out_dir, args.hdf5file)
    else:
        if not os.path.isdir(args.outputdir):
            warn("Output dir {} does not exist: creating it".format(args.outputdir))
            os.mkdir(args.outputdir)
        if len(os.listdir(args.outputdir)) > 0:
            warn("Output dir {} already contains files: deleting them".format(args.outputdir))
            import shutil
            shutil.rmtree(args.outputdir) 
            os.mkdir(args.outputdir)
        main(args.outputdir, args.hdf5file)
