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

def run_ARGweaver(Ne, mut_rate, recomb_rate, executable, ARGweaver_in_filehandle, ARGweaver_out_dir, out_prefix="aw", seed=None, iterations=100, status_to=sys.stdout, quiet=False):
    import os
    from subprocess import call
    ARGweaver_out_dir = os.path.join(ARGweaver_out_dir,out_prefix)
    if status_to:
        print("== Running ARGweaver ==", file=status_to)
    exe = [str(executable), '--output', ARGweaver_out_dir, '--popsize', str(Ne), '--mutrate', str(mut_rate), '--recombrate', str(recomb_rate), '--iters', str(iterations), '--sites']
    if quiet:
        exe.insert(1,'--quiet')
    call(exe + [ARGweaver_in_filehandle.name])
    
    
def ARGweaver_smc_to_msprime_txts(smc2bin_executable, prefix, tree_filehandle, status_to=sys.stdout):
    """
    convert the ARGweaver smc representation to coalescence records format
    """
    from subprocess import call
    if status_to:
        print("Converting the ARGweaver smc output file '{}' to .arg format using '{}'".format(prefix + ".smc.gz", smc2bin_executable), file=status_to)
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
            if node_name in temporary_marks:
                raise LookupError('ARG is not acyclic!')
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
            shutil.copyfile(ARGweaver_arg_filehandle.name, "bad.arg")
            raise type(e)(str(e) + "\n" +
                      "A copy of the file '{}' is being saved to 'bad.arg' for debugging".format(ARGweaver_arg_filehandle.name)).with_traceback(sys.exc_info()[2])
    tree_filehandle.flush()
    
    return(node_names)
    
   
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
    full_prefix = os.path.join(directory, os.path.splitext(os.path.basename(hdf5_filehandle.name))[0])
    with open(full_prefix+".sites", "w+") as aw_in:
        params = msprime_hdf5_to_ARGweaver_in(hdf5_filehandle, aw_in)
        run_ARGweaver(Ne=params['Ne'], 
                      mut_rate=params['mutation_rate'],
                      recomb_rate=params['recombination_rate'],
                      executable= os.path.join(args.ARGweaver_executable_dir, args.ARGweaver_sample_executable),
                      ARGweaver_in_filehandle=aw_in,
                      ARGweaver_out_dir=directory)
        if os.stat(aw_in.name).st_size == 0:
            warn("Initial .sites file is empty")            
        for smc_file in os.listdir(directory):
            if smc_file.endswith(".smc.gz"):
                prefix = os.path.join(directory,smc_file.replace(".smc.gz", ""))
                print(prefix)
                with open(prefix + ".msprime", "w+") as tree, \
                     open(prefix + ".msmuts", "w+") as muts:
                    ARGweaver_smc_to_msprime_txts(os.path.join(args.ARGweaver_executable_dir, args.ARGweaver_smc2arg_executable), prefix, tree)
                    msprime_txts_to_hdf5(tree, prefix + ".hdf5")
                    #NB, the ARGweaver output does not specify where mutations occur on the ARG, so we cannot
                    #reconstruct the sequences implied by this ARG for testing purposes, and thus cannot compare
                    #the original sequences with the reconstructed ones
        

if __name__ == "__main__":
    import argparse
    import filecmp
    import os
    from warnings import warn
    parser = argparse.ArgumentParser(description='Check ARGweaver imports by running msprime simulation through it and back out')
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
