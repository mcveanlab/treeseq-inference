#!/usr/bin/env python3
"""Various functions to convert msprime simulation file to ARGweaver input format, and from .arg files to msprime input.

When run as a script, takes an msprime simulation in hdf5 format, saves to ARGweaver input format (haplotype sequences), runs ARGweaver inference on it, converts the ARGweaver output file (arg and mutations) to msprime input files, reads these into msprime outputs the haplotype sequences again, and checks that the haplotype sequences are the same.

E.g. ./msprime_ARGweaver.py ../test_files/4000.hdf5 -x ../argweaver/bin/****

"""
import sys
import os.path
sys.path.insert(1,os.path.join(sys.path[0],'..','msprime')) # use the local copy of msprime in preference to the global one
import msprime
from warnings import warn

def msprime_hdf5_to_ARGweaver_in(msprime_hdf5, ARGweaver_filehandle):
    print("== Saving to ARGweaver input format ==")
    ts = msprime.load(msprime_hdf5.name)
    msprime_to_ARGweaver_in(ts, ARGweaver_filehandle)
    
def msprime_to_ARGweaver_in(ts, ARGweaver_filehandle):
    simple_ts = ts.simplify()
    for j, v in enumerate(simple_ts.variants(as_bytes=True)):
        print(j, v.genotypes.decode(), sep="\t", file=ARGweaver_filehandle)
    ARGweaver_filehandle.flush()
    ARGweaver_filehandle.seek(0)

def run_ARGweaver(executable, ARGweaver_in_filehandle, ARGweaver_out_filehandle, seed=None):
    print("== Running ARGweaver ==")
    from subprocess import call
    exe = [str(executable), 'build']
    if seed:
        exe += ['-s', str(int(seed))]
    call(exe + [fastARG_in_filehandle.name], stdout=fastARG_out_filehandle)
    fastARG_out_filehandle.flush()
    fastARG_out_filehandle.seek(0)
    
    
def ARGweaver_arg_to_msprime_txt(ARGweaver_arg_filehandle, tree_filehandle):
    """
    convert the ARGweaver arg representation to coalescence records format
    
    We need to split ARGweaver records that extend over the whole genome into sections that cover just that 
    coalescence point.
    
    returns the mapping of ARGweaver node names to msprime node names
    """
    import re
    import csv
    print("== Converting .arg output to msprime ==")
    ARG_nodes={} #cr[X] = child1:[left,right], child2:[left,right],... : serves as intermediate ARG storage 
    ARG_node_times={} #node_name => time
    node_names={} #map of ARGweaver names -> numbers
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
    
    #now relabel the nodes
    n_tips = len(node_names)
    for key in ARG_nodes:
        node_names[key]=int(key)+n_tips
    
    
    #recursive hack to make times strictly decreasing
    ARG_node_time_epsilon= {name:0 for name in ARG_node_times}
    def set_child_times(node_name, epsilon=1):
        ARG_node_time_epsilon[node_name] += epsilon
        for child_name in ARG_nodes[node_name]:
            try:
                set_child_times(child_name, epsilon+1)
            except:
                print(child_name)
                pass    
    set_child_times(root_node)
    max_epsilon = max(ARG_node_time_epsilon.values())
    print(max_epsilon, ARG_node_time_epsilon)
    for nm in ARG_node_times:
        ARG_node_times[nm] += (max_epsilon - ARG_node_time_epsilon[nm])/10000
    
    
    for node_name in sorted(ARG_node_times, key=ARG_node_times.get): #sort by time
        #look at the break points for all the child sequences, and break up into that number of records
        try:
            children = ARG_nodes[node_name]
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
            assert node_name.startswith('n'), "The node {} is not a parent of any other node".format(node_name)
            
    tree_filehandle.flush()
    
    return(node_names)
    
   
def fastARG_root_seq(fastARG_out_filehandle):
    for line in fastARG_out_filehandle:
        spl = line.split(None,3)
        if spl[0]=="S":
            root_seq = spl[2]
            fastARG_out_filehandle.seek(0)
            return([seq!="0" for seq in root_seq])
    warn("No root sequence found in '{}'".format(fastARG_out_filehandle.name))
    fastARG_out_filehandle.seek(0)
    return([])

def msprime_txts_to_fastARG_in_revised(tree_filehandle, mutations_filehandle, root_seq, fastARG_filehandle, hdf5_outname=None):
    if hdf5_outname:
        print("== Saving new msprime ARG as hdf5 and also as input format for fastARG ==")
    else:
        print("== Saving new msprime ARG as input format for fastARG ==")
    
    ts = msprime.load_txt(tree_filehandle.name, mutations_filehandle.name)
    try:
        simple_ts = ts.simplify()
    except:
        ts.dump("bad.hdf5")
        raise
    if hdf5_outname:
        simple_ts.dump(hdf5_outname)
    for j, v in enumerate(simple_ts.variants(as_bytes=True)):
        if root_seq[j]:
            print(j, v.genotypes.decode().translate(str.maketrans('01','10')), sep="\t", file=fastARG_filehandle)
        else:
            print(j, v.genotypes.decode(), sep="\t", file=fastARG_filehandle)
    fastARG_filehandle.flush()
    fastARG_filehandle.seek(0)
    return(simple_ts)

if __name__ == "__main__":
    import argparse
    import filecmp
    import os
    parser = argparse.ArgumentParser(description='Check fastARG imports by running msprime simulation through it and back out')
    parser.add_argument('hdf5file', type=argparse.FileType('r', encoding='UTF-8'), help='an msprime hdf5 file')
    parser.add_argument('--fastARG_executable', '-x', default="fastARG", help='the path & name of the fastARG executable')
    parser.add_argument('outputdir', nargs="?", default=None, help='the directory in which to store the intermediate files. If None, files are saved under temporary names')
    args = parser.parse_args()
    
    
    if args.outputdir == None:
        print("Saving everything to temporary files")
        from tempfile import NamedTemporaryFile
        with NamedTemporaryFile("w+") as fa_in, \
             NamedTemporaryFile("w+") as fa_out, \
             NamedTemporaryFile("w+") as tree, \
             NamedTemporaryFile("w+") as muts, \
             NamedTemporaryFile("w+") as fa_revised: 
    
            msprime_hdf5_to_fastARG_in(args.hdf5file, fa_in)
            run_fastARG(args.fastARG_executable, fa_in, fa_out)
            root_seq = fastARG_root_seq(fa_out)
            fastARG_out_to_msprime_txts(fa_out, tree, muts)
            msprime_txts_to_fastARG_in_revised(tree, muts, root_seq, fa_revised)
            if os.stat(fa_in.name).st_size == 0:
                warn("Initial fastARG input file is empty")            
            elif filecmp.cmp(fa_in.name, fa_revised.name, shallow=False) == False:
                warn("Initial fastARG input file differs from processed fastARG file")
            else:
                print("Conversion via fastARG has worked! Input and output files are identical")
    else:
        prefix = os.path.splitext(os.path.basename(args.hdf5file.name))[0]
        full_prefix = os.path.join(args.outputdir, prefix)
        with open(full_prefix + ".haps", "w+") as fa_in, \
             open(full_prefix + ".fastarg", "w+") as fa_out, \
             open(full_prefix + ".msprime", "w+") as tree, \
             open(full_prefix + ".msmuts", "w+") as muts, \
             open(full_prefix + ".haps_revised", "w+") as fa_revised:
            msprime_hdf5_to_fastARG_in(args.hdf5file, fa_in)
            run_fastARG(args.fastARG_executable, fa_in, fa_out)
            root_seq = fastARG_root_seq(fa_out)
            fastARG_out_to_msprime_txts(fa_out, tree, muts)
            msprime_txts_to_fastARG_in_revised(tree, muts, root_seq, fa_revised, full_prefix + ".hdf5_revised")
            if os.stat(fa_in.name).st_size == 0:
                warn("Initial fastARG input file is empty")            
            elif filecmp.cmp(fa_in.name, fa_revised.name, shallow=False) == False:
                warn("Initial fastARG input file differs from processed fastARG file")
            else:
                print("Conversion via fastARG has worked! Input and output files are identical")
