#!/usr/bin/env python3
"""
Various functions to convert msprime simulation file to run in Heng Li's fastARG program 
and back again.

When run as a script, takes an msprime simulation in .trees format, saves to fastARG
input format (haplotype sequences), runs fastARG on it, converts fastARG output to
msprime input (2 files), reads these into msprime outputs the haplotype sequences again,
and checks that they agree.

E.g. ./ts_fastARG.py ../test_files/4000.trees -x ../fastARG/fastARG

"""
import os
import sys
import logging
import tempfile

import numpy as np

import msprime

def ts_to_fastARG_in(ts, fastARG_filehandle):
    # There is an odd intermittent bug in FastARG which fails unpredictably when there
    # is a trailing newline in the file, so print newlines *before* all but 1st line
    line_prefix = ""
    for v in ts.variants(as_bytes=True):
        #do not print out non-variable sites
        if len(v.alleles) > 1:
            print(
                line_prefix + v.site.position, v.genotypes.decode(), 
                sep="\t", end="", file=fastARG_filehandle)
            line_prefix = "\n"
    fastARG_filehandle.flush()
    fastARG_filehandle.seek(0)

def samples_to_fastARG_in(sample_data, fastARG_filehandle):
    # There is an odd intermittent bug in FastARG which fails unpredictably when there
    # is a trailing newline in the file, so print newlines *before* all but 1st line
    line_prefix = ""
    position = sample_data.sites_position[:] #decompress all in one go to avoid sequential unpacking
    for id, genotypes in sample_data.genotypes():
        #do not print out non-variable sites
        if np.any(genotypes) and not np.all(genotypes):
            s = (genotypes + ord('0')).tobytes().decode() #convert to string
            print(
                line_prefix + str(position[id]), s, 
                sep= "\t", end="", file=fastARG_filehandle)
            line_prefix = "\n"
    fastARG_filehandle.flush()
    fastARG_filehandle.seek(0)

def get_cmd(executable, fastARG_in, seed):
    files = [fastARG_in.name if hasattr(fastARG_in, "name") else fastARG_in]
    exe = [str(executable), 'build'] + files
    if seed is not None:
        exe += ['-s', str(int(seed))]
    return exe

def variant_positions_from_fastARGin_name(fastARG_in):
    with open(fastARG_in) as fastARG_in_filehandle:
        return(variant_positions_from_fastARGin(fastARG_in_filehandle))

def variant_positions_from_fastARGin(fastARG_in_filehandle):
    import numpy as np
    fastARG_in_filehandle.seek(0) #make sure we reset to the start of the infile
    vp=[]
    for line in fastARG_in_filehandle:
        try:
            pos = line.split()[0]
            vp.append(float(pos))
        except:
            logging.warning("Could not convert the title on the following line to a floating point value:\n {}".format(line))
    return(np.array(vp))

def fastARG_out_to_ts_txts(
        fastARG_out_filehandle, variant_positions,
        nodes_filehandle, edges_filehandle, sites_filehandle, mutations_filehandle,
        seq_len=None):
    """
    convert the fastARG output format (plus a list of positions) to 4 text
    files (edges, nodes, mutations, and sites) which can be read in to a tree sequence.
    We need to split fastARG records that extend over the whole genome into
    sections that cover just that coalescence point.
    """
    import csv
    import numpy as np
    logging.debug("== Converting fastARG output to ts ==")
    fastARG_out_filehandle.seek(0) #make sure we reset to the start of the infile
    ARG_nodes={} #An[X] = child1:[left,right], child2:[left,right],... : serves as intermediate ARG storage
    mutations={}
    mutation_nodes=set() #to check there aren't duplicate nodes with the same mutations - a fastARG restriction
    seq_len = seq_len if seq_len is not None else  variant_positions[-1]+1 #if not given, hack seq length to max variant pos +1
    try:
        breakpoints = np.concatenate([[0],np.diff(variant_positions)/2 + variant_positions[:-1], [seq_len]])
    except IndexError:
        raise ValueError(
            "Some variant positions seem to lie outside the sequence length "
            "(l={}):\n{}".format(seq_len, variant_positions))
    for line_num, fields in enumerate(csv.reader(fastARG_out_filehandle, delimiter='\t')):
        if   fields[0]=='E':
            srand48seed = int(fields[1])

        elif fields[0]=='N':
            haplotypes, loci = int(fields[1]), int(fields[2])
            csv.field_size_limit(max(csv.field_size_limit(), loci)) #needed because the last line has a field as long as loci


        elif fields[0]=='C' or fields[0]=='R':
            #coalescence or recombination events -
            #coalescence events should come in pairs
            #format is curr_node, child_node, seq_start_inclusive, seq_end_noninclusive, num_mutations, mut_loc1, mut_loc2, ....
            curr_node, child_node, left, right, n_mutations = [int(i) for i in fields[1:6]]
            if curr_node<child_node:
                raise ValueError(
                    "Line {} has an ancestor node_id less than the id of its "
                    "children, so node ID cannot be used as a proxy for age"
                    .format(line_num))
            #one record for each node
            if curr_node not in ARG_nodes:
                ARG_nodes[curr_node] = {}
            if child_node in ARG_nodes[curr_node]:
               raise ValueError(
                   "Child node {} already exists for node {} (line {})"
                   .format(child_node, curr_node, line_num))
            if child_node not in ARG_nodes:
                ARG_nodes[child_node]={}
            ARG_nodes[curr_node][child_node]=[breakpoints[left], breakpoints[right]]
            #mutations in tree seqs are placed as ancestral to a target node, rather than descending from a child node (as in fastARG)
            #we should check that (a) mutation at the same locus does not occur twice
            # (b) the same child node is not used more than once (NB this should be a fastARG restriction)
            if n_mutations:
                if child_node in mutation_nodes:
                    logging.warning("Node {} already has some mutations: more are being added from line {}.".format(child_node, line_num))
                mutation_nodes.add(child_node)
                for pos in fields[6:(6+n_mutations)]:
                    p = int(pos)
                    if p in mutations:
                        logging.warning("Duplicate mutations at a single locus: {}. One is from node {}.".format(m, child_node))
                    else:
                        mutations[p]=child_node


        elif fields[0]=='S':
            #sequence at root
            root_node = int(fields[1])
            base_sequence = fields[2]
        else:
            raise ValueError(
                "Bad line format - the fastARG file has a line that does not "
                "begin with E, N, C, R, or S:\n{}"
                .format("\t".join(fields)))

    #Done reading in

    if len([1 for k,v in ARG_nodes.items() if len(v)==0]) != haplotypes:
        raise ValueError(
            "We expect the same number of childless nodes as haplotypes "
            "but they are different ({} vs {})"
            .format(len([1 for k,v in ARG_nodes.items() if len(v)==0]), haplotypes))

    cols = {
        'nodes':     ["id","is_sample","time"],
        'edges':  ["left","right","parent","child"],
        'sites':     ["position","ancestral_state"],
        'mutations': ["site","node","derived_state"]
    }
    lines = {k:"\t".join(["{%s}" % s for s in v]) for k,v in cols.items()}
    for k,v in cols.items():
        #print the header lines
        print("\t".join(v), file=locals()[k+'_filehandle'])
    for node,children in sorted(ARG_nodes.items()): #sort by node number == time
        #look at the break points for all the child sequences, and break up into that number of records

        print(lines['nodes'].format(id=node, is_sample=int(len(children)==0), time=0 if len(children)==0 else node),
            file=nodes_filehandle)
        breaks = set()
        if len(children):
            for leftright in children.values():
                breaks.update(leftright)
            breaks = sorted(breaks)
            for i in range(1,len(breaks)):
                leftbreak = breaks[i-1]
                rightbreak = breaks[i]
                #NB - here we could try to input the values straight into an tree seq python structure,
                target_children = [str(cnode) for cnode, cspan in sorted(children.items()) if cspan[0]<rightbreak and cspan[1]>leftbreak]
                print(lines['edges'].format(left=leftbreak, right=rightbreak,
                    parent=node, child=",".join(target_children)), file=edges_filehandle)

    for i,base in enumerate(base_sequence):
        if base not in ['1','0']:
            raise ValueError('The ancestral sequence is not 0 or 1')
        print(lines['sites'].format(
            position=variant_positions[i], ancestral_state=base),
            file=sites_filehandle)

    for pos, node in sorted(mutations.items()):
        print(lines['mutations'].format(
            site=pos, node=node, derived_state=1-int(base_sequence[pos])),
            file=mutations_filehandle)

    for k in cols.keys():
        locals()[k+'_filehandle'].flush()
        locals()[k+'_filehandle'].seek(0)

def fastARG_out_to_ts(fastARG_out_filehandle, variant_positions, seq_len=None):
    """
    The same as fastARG_out_to_msprime_txts, but use temporary files and return a ts.
    """
    with tempfile.NamedTemporaryFile("w+") as nodes, \
        tempfile.NamedTemporaryFile("w+") as edges, \
        tempfile.NamedTemporaryFile("w+") as sites, \
        tempfile.NamedTemporaryFile("w+") as mutations:
        fastARG_out_to_ts_txts(fastARG_out_filehandle, variant_positions,
            nodes, edges, sites, mutations, seq_len=seq_len)
        ts = msprime.load_text(nodes=nodes, edges=edges, sites=sites, mutations=mutations).simplify()
        return ts


def main(ts, fastARG_executable, fa_in, fa_out, nodes_fh, edges_fh, sites_fh, muts_fh):
    """
    This is just to test if fastarg produces the same haplotypes
    """
    import subprocess
    seq_len = ts.get_sequence_length()
    ts_to_fastARG_in(ts, fa_in)
    subprocess.call([fastARG_executable, 'build', fa_in.name], stdout=fa_out)
    fastARG_out_to_ts_txts(fa_out, variant_positions_from_fastARGin(fa_in),
        nodes_fh, edges_fh, sites_fh, muts_fh, seq_len=seq_len)

    new_ts = msprime.load_text(nodes=nodes_fh, edges=edges_fh, sites=sites_fh, mutations=muts_fh)
    simple_ts = new_ts.simplify()
    logging.debug("Simplified num_records should always be < unsimplified num_records.\n"
        "For low mutationRate:recombinationRate ratio,"
        " the initial num records will probably be higher than the"
        " fastarg num_records, as the original simulation will have records"
        " which leave no mutational trace. As the mutation rate increases,"
        " we expect the fastarg num_records to equal, then exceed the original"
        " as fastarg starts inferring the wrong (less parsimonious) set of trees")
    logging.debug(
        "Initial num records = {}, fastARG (simplified) = {}, fastARG (unsimplified) = {}".format(
        ts.get_num_records(), simple_ts.get_num_records(), new_ts.get_num_records()))


if __name__ == "__main__":
    import argparse
    import os
    default_sim={'sample_size':20, 'Ne':1e4, 'length':10000, 'mutation_rate':2e-8, 'recombination_rate':2e-8}
    parser = argparse.ArgumentParser(description='Check fastARG imports by running msprime simulation through it and back out')
    parser.add_argument('--trees_file', '-f', type=argparse.FileType('r', encoding='UTF-8'), help='an msprime .trees file. If none, a simulation is created with {}'.format(default_sim))
    parser.add_argument('--fastARG_executable', '-x', default="fastARG", help='the path & name of the fastARG executable')
    parser.add_argument('--outputdir', '-o', help='the directory in which to store the intermediate files. If None, files are saved under temporary names')
    args = parser.parse_args()
    logging.basicConfig(
        format='%(asctime)s %(message)s', level=logging.DEBUG, stream=sys.stdout)
    if args.trees_file:
        ts = msprime.load(args.trees_file)
        prefix = os.path.splitext(os.path.basename(args.trees_file.name))[0]
    else:
        logging.debug("Creating new simulation with {}".format(default_sim))
        ts = msprime.simulate(**default_sim)
        prefix = "sim"
        logging.debug("mutationRate:recombinationRate ratio is {}".format(default_sim['mutation_rate']/default_sim['recombination_rate']))
    if args.outputdir == None:
        logging.info("Saving everything to temporary files")
        with tempfile.NamedTemporaryFile("w+") as fa_in, \
             tempfile.NamedTemporaryFile("w+") as fa_out, \
             tempfile.NamedTemporaryFile("w+") as nodes, \
             tempfile.NamedTemporaryFile("w+") as edges, \
             tempfile.NamedTemporaryFile("w+") as sites, \
             tempfile.NamedTemporaryFile("w+") as mutations:
            main(ts, args.fastARG_executable, fa_in, fa_out, nodes, edges, sites, mutations)
    else:
        full_prefix = os.path.join(args.outputdir, prefix)
        with open(full_prefix + ".haps", "w+") as fa_in, \
             open(full_prefix + ".fastarg", "w+") as fa_out, \
             open(full_prefix + ".TSnodes", "w+") as nodes, \
             open(full_prefix + ".TSedges", "w+") as edges, \
             open(full_prefix + ".TSsites", "w+") as sites, \
             open(full_prefix + ".TSmuts", "w+") as mutations:
            main(ts, args.fastARG_executable, fa_in, fa_out, nodes, edges, sites, mutations)
