import requests
import re
import json
import math
import os
import sys
import base64

import subprocess

import numpy as np

import msprime
sys.path.insert(1,os.path.join(sys.path[0],'..','tsinfer'))
import tsinfer

#tsinfer_executable = os.path.join(sys.path[0],'..','src','run_old_tsinfer.py') #use with e.g. `git checkout b1fa4ed8`
tsinfer_executable = os.path.join(sys.path[0],'..','src','run_tsinfer.py')

print("Running using")
subprocess.call(["python3", tsinfer_executable, "--version"])

if 'run_old_tsinfer' in tsinfer_executable:
    import dbm
    #pre-release tsinfer versions don't have sensible version numbers set
    #and have to be called with SampleData.initialise
    #monkey patch so that old versions of tsinfer (which use add_variant not add_site) work.
    #This can be deleted unless you want to use old pre-release (e.g. Feb 2018) versions of tsinfer
    tsinfer.SampleData.add_site = lambda a,b,c: tsinfer.SampleData.add_variant(a,c,b)
    tsinfer.SampleData.sites_position = tsinfer.SampleData.position
    NoSuchTsinferFileError = dbm.error
    def make_sample_data(path, num_samples):
        return tsinfer.SampleData.initialise(
            filename=path, sequence_length=34025983, num_samples = num_samples)
else:
    NoSuchTsinferFileError = tsinfer.exceptions.FileFormatError
    def make_sample_data(path, num_samples):
        return tsinfer.SampleData(path=path)



def treestring(name, tree):
    return "TREE " + name + " = [&R] " + tree.newick(precision=14)[:-1] + ":0;\n"

def header(n_tips, node_labels):
    #convert back to 0-based tip labels, or use string node labels (escaped with single quotes)
    tip_map = [
        str(i+1) + " " + ("'{}'".format(node_labels[i].replace("'","''")) if i in node_labels else str(i))
        for i in range(n_tips)]
    return "#NEXUS\nBEGIN TREES;\nTRANSLATE\n{};\n".format(",\n".join(tip_map))
        
def footer():
    return "END;\n"

def write_nexus_tree(tree, treefile, node_labels={}):
    """
    Writes out a single SparseTree from a tree sequence to a nexus tree file
    which allows tips to be relabelled so that they either correspond to the
    tskit numbers (0 based) or to a user specified set of names
    """
    treefile.write(header(tree.num_samples(), node_labels))
    treefile.write(treestring(str(tree.interval[1]), tree))
    treefile.write(footer())
    

def save_nexus_tree(tree, fn, **kwargs):
    """
    Same as write_nexus_tree() only use a file name not a file handle
    """
    with open(fn, "w+") as out:
        write_nexus_tree(tree, out, **kwargs)


focal_10 = set([
    'Oceania_PapuaNewGuinea.2451.2.SGDP',
    'WestEurasia_Finland.3014.2.SGDP',
    'MSL.2378.1.1000G',
    'ASM.Europe_EasternEurope.303.1.EGDP',
    'Africa_Congo.107.1.SGDP',
    'CHB.1027.2.1000G',
    'FIN.1452.1.1000G',
    'ITU.1971.1.1000G',
    'YRI.3179.1.1000G',
    'CHB.1014.1.1000G'])

google_files = {
    "filtered": "1jjRkDBO2QmAmsCZ6cgMjMY6g17pecp2G",
    "full"    : "1NZAp11dUnfq9VmXs4gN5qZ3d6wJiHt-Y", 
    "allsnps" : "1AF4lZHOO2IrtZJVZGPe4QaEKlk6ttq7B",
    }
for dataset, doc in google_files.items():
    a=[]
    subsample=set()
    url = "https://drive.google.com/uc?export=download&id={}".format(doc)
    chr = None
    l = 0
    r = requests.get(url, stream=True)
    l_names = []
    s_names = []
    for line in r.iter_lines(decode_unicode=True):
        # filter out keep-alive new lines
        if line:
            l += 1
            data = line.split() #split on whitespace
            if l==1: #header line
                lsd_fn = "{}.samples".format(dataset)
                ssd_fn = "{}_small.samples".format(dataset)
                large_sample_data = make_sample_data(lsd_fn, num_samples = len(data[1:]))
                small_sample_data = make_sample_data(ssd_fn, num_samples = len(focal_10))
                for i, name in enumerate([x.strip("\"") for x in data[1:]]): #first is position
                    if 'run_old_tsinfer' in tsinfer_executable:
                        l_names.append(name)
                    else:
                        #the data is present as haploid samples
                        large_sample_data.add_individual(ploidy=1, metadata={'name':name})
                    if name in focal_10:
                        subsample.add(i)
                        if 'run_old_tsinfer' in tsinfer_executable:
                            s_names.append(name)
                        else:
                            small_sample_data.add_individual(ploidy=1, metadata={'name':name})

            else:
                name = data[0].strip("\"")
                #match e.g. chr20:33896756-33896757
                snp_position = re.match(r'([^:]+):(\d+)-(\d+)', data[1])
                #check we have a sensible position
                assert snp_position, "SNP position is {}".format(data[1])
                #check all are on the same chromosome
                assert (chr == None) or (chr == snp_position.group(1)), \
                    "Different chromosomes, {} vs {}".format(chr, snp_position.group(1))
                chr = snp_position.group(1)
                #check these are single SNPs
                assert int(snp_position.group(2))+1 == int(snp_position.group(3)), \
                    "start+1 != end ({} vs {})".format(snp_position.group(2), snp_position.group(3))
                #read in
                a.append(int(snp_position.group(2)))
                if all([0<=int(i)<=1 for i in data[2:]]):
                    large_sample_data.add_site(int(snp_position.group(2)), [int(x) for x in data[2:]], ["0", "1"])
                    small_sample_data.add_site(int(snp_position.group(2)), [int(x) for i, x in enumerate(data[2:]) if i in subsample], ["0", "1"])
                else:
                    print("Problem in line {} of {} (SNP {} {})".format(l, doc, data[0], data[1]))
    large_sample_data.finalise()
    small_sample_data.finalise()
    
    extra_params = ["--recombination-rate", str(1e-8), "--length", str(34025983)]
    l_inferred_ts = subprocess.call(["python3", tsinfer_executable, lsd_fn, "{}.trees".format(dataset)] + extra_params)
    s_inferred_ts = subprocess.call(["python3", tsinfer_executable, ssd_fn, "{}_small.trees".format(dataset)] + extra_params)
    l_inferred_ts = msprime.load("{}.trees".format(dataset))
    s_inferred_ts = msprime.load("{}_small.trees".format(dataset))
    print("For {} dataset, pos {} to {}  ({} SNPs):  inferred {} trees using all data & {} trees using {} samples".format(
        dataset, min(a), max(a), l-1, l_inferred_ts.num_trees, s_inferred_ts.num_trees,len(focal_10)))
    
    if l_names:
        node_labels={k:v for k,v in enumerate(l_names)}
    else:
        node_labels={n.id:json.loads(l_inferred_ts.individual(n.individual).metadata.decode('utf-8'))['name'] for n in l_inferred_ts.nodes() if n.is_sample()}
    for tree in l_inferred_ts.trees():
        low, high = tree.interval
        if low <= 33952619 < high:
            #only output the tree at the locus of interest
            identifier = "{}-{}_{}".format(low, high, dataset)
            svg = tree.draw(path=identifier + ".svg", format="svg", width = 1000*math.log(l_inferred_ts.num_samples), height=1000, node_labels=node_labels)
            save_nexus_tree(tree, identifier + ".nex", node_labels=node_labels)

    #take the huge inferred tree and subsample down to only the 10 selected tips
    #NB - a hack in the simplify() method means we need to reconstruct the 
    tmp_tables = l_inferred_ts.dump_tables()
    tmp_tables.nodes.clear()
    tmp_tables.individuals.clear()
    for node in l_inferred_ts.nodes():
        if node.is_sample():
            metadata = l_inferred_ts.individual(node.individual).metadata
        else:
            metadata = None
        tmp_tables.nodes.add_row(
            flags=node.flags, time=node.time, metadata=metadata)
    tmp_tables.nodes.set_columns(
       flags=tmp_tables.nodes.flags,
       time=tmp_tables.nodes.time,
       metadata=tmp_tables.nodes.metadata,
       metadata_offset=tmp_tables.nodes.metadata_offset)
    ls_inferred_ts = tmp_tables.tree_sequence().simplify(sorted([s for s in subsample]))
    
    if l_names:
        node_labels={k:v for k,v in enumerate(s_names)}
    else:
        node_labels={n.id:json.loads(n.metadata.decode('utf-8'))['name'] for n in ls_inferred_ts.nodes() if n.is_sample()}
    for tree in ls_inferred_ts.trees():
        low, high = tree.interval
        if low <= 33952619 < high:
            #only output the tree at the locus of interest
            identifier = "{}-{}_{}".format(low, high, dataset)
            svg = tree.draw(path=identifier + "_subs.svg", format="svg", width = 1000*math.log(ls_inferred_ts.num_samples), height=1000, node_labels=node_labels)
            save_nexus_tree(tree, identifier + "_subs.nex", node_labels=node_labels)

    #look at the inference over only 10 inferred tips
    if l_names:
        node_labels={k:v for k,v in enumerate(s_names)}
    else:
        node_labels={n.id:json.loads(s_inferred_ts.individual(n.individual).metadata.decode('utf-8'))['name'] for n in s_inferred_ts.nodes() if n.is_sample()}
    for tree in s_inferred_ts.trees():
        low, high = tree.interval
        if low <= 33952619 < high:
            #only output the tree at the locus of interest
            identifier = "{}-{}_{}".format(low, high, dataset)
            svg = tree.draw(path=identifier + "_small.svg", format="svg", width = 1000*math.log(s_inferred_ts.num_samples), height=1000, node_labels=node_labels)
            save_nexus_tree(tree, identifier + "_small.nex", node_labels=node_labels)
            
"""
#R code for analysis

library(treespace)
library(phytools)

names <- c("filtered", "full", "allsnps")

RaxML.tree = read.newick(text="(((Africa_Congo.107.1.SGDP:0.5345217134,((CHB.1014.1.1000G:0.0428432624,CHB.1027.2.1000G:0.2913519629):0.1366214608,ASM.Europe_EasternEurope.303.1.EGDP:0.05531952005):0.2111176833):1.168343643,(WestEurasia_Finland.3014.2.SGDP:0.09731705206,((MSL.2378.1.1000G:0.1389501155,(ITU.1971.1.1000G:0.3892374349,FIN.1452.1.1000G:0.1465821786):2.276752947):0.02793908635,YRI.3179.1.1000G:0.6379600936):0.7898707392):1.50000075e-05):0.1321460672,Oceania_PapuaNewGuinea.2451.2.SGDP:0.0546869088);")

#Check tips
n1 = read.nexus(Sys.glob(sprintf("*_%s_small.nex", "full")))
tip.match = sapply(n1$tip.label, function(n) {n %in% RaxML.tree$tip.label})
if (all(tip.match)) {
    sapply(names, function(n) {
        c(origsize10_with_poly  = treeDist(RaxML.tree, read.nexus(Sys.glob(sprintf("*_%s_small.nex", n)))),
          origsize10_break_poly = mean(replicate(50, treeDist(RaxML.tree, multi2di(read.nexus(Sys.glob(sprintf("*_%s_small.nex", n))))))),
          subsamp6k_with_poly   = treeDist(RaxML.tree, read.nexus(Sys.glob(sprintf("*_%s_subs.nex", n)))),
          subsamp6k_break_poly  = mean(replicate(50, treeDist(RaxML.tree, multi2di(read.nexus(Sys.glob(sprintf("*_%s_subs.nex", n)))))))
        )
    })
} else {
    tip.match
}


"""

"""
Gives the following for the new (MASTER) codebase (draconian stopping rule)

filtered      full  allsnps                       
7.810250  7.810250 11.09054  small_with_polytomies  
7.297314  7.217913 10.80740  small_break_polytomies 
7.745967 11.832160 11.48913  subs_with_polytomies   
7.657892 11.815536 11.15343  subs_break_polytomies  


And for the old codebase (very liberal stopping rule, nearly all ancestors v long):

filtered     full  allsnps                         
7.615773 7.615773 9.486833  small_with_polytomies  
7.136824 7.225302 9.812799  small_break_polytomies 
7.071068 7.071068 9.746794  subs_with_polytomies   
7.071068 7.071068 9.756919  subs_break_polytomies  


And for Jerome's recent modified ancestor construction in 
jeromekelleher/update-ancestors-alg

#revision c8f3abf (less draconian stopping rule than MASTER)
                       filtered     full  allsnps
small_with_polytomies  7.000000 7.615773 9.949874
small_break_polytomies 6.719144 7.181756 9.961036
subs_with_polytomies   7.071068 7.483315 9.643651
subs_break_polytomies  6.799543 7.342762 9.916682


#revision 00ebf1d (less draconian stopping rule, allow for some error, no trailing zeros)
                       filtered     full   allsnps
small_with_polytomies  7.549834 7.745967 10.000000
small_break_polytomies 7.549834 7.745967 10.000000
subs_with_polytomies   7.280110 6.855655  9.539392
subs_break_polytomies  7.171893 6.661915  9.866577
"""
