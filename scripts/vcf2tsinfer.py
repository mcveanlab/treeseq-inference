#!/usr/bin/env python3
description="""
Take a vcf file with ancestral allele information (e.g. from 1000 genomes phase 1, such as
ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
and convert it to the huge samples _times_ sites array used in tsinfer. Store this
together with the variant positions and the sample names in hdf5 output files, one per
chromosome.

This script deals with samples with different ploidies (assuming the data is phased) by
adding '#a', '#b' etc to the sample names. For haploid samples, the '#x' suffix is removed.

Since we require complete data, we need to prune out mising values. At the moment we do
this by pruning samples with missing values, but we could equally do so by pruning sites
instead, or perhaps optimally remove both sites and samples (see 
https://stackoverflow.com/questions/48355644/optimally-remove-rows-or-columns-in-numpy-to-remove-missing-data)
"""

"""
NB - a slightly odd observation: all the alleles for which the ancestral state is not present in the dataset are SNPs

read using
with h5py.File(filename, 'r') as f:
    print(f['data']['variants'][()])
"""
import collections
import argparse
import string

import numpy as np
import h5py
import pysam

parser = argparse.ArgumentParser(description=description)
parser.add_argument('infile', 
                    help='a vcf or vcf.gz file')
parser.add_argument('outfile',
                    help='the output file prefix: data will be saved in hdf5 format as PREFIXchr.hdf5')
parser.add_argument('--only_use_n_variants', '-n', type=int, default=None,
                    help='For testing purposes, only use the first n variants')

args = parser.parse_args()

vcf_in = pysam.VariantFile(args.infile)

def process_variant(rec, rows, max_ploidy):
    """
    Return the true position of this variant, or None if we wish to omit
    Also fills out the site_data array
    """
    #restrict to cases where ancestral state contains only letters ATCG
    # i.e. where the ancestral state is certain (see ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README
    if "AA" in rec.info and all(letter in "ATCG" for letter in rec.info["AA"]):
        #only use biallelic variants 
        if len(rec.alleles) == 2:
            if rec.info["AA"] not in rec.alleles:
                print("Ancestral state {} not in allele list ({}) for position {}".format(\
                    rec.info["AA"], rec.alleles, rec.pos))
                pass
            else:
                omit=False
                #check for duplicate positions, often caused e.g. by C/CAT as opposed to -/AT
                #if the first x letters are the same, we have an intermediate position, e.g.
                #if C is at position 123, we can place the indel AT at position 123.5
                allele_start = 0
                for i in range(min([len(a) for a in rec.alleles])):
                    if len(set([a[i] for a in rec.alleles])) > 1:
                        #alleles differ here
                        break
                    allele_start += 1
                if allele_start != 0:
                    pos = rec.pos+allele_start
                    if len(set([len(a) for a in rec.alleles])) == 1:
                        #all alleles are the same length, => this is not an indel
                        print("The variants at {} share sequence, but are not an indel: {}".format({rec.id:rec.pos}, rec.alleles))
                    else:
                        pos-=0.5
                    if allele_start > 1:
                        print("The variants at {} share more than one starting letter: {}".format({rec.id:rec.pos}, rec.alleles))
                        
                    #print("Starting allele at an incremented position ({} not {}) for {} (alleles:{})".format(
                    #    pos, rec.pos, rec.id, rec.alleles))
                else:
                    allele_start=0
                    pos = rec.pos

                site_data = -np.ones((len(rows),), dtype="i1")
                for label, samp in rec.samples.items():
                    for i in range(min(max_ploidy, len(samp.alleles))):
                        #print("allele {} (pos {}, sample {}, ancestral state {}, alleles {})".format( \
                        #    rec.alleles[i], rec.pos, label+suffix, rec.info["AA"], rec.alleles))
                        if samp.alleles[i] not in rec.alleles:
                            continue
                        suffix = string.ascii_lowercase[i]
                        site_data[rows[label+'#'+suffix]] = samp.alleles[i][allele_start:]!=rec.info["AA"][allele_start:]
                return rec.chrom, pos, site_data
    return rec.chrom, None, None

allele_count = {}
rows, row = {}, 0
max_ploidy = 2
suffixes = string.ascii_lowercase[:max_ploidy] #use 'a' and 'b' for ploidy suffix

for sample_name in vcf_in.header.samples:
    #assume each sample can have a maximum of 2 variants (maternal & paternal)
    #i.e. max is diploid. For haploid samples, we just use the suffix 'a'
    #and hack to remove it afterwards
    for suffix in suffixes:
        rows[sample_name+'#'+suffix]=row
        row+=1
extend_amount = 10000 #numpy storage array will be extended by this much at a time
chromosomes = {} #most data stored here
output_freq_variants = 1e3 #output status after multiples of this many variants read


for j, variant in enumerate(vcf_in.fetch()):
    if j==args.only_use_n_variants:
        break
    #keep track of allele numbers
    allele_count[len(variant.alleles)] = allele_count.get(len(variant.alleles),0) + 1

    c, position, locus_data = process_variant(variant, rows, max_ploidy)
    #have we stored this chromosome yet
    if c not in chromosomes:
        chromosomes[c] = {
            'position':collections.OrderedDict(),
            'previous_position': None,
            'sites_by_samples': -np.ones((extend_amount, len(rows)), dtype="i1")}
    #check if we now need more storage
    if chromosomes[c]['sites_by_samples'].shape[0] <= len(chromosomes[c]['position']):
        chromosomes[c]['sites_by_samples'] = np.append(chromosomes[c]['sites_by_samples'],
            -np.ones((extend_amount, len(rows)), dtype="i1"), axis=0)
    chromosomes[c]['previous_position'] = -1
    if position is not None:
        if position in chromosomes[c]['position'] or position == chromosomes[c]['previous_position']:
            print("Trying to store more than one set of variants at position {} of chr {}. ".format(
                position, c))
            print("Previous was {}, this is {}.".format(
                position.get(position, chromosomes[c]['previous_position']), {variant.id:variant.alleles}))
            if position in chromosomes[c]['position']:
                chromosomes[c]['previous_position'], prev_data = chromosomes[c]['position'].popitem()
                print("Deleting original at position {} as well as duplicates.".format(
                    chromosomes[c]['previous_position']))
        elif position < chromosomes[c]['previous_position']:
            print("A position is out of order. We require a vcf file with all positions in strictly increasing order")
        else:
            chromosomes[c]['sites_by_samples'][len(chromosomes[c]['position']),:]=locus_data
            chromosomes[c]['position'][position]={variant.id: variant.alleles}
            chromosomes[c]['previous_position']=position
    
    if j % output_freq_variants == 0:
        print("{} variants read ({} saved). Base position {} Mb chr {} (alleles per site: {})".format(
            j+1, len(chromosomes[c]['position']), variant.pos/1e6, c, 
            [(k, allele_count[k]) for k in sorted(allele_count.keys())]), 
            flush=True)

#Finished reading

for c, dat in chromosomes.items():
    print("Processing chromosome {}".format(c))

    #check for samples with entirely missing data (e.g. if haploid)
    use_samples = np.ones((dat['sites_by_samples'].shape[1],), np.bool)
    for colnum in range(dat['sites_by_samples'].shape[1]):
        if np.all(dat['sites_by_samples'][:,colnum]==-1):
            use_samples[colnum]=False
    sites_by_samples = dat['sites_by_samples'][:,use_samples]
    reduced_rows = [name for i,name in enumerate(sorted(rows, key=rows.get)) if use_samples[i]]
    #simplify the names for haploids (remove e.g. terminal 'a' from the key 'XXXXa'
    # if there is not an XXXXb in the list)
    reduced_rows = [n if any((n[:-1]+suffix) in reduced_rows for suffix in suffixes if suffix != n[-1]) else n[:-2] for n in reduced_rows]
    print(" removed {}/{} unused sample slots (if this is haploid, half should be removed)".format(
        sum(~use_samples), use_samples.shape[0]))

    #check for missing data (all -1's for a site)
    use_sites = np.zeros((sites_by_samples.shape[0],), np.bool)
    for rownum in range(len(dat['position'])):
        if np.any(sites_by_samples[rownum]):
            use_sites[rownum] = True
    sites_by_samples = sites_by_samples[use_sites,:]
    print("Pruned down to {} sites".format(sites_by_samples.shape[0]))
    
    #check for samples with partial missing data
    use_samples = np.ones((sites_by_samples.shape[1],), np.bool)
    for colnum in range(sites_by_samples.shape[1]):
        if np.any(sites_by_samples[:,colnum]==-1):
            use_samples[colnum]=False
    sites_by_samples = sites_by_samples[:,use_samples]
    reduced_rows = [name for i,name in enumerate(reduced_rows) if use_samples[i]]
    print("Removed {}/{} samples which have incomplete data".format(
        sum(~use_samples), use_samples.shape[0]))


    outfile = args.outfile + str(c) + '.hdf5'
    with h5py.File(outfile, "w") as f:
        g = f.create_group("data")
        g.create_dataset("position", data=list(dat['position'].keys()))
        g.create_dataset("samples", data=[s.encode() for s in reduced_rows])
        g.create_dataset("variants", data=np.transpose(sites_by_samples), compression='gzip', compression_opts=9)
    print("Saved {} biallelic loci for {} samples into {}".format(len(dat['position']), len(reduced_rows), outfile))
