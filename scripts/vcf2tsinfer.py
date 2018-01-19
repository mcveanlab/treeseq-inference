#!/usr/bin/env python3
"""
Take a vcf file with ancestral information (e.g. from 1000 genomes phase 1, such as
ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
and convert it to the huge samples X sites array used in tsinfer. Store this
in the hdf5 output file, and also store the variant positions and the sample names.

vcf2tsinfer.py 1000G_chr22.vcf.gz outputarray

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

parser = argparse.ArgumentParser(description='Take a vcf file with ancestral states' \
    ' and save it into a format suitable for run_vcf.py')
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
                            print("@ position {}{}, sample allele {} is not in {} - could be missing data. Omitting this site".format(
                                rec.pos, "+{}".format(allele_start) if allele_start else "", samp.alleles[i], rec.alleles))
                            return None
                        suffix = string.ascii_lowercase[i]
                        site_data[rows[label+suffix]] = samp.alleles[i][allele_start:]!=rec.info["AA"][allele_start:]
                return pos, site_data
    return None, None

allele_count = {}
rows, row = {}, 0
max_ploidy = 2
suffixes = string.ascii_lowercase[:max_ploidy] #use 'a' and 'b' for ploidy suffix

for sample_name in vcf_in.header.samples:
    #assume each sample can have a maximum of 2 variants (maternal & paternal)
    #i.e. max is diploid. For haploid samples, we just use the suffix 'a'
    #and hack to remove it afterwards
    for suffix in suffixes:
        rows[sample_name+suffix]=row
        row+=1
position = collections.OrderedDict()
extend_amount = 10000 #numpy storage array will be extended by this much at a time
sites_by_samples = -np.ones((extend_amount, len(rows)), dtype="i1")
output_freq_variants = 1e3 #output status after multiples of this many variants read


for j, variant in enumerate(vcf_in.fetch()):
    if j==args.only_use_n_variants:
        break
    #keep track of allele numbers
    allele_count[len(variant.alleles)] = allele_count.get(len(variant.alleles),0) + 1
    #check if we now need more storage
    if sites_by_samples.shape[0] <= len(position):
        sites_by_samples = np.append(sites_by_samples, np.zeros((extend_amount, len(rows)), dtype="i1"), axis=0)

    previous_position = -1
    actual_position, locus_data = process_variant(variant, rows, max_ploidy)
    if actual_position is not None:
        if actual_position in position or actual_position == previous_position:
            print("Trying to store more than one set of variants at position {}. ".format(actual_position))
            print("Previous was {}, this is {}.".format(
                position.get(actual_position, previous_position), {variant.id:variant.alleles}))
            if actual_position in position:
                previous_position, prev_data = position.popitem()
                print("Deleting original at position {} as well as duplicates.".format(previous_position))
        else:
            sites_by_samples[len(position),:]=locus_data
            position[actual_position]={variant.id: variant.alleles}
            previous_postition=actual_position
    
    if j % output_freq_variants == 0:
        print("{} variants read ({} saved). Base position {} Mb (alleles per site: {})".format(
            j+1, len(position), variant.pos/1e6, [(k, allele_count[k]) for k in sorted(allele_count.keys())]), 
            flush=True)

#Finished reading

#check for unused rows (e.g. if haploid)
use_samples = np.ones((sites_by_samples.shape[1],), np.bool)
for colnum in range(sites_by_samples.shape[1]):
    if np.all(sites_by_samples[:,colnum]==-1):
        use_samples[colnum]=False

sites_by_samples = sites_by_samples[:,use_samples]
rows = {name:order for name,order in rows.items() if use_samples[order]}
#simplify the names for haploids (remove e.g. terminal 'a' from the key 'XXXXa'
# if there is not an XXXXb in the list)
rows = {n:o if any((n[:-1]+suffix) in rows for suffix in suffixes if suffix != n[-1]) else n[:-1] for n,o in rows.items()}

#check for missing data (-1's)
use_sites = np.zeros((sites_by_samples.shape[0],), np.bool)
for rownum in range(len(position)):
    if np.all(sites_by_samples[rownum]!=-1):
        use_sites[rownum] = True
sites_by_samples = sites_by_samples[use_sites,:]


with h5py.File(args.outfile + '.hdf5', "w") as f:
    g = f.create_group("data")
    g.create_dataset("position", data=list(position.keys()))
    g.create_dataset("samples", data=[s.encode() for s in sorted(rows, key=rows.get)])
    g.create_dataset("variants", data=np.transpose(sites_by_samples), compression='gzip', compression_opts=9)
print("Saved {} biallelic loci for {} samples into {}".format(len(position), len(rows), args.outfile))
