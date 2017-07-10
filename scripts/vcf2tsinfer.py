#!/usr/bin/env python3
"""
Take a vcf file with ancestral information (e.g. from 1000 genomes phase 1, such as
ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
and convert it to the huge samples X sites array used in tsinfer. Store this
in the hdf5 output file, and also store the variant positions and the sample names.

vcf2tsinfer.py 1000G_chr22.vcf.gz outputarrays.hdf5

NB - a slightly odd observation: all the alleles for which the ancestral state is not present in the dataset are SNPs

read using
with h5py.File(filename, 'r') as f:
    print(f['data']['variants'][()])
"""
import sys
import collections

import numpy as np
import h5py
import pysam
#To do - use argparse module 

vcf_in = pysam.VariantFile(sys.argv[1])

def process_variant(rec, data):
    """
    return True only when
    """
    data['allele_count'][len(rec.alleles)] = data['allele_count'].get(len(rec.alleles),0) + 1
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
                    pos = rec.pos+allele_start-0.5
                    #print("Starting allele at an incremented position ({} not {}) for {} (alleles:{})".format(
                    #    pos, rec.pos, rec.id, rec.alleles))
                else:
                    allele_start=0
                    pos = rec.pos

                if pos in data['position']:                    
                    print("More than one set of variants at position {}. ".format(pos))
                    print("Previous was {}. Omitting a subsequent duplicate (id {})".format(
                        data['position'][pos], {rec.id:rec.alleles}))
                    return False

                column = np.zeros((len(data['rows']),), dtype="i1")
                for label, samp in rec.samples.items():
                    for i,suffix in enumerate(('a','b')):
                        #print("allele {} (pos {}, sample {}, ancestral state {}, alleles {})".format( \
                        #    rec.alleles[i], rec.pos, label+suffix, rec.info["AA"], rec.alleles))
                        if samp.alleles[i] not in rec.alleles:
                            print("@ position {}{}, sample allele {} is not in {} - could be missing data. Omitting this row".format(
                                rec.pos, "+{}".format(allele_start) if allele_start else "", samp.alleles[i], rec.alleles))
                            return False
                        column[data['rows'][label+suffix]] = samp.alleles[i][allele_start:]!=rec.info["AA"][allele_start:]
                if data['sites_by_samples'].shape[0] <= len(data['position']):
                    #need more storage
                    data['sites_by_samples'] = np.append(data['sites_by_samples'], np.zeros((data['extend_amount'], len(data['rows'])), dtype="i1"), axis=0)
                data['sites_by_samples'][len(data['position'])]=column
                data['position'][pos]={rec.id:rec.alleles}
                return True
    return False

store={}
store['allele_count'] = {}
store['rows'], row = {}, 0
for sample_name in vcf_in.header.samples:
    for suffix in ('a','b'):
        store['rows'][sample_name+suffix]=row
        row+=1
store['position'] = collections.OrderedDict()
store['extend_amount'] = 10000 #numpy storage array will be extended by this much at a time
store['sites_by_samples'] = np.zeros((store['extend_amount'], len(store['rows'])), dtype="i1")
output_freq_variants = 1e3 #output status after multiples of this many variants read

for j, variant in enumerate(vcf_in.fetch()):
    process_variant(variant, store)
    if j % output_freq_variants == 0:
        print("{} variants read ({} saved). Base position {} Mb (alleles per site: {})".format(
            j+1, len(store['position']), variant.pos/1e6, [(k, store['allele_count'][k]) for k in sorted(store['allele_count'].keys())]), 
            flush=True)

with h5py.File(sys.argv[2], "w") as f:
    g = f.create_group("data")
    g.create_dataset("position", data=position.keys())
    g.create_dataset("samples", data=[s.encode() for s in sorted(store['rows'], key= rows.get)])
    g.create_dataset("variants", data=np.transpose(sites_by_samples[:len(store['position'])]))
print("Saved {} biallelic loci for {} samples into {}".format(len(position), len(rows), sys.argv[2]))
