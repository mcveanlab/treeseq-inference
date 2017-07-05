#!/usr/bin/env python3
"""
Take a vcf file with ancestral information (e.g. from 1000 genomes phase 1, such as
ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
and convert it to a huge numpy array. Also output the variant positions and the sample names. Use as

vcf2tsinfer.py 1000G_chr22.vcf.gz outputarrays.npz
"""
import sys

import pysam
import numpy as np

vcf_in = pysam.VariantFile(sys.argv[1])

allele_count = {1:0,2:0,3:0,4:0}
rows, row = {}, 0
for sample_name in vcf_in.header.samples:
    for suffix in ('a','b'):
        rows[sample_name+suffix]=row
        row+=1
cols=0
positions = []
data = []
output_status_bases = 1e6 #how many bases between status outputs
curr_output_after = output_status_bases

print(len(vcf_in.header.records))

for rec in vcf_in.fetch():
    #restrict to cases where ancestral state contains only letters ATCG
    # i.e. where the ancestral state is certain (see ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README
    if "AA" in rec.info and all(letter in "ATCG" for letter in rec.info["AA"]):
        allele_count[len(rec.alleles)] += 1
        #only use biallelic variants
        if len(rec.alleles) == 2:
            if rec.info["AA"] not in rec.alleles:
                print("Ancestral state {} not in allele list ({})".format(\
                    rec.info["AA"], rec.alleles))
            else:
                column = np.zeros((len(rows),), dtype=np.bool)
                for label, samp in rec.samples.items():
                    for i,suffix in enumerate(('a','b')):
                        #print("allele {} (pos {}, sample {}, ancestral state {}, alleles {})".format( \
                        #    data.alleles[i], rec.pos, label+suffix, rec.info["AA"], rec.alleles))
                        column[rows[label+suffix]] = samp.alleles[i]!=rec.info["AA"]
                positions.append(rec.pos)
                data.append(column)
        if rec.pos > curr_output_after:
            print("@ base position {} Mb (alleles per site: {})".format(rec.pos/1e6, allele_count))
            while curr_output_after < rec.pos:
                curr_output_after += output_status_bases
np.savez(sys.argv[2], data=np.array(data), samples=np.array(sorted(rows, key= rows.get),dtype='|U{}'.format(max([len(n) for n in rows.keys()]))), positions=np.array(positions))