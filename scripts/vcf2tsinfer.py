#!/usr/bin/env python3
"""
Take a vcf file with ancestral information (e.g. from 1000 genomes phase 1, such as
ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf
and convert it to a huge numpy array. Also output the variant positions and the sample names
"""


import pysam
import numpy as np

vcf_in = pysam.VariantFile("/Volumes/SDdisk/1000G/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz")

allele_count = {0:0,1:0,2:0,3:0,4:0}
rows, row = {}, 0
for sample_name in vcf_in.header.samples:
    for suffix in ('a','b'):
        rows[sample_name+suffix]=row
        row+=1
cols=0
positions = []


for rec in vcf_in.fetch():
    #only print SNPs where the ancestral state is certain (see ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README
    if "AA" in rec.info and all(letter in "ATCG" for letter in rec.info["AA"]):
        #only use biallelic variants
        allele_count[len(rec.alleles)] += 1
        if len(rec.alleles) == 2:
            if rec.info["AA"] not in rec.alleles:
                print("alleles: {}, ancestral state: {} (currently at saved variant {})".format(rec.alleles, rec.info["AA"], cols))
            else:
                column = np.zeros((len(rows),), dtype=np.bool)
                for label, data in rec.samples.items():
                    for i,suffix in enumerate(('a','b')):
                        #print("allele {} (pos {}, sample {}, ancestral state {}, alleles {})".format( \
                        #    data.alleles[i], rec.pos, label+suffix, rec.info["AA"], rec.alleles))
                        column[rows[label+suffix]] = data.alleles[i]!=rec.info["AA"]
                positions.append(rec.pos)
                cols+=1
print(cols)