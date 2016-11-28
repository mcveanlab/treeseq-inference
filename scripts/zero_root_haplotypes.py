#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage: ../scripts/zero_root_haplotypes.py 4000.haps 4000.fastarg 4000.haps_zeroed
"""

#find line beginning with S
import argparse
parser = argparse.ArgumentParser(description='Convert haplotype sequences so that all mutations are from 0->1, using a base sequence taken from a fastARG output file')
parser.add_argument('seqfile', type=argparse.FileType('r', encoding='UTF-8'), help='a set of sequences (one line per locus) as using in fastARG input')
parser.add_argument('fastARGfile', type=argparse.FileType('r', encoding='UTF-8'), help='an output file from fastARG (https://github.com/lh3/fastARG ), containing a root sequence line beginning with S')
parser.add_argument('outfile',  type=argparse.FileType('w', encoding='UTF-8'), help='an output file of the same format as the seqfile input')
args = parser.parse_args()

for line in args.fastARGfile:
    spl = line.split(None,3)
    if spl[0]=="S":
        base_seq = spl[2]
        print("sequence length: {}".format(len(base_seq)))
        break
    
for i, line in enumerate(args.seqfile):    
    spl = line.split(None,2)
    args.outfile.write(spl[0])
    args.outfile.write("\t")
    if base_seq[i] == '1':
        args.outfile.write(spl[1].translate(str.maketrans('01','10')))
    else:
        args.outfile.write(spl[1])
    args.outfile.write("\n")
