"""
Convert UK-biobank phased BGEN files to tsinfer sample input
"""

import csv
from datetime import datetime as dt
import argparse

import numpy as np
from tqdm import tqdm, trange

import tsinfer
from bgen_reader import read_bgen
"""We assume that although bgen_reader.read_bgen is no coded to work on phased BGEN files
the results stored in the ['genotypes'] array can be sensibly interpreted as phased
see https://github.com/limix/bgen-reader-py/issues/8
"""


def main():
    parser = argparse.ArgumentParser(
        description=__doc__)
    parser.add_argument("-chr", "--chromosome", type=int, default=19, help="Which chromosome number (integer)")
    args = parser.parse_args()
    chr = args.chromosome


    ancestral_alleles = {}
    for ref_alt_file, aa_allele_pos in {'ref':0, 'alt':1}.items():
        with open("/well/mcvean/pkalbers/ancestral/aa.bi.{}.snv.1000G.chr{}.txt".format(ref_alt_file, chr), newline='') as aa_file:
            print("Reading ancestral allele states for {} from {} (start time {})".format(ref_alt_file, aa_file.name, dt.now()))
            row_count = sum(1 for row in aa_file)
            dummy = aa_file.seek(0)
            aa_reader = csv.DictReader(aa_file, delimiter=' ')
            for row in tqdm(aa_reader, total=row_count):
                ancestral_alleles[row['rsID']]=row['Alleles'].split(",")[aa_allele_pos].upper()
    
    
    sample_data =  tsinfer.SampleData(path="chrom{}.samples".format(chr))
    bgen_file = "/well/ukbb-wtchg/v2/haplotypes/ukb_hap_chr{}_v2.bgen".format(chr)
    bgen = read_bgen(bgen_file, verbose=False, size = 200)
    sample_df = bgen['samples']
    print("Inputting sample names from {} (start time {})".format(chr, dt.now()))
    for i in tqdm(range(len(sample_df))):
        samp_name = sample_df.iat[i,0]
        sample_data.add_sample({'name':samp_name+"a"})
        sample_data.add_sample({'name':samp_name+"b"})
    
    
    
    
    ref_is_ancestral = alt_is_ancestral = 0
    v_not_in_aa = []
    a_not_in_aa = []
    variant_df = bgen['variants']
    derived_freq = {}
    print("Inputting haplotype data from chromosome {} (start time {})".format(chr, dt.now()))
    for idx in trange(len(bgen['genotype'])):
        if variant_df.at[idx,'nalleles'] == 2:
            rsID = variant_df.at[idx,'rsid']
            pos = variant_df.at[idx,'pos']
            alleles = variant_df.at[idx,'allele_ids'].upper().split(",")
            assert len(alleles) == 2
            if rsID not in ancestral_alleles:
                v_not_in_aa.append(rsID)
                try:
                    sample_data.add_site(pos, alleles, np.ravel(bgen['genotype'][idx,:,0:2].compute(num_workers=40).astype(np.uint8)), inference=False)
                except:
                    debug_file="debug.npy"
                    print("Error for pos {} (index {}) with alleles {}. Saving genotype data for debug to {}".format(pos, idx, alleles, debug_file))
                    np.save(debug_file, np.ravel(bgen['genotype'][idx,:,0:2].compute(num_workers=40).astype(np.uint8)))
                    raise
            elif ancestral_alleles[rsID] not in alleles:
                a_not_in_aa.append(rsID)
                sample_data.add_site(pos, alleles, np.ravel(bgen['genotype'][idx,:,0:2].compute(num_workers=40).astype(np.uint8)), inference=False)
            elif ancestral_alleles[rsID]==alleles[1]:
                sample_data.add_site(pos, alleles[::-1], 1-np.ravel(bgen['genotype'][idx,:,0:2].compute(num_workers=40).astype(np.uint8)))
                alt_is_ancestral += 1
                derived_freq[rsID] = np.mean(1-bgen['genotype'][idx,:,0:2].compute(num_workers=40))
            elif ancestral_alleles[rsID]==alleles[0]:
                sample_data.add_site(pos, alleles, np.ravel(bgen['genotype'][idx,:,0:2].compute(num_workers=40).astype(np.uint8)))
                ref_is_ancestral += 1
                derived_freq[rsID] = np.mean(bgen['genotype'][idx,:,0:2].compute(num_workers=40))
            else:
                print("Unspecified error for {} on chr{} on line {} of {}".format(rsID, chr, idx, bgen_file))
    
    sample_data.finalise()
    
    ff_name = 'UKBB_freqs_chr{}.csv'.format(chr)
    print("Writing frequencies to " + ff_name)
    with open(ff_name,'w') as f:
        for k,v in tqdm(derived_freq.items(), total=len(derived_freq)):
            print("{},{}".format(k,v))
    print("frequencies of derived alleles saved to freqs.csv")
    print("{} variants marked NOT for inference ({} with no variant of that name, {} with no matching ancestral allele)".format(
        len(v_not_in_aa) + len(a_not_in_aa), len(v_not_in_aa), len(a_not_in_aa)))
    print("{} variants marked for inference ({} with ref ancestral, {} with alt ancestral)".format(
        ref_is_ancestral + alt_is_ancestral, ref_is_ancestral, alt_is_ancestral))

if __name__ == "__main__":
    main()



"""

#To check frequencies

from datetime import datetime as dt
import csv

from tqdm import tqdm, trange

chr = 1
ancestral_alleles = {}
for ref_alt_file, aa_allele_pos in {'ref':0, 'alt':1}.items():
    with open("/well/mcvean/pkalbers/ancestral/aa.bi.{}.snv.1000G.chr{}.txt".format(ref_alt_file, chr)) as aa_file:
        print("Reading ancestral allele states for {} from {} (start time {})".format(ref_alt_file, aa_file.name, dt.now()))
        row_count = sum(1 for row in aa_file)
        dummy = aa_file.seek(0)
        aa_reader = csv.DictReader(aa_file, delimiter=' ')
        for row in tqdm(aa_reader, total=row_count):
            ancestral_alleles[row['rsID']]=row['Alleles'].split(",")[aa_allele_pos].upper()

derived_freq = {}
num_samples = 0
with open("/well/mcvean/ukbb12788/1000G/1000G.chr{}.marker.txt".format(chr)) as thousand_G:
    row_count = sum(1 for row in thousand_G)
    dummy = thousand_G.seek(0)
    tg_reader = csv.DictReader(thousand_G, delimiter=' ')
    for row in tqdm(tg_reader, total=row_count):
        rsID = row['Label']
        num_samples = max(num_samples, int(row['AlleleCount0']), int(row['AlleleCount0']), int(row['AlleleCountX']))
        alleles = row['Alleles'].split(",")
        if len(alleles) == 2:
            if rsID in ancestral_alleles and ancestral_alleles[rsID] in alleles:
                if ancestral_alleles[rsID] == alleles[0]:
                    derived_freq[rsID] = int(row['AlleleCount1'])
                elif ancestral_alleles[rsID] == alleles[1]:
                    derived_freq[rsID] = int(row['AlleleCount0'])
                else:
                    assert False

ff_name = '1000G_freqs_chr{}.csv'.format(chr)
print("Writing frequencies to " + ff_name)
with open(ff_name,'w') as f:
    for k,v in tqdm(derived_freq.items(), total=len(derived_freq)):
        print("{},{}".format(k,v/num_samples), file=f)
"""