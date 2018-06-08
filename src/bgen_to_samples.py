#assumes that although bgen_reader.read_bgen is no coded to work on phased BGEN files
# the results stored in the ['genotypes'] array can be sensibly interpreted as phased
# see https://github.com/limix/bgen-reader-py/issues/8

import csv
from datetime import datetime as dt

import numpy as np
from tqdm import tqdm

import tsinfer
from bgen_reader import read_bgen

chr = 20

ancestral_alleles = {}
for ref_alt_file, aa_allele_pos in {'ref':0, 'alt':1}.items():
    with open("/well/mcvean/pkalbers/ancestral/aa.bi.{}.snv.1000G.chr{}.txt".format(ref_alt_file, chr), newline='') as aa_file:
        print("Reading ancestral allele states for {} from {} (start time {})".format(ref_alt_file, aa_file.name, dt.now()))
        row_count = sum(1 for row in aa_file)
        dummy = aa_file.seek(0)
        aa_reader = csv.DictReader(aa_file, delimiter=' ')
        for row in tqdm(aa_reader, total=row_count):
            ancestral_alleles[row['rsID']]=row['Alleles'].split(",")[aa_allele_pos].upper()


sample_data =  tsinfer.SampleData(path="chr{}.samples".format(chr))
bgen_file = "/well/ukbb-wtchg/v2/haplotypes/ukb_hap_chr{}_v2.bgen".format(chr)
bgen = read_bgen(bgen_file, verbose=False, metadata_file=False)

sample_df = bgen['samples']
print("Inputting sample names from {} (start time {})".format(chr, dt.now()))
for i in tqdm(range(len(sample_df))):
    samp_name = sample_df.iat[i,0]
    sample_data.add_sample({'name':samp_name+"a"})
    sample_data.add_sample({'name':samp_name+"b"})


def chunkedRange(start,end,step):
    i = start
    j = i+step
    while j < end:
        yield (i,j)
        i, j = j, j+step
    yield (i,end)

n_var = len(bgen['genotype'])

print("Inputting haplotype data from {} (start time {})".format(chr, dt.now()))

batch_size = 1000
ref_is_ancestral = alt_is_ancestral = 0
v_not_in_aa = []
a_not_in_aa = []
variant_df = bgen['variants']
with tqdm(total=n_var) as pbar:
    for start,end in chunkedRange(0,len(bgen['genotype']), batch_size):
        genotypes = bgen['genotype'][start:end,:,0:2].compute(num_workers=40).astype(np.uint8)
        for local_idx, idx in enumerate(range(start,end)): #genotypes is a subsetted np array which need local indexing
            if variant_df.at[idx,'nalleles'] == 2:
                rsID = variant_df.at[idx,'rsid']
                pos = variant_df.at[idx,'pos']
                alleles = variant_df.at[idx,'allele_ids'].upper().split(",")
                assert len(alleles) == 2
                if rsID not in ancestral_alleles:
                    v_not_in_aa.append(rsID)
                    sample_data.add_site(pos, alleles, np.ravel(genotypes[local_idx,:,:]), inference=False)
                elif ancestral_alleles[rsID] not in alleles:
                    a_not_in_aa.append(rsID)
                    sample_data.add_site(pos, alleles, np.ravel(genotypes[local_idx,:,:]), inference=False)
                elif ancestral_alleles[rsID]==alleles[1]:
                    sample_data.add_site(pos, alleles[::-1], 1-np.ravel(genotypes[local_idx,:,:]))
                    alt_is_ancestral += 1
                elif ancestral_alleles[rsID]==alleles[0]:
                    sample_data.add_site(pos, alleles, np.ravel(genotypes[local_idx,:,:]))
                    ref_is_ancestral += 1
                else:
                    print("Unspecified error for {} on chr{} on line {} of {}".format(rsID, chr, idx, bgen_file))
                
            pbar.update()

sample_data.finalise()

print("{} variants marked NOT for inference ({} with no variant of that name, {} with no matching ancestral allele)".format(
    len(v_not_in_aa), len(a_not_in_aa))
print("{} variants marked for inference ({} with ref ancestral, {} with alt ancestral)".format(
    ref_is_ancestral + alt_is_ancestral, ref_is_ancestral, alt_is_ancestral))
