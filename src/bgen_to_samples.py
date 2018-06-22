"""
Convert UK-biobank phased BGEN files to tsinfer sample input
"""
import numpy as np
import tsinfer
import bgen_reader
from tqdm import tqdm

for chr in range(20,21):
    sd =  tsinfer.load("1000G_GRCh38/H_sap_chr{}.samples".format(chr))
    bgen_file = "/well/ukbb-wtchg/v2/haplotypes/ukb_hap_chr{}_v2.bgen".format(chr)
    bgen = bgen_reader.read_bgen(bgen_file, verbose=False, size = 500)
    
    GRCh38_1000G = [s['ID'] for s in sd.sites_metadata[:]]
    GRCh37_UK_BB = [s for s in bgen['variants']['rsid']]
    GRCh37_UK_BB_info = {k['rsid']:[i,k['allele_ids'].split(",")] for i,k in bgen['variants'].iterrows()}
    
    s1000G = frozenset(GRCh38_1000G)
    sUK_BB = frozenset(GRCh37_UK_BB)
    intersection_order1000G = [x for x in GRCh38_1000G if x in sUK_BB]
    intersection_orderUK_BB = [x for x in GRCh37_UK_BB if x in s1000G]
    assert len(intersection_order1000G) == len(intersection_orderUK_BB)
    
    print("Chromosome {}, {} SNPs shared between 1000G GRCh38 ({} with aa) and UK-BB GRCh37 ({} total). Shared order = {}".format(
        chr, len(intersection_order1000G), len(GRCh38_1000G), len(GRCh37_UK_BB), 
        intersection_order1000G == intersection_orderUK_BB))
    
    with tsinfer.SampleData(path="UK_BB/UK-BB_chr{}.samples") as sample_data:
        sample_df = bgen['samples']
        for i in tqdm(range(len(sample_df)), desc="Read sample names"):
            samp_name = sample_df.iat[i,0]
            dummy = sample_data.add_individual(metadata={'name':samp_name}, ploidy=2)
        
        
        ref_is_ancestral = alt_is_ancestral = alleles_differ = nomatch = 0
        allele1_index = np.array([0], dtype=np.uint8) #the index of allele 1 in the UKBB allele list. 
        allele2_index = np.array([1], dtype=np.uint8) #the index of allele 2 in the UKBB allele list
        allele1_indices = np.array([0,2]) #the bgen indexes of the prob of allele1 in maternal & paternal chrs
        
        
        for pos, md, alleles in tqdm(
            zip(sd.sites_position[:], sd.sites_metadata[:], sd.sites_alleles[:]),
            desc="Read bgen", total=len(bgen['genotype'])):
            try:
                idx, UKBBalleles = GRCh37_UK_BB_info[md['ID']]
                if UKBBalleles==alleles:
                    ref_is_ancestral += 1
                    genotypes = np.where(
                        #NB: np.ravel does this in the right way, [[1,2],[3,4],[5,6]] -> [1,2,3,4,5,6]
                        np.ravel(bgen['genotype'][idx,:,:].compute(num_workers=40)[:,allele1_indices] > 0.5),
                        allele1_index, allele2_index)
                    dummy = sample_data.add_site(pos, genotypes, alleles)
                elif list(reversed(UKBBalleles))==alleles:
                    alt_is_ancestral += 1
                    genotypes = np.where(
                        np.ravel(bgen['genotype'][idx,:,:].compute(num_workers=40)[:,allele1_indices] > 0.5),
                        allele2_index, allele1_index) #reverse alt & ref 
                    dummy = sample_data.add_site(pos, genotypes, alleles)
                else:
                    alleles_differ += 1
                    print("Alleles differ {}, {}".format(UKBBalleles,alleles))
            except KeyError:
                #can omit any rsIDs not in GRCh37_UK_BB_index
                nomatch += 1
    
    print("Finished: {} with ref ancestral and {} with alt ancestral. {} with neither, and {} omitted".format(
        ref_is_ancestral, alt_is_ancestral, alleles_differ, nomatch))
