NUM_THREADS ?= 0

# Requirements: bcftools, tabix, python3.
# See requirements.txt for Python package requirements.

all: 1kg_chr20.samples sgdp_chr20.samples

simplebgen: setup.py simplebgenmodule.c
	python3 setup.py build_ext --inplace

%.bcf.csi: %.bcf
	bcftools index $(patsubst %.bcf.csi,%.bcf,$@)

%.vcf.gz.csi: %.vcf.gz
	bcftools index $(patsubst %.vcf.gz.csi,%.vcf.gz,$@)

#############################################
# Ancestral states from Ensembl
#############################################

# SGDP and 1000G are aligned to hs37d5, which is GRCh37 plus extra decoy sequences. 
# So we download the ancestral states for GRCh37. 

# Recorded in the sample file provenance.
REFERENCE_NAME=GRCh37

# We're using Ensembl release 92. 
# I tried using 93 but the VCF header said it was release 92.
#
ANCESTRAL_STATES_URL=ftp://ftp.ensembl.org/pub/grch37/release-92/variation/vcf/homo_sapiens/homo_sapiens.vcf.gz 

ancestral_states.vcf.gz.tbi:
	curl ${ANCESTRAL_STATES_URL}.tbi -o $@

ancestral_states.vcf.gz:
	curl ${ANCESTRAL_STATES_URL} -o $@

# chr20
chr20_ancestral_states.bcf: ancestral_states.vcf.gz ancestral_states.vcf.gz.tbi
	# Filter out multiallelic sites and indels to save time.
	bcftools view --max-alleles 2 --exclude-types indels ancestral_states.vcf.gz --regions 20 -Ob -o $@

# chr1
chr1_ancestral_states.bcf: ancestral_states.vcf.gz ancestral_states.vcf.gz.tbi
	# Filter out multiallelic sites and indels to save time.
	bcftools view --max-alleles 2 --exclude-types indels ancestral_states.vcf.gz --regions 1 -Ob -o $@


#############################################
# 1000 Genomes data.
#############################################

GENOTYPES_BASE=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/

1kg_samples.ped:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped \
		-o $@
# chr20
1kg_chr20_genotypes.bcf:
	curl ${GENOTYPES_BASE}/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf -o $@

1kg_chr20.samples: 1kg_chr20_genotypes.bcf.csi chr20_ancestral_states.bcf.csi 1kg_samples.ped
	python3 convert.py 1kg -p \
		1kg_chr20_genotypes.bcf \
		chr20_ancestral_states.bcf \
		-m 1kg_samples.ped \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@ 

1kg_chr20.ancestors: 1kg_chr20.samples
	python3 -m tsinfer ga -vp -t ${NUM_THREADS} $^

1kg_chr20.ancestors.trees: 1kg_chr20.ancestors
	python3 -m tsinfer ma -vp -t ${NUM_THREADS} 1kg_chr20.samples

1kg_chr20.nosimplify.trees: 1kg_chr20.ancestors.trees
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} 1kg_chr20.samples -O $@

# chr1
1kg_chr1_genotypes.bcf:
	curl ${GENOTYPES_BASE}/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf -o $@

1kg_chr1.samples: 1kg_chr1_genotypes.bcf.csi chr1_ancestral_states.bcf.csi 1kg_samples.ped
	python3 convert.py 1kg -p \
		1kg_chr1_genotypes.bcf \
		chr1_ancestral_states.bcf \
		-m 1kg_samples.ped \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@ 

1kg_chr1.ancestors: 1kg_chr1.samples
	python3 -m tsinfer ga -vp -t ${NUM_THREADS} $^

1kg_chr1.ancestors.trees: 1kg_chr1.ancestors
	python3 -m tsinfer ma -vp -t ${NUM_THREADS} 1kg_chr1.samples

1kg_chr1.nosimplify.trees: 1kg_chr1.ancestors.trees
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} 1kg_chr1.samples -O $@

#############################################
# SGDP data.
#############################################

sgdp_samples.txt:
	curl https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.samples.txt -o $@
	 
# TODO download the full tarball for the cteam_extended files?

sgdp_chr20_genotypes.vcf.gz: cteam_extended.v4.PS2_phase.public.chr20.vcf.gz
	ln cteam_extended.v4.PS2_phase.public.chr20.vcf.gz $@

sgdp_chr20.samples: sgdp_chr20_genotypes.vcf.gz.csi chr20_ancestral_states.bcf.csi sgdp_samples.txt
	python3 convert.py sgdp -p \
		sgdp_chr20_genotypes.vcf.gz  \
		chr20_ancestral_states.bcf \
		-m sgdp_samples.txt \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@ 

sgdp_chr20.ancestors: sgdp_chr20.samples
	python3 -m tsinfer ga -vp -t ${NUM_THREADS} $^

sgdp_chr20.ancestors.trees: sgdp_chr20.ancestors
	python3 -m tsinfer ma -vp -t ${NUM_THREADS} sgdp_chr20.samples

sgdp_chr20.nosimplify.trees: sgdp_chr20.ancestors.trees
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} sgdp_chr20.samples -O $@


#############################################
# UKBB
#############################################

ukbb_chr20_genotypes.bgen: 
	ln -s /gpfs2/well/ukbb-wtchg/v2/haplotypes/ukb_hap_chr20_v2.bgen $@

ukbb_chr20.samples: simplebgen ukbb_chr20_genotypes.bgen chr20_ancestral_states.bcf.csi 
	python3 convert.py ukbb -p \
		ukbb_chr20_genotypes.bgen  \
		chr20_ancestral_states.bcf \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@ 
		 
ukbb_chr20.ancestors: ukbb_chr20.samples
	python3 -m tsinfer ga -vp -t ${NUM_THREADS} $^

ukbb_chr20.ancestors.trees: ukbb_chr20.ancestors
	python3 -m tsinfer ma -vp -t ${NUM_THREADS} ukbb_chr20.samples

ukbb_chr20.nosimplify.trees: ukbb_chr20.ancestors.trees
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} ukbb_chr20.samples -O $@



clean:
	rm -f 1kg_samples.ped sgdp_samples.txt *.vcf* *.samples
