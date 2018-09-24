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

.PRECIOUS: %.trees

#############################################
# Standard pipeline for samples file to .trees
#############################################

%.ancestors: %.samples
	python3 -m tsinfer ga -vp -t ${NUM_THREADS} $^

%.ancestors.trees: %.ancestors
	python3 -m tsinfer ma -vp -t ${NUM_THREADS} $*.samples

%.nosimplify.trees: %.ancestors.trees
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} $*.samples -O $@

%.trees: %.nosimplify.trees
	python3 tsutil.py simplify $^ $@

%.trees.gz: %.trees
	gzip -k $^


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

chr%_ancestral_states.bcf: ancestral_states.vcf.gz ancestral_states.vcf.gz.tbi
	# Filter out multiallelic sites and indels to save time.
	bcftools view --max-alleles 2 --exclude-types indels ancestral_states.vcf.gz --regions $* -Ob -o $@


#############################################
# 1000 Genomes data.
#############################################

GENOTYPES_BASE=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files

1kg_samples.ped:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped \
		-o $@
1kg_%_genotypes.bcf:
	curl ${GENOTYPES_BASE}/ALL.$*.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf -o $@

1kg_%.samples: 1kg_%_genotypes.bcf.csi %_ancestral_states.bcf.csi 1kg_samples.ped
	python3 convert.py 1kg -p \
		1kg_$*_genotypes.bcf \
		$*_ancestral_states.bcf \
		-m 1kg_samples.ped \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@ 

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

#############################################
# UKBB
#############################################

ukbb_chr20_genotypes.bgen: 
	ln -s /gpfs2/well/ukbb-wtchg/v2/haplotypes/ukb_hap_chr20_v2.bgen $@

# Also requires the simplebgen module above. Don't want to require it here as 
# any changes will cascade a full rebuild on UKBB.
ukbb_chr20.samples: ukbb_chr20_genotypes.bgen chr20_ancestral_states.bcf.csi 
	python3 convert.py ukbb -p \
		ukbb_chr20_genotypes.bgen  \
		chr20_ancestral_states.bcf \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@ 

ukbb-chr20.augmented_131072.nosimplify.trees: ukbb_chr20.samples ukbb_chr20.ancestors.trees
	python3 tsutil.py sequential-augment ukbb_chr20.samples --num-threads=${NUM_THREADS}

clean:
	rm -f 1kg_samples.ped sgdp_samples.txt *.vcf* *.samples
