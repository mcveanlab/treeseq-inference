NUM_THREADS ?= 0

# Requirements: bcftools, tabix, python3.
# See requirements.txt for Python package requirements.
#
help:
	@echo WRITE SOME HELP

all: 1kg_chr20.samples sgdp_chr20.samples


simplebgen: setup.py simplebgenmodule.c
	python3 setup.py build_ext --inplace

%.bcf.csi: %.bcf
	bcftools index $(patsubst %.bcf.csi,%.bcf,$@)

%.vcf.gz.csi: %.vcf.gz
	bcftools index $(patsubst %.vcf.gz.csi,%.vcf.gz,$@)

# Any implicit targets that are built will be deleted unless we make them as 'precious'.
.PRECIOUS: \
	%.ancestors %.ancestors.trees %.nosimplify.trees %.trees %.trees.gz\
	%.bcf.csi %.vcf.gz.csi \
	chr%_ancestral_states.fa \
	1kg_%_genotypes.bcf 1kg_%.samples

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

ANCESTRAL_STATES_PREFIX=homo_sapiens_ancestor_GRCh37_e71
ANCESTRAL_STATES_TARBALL=${ANCESTRAL_STATES_PREFIX}.tar.bz2
ANCESTRAL_STATES_URL=ftp://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/${ANCESTRAL_STATES_TARBALL}

${ANCESTRAL_STATES_TARBALL}:
	curl ${ANCESTRAL_STATES_URL} -o ${ANCESTRAL_STATES_TARBALL}

${ANCESTRAL_STATES_PREFIX}/README: ${ANCESTRAL_STATES_TARBALL}
	rm -fR ${ANCESTRAL_STATES_PREFIX}
	tar -jxvf ${ANCESTRAL_STATES_TARBALL}
	# Update access times or we'll keep rebuilding this rule
	touch ${ANCESTRAL_STATES_PREFIX}/*

chr%_ancestral_states.fa: ${ANCESTRAL_STATES_PREFIX}/README
	ln -f ${ANCESTRAL_STATES_PREFIX}/homo_sapiens_ancestor_$*.fa $@


#############################################
# 1000 Genomes data.
#############################################

GENOTYPES_BASE=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files

1kg_samples.ped:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped \
		-o $@
1kg_%_genotypes.bcf:
	curl ${GENOTYPES_BASE}/ALL.$*.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf -o $@

1kg_%.samples: 1kg_%_genotypes.bcf.csi %_ancestral_states.fa 1kg_samples.ped
	python3 convert.py 1kg -p \
		1kg_$*_genotypes.bcf \
		$*_ancestral_states.fa \
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

# TODO we need to make sure that the source VCF is indexed here or the sample
# removal will fail.
#
sgdp_chr20_genotypes.bcf: cteam_extended.v4.PS2_phase.public.chr20.vcf.gz 
	# Remove the S_Naxi-2 individual because (a) it doesn't have any metadata in the 
	# file we're using and (b) it has a massively elevated sample edge count if we 
	# leave it in.
	bcftools view -s '^S_Naxi-2' $^ -O b -o $@

sgdp_chr20.samples: sgdp_chr20_genotypes.bcf.csi chr20_ancestral_states.fa sgdp_samples.txt
	python3 convert.py sgdp -p \
		sgdp_chr20_genotypes.bcf \
		chr20_ancestral_states.fa \
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
ukbb_chr20.samples: ukbb_chr20_genotypes.bgen chr20_ancestral_states.fa
	python3 convert.py ukbb -p \
		ukbb_chr20_genotypes.bgen  \
		chr20_ancestral_states.fa \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@ 

ukbb-chr20.augmented_131072.nosimplify.trees: ukbb_chr20.samples ukbb_chr20.ancestors.trees
	python3 tsutil.py sequential-augment ukbb_chr20.samples --num-threads=${NUM_THREADS}

clean:
	rm -f 1kg_samples.ped sgdp_samples.txt *.vcf* *.samples
