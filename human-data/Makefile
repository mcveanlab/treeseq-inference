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
	%.ancestors %.ancestors.trees %.nosimplify.trees %.nosimplify.nopc.trees %.trees %.trees.gz\
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
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} $*.samples -O $@ --no-simplify

%.nosimplify.nopc.trees: %.ancestors.trees
	python3 -m tsinfer ms -vp -t ${NUM_THREADS} $*.samples -O $@ --no-simplify --no-path-compression

%.trees: %.nosimplify.trees
	python3 tsutil.py simplify $^ $@

%.trees.gz: %.trees
	gzip -c $^ > $@

%.trees.tsz: %.trees
	tszip -k $^ 

%.trees.bcf: %.trees
	msp vcf -P 2 $^ | bcftools view - -O b -o $@

%.snipped.trees: %.trees ${CENTROMERES_CSV}
	python3 tsutil.py snip-centromere $< $@ $* ${CENTROMERES_CSV}

#############################################
# Centromere locations for GRCh37 (aka hg19) from UCSC
# See https://www.biostars.org/p/2349/
#############################################
CENTROMERES_CSV=centromeres.csv
${CENTROMERES_CSV}:
	 curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz > cytoband.txt.gz
	 echo "chrom,start,end" > ${CENTROMERES_CSV}
	 # Start and end coordinates are on different lines, so we merge them.
	 zcat cytoband.txt.gz | grep acen | sort | paste -d " " - - \
		 | cut -f 1,2,7 --output-delim="," >> ${CENTROMERES_CSV}

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
	# Update access times or we'll keep rebuilding this rule. Have to make sure 
	# that the README we touch is older than the actual fa files.
	touch $@
	touch ${ANCESTRAL_STATES_PREFIX}/*.fa

chr%_ancestral_states.fa: ${ANCESTRAL_STATES_PREFIX}/README
	ln -sf ${ANCESTRAL_STATES_PREFIX}/homo_sapiens_ancestor_$*.fa $@


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
		$@  > $@.report

1kg_chr20.snipped.trees.gnn.csv: 1kg_chr20.snipped.trees
	python3 tsutil.py compute-1kg-gnn $^ $@ --num-threads=16


#############################################
# SGDP data.
#############################################


SGDP_GENOTYPES_BASE=https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data/PS2_multisample_public

sgdp_samples.txt:
	curl https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.samples.txt -o $@

sgdp_%_genotypes.vcf.gz:
	curl ${SGDP_GENOTYPES_BASE}/cteam_extended.v4.PS2_phase.public.$*.vcf.gz -o $@
	curl ${SGDP_GENOTYPES_BASE}/cteam_extended.v4.PS2_phase.public.$*.vcf.gz.csi -o $@.csi

sgdp_%_genotypes.bcf: sgdp_%_genotypes.vcf.gz
	# Remove the S_Naxi-2 individual because (a) it doesn't have any metadata in the 
	# file we're using and (b) it has a massively elevated sample edge count if we 
	# leave it in.
	bcftools view -s '^S_Naxi-2' $^ -O b -o $@

sgdp_%.samples: sgdp_%_genotypes.bcf.csi %_ancestral_states.fa sgdp_samples.txt
	python3 convert.py sgdp -p \
		sgdp_$*_genotypes.bcf \
		$*_ancestral_states.fa \
		-m sgdp_samples.txt \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@  > $@.report

sgdp_%.snipped.trees.gnn.csv: sgdp_%.snipped.trees
	python3 tsutil.py compute-sgdp-gnn $^ $@ --num-threads=16

#############################################
# UKBB
#############################################

ukbb_chr20_genotypes.bgen: 
	ln -s /gpfs2/well/ukbb-wtchg/v2/haplotypes/ukb_hap_chr20_v2.bgen $@

# Also requires the simplebgen module above. Don't want to require it here as 
# any changes will cascade a full rebuild on UKBB.
ukbb_chr20.samples: ukbb_chr20_genotypes.bgen chr20_ancestral_states.fa ukbb_metadata.csv
	python3 convert.py ukbb -p \
		ukbb_chr20_genotypes.bgen  \
		chr20_ancestral_states.fa \
		-m ukbb_metadata.csv \
		--ancestral-states-url=${ANCESTRAL_STATES_URL}\
		--reference-name=${REFERENCE_NAME}\
		$@ > $@.report

# TODO change to ukbb_chr20
ukbb_chr20.augmented_131072.nosimplify.trees: ukbb_chr20.samples ukbb_chr20.ancestors.trees
	python3 tsutil.py sequential-augment ukbb_chr20.samples --num-threads=${NUM_THREADS}

1kg_ukbb_chr20.samples 1kg_ukbb_chr20.ancestors.trees: ukbb_chr20.samples 1kg_chr20.trees
	# FIXME: to make things run a bit more quickly use 100K samples. If we try
	# any less, we end up having no variation at some sites.
	#python3 tsutil.py combine-ukbb-1kg chr20 --num-individuals=50000
	python3 tsutil.py combine-ukbb-1kg chr20 

1kg_ukbb_chr20.snipped.trees: 1kg_ukbb_chr20.trees ${CENTROMERES_CSV}
	python3 tsutil.py snip-centromere $< $@ chr20 ${CENTROMERES_CSV}
	 
1kg_ukbb_chr20.snipped.trees.gnn.csv: 1kg_ukbb_chr20.snipped.trees 
	python3 tsutil.py compute-1kg-ukbb-gnn $< $@

ukbb_chr20.augmented_131072.snipped.trees: ukbb_chr20.augmented_131072.trees ${CENTROMERES_CSV}
	python3 tsutil.py snip-centromere $< $@ chr20 ${CENTROMERES_CSV}
	 
ukbb_chr20.augmented_131072.snipped.trees.gnn.csv: ukbb_chr20.augmented_131072.snipped.trees
	python3 tsutil.py compute-ukbb-gnn $< $@

clean:
	rm -f 1kg_samples.ped sgdp_samples.txt *.vcf* *.samples
