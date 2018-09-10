
GENOTYPES_BASE=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/

all: 1kg_chr20.samples

1kg_samples.ped:
	curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped \
		-o $@

1kg_ancestral_states.vcf.gz.tbi:
	curl http://ftp.ensembl.org/pub/release-91/variation/vcf/homo_sapiens/homo_sapiens.vcf.gz.tbi -o $@

1kg_ancestral_states.vcf.gz:
	curl http://ftp.ensembl.org/pub/release-91/variation/vcf/homo_sapiens/homo_sapiens.vcf.gz -o $@

%.vcf.gz.tbi: %.vcf.gz
	tabix $(patsubst %.vcf.gz.tbi,%.vcf.gz,$@)

# chr20
1kg_chr20_ancestral_states.vcf.gz: 1kg_ancestral_states.vcf.gz 1kg_ancestral_states.vcf.gz.tbi
	bcftools view 1kg_ancestral_states.vcf.gz --regions 20 -Oz -o $@

1kg_chr20_genotypes.vcf.gz:
	curl ${GENOTYPES_BASE}ALL.chr20_GRCh38.genotypes.20170504.vcf.gz -o $@

1kg_chr20.samples: 1kg_chr20_genotypes.vcf.gz.tbi 1kg_chr20_ancestral_states.vcf.gz.tbi 1kg_samples.ped
	python3 convert_1kg.py -p \
		1kg_chr20_genotypes.vcf.gz  \
		1kg_chr20_ancestral_states.vcf.gz \
		1kg_samples.ped \
		$@

# chr1
1kg_chr1_ancestral_states.vcf.gz: 1kg_ancestral_states.vcf.gz 1kg_ancestral_states.vcf.gz.tbi
	bcftools view 1kg_ancestral_states.vcf.gz --regions 1 -Oz -o $@

1kg_chr1_genotypes.vcf.gz:
	curl ${GENOTYPES_BASE}ALL.chr1_GRCh38.genotypes.1170504.vcf.gz -o $@

1kg_chr1.samples: 1kg_chr1_genotypes.vcf.gz.tbi 1kg_chr1_ancestral_states.vcf.gz.tbi 1kg_samples.ped
	python3 convert_1kg.py -p \
		1kg_chr1_genotypes.vcf.gz  \
		1kg_chr1_ancestral_states.vcf.gz \
		1kg_samples.ped \
		$@


clean:
	rm -f 1kg_samples.ped *.vcf* *.samples
