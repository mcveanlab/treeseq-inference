
FIGURES=\
	figures/storing_everyone.pdf\
	figures/sample_edges.pdf


help:
	echo WRITE SOME HELP

.PRECIOUS: data/%.csv

figures/%.pdf: data/%.csv
	python3 src/plot.py $*

figs: ${FIGURES}


# Storing everyone plot
STORING_EVERYONE_PATH=data/raw__NOBACKUP__/storing_everyone
STORING_EVERYONE_VCF=${STORING_EVERYONE_PATH}/1000000.vcf.gz
STORING_EVERYONE_TREES=${STORING_EVERYONE_PATH}/1000000.trees

${STORING_EVERYONE_TREES}: 
	python3 src/storing_everyone.py simulate

${STORING_EVERYONE_VCF}: ${STORING_EVERYONE_TREES}
	python3 src/storing_everyone.py convert-files

data/storing_everyone.csv: ${STORING_EVERYONE_VCF}
	python3 src/storing_everyone.py make-data

deps:
	make -C src
	# this should download and compile fastARG, ARGweaver, RentPlus, ftprime, etc
	make -C tools
