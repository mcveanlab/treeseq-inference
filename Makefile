
FIGURES=\
	figures/storing_everyone.pdf\
	figures/sample_edges.pdf


help:
	echo WRITE SOME HELP

.PRECIOUS: data/%.csv

# TODO add rules for the making the data files.

figures/%.pdf: data/%.csv
	python3 src/plot.py $*

figs: ${FIGURES}

deps:
	make -C src
	# this should download and compile fastARG, ARGweaver, RentPlus, ftprime, etc
	make -C tools
