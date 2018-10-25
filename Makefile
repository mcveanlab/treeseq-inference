
FIGURES=\
	figures/storing_everyone.pdf\
	figures/metric_all_tools.pdf\
	figures/sample_edges.pdf\
	figures/metrics_all_tools_accuracy.pdf\
	figures/metrics_all_tools_demography.pdf\
	figures/metric_subsampling.pdf\
	figures/metric_all_tools_accuracy_sweep.pdf\
	figures/cputime_all_tools_by_sample_size.pdf\
	figures/fastarg_tsinfer_comparison_memory.pdf\
	figures/fastarg_tsinfer_comparison_time.pdf\
	figures/tsinfer_compression_ln.pdf\
	figures/tsinfer_edges_ln.pdf \
	figures/ukbb_structure.pdf

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

data/storing_everyone.csv: 
	python3 src/storing_everyone.py make-data

data/sample_edges.csv:
	python3 src/analyse_human_data.py sample_edges

figures/ukbb_structure.pdf: data/1kg_ukbb_british_centre.csv data/ukbb_ukbb_british_centre.csv
	python3 src/plot.py ukbb_structure

data/1kg_ukbb_british_centre.csv:
	python3 src/analyse_human_data.py 1kg_ukbb_gnn
	 
data/ukbb_ukbb_british_centre.csv:
	python3 src/analyse_human_data.py ukbb_ukbb_gnn

deps:
	make -C src
	# this should download and compile fastARG, ARGweaver, RentPlus, ftprime, etc
	make -C tools
