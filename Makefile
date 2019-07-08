
FIGURES=\
	figures/storing_everyone.pdf\
	figures/metric_all_tools.pdf\
	figures/metrics_all_tools_accuracy.pdf\
	figures/metric_all_tools_accuracy_demography.pdf\
	figures/metric_subsampling.pdf\
	figures/metric_all_tools_accuracy_sweep.pdf\
	figures/cputime_all_tools_by_sample_size.pdf\
	figures/mem_time_fastarg_tsinfer.pdf\
	figures/tsinfer_vcf_compression_ln.pdf\
	figures/tsinfer_ts_filesize_ln.pdf\
	figures/ukbb_structure.pdf \
	figures/sample_edges_1kg.pdf \
	figures/global_structure.pdf 

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

data/1kg_ukbb_british_centre.csv:
	python3 src/analyse_human_data.py 1kg_ukbb_gnn
	 
data/ukbb_ukbb_british_centre.csv:
	python3 src/analyse_human_data.py ukbb_ukbb_gnn

data/1kg_gnn.csv:
	cp human-data/1kg_chr20.snipped.trees.gnn.csv $@

data/sgdp_gnn.csv:
	cp human-data/sgdp_chr20.snipped.trees.gnn.csv $@

data/HG01933_local_gnn.csv:
	python3 src/analyse_human_data.py hg01933_local_gnn

figures/sample_edges_1kg.pdf: data/sample_edges.csv
	python3 src/plot.py sample_edges

figures/ukbb_structure.pdf: data/1kg_ukbb_british_centre.csv data/ukbb_ukbb_british_centre.csv
	python3 src/plot.py ukbb_structure

figures/global_structure.pdf: data/sample_edges.csv data/1kg_gnn.csv data/HG01933_local_gnn_0.csv
	python3 src/plot.py global_structure

deps:
	make -C src
	# this should download and compile fastARG, ARGweaver, RentPlus, ftprime, etc
	make -C tools
